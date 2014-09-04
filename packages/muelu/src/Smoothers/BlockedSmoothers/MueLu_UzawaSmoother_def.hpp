// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
/*
 * MueLu_UzawaSmoother_def.hpp
 *
 *  Created on: 13 May 2014
 *      Author: wiesner
 */

#ifndef MUELU_UZAWASMOOTHER_DEF_HPP_
#define MUELU_UZAWASMOOTHER_DEF_HPP_

#include "Teuchos_ArrayViewDecl.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "MueLu_ConfigDefs.hpp"

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_UzawaSmoother_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_HierarchyHelpers.hpp"
#include "MueLu_SmootherBase.hpp"
#include "MueLu_SubBlockAFactory.hpp"

// include files for default FactoryManager
#include "MueLu_SchurComplementFactory.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_FactoryManager.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  UzawaSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::UzawaSmoother()
    : type_("Uzawa"), A_(Teuchos::null)
  {
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  UzawaSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~UzawaSmoother() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const ParameterList> UzawaSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< RCP<const FactoryBase> >("A",                  Teuchos::null, "Generating factory of the matrix A");
    validParamList->set< Scalar >                ("Damping factor",     1.0, "Damping/Scaling factor in SIMPLE");
    validParamList->set< LocalOrdinal >          ("Sweeps",             1, "Number of SIMPLE sweeps (default = 1)");

    return validParamList;
  }

  //! Add a factory manager at a specific position
  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void UzawaSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AddFactoryManager(RCP<const FactoryManagerBase> FactManager, int pos) {
    TEUCHOS_TEST_FOR_EXCEPTION(pos < 0, Exceptions::RuntimeError, "MueLu::UzawaSmoother::AddFactoryManager: parameter \'pos\' must not be negative! error.");

    size_t myPos = Teuchos::as<size_t>(pos);

    if (myPos < FactManager_.size()) {
      // replace existing entris in FactManager_ vector
      FactManager_.at(myPos) = FactManager;
    } else if( myPos == FactManager_.size()) {
      // add new Factory manager in the end of the vector
      FactManager_.push_back(FactManager);
    } else { // if(myPos > FactManager_.size())
      RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      *out << "Warning: cannot add new FactoryManager at proper position " << pos << ". The FactoryManager is just appended to the end. Check this!" << std::endl;

      // add new Factory manager in the end of the vector
      FactManager_.push_back(FactManager);
    }

  }


  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void UzawaSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetVelocityPredictionFactoryManager(RCP<FactoryManager> FactManager) {
    AddFactoryManager(FactManager, 0); // overwrite factory manager for predicting the primary variable
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void UzawaSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetSchurCompFactoryManager(RCP<FactoryManager> FactManager) {
    AddFactoryManager(FactManager, 1); // overwrite factory manager for SchurComplement
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void UzawaSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    currentLevel.DeclareInput("A",this->GetFactory("A").get());

    TEUCHOS_TEST_FOR_EXCEPTION(FactManager_.size() != 2, Exceptions::RuntimeError,"MueLu::UzawaSmoother::DeclareInput: You have to declare two FactoryManagers with a \"Smoother\" object: One for predicting the primary variable and one for the SchurComplement system. The smoother for the SchurComplement system needs a SchurComplementFactory as input for variable \"A\". make sure that you use the same proper damping factors for omega both in the SchurComplementFactory and in the SIMPLE smoother!");

    // loop over all factory managers for the subblocks of blocked operator A
    std::vector<Teuchos::RCP<const FactoryManagerBase> >::const_iterator it;
    for(it = FactManager_.begin(); it!=FactManager_.end(); ++it) {
      SetFactoryManager currentSFM  (rcpFromRef(currentLevel),   *it);

      // request "Smoother" for current subblock row.
      currentLevel.DeclareInput("PreSmoother",(*it)->GetFactory("Smoother").get());
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void UzawaSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Setup(Level &currentLevel) {

    FactoryMonitor m(*this, "Setup blocked Uzawa Smoother", currentLevel);

    if (SmootherPrototype::IsSetup() == true)
      this->GetOStream(Warnings0) << "MueLu::UzawaSmoother::Setup(): Setup() has already been called";

    // extract blocked operator A from current level
    A_ = Factory::Get<RCP<Matrix> > (currentLevel, "A");

    RCP<BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A_);
    TEUCHOS_TEST_FOR_EXCEPTION(bA == Teuchos::null, Exceptions::BadCast, "MueLu::UzawaSmoother::Setup: input matrix A is not of type BlockedCrsMatrix! error.");

    // store map extractors
    rangeMapExtractor_ = bA->getRangeMapExtractor();
    domainMapExtractor_ = bA->getDomainMapExtractor();

    // Store the blocks in local member variables
    Teuchos::RCP<CrsMatrix> A00 = bA->getMatrix(0, 0);
    Teuchos::RCP<CrsMatrix> A01 = bA->getMatrix(0, 1);
    Teuchos::RCP<CrsMatrix> A10 = bA->getMatrix(1, 0);
    Teuchos::RCP<CrsMatrix> A11 = bA->getMatrix(1, 1);

    Teuchos::RCP<CrsMatrixWrap> Op00 = Teuchos::rcp(new CrsMatrixWrap(A00));
    Teuchos::RCP<CrsMatrixWrap> Op01 = Teuchos::rcp(new CrsMatrixWrap(A01));
    Teuchos::RCP<CrsMatrixWrap> Op10 = Teuchos::rcp(new CrsMatrixWrap(A10));
    Teuchos::RCP<CrsMatrixWrap> Op11 = Teuchos::rcp(new CrsMatrixWrap(A11));

    F_ = Teuchos::rcp_dynamic_cast<Matrix>(Op00);
    G_ = Teuchos::rcp_dynamic_cast<Matrix>(Op01);
    D_ = Teuchos::rcp_dynamic_cast<Matrix>(Op10);
    Z_ = Teuchos::rcp_dynamic_cast<Matrix>(Op11);

    // TODO move this to BlockedCrsMatrix->getMatrix routine...
    F_->CreateView("stridedMaps", bA->getRangeMap(0), bA->getDomainMap(0));
    G_->CreateView("stridedMaps", bA->getRangeMap(0), bA->getDomainMap(1));
    D_->CreateView("stridedMaps", bA->getRangeMap(1), bA->getDomainMap(0));
    Z_->CreateView("stridedMaps", bA->getRangeMap(1), bA->getDomainMap(1));

    // Set the Smoother
    // carefully switch to the SubFactoryManagers (defined by the users)
    {
      RCP<const FactoryManagerBase> velpredictFactManager = FactManager_.at(0);
      SetFactoryManager currentSFM  (rcpFromRef(currentLevel), velpredictFactManager);
      velPredictSmoo_ = currentLevel.Get< RCP<SmootherBase> >("PreSmoother",velpredictFactManager->GetFactory("Smoother").get());
    }
    {
      RCP<const FactoryManagerBase> schurFactManager = FactManager_.at(1);
      SetFactoryManager currentSFM  (rcpFromRef(currentLevel), schurFactManager);
      schurCompSmoo_ = currentLevel.Get< RCP<SmootherBase> >("PreSmoother", schurFactManager->GetFactory("Smoother").get());
    }

    SmootherPrototype::IsSetup(true);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void UzawaSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::UzawaSmoother::Apply(): Setup() has not been called");

    Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));

    SC zero = Teuchos::ScalarTraits<SC>::zero(), one = Teuchos::ScalarTraits<SC>::one();

    // extract parameters from internal parameter list
    const ParameterList & pL = Factory::GetParameterList();
    LocalOrdinal nSweeps = pL.get<LocalOrdinal>("Sweeps");
    Scalar omega = pL.get<Scalar>("Damping factor");

    // wrap current solution vector in RCP
    RCP<MultiVector> rcpX = Teuchos::rcpFromRef(X);

    // create residual vector
    // contains current residual of current solution X with rhs B
    RCP<MultiVector> residual = MultiVectorFactory::Build(B.getMap(), B.getNumVectors());

    // incrementally improve solution vector X
    for (LocalOrdinal run = 0; run < nSweeps; ++run) {
      // 1) calculate current residual
      residual->update(one,B,zero); // residual = B
      A_->apply(*rcpX, *residual, Teuchos::NO_TRANS, -one, one);

      // split residual vector
      Teuchos::RCP<MultiVector> r1 = rangeMapExtractor_->ExtractVector(residual, 0);
      Teuchos::RCP<MultiVector> r2 = rangeMapExtractor_->ExtractVector(residual, 1);

      // 2) solve F * \Delta \tilde{x}_1 = r_1
      //    start with zero guess \Delta \tilde{x}_1
      RCP<MultiVector> xtilde1 = MultiVectorFactory::Build(F_->getRowMap(),1);
      xtilde1->putScalar(zero);
      velPredictSmoo_->Apply(*xtilde1,*r1);

      // 3) calculate rhs for SchurComp equation
      //    r_2 - D \Delta \tilde{x}_1
      RCP<MultiVector> schurCompRHS = MultiVectorFactory::Build(Z_->getRowMap(),1);
      D_->apply(*xtilde1,*schurCompRHS);
      schurCompRHS->update(one,*r2,-one);

      // 4) solve SchurComp equation
      //    start with zero guess \Delta \tilde{x}_2
      RCP<MultiVector> xtilde2 = MultiVectorFactory::Build(Z_->getRowMap(),1);
      xtilde2->putScalar(zero);
      schurCompSmoo_->Apply(*xtilde2,*schurCompRHS);

      // 5) extract parts of solution vector X
      Teuchos::RCP<MultiVector> x1 = domainMapExtractor_->ExtractVector(rcpX, 0);
      Teuchos::RCP<MultiVector> x2 = domainMapExtractor_->ExtractVector(rcpX, 1);

      // 6) update solution vector with increments xhat1 and xhat2
      //    rescale increment for x2 with omega_
      x1->update(omega,*xtilde1,one); // x1 = x1_old + omega xtilde1
      x2->update(omega,*xtilde2,one); // x2 = x2_old + omega xtilde2

      // write back solution in global vector X
      domainMapExtractor_->InsertVector(x1, 0, rcpX);
      domainMapExtractor_->InsertVector(x2, 1, rcpX);

    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  UzawaSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Copy() const {
    return rcp( new UzawaSmoother(*this) );
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string UzawaSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::description() const {
    std::ostringstream out;
    out << SmootherPrototype::description();
    out << "{type = " << type_ << "}";
    return out.str();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void UzawaSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0) {
      out0 << "Prec. type: " << type_ << /*" Sweeps: " << nSweeps_ << " damping: " << omega_ <<*/ std::endl;
    }

    if (verbLevel & Debug) {
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl;
    }
  }

} // namespace MueLu


#endif /* MUELU_UZAWASMOOTHER_DEF_HPP_ */
