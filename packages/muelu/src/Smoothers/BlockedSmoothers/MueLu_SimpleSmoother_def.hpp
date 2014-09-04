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
 * MueLu_SimpleSmoother_def.hpp
 *
 *  Created on: 19.03.2013
 *      Author: wiesner
 */

#ifndef MUELU_SIMPLESMOOTHER_DEF_HPP_
#define MUELU_SIMPLESMOOTHER_DEF_HPP_

#include "Teuchos_ArrayViewDecl.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "MueLu_ConfigDefs.hpp"

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_SimpleSmoother_decl.hpp"
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
  SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SimpleSmoother()
    : type_("SIMPLE"), A_(Teuchos::null)
  {
    //Factory::SetParameter("Sweeps", Teuchos::ParameterEntry(sweeps));
    //Factory::SetParameter("Damping factor",Teuchos::ParameterEntry(omega));
    //Factory::SetParameter("UseSIMPLEC", Teuchos::ParameterEntry(SIMPLEC));

#if 0
    // when declaring default factories without overwriting them leads to a multipleCallCheck exception
    // TODO: debug into this
    // workaround: always define your factory managers outside either using the C++ API or the XML files
    RCP<SchurComplementFactory> SchurFact = Teuchos::rcp(new SchurComplementFactory());
    SchurFact->SetParameter("omega",Teuchos::ParameterEntry(omega));
    SchurFact->SetParameter("lumping",Teuchos::ParameterEntry(SIMPLEC));
    SchurFact->SetFactory("A", this->GetFactory("A"));

    // define smoother/solver for SchurComplement equation
    Teuchos::ParameterList SCparams;
    std::string SCtype;
    RCP<SmootherPrototype> smoProtoSC     = rcp( new DirectSolver(SCtype,SCparams) );
    smoProtoSC->SetFactory("A", SchurFact);

    RCP<SmootherFactory> SmooSCFact = rcp( new SmootherFactory(smoProtoSC) );

    RCP<FactoryManager> schurFactManager = rcp(new FactoryManager());
    schurFactManager->SetFactory("A", SchurFact);
    schurFactManager->SetFactory("Smoother", SmooSCFact);
    schurFactManager->SetIgnoreUserData(true);

    // define smoother/solver for velocity prediction
    RCP<SubBlockAFactory> A00Fact = Teuchos::rcp(new SubBlockAFactory(/*this->GetFactory("A"), 0, 0*/));
    A00Fact->SetFactory("A",this->GetFactory("A"));
    A00Fact->SetParameter("block row",ParameterEntry(0));
    A00Fact->SetParameter("block col",ParameterEntry(0));
    Teuchos::ParameterList velpredictParams;
    std::string velpredictType;
    RCP<SmootherPrototype> smoProtoPredict     = rcp( new DirectSolver(velpredictType,velpredictParams) );
    smoProtoPredict->SetFactory("A", A00Fact);
    RCP<SmootherFactory> SmooPredictFact = rcp( new SmootherFactory(smoProtoPredict) );

    RCP<FactoryManager> velpredictFactManager = rcp(new FactoryManager());
    velpredictFactManager->SetFactory("A", A00Fact);
    velpredictFactManager->SetFactory("Smoother", SmooPredictFact);
    velpredictFactManager->SetIgnoreUserData(true);

    AddFactoryManager(velpredictFactManager, 0);
    AddFactoryManager(schurFactManager, 1);
#endif
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~SimpleSmoother() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const ParameterList> SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< RCP<const FactoryBase> >("A",                  Teuchos::null, "Generating factory of the matrix A");
    validParamList->set< Scalar >                ("Damping factor",     1.0, "Damping/Scaling factor in SIMPLE");
    validParamList->set< LocalOrdinal >          ("Sweeps",             1, "Number of SIMPLE sweeps (default = 1)");
    validParamList->set< bool >                  ("UseSIMPLEC",         false, "Use SIMPLEC instead of SIMPLE (default = false)");

    return validParamList;
  }

  //! Add a factory manager at a specific position
  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AddFactoryManager(RCP<const FactoryManagerBase> FactManager, int pos) {
    TEUCHOS_TEST_FOR_EXCEPTION(pos < 0, Exceptions::RuntimeError, "MueLu::SimpleSmoother::AddFactoryManager: parameter \'pos\' must not be negative! error.");

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
  void SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetVelocityPredictionFactoryManager(RCP<FactoryManager> FactManager) {
    AddFactoryManager(FactManager, 0); // overwrite factory manager for predicting the primary variable
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetSchurCompFactoryManager(RCP<FactoryManager> FactManager) {
    AddFactoryManager(FactManager, 1); // overwrite factory manager for SchurComplement
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    currentLevel.DeclareInput("A",this->GetFactory("A").get());

    TEUCHOS_TEST_FOR_EXCEPTION(FactManager_.size() != 2, Exceptions::RuntimeError,"MueLu::SimpleSmoother::DeclareInput: You have to declare two FactoryManagers with a \"Smoother\" object: One for predicting the primary variable and one for the SchurComplement system. The smoother for the SchurComplement system needs a SchurComplementFactory as input for variable \"A\". make sure that you use the same proper damping factors for omega both in the SchurComplementFactory and in the SIMPLE smoother!");

    // loop over all factory managers for the subblocks of blocked operator A
    std::vector<Teuchos::RCP<const FactoryManagerBase> >::const_iterator it;
    for(it = FactManager_.begin(); it!=FactManager_.end(); ++it) {
      SetFactoryManager currentSFM  (rcpFromRef(currentLevel),   *it);

      // request "Smoother" for current subblock row.
      currentLevel.DeclareInput("PreSmoother",(*it)->GetFactory("Smoother").get());
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Setup(Level &currentLevel) {
    //*********************************************
    // Setup routine can be summarized in 4 steps:
    // - Set the map extractors
    // - Set the blocks
    // - Create and set the inverse of the diagonal of F
    // - Set the smoother for the Schur Complement

    FactoryMonitor m(*this, "Setup blocked SIMPLE Smoother", currentLevel);

    if (SmootherPrototype::IsSetup() == true)
      this->GetOStream(Warnings0) << "MueLu::SimpleSmoother::Setup(): Setup() has already been called";

    // extract blocked operator A from current level
    A_ = Factory::Get<RCP<Matrix> > (currentLevel, "A");

    RCP<BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A_);
    TEUCHOS_TEST_FOR_EXCEPTION(bA == Teuchos::null, Exceptions::BadCast, "MueLu::SimpleSmoother::Setup: input matrix A is not of type BlockedCrsMatrix! error.");

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

    const ParameterList & pL = Factory::GetParameterList();
    bool bSIMPLEC = pL.get<bool>("UseSIMPLEC");

    // Create the inverse of the diagonal of F
    RCP<Vector> diagFVector = VectorFactory::Build(F_->getRowMap());
    if(!bSIMPLEC) {
      F_->getLocalDiagCopy(*diagFVector);       // extract diagonal of F
      diagFVector->reciprocal(*diagFVector);    // build reciprocal
    } else {
      const RCP<const Map> rowmap = F_->getRowMap();
      size_t locSize = rowmap->getNodeNumElements();
      Teuchos::ArrayRCP<SC> diag = diagFVector->getDataNonConst(0);
      Teuchos::ArrayView<const LO> cols;
      Teuchos::ArrayView<const SC> vals;
      for (size_t i=0; i<locSize; ++i) { // loop over rows
        F_->getLocalRowView(i,cols,vals);
        Scalar absRowSum = Teuchos::ScalarTraits<Scalar>::zero();
        for (LO j=0; j<cols.size(); ++j) { // loop over cols
          absRowSum += Teuchos::ScalarTraits<Scalar>::magnitude(vals[j]);
        }
        diag[i] = absRowSum;
      }
      diagFVector->reciprocal(*diagFVector);    // build reciprocal
    }
    diagFinv_ = diagFVector;

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
  void SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::SimpleSmoother::Apply(): Setup() has not been called");

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

      // 5) scale xtilde2 with omega
      //    store this in xhat2
      RCP<MultiVector> xhat2 = MultiVectorFactory::Build(Z_->getRowMap(),1);
      xhat2->update(omega,*xtilde2,zero);

      // 6) calculate xhat1
      RCP<MultiVector> xhat1      = MultiVectorFactory::Build(F_->getRowMap(),1);
      RCP<MultiVector> xhat1_temp = MultiVectorFactory::Build(F_->getRowMap(),1);
      G_->apply(*xhat2,*xhat1_temp); // store result temporarely in xtilde1_temp
      xhat1->elementWiseMultiply(one/*/omega*/,*diagFinv_,*xhat1_temp,zero);
      xhat1->update(one,*xtilde1,-one);

      // 7) extract parts of solution vector X
      Teuchos::RCP<MultiVector> x1 = domainMapExtractor_->ExtractVector(rcpX, 0);
      Teuchos::RCP<MultiVector> x2 = domainMapExtractor_->ExtractVector(rcpX, 1);

      // 8) update solution vector with increments xhat1 and xhat2
      //    rescale increment for x2 with omega_
      x1->update(one,*xhat1,one);    // x1 = x1_old + xhat1
      x2->update(/*omega*/ one,*xhat2,one); // x2 = x2_old + omega xhat2

      // write back solution in global vector X
      domainMapExtractor_->InsertVector(x1, 0, rcpX);
      domainMapExtractor_->InsertVector(x2, 1, rcpX);

    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Copy () const {
    return rcp( new SimpleSmoother(*this) );
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::description() const {
    std::ostringstream out;
    out << SmootherPrototype::description();
    out << "{type = " << type_ << "}";
    return out.str();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0) {
      out0 << "Prec. type: " << type_ << /*" Sweeps: " << nSweeps_ << " damping: " << omega_ <<*/ std::endl;
    }

    if (verbLevel & Debug) {
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl;
    }
  }

} // namespace MueLu


#endif /* MUELU_SIMPLESMOOTHER_DEF_HPP_ */
