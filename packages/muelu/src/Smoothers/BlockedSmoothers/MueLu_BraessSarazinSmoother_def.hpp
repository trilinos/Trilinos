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
 * MueLu_BraessSarazinSmoother_def.hpp
 *
 *  Created on: Apr 16, 2012
 *      Author: wiesner
 */

#ifndef MUELU_BRAESSSARAZINSMOOTHER_DEF_HPP_
#define MUELU_BRAESSSARAZINSMOOTHER_DEF_HPP_

#include "Teuchos_ArrayViewDecl.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "MueLu_ConfigDefs.hpp"

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_BraessSarazinSmoother_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_HierarchyUtils.hpp"
#include "MueLu_SmootherBase.hpp"

// include files for default FactoryManager
#include "MueLu_SchurComplementFactory.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_FactoryManager.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BraessSarazinSmoother()
    : type_("Braess Sarazin"), A_(null)
  {
    //Factory::SetParameter("Sweeps", ParameterEntry(sweeps));
    //Factory::SetParameter("Damping factor",ParameterEntry(omega));

#if 0
    // when declaring default factories without overwriting them leads to a multipleCallCheck exception
    // TODO: debug into this
    // workaround: always define your factory managers outside either using the C++ API or the XML files
    RCP<SchurComplementFactory> SchurFact = rcp(new SchurComplementFactory());
    SchurFact->SetParameter("omega",ParameterEntry(omega));
    SchurFact->SetFactory("A", this->GetFactory("A"));

    // define smoother/solver for BraessSarazin
    ParameterList SCparams;
    std::string SCtype;
    RCP<SmootherPrototype> smoProtoSC     = rcp( new DirectSolver(SCtype,SCparams) );
    smoProtoSC->SetFactory("A", SchurFact);

    RCP<SmootherFactory> SmooSCFact = rcp( new SmootherFactory(smoProtoSC) );

    RCP<FactoryManager> FactManager = rcp(new FactoryManager());
    FactManager->SetFactory("A", SchurFact);
    FactManager->SetFactory("Smoother", SmooSCFact);
    FactManager->SetIgnoreUserData(true);

    AddFactoryManager(FactManager,0);
#endif
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~BraessSarazinSmoother() {}

  //! Add a factory manager at a specific position
  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddFactoryManager(RCP<const FactoryManagerBase> FactManager, int pos) {
    TEUCHOS_TEST_FOR_EXCEPTION(pos != 0, Exceptions::RuntimeError, "MueLu::BraessSarazinSmoother::AddFactoryManager: parameter \'pos\' must be zero! error.");
    FactManager_ = FactManager;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    SC one = Teuchos::ScalarTraits<SC>::one();

    validParamList->set<RCP<const FactoryBase> >("A",                  null, "Generating factory of the matrix A");

    validParamList->set<bool>                   ("lumping",           false, "Use lumping to construct diag(A(0,0)), i.e. use row sum of the abs values on the diagonal "
                                                                             "as approximation of A00 (and A00^{-1})");
    validParamList->set<SC>                     ("Damping factor",      one, "Damping/Scaling factor in BraessSarazin (usually has to be chosen > 1, default = 1.0)");
    validParamList->set<LO>                     ("Sweeps",                1, "Number of BraessSarazin sweeps (default = 1)");
    validParamList->set<ParameterList>          ("block1",  ParameterList(), "Sublist for parameters for SchurComplement block. At least has to contain some information about a smoother \"Smoother\" for variable \"A\" which is generated by a SchurComplementFactory.");
    validParamList->set<bool>                   ("q2q1 mode",         false, "Run in the mode matching Q2Q1 matlab code");

    return validParamList;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
    this->Input(currentLevel, "A");

    TEUCHOS_TEST_FOR_EXCEPTION(FactManager_.is_null(), Exceptions::RuntimeError,
                               "MueLu::BraessSarazinSmoother::DeclareInput: FactManager_ must not be null! "
                               "Introduce a FactoryManager for the SchurComplement equation.");

    // carefully call DeclareInput after switching to internal FactoryManager
    {
      SetFactoryManager currentSFM(rcpFromRef(currentLevel), FactManager_);

      // request "Smoother" for current subblock row.
      currentLevel.DeclareInput("PreSmoother", FactManager_->GetFactory("Smoother").get());

      // request Schur matrix just in case
      currentLevel.DeclareInput("A", FactManager_->GetFactory("A").get());
    }
  }

  // Setup routine can be summarized in 4 steps:
  // - set the map extractors
  // - set the blocks
  // - create and set the inverse of the diagonal of F
  // - set the smoother for the Schur Complement
  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level& currentLevel) {
    FactoryMonitor m(*this, "Setup Smoother", currentLevel);

    if (SmootherPrototype::IsSetup() == true)
      this->GetOStream(Warnings0) << "MueLu::BreaessSarazinSmoother::Setup(): Setup() has already been called";

    // Extract blocked operator A from current level
    A_ = Factory::Get<RCP<Matrix> > (currentLevel, "A");
    RCP<BlockedCrsMatrix> bA = rcp_dynamic_cast<BlockedCrsMatrix>(A_);
    TEUCHOS_TEST_FOR_EXCEPTION(bA.is_null(), Exceptions::BadCast,
                               "MueLu::BraessSarazinSmoother::Setup: input matrix A is not of type BlockedCrsMatrix! error.");

    // Store map extractors
    rangeMapExtractor_  = bA->getRangeMapExtractor();
    domainMapExtractor_ = bA->getDomainMapExtractor();

    // Store the blocks in local member variables
    A00_ = MueLu::Utilities<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Crs2Op(bA->getMatrix(0,0));
    A01_ = MueLu::Utilities<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Crs2Op(bA->getMatrix(0,1));
    A10_ = MueLu::Utilities<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Crs2Op(bA->getMatrix(1,0));
    A11_ = MueLu::Utilities<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Crs2Op(bA->getMatrix(1,1));

    // TODO move this to BlockedCrsMatrix->getMatrix routine...
    A00_->CreateView("stridedMaps", bA->getRangeMap(0), bA->getDomainMap(0));
    A01_->CreateView("stridedMaps", bA->getRangeMap(0), bA->getDomainMap(1));
    A10_->CreateView("stridedMaps", bA->getRangeMap(1), bA->getDomainMap(0));
    if (!A11_.is_null())
      A11_->CreateView("stridedMaps", bA->getRangeMap(1), bA->getDomainMap(1));

    const ParameterList& pL = Factory::GetParameterList();
    SC omega = pL.get<SC>("Damping factor");

    // Create the inverse of the diagonal of F
    D_ = VectorFactory::Build(A00_->getRowMap());

    ArrayRCP<SC> diag;
    if (pL.get<bool>("lumping") == false)
      diag = Utilities::GetMatrixDiagonal      (*A00_);
    else
      diag = Utilities::GetLumpedMatrixDiagonal(*A00_);

    SC one = Teuchos::ScalarTraits<SC>::one();

    ArrayRCP<SC> Ddata = D_->getDataNonConst(0);
    for (GO row = 0; row < Ddata.size(); row++)
      Ddata[row] = one / (diag[row]*omega);

    // Set the Smoother
    // carefully switch to the SubFactoryManagers (defined by the users)
    {
      SetFactoryManager currentSFM(rcpFromRef(currentLevel), FactManager_);
      smoo_ = currentLevel.Get<RCP<SmootherBase> >("PreSmoother", FactManager_->GetFactory("Smoother").get());
      S_    = currentLevel.Get<RCP<Matrix> >      ("A",           FactManager_->GetFactory("A").get());
    }

    SmootherPrototype::IsSetup(true);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError,
                               "MueLu::BraessSarazinSmoother::Apply(): Setup() has not been called");

    RCP<MultiVector> rcpX    = rcpFromRef(X);
    RCP<MultiVector> deltaX0 = MultiVectorFactory::Build(A00_->getRowMap(), 1);
    RCP<MultiVector> deltaX1 = MultiVectorFactory::Build(A10_->getRowMap(), 1);
    RCP<MultiVector> Rtmp    = MultiVectorFactory::Build(A10_->getRowMap(), 1);

    typedef Teuchos::ScalarTraits<SC> STS;
    SC one = STS::one(), zero = STS::zero();

    // extract parameters from internal parameter list
    const ParameterList& pL = Factory::GetParameterList();
    LO nSweeps = pL.get<LO>("Sweeps");

    RCP<MultiVector> R;
    if (InitialGuessIsZero)  {
      R = MultiVectorFactory::Build(B.getMap(), B.getNumVectors());
      R->update(one, B, zero);
    } else {
      R = Utilities::Residual(*A_, X, B);
    }

    ArrayRCP<SC> Sdiag = Utilities::GetMatrixDiagonal(*S_);

    for (LO run = 0; run < nSweeps; ++run) {
      // Extract corresponding subvectors from X and R
      RCP<MultiVector> R0 = rangeMapExtractor_ ->ExtractVector(R, 0);
      RCP<MultiVector> R1 = rangeMapExtractor_ ->ExtractVector(R, 1);

      RCP<MultiVector> X0 = domainMapExtractor_->ExtractVector(rcpX, 0);
      RCP<MultiVector> X1 = domainMapExtractor_->ExtractVector(rcpX, 1);

      // Calculate Rtmp = R1 - D * deltaX0 (equation 8.14)
      deltaX0->putScalar(zero);
      deltaX0->elementWiseMultiply(one, *D_, *R0, zero);    // deltaX0 = D * R0 (equation 8.13)
      A10_->apply(*deltaX0, *Rtmp);                         // Rtmp    = A10*D*deltaX0 (intermediate step)
      Rtmp->update(one, *R1, -one);                         // Rtmp    = R1 - A10*D*deltaX0

      if (!pL.get<bool>("q2q1 mode")) {
        deltaX1->putScalar(zero);
      } else {
        ArrayRCP<SC> deltaX1data = deltaX1->getDataNonConst(0);
        ArrayRCP<SC> Rtmpdata    = Rtmp->getDataNonConst(0);
        for (GO row = 0; row < deltaX1data.size(); row++)
          deltaX1data[row] = 1.1*Rtmpdata[row] / Sdiag[row];
      }

      // Compute deltaX1 (pressure correction)
      // We use user provided preconditioner
      smoo_->Apply(*deltaX1, *Rtmp);

      // Compute deltaX0
      deltaX0->putScalar(zero);                             // just for safety
      A01_->apply(*deltaX1, *deltaX0);                      // deltaX0 = A01*deltaX1
      deltaX0->update(one, *R0, -one);                      // deltaX0 = R0 - A01*deltaX1
      R0.swap(deltaX0);
      deltaX0->elementWiseMultiply(one, *D_, *R0, zero);    // deltaX0 = D*(R0 - A01*deltaX1)

      // Update solution
      X0->update(one, *deltaX0, one);
      X1->update(one, *deltaX1, one);

      domainMapExtractor_->InsertVector(X0, 0, rcpX);
      domainMapExtractor_->InsertVector(X1, 1, rcpX);

      if (run < nSweeps-1)
        R = Utilities::Residual(*A_, X, B);
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Copy () const {
    return rcp (new BraessSarazinSmoother (*this));
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  description () const {
    std::ostringstream out;
    out << SmootherPrototype::description();
    out << "{type = " << type_ << "}";
    return out.str();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0) {
      out0 << "Prec. type: " << type_ << /*" Sweeps: " << nSweeps_ << " damping: " << omega_ <<*/ std::endl;
    }

    if (verbLevel & Debug) {
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl;
    }
  }

} // namespace MueLu

#endif /* MUELU_BRAESSSARAZINSMOOTHER_DEF_HPP_ */
