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
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_ReorderedBlockedCrsMatrix.hpp>

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
    //validParamList->set<ParameterList>          ("block1",  ParameterList(), "Sublist for parameters for SchurComplement block. At least has to contain some information about a smoother \"Smoother\" for variable \"A\" which is generated by a SchurComplementFactory.");
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
    A00_ = bA->getMatrix(0,0);
    A01_ = bA->getMatrix(0,1);
    A10_ = bA->getMatrix(1,0);
    A11_ = bA->getMatrix(1,1);

    const ParameterList& pL = Factory::GetParameterList();
    SC omega = pL.get<SC>("Damping factor");

#if 0 // old code
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
      Ddata[row] = one / (diag[row]*omega);*/
#else
    // Create the inverse of the diagonal of F
    // TODO add safety check for zeros on diagonal of F!
    RCP<Vector> diagFVector = VectorFactory::Build(A00_->getRowMap());
    if (pL.get<bool>("lumping") == false) {
      A00_->getLocalDiagCopy(*diagFVector);       // extract diagonal of F
    } else {
      diagFVector = Utilities::GetLumpedMatrixDiagonal(A00_);
    }
    diagFVector->scale(omega);
    D_ = Utilities::GetInverse(diagFVector);

    // check whether D_ is a blocked vector with only 1 block
    RCP<BlockedVector> bD = Teuchos::rcp_dynamic_cast<BlockedVector>(D_);
    if(bD.is_null() == false && bD->getBlockedMap()->getNumMaps() == 1) {
      RCP<Vector> nestedVec = bD->getMultiVector(0,bD->getBlockedMap()->getThyraMode())->getVectorNonConst(0);
      D_.swap(nestedVec);
    }
#endif

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

#if 0 //def HAVE_MUELU_DEBUG
    // TODO simplify this debug check
    RCP<MultiVector> rcpDebugX = Teuchos::rcpFromRef(X);
    RCP<const MultiVector> rcpDebugB = Teuchos::rcpFromRef(B);
    RCP<BlockedMultiVector> rcpBDebugX = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(rcpDebugX);
    RCP<const BlockedMultiVector> rcpBDebugB = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(rcpDebugB);
    RCP<BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A_);
    if(rcpBDebugB.is_null() == false) {
      //TEUCHOS_TEST_FOR_EXCEPTION(A_->getRangeMap()->isSameAs(*(B.getMap())) == false, Exceptions::RuntimeError, "MueLu::BlockedGaussSeidelSmoother::Apply(): The map of RHS vector B is not the same as range map of the blocked operator A. Please check the map of B and A.");
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(bA->getFullRangeMap()->isSameAs(*(B.getMap())) == false, Exceptions::RuntimeError, "MueLu::BlockedGaussSeidelSmoother::Apply(): The map of RHS vector B is not the same as range map of the blocked operator A. Please check the map of B and A.");
    }
    if(rcpBDebugX.is_null() == false) {
      //TEUCHOS_TEST_FOR_EXCEPTION(A_->getDomainMap()->isSameAs(*(X.getMap())) == false, Exceptions::RuntimeError, "MueLu::BlockedGaussSeidelSmoother::Apply(): The map of the solution vector X is not the same as domain map of the blocked operator A. Please check the map of X and A.");
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(bA->getFullDomainMap()->isSameAs(*(X.getMap())) == false, Exceptions::RuntimeError, "MueLu::BlockedGaussSeidelSmoother::Apply(): The map of the solution vector X is not the same as domain map of the blocked operator A. Please check the map of X and A.");
    }
#endif


    RCP<MultiVector> rcpX       = rcpFromRef(X);
    RCP<const MultiVector> rcpB = rcpFromRef(B);

    // make sure that both rcpX and rcpB are BlockedMultiVector objects
    bool bCopyResultX = false;
    RCP<BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A_);
    MUELU_TEST_FOR_EXCEPTION(bA.is_null() == true, Exceptions::RuntimeError, "MueLu::BlockedGaussSeidelSmoother::Apply(): A_ must be a BlockedCrsMatrix");
    RCP<BlockedMultiVector> bX = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(rcpX);
    RCP<const BlockedMultiVector> bB = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(rcpB);

    if(bX.is_null() == true) {
      RCP<MultiVector> test = Teuchos::rcp(new BlockedMultiVector(bA->getBlockedDomainMap(),rcpX));
      rcpX.swap(test);
      bCopyResultX = true;
    }

    if(bB.is_null() == true) {
      RCP<const MultiVector> test = Teuchos::rcp(new BlockedMultiVector(bA->getBlockedRangeMap(),rcpB));
      rcpB.swap(test);
    }

    // we now can guarantee that X and B are blocked multi vectors
    bX = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(rcpX);
    bB = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(rcpB);

    // check the type of operator
    RCP<Xpetra::ReorderedBlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > rbA = Teuchos::rcp_dynamic_cast<Xpetra::ReorderedBlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(bA);
    if(rbA.is_null() == false) {
      // A is a ReorderedBlockedCrsMatrix
      Teuchos::RCP<const Xpetra::BlockReorderManager > brm = rbA->getBlockReorderManager();

      // check type of X vector
      if(bX->getBlockedMap()->getNumMaps() != bA->getDomainMapExtractor()->NumMaps()) {
        // X is a blocked multi vector but incompatible to the reordered blocked operator A
        Teuchos::RCP<MultiVector> test =
            buildReorderedBlockedMultiVector(brm, bX);
        rcpX.swap(test);
      }
      if(bB->getBlockedMap()->getNumMaps() != bA->getRangeMapExtractor()->NumMaps()) {
        // B is a blocked multi vector but incompatible to the reordered blocked operator A
        Teuchos::RCP<const MultiVector> test =
            buildReorderedBlockedMultiVector(brm, bB);
        rcpB.swap(test);
      }
    }

    // use the GIDs of the sub blocks
    // This is valid as the subblocks actually represent the GIDs (either Thyra, Xpetra or pseudo Xpetra)

    bool bRangeThyraMode  = rangeMapExtractor_->getThyraMode();
    bool bDomainThyraMode = domainMapExtractor_->getThyraMode();

    RCP<MultiVector> deltaX     = MultiVectorFactory::Build(rcpX->getMap(), rcpX->getNumVectors());
    RCP<BlockedMultiVector> bdeltaX = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(deltaX);
    RCP<MultiVector> deltaX0 = bdeltaX->getMultiVector(0,bDomainThyraMode);
    RCP<MultiVector> deltaX1 = bdeltaX->getMultiVector(1,bDomainThyraMode);

    RCP<MultiVector> Rtmp    = rangeMapExtractor_->getVector(1, rcpB->getNumVectors(), bRangeThyraMode);


    typedef Teuchos::ScalarTraits<SC> STS;
    SC one = STS::one(), zero = STS::zero();

    // extract parameters from internal parameter list
    const ParameterList& pL = Factory::GetParameterList();
    LO nSweeps = pL.get<LO>("Sweeps");

    RCP<MultiVector> R;// = MultiVectorFactory::Build(rcpB->getMap(), rcpB->getNumVectors());;
    if (InitialGuessIsZero)  {
      R = MultiVectorFactory::Build(rcpB->getMap(), rcpB->getNumVectors());
      R->update(one, *rcpB, zero);
    } else {
      R = Utilities::Residual(*A_, *rcpX, *rcpB);
    }

    // extract diagonal of Schur complement operator
    RCP<Vector> diagSVector = VectorFactory::Build(S_->getRowMap());
    S_->getLocalDiagCopy(*diagSVector);

    for (LO run = 0; run < nSweeps; ++run) {
      // Extract corresponding subvectors from X and R
      // Automatically detect whether we use Thyra or Xpetra GIDs
      // The GIDs should be always compatible with the GIDs of A00, A01, etc...
      RCP<MultiVector> R0 = rangeMapExtractor_ ->ExtractVector(R, 0, bRangeThyraMode);
      RCP<MultiVector> R1 = rangeMapExtractor_ ->ExtractVector(R, 1, bRangeThyraMode);

      // Calculate Rtmp = R1 - D * deltaX0 (equation 8.14)
      deltaX0->putScalar(zero);
      deltaX0->elementWiseMultiply(one, *D_, *R0, zero);    // deltaX0 = D * R0 (equation 8.13)
      A10_->apply(*deltaX0, *Rtmp);                         // Rtmp    = A10*D*deltaX0 (intermediate step)
      Rtmp->update(one, *R1, -one);                         // Rtmp    = R1 - A10*D*deltaX0

      if (!pL.get<bool>("q2q1 mode")) {
        deltaX1->putScalar(zero);
      } else {
        // special code for q2q1
        if(Teuchos::rcp_dynamic_cast<BlockedVector>(diagSVector) == Teuchos::null) {
          ArrayRCP<SC> Sdiag = diagSVector->getDataNonConst(0);
          ArrayRCP<SC> deltaX1data = deltaX1->getDataNonConst(0);
          ArrayRCP<SC> Rtmpdata    = Rtmp->getDataNonConst(0);
          for (GO row = 0; row < deltaX1data.size(); row++)
            deltaX1data[row] = Teuchos::as<SC>(1.1)*Rtmpdata[row] / Sdiag[row];
        } else {
          TEUCHOS_TEST_FOR_EXCEPTION(true,MueLu::Exceptions::RuntimeError,"MueLu::BraessSarazinSmoother: q2q1 mode only supported for non-blocked operators.")
        }
      }

      smoo_->Apply(*deltaX1,*Rtmp);

      // Compute deltaX0
      deltaX0->putScalar(zero);                             // just for safety
      A01_->apply(*deltaX1, *deltaX0);                      // deltaX0 = A01*deltaX1
      deltaX0->update(one, *R0, -one);                      // deltaX0 = R0 - A01*deltaX1
      R0.swap(deltaX0);
      deltaX0->elementWiseMultiply(one, *D_, *R0, zero);    // deltaX0 = D*(R0 - A01*deltaX1)

      RCP<MultiVector> X0 = domainMapExtractor_->ExtractVector(rcpX, 0, bDomainThyraMode);
      RCP<MultiVector> X1 = domainMapExtractor_->ExtractVector(rcpX, 1, bDomainThyraMode);

      // Update solution
      X0->update(one, *deltaX0, one);
      X1->update(one, *deltaX1, one);

      domainMapExtractor_->InsertVector(X0, 0, rcpX, bDomainThyraMode);
      domainMapExtractor_->InsertVector(X1, 1, rcpX, bDomainThyraMode);

      if (run < nSweeps-1) {
        R = Utilities::Residual(*A_, *rcpX, *rcpB);
      }

    }

    if (bCopyResultX == true) {
      RCP<MultiVector> Xmerged = bX->Merge();
      X.update(one, *Xmerged, zero);
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

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNodeSmootherComplexity() const {
    // FIXME: This is a placeholder
    return Teuchos::OrdinalTraits<size_t>::invalid();
  }
  

} // namespace MueLu

#endif /* MUELU_BRAESSSARAZINSMOOTHER_DEF_HPP_ */
