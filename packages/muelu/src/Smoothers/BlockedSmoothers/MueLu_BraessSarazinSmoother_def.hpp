// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_BRAESSSARAZINSMOOTHER_DEF_HPP_
#define MUELU_BRAESSSARAZINSMOOTHER_DEF_HPP_

#include <Teuchos_ArrayViewDecl.hpp>
#include <Teuchos_ScalarTraits.hpp>

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

#include "MueLu_FactoryManager.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddFactoryManager(RCP<const FactoryManagerBase> FactManager, int pos) {
  TEUCHOS_TEST_FOR_EXCEPTION(pos != 0, Exceptions::RuntimeError, "MueLu::BraessSarazinSmoother::AddFactoryManager: parameter \'pos\' must be zero! error.");
  FactManager_ = FactManager;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  SC one = Teuchos::ScalarTraits<SC>::one();

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory of the matrix A");

  validParamList->set<bool>("lumping", false,
                            "Use lumping to construct diag(A(0,0)), "
                            "i.e. use row sum of the abs values on the diagonal as approximation of A00 (and A00^{-1})");
  validParamList->set<SC>("Damping factor", one, "Damping/Scaling factor in BraessSarazin (usually has to be chosen > 1, default = 1.0)");
  validParamList->set<LO>("Sweeps", 1, "Number of BraessSarazin sweeps (default = 1)");
  validParamList->set<bool>("q2q1 mode", false, "Run in the mode matching Q2Q1 matlab code");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
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

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level& currentLevel) {
  FactoryMonitor m(*this, "Setup Smoother", currentLevel);

  if (SmootherPrototype::IsSetup() == true)
    this->GetOStream(Warnings0) << "MueLu::BreaessSarazinSmoother::Setup(): Setup() has already been called";

  // Extract blocked operator A from current level
  A_                       = Factory::Get<RCP<Matrix> >(currentLevel, "A");
  RCP<BlockedCrsMatrix> bA = rcp_dynamic_cast<BlockedCrsMatrix>(A_);
  TEUCHOS_TEST_FOR_EXCEPTION(bA.is_null(), Exceptions::BadCast,
                             "MueLu::BraessSarazinSmoother::Setup: input matrix A is not of type BlockedCrsMatrix! error.");

  // Store map extractors
  rangeMapExtractor_  = bA->getRangeMapExtractor();
  domainMapExtractor_ = bA->getDomainMapExtractor();

  // Store the blocks in local member variables
  A00_ = bA->getMatrix(0, 0);
  A01_ = bA->getMatrix(0, 1);
  A10_ = bA->getMatrix(1, 0);
  A11_ = bA->getMatrix(1, 1);

  const ParameterList& pL = Factory::GetParameterList();
  SC omega                = pL.get<SC>("Damping factor");

  // Create the inverse of the diagonal of F
  // TODO add safety check for zeros on diagonal of F!
  RCP<Vector> diagFVector = VectorFactory::Build(A00_->getRowMap());
  if (pL.get<bool>("lumping") == false) {
    A00_->getLocalDiagCopy(*diagFVector);  // extract diagonal of F
  } else {
    diagFVector = Utilities::GetLumpedMatrixDiagonal(*A00_);
  }
  diagFVector->scale(omega);
  D_ = Utilities::GetInverse(diagFVector);

  // check whether D_ is a blocked vector with only 1 block
  RCP<BlockedVector> bD = Teuchos::rcp_dynamic_cast<BlockedVector>(D_);
  if (bD.is_null() == false && bD->getBlockedMap()->getNumMaps() == 1) {
    RCP<Vector> nestedVec = bD->getMultiVector(0, bD->getBlockedMap()->getThyraMode())->getVectorNonConst(0);
    D_.swap(nestedVec);
  }

  // Set the Smoother
  // carefully switch to the SubFactoryManagers (defined by the users)
  {
    SetFactoryManager currentSFM(rcpFromRef(currentLevel), FactManager_);
    smoo_ = currentLevel.Get<RCP<SmootherBase> >("PreSmoother", FactManager_->GetFactory("Smoother").get());
    S_    = currentLevel.Get<RCP<Matrix> >("A", FactManager_->GetFactory("A").get());
  }

  SmootherPrototype::IsSetup(true);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const {
  TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError,
                             "MueLu::BraessSarazinSmoother::Apply(): Setup() has not been called");

  RCP<MultiVector> rcpX       = rcpFromRef(X);
  RCP<const MultiVector> rcpB = rcpFromRef(B);

  // make sure that both rcpX and rcpB are BlockedMultiVector objects
  bool bCopyResultX        = false;
  bool bReorderX           = false;
  bool bReorderB           = false;
  RCP<BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A_);
  MUELU_TEST_FOR_EXCEPTION(bA.is_null() == true, Exceptions::RuntimeError, "MueLu::BlockedGaussSeidelSmoother::Apply(): A_ must be a BlockedCrsMatrix");
  RCP<BlockedMultiVector> bX       = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(rcpX);
  RCP<const BlockedMultiVector> bB = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(rcpB);

  // check the type of operator
  RCP<Xpetra::ReorderedBlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rbA =
      Teuchos::rcp_dynamic_cast<Xpetra::ReorderedBlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(bA);

  if (rbA.is_null() == false) {
    // A is a ReorderedBlockedCrsMatrix and thus uses nested maps, retrieve BlockedCrsMatrix and use plain blocked
    // map for the construction of the blocked multivectors

    // check type of X vector
    if (bX.is_null() == true) {
      RCP<MultiVector> vectorWithBlockedMap = Teuchos::rcp(new BlockedMultiVector(rbA->getBlockedCrsMatrix()->getBlockedDomainMap(), rcpX));
      rcpX.swap(vectorWithBlockedMap);
      bCopyResultX = true;
      bReorderX    = true;
    }

    // check type of B vector
    if (bB.is_null() == true) {
      RCP<const MultiVector> vectorWithBlockedMap = Teuchos::rcp(new BlockedMultiVector(rbA->getBlockedCrsMatrix()->getBlockedRangeMap(), rcpB));
      rcpB.swap(vectorWithBlockedMap);
      bReorderB = true;
    }
  } else {
    // A is a BlockedCrsMatrix and uses a plain blocked map
    if (bX.is_null() == true) {
      RCP<MultiVector> vectorWithBlockedMap = Teuchos::rcp(new BlockedMultiVector(bA->getBlockedDomainMap(), rcpX));
      rcpX.swap(vectorWithBlockedMap);
      bCopyResultX = true;
    }

    if (bB.is_null() == true) {
      RCP<const MultiVector> vectorWithBlockedMap = Teuchos::rcp(new BlockedMultiVector(bA->getBlockedRangeMap(), rcpB));
      rcpB.swap(vectorWithBlockedMap);
    }
  }

  // we now can guarantee that X and B are blocked multi vectors
  bX = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(rcpX);
  bB = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(rcpB);

  // Finally we can do a reordering of the blocked multivectors if A is a ReorderedBlockedCrsMatrix
  if (rbA.is_null() == false) {
    Teuchos::RCP<const Xpetra::BlockReorderManager> brm = rbA->getBlockReorderManager();

    // check type of X vector
    if (bX->getBlockedMap()->getNumMaps() != bA->getDomainMapExtractor()->NumMaps() || bReorderX) {
      // X is a blocked multi vector but incompatible to the reordered blocked operator A
      Teuchos::RCP<MultiVector> reorderedBlockedVector = buildReorderedBlockedMultiVector(brm, bX);
      rcpX.swap(reorderedBlockedVector);
    }

    if (bB->getBlockedMap()->getNumMaps() != bA->getRangeMapExtractor()->NumMaps() || bReorderB) {
      // B is a blocked multi vector but incompatible to the reordered blocked operator A
      Teuchos::RCP<const MultiVector> reorderedBlockedVector = buildReorderedBlockedMultiVector(brm, bB);
      rcpB.swap(reorderedBlockedVector);
    }
  }

  // use the GIDs of the sub blocks
  // This is valid as the subblocks actually represent the GIDs (either Thyra, Xpetra or pseudo Xpetra)

  bool bRangeThyraMode  = rangeMapExtractor_->getThyraMode();
  bool bDomainThyraMode = domainMapExtractor_->getThyraMode();

  RCP<MultiVector> deltaX         = MultiVectorFactory::Build(rcpX->getMap(), rcpX->getNumVectors());
  RCP<BlockedMultiVector> bdeltaX = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(deltaX);
  RCP<MultiVector> deltaX0        = bdeltaX->getMultiVector(0, bDomainThyraMode);
  RCP<MultiVector> deltaX1        = bdeltaX->getMultiVector(1, bDomainThyraMode);

  RCP<MultiVector> Rtmp = rangeMapExtractor_->getVector(1, rcpB->getNumVectors(), bRangeThyraMode);

  typedef Teuchos::ScalarTraits<SC> STS;
  SC one = STS::one(), zero = STS::zero();

  // extract parameters from internal parameter list
  const ParameterList& pL = Factory::GetParameterList();
  LO nSweeps              = pL.get<LO>("Sweeps");

  RCP<MultiVector> R;  // = MultiVectorFactory::Build(rcpB->getMap(), rcpB->getNumVectors());;
  if (InitialGuessIsZero) {
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
    RCP<MultiVector> R0 = rangeMapExtractor_->ExtractVector(R, 0, bRangeThyraMode);
    RCP<MultiVector> R1 = rangeMapExtractor_->ExtractVector(R, 1, bRangeThyraMode);

    // Calculate Rtmp = R1 - D * deltaX0 (equation 8.14)
    deltaX0->putScalar(zero);
    deltaX0->elementWiseMultiply(one, *D_, *R0, zero);  // deltaX0 = D * R0 (equation 8.13)
    A10_->apply(*deltaX0, *Rtmp);                       // Rtmp    = A10*D*deltaX0 (intermediate step)
    Rtmp->update(one, *R1, -one);                       // Rtmp    = R1 - A10*D*deltaX0

    if (!pL.get<bool>("q2q1 mode")) {
      deltaX1->putScalar(zero);
    } else {
      // special code for q2q1
      if (Teuchos::rcp_dynamic_cast<BlockedVector>(diagSVector) == Teuchos::null) {
        ArrayRCP<SC> Sdiag       = diagSVector->getDataNonConst(0);
        ArrayRCP<SC> deltaX1data = deltaX1->getDataNonConst(0);
        ArrayRCP<SC> Rtmpdata    = Rtmp->getDataNonConst(0);
        for (GO row = 0; row < deltaX1data.size(); row++)
          deltaX1data[row] = Teuchos::as<SC>(1.1) * Rtmpdata[row] / Sdiag[row];
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "MueLu::BraessSarazinSmoother: q2q1 mode only supported for non-blocked operators.")
      }
    }

    smoo_->Apply(*deltaX1, *Rtmp);

    // Compute deltaX0
    deltaX0->putScalar(zero);         // just for safety
    A01_->apply(*deltaX1, *deltaX0);  // deltaX0 = A01*deltaX1
    deltaX0->update(one, *R0, -one);  // deltaX0 = R0 - A01*deltaX1
    R0.swap(deltaX0);
    deltaX0->elementWiseMultiply(one, *D_, *R0, zero);  // deltaX0 = D*(R0 - A01*deltaX1)

    RCP<MultiVector> X0 = domainMapExtractor_->ExtractVector(rcpX, 0, bDomainThyraMode);
    RCP<MultiVector> X1 = domainMapExtractor_->ExtractVector(rcpX, 1, bDomainThyraMode);

    // Update solution
    X0->update(one, *deltaX0, one);
    X1->update(one, *deltaX1, one);

    domainMapExtractor_->InsertVector(X0, 0, rcpX, bDomainThyraMode);
    domainMapExtractor_->InsertVector(X1, 1, rcpX, bDomainThyraMode);

    if (run < nSweeps - 1) {
      R = Utilities::Residual(*A_, *rcpX, *rcpB);
    }
  }

  if (bCopyResultX == true) {
    RCP<MultiVector> Xmerged = bX->Merge();
    X.update(one, *Xmerged, zero);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Copy() const {
  return rcp(new BraessSarazinSmoother(*this));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    description() const {
  std::ostringstream out;
  out << SmootherPrototype::description();
  out << "{type = " << type_ << "}";
  return out.str();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream& out, const VerbLevel verbLevel) const {
  MUELU_DESCRIBE;

  if (verbLevel & Parameters0) {
    out0 << "Prec. type: " << type_ << /*" Sweeps: " << nSweeps_ << " damping: " << omega_ <<*/ std::endl;
  }

  if (verbLevel & Debug) {
    out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl;
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNodeSmootherComplexity() const {
  // FIXME: This is a placeholder
  return Teuchos::OrdinalTraits<size_t>::invalid();
}

}  // namespace MueLu

#endif /* MUELU_BRAESSSARAZINSMOOTHER_DEF_HPP_ */
