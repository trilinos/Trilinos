// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
#include <Xpetra_ReorderedBlockedCrsMatrix.hpp>

#include "MueLu_SimpleSmoother_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_HierarchyUtils.hpp"
#include "MueLu_SmootherBase.hpp"

// include files for default FactoryManager
#include "MueLu_SchurComplementFactory.hpp"
#include "MueLu_FactoryManager.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SimpleSmoother()
  : type_("SIMPLE")
  , A_(Teuchos::null) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~SimpleSmoother() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase>>("A", Teuchos::null, "Generating factory of the matrix A");
  validParamList->set<Scalar>("Damping factor", 1.0, "Damping/Scaling factor in SIMPLE");
  validParamList->set<LocalOrdinal>("Sweeps", 1, "Number of SIMPLE sweeps (default = 1)");
  validParamList->set<bool>("UseSIMPLEC", false, "Use SIMPLEC instead of SIMPLE (default = false)");

  return validParamList;
}

//! Add a factory manager at a specific position
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddFactoryManager(RCP<const FactoryManagerBase> FactManager, int pos) {
  TEUCHOS_TEST_FOR_EXCEPTION(pos < 0, Exceptions::RuntimeError, "MueLu::SimpleSmoother::AddFactoryManager: parameter \'pos\' must not be negative! error.");

  size_t myPos = Teuchos::as<size_t>(pos);

  if (myPos < FactManager_.size()) {
    // replace existing entris in FactManager_ vector
    FactManager_.at(myPos) = FactManager;
  } else if (myPos == FactManager_.size()) {
    // add new Factory manager in the end of the vector
    FactManager_.push_back(FactManager);
  } else {  // if(myPos > FactManager_.size())
    RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    *out << "Warning: cannot add new FactoryManager at proper position " << pos << ". The FactoryManager is just appended to the end. Check this!" << std::endl;

    // add new Factory manager in the end of the vector
    FactManager_.push_back(FactManager);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetVelocityPredictionFactoryManager(RCP<FactoryManager> FactManager) {
  AddFactoryManager(FactManager, 0);  // overwrite factory manager for predicting the primary variable
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetSchurCompFactoryManager(RCP<FactoryManager> FactManager) {
  AddFactoryManager(FactManager, 1);  // overwrite factory manager for SchurComplement
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
  currentLevel.DeclareInput("A", this->GetFactory("A").get());

  TEUCHOS_TEST_FOR_EXCEPTION(FactManager_.size() != 2, Exceptions::RuntimeError, "MueLu::SimpleSmoother::DeclareInput: You have to declare two FactoryManagers with a \"Smoother\" object: One for predicting the primary variable and one for the SchurComplement system. The smoother for the SchurComplement system needs a SchurComplementFactory as input for variable \"A\". make sure that you use the same proper damping factors for omega both in the SchurComplementFactory and in the SIMPLE smoother!");

  // loop over all factory managers for the subblocks of blocked operator A
  std::vector<Teuchos::RCP<const FactoryManagerBase>>::const_iterator it;
  for (it = FactManager_.begin(); it != FactManager_.end(); ++it) {
    SetFactoryManager currentSFM(rcpFromRef(currentLevel), *it);

    // request "Smoother" for current subblock row.
    currentLevel.DeclareInput("PreSmoother", (*it)->GetFactory("Smoother").get());
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level &currentLevel) {
  FactoryMonitor m(*this, "Setup blocked SIMPLE Smoother", currentLevel);

  if (SmootherPrototype::IsSetup() == true)
    this->GetOStream(Warnings0) << "MueLu::SimpleSmoother::Setup(): Setup() has already been called";

  // extract blocked operator A from current level
  A_ = Factory::Get<RCP<Matrix>>(currentLevel, "A");

  RCP<BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A_);
  TEUCHOS_TEST_FOR_EXCEPTION(bA == Teuchos::null, Exceptions::BadCast, "MueLu::SimpleSmoother::Setup: input matrix A is not of type BlockedCrsMatrix! error.");

  // store map extractors
  rangeMapExtractor_  = bA->getRangeMapExtractor();
  domainMapExtractor_ = bA->getDomainMapExtractor();

  // Store the blocks in local member variables
  F_ = bA->getMatrix(0, 0);
  G_ = bA->getMatrix(0, 1);
  D_ = bA->getMatrix(1, 0);
  Z_ = bA->getMatrix(1, 1);

  const ParameterList &pL = Factory::GetParameterList();
  bool bSIMPLEC           = pL.get<bool>("UseSIMPLEC");

  // Create the inverse of the diagonal of F
  // TODO add safety check for zeros on diagonal of F!
  RCP<Vector> diagFVector = VectorFactory::Build(F_->getRowMap());
  if (!bSIMPLEC) {
    F_->getLocalDiagCopy(*diagFVector);  // extract diagonal of F
  } else {
    /*const RCP<const Map> rowmap = F_->getRowMap();
    size_t locSize = rowmap->getLocalNumElements();
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
    }*/
    // TODO this does not work if F_ is nested!
    diagFVector = Utilities::GetLumpedMatrixDiagonal(*F_);
  }
  diagFinv_ = Utilities::GetInverse(diagFVector);

  // check whether diagFinv_ is a blocked vector with only 1 block
  RCP<BlockedVector> bdiagFinv = Teuchos::rcp_dynamic_cast<BlockedVector>(diagFinv_);
  if (bdiagFinv.is_null() == false && bdiagFinv->getBlockedMap()->getNumMaps() == 1) {
    RCP<Vector> nestedVec = bdiagFinv->getMultiVector(0, bdiagFinv->getBlockedMap()->getThyraMode())->getVectorNonConst(0);
    diagFinv_.swap(nestedVec);
  }

  // Set the Smoother
  // carefully switch to the SubFactoryManagers (defined by the users)
  {
    RCP<const FactoryManagerBase> velpredictFactManager = FactManager_.at(0);
    SetFactoryManager currentSFM(rcpFromRef(currentLevel), velpredictFactManager);
    velPredictSmoo_ = currentLevel.Get<RCP<SmootherBase>>("PreSmoother", velpredictFactManager->GetFactory("Smoother").get());
  }
  {
    RCP<const FactoryManagerBase> schurFactManager = FactManager_.at(1);
    SetFactoryManager currentSFM(rcpFromRef(currentLevel), schurFactManager);
    schurCompSmoo_ = currentLevel.Get<RCP<SmootherBase>>("PreSmoother", schurFactManager->GetFactory("Smoother").get());
  }

  SmootherPrototype::IsSetup(true);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Apply(MultiVector &X, const MultiVector &B, bool InitialGuessIsZero) const {
  TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError,
                             "MueLu::SimpleSmoother::Apply(): Setup() has not been called");
#if 0
    // TODO simplify this debug check
    RCP<MultiVector> rcpDebugX = Teuchos::rcpFromRef(X);
    RCP<const MultiVector> rcpDebugB = Teuchos::rcpFromRef(B);
    RCP<BlockedMultiVector> rcpBDebugX = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(rcpDebugX);
    RCP<const BlockedMultiVector> rcpBDebugB = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(rcpDebugB);
    if(rcpBDebugB.is_null() == false) {
      //TEUCHOS_TEST_FOR_EXCEPTION(A_->getRangeMap()->isSameAs(*(B.getMap())) == false, Exceptions::RuntimeError, "MueLu::SimpleSmoother::Apply(): The map of RHS vector B is not the same as range map of the blocked operator A. Please check the map of B and A.");
    } else {
      //TEUCHOS_TEST_FOR_EXCEPTION(bA->getFullRangeMap()->isSameAs(*(B.getMap())) == false, Exceptions::RuntimeError, "MueLu::SimpleSmoother::Apply(): The map of RHS vector B is not the same as range map of the blocked operator A. Please check the map of B and A.");
    }
    if(rcpBDebugX.is_null() == false) {
      //TEUCHOS_TEST_FOR_EXCEPTION(A_->getDomainMap()->isSameAs(*(X.getMap())) == false, Exceptions::RuntimeError, "MueLu::SimpleSmoother::Apply(): The map of the solution vector X is not the same as domain map of the blocked operator A. Please check the map of X and A.");
    } else {
      //TEUCHOS_TEST_FOR_EXCEPTION(bA->getFullDomainMap()->isSameAs(*(X.getMap())) == false, Exceptions::RuntimeError, "MueLu::SimpleSmoother::Apply(): The map of the solution vector X is not the same as domain map of the blocked operator A. Please check the map of X and A.");
    }
#endif

  const SC zero = Teuchos::ScalarTraits<SC>::zero();
  const SC one  = Teuchos::ScalarTraits<SC>::one();

  // extract parameters from internal parameter list
  const ParameterList &pL = Factory::GetParameterList();
  LocalOrdinal nSweeps    = pL.get<LocalOrdinal>("Sweeps");
  Scalar omega            = pL.get<Scalar>("Damping factor");

  // The boolean flags check whether we use Thyra or Xpetra style GIDs
  bool bRangeThyraMode  = rangeMapExtractor_->getThyraMode();   //  && (Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(F_) == Teuchos::null);
  bool bDomainThyraMode = domainMapExtractor_->getThyraMode();  // && (Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(F_) == Teuchos::null);

  // RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  // wrap current solution vector in RCP
  RCP<MultiVector> rcpX       = Teuchos::rcpFromRef(X);
  RCP<const MultiVector> rcpB = Teuchos::rcpFromRef(B);

  // make sure that both rcpX and rcpB are BlockedMultiVector objects
  bool bCopyResultX        = false;
  bool bReorderX           = false;
  bool bReorderB           = false;
  RCP<BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A_);
  MUELU_TEST_FOR_EXCEPTION(bA.is_null() == true, Exceptions::RuntimeError, "MueLu::BlockedGaussSeidelSmoother::Apply(): A_ must be a BlockedCrsMatrix");
  RCP<BlockedMultiVector> bX       = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(rcpX);
  RCP<const BlockedMultiVector> bB = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(rcpB);

  // check the type of operator
  RCP<Xpetra::ReorderedBlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> rbA =
      Teuchos::rcp_dynamic_cast<Xpetra::ReorderedBlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(bA);

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

  // Throughout the rest of the algorithm rcpX and rcpB are used for solution vector and RHS

  // create residual vector
  // contains current residual of current solution X with rhs B
  RCP<MultiVector> residual         = MultiVectorFactory::Build(rcpB->getMap(), rcpB->getNumVectors());
  RCP<BlockedMultiVector> bresidual = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(residual);
  RCP<MultiVector> r1               = bresidual->getMultiVector(0, bRangeThyraMode);
  RCP<MultiVector> r2               = bresidual->getMultiVector(1, bRangeThyraMode);

  // helper vector 1
  RCP<MultiVector> xtilde         = MultiVectorFactory::Build(rcpX->getMap(), rcpX->getNumVectors());
  RCP<BlockedMultiVector> bxtilde = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(xtilde);
  RCP<MultiVector> xtilde1        = bxtilde->getMultiVector(0, bDomainThyraMode);
  RCP<MultiVector> xtilde2        = bxtilde->getMultiVector(1, bDomainThyraMode);

  // helper vector 2
  RCP<MultiVector> xhat         = MultiVectorFactory::Build(rcpX->getMap(), rcpX->getNumVectors());
  RCP<BlockedMultiVector> bxhat = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(xhat);
  RCP<MultiVector> xhat1        = bxhat->getMultiVector(0, bDomainThyraMode);
  RCP<MultiVector> xhat2        = bxhat->getMultiVector(1, bDomainThyraMode);

  // incrementally improve solution vector X
  for (LocalOrdinal run = 0; run < nSweeps; ++run) {
    // 1) calculate current residual
    residual->update(one, *rcpB, zero);  // residual = B
    if (InitialGuessIsZero == false || run > 0)
      A_->apply(*rcpX, *residual, Teuchos::NO_TRANS, -one, one);

    // 2) solve F * \Delta \tilde{x}_1 = r_1
    //    start with zero guess \Delta \tilde{x}_1
    xtilde1->putScalar(zero);
    xtilde2->putScalar(zero);
    velPredictSmoo_->Apply(*xtilde1, *r1);

    // 3) calculate rhs for SchurComp equation
    //    r_2 - D \Delta \tilde{x}_1
    RCP<MultiVector> schurCompRHS = rangeMapExtractor_->getVector(1, rcpB->getNumVectors(), bRangeThyraMode);
    D_->apply(*xtilde1, *schurCompRHS);

    schurCompRHS->update(one, *r2, -one);

    // 4) solve SchurComp equation
    //    start with zero guess \Delta \tilde{x}_2
    schurCompSmoo_->Apply(*xtilde2, *schurCompRHS);

    // 5) scale xtilde2 with omega
    //    store this in xhat2
    xhat2->update(omega, *xtilde2, zero);

    // 6) calculate xhat1
    RCP<MultiVector> xhat1_temp = domainMapExtractor_->getVector(0, rcpX->getNumVectors(), bDomainThyraMode);
    G_->apply(*xhat2, *xhat1_temp);  // store result temporarely in xtilde1_temp

    xhat1->elementWiseMultiply(one /*/omega*/, *diagFinv_, *xhat1_temp, zero);
    xhat1->update(one, *xtilde1, -one);

    rcpX->update(one, *bxhat, one);
  }

  if (bCopyResultX == true) {
    RCP<const MultiVector> Xmerged = bX->Merge();
    X.update(one, *Xmerged, zero);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Copy() const {
  return rcp(new SimpleSmoother(*this));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
  std::ostringstream out;
  out << SmootherPrototype::description();
  out << "{type = " << type_ << "}";
  return out.str();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
  MUELU_DESCRIBE;

  if (verbLevel & Parameters0) {
    out0 << "Prec. type: " << type_ << /*" Sweeps: " << nSweeps_ << " damping: " << omega_ <<*/ std::endl;
  }

  if (verbLevel & Debug) {
    out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl;
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t SimpleSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNodeSmootherComplexity() const {
  // FIXME: This is a placeholder
  return Teuchos::OrdinalTraits<size_t>::invalid();
}

}  // namespace MueLu

#endif /* MUELU_SIMPLESMOOTHER_DEF_HPP_ */