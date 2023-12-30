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

#ifndef MUELU_UZAWASMOOTHER_DEF_HPP_
#define MUELU_UZAWASMOOTHER_DEF_HPP_

#include <Teuchos_ArrayViewDecl.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include "MueLu_ConfigDefs.hpp"

#include <Xpetra_Matrix.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_ReorderedBlockedCrsMatrix.hpp>

#include "MueLu_UzawaSmoother_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_HierarchyUtils.hpp"
#include "MueLu_SmootherBase.hpp"

#include "MueLu_FactoryManager.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> UzawaSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase>>("A", Teuchos::null, "Generating factory of the matrix A");
  validParamList->set<Scalar>("Damping factor", 1.0, "Damping/Scaling factor in SIMPLE");
  validParamList->set<LocalOrdinal>("Sweeps", 1, "Number of SIMPLE sweeps (default = 1)");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UzawaSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddFactoryManager(RCP<const FactoryManagerBase> FactManager, int pos) {
  TEUCHOS_TEST_FOR_EXCEPTION(pos < 0, Exceptions::RuntimeError, "MueLu::UzawaSmoother::AddFactoryManager: parameter \'pos\' must not be negative! error.");

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
void UzawaSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetVelocityPredictionFactoryManager(RCP<FactoryManager> FactManager) {
  AddFactoryManager(FactManager, 0);  // overwrite factory manager for predicting the primary variable
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UzawaSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetSchurCompFactoryManager(RCP<FactoryManager> FactManager) {
  AddFactoryManager(FactManager, 1);  // overwrite factory manager for SchurComplement
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UzawaSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
  currentLevel.DeclareInput("A", this->GetFactory("A").get());

  TEUCHOS_TEST_FOR_EXCEPTION(FactManager_.size() != 2, Exceptions::RuntimeError, "MueLu::UzawaSmoother::DeclareInput: You have to declare two FactoryManagers with a \"Smoother\" object: One for predicting the primary variable and one for the SchurComplement system. The smoother for the SchurComplement system needs a SchurComplementFactory as input for variable \"A\". make sure that you use the same proper damping factors for omega both in the SchurComplementFactory and in the SIMPLE smoother!");

  // loop over all factory managers for the subblocks of blocked operator A
  std::vector<Teuchos::RCP<const FactoryManagerBase>>::const_iterator it;
  for (it = FactManager_.begin(); it != FactManager_.end(); ++it) {
    SetFactoryManager currentSFM(rcpFromRef(currentLevel), *it);

    // request "Smoother" for current subblock row.
    currentLevel.DeclareInput("PreSmoother", (*it)->GetFactory("Smoother").get());
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UzawaSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level &currentLevel) {
  FactoryMonitor m(*this, "Setup blocked Uzawa Smoother", currentLevel);

  if (SmootherPrototype::IsSetup() == true)
    this->GetOStream(Warnings0) << "MueLu::UzawaSmoother::Setup(): Setup() has already been called";

  // extract blocked operator A from current level
  A_ = Factory::Get<RCP<Matrix>>(currentLevel, "A");

  RCP<BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A_);
  TEUCHOS_TEST_FOR_EXCEPTION(bA == Teuchos::null, Exceptions::BadCast, "MueLu::UzawaSmoother::Setup: input matrix A is not of type BlockedCrsMatrix! error.");

  // store map extractors
  rangeMapExtractor_  = bA->getRangeMapExtractor();
  domainMapExtractor_ = bA->getDomainMapExtractor();

  // Store the blocks in local member variables
  F_ = bA->getMatrix(0, 0);
  G_ = bA->getMatrix(0, 1);
  D_ = bA->getMatrix(1, 0);
  Z_ = bA->getMatrix(1, 1);

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
void UzawaSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Apply(MultiVector &X, const MultiVector &B, bool /* InitialGuessIsZero */) const {
  TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::UzawaSmoother::Apply(): Setup() has not been called");

  Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));

  SC zero = Teuchos::ScalarTraits<SC>::zero(), one = Teuchos::ScalarTraits<SC>::one();

  bool bRangeThyraMode  = rangeMapExtractor_->getThyraMode();
  bool bDomainThyraMode = domainMapExtractor_->getThyraMode();

  // extract parameters from internal parameter list
  const ParameterList &pL = Factory::GetParameterList();
  LocalOrdinal nSweeps    = pL.get<LocalOrdinal>("Sweeps");
  Scalar omega            = pL.get<Scalar>("Damping factor");

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
  Teuchos::RCP<MultiVector> r1      = bresidual->getMultiVector(0, bRangeThyraMode);
  Teuchos::RCP<MultiVector> r2      = bresidual->getMultiVector(1, bRangeThyraMode);

  // helper vector 1
  RCP<MultiVector> xtilde         = MultiVectorFactory::Build(rcpX->getMap(), rcpX->getNumVectors());
  RCP<BlockedMultiVector> bxtilde = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(xtilde);
  RCP<MultiVector> xtilde1        = bxtilde->getMultiVector(0, bDomainThyraMode);
  RCP<MultiVector> xtilde2        = bxtilde->getMultiVector(1, bDomainThyraMode);

  // incrementally improve solution vector X
  for (LocalOrdinal run = 0; run < nSweeps; ++run) {
    // 1) calculate current residual
    residual->update(one, *rcpB, zero);  // residual = B
    A_->apply(*rcpX, *residual, Teuchos::NO_TRANS, -one, one);

    // 2) solve F * \Delta \tilde{x}_1 = r_1
    //    start with zero guess \Delta \tilde{x}_1
    bxtilde->putScalar(zero);
    velPredictSmoo_->Apply(*xtilde1, *r1);

    // 3) calculate rhs for SchurComp equation
    //    r_2 - D \Delta \tilde{x}_1
    RCP<MultiVector> schurCompRHS = MultiVectorFactory::Build(r2->getMap(), rcpB->getNumVectors());
    D_->apply(*xtilde1, *schurCompRHS);
    schurCompRHS->update(one, *r2, -one);

    // 4) solve SchurComp equation
    //    start with zero guess \Delta \tilde{x}_2
    schurCompSmoo_->Apply(*xtilde2, *schurCompRHS);

    rcpX->update(omega, *bxtilde, one);
  }

  if (bCopyResultX == true) {
    RCP<MultiVector> Xmerged = bX->Merge();
    X.update(one, *Xmerged, zero);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
UzawaSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Copy() const {
  return rcp(new UzawaSmoother(*this));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string UzawaSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
  std::ostringstream out;
  out << SmootherPrototype::description();
  out << "{type = " << type_ << "}";
  return out.str();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void UzawaSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
  MUELU_DESCRIBE;

  if (verbLevel & Parameters0) {
    out0 << "Prec. type: " << type_ << /*" Sweeps: " << nSweeps_ << " damping: " << omega_ <<*/ std::endl;
  }

  if (verbLevel & Debug) {
    out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl;
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t UzawaSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNodeSmootherComplexity() const {
  // FIXME: This is a placeholder
  return Teuchos::OrdinalTraits<size_t>::invalid();
}

}  // namespace MueLu

#endif /* MUELU_UZAWASMOOTHER_DEF_HPP_ */
