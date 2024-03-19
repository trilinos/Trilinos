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
 * MueLu_BlockedDirectSolver_def.hpp
 *
 *  Created on: 09.02.2014
 *      Author: tobias
 */

#ifndef MUELU_BLOCKEDDIRECTSOLVER_DEF_HPP_
#define MUELU_BLOCKEDDIRECTSOLVER_DEF_HPP_

#include "Teuchos_ArrayViewDecl.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "MueLu_ConfigDefs.hpp"

#include <Xpetra_Matrix.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_BlockedDirectSolver_decl.hpp"
#include "MueLu_MergedBlockedMatrixFactory.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_DirectSolver.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedDirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BlockedDirectSolver(const std::string &type, const Teuchos::ParameterList &paramList) {
  MergedAFact_ = Teuchos::rcp(new MergedBlockedMatrixFactory());
  s_           = Teuchos::rcp(new DirectSolver(type, paramList));
  type_        = "blocked direct solver (" + type + ")";
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> BlockedDirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase> >("A", null, "Generating factory of the matrix A");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedDirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
  // Note that we have a nested smoother/solver object (of type DirectSolver), so we have to declare the dependencies by hand
  // call DeclareInput by hand, since this->Input(currentLevel, "A") would not properly free A in the release mode (dependencies)
  // We need the blocked version of A as input for the MergedAFact_
  currentLevel.DeclareInput("A", this->GetFactory("A").get());

  // syncronize input factory for "A" and nested representation for "A"
  MergedAFact_->SetFactory("A", this->GetFactory("A"));

  // declare input factories for nested direct solver
  s_->SetFactory("A", MergedAFact_);
  s_->DeclareInput(currentLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedDirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level &currentLevel) {
  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  FactoryMonitor m(*this, "Setup BlockedDirectSolver", currentLevel);
  if (this->IsSetup() == true)
    this->GetOStream(Warnings0) << "MueLu::BlockedDirectSolver::Setup(): Setup() has already been called";

  // extract blocked operator A from current level
  A_                       = Factory::Get<RCP<Matrix> >(currentLevel, "A");  // A needed for extracting map extractors
  RCP<BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A_);
  TEUCHOS_TEST_FOR_EXCEPTION(bA.is_null(), Exceptions::BadCast,
                             "MueLu::BlockedDirectSolver::Build: input matrix A is not of type BlockedCrsMatrix.");

  s_->Setup(currentLevel);

  this->IsSetup(true);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedDirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Apply(MultiVector &X, const MultiVector &B, bool InitialGuessIsZero) const {
  TEUCHOS_TEST_FOR_EXCEPTION(this->IsSetup() == false, Exceptions::RuntimeError,
                             "MueLu::BlockedDirectSolver::Apply(): Setup() has not been called");

  RCP<MultiVector> rcpX            = Teuchos::rcpFromRef(X);
  RCP<const MultiVector> rcpB      = Teuchos::rcpFromRef(B);
  RCP<BlockedMultiVector> bX       = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(rcpX);
  RCP<const BlockedMultiVector> bB = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(rcpB);

#ifdef HAVE_MUELU_DEBUG
  RCP<BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A_);
  if (bB.is_null() == false) {
    // TEUCHOS_TEST_FOR_EXCEPTION(A_->getRangeMap()->isSameAs(*(B.getMap())) == false, Exceptions::RuntimeError, "MueLu::BlockedDirectSolver::Apply(): The map of RHS vector B is not the same as range map of the blocked operator A. Please check the map of B and A.");
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(bA->getFullRangeMap()->isSameAs(*(B.getMap())) == false, Exceptions::RuntimeError, "MueLu::BlockedDirectSolver::Apply(): The map of RHS vector B is not the same as range map of the blocked operator A. Please check the map of B and A.");
  }
  if (bX.is_null() == false) {
    // TEUCHOS_TEST_FOR_EXCEPTION(A_->getDomainMap()->isSameAs(*(X.getMap())) == false, Exceptions::RuntimeError, "MueLu::BlockedDirectSolver::Apply(): The map of the solution vector X is not the same as domain map of the blocked operator A. Please check the map of X and A.");
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(bA->getFullDomainMap()->isSameAs(*(X.getMap())) == false, Exceptions::RuntimeError, "MueLu::BlockedDirectSolver::Apply(): The map of the solution vector X is not the same as domain map of the blocked operator A. Please check the map of X and A.");
  }
#endif

  if (bB.is_null() == true && bX.is_null() == true) {
    // standard case (neither B nor X are blocked)
    s_->Apply(X, B, InitialGuessIsZero);
  } else if (bB.is_null() == false && bX.is_null() == false) {
    // both B and X are blocked
    RCP<MultiVector> mergedX       = bX->Merge();
    RCP<const MultiVector> mergedB = bB->Merge();
    s_->Apply(*mergedX, *mergedB, InitialGuessIsZero);
    RCP<MultiVector> xx = Teuchos::rcp(new BlockedMultiVector(bX->getBlockedMap(), mergedX));
    SC zero = Teuchos::ScalarTraits<SC>::zero(), one = Teuchos::ScalarTraits<SC>::one();
    X.update(one, *xx, zero);
  } else if (bB.is_null() == true && bX.is_null() == false) {
    // the solution vector is blocked
    RCP<MultiVector> mergedX = bX->Merge();
    s_->Apply(*mergedX, B, InitialGuessIsZero);
    RCP<MultiVector> xx = Teuchos::rcp(new BlockedMultiVector(bX->getBlockedMap(), mergedX));
    SC zero = Teuchos::ScalarTraits<SC>::zero(), one = Teuchos::ScalarTraits<SC>::one();
    X.update(one, *xx, zero);
  } else if (bB.is_null() == false && bX.is_null() == true) {
    // only the RHS vector is blocked
    RCP<const MultiVector> mergedB = bB->Merge();
    s_->Apply(X, *mergedB, InitialGuessIsZero);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
BlockedDirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Copy() const {
  return rcp(new BlockedDirectSolver(*this));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string BlockedDirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
  std::ostringstream out;
  out << SmootherPrototype::description();
  out << "{type = " << type_ << "}";
  return out.str();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedDirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
  MUELU_DESCRIBE;

  if (verbLevel & Parameters0)
    out0 << "Prec. type: " << type_ << std::endl;

  if (verbLevel & Debug)
    out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t BlockedDirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNodeSmootherComplexity() const {
  // FIXME: This is a placeholder
  return Teuchos::OrdinalTraits<size_t>::invalid();
}

}  // namespace MueLu

#endif /* MUELU_BLOCKEDDIRECTSOLVER_DEF_HPP_ */
