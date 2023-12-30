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
 * MueLu_PermutingSmoother_def.hpp
 *
 *  Created on: Nov 28, 2012
 *      Author: wiesner
 */

#ifndef MUELU_PERMUTINGSMOOTHER_DEF_HPP
#define MUELU_PERMUTINGSMOOTHER_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_PermutingSmoother_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_IfpackSmoother.hpp"
#include "MueLu_PermutationFactory.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
PermutingSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::PermutingSmoother(const std::string& mapName, const RCP<const FactoryBase>& mapFact, const std::string& type, const Teuchos::ParameterList& paramList, const LO& overlap, RCP<FactoryBase> permFact)
  : type_(type)
  , overlap_(overlap)
  , permQT_(Teuchos::null)
  , permP_(Teuchos::null)
  , diagScalingOp_(Teuchos::null) {
  this->SetParameterList(paramList);

  permFact_ = permFact;
  if (permFact_ == Teuchos::null) {
    RCP<PermutationFactory> newPermFact = Teuchos::rcp(new PermutationFactory());
    newPermFact->SetParameter("PermutationRowMapName", Teuchos::ParameterEntry(mapName));
    newPermFact->SetFactory("PermutationRowMapFactory", mapFact);
    permFact_ = newPermFact;
  }

  // create internal smoother
  if (type_ == "ILU") {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_IFPACK)
    s_ = MueLu::GetIfpackSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>(type_, this->GetParameterList(), overlap_);
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::PermutingSmoother requires Epetra and Ifpack.");
#endif
  } else {
    s_ = Teuchos::rcp(new TrilinosSmoother(type_, this->GetParameterList(), overlap_));
  }
  TEUCHOS_TEST_FOR_EXCEPTION(s_ == Teuchos::null, Exceptions::RuntimeError, "");

  // Use permuted matrix A
  s_->SetFactory("A", permFact_);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
PermutingSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~PermutingSmoother() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PermutingSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  currentLevel.DeclareInput("permP", permFact_.get());
  currentLevel.DeclareInput("permQT", permFact_.get());
  currentLevel.DeclareInput("permScaling", permFact_.get());

  s_->DeclareInput(currentLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PermutingSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level& currentLevel) {
  FactoryMonitor monitor(*this, "Permuting Smoother", currentLevel);

  if (SmootherPrototype::IsSetup() == true)
    this->GetOStream(Warnings0) << "MueLu::PermutingSmoother::Setup(): Setup() has already been called" << std::endl;

  // extract information from level class
  permP_         = currentLevel.Get<RCP<Matrix> >("permP", permFact_.get());
  permQT_        = currentLevel.Get<RCP<Matrix> >("permQT", permFact_.get());
  diagScalingOp_ = currentLevel.Get<RCP<Matrix> >("permScaling", permFact_.get());

  s_->Setup(currentLevel);

  SmootherPrototype::IsSetup(true);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PermutingSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const {
  TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::PermutingSmoother::Apply(): Setup() has not been called");

  typedef Teuchos::ScalarTraits<Scalar> STS;

  Teuchos::RCP<MultiVector> Xtemp = MultiVectorFactory::Build(X.getMap(), 1, true);
  Xtemp->update(STS::one(), X, STS::zero());

  // TODO: unify scaling and left permutation operator
  Teuchos::RCP<MultiVector> Btemp  = MultiVectorFactory::Build(B.getMap(), 1, true);
  Teuchos::RCP<MultiVector> Btemp2 = MultiVectorFactory::Build(B.getMap(), 1, true);
  permP_->apply(B, *Btemp, Teuchos::NO_TRANS);                // apply permutation operator to rhs
  diagScalingOp_->apply(*Btemp, *Btemp2, Teuchos::NO_TRANS);  // apply scaling operator to rhs

  // apply smoother to permuted linear system
  s_->Apply(*Xtemp, *Btemp2, InitialGuessIsZero);

  // retransform smooth solution
  permQT_->apply(*Xtemp, X, Teuchos::NO_TRANS);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
PermutingSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Copy() const {
  return rcp(new PermutingSmoother(*this));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string PermutingSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
  std::ostringstream out;
  out << SmootherPrototype::description();
  return out.str();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PermutingSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream& out, const VerbLevel verbLevel) const {
  MUELU_DESCRIBE;
  out0 << "";  // avoid warning
}

}  // namespace MueLu

#endif /* MUELU_PERMUTINGSMOOTHER_DEF_HPP */
