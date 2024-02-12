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
#ifndef MUELU_PROJECTORSMOOTHER_DEF_HPP
#define MUELU_PROJECTORSMOOTHER_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_ProjectorSmoother_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ProjectorSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ProjectorSmoother(RCP<SmootherPrototype> coarseSolver)
  : coarseSolver_(coarseSolver) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ProjectorSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~ProjectorSmoother() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ProjectorSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
  Factory::Input(currentLevel, "A");
  Factory::Input(currentLevel, "Nullspace");

  coarseSolver_->DeclareInput(currentLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ProjectorSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level &currentLevel) {
  FactoryMonitor monitor(*this, "Projector Smoother", currentLevel);

  coarseSolver_->Setup(currentLevel);

  if (SmootherPrototype::IsSetup() == true)
    this->GetOStream(Warnings0) << "MueLu::ProjectorSmoother::Setup(): Setup() has already been called" << std::endl;

  RCP<Matrix> A      = Factory::Get<RCP<Matrix> >(currentLevel, "A");
  RCP<MultiVector> B = Factory::Get<RCP<MultiVector> >(currentLevel, "Nullspace");

  int m = B->getNumVectors();

  // Find which vectors we want to keep
  Array<Scalar> br(m), bb(m);
  RCP<MultiVector> R = MultiVectorFactory::Build(B->getMap(), m);
  A->apply(*B, *R);
  B->dot(*R, br);
  B->dot(*B, bb);

  Array<size_t> selectedIndices;
  for (int i = 0; i < m; i++) {
    Scalar rayleigh = br[i] / bb[i];

    if (Teuchos::ScalarTraits<Scalar>::magnitude(rayleigh) < 1e-12)
      selectedIndices.push_back(i);
  }
  this->GetOStream(Runtime0) << "Coarse level orth indices: " << selectedIndices << std::endl;

#if defined(HAVE_XPETRA_TPETRA)
#ifdef HAVE_MUELU_TPETRA_INST_INT_INT
  // Orthonormalize
  RCP<const Tpetra::MultiVector<SC, LO, GO, NO> > B_ = Utilities::MV2TpetraMV(B);
  // TAW: Oct 16 2015: subCopy is not part of Xpetra. One should either add it to Xpetra (with an emulator for Epetra)
  //                   or replace this call by a local loop. I'm not motivated to do this now...
  RCP<Tpetra::MultiVector<SC, LO, GO, NO> > Borth = B_->subCopy(selectedIndices);  // copy
  for (int i = 0; i < selectedIndices.size(); i++) {
    RCP<Tpetra::Vector<SC, LO, GO, NO> > Bi = Borth->getVectorNonConst(i);

    Scalar dot;
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm;
    for (int k = 0; k < i; k++) {  // orthogonalize
      RCP<const Tpetra::Vector<SC, LO, GO, NO> > Bk = Borth->getVector(k);

      dot = Bi->dot(*Bk);
      Bi->update(-dot, *Bk, Teuchos::ScalarTraits<Scalar>::one());
    }

    norm = Bi->norm2();
    Bi->scale(Teuchos::ScalarTraits<Scalar>::one() / norm);  // normalize
  }

  Borth_ = rcp(static_cast<MultiVector *>(new TpetraMultiVector(Borth)));
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Tpetra with GO=int not available. The code in ProjectorSmoother should be rewritten!");
#endif
#endif

  SmootherPrototype::IsSetup(true);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ProjectorSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Apply(MultiVector &X, const MultiVector &B, bool InitialGuessIsZero) const {
  coarseSolver_->Apply(X, B, InitialGuessIsZero);

  int m = Borth_->getNumVectors();
  int n = X.getNumVectors();

  RCP<Xpetra::MultiVector<SC, LO, GO, NO> > X_ = Teuchos::rcpFromRef(X);
  for (int i = 0; i < n; i++) {
    RCP<Xpetra::Vector<SC, LO, GO, NO> > Xi = X_->getVectorNonConst(i);

    Array<Scalar> dot(1);
    for (int k = 0; k < m; k++) {  // orthogonalize
      RCP<const Xpetra::Vector<SC, LO, GO, NO> > Bk = Borth_->getVector(k);

      Xi->dot(*Bk, dot());
      Xi->update(-dot[0], *Bk, Teuchos::ScalarTraits<Scalar>::one());
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> > ProjectorSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Copy() const {
  return rcp(new ProjectorSmoother(*this));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string ProjectorSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
  std::ostringstream out;
  out << SmootherPrototype::description();
  return out.str();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ProjectorSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
  MUELU_DESCRIBE;
  out0 << "";
}

}  // namespace MueLu

#endif  // MUELU_PROJECTORSMOOTHER_DEF_HPP
