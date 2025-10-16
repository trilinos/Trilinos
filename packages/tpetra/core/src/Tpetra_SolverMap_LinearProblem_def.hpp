// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_SOLVERMAP_LINEARPROBLEM_DEF_HPP
#define TPETRA_SOLVERMAP_LINEARPROBLEM_DEF_HPP

/// \file Tpetra_SolverMap_LinearProblem_def.hpp
/// \brief Definition of the Tpetra::SolverMap_LinearProblem class
///
/// If you want to use Tpetra::SolverMap_LinearProblem, include
/// "Tpetra_SolverMap_LinearProblem.hpp", a file which CMake generates
/// and installs for you.
///
/// If you only want the declaration of Tpetra::SolverMap_LinearProblem,
/// include "Tpetra_SolverMap_LinearProblem_decl.hpp".

#include <Tpetra_SolverMap_LinearProblem_decl.hpp>

namespace Tpetra {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
SolverMap_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SolverMap_LinearProblem()
  : StructuralSameTypeTransform<LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node> >()
  , solverMapCrsMatrixTrans_() {
  // Nothing to do
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
SolverMap_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~SolverMap_LinearProblem() {
  // Nothing to do
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
typename SolverMap_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NewType
SolverMap_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
operator()(OriginalType const &origProblem) {
  using mv_t = MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using cm_t = CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using lp_t = LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  this->origObj_ = origProblem;

  cm_t *oldMatrix              = dynamic_cast<cm_t *>(origProblem->getMatrix().get());
  Teuchos::RCP<mv_t> oldRHS    = origProblem->getRHS();
  Teuchos::RCP<mv_t> oldLHS    = origProblem->getLHS();
  Teuchos::RCP<cm_t> newMatrix = solverMapCrsMatrixTrans_(Teuchos::rcp<cm_t>(oldMatrix, false));

  if (newMatrix.get() == oldMatrix) {
    // Same matrix, so use same problem
    this->newObj_ = this->origObj_;
  } else {
    this->newObj_ = Teuchos::rcp<lp_t>(new lp_t(newMatrix, oldLHS, oldRHS));
  }

  return this->newObj_;
}

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_SOLVERMAPLINEARPROBLEM_INSTANT(SCALAR, LO, GO, NODE) \
  template class SolverMap_LinearProblem<SCALAR, LO, GO, NODE>;

}  // namespace Tpetra

#endif  // TPETRA_SOLVERMAP_LINEARPROBLEM_DEF_HPP
