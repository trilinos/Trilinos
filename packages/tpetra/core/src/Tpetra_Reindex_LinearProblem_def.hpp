// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_REINDEX_LINEARPROBLEM_DEF_HPP
#define TPETRA_REINDEX_LINEARPROBLEM_DEF_HPP

/// \file Tpetra_Reindex_LinearProblem_def.hpp
/// \brief Definition of the Tpetra::Reindex_LinearProblem class
///
/// If you want to use Tpetra::Reindex_LinearProblem, include
/// "Tpetra_Reindex_LinearProblem.hpp", a file which CMake generates
/// and installs for you.
///
/// If you only want the declaration of Tpetra::Reindex_LinearProblem,
/// include "Tpetra_Reindex_LinearProblem_decl.hpp".

#include <Tpetra_Reindex_LinearProblem_decl.hpp>

namespace Tpetra {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Reindex_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Reindex_LinearProblem(Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> const> newRowMap)
  : ViewTransform<LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node> >()
  , newRowMap_(newRowMap)
  , matTrans_(Teuchos::null)
  , lhsTrans_(Teuchos::null)
  , rhsTrans_(Teuchos::null) {
  // Nothing to do
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Reindex_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~Reindex_LinearProblem() {
  // Nothing to do
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
typename Reindex_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NewType
Reindex_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::operator()(OriginalType const &origProblem) {
  using map_t = Map<LocalOrdinal, GlobalOrdinal, Node>;
  using mv_t  = MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using cm_t  = CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using lp_t  = LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  // Save original object
  this->origObj_ = origProblem;

  Teuchos::RCP<cm_t> origMatrix        = Teuchos::rcp<cm_t>(dynamic_cast<cm_t *>(origProblem->getMatrix().get()), false);
  Teuchos::RCP<mv_t> origRHS           = origProblem->getRHS();
  Teuchos::RCP<mv_t> origLHS           = origProblem->getLHS();
  Teuchos::RCP<map_t const> origRowMap = origMatrix->getMap();

  // If no new map has been passed in, create one
  if (newRowMap_.get() == nullptr) {
    newRowMap_ = Teuchos::rcp<map_t const>(new map_t(origRowMap->getGlobalNumElements(), origRowMap->getLocalNumElements(), 0, origRowMap->getComm()));
  }

  using r_cm_t = Reindex_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using r_mv_t = Reindex_MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  matTrans_ = Teuchos::rcp<r_cm_t>(new r_cm_t(newRowMap_));
  lhsTrans_ = Teuchos::rcp<r_mv_t>(new r_mv_t(newRowMap_));
  rhsTrans_ = Teuchos::rcp<r_mv_t>(new r_mv_t(newRowMap_));

  Teuchos::RCP<cm_t> newMatrix = ((*matTrans_)(origMatrix));
  Teuchos::RCP<mv_t> newLHS    = ((*lhsTrans_)(origLHS));
  Teuchos::RCP<mv_t> newRHS    = ((*rhsTrans_)(origRHS));

  this->newObj_ = Teuchos::rcp<lp_t>(new lp_t(newMatrix, newLHS, newRHS));

  return this->newObj_;
}

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_REINDEXLINEARPROBLEM_INSTANT(SCALAR, LO, GO, NODE) \
  template class Reindex_LinearProblem<SCALAR, LO, GO, NODE>;

}  // namespace Tpetra

#endif  // TPETRA_REINDEX_LINEARPROBLEM_DEF_HPP
