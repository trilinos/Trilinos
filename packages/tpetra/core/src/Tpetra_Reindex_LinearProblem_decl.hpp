// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_REINDEX_LINEARPROBLEM_DECL_HPP
#define TPETRA_REINDEX_LINEARPROBLEM_DECL_HPP

/// \file Tpetra_Reindex_LinearProblem_decl.hpp
/// \brief Declaration of the Tpetra::Reindex_LinearProblem class

#include <Tpetra_Transform.hpp>
#include <Tpetra_LinearProblem.hpp>
#include <Tpetra_Reindex_CrsMatrix.hpp>
#include <Tpetra_Reindex_MultiVector.hpp>

namespace Tpetra {

///
/** Given and input Tpetra LinearProblem, a "reindexed" version will be returned
 *  using the given NewRowMap. If a null map is given, a lexigraphically indexed
 *  LP will be returned. The data in the new T_LP is a "reindexed" view of the
 *  original.
 */
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class Reindex_LinearProblem : public ViewTransform<LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node> > {
 public:
  using NewType      = typename ViewTransform<LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::NewType;
  using OriginalType = typename ViewTransform<LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::OriginalType;

  ///
  /** Constructor
   */
  Reindex_LinearProblem(Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> const> newRowMap);

  ///
  /** Destructor
   */
  ~Reindex_LinearProblem();

  ///
  /** Constructs a new view the original LP, "reindexed" using the given NewRowMap.
   */
  NewType operator()(OriginalType const& origProblem);

 private:
  Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> const> newRowMap_;

  Teuchos::RCP<Reindex_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > matTrans_;
  Teuchos::RCP<Reindex_MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > lhsTrans_;
  Teuchos::RCP<Reindex_MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rhsTrans_;
};

}  // namespace Tpetra

#endif  // TPETRA_REINDEX_LINEARPROBLEM_DECL_HPP
