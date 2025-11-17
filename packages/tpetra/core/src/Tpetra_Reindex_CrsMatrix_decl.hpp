// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_REINDEX_CRSMATRIX_DECL_HPP
#define TPETRA_REINDEX_CRSMATRIX_DECL_HPP

/// \file Tpetra_Reindex_CrsMatrix_decl.hpp
/// \brief Declaration of the Tpetra::Reindex_CrsMatrix class

#include <Tpetra_Transform.hpp>
#include <Tpetra_CrsMatrix.hpp>

namespace Tpetra {

///
/** Given a Tpetra CrsMatrix, a "reindexed" version is returned based on
 *  the new row map. The row map must be conformal to the original. The
 *  Matrix data will be shared by the new Matrix using the new indexing
 */
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class Reindex_CrsMatrix : public ViewTransform<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > {
 public:
  using NewType      = typename ViewTransform<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::NewType;
  using OriginalType = typename ViewTransform<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::OriginalType;

  ///
  /** Constructor
   */
  Reindex_CrsMatrix(Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> const> newRowMap);

  ///
  /** Destructor
   */
  ~Reindex_CrsMatrix();

  ///
  /** Constructs "reindexed" Matrix
   */
  NewType operator()(OriginalType const& origMatrix);

 private:
  Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> const> newRowMap_;
  Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> > newColMap_;
};

}  // namespace Tpetra

#endif  // TPETRA_REINDEX_CRSMATRIX_DECL_HPP
