// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_SOLVERMAP_CRSMATRIX_DECL_HPP
#define TPETRA_SOLVERMAP_CRSMATRIX_DECL_HPP

/// \file Tpetra_SolverMap_CrsMatrix_decl.hpp
/// \brief Declaration of the Tpetra::SolverMap_CrsMatrix class

#include <Tpetra_Transform.hpp>
#include <Tpetra_CrsMatrix.hpp>

namespace Tpetra {

///
/** Given an input CrsMatrix, the column map is checked for missing indices
 *  associated with the local rows. If found, a new CrsMatrix is formed using
 *  a new CrsGraph with a fixed column mapping including all local row indices.
 */
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class SolverMap_CrsMatrix : public StructuralSameTypeTransform<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > {
 public:
  using NewType      = typename StructuralSameTypeTransform<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::NewType;
  using OriginalType = typename StructuralSameTypeTransform<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::OriginalType;

  ///
  /** Constructor
   */
  SolverMap_CrsMatrix();

  ///
  /** Destructor
   */
  ~SolverMap_CrsMatrix();

  ///
  /** Constructs fixed view of CrsMatrix as necessary.
   */
  NewType operator()(OriginalType const& origMatrix);

 private:
  Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> > newColMap_;
  Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > newGraph_;
};

}  // namespace Tpetra

#endif  // TPETRA_SOLVERMAP_CRSMATRIX_DECL_HPP
