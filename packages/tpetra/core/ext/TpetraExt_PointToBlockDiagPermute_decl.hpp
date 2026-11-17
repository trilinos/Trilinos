// @HEADER
// *****************************************************************************
//      TpetraExt: Tpetra Extended - Linear Algebra Services Package
//
// Copyright 2025 NTESS
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRAEXT_POINTTOBLOCKDIAGPERMUTE_DECL_HPP
#define TPETRAEXT_POINTTOBLOCKDIAGPERMUTE_DECL_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"

namespace Tpetra::Ext {

/** \brief Tpetra analogue of EpetraExt_PointToBlockDiagPermute
 *
 * This class extracts a point-to-block-diagonal approximation from a Tpetra::CrsMatrix.
 *
 * Supported modes:
 *   - contiguous block size
 *   - explicit local block starts + local/global block entry ids
 *
 * Blocks are purely local: each rank only keeps the entries of a block that it
 * owns. compute() inverts every block and stores the result as a
 * Tpetra::CrsMatrix (see getInverseBlockDiagonal()), so applying the inverse is
 * a sparse matrix-vector product on the matrix's execution space.
 */
template <class Scalar        = ::Tpetra::Details::DefaultTypes::scalar_type,
          class LocalOrdinal  = ::Tpetra::Details::DefaultTypes::local_ordinal_type,
          class GlobalOrdinal = ::Tpetra::Details::DefaultTypes::global_ordinal_type,
          class Node          = ::Tpetra::Details::DefaultTypes::node_type>
class PointToBlockDiagPermute {
 public:
  using map_type    = Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
  using crs_type    = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using mv_type     = Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using vector_type = Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  explicit PointToBlockDiagPermute(const crs_type& A);

  virtual ~PointToBlockDiagPermute() = default;

  int setParameters(Teuchos::ParameterList& list);
  int compute();

  /** \brief Apply the inverse block diagonal to X, returning Y.
   *
   * This applies the inverse of the extracted local block diagonal, not the
   * original matrix. Equivalent to <tt>getInverseBlockDiagonal()->apply(X, Y)</tt>.
   */
  int applyInverse(const mv_type& X, mv_type& Y) const;

  /** \brief Return the inverse of the extracted block diagonal as a sparse matrix.
   *
   * The matrix is computed by compute() and shared, not copied.
   */
  Teuchos::RCP<crs_type> getInverseBlockDiagonal() const { return invBlockDiagMatrix_; }

#if defined(TPETRA_ENABLE_DEPRECATED_CODE)
  /** \brief Return the extracted (non-inverted) block diagonal as a sparse matrix.
   *
   * \deprecated This is the Trilinos 17.x behavior. For the inverse, which is
   * what a block-diagonal preconditioner needs, use getInverseBlockDiagonal().
   */
  TPETRA_DEPRECATED_MSG("Use getInverseBlockDiagonal() for the inverse block diagonal")
  Teuchos::RCP<crs_type> createCrsMatrix() const;
#endif

 private:
  const crs_type* matrix_;
  Teuchos::ParameterList list_;

  bool purelyLocalMode_;
  bool contiguousBlockMode_;
  int contiguousBlockSize_;

  int numBlocks_;
  Teuchos::Array<int> blockStarts_;
  Teuchos::Array<GlobalOrdinal> blockGids_;

  Teuchos::RCP<const map_type> compatibleMap_;
  Teuchos::RCP<crs_type> invBlockDiagMatrix_;

  // Local row -> local block number (-1 if the row is in no block).
  std::vector<int> localToBlock_;

  void cleanupBlockInfo();
  int setupContiguousMode();
  int extractBlockDiagonal();
};

}  // namespace Tpetra::Ext

#endif
