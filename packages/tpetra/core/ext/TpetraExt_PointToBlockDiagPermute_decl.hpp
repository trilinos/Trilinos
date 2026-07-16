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
#include "Teuchos_SerialDenseMatrix.hpp"
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
 * The extracted block diagonal may be materialized as a Tpetra::CrsMatrix.
 * The materialized matrix represents the explicit inverse of the extracted
 * block diagonal, matching the semantics expected by Teko's BlkDiag path.
 */
template <class Scalar        = ::Tpetra::Details::DefaultTypes::scalar_type,
          class LocalOrdinal  = ::Tpetra::Details::DefaultTypes::local_ordinal_type,
          class GlobalOrdinal = ::Tpetra::Details::DefaultTypes::global_ordinal_type,
          class Node          = ::Tpetra::Details::DefaultTypes::node_type>
class PointToBlockDiagPermute {
 public:
  using map_type       = Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
  using crs_type       = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using mv_type        = Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using vector_type    = Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using dense_mat_type = Teuchos::SerialDenseMatrix<int, Scalar>;

  explicit PointToBlockDiagPermute(const crs_type& A);

  virtual ~PointToBlockDiagPermute() = default;

  int setParameters(Teuchos::ParameterList& list);
  int compute();

  /** \brief Apply the inverse block diagonal to X, returning Y.
   *
   * This applies the inverse of the extracted local block diagonal, not the
   * original matrix.
   */
  int applyInverse(const mv_type& X, mv_type& Y) const;

  /** \brief Return an explicit sparse matrix representation of the inverse
   * block diagonal.
   */
  Teuchos::RCP<crs_type> createCrsMatrix() const { return blockDiagMatrix_; }

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
  Teuchos::RCP<crs_type> blockDiagMatrix_;

  // Inverse dense blocks and mapping from block row -> local matrix row.
  std::vector<dense_mat_type> invBlocks_;
  std::vector<std::vector<LocalOrdinal> > blockRows_;

  void cleanupBlockInfo();
  int setupContiguousMode();
  int extractBlockDiagonal();
};

}  // namespace Tpetra::Ext

#endif
