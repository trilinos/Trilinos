// @HEADER
// *****************************************************************************
//      TpetraExt: Tpetra Extended - Linear Algebra Services Package
//
// Copyright 2025 NTESS
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRAEXT_POINTTOBLOCKDIAGPERMUTE_DEF_HPP
#define TPETRAEXT_POINTTOBLOCKDIAGPERMUTE_DEF_HPP

#include "TpetraExt_PointToBlockDiagPermute_decl.hpp"

#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"

#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace Tpetra::Ext {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
PointToBlockDiagPermute<Scalar, LocalOrdinal, GlobalOrdinal, Node>::PointToBlockDiagPermute(
    const crs_type& A)
  : matrix_(&A)
  , list_()
  , purelyLocalMode_(true)
  , contiguousBlockMode_(false)
  , contiguousBlockSize_(0)
  , numBlocks_(0)
  , blockStarts_()
  , blockGids_()
  , compatibleMap_(Teuchos::null)
  , invBlockDiagMatrix_(Teuchos::null)
  , localToBlock_() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PointToBlockDiagPermute<Scalar, LocalOrdinal, GlobalOrdinal, Node>::cleanupBlockInfo() {
  blockStarts_.clear();
  blockGids_.clear();
  numBlocks_          = 0;
  compatibleMap_      = Teuchos::null;
  invBlockDiagMatrix_ = Teuchos::null;
  localToBlock_.clear();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int PointToBlockDiagPermute<Scalar, LocalOrdinal, GlobalOrdinal, Node>::setParameters(
    Teuchos::ParameterList& list) {
  cleanupBlockInfo();
  list_ = list;

  contiguousBlockSize_ = list_.get("contiguous block size", 0);
  TEUCHOS_TEST_FOR_EXCEPTION(contiguousBlockSize_ < 0, std::runtime_error,
                             "PointToBlockDiagPermute: contiguous block size must be non-negative");
  contiguousBlockMode_ = (contiguousBlockSize_ != 0);
  purelyLocalMode_     = true;

  if (contiguousBlockMode_) {
    return setupContiguousMode();
  }

  numBlocks_ = list_.get("number of local blocks", 0);
  TEUCHOS_TEST_FOR_EXCEPTION(numBlocks_ < 0, std::runtime_error,
                             "PointToBlockDiagPermute: invalid number of local blocks");

  if (numBlocks_ == 0) {
    blockStarts_.assign(1, 0);
    blockGids_.clear();
    return 0;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(!list_.isParameter("block start index"), std::runtime_error,
                             "PointToBlockDiagPermute: missing block start index");
  TEUCHOS_TEST_FOR_EXCEPTION(!list_.isParameter("block entry gids"), std::runtime_error,
                             "PointToBlockDiagPermute: missing block entry gids");
  blockStarts_ = list_.get<Teuchos::Array<int> >("block start index");
  blockGids_   = list_.get<Teuchos::Array<GlobalOrdinal> >("block entry gids");

  TEUCHOS_TEST_FOR_EXCEPTION(blockStarts_.size() < numBlocks_ + 1, std::runtime_error,
                             "PointToBlockDiagPermute: block start index array size is too small");
  TEUCHOS_TEST_FOR_EXCEPTION(blockGids_.size() < blockStarts_[numBlocks_], std::runtime_error,
                             "PointToBlockDiagPermute: block entry gids array size is too small");

  return 0;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int PointToBlockDiagPermute<Scalar, LocalOrdinal, GlobalOrdinal, Node>::setupContiguousMode() {
  if (!contiguousBlockMode_) return 0;

  const auto rowMap = matrix_->getRowMap();
  if (rowMap->getLocalNumElements() == 0) {
    numBlocks_ = 0;
    blockStarts_.assign(1, 0);
    blockGids_.clear();
    return 0;
  }

  const GlobalOrdinal minMyGID = rowMap->getMinGlobalIndex();
  const GlobalOrdinal maxMyGID = rowMap->getMaxGlobalIndex();
  const GlobalOrdinal base     = rowMap->getIndexBase();

  const GlobalOrdinal myFirstBlockGID =
      static_cast<GlobalOrdinal>(contiguousBlockSize_ *
                                     std::floor(static_cast<double>(minMyGID - base) /
                                                static_cast<double>(contiguousBlockSize_)) +
                                 base);

  numBlocks_ = static_cast<int>(
      std::ceil(static_cast<double>(maxMyGID - myFirstBlockGID + 1.0) /
                static_cast<double>(contiguousBlockSize_)));

  blockStarts_.resize(numBlocks_ + 1);

  const size_t numBlockEntries =
      static_cast<size_t>(numBlocks_) * static_cast<size_t>(contiguousBlockSize_);
  blockGids_.resize(numBlockEntries);

  blockStarts_[numBlocks_] = static_cast<int>(numBlockEntries);

  for (int i = 0, ct = 0; i < numBlocks_; i++) {
    blockStarts_[i] = ct;
    for (int j = 0; j < contiguousBlockSize_; j++, ct++) {
      blockGids_[ct] = myFirstBlockGID + ct;
    }
  }

  return 0;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int PointToBlockDiagPermute<Scalar, LocalOrdinal, GlobalOrdinal, Node>::compute() {
  return extractBlockDiagonal();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int PointToBlockDiagPermute<Scalar, LocalOrdinal, GlobalOrdinal, Node>::extractBlockDiagonal() {
  TEUCHOS_TEST_FOR_EXCEPTION(matrix_ == nullptr, std::runtime_error,
                             "PointToBlockDiagPermute: null matrix");

  using dense_mat_type = Teuchos::SerialDenseMatrix<int, Scalar>;

  const auto rowMap = matrix_->getRowMap();
  const auto colMap = matrix_->getColMap();

  const LocalOrdinal numMyRows = static_cast<LocalOrdinal>(rowMap->getLocalNumElements());

  localToBlock_.assign(numMyRows, -1);
  std::vector<int> blockOffset(numMyRows, -1);
  std::vector<std::vector<LocalOrdinal> > blockRows(numBlocks_);
  std::vector<dense_mat_type> invBlocks(numBlocks_);

  // Purely local mode: each block only contains the rows owned by this rank.
  // Block entries that are not locally owned (e.g., a contiguous block that
  // straddles a rank boundary) are dropped, as is any row already assigned to
  // an earlier block.
  for (int b = 0; b < numBlocks_; ++b) {
    for (int j = blockStarts_[b]; j < blockStarts_[b + 1]; ++j) {
      GlobalOrdinal gid = blockGids_[j];
      LocalOrdinal lid  = rowMap->getLocalElement(gid);
      if (lid == Teuchos::OrdinalTraits<LocalOrdinal>::invalid()) continue;
      if (localToBlock_[lid] >= 0) continue;
      localToBlock_[lid] = b;
      blockOffset[lid]   = static_cast<int>(blockRows[b].size());
      blockRows[b].push_back(lid);
    }

    // shape() zero-initializes the block
    const int bsz = static_cast<int>(blockRows[b].size());
    invBlocks[b].shape(bsz, bsz);
  }

  typename crs_type::nonconst_local_inds_host_view_type inds(
      "inds", matrix_->getLocalMaxNumRowEntries());
  typename crs_type::nonconst_values_host_view_type vals(
      "vals", matrix_->getLocalMaxNumRowEntries());

  for (LocalOrdinal lrow = 0; lrow < numMyRows; ++lrow) {
    int blockNum = localToBlock_[lrow];
    if (blockNum < 0) continue;

    const int rowInBlock = blockOffset[lrow];

    size_t numEntries = Teuchos::OrdinalTraits<size_t>::invalid();
    matrix_->getLocalRowCopy(lrow, inds, vals, numEntries);

    for (size_t k = 0; k < numEntries; ++k) {
      LocalOrdinal lcol = inds(k);
      if (lcol == Teuchos::OrdinalTraits<LocalOrdinal>::invalid()) continue;

      GlobalOrdinal colGid = colMap->getGlobalElement(lcol);
      if (colGid == Teuchos::OrdinalTraits<GlobalOrdinal>::invalid()) continue;

      LocalOrdinal rowLid = rowMap->getLocalElement(colGid);
      if (rowLid == Teuchos::OrdinalTraits<LocalOrdinal>::invalid()) continue;
      if (localToBlock_[rowLid] != blockNum) continue;

      const int colInBlock                        = blockOffset[rowLid];
      invBlocks[blockNum](rowInBlock, colInBlock) = vals(k);
    }
  }

  // Invert dense blocks in-place
  Teuchos::LAPACK<int, Scalar> lapack;
  for (int b = 0; b < numBlocks_; ++b) {
    const int n = invBlocks[b].numRows();
    if (n == 0) continue;
    std::vector<int> ipiv(n);
    int info = 0;

    lapack.GETRF(n, n, invBlocks[b].values(), invBlocks[b].stride(), ipiv.data(), &info);
    TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::runtime_error,
                               "PointToBlockDiagPermute: GETRF failed for block " << b
                                                                                  << " with info = " << info);

    std::vector<Scalar> work(std::max(1, n * n));
    lapack.GETRI(n, invBlocks[b].values(), invBlocks[b].stride(), ipiv.data(),
                 work.data(), static_cast<int>(work.size()), &info);
    TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::runtime_error,
                               "PointToBlockDiagPermute: GETRI failed for block " << b
                                                                                  << " with info = " << info);
  }

  // Assemble explicit inverse block-diagonal CRS matrix
  auto invBlockDiag =
      Teuchos::rcp(new crs_type(rowMap, contiguousBlockSize_ > 0 ? contiguousBlockSize_ : 8));

  std::vector<GlobalOrdinal> outCols;
  std::vector<Scalar> outVals;
  for (int b = 0; b < numBlocks_; ++b) {
    const int bsz = invBlocks[b].numRows();

    for (int i = 0; i < bsz; ++i) {
      const GlobalOrdinal rowGid = rowMap->getGlobalElement(blockRows[b][i]);

      outCols.clear();
      outVals.clear();
      for (int j = 0; j < bsz; ++j) {
        const Scalar val = invBlocks[b](i, j);
        if (val != Teuchos::ScalarTraits<Scalar>::zero()) {
          outCols.push_back(rowMap->getGlobalElement(blockRows[b][j]));
          outVals.push_back(val);
        }
      }

      invBlockDiag->insertGlobalValues(rowGid,
                                       Teuchos::ArrayView<const GlobalOrdinal>(outCols),
                                       Teuchos::ArrayView<const Scalar>(outVals));
    }
  }

  invBlockDiag->fillComplete(matrix_->getDomainMap(), matrix_->getRangeMap());
  compatibleMap_      = rowMap;
  invBlockDiagMatrix_ = invBlockDiag;

  return 0;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int PointToBlockDiagPermute<Scalar, LocalOrdinal, GlobalOrdinal, Node>::applyInverse(
    const mv_type& X, mv_type& Y) const {
  TEUCHOS_TEST_FOR_EXCEPTION(invBlockDiagMatrix_.is_null(), std::runtime_error,
                             "PointToBlockDiagPermute::applyInverse: compute() must be called first");
  invBlockDiagMatrix_->apply(X, Y);
  return 0;
}

#if defined(TPETRA_ENABLE_DEPRECATED_CODE)
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename PointToBlockDiagPermute<Scalar, LocalOrdinal, GlobalOrdinal, Node>::crs_type>
PointToBlockDiagPermute<Scalar, LocalOrdinal, GlobalOrdinal, Node>::createCrsMatrix() const {
  if (invBlockDiagMatrix_.is_null()) return Teuchos::null;

  const auto rowMap = matrix_->getRowMap();
  const auto colMap = matrix_->getColMap();

  const LocalOrdinal numMyRows = static_cast<LocalOrdinal>(rowMap->getLocalNumElements());

  auto blockDiag =
      Teuchos::rcp(new crs_type(rowMap, contiguousBlockSize_ > 0 ? contiguousBlockSize_ : 8));

  typename crs_type::nonconst_local_inds_host_view_type inds(
      "inds", matrix_->getLocalMaxNumRowEntries());
  typename crs_type::nonconst_values_host_view_type vals(
      "vals", matrix_->getLocalMaxNumRowEntries());

  std::vector<GlobalOrdinal> outCols;
  std::vector<Scalar> outVals;
  for (LocalOrdinal lrow = 0; lrow < numMyRows; ++lrow) {
    const int blockNum = localToBlock_[lrow];
    if (blockNum < 0) continue;

    size_t numEntries = Teuchos::OrdinalTraits<size_t>::invalid();
    matrix_->getLocalRowCopy(lrow, inds, vals, numEntries);

    outCols.clear();
    outVals.clear();
    for (size_t k = 0; k < numEntries; ++k) {
      const GlobalOrdinal colGid = colMap->getGlobalElement(inds(k));
      const LocalOrdinal rowLid  = rowMap->getLocalElement(colGid);
      if (rowLid == Teuchos::OrdinalTraits<LocalOrdinal>::invalid()) continue;
      if (localToBlock_[rowLid] != blockNum) continue;

      outCols.push_back(colGid);
      outVals.push_back(vals(k));
    }

    if (outCols.empty()) {
      outCols.push_back(rowMap->getGlobalElement(lrow));
      outVals.push_back(Teuchos::ScalarTraits<Scalar>::one());
    }

    blockDiag->insertGlobalValues(rowMap->getGlobalElement(lrow),
                                  Teuchos::ArrayView<const GlobalOrdinal>(outCols),
                                  Teuchos::ArrayView<const Scalar>(outVals));
  }

  blockDiag->fillComplete(matrix_->getDomainMap(), matrix_->getRangeMap());
  return blockDiag;
}
#endif

}  // namespace Tpetra::Ext

#endif
