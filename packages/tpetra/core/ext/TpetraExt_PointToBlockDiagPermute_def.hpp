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
  , purelyLocalMode_(true)
  , contiguousBlockMode_(false)
  , contiguousBlockSize_(0)
  , numBlocks_(0) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PointToBlockDiagPermute<Scalar, LocalOrdinal, GlobalOrdinal, Node>::cleanupBlockInfo() {
  blockStarts_.clear();
  blockGids_.clear();
  numBlocks_ = 0;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int PointToBlockDiagPermute<Scalar, LocalOrdinal, GlobalOrdinal, Node>::setParameters(
    Teuchos::ParameterList& list) {
  cleanupBlockInfo();
  list_ = list;

  contiguousBlockSize_ = list_.get("contiguous block size", 0);
  contiguousBlockMode_ = (contiguousBlockSize_ != 0);
  purelyLocalMode_     = true;

  if (contiguousBlockMode_) {
    return setupContiguousMode();
  }

  numBlocks_ = list_.get("number of local blocks", 0);

  TEUCHOS_TEST_FOR_EXCEPTION(!list_.isParameter("block start index"), std::runtime_error,
                             "PointToBlockDiagPermute: missing block start index");
  TEUCHOS_TEST_FOR_EXCEPTION(!list_.isParameter("block entry gids"), std::runtime_error,
                             "PointToBlockDiagPermute: missing block entry gids");
  blockStarts_ = list_.get<Teuchos::Array<int>>("block start index");
  blockGids_   = list_.get<Teuchos::Array<GlobalOrdinal>>("block entry gids");

  TEUCHOS_TEST_FOR_EXCEPTION(numBlocks_ <= 0, std::runtime_error,
                             "PointToBlockDiagPermute: invalid number of local blocks");

  return 0;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int PointToBlockDiagPermute<Scalar, LocalOrdinal, GlobalOrdinal, Node>::setupContiguousMode() {
  if (!contiguousBlockMode_) return 0;

  const auto rowMap            = matrix_->getRowMap();
  const GlobalOrdinal minMyGID = rowMap->getMinGlobalIndex();
  const GlobalOrdinal maxMyGID = rowMap->getMaxGlobalIndex();
  const GlobalOrdinal base     = rowMap->getIndexBase();

  const GlobalOrdinal myFirstBlockGID =
      static_cast<GlobalOrdinal>(contiguousBlockSize_ *
                                     std::ceil(static_cast<double>(minMyGID - base) /
                                               static_cast<double>(contiguousBlockSize_)) +
                                 base);

  numBlocks_ = static_cast<int>(
      std::ceil(static_cast<double>(maxMyGID - myFirstBlockGID + 1.0) /
                static_cast<double>(contiguousBlockSize_)));

  blockStarts_.resize(numBlocks_ + 1);
  blockGids_.resize(numBlocks_ * contiguousBlockSize_);

  blockStarts_[numBlocks_] = numBlocks_ * contiguousBlockSize_;

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

  const auto rowMap = matrix_->getRowMap();
  const auto colMap = matrix_->getColMap();

  const LocalOrdinal numMyRows = static_cast<LocalOrdinal>(rowMap->getLocalNumElements());

  std::vector<int> localToBlock(numMyRows, -1);
  std::vector<int> blockOffset(numMyRows, -1);

  for (int b = 0; b < numBlocks_; ++b) {
    for (int j = blockStarts_[b]; j < blockStarts_[b + 1]; ++j) {
      GlobalOrdinal gid = blockGids_[j];
      LocalOrdinal lid  = rowMap->getLocalElement(gid);
      if (lid != Teuchos::OrdinalTraits<LocalOrdinal>::invalid()) {
        localToBlock[lid] = b;
        blockOffset[lid]  = j - blockStarts_[b];
      }
    }
  }

  auto blockDiag =
      Teuchos::rcp(new crs_type(rowMap, contiguousBlockSize_ > 0 ? contiguousBlockSize_ : 8));

  typename crs_type::nonconst_local_inds_host_view_type inds(
      "inds", matrix_->getGlobalMaxNumRowEntries());
  typename crs_type::nonconst_values_host_view_type vals(
      "vals", matrix_->getGlobalMaxNumRowEntries());

  for (LocalOrdinal lrow = 0; lrow < numMyRows; ++lrow) {
    int blockNum = localToBlock[lrow];
    if (blockNum < 0) continue;

    GlobalOrdinal rowGid = rowMap->getGlobalElement(lrow);

    size_t numEntries = Teuchos::OrdinalTraits<size_t>::invalid();
    matrix_->getLocalRowCopy(lrow, inds, vals, numEntries);

    std::vector<GlobalOrdinal> outCols;
    std::vector<Scalar> outVals;

    for (size_t k = 0; k < numEntries; ++k) {
      LocalOrdinal lcol = inds(k);
      if (lcol == Teuchos::OrdinalTraits<LocalOrdinal>::invalid()) continue;
      if (lcol >= static_cast<LocalOrdinal>(localToBlock.size())) continue;
      if (localToBlock[lcol] != blockNum) continue;

      GlobalOrdinal colGid = colMap->getGlobalElement(lcol);
      outCols.push_back(colGid);
      outVals.push_back(vals(k));
    }

    if (outCols.empty()) {
      outCols.push_back(rowGid);
      outVals.push_back(Teuchos::ScalarTraits<Scalar>::one());
    }

    blockDiag->insertGlobalValues(rowGid, Teuchos::ArrayView<const GlobalOrdinal>(outCols),
                                  Teuchos::ArrayView<const Scalar>(outVals));
  }

  blockDiag->fillComplete(matrix_->getDomainMap(), matrix_->getRangeMap());
  compatibleMap_   = rowMap;
  blockDiagMatrix_ = blockDiag;

  return 0;
}

}  // namespace Tpetra::Ext

#endif
