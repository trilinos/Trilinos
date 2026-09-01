// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_BLOCKEDCRSMATRIX_DEF_HPP
#define TPETRA_BLOCKEDCRSMATRIX_DEF_HPP

#include <algorithm>
#include <sstream>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_OrdinalTraits.hpp>

#include "Tpetra_BlockedCrsMatrix_decl.hpp"

#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include "Tpetra_BlockedMap_def.hpp"
#include "Tpetra_MapExtractor_def.hpp"
#include "Tpetra_BlockedMultiVector_def.hpp"
#include "Tpetra_BlockedVector_def.hpp"

namespace Tpetra {

// ---------------------------------------------------------------------------
// Constructors / destructor
// ---------------------------------------------------------------------------

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BlockedCrsMatrix(const Teuchos::RCP<const BlockedMap>& rangeMaps,
                                                                              const Teuchos::RCP<const BlockedMap>& domainMaps,
                                                                              size_t numEntriesPerRow)
  : is_diagonal_(true) {
  domainmaps_       = Teuchos::rcp(new MapExtractor(domainMaps));
  rangemaps_        = Teuchos::rcp(new MapExtractor(rangeMaps));
  bRangeThyraMode_  = rangeMaps->getThyraMode();
  bDomainThyraMode_ = domainMaps->getThyraMode();

  blocks_.reserve(Rows() * Cols());

  // add CrsMatrix objects in row,column order
  for (size_t r = 0; r < Rows(); ++r)
    for (size_t c = 0; c < Cols(); ++c)
      blocks_.push_back(Teuchos::rcp(new CrsMatrix(getRangeMap(r, bRangeThyraMode_), numEntriesPerRow)));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BlockedCrsMatrix(Teuchos::RCP<const MapExtractor>& rangeMapExtractor,
                                                                              Teuchos::RCP<const MapExtractor>& domainMapExtractor,
                                                                              size_t numEntriesPerRow)
  : is_diagonal_(true)
  , domainmaps_(domainMapExtractor)
  , rangemaps_(rangeMapExtractor) {
  bRangeThyraMode_  = rangeMapExtractor->getThyraMode();
  bDomainThyraMode_ = domainMapExtractor->getThyraMode();

  blocks_.reserve(Rows() * Cols());

  // add CrsMatrix objects in row,column order
  for (size_t r = 0; r < Rows(); ++r)
    for (size_t c = 0; c < Cols(); ++c)
      blocks_.push_back(Teuchos::rcp(new CrsMatrix(getRangeMap(r, bRangeThyraMode_), numEntriesPerRow)));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~BlockedCrsMatrix() = default;

// ---------------------------------------------------------------------------
// Insertion / removal
// ---------------------------------------------------------------------------

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::insertGlobalValues(GlobalOrdinal globalRow, const Teuchos::ArrayView<const GlobalOrdinal>& cols, const Teuchos::ArrayView<const Scalar>& vals) {
  if (Rows() == 1 && Cols() == 1) {
    Teuchos::RCP<RowMatrix> block          = getMatrix(0, 0);
    Teuchos::RCP<BlockedCrsMatrix> bblock  = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(block);
    if (!bblock.is_null()) {
      bblock->insertGlobalValues(globalRow, cols, vals);
      return;
    }
    Teuchos::RCP<CrsMatrix> cblock = Teuchos::rcp_dynamic_cast<CrsMatrix>(block, true);
    cblock->insertGlobalValues(globalRow, cols, vals);
    return;
  }
  throw std::runtime_error("insertGlobalValues not supported by BlockedCrsMatrix");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::insertLocalValues(LocalOrdinal localRow, const Teuchos::ArrayView<const LocalOrdinal>& cols, const Teuchos::ArrayView<const Scalar>& vals) {
  if (Rows() == 1 && Cols() == 1) {
    Teuchos::RCP<RowMatrix> block         = getMatrix(0, 0);
    Teuchos::RCP<BlockedCrsMatrix> bblock = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(block);
    if (!bblock.is_null()) {
      bblock->insertLocalValues(localRow, cols, vals);
      return;
    }
    Teuchos::RCP<CrsMatrix> cblock = Teuchos::rcp_dynamic_cast<CrsMatrix>(block, true);
    cblock->insertLocalValues(localRow, cols, vals);
    return;
  }
  throw std::runtime_error("insertLocalValues not supported by BlockedCrsMatrix");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::removeEmptyProcessesInPlace(const Teuchos::RCP<const Map>& newMap) {
  if (Rows() == 1 && Cols() == 1) {
    Teuchos::RCP<CrsMatrix> cblock = Teuchos::rcp_dynamic_cast<CrsMatrix>(getMatrix(0, 0), true);
    cblock->removeEmptyProcessesInPlace(newMap);
    return;
  }
  throw std::runtime_error("removeEmptyProcesses not supported by BlockedCrsMatrix");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceGlobalValues(GlobalOrdinal globalRow,
                                                                                      const Teuchos::ArrayView<const GlobalOrdinal>& cols,
                                                                                      const Teuchos::ArrayView<const Scalar>& vals) {
  if (Rows() == 1 && Cols() == 1) {
    Teuchos::RCP<RowMatrix> block         = getMatrix(0, 0);
    Teuchos::RCP<BlockedCrsMatrix> bblock = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(block);
    if (!bblock.is_null()) {
      bblock->replaceGlobalValues(globalRow, cols, vals);
      return;
    }
    Teuchos::RCP<CrsMatrix> cblock = Teuchos::rcp_dynamic_cast<CrsMatrix>(block, true);
    cblock->replaceGlobalValues(globalRow, cols, vals);
    return;
  }
  throw std::runtime_error("replaceGlobalValues not supported by BlockedCrsMatrix");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceLocalValues(LocalOrdinal localRow,
                                                                                     const Teuchos::ArrayView<const LocalOrdinal>& cols,
                                                                                     const Teuchos::ArrayView<const Scalar>& vals) {
  if (Rows() == 1 && Cols() == 1) {
    Teuchos::RCP<RowMatrix> block         = getMatrix(0, 0);
    Teuchos::RCP<BlockedCrsMatrix> bblock = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(block);
    if (!bblock.is_null()) {
      bblock->replaceLocalValues(localRow, cols, vals);
      return;
    }
    Teuchos::RCP<CrsMatrix> cblock = Teuchos::rcp_dynamic_cast<CrsMatrix>(block, true);
    cblock->replaceLocalValues(localRow, cols, vals);
    return;
  }
  throw std::runtime_error("replaceLocalValues not supported by BlockedCrsMatrix");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::setAllToScalar(const Scalar& alpha) {
  for (size_t row = 0; row < Rows(); row++) {
    for (size_t col = 0; col < Cols(); col++) {
      Teuchos::RCP<RowMatrix> block = getMatrix(row, col);
      if (block.is_null()) continue;
      Teuchos::RCP<BlockedCrsMatrix> bblock = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(block);
      if (!bblock.is_null()) {
        bblock->setAllToScalar(alpha);
      } else {
        Teuchos::RCP<CrsMatrix> cblock = Teuchos::rcp_dynamic_cast<CrsMatrix>(block, true);
        cblock->setAllToScalar(alpha);
      }
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::scale(const Scalar& alpha) {
  for (size_t row = 0; row < Rows(); row++) {
    for (size_t col = 0; col < Cols(); col++) {
      Teuchos::RCP<RowMatrix> block = getMatrix(row, col);
      if (block.is_null()) continue;
      Teuchos::RCP<BlockedCrsMatrix> bblock = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(block);
      if (!bblock.is_null()) {
        bblock->scale(alpha);
      } else {
        Teuchos::RCP<CrsMatrix> cblock = Teuchos::rcp_dynamic_cast<CrsMatrix>(block, true);
        cblock->scale(alpha);
      }
    }
  }
}

// ---------------------------------------------------------------------------
// Transformational
// ---------------------------------------------------------------------------

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::resumeFill(const Teuchos::RCP<Teuchos::ParameterList>& params) {
  for (size_t row = 0; row < Rows(); row++) {
    for (size_t col = 0; col < Cols(); col++) {
      Teuchos::RCP<RowMatrix> block = getMatrix(row, col);
      if (block.is_null()) continue;
      Teuchos::RCP<BlockedCrsMatrix> bblock = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(block);
      if (!bblock.is_null()) {
        bblock->resumeFill(params);
      } else {
        Teuchos::RCP<CrsMatrix> cblock = Teuchos::rcp_dynamic_cast<CrsMatrix>(block, true);
        cblock->resumeFill(params);
      }
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::fillComplete(const Teuchos::RCP<const Map>& domainMap, const Teuchos::RCP<const Map>& rangeMap, const Teuchos::RCP<Teuchos::ParameterList>& params) {
  if (Rows() == 1 && Cols() == 1) {
    Teuchos::RCP<RowMatrix> block         = getMatrix(0, 0);
    Teuchos::RCP<BlockedCrsMatrix> bblock = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(block);
    if (!bblock.is_null()) {
      bblock->fillComplete(domainMap, rangeMap, params);
      return;
    }
    Teuchos::RCP<CrsMatrix> cblock = Teuchos::rcp_dynamic_cast<CrsMatrix>(block, true);
    cblock->fillComplete(domainMap, rangeMap, params);
    return;
  }
  fillComplete(params);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::fillComplete(const Teuchos::RCP<Teuchos::ParameterList>& params) {
  TEUCHOS_TEST_FOR_EXCEPTION(rangemaps_ == Teuchos::null, std::runtime_error, "BlockedCrsMatrix::fillComplete: rangemaps_ is not set. Error.");

  for (size_t r = 0; r < Rows(); ++r)
    for (size_t c = 0; c < Cols(); ++c) {
      Teuchos::RCP<RowMatrix> Ablock = getMatrix(r, c);
      if (Ablock.is_null()) continue;
      if (r != c) is_diagonal_ = false;
      if (!Ablock->isFillComplete()) {
        Teuchos::RCP<BlockedCrsMatrix> bAblock = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(Ablock);
        if (!bAblock.is_null()) {
          bAblock->fillComplete(getDomainMap(c, bDomainThyraMode_), getRangeMap(r, bRangeThyraMode_), params);
        } else {
          Teuchos::RCP<CrsMatrix> cAblock = Teuchos::rcp_dynamic_cast<CrsMatrix>(Ablock, true);
          cAblock->fillComplete(getDomainMap(c, bDomainThyraMode_), getRangeMap(r, bRangeThyraMode_), params);
        }
      }
    }
}

// ---------------------------------------------------------------------------
// RowMatrix query methods (aggregate over blocks)
// ---------------------------------------------------------------------------

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Teuchos::Comm<int>> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getComm() const {
  return rangemaps_->getFullMap()->getComm();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Map> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getRowMap() const {
  return rangemaps_->getFullMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Map> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getColMap() const {
  if (Rows() == 1 && Cols() == 1)
    return getMatrix(0, 0)->getColMap();
  return domainmaps_->getFullMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::RowGraph> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getGraph() const {
  if (Rows() == 1 && Cols() == 1)
    return getMatrix(0, 0)->getGraph();
  throw std::runtime_error("getGraph() not supported by BlockedCrsMatrix");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
global_size_t BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getGlobalNumRows() const {
  global_size_t globalNumRows = 0;
  for (size_t row = 0; row < Rows(); row++)
    for (size_t col = 0; col < Cols(); col++)
      if (!getMatrix(row, col).is_null()) {
        globalNumRows += getMatrix(row, col)->getGlobalNumRows();
        break;  // we need only one non-null matrix in a row
      }
  return globalNumRows;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
global_size_t BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getGlobalNumCols() const {
  global_size_t globalNumCols = 0;
  for (size_t col = 0; col < Cols(); col++)
    for (size_t row = 0; row < Rows(); row++)
      if (!getMatrix(row, col).is_null()) {
        globalNumCols += getMatrix(row, col)->getGlobalNumCols();
        break;  // we need only one non-null matrix in a col
      }
  return globalNumCols;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getLocalNumRows() const {
  size_t nodeNumRows = 0;
  for (size_t row = 0; row < Rows(); ++row)
    for (size_t col = 0; col < Cols(); col++)
      if (!getMatrix(row, col).is_null()) {
        nodeNumRows += getMatrix(row, col)->getLocalNumRows();
        break;  // we need only one non-null matrix in a row
      }
  return nodeNumRows;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getLocalNumCols() const {
  return getColMap()->getLocalNumElements();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getIndexBase() const {
  return getRowMap()->getIndexBase();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
global_size_t BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getGlobalNumEntries() const {
  global_size_t globalNumEntries = 0;
  for (size_t row = 0; row < Rows(); ++row)
    for (size_t col = 0; col < Cols(); ++col)
      if (!getMatrix(row, col).is_null())
        globalNumEntries += getMatrix(row, col)->getGlobalNumEntries();
  return globalNumEntries;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getLocalNumEntries() const {
  size_t nodeNumEntries = 0;
  for (size_t row = 0; row < Rows(); ++row)
    for (size_t col = 0; col < Cols(); ++col)
      if (!getMatrix(row, col).is_null())
        nodeNumEntries += getMatrix(row, col)->getLocalNumEntries();
  return nodeNumEntries;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNumEntriesInLocalRow(LocalOrdinal localRow) const {
  GlobalOrdinal gid           = this->getRowMap()->getGlobalElement(localRow);
  size_t row                  = getBlockedRangeMap()->getMapIndexForGID(gid);
  LocalOrdinal lid            = getBlockedRangeMap()->getMap(row)->getLocalElement(gid);
  size_t numEntriesInLocalRow = 0;
  for (size_t col = 0; col < Cols(); ++col)
    if (!getMatrix(row, col).is_null())
      numEntriesInLocalRow += getMatrix(row, col)->getNumEntriesInLocalRow(lid);
  return numEntriesInLocalRow;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const {
  size_t row                   = getBlockedRangeMap()->getMapIndexForGID(globalRow);
  size_t numEntriesInGlobalRow = 0;
  for (size_t col = 0; col < Cols(); ++col)
    if (!getMatrix(row, col).is_null())
      numEntriesInGlobalRow += getMatrix(row, col)->getNumEntriesInGlobalRow(globalRow);
  return numEntriesInGlobalRow;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getGlobalMaxNumRowEntries() const {
  global_size_t globalMaxEntries = 0;
  for (size_t row = 0; row < Rows(); row++) {
    global_size_t globalMaxEntriesBlockRows = 0;
    for (size_t col = 0; col < Cols(); col++) {
      if (!getMatrix(row, col).is_null()) {
        globalMaxEntriesBlockRows += getMatrix(row, col)->getGlobalMaxNumRowEntries();
      }
    }
    if (globalMaxEntriesBlockRows > globalMaxEntries)
      globalMaxEntries = globalMaxEntriesBlockRows;
  }
  return globalMaxEntries;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getBlockSize() const {
  return Teuchos::as<LocalOrdinal>(1);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getLocalMaxNumRowEntries() const {
  size_t localMaxEntries = 0;
  for (size_t row = 0; row < Rows(); row++) {
    size_t localMaxEntriesBlockRows = 0;
    for (size_t col = 0; col < Cols(); col++) {
      if (!getMatrix(row, col).is_null()) {
        localMaxEntriesBlockRows += getMatrix(row, col)->getLocalMaxNumRowEntries();
      }
    }
    if (localMaxEntriesBlockRows > localMaxEntries)
      localMaxEntries = localMaxEntriesBlockRows;
  }
  return localMaxEntries;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::hasColMap() const {
  if (Rows() == 1 && Cols() == 1)
    return getMatrix(0, 0)->hasColMap();
  return false;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isLocallyIndexed() const {
  for (size_t i = 0; i < blocks_.size(); ++i)
    if (blocks_[i] != Teuchos::null && !blocks_[i]->isLocallyIndexed())
      return false;
  return true;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isGloballyIndexed() const {
  for (size_t i = 0; i < blocks_.size(); i++)
    if (blocks_[i] != Teuchos::null && !blocks_[i]->isGloballyIndexed())
      return false;
  return true;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isFillComplete() const {
  for (size_t i = 0; i < blocks_.size(); i++)
    if (blocks_[i] != Teuchos::null && !blocks_[i]->isFillComplete())
      return false;
  return true;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::supportsRowViews() const {
  if (Rows() == 1 && Cols() == 1)
    return getMatrix(0, 0)->supportsRowViews();
  return false;
}

// ---------------------------------------------------------------------------
// RowMatrix extraction methods
// ---------------------------------------------------------------------------

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getGlobalRowCopy(GlobalOrdinal GlobalRow,
                                                                                   nonconst_global_inds_host_view_type& Indices,
                                                                                   nonconst_values_host_view_type& Values,
                                                                                   size_t& NumEntries) const {
  if (Rows() == 1 && Cols() == 1) {
    getMatrix(0, 0)->getGlobalRowCopy(GlobalRow, Indices, Values, NumEntries);
    return;
  }
  throw std::runtime_error("getGlobalRowCopy not supported by BlockedCrsMatrix");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getLocalRowCopy(LocalOrdinal LocalRow,
                                                                                  nonconst_local_inds_host_view_type& Indices,
                                                                                  nonconst_values_host_view_type& Values,
                                                                                  size_t& NumEntries) const {
  if (Rows() == 1 && Cols() == 1) {
    getMatrix(0, 0)->getLocalRowCopy(LocalRow, Indices, Values, NumEntries);
    return;
  }
  throw std::runtime_error("getLocalRowCopy not supported by BlockedCrsMatrix");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getGlobalRowView(GlobalOrdinal GlobalRow,
                                                                                   global_inds_host_view_type& indices,
                                                                                   values_host_view_type& values) const {
  if (Rows() == 1 && Cols() == 1) {
    getMatrix(0, 0)->getGlobalRowView(GlobalRow, indices, values);
    return;
  }
  throw std::runtime_error("getGlobalRowView not supported by BlockedCrsMatrix");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getLocalRowView(LocalOrdinal LocalRow,
                                                                                  local_inds_host_view_type& indices,
                                                                                  values_host_view_type& values) const {
  if (Rows() == 1 && Cols() == 1) {
    getMatrix(0, 0)->getLocalRowView(LocalRow, indices, values);
    return;
  } else if (is_diagonal_) {
    GlobalOrdinal gid = this->getRowMap()->getGlobalElement(LocalRow);
    size_t row        = getBlockedRangeMap()->getMapIndexForGID(gid);
    getMatrix(row, row)->getLocalRowView(getMatrix(row, row)->getRowMap()->getLocalElement(gid), indices, values);
    return;
  }
  throw std::runtime_error("getLocalRowView not supported by BlockedCrsMatrix");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getLocalDiagCopy(Vector& diag) const {
  Teuchos::RCP<Vector> rcpdiag      = Teuchos::rcpFromRef(diag);
  Teuchos::RCP<BlockedVector> bdiag = Teuchos::rcp_dynamic_cast<BlockedVector>(rcpdiag);

  // special treatment for 1x1 block matrices: leaf blocks store plain Vectors.
  if (Rows() == 1 && Cols() == 1 && bdiag.is_null() == true) {
    Teuchos::RCP<const RowMatrix> rm = getMatrix(0, 0);
    rm->getLocalDiagCopy(diag);
    return;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(bdiag.is_null() == true, std::runtime_error, "BlockedCrsMatrix::getLocalDiagCopy: diag must be a Blocked(Multi)Vector.");
  TEUCHOS_TEST_FOR_EXCEPTION(bdiag->getNumVectors() != 1, std::runtime_error, "BlockedCrsMatrix::getLocalDiagCopy: diag must be a Blocked(Multi)Vector with exactly one vector. However, the number of stored vectors is " << bdiag->getNumVectors());
  TEUCHOS_TEST_FOR_EXCEPTION(bdiag->getBlockedMap()->getNumMaps() != this->Rows(), std::runtime_error,
                             "BlockedCrsMatrix::getLocalDiagCopy(): the number of blocks in diag differ from the number of blocks in this operator.");

  for (size_t row = 0; row < Rows(); row++) {
    Teuchos::RCP<const RowMatrix> rm = getMatrix(row, row);
    if (!rm.is_null()) {
      Teuchos::RCP<Vector> rv = Teuchos::rcp(new Vector(bdiag->getBlockedMap()->getMap(row, bdiag->getBlockedMap()->getThyraMode())));
      rm->getLocalDiagCopy(*rv);
      bdiag->setMultiVector(row, rv, bdiag->getBlockedMap()->getThyraMode());
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::leftScale(const Vector& x) {
  Teuchos::RCP<const Vector> rcpx      = Teuchos::rcpFromRef(x);
  Teuchos::RCP<const BlockedVector> bx = Teuchos::rcp_dynamic_cast<const BlockedVector>(rcpx);

  // special treatment for 1xn block matrices: leaf blocks store plain Vectors.
  if (Rows() == 1 && bx.is_null() == true) {
    for (size_t col = 0; col < Cols(); ++col) {
      Teuchos::RCP<RowMatrix> rm = getMatrix(0, col);
      rm->leftScale(x);
    }
    return;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(bx.is_null() == true, std::runtime_error, "BlockedCrsMatrix::leftScale: x must be a Blocked(Multi)Vector.");
  TEUCHOS_TEST_FOR_EXCEPTION(bx->getNumVectors() != 1, std::runtime_error, "BlockedCrsMatrix::leftScale: x must be a Blocked(Multi)Vector with exactly one vector. However, the number of stored vectors is " << bx->getNumVectors());
  TEUCHOS_TEST_FOR_EXCEPTION(bx->getBlockedMap()->getNumMaps() != this->Rows(), std::runtime_error,
                             "BlockedCrsMatrix::leftScale(): the number of blocks in diag differ from the number of blocks in this operator.");

  for (size_t row = 0; row < Rows(); row++) {
    Teuchos::RCP<const MultiVector> rmv = bx->getMultiVector(row);
    Teuchos::RCP<const Vector> rscale   = rmv->getVector(0);
    TEUCHOS_TEST_FOR_EXCEPTION(rscale.is_null() == true, std::runtime_error, "BlockedCrsMatrix::leftScale: x must be a Vector.");
    for (size_t col = 0; col < Cols(); ++col) {
      Teuchos::RCP<RowMatrix> rm = getMatrix(row, col);
      if (!rm.is_null()) {
        rm->leftScale(*rscale);
      }
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::rightScale(const Vector& x) {
  Teuchos::RCP<const Vector> rcpx      = Teuchos::rcpFromRef(x);
  Teuchos::RCP<const BlockedVector> bx = Teuchos::rcp_dynamic_cast<const BlockedVector>(rcpx);

  // special treatment for nx1 block matrices: leaf blocks store plain Vectors.
  if (Cols() == 1 && bx.is_null() == true) {
    for (size_t row = 0; row < Rows(); ++row) {
      Teuchos::RCP<RowMatrix> rm = getMatrix(row, 0);
      rm->rightScale(x);
    }
    return;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(bx.is_null() == true, std::runtime_error, "BlockedCrsMatrix::rightScale: x must be a Blocked(Multi)Vector.");
  TEUCHOS_TEST_FOR_EXCEPTION(bx->getNumVectors() != 1, std::runtime_error, "BlockedCrsMatrix::rightScale: x must be a Blocked(Multi)Vector with exactly one vector. However, the number of stored vectors is " << bx->getNumVectors());
  TEUCHOS_TEST_FOR_EXCEPTION(bx->getBlockedMap()->getNumMaps() != this->Cols(), std::runtime_error,
                             "BlockedCrsMatrix::rightScale(): the number of blocks in diag differ from the number of blocks in this operator.");

  for (size_t col = 0; col < Cols(); ++col) {
    Teuchos::RCP<const MultiVector> rmv = bx->getMultiVector(col);
    Teuchos::RCP<const Vector> rscale   = rmv->getVector(0);
    TEUCHOS_TEST_FOR_EXCEPTION(rscale.is_null() == true, std::runtime_error, "BlockedCrsMatrix::rightScale: x must be a Vector.");
    for (size_t row = 0; row < Rows(); row++) {
      Teuchos::RCP<RowMatrix> rm = getMatrix(row, col);
      if (!rm.is_null()) {
        rm->rightScale(*rscale);
      }
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::mag_type BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getFrobeniusNorm() const {
  mag_type ret = Teuchos::ScalarTraits<mag_type>::zero();
  for (size_t col = 0; col < Cols(); ++col) {
    for (size_t row = 0; row < Rows(); ++row) {
      if (getMatrix(row, col) != Teuchos::null) {
        mag_type n = getMatrix(row, col)->getFrobeniusNorm();
        ret += n * n;
      }
    }
  }
  return Teuchos::ScalarTraits<mag_type>::squareroot(ret);
}

// ---------------------------------------------------------------------------
// Operator interface (apply)
// ---------------------------------------------------------------------------

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::apply(const MultiVector& /* X */, MultiVector& /* Y */, Teuchos::ETransp /* mode */, Scalar /* alpha */, Scalar /* beta */, bool /* sumInterfaceValues */,
                                                                        const Teuchos::RCP<Import>& /* regionInterfaceImporter */,
                                                                        const Teuchos::ArrayRCP<LocalOrdinal>& /* regionInterfaceLIDs */) const {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::apply(const MultiVector& X, MultiVector& Y,
                                                                        Teuchos::ETransp mode,
                                                                        Scalar alpha,
                                                                        Scalar beta) const {
  TEUCHOS_TEST_FOR_EXCEPTION(mode != Teuchos::NO_TRANS && mode != Teuchos::TRANS, std::runtime_error,
                             "apply() only supports the following modes: NO_TRANS and TRANS.");

  // check whether input parameters are blocked or not
  Teuchos::RCP<const MultiVector> refX         = Teuchos::rcpFromRef(X);
  Teuchos::RCP<const BlockedMultiVector> refbX = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(refX);

  bool bBlockedX = (refbX != Teuchos::null) ? true : false;

  // NOTE (Tpetra port): unlike Xpetra::MultiVector, Tpetra::MultiVector::getNumVectors(),
  // getMap() and update() are NOT virtual.  When Y is actually a BlockedMultiVector but is
  // referenced here as a plain MultiVector&, those accessors would bind to the (empty)
  // MultiVector base sub-object and report 0 columns / a null map.  Recover the blocked
  // type and query/operate through it so the output vectors are sized correctly.
  Teuchos::RCP<MultiVector> refY         = Teuchos::rcpFromRef(Y);
  Teuchos::RCP<BlockedMultiVector> refbY = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(refY);
  const bool bBlockedY                   = !refbY.is_null();
  const size_t numVecsY                  = bBlockedY ? refbY->getNumVectors() : Y.getNumVectors();

  // create (temporary) vectors for output; in the end we call Y.update(alpha, *tmpY, beta).
  Teuchos::RCP<MultiVector> tmpY;
  if (bBlockedY)
    tmpY = Teuchos::rcp(new BlockedMultiVector(refbY->getBlockedMap(), numVecsY, true));
  else
    tmpY = Teuchos::rcp(new MultiVector(Y.getMap(), numVecsY, true));

  Scalar one = Teuchos::ScalarTraits<Scalar>::one();

  if (mode == Teuchos::NO_TRANS) {
    for (size_t row = 0; row < Rows(); row++) {
      Teuchos::RCP<MultiVector> Yblock = rangemaps_->getVector(row, numVecsY, bRangeThyraMode_, true);
      for (size_t col = 0; col < Cols(); col++) {
        Teuchos::RCP<RowMatrix> Ablock = getMatrix(row, col);

        if (Ablock.is_null())
          continue;

        // check whether Ablock is itself a blocked operator
        bool bBlockedSubMatrix = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(Ablock) == Teuchos::null ? false : true;

        Teuchos::RCP<const MultiVector> Xblock = Teuchos::null;

        // extract sub part of X using Xpetra or Thyra GIDs
        if (bBlockedX)
          Xblock = domainmaps_->ExtractVector(refbX, col, bDomainThyraMode_);
        else
          Xblock = domainmaps_->ExtractVector(refX, col, bBlockedSubMatrix == true ? false : bDomainThyraMode_);

        Teuchos::RCP<MultiVector> tmpYblock = rangemaps_->getVector(row, numVecsY, bRangeThyraMode_, false);
        Ablock->apply(*Xblock, *tmpYblock);

        Yblock->update(one, *tmpYblock, one);
      }
      rangemaps_->InsertVector(Yblock, row, tmpY, bRangeThyraMode_);
    }

  } else if (mode == Teuchos::TRANS) {
    for (size_t col = 0; col < Cols(); col++) {
      Teuchos::RCP<MultiVector> Yblock = domainmaps_->getVector(col, numVecsY, bDomainThyraMode_, true);

      for (size_t row = 0; row < Rows(); row++) {
        Teuchos::RCP<RowMatrix> Ablock = getMatrix(row, col);

        if (Ablock.is_null())
          continue;

        bool bBlockedSubMatrix = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(Ablock) == Teuchos::null ? false : true;

        Teuchos::RCP<const MultiVector> Xblock = Teuchos::null;

        if (bBlockedX)
          Xblock = rangemaps_->ExtractVector(refbX, row, bRangeThyraMode_);
        else
          Xblock = rangemaps_->ExtractVector(refX, row, bBlockedSubMatrix == true ? false : bRangeThyraMode_);
        Teuchos::RCP<MultiVector> tmpYblock = domainmaps_->getVector(col, numVecsY, bDomainThyraMode_, false);
        Ablock->apply(*Xblock, *tmpYblock, Teuchos::TRANS);

        Yblock->update(one, *tmpYblock, one);
      }
      domainmaps_->InsertVector(Yblock, col, tmpY, bDomainThyraMode_);
    }
  }
  // NOTE (Tpetra port): dispatch update() through the blocked type when Y is blocked
  // (base MultiVector::update() is non-virtual and would no-op on the empty base).
  if (bBlockedY)
    refbY->update(alpha, *tmpY, beta);
  else
    Y.update(alpha, *tmpY, beta);
}

// ---------------------------------------------------------------------------
// Domain / range map accessors
// ---------------------------------------------------------------------------

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Map> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getDomainMap() const {
  return domainmaps_->getFullMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Map> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getRangeMap() const {
  return rangemaps_->getFullMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Map> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getFullDomainMap() const {
  return domainmaps_->getFullMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BlockedMap> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getBlockedDomainMap() const {
  return domainmaps_->getBlockedMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Map> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getDomainMap(size_t i) const {
  return domainmaps_->getMap(i, bDomainThyraMode_);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Map> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getDomainMap(size_t i, bool bThyraMode) const {
  return domainmaps_->getMap(i, bThyraMode);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Map> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getFullRangeMap() const {
  return rangemaps_->getFullMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BlockedMap> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getBlockedRangeMap() const {
  return rangemaps_->getBlockedMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Map> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getRangeMap(size_t i) const {
  return rangemaps_->getMap(i, bRangeThyraMode_);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Map> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getRangeMap(size_t i, bool bThyraMode) const {
  return rangemaps_->getMap(i, bThyraMode);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MapExtractor> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getRangeMapExtractor() const {
  return rangemaps_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MapExtractor> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getDomainMapExtractor() const {
  return domainmaps_;
}

// ---------------------------------------------------------------------------
// bgs_apply
// ---------------------------------------------------------------------------

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::bgs_apply(const MultiVector& X,
                                                                            MultiVector& Y,
                                                                            size_t row,
                                                                            Teuchos::ETransp mode,
                                                                            Scalar alpha,
                                                                            Scalar beta) const {
  TEUCHOS_TEST_FOR_EXCEPTION(mode != Teuchos::NO_TRANS && mode != Teuchos::TRANS, std::runtime_error,
                             "apply() only supports the following modes: NO_TRANS and TRANS.");

  Teuchos::RCP<const MultiVector> refX         = Teuchos::rcpFromRef(X);
  Teuchos::RCP<const BlockedMultiVector> refbX = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(refX);

  bool bBlockedX = (refbX != Teuchos::null) ? true : false;

  // NOTE (Tpetra port): see apply() above. getNumVectors()/getMap()/update() are non-virtual
  // on the MultiVector base, so recover the blocked type of Y and work through it.
  Teuchos::RCP<MultiVector> refY         = Teuchos::rcpFromRef(Y);
  Teuchos::RCP<BlockedMultiVector> refbY = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(refY);
  const bool bBlockedY                   = !refbY.is_null();
  const size_t numVecsY                  = bBlockedY ? refbY->getNumVectors() : Y.getNumVectors();

  Teuchos::RCP<MultiVector> tmpY;
  if (bBlockedY)
    tmpY = Teuchos::rcp(new BlockedMultiVector(refbY->getBlockedMap(), numVecsY, true));
  else
    tmpY = Teuchos::rcp(new MultiVector(Y.getMap(), numVecsY, true));

  Scalar one = Teuchos::ScalarTraits<Scalar>::one();

  if (mode == Teuchos::NO_TRANS) {
    Teuchos::RCP<MultiVector> Yblock = rangemaps_->getVector(row, numVecsY, bRangeThyraMode_, true);
    for (size_t col = 0; col < Cols(); col++) {
      Teuchos::RCP<RowMatrix> Ablock = getMatrix(row, col);

      if (Ablock.is_null())
        continue;

      bool bBlockedSubMatrix = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(Ablock) == Teuchos::null ? false : true;

      Teuchos::RCP<const MultiVector> Xblock = Teuchos::null;

      if (bBlockedX)
        Xblock = domainmaps_->ExtractVector(refbX, col, bDomainThyraMode_);
      else
        Xblock = domainmaps_->ExtractVector(refX, col, bBlockedSubMatrix == true ? false : bDomainThyraMode_);

      Teuchos::RCP<MultiVector> tmpYblock = rangemaps_->getVector(row, numVecsY, bRangeThyraMode_, false);
      Ablock->apply(*Xblock, *tmpYblock);

      Yblock->update(one, *tmpYblock, one);
    }
    rangemaps_->InsertVector(Yblock, row, tmpY, bRangeThyraMode_);
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Tpetra::BlockedCrsMatrix::bgs_apply: not implemented for transpose case.");
  }
  if (bBlockedY)
    refbY->update(alpha, *tmpY, beta);
  else
    Y.update(alpha, *tmpY, beta);
}

// ---------------------------------------------------------------------------
// DistObject-like interface (delegate-or-throw)
// ---------------------------------------------------------------------------

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doImport(const RowMatrix& source, const Import& importer, CombineMode CM) {
  if (Rows() == 1 && Cols() == 1) {
    Teuchos::RCP<CrsMatrix> cblock         = Teuchos::rcp_dynamic_cast<CrsMatrix>(getMatrix(0, 0), true);
    const CrsMatrix* csource               = dynamic_cast<const CrsMatrix*>(&source);
    TEUCHOS_TEST_FOR_EXCEPTION(csource == nullptr, std::runtime_error, "BlockedCrsMatrix::doImport(): source must be a Tpetra::CrsMatrix.");
    cblock->doImport(*csource, importer, CM);
    return;
  }
  throw std::runtime_error("BlockedCrsMatrix::doImport(): operation not supported.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doExport(const RowMatrix& dest, const Import& importer, CombineMode CM) {
  if (Rows() == 1 && Cols() == 1) {
    Teuchos::RCP<CrsMatrix> cblock = Teuchos::rcp_dynamic_cast<CrsMatrix>(getMatrix(0, 0), true);
    const CrsMatrix* cdest         = dynamic_cast<const CrsMatrix*>(&dest);
    TEUCHOS_TEST_FOR_EXCEPTION(cdest == nullptr, std::runtime_error, "BlockedCrsMatrix::doExport(): dest must be a Tpetra::CrsMatrix.");
    cblock->doExport(*cdest, importer, CM);
    return;
  }
  throw std::runtime_error("BlockedCrsMatrix::doExport(): operation not supported.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doImport(const RowMatrix& source, const Export& exporter, CombineMode CM) {
  if (Rows() == 1 && Cols() == 1) {
    Teuchos::RCP<CrsMatrix> cblock = Teuchos::rcp_dynamic_cast<CrsMatrix>(getMatrix(0, 0), true);
    const CrsMatrix* csource       = dynamic_cast<const CrsMatrix*>(&source);
    TEUCHOS_TEST_FOR_EXCEPTION(csource == nullptr, std::runtime_error, "BlockedCrsMatrix::doImport(): source must be a Tpetra::CrsMatrix.");
    cblock->doImport(*csource, exporter, CM);
    return;
  }
  throw std::runtime_error("BlockedCrsMatrix::doImport(): operation not supported.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doExport(const RowMatrix& dest, const Export& exporter, CombineMode CM) {
  if (Rows() == 1 && Cols() == 1) {
    Teuchos::RCP<CrsMatrix> cblock = Teuchos::rcp_dynamic_cast<CrsMatrix>(getMatrix(0, 0), true);
    const CrsMatrix* cdest         = dynamic_cast<const CrsMatrix*>(&dest);
    TEUCHOS_TEST_FOR_EXCEPTION(cdest == nullptr, std::runtime_error, "BlockedCrsMatrix::doExport(): dest must be a Tpetra::CrsMatrix.");
    cblock->doExport(*cdest, exporter, CM);
    return;
  }
  throw std::runtime_error("BlockedCrsMatrix::doExport(): operation not supported.");
}

// ---------------------------------------------------------------------------
// Describable / LabeledObject
// ---------------------------------------------------------------------------

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const { return "Tpetra::BlockedCrsMatrix"; }

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const {
  out << "Tpetra::BlockedCrsMatrix: " << Rows() << " x " << Cols() << std::endl;

  if (isFillComplete()) {
    out << "BlockMatrix is fillComplete" << std::endl;
  } else {
    out << "BlockMatrix is NOT fillComplete" << std::endl;
  }

  for (size_t r = 0; r < Rows(); ++r)
    for (size_t c = 0; c < Cols(); ++c) {
      if (getMatrix(r, c) != Teuchos::null) {
        out << "Block(" << r << "," << c << ")" << std::endl;
        getMatrix(r, c)->describe(out, verbLevel);
      } else
        out << "Block(" << r << "," << c << ") = null" << std::endl;
    }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::setObjectLabel(const std::string& objectLabel) {
  for (size_t r = 0; r < Rows(); ++r)
    for (size_t c = 0; c < Cols(); ++c) {
      if (getMatrix(r, c) != Teuchos::null) {
        std::ostringstream oss;
        oss << objectLabel << "(" << r << "," << c << ")";
        getMatrix(r, c)->setObjectLabel(oss.str());
      }
    }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::hasCrsGraph() const {
  return (Rows() == 1 && Cols() == 1);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CrsGraph> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getCrsGraph() const {
  if (Rows() == 1 && Cols() == 1) {
    Teuchos::RCP<CrsMatrix> cblock = Teuchos::rcp_dynamic_cast<CrsMatrix>(getMatrix(0, 0), true);
    return cblock->getCrsGraph();
  }
  throw std::runtime_error("getCrsGraph() not supported by BlockedCrsMatrix");
}

// ---------------------------------------------------------------------------
// Block matrix access
// ---------------------------------------------------------------------------

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isDiagonal() const { return is_diagonal_; }

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Rows() const {
  return rangemaps_->NumMaps();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Cols() const {
  return domainmaps_->NumMaps();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::RowMatrix> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getCrsMatrix() const {
  TEUCHOS_TEST_FOR_EXCEPTION(Rows() != 1, std::out_of_range, "Can only unwrap a 1x1 blocked matrix. The matrix has " << Rows() << " block rows, though.");
  TEUCHOS_TEST_FOR_EXCEPTION(Cols() != 1, std::out_of_range, "Can only unwrap a 1x1 blocked matrix. The matrix has " << Cols() << " block columns, though.");

  Teuchos::RCP<RowMatrix> mat            = getMatrix(0, 0);
  Teuchos::RCP<BlockedCrsMatrix> bmat    = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(mat);
  if (bmat == Teuchos::null) return mat;
  return bmat->getCrsMatrix();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::RowMatrix> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getInnermostCrsMatrix() {
  size_t row = Rows() + 1, col = Cols() + 1;
  for (size_t r = 0; r < Rows(); ++r)
    for (size_t c = 0; c < Cols(); ++c)
      if (getMatrix(r, c) != Teuchos::null) {
        row = r;
        col = c;
        break;
      }
  TEUCHOS_TEST_FOR_EXCEPTION(row == Rows() + 1 || col == Cols() + 1, std::runtime_error, "Tpetra::BlockedCrsMatrix::getInnermostCrsMatrix: Could not find a non-zero sub-block in blocked operator.");
  Teuchos::RCP<RowMatrix> mm          = getMatrix(row, col);
  Teuchos::RCP<BlockedCrsMatrix> bmat = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(mm);
  if (bmat == Teuchos::null) return mm;
  return bmat->getInnermostCrsMatrix();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::RowMatrix> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getMatrix(size_t r, size_t c) const {
  TEUCHOS_TEST_FOR_EXCEPTION(r > Rows(), std::out_of_range, "Error, r = " << Rows() << " is too big");
  TEUCHOS_TEST_FOR_EXCEPTION(c > Cols(), std::out_of_range, "Error, c = " << Cols() << " is too big");
  return blocks_[r * Cols() + c];
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::setMatrix(size_t r, size_t c, Teuchos::RCP<RowMatrix> mat) {
  TEUCHOS_TEST_FOR_EXCEPTION(r > Rows(), std::out_of_range, "Error, r = " << Rows() << " is too big");
  TEUCHOS_TEST_FOR_EXCEPTION(c > Cols(), std::out_of_range, "Error, c = " << Cols() << " is too big");
  if (!mat.is_null() && r != c) is_diagonal_ = false;
  blocks_[r * Cols() + c] = mat;
}

// ---------------------------------------------------------------------------
// Merge
// ---------------------------------------------------------------------------

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CrsMatrix> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Merge() const {
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  Scalar one = Teuchos::ScalarTraits<Scalar>::one();

  TEUCHOS_TEST_FOR_EXCEPTION(bRangeThyraMode_ != bDomainThyraMode_, std::runtime_error,
                             "BlockedCrsMatrix::Merge: only implemented for Xpetra-style or Thyra-style numbering. No mixup allowed!");

  TEUCHOS_TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
                             "BlockedCrsMatrix::Merge: BlockMatrix must be fill-completed.");

  LocalOrdinal lclNumRows = static_cast<LocalOrdinal>(getFullRangeMap()->getLocalNumElements());
  Teuchos::ArrayRCP<size_t> numEntPerRow(lclNumRows);
  for (LocalOrdinal lclRow = 0; lclRow < lclNumRows; ++lclRow)
    numEntPerRow[lclRow] = getNumEntriesInLocalRow(lclRow);

  RCP<CrsMatrix> sparse = Teuchos::rcp(new CrsMatrix(getFullRangeMap(), numEntPerRow()));

  if (bRangeThyraMode_ == false) {
    // Xpetra mode
    for (size_t i = 0; i < Rows(); i++) {
      for (size_t j = 0; j < Cols(); j++) {
        if (getMatrix(i, j) != Teuchos::null) {
          RCP<const RowMatrix> mat = getMatrix(i, j);

          // recursively call Merge routine
          RCP<const BlockedCrsMatrix> bMat = rcp_dynamic_cast<const BlockedCrsMatrix>(mat);
          if (bMat != Teuchos::null) mat = bMat->Merge();

          bMat = rcp_dynamic_cast<const BlockedCrsMatrix>(mat);
          TEUCHOS_TEST_FOR_EXCEPTION(bMat != Teuchos::null, std::runtime_error,
                                     "BlockedCrsMatrix::Merge: Merging of blocked sub-operators failed?!");

          // jump over empty blocks
          if (mat->getLocalNumEntries() == 0) continue;

          this->Add(*mat, one, *sparse, one);
        }
      }
    }
  } else {
    // Thyra mode
    for (size_t i = 0; i < Rows(); i++) {
      for (size_t j = 0; j < Cols(); j++) {
        if (getMatrix(i, j) != Teuchos::null) {
          RCP<const RowMatrix> mat = getMatrix(i, j);
          // recursively call Merge routine
          RCP<const BlockedCrsMatrix> bMat = rcp_dynamic_cast<const BlockedCrsMatrix>(mat);
          if (bMat != Teuchos::null) mat = bMat->Merge();

          bMat = rcp_dynamic_cast<const BlockedCrsMatrix>(mat);
          TEUCHOS_TEST_FOR_EXCEPTION(bMat != Teuchos::null, std::runtime_error,
                                     "BlockedCrsMatrix::Merge: Merging of blocked sub-operators failed?!");

          // these are the thyra style maps of the matrix
          RCP<const Map> trowMap = mat->getRowMap();
          RCP<const Map> tcolMap = mat->getColMap();
          RCP<const Map> tdomMap = mat->getDomainMap();

          // get Xpetra maps
          RCP<const Map> xrowMap = getRangeMapExtractor()->getMap(i, false);
          RCP<const Map> xdomMap = getDomainMapExtractor()->getMap(j, false);

          // generate column map with Xpetra GIDs
          Teuchos::RCP<Map> xcolMap = transformThyra2XpetraGIDs(*tcolMap, *tdomMap, *xdomMap);

          // jump over empty blocks
          if (mat->getLocalNumEntries() == 0) continue;

          size_t maxNumEntries = mat->getLocalMaxNumRowEntries();

          size_t numEntries;
          nonconst_global_inds_host_view_type inds("inds", maxNumEntries);
          Teuchos::Array<GlobalOrdinal> inds2(maxNumEntries);
          nonconst_values_host_view_type vals("vals", maxNumEntries);

          // loop over all rows and add entries
          for (size_t k = 0; k < mat->getLocalNumRows(); k++) {
            GlobalOrdinal rowTGID = trowMap->getGlobalElement(k);
            mat->getGlobalRowCopy(rowTGID, inds, vals, numEntries);

            // create new indices array
            for (size_t l = 0; l < numEntries; ++l) {
              LocalOrdinal lid = tcolMap->getLocalElement(inds[l]);
              inds2[l]         = xcolMap->getGlobalElement(lid);
            }

            GlobalOrdinal rowXGID = xrowMap->getGlobalElement(k);
            sparse->insertGlobalValues(
                rowXGID, inds2(0, numEntries),
                Teuchos::ArrayView<const Scalar>(reinterpret_cast<const Scalar*>(vals.data()), numEntries)(0, numEntries));
          }
        }
      }
    }
  }

  sparse->fillComplete(getFullDomainMap(), getFullRangeMap());

  TEUCHOS_TEST_FOR_EXCEPTION(sparse->getLocalNumEntries() != getLocalNumEntries(), std::runtime_error,
                             "BlockedCrsMatrix::Merge: Local number of entries of merged matrix does not coincide with local number of entries of blocked operator.");

  TEUCHOS_TEST_FOR_EXCEPTION(sparse->getGlobalNumEntries() != getGlobalNumEntries(), std::runtime_error,
                             "BlockedCrsMatrix::Merge: Global number of entries of merged matrix does not coincide with global number of entries of blocked operator.");

  return sparse;
}

// ---------------------------------------------------------------------------
// Misc
// ---------------------------------------------------------------------------

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetStorageBlockSize() const { return 1; }

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::residual(const MultiVector& X,
                                                                           const MultiVector& B,
                                                                           MultiVector& R) const {
  using STS = Teuchos::ScalarTraits<Scalar>;
  R.update(STS::one(), B, STS::zero());
  this->apply(X, R, Teuchos::NO_TRANS, -STS::one(), STS::one());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Add(const RowMatrix& A, const Scalar scalarA, CrsMatrix& B, const Scalar scalarB) const {
  TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), std::runtime_error, "Matrix A is not completed");

  B.scale(scalarB);

  Scalar one  = Teuchos::ScalarTraits<Scalar>::one();
  Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();

  if (scalarA == zero)
    return;

  const impl_scalar_type implScalarA = static_cast<impl_scalar_type>(scalarA);

  size_t maxNumEntries = A.getLocalMaxNumRowEntries();

  size_t numEntries;
  nonconst_global_inds_host_view_type inds("inds", maxNumEntries);
  nonconst_values_host_view_type vals("vals", maxNumEntries);

  Teuchos::RCP<const Map> rowMap = A.getRowMap();
  for (size_t i = 0; i < A.getLocalNumRows(); i++) {
    GlobalOrdinal row = rowMap->getGlobalElement(static_cast<LocalOrdinal>(i));
    A.getGlobalRowCopy(row, inds, vals, numEntries);

    if (scalarA != one)
      for (size_t j = 0; j < numEntries; ++j)
        vals[j] *= implScalarA;

    // insert should be ok, since blocks in a BlockedCrsMatrix do not overlap!
    B.insertGlobalValues(row,
                         Teuchos::ArrayView<const GlobalOrdinal>(inds.data(), numEntries),
                         Teuchos::ArrayView<const Scalar>(reinterpret_cast<const Scalar*>(vals.data()), numEntries));
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Map> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::transformThyra2XpetraGIDs(const Map& input,
                                                                                                                                                                            const Map& nonOvlInput,
                                                                                                                                                                            const Map& nonOvlReferenceInput) {
  TEUCHOS_TEST_FOR_EXCEPTION(nonOvlInput.getLocalNumElements() != nonOvlReferenceInput.getLocalNumElements(), std::runtime_error,
                             "Tpetra::BlockedCrsMatrix::transformThyra2XpetraGIDs: the number of local Xpetra reference GIDs and local Thyra GIDs of the non-overlapping maps must be the same!");

  Teuchos::RCP<const Teuchos::Comm<int>> comm = input.getComm();

  // fill translation map as far as possible
  std::map<const GlobalOrdinal, GlobalOrdinal> thyra2xpetraGID;
  for (size_t i = 0; i < nonOvlInput.getLocalNumElements(); i++) {
    thyra2xpetraGID[nonOvlInput.getGlobalElement(static_cast<LocalOrdinal>(i))] =
        nonOvlReferenceInput.getGlobalElement(static_cast<LocalOrdinal>(i));
  }

  // find all GIDs of the overlapping Thyra map which are not owned by this proc
  Teuchos::Array<GlobalOrdinal> ovlUnknownStatusGids;
  for (size_t i = 0; i < input.getLocalNumElements(); i++) {
    GlobalOrdinal gcid = input.getGlobalElement(static_cast<LocalOrdinal>(i));
    if (nonOvlInput.isNodeGlobalElement(gcid) == false) {
      ovlUnknownStatusGids.push_back(gcid);
    }
  }

  // Communicate the number of DOFs on each processor
  std::vector<int> myUnknownDofGIDs(comm->getSize(), 0);
  std::vector<int> numUnknownDofGIDs(comm->getSize(), 0);
  myUnknownDofGIDs[comm->getRank()] = static_cast<int>(ovlUnknownStatusGids.size());
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, comm->getSize(), &myUnknownDofGIDs[0], &numUnknownDofGIDs[0]);

  // create array containing all DOF GIDs
  size_t cntUnknownDofGIDs = 0;
  for (int p = 0; p < comm->getSize(); p++) cntUnknownDofGIDs += numUnknownDofGIDs[p];
  std::vector<GlobalOrdinal> lUnknownDofGIDs(cntUnknownDofGIDs, 0);
  std::vector<GlobalOrdinal> gUnknownDofGIDs(cntUnknownDofGIDs, 0);
  size_t cntUnknownOffset = 0;
  for (int p = 0; p < comm->getRank(); p++) cntUnknownOffset += numUnknownDofGIDs[p];
  for (size_t k = 0; k < Teuchos::as<size_t>(ovlUnknownStatusGids.size()); k++) {
    lUnknownDofGIDs[k + cntUnknownOffset] = ovlUnknownStatusGids[k];
  }
  if (cntUnknownDofGIDs > 0)
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, Teuchos::as<int>(cntUnknownDofGIDs), &lUnknownDofGIDs[0], &gUnknownDofGIDs[0]);
  std::vector<GlobalOrdinal> lTranslatedDofGIDs(cntUnknownDofGIDs, 0);
  std::vector<GlobalOrdinal> gTranslatedDofGIDs(cntUnknownDofGIDs, 0);
  for (size_t k = 0; k < gUnknownDofGIDs.size(); k++) {
    GlobalOrdinal curgid = gUnknownDofGIDs[k];
    if (nonOvlInput.isNodeGlobalElement(curgid)) {
      lTranslatedDofGIDs[k] = thyra2xpetraGID[curgid];
    }
  }
  if (cntUnknownDofGIDs > 0)
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, Teuchos::as<int>(cntUnknownDofGIDs), &lTranslatedDofGIDs[0], &gTranslatedDofGIDs[0]);

  for (size_t k = 0; k < Teuchos::as<size_t>(ovlUnknownStatusGids.size()); k++) {
    thyra2xpetraGID[ovlUnknownStatusGids[k]] = gTranslatedDofGIDs[k + cntUnknownOffset];
  }
  Teuchos::Array<GlobalOrdinal> ovlDomainMapArray;
  for (size_t i = 0; i < input.getLocalNumElements(); i++) {
    GlobalOrdinal gcid = input.getGlobalElement(static_cast<LocalOrdinal>(i));
    ovlDomainMapArray.push_back(thyra2xpetraGID[gcid]);
  }
  const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid();
  Teuchos::RCP<Map> ovlDomainMap = Teuchos::rcp(new Map(INVALID, ovlDomainMapArray(), 0, comm));

  TEUCHOS_TEST_FOR_EXCEPTION(input.getLocalNumElements() != ovlDomainMap->getLocalNumElements(), std::runtime_error,
                             "Tpetra::BlockedCrsMatrix::transformThyra2XpetraGIDs: the number of local Thyra reference GIDs (overlapping) and local Xpetra GIDs (overlapping) must be the same!");

  return ovlDomainMap;
}

}  // namespace Tpetra

//
// Explicit instantiation macro
//
#define TPETRA_BLOCKEDCRSMATRIX_INSTANT(SCALAR, LO, GO, NODE) \
  template class BlockedCrsMatrix<SCALAR, LO, GO, NODE>;

#endif  // TPETRA_BLOCKEDCRSMATRIX_DEF_HPP
