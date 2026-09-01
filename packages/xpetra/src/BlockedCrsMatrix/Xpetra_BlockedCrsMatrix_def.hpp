// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_BLOCKEDCRSMATRIX_DEF_HPP
#define XPETRA_BLOCKEDCRSMATRIX_DEF_HPP

#include <Tpetra_KokkosCompat_DefaultNode.hpp>

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_Hashtable.hpp>

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Exceptions.hpp"

#include "Xpetra_BlockedCrsMatrix_decl.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_MultiVector.hpp"
#include "Xpetra_BlockedMultiVector.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#include "Xpetra_BlockedVector.hpp"
#include "Xpetra_CrsGraph.hpp"
#include "Xpetra_CrsMatrix.hpp"
#include "Xpetra_CrsMatrixFactory.hpp"

#include "Xpetra_MapExtractor.hpp"
#include "Xpetra_MapExtractorFactory.hpp"

#include "Xpetra_Matrix.hpp"
#include "Xpetra_MatrixFactory.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"

// Helpers::Op2NonConstTpetraCrs unwraps an Xpetra::Matrix (CrsMatrixWrap) to its
// underlying Tpetra::CrsMatrix, throwing Xpetra::Exceptions::BadCast on a non-Tpetra
// (e.g. Epetra) block -- this is the Epetra guard for the Tpetra-only delegation.
#include "Xpetra_Helpers.hpp"
// TpetraCrsMatrix is used to re-wrap the Tpetra::CrsMatrix returned by op_->Merge().
#include "Xpetra_TpetraCrsMatrix.hpp"

#ifdef HAVE_XPETRA_THYRA
#include <Thyra_ProductVectorSpaceBase.hpp>
#include <Thyra_VectorSpaceBase.hpp>
#include <Thyra_LinearOpBase.hpp>
#include <Thyra_BlockedLinearOpBase.hpp>
#include <Thyra_PhysicallyBlockedLinearOpBase.hpp>
#include "Xpetra_ThyraUtils.hpp"
#endif

#include "Xpetra_VectorFactory.hpp"

/** \file Xpetra_BlockedCrsMatrix.hpp

  Declarations for the class Xpetra::BlockedCrsMatrix.
*/
namespace Xpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BlockedCrsMatrix(const Teuchos::RCP<const BlockedMap>& rangeMaps,
                                                                              const Teuchos::RCP<const BlockedMap>& domainMaps,
                                                                              size_t numEntriesPerRow)
  : is_diagonal_(true) {
  domainmaps_       = MapExtractorFactory::Build(domainMaps);
  rangemaps_        = MapExtractorFactory::Build(rangeMaps);
  bRangeThyraMode_  = rangeMaps->getThyraMode();
  bDomainThyraMode_ = domainMaps->getThyraMode();

  blocks_.reserve(Rows() * Cols());

  // add CrsMatrix objects in row,column order
  for (size_t r = 0; r < Rows(); ++r)
    for (size_t c = 0; c < Cols(); ++c)
      blocks_.push_back(MatrixFactory::Build(getRangeMap(r, bRangeThyraMode_), numEntriesPerRow));

  // Build and sync the wrapped Tpetra core object
  buildAndSyncOp();

  // Default view
  CreateDefaultView();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BlockedCrsMatrix(Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& rangeMapExtractor,
                                                                              Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& domainMapExtractor,
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
      blocks_.push_back(MatrixFactory::Build(getRangeMap(r, bRangeThyraMode_), numEntriesPerRow));

  // Build and sync the wrapped Tpetra core object
  buildAndSyncOp();

  // Default view
  CreateDefaultView();
}

#ifdef HAVE_XPETRA_THYRA
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BlockedCrsMatrix(const Teuchos::RCP<const Thyra::BlockedLinearOpBase<Scalar>>& thyraOp, const Teuchos::RCP<const Teuchos::Comm<int>>& /* comm */)
  : is_diagonal_(true)
  , thyraOp_(thyraOp) {
  // extract information from Thyra blocked operator and rebuilt information
  const Teuchos::RCP<const Thyra::ProductVectorSpaceBase<Scalar>> productRangeSpace  = thyraOp->productRange();
  const Teuchos::RCP<const Thyra::ProductVectorSpaceBase<Scalar>> productDomainSpace = thyraOp->productDomain();

  int numRangeBlocks  = productRangeSpace->numBlocks();
  int numDomainBlocks = productDomainSpace->numBlocks();

  // build range map extractor from Thyra::BlockedLinearOpBase object
  std::vector<Teuchos::RCP<const Map>> subRangeMaps(numRangeBlocks);
  for (size_t r = 0; r < Teuchos::as<size_t>(numRangeBlocks); ++r) {
    for (size_t c = 0; c < Teuchos::as<size_t>(numDomainBlocks); ++c) {
      if (thyraOp->blockExists(r, c)) {
        // we only need at least one block in each block row to extract the range map
        Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> const_op = thyraOp->getBlock(r, c);  // nonConst access is not allowed.
        Teuchos::RCP<const Xpetra::Matrix<Scalar, LO, GO, Node>> xop =
            Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toXpetra(const_op);
        subRangeMaps[r] = xop->getRangeMap();
        if (r != c) is_diagonal_ = false;
        break;
      }
    }
  }
  Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> fullRangeMap = mergeMaps(subRangeMaps);

  // check whether the underlying Thyra operator uses Thyra-style numbering (default for most applications) or
  // Xpetra style numbering
  bool bRangeUseThyraStyleNumbering = false;
  size_t numAllElements             = 0;
  for (size_t v = 0; v < subRangeMaps.size(); ++v) {
    numAllElements += subRangeMaps[v]->getGlobalNumElements();
  }
  if (fullRangeMap->getGlobalNumElements() != numAllElements) bRangeUseThyraStyleNumbering = true;
  rangemaps_ = Teuchos::rcp(new Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>(fullRangeMap, subRangeMaps, bRangeUseThyraStyleNumbering));

  // build domain map extractor from Thyra::BlockedLinearOpBase object
  std::vector<Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>> subDomainMaps(numDomainBlocks);
  for (size_t c = 0; c < Teuchos::as<size_t>(numDomainBlocks); ++c) {
    for (size_t r = 0; r < Teuchos::as<size_t>(numRangeBlocks); ++r) {
      if (thyraOp->blockExists(r, c)) {
        // we only need at least one block in each block row to extract the range map
        Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> const_op = thyraOp->getBlock(r, c);  // nonConst access is not allowed.
        Teuchos::RCP<const Xpetra::Matrix<Scalar, LO, GO, Node>> xop =
            Xpetra::ThyraUtils<Scalar, LO, GO, Node>::toXpetra(const_op);
        subDomainMaps[c] = xop->getDomainMap();
        break;
      }
    }
  }
  Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> fullDomainMap = mergeMaps(subDomainMaps);
  // plausibility check for numbering style (Xpetra or Thyra)
  bool bDomainUseThyraStyleNumbering = false;
  numAllElements                     = 0;
  for (size_t v = 0; v < subDomainMaps.size(); ++v) {
    numAllElements += subDomainMaps[v]->getGlobalNumElements();
  }
  if (fullDomainMap->getGlobalNumElements() != numAllElements) bDomainUseThyraStyleNumbering = true;
  domainmaps_ = Teuchos::rcp(new Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>(fullDomainMap, subDomainMaps, bDomainUseThyraStyleNumbering));

  // store numbering mode
  bRangeThyraMode_  = bRangeUseThyraStyleNumbering;
  bDomainThyraMode_ = bDomainUseThyraStyleNumbering;

  // add CrsMatrix objects in row,column order
  blocks_.reserve(Rows() * Cols());
  for (size_t r = 0; r < Rows(); ++r) {
    for (size_t c = 0; c < Cols(); ++c) {
      if (thyraOp->blockExists(r, c)) {
        // TODO we do not support nested Thyra operators here!
        Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> const_op = thyraOp->getBlock(r, c);                                         // nonConst access is not allowed.
        Teuchos::RCP<Thyra::LinearOpBase<Scalar>> op             = Teuchos::rcp_const_cast<Thyra::LinearOpBase<Scalar>>(const_op);  // cast away const
        Teuchos::RCP<Xpetra::Matrix<Scalar, LO, GO, Node>> xop =
            Xpetra::ThyraUtils<Scalar, LO, GO, Node>::toXpetra(op);
        Teuchos::RCP<Xpetra::CrsMatrixWrap<Scalar, LO, GO, Node>> xwrap =
            Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrixWrap<Scalar, LO, GO, Node>>(xop, true);
        blocks_.push_back(xwrap);
      } else {
        // add empty block
        blocks_.push_back(MatrixFactory::Build(getRangeMap(r, bRangeThyraMode_), 0));
      }
    }
  }

  // Build and sync the wrapped Tpetra core object
  buildAndSyncOp();

  // Default view
  CreateDefaultView();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::mergeMaps(std::vector<Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>>& subMaps) {
  // TODO merging for Thyra mode is missing (similar to what we do in constructor of MapExtractor

  // merge submaps to global map
  std::vector<GlobalOrdinal> gids;
  for (size_t tt = 0; tt < subMaps.size(); ++tt) {
    Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> subMap = subMaps[tt];
#if 1
    Teuchos::ArrayView<const GlobalOrdinal> subMapGids = subMap->getLocalElementList();
    gids.insert(gids.end(), subMapGids.begin(), subMapGids.end());
#else
    size_t myNumElements = subMap->getLocalNumElements();
    for (LocalOrdinal l = 0; l < Teuchos::as<LocalOrdinal>(myNumElements); ++l) {
      GlobalOrdinal gid = subMap->getGlobalElement(l);
      gids.push_back(gid);
    }
#endif
  }

  // we have to sort the matrix entries and get rid of the double entries
  // since we use this to detect Thyra-style numbering or Xpetra-style
  // numbering. In Thyra-style numbering mode, the Xpetra::MapExtractor builds
  // the correct row maps.
  const GO INVALID = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
  std::sort(gids.begin(), gids.end());
  gids.erase(std::unique(gids.begin(), gids.end()), gids.end());
  Teuchos::ArrayView<GO> gidsView(&gids[0], gids.size());
  Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> fullMap = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(subMaps[0]->lib(), INVALID, gidsView, subMaps[0]->getIndexBase(), subMaps[0]->getComm());
  return fullMap;
}

#endif

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~BlockedCrsMatrix() = default;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<::Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::unwrapBlock(const Teuchos::RCP<Matrix>& mat) const {
  if (mat.is_null())
    return Teuchos::null;

  // A nested blocked block unwraps to its own wrapped Tpetra::BlockedCrsMatrix.
  Teuchos::RCP<BlockedCrsMatrix> bmat = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(mat);
  if (!bmat.is_null())
    return bmat->getTpetra_BlockedCrsMatrix();

  // A leaf CrsMatrixWrap unwraps to its Tpetra::CrsMatrix.  Op2NonConstTpetraCrs throws
  // Xpetra::Exceptions::BadCast if the block is not Tpetra-backed (e.g. Epetra).
  return Xpetra::Helpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Op2NonConstTpetraCrs(mat);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::buildAndSyncOp() {
  // Build the Tpetra core from the same map extractors this wrapper holds.  The core
  // constructor auto-builds its own (empty) Tpetra blocks; we immediately overwrite them
  // with the unwrapped Xpetra cache blocks so that op_ and the cache share the same
  // underlying Tpetra objects.  The numEntriesPerRow argument is irrelevant here since the
  // auto-built blocks are replaced; pass 0.
  Teuchos::RCP<const ::Tpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tRange =
      rangemaps_->getTpetra_MapExtractor();
  Teuchos::RCP<const ::Tpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tDomain =
      domainmaps_->getTpetra_MapExtractor();
  op_ = Teuchos::rcp(new tpetra_blockedcrsmatrix_type(tRange, tDomain, 0));

  for (size_t r = 0; r < Rows(); ++r)
    for (size_t c = 0; c < Cols(); ++c)
      op_->setMatrix(r, c, unwrapBlock(blocks_[r * Cols() + c]));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::insertGlobalValues(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal>& cols, const ArrayView<const Scalar>& vals) {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::insertGlobalValues");
  if (Rows() == 1 && Cols() == 1) {
    getMatrix(0, 0)->insertGlobalValues(globalRow, cols, vals);
    return;
  }
  throw Xpetra::Exceptions::RuntimeError("insertGlobalValues not supported by BlockedCrsMatrix");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::insertLocalValues(LocalOrdinal localRow, const ArrayView<const LocalOrdinal>& cols, const ArrayView<const Scalar>& vals) {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::insertLocalValues");
  if (Rows() == 1 && Cols() == 1) {
    getMatrix(0, 0)->insertLocalValues(localRow, cols, vals);
    return;
  }
  throw Xpetra::Exceptions::RuntimeError("insertLocalValues not supported by BlockedCrsMatrix");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::removeEmptyProcessesInPlace(const Teuchos::RCP<const Map>& newMap) {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::removeEmptyProcessesInPlace");
  if (Rows() == 1 && Cols() == 1) {
    getMatrix(0, 0)->removeEmptyProcessesInPlace(newMap);
    return;
  }
  throw Xpetra::Exceptions::RuntimeError("removeEmptyProcesses not supported by BlockedCrsMatrix");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceGlobalValues(GlobalOrdinal globalRow,
                                                                                      const ArrayView<const GlobalOrdinal>& cols,
                                                                                      const ArrayView<const Scalar>& vals) {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::replaceGlobalValues");
  if (Rows() == 1 && Cols() == 1) {
    getMatrix(0, 0)->replaceGlobalValues(globalRow, cols, vals);
    return;
  }
  throw Xpetra::Exceptions::RuntimeError("replaceGlobalValues not supported by BlockedCrsMatrix");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceLocalValues(LocalOrdinal localRow,
                                                                                     const ArrayView<const LocalOrdinal>& cols,
                                                                                     const ArrayView<const Scalar>& vals) {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::replaceLocalValues");
  if (Rows() == 1 && Cols() == 1) {
    getMatrix(0, 0)->replaceLocalValues(localRow, cols, vals);
    return;
  }
  throw Xpetra::Exceptions::RuntimeError("replaceLocalValues not supported by BlockedCrsMatrix");
}

//! Set all matrix entries equal to scalar
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::setAllToScalar(const Scalar& alpha) {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::setAllToScalar");
  for (size_t row = 0; row < Rows(); row++) {
    for (size_t col = 0; col < Cols(); col++) {
      if (!getMatrix(row, col).is_null()) {
        getMatrix(row, col)->setAllToScalar(alpha);
      }
    }
  }
}

//! Scale the current values of a matrix, this = alpha*this.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::scale(const Scalar& alpha) {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::scale");
  for (size_t row = 0; row < Rows(); row++) {
    for (size_t col = 0; col < Cols(); col++) {
      if (!getMatrix(row, col).is_null()) {
        getMatrix(row, col)->scale(alpha);
      }
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::resumeFill(const RCP<ParameterList>& params) {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::resumeFill");
  for (size_t row = 0; row < Rows(); row++) {
    for (size_t col = 0; col < Cols(); col++) {
      if (!getMatrix(row, col).is_null()) {
        getMatrix(row, col)->resumeFill(params);
      }
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::fillComplete(const RCP<const Map>& domainMap, const RCP<const Map>& rangeMap, const RCP<ParameterList>& params) {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::fillComplete");
  if (Rows() == 1 && Cols() == 1) {
    getMatrix(0, 0)->fillComplete(domainMap, rangeMap, params);
    return;
  }
  fillComplete(params);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::fillComplete(const RCP<ParameterList>& params) {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::fillComplete");
  TEUCHOS_TEST_FOR_EXCEPTION(rangemaps_ == Teuchos::null, Xpetra::Exceptions::RuntimeError, "BlockedCrsMatrix::fillComplete: rangemaps_ is not set. Error.");

  for (size_t r = 0; r < Rows(); ++r)
    for (size_t c = 0; c < Cols(); ++c) {
      if (getMatrix(r, c) != Teuchos::null) {
        Teuchos::RCP<Matrix> Ablock = getMatrix(r, c);
        if (r != c) is_diagonal_ = false;
        if (Ablock != Teuchos::null && !Ablock->isFillComplete())
          Ablock->fillComplete(getDomainMap(c, bDomainThyraMode_), getRangeMap(r, bRangeThyraMode_), params);
      }
    }

#if 0
      // get full row map
      RCP<const Map> rangeMap = rangemaps_->getFullMap();
      fullrowmap_ = MapFactory::Build(rangeMap()->lib(), rangeMap()->getGlobalNumElements(), rangeMap()->getLocalElementList(), rangeMap()->getIndexBase(), rangeMap()->getComm());

      // build full col map
      fullcolmap_ = Teuchos::null; // delete old full column map

      std::vector<GO> colmapentries;
      for (size_t c = 0; c < Cols(); ++c) {
        // copy all local column lids of all block rows to colset
        std::set<GO> colset;
        for (size_t r = 0; r < Rows(); ++r) {
          Teuchos::RCP<CrsMatrix> Ablock = getMatrix(r,c);

          if (Ablock != Teuchos::null) {
            Teuchos::ArrayView<const GO> colElements = Ablock->getColMap()->getLocalElementList();
            Teuchos::RCP<const Map>      colmap      = Ablock->getColMap();
            copy(colElements.getRawPtr(), colElements.getRawPtr() + colElements.size(), inserter(colset, colset.begin()));
          }
        }

        // remove duplicates (entries which are in column maps of more than one block row)
        colmapentries.reserve(colmapentries.size() + colset.size());
        copy(colset.begin(), colset.end(), back_inserter(colmapentries));
        sort(colmapentries.begin(), colmapentries.end());
        typename std::vector<GO>::iterator gendLocation;
        gendLocation = std::unique(colmapentries.begin(), colmapentries.end());
        colmapentries.erase(gendLocation,colmapentries.end());
      }

      // sum up number of local elements
      size_t numGlobalElements = 0;
      Teuchos::reduceAll(*(rangeMap->getComm()), Teuchos::REDUCE_SUM, colmapentries.size(), Teuchos::outArg(numGlobalElements));

      // store global full column map
      const Teuchos::ArrayView<const GO> aView = Teuchos::ArrayView<const GO>(colmapentries);
      fullcolmap_ = MapFactory::Build(rangeMap->lib(), numGlobalElements, aView, 0, rangeMap->getComm());
#endif
}

// The RowMatrix aggregate queries below (counts, indexing state, fill state) are the
// real block-reduction logic and live in the Tpetra core; the wrapper simply forwards.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
global_size_t BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getGlobalNumRows() const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getGlobalNumRows");
  return op_->getGlobalNumRows();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
global_size_t BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getGlobalNumCols() const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getGlobalNumCols");
  return op_->getGlobalNumCols();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getLocalNumRows() const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getLocalNumRows");
  return op_->getLocalNumRows();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
global_size_t BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getGlobalNumEntries() const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getGlobalNumEntries");
  return op_->getGlobalNumEntries();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getLocalNumEntries() const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getLocalNumEntries");
  return op_->getLocalNumEntries();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNumEntriesInLocalRow(LocalOrdinal localRow) const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getNumEntriesInLocalRow");
  return op_->getNumEntriesInLocalRow(localRow);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getNumEntriesInGlobalRow");
  return op_->getNumEntriesInGlobalRow(globalRow);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getGlobalMaxNumRowEntries() const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getGlobalMaxNumRowEntries");
  return op_->getGlobalMaxNumRowEntries();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getLocalMaxNumRowEntries() const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getLocalMaxNumRowEntries");
  return op_->getLocalMaxNumRowEntries();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isLocallyIndexed() const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::isLocallyIndexed");
  return op_->isLocallyIndexed();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isGloballyIndexed() const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::isGloballyIndexed");
  return op_->isGloballyIndexed();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isFillComplete() const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::isFillComplete");
  return op_->isFillComplete();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getLocalRowCopy(LocalOrdinal LocalRow,
                                                                                  const ArrayView<LocalOrdinal>& Indices,
                                                                                  const ArrayView<Scalar>& Values,
                                                                                  size_t& NumEntries) const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getLocalRowCopy");
  if (Rows() == 1 && Cols() == 1) {
    getMatrix(0, 0)->getLocalRowCopy(LocalRow, Indices, Values, NumEntries);
    return;
  }
  throw Xpetra::Exceptions::RuntimeError("getLocalRowCopy not supported by BlockedCrsMatrix");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal>& indices, ArrayView<const Scalar>& values) const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getGlobalRowView");
  if (Rows() == 1 && Cols() == 1) {
    getMatrix(0, 0)->getGlobalRowView(GlobalRow, indices, values);
    return;
  }
  throw Xpetra::Exceptions::RuntimeError("getGlobalRowView not supported by BlockedCrsMatrix");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getLocalRowView(LocalOrdinal LocalRow, ArrayView<const LocalOrdinal>& indices, ArrayView<const Scalar>& values) const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getLocalRowView");
  if (Rows() == 1 && Cols() == 1) {
    getMatrix(0, 0)->getLocalRowView(LocalRow, indices, values);
    return;
  } else if (is_diagonal_) {
    GlobalOrdinal gid = this->getRowMap()->getGlobalElement(LocalRow);
    size_t row        = getBlockedRangeMap()->getMapIndexForGID(gid);
    getMatrix(row, row)->getLocalRowView(getMatrix(row, row)->getRowMap()->getLocalElement(gid), indices, values);
    return;
  }
  throw Xpetra::Exceptions::RuntimeError("getLocalRowView not supported by BlockedCrsMatrix");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getLocalDiagCopy(Vector& diag) const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getLocalDiagCopy");

  // RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  Teuchos::RCP<Vector> rcpdiag      = Teuchos::rcpFromRef(diag);
  Teuchos::RCP<BlockedVector> bdiag = Teuchos::rcp_dynamic_cast<BlockedVector>(rcpdiag);

  // special treatment for 1x1 block matrices
  // ReorderedBlockedCrsMatrix object encapsulate single blocks in ReorderedBlocks
  // BlockedVectors have Vector objects as Leaf objects.
  if (Rows() == 1 && Cols() == 1 && bdiag.is_null() == true) {
    Teuchos::RCP<const Matrix> rm = getMatrix(0, 0);
    rm->getLocalDiagCopy(diag);
    return;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(bdiag.is_null() == true, Xpetra::Exceptions::RuntimeError, "BlockedCrsMatrix::getLocalDiagCopy: diag must be a Blocked(Multi)Vector.");
  TEUCHOS_TEST_FOR_EXCEPTION(bdiag->getNumVectors() != 1, Xpetra::Exceptions::RuntimeError, "BlockedCrsMatrix::getLocalDiagCopy: diag must be a Blocked(Multi)Vector with exactly one vector. However, the number of stored vectors is " << bdiag->getNumVectors());
  TEUCHOS_TEST_FOR_EXCEPTION(bdiag->getBlockedMap()->getNumMaps() != this->Rows(), Xpetra::Exceptions::RuntimeError,
                             "BlockedCrsMatrix::getLocalDiagCopy(): the number of blocks in diag differ from the number of blocks in this operator.");
  // XPETRA_TEST_FOR_EXCEPTION(bdiag->getMap()->isSameAs(*(getMap())) == false, Xpetra::Exceptions::RuntimeError,
  //   "BlockedCrsMatrix::getLocalDiagCopy(): the map of the vector diag is not compatible with the map of the blocked operator." );

  for (size_t row = 0; row < Rows(); row++) {
    Teuchos::RCP<const Matrix> rm = getMatrix(row, row);
    if (!rm.is_null()) {
      Teuchos::RCP<Vector> rv = VectorFactory::Build(bdiag->getBlockedMap()->getMap(row, bdiag->getBlockedMap()->getThyraMode()));
      rm->getLocalDiagCopy(*rv);
      bdiag->setMultiVector(row, rv, bdiag->getBlockedMap()->getThyraMode());
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::leftScale(const Vector& x) {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::leftScale");

  Teuchos::RCP<const Vector> rcpx      = Teuchos::rcpFromRef(x);
  Teuchos::RCP<const BlockedVector> bx = Teuchos::rcp_dynamic_cast<const BlockedVector>(rcpx);

  // RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  // special treatment for 1xn block matrices
  // ReorderedBlockedCrsMatrix object encapsulate single blocks in ReorderedBlocks
  // BlockedVectors have Vector objects as Leaf objects.
  if (Rows() == 1 && bx.is_null() == true) {
    for (size_t col = 0; col < Cols(); ++col) {
      Teuchos::RCP<Matrix> rm = getMatrix(0, col);
      rm->leftScale(x);
    }
    return;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(bx.is_null() == true, Xpetra::Exceptions::RuntimeError, "BlockedCrsMatrix::leftScale: x must be a Blocked(Multi)Vector.");
  TEUCHOS_TEST_FOR_EXCEPTION(bx->getNumVectors() != 1, Xpetra::Exceptions::RuntimeError, "BlockedCrsMatrix::leftScale: x must be a Blocked(Multi)Vector with exactly one vector. However, the number of stored vectors is " << bx->getNumVectors());
  TEUCHOS_TEST_FOR_EXCEPTION(bx->getBlockedMap()->getNumMaps() != this->Rows(), Xpetra::Exceptions::RuntimeError,
                             "BlockedCrsMatrix::leftScale(): the number of blocks in diag differ from the number of blocks in this operator.");

  for (size_t row = 0; row < Rows(); row++) {
    Teuchos::RCP<const MultiVector> rmv = bx->getMultiVector(row);
    Teuchos::RCP<const Vector> rscale   = rmv->getVector(0);
    XPETRA_TEST_FOR_EXCEPTION(rscale.is_null() == true, Xpetra::Exceptions::RuntimeError, "BlockedCrsMatrix::leftScale: x must be a Vector.");
    for (size_t col = 0; col < Cols(); ++col) {
      Teuchos::RCP<Matrix> rm = getMatrix(row, col);
      if (!rm.is_null()) {
        rm->leftScale(*rscale);
      }
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::rightScale(const Vector& x) {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::rightScale");

  Teuchos::RCP<const Vector> rcpx      = Teuchos::rcpFromRef(x);
  Teuchos::RCP<const BlockedVector> bx = Teuchos::rcp_dynamic_cast<const BlockedVector>(rcpx);

  // RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  // special treatment for nx1 block matrices
  // ReorderedBlockedCrsMatrix object encapsulate single blocks in ReorderedBlocks
  // BlockedVectors have Vector objects as Leaf objects.
  if (Cols() == 1 && bx.is_null() == true) {
    for (size_t row = 0; row < Rows(); ++row) {
      Teuchos::RCP<Matrix> rm = getMatrix(row, 0);
      rm->rightScale(x);
    }
    return;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(bx.is_null() == true, Xpetra::Exceptions::RuntimeError, "BlockedCrsMatrix::rightScale: x must be a Blocked(Multi)Vector.");
  TEUCHOS_TEST_FOR_EXCEPTION(bx->getNumVectors() != 1, Xpetra::Exceptions::RuntimeError, "BlockedCrsMatrix::rightScale: x must be a Blocked(Multi)Vector with exactly one vector. However, the number of stored vectors is " << bx->getNumVectors());
  TEUCHOS_TEST_FOR_EXCEPTION(bx->getBlockedMap()->getNumMaps() != this->Cols(), Xpetra::Exceptions::RuntimeError,
                             "BlockedCrsMatrix::rightScale(): the number of blocks in diag differ from the number of blocks in this operator.");

  for (size_t col = 0; col < Cols(); ++col) {
    Teuchos::RCP<const MultiVector> rmv = bx->getMultiVector(col);
    Teuchos::RCP<const Vector> rscale   = rmv->getVector(0);
    XPETRA_TEST_FOR_EXCEPTION(rscale.is_null() == true, Xpetra::Exceptions::RuntimeError, "BlockedCrsMatrix::leftScale: x must be a Vector.");
    for (size_t row = 0; row < Rows(); row++) {
      Teuchos::RCP<Matrix> rm = getMatrix(row, col);
      if (!rm.is_null()) {
        rm->rightScale(*rscale);
      }
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename ScalarTraits<Scalar>::magnitudeType BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getFrobeniusNorm() const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getFrobeniusNorm");
  return op_->getFrobeniusNorm();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::haveGlobalConstants() const { return true; }

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::apply(const MultiVector& X, MultiVector& Y, Teuchos::ETransp mode, Scalar alpha, Scalar beta, bool sumInterfaceValues,
                                                                        const RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>>& regionInterfaceImporter,
                                                                        const Teuchos::ArrayRCP<LocalOrdinal>& regionInterfaceLIDs) const {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::apply(const MultiVector& X, MultiVector& Y,
                                                                        Teuchos::ETransp mode,
                                                                        Scalar alpha,
                                                                        Scalar beta) const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::apply");

  // The full block-apply algorithm (blocked/non-blocked detection, Thyra/Xpetra GID
  // handling, per-block dispatch) lives in the Tpetra core.  Unwrap X and Y to their
  // underlying Tpetra objects -- blocked Xpetra vectors unwrap to Tpetra::BlockedMultiVector
  // so op_ sees the same blocked structure -- and forward.  Bind explicitly-typed RCP
  // intermediates first to avoid an unwrapMultiVector overload ambiguity.
  Teuchos::RCP<const MultiVector> refX = Teuchos::rcpFromRef(X);
  Teuchos::RCP<MultiVector> refY       = Teuchos::rcpFromRef(Y);
  Teuchos::RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tX =
      BlockedMultiVectorDetails::unwrapMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(refX);
  Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tY =
      BlockedMultiVectorDetails::unwrapMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(refY);
  op_->apply(*tX, *tY, mode, alpha, beta);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getFullDomainMap() const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getFullDomainMap()");
  return domainmaps_->getFullMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getBlockedDomainMap() const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getBlockedDomainMap()");
  return domainmaps_->getBlockedMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getDomainMap() const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getDomainMap()");
  return domainmaps_->getMap(); /*domainmaps_->getFullMap();*/
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getDomainMap(size_t i) const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getDomainMap(size_t)");
  return domainmaps_->getMap(i, bDomainThyraMode_);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getDomainMap(size_t i, bool bThyraMode) const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getDomainMap(size_t,bool)");
  return domainmaps_->getMap(i, bThyraMode);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getFullRangeMap() const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getRangeMap()");
  return rangemaps_->getFullMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getBlockedRangeMap() const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getBlockedRangeMap()");
  return rangemaps_->getBlockedMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getRangeMap() const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getRangeMap()");
  return rangemaps_->getMap(); /*rangemaps_->getFullMap();*/
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getRangeMap(size_t i) const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getRangeMap(size_t)");
  return rangemaps_->getMap(i, bRangeThyraMode_);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getRangeMap(size_t i, bool bThyraMode) const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getRangeMap(size_t,bool)");
  return rangemaps_->getMap(i, bThyraMode);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getRangeMapExtractor() const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getRangeMapExtractor()");
  return rangemaps_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getDomainMapExtractor() const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getDomainMapExtractor()");
  return domainmaps_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::bgs_apply(
    const MultiVector& X,   ///< Vector to be multiplied by matrix (input)
    MultiVector& Y,         ///< result vector
    size_t row,             ///< Index of block row to be treated
    Teuchos::ETransp mode,  ///< Transpose mode
    Scalar alpha,           ///< scaling factor for result of matrix-vector product
    Scalar beta             ///< scaling factor for linear combination with result vector
) const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::bgs_apply");

  // Block-Gauss-Seidel single-block-row apply lives in the Tpetra core.  Unwrap X/Y as in
  // apply() (blocked Xpetra vectors unwrap to Tpetra::BlockedMultiVector) and forward.
  Teuchos::RCP<const MultiVector> refX = Teuchos::rcpFromRef(X);
  Teuchos::RCP<MultiVector> refY       = Teuchos::rcpFromRef(Y);
  Teuchos::RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tX =
      BlockedMultiVectorDetails::unwrapMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(refX);
  Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tY =
      BlockedMultiVectorDetails::unwrapMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(refY);
  op_->bgs_apply(*tX, *tY, row, mode, alpha, beta);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getMap() const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getMap");
  if (Rows() == 1 && Cols() == 1) {
    return getMatrix(0, 0)->getMap();
  }
  throw Xpetra::Exceptions::RuntimeError("BlockedCrsMatrix::getMap(): operation not supported.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doImport(const Matrix& source, const Import& importer, CombineMode CM) {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::doImport");
  if (Rows() == 1 && Cols() == 1) {
    getMatrix(0, 0)->doImport(source, importer, CM);
    return;
  }
  throw Xpetra::Exceptions::RuntimeError("BlockedCrsMatrix::doImport(): operation not supported.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doExport(const Matrix& dest, const Import& importer, CombineMode CM) {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::doExport");
  if (Rows() == 1 && Cols() == 1) {
    getMatrix(0, 0)->doExport(dest, importer, CM);
    return;
  }
  throw Xpetra::Exceptions::RuntimeError("BlockedCrsMatrix::doExport(): operation not supported.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doImport(const Matrix& source, const Export& exporter, CombineMode CM) {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::doImport");
  if (Rows() == 1 && Cols() == 1) {
    getMatrix(0, 0)->doImport(source, exporter, CM);
    return;
  }
  throw Xpetra::Exceptions::RuntimeError("BlockedCrsMatrix::doImport(): operation not supported.");
}

//! Export (using an Importer).
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doExport(const Matrix& dest, const Export& exporter, CombineMode CM) {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::doExport");
  if (Rows() == 1 && Cols() == 1) {
    getMatrix(0, 0)->doExport(dest, exporter, CM);
    return;
  }
  throw Xpetra::Exceptions::RuntimeError("BlockedCrsMatrix::doExport(): operation not supported.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const { return "Xpetra_BlockedCrsMatrix.description()"; }

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const {
  out << "Xpetra::BlockedCrsMatrix: " << Rows() << " x " << Cols() << std::endl;

  if (isFillComplete()) {
    out << "BlockMatrix is fillComplete" << std::endl;

    /*if(fullrowmap_ != Teuchos::null) {
      out << "fullRowMap" << std::endl;
      fullrowmap_->describe(out,verbLevel);
    } else {
      out << "fullRowMap not set. Check whether block matrix is properly fillCompleted!" << std::endl;
    }*/

    // out << "fullColMap" << std::endl;
    // fullcolmap_->describe(out,verbLevel);

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
  XPETRA_MONITOR("TpetraBlockedCrsMatrix::setObjectLabel");
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
  if (Rows() == 1 && Cols() == 1)
    return true;
  else
    return false;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getCrsGraph() const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getCrsGraph");
  if (Rows() == 1 && Cols() == 1) {
    return getMatrix(0, 0)->getCrsGraph();
  }
  throw Xpetra::Exceptions::RuntimeError("getCrsGraph() not supported by BlockedCrsMatrix");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isDiagonal() const { return is_diagonal_; }

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Rows() const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::Rows");
  return rangemaps_->NumMaps();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Cols() const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::Cols");
  return domainmaps_->NumMaps();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getCrsMatrix() const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getCrsMatrix");
  TEUCHOS_TEST_FOR_EXCEPTION(Rows() != 1, std::out_of_range, "Can only unwrap a 1x1 blocked matrix. The matrix has " << Rows() << " block rows, though.");
  TEUCHOS_TEST_FOR_EXCEPTION(Cols() != 1, std::out_of_range, "Can only unwrap a 1x1 blocked matrix. The matrix has " << Cols() << " block columns, though.");

  RCP<Matrix> mat            = getMatrix(0, 0);
  RCP<BlockedCrsMatrix> bmat = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(mat);
  if (bmat == Teuchos::null) return mat;
  return bmat->getCrsMatrix();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getInnermostCrsMatrix() {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getInnermostCrsMatrix");
  size_t row = Rows() + 1, col = Cols() + 1;
  for (size_t r = 0; r < Rows(); ++r)
    for (size_t c = 0; c < Cols(); ++c)
      if (getMatrix(r, c) != Teuchos::null) {
        row = r;
        col = c;
        break;
      }
  TEUCHOS_TEST_FOR_EXCEPTION(row == Rows() + 1 || col == Cols() + 1, Xpetra::Exceptions::Incompatible, "Xpetra::BlockedCrsMatrix::getInnermostCrsMatrix: Could not find a non-zero sub-block in blocked operator.")
  RCP<Matrix> mm             = getMatrix(row, col);
  RCP<BlockedCrsMatrix> bmat = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(mm);
  if (bmat == Teuchos::null) return mm;
  return bmat->getInnermostCrsMatrix();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getMatrix(size_t r, size_t c) const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::getMatrix");
  TEUCHOS_TEST_FOR_EXCEPTION(r > Rows(), std::out_of_range, "Error, r = " << Rows() << " is too big");
  TEUCHOS_TEST_FOR_EXCEPTION(c > Cols(), std::out_of_range, "Error, c = " << Cols() << " is too big");

  // transfer strided/blocked map information
  /*      if (blocks_[r*Cols()+c] != Teuchos::null &&
      blocks_[r*Cols()+c]->IsView("stridedMaps") == false)
      blocks_[r*Cols()+c]->CreateView("stridedMaps", getRangeMap(r,bRangeThyraMode_), getDomainMap(c,bDomainThyraMode_));*/
  return blocks_[r * Cols() + c];
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::setMatrix(size_t r, size_t c, Teuchos::RCP<Matrix> mat) {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::setMatrix");
  // TODO: if filled -> return error

  TEUCHOS_TEST_FOR_EXCEPTION(r > Rows(), std::out_of_range, "Error, r = " << Rows() << " is too big");
  TEUCHOS_TEST_FOR_EXCEPTION(c > Cols(), std::out_of_range, "Error, c = " << Cols() << " is too big");
  if (!mat.is_null() && r != c) is_diagonal_ = false;
  // set matrix in the Xpetra identity cache ...
  blocks_[r * Cols() + c] = mat;
  // ... and, unwrapped, in the Tpetra core (throws BadCast on an Epetra-backed block).
  op_->setMatrix(r, c, unwrapBlock(mat));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Merge() const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::Merge");

  // The full merge (Xpetra- and Thyra-style GID handling, recursive merge of nested
  // blocked sub-operators, assembly into a single fill-completed CrsMatrix) lives in the
  // Tpetra core.  Forward to it and re-wrap the returned Tpetra::CrsMatrix as an
  // Xpetra::Matrix (CrsMatrixWrap around a TpetraCrsMatrix), matching the legacy return type.
  Teuchos::RCP<::Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tmerged = op_->Merge();
  Teuchos::RCP<CrsMatrix> xmerged = Teuchos::rcp(new TpetraCrsMatrix(tmerged));
  return Teuchos::rcp(new CrsMatrixWrap(xmerged));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getLocalMatrixDevice() const {
  if (Rows() == 1 && Cols() == 1) {
    return getMatrix(0, 0)->getLocalMatrixDevice();
  }
  throw Xpetra::Exceptions::RuntimeError("BlockedCrsMatrix::getLocalMatrix(): operation not supported.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type::host_mirror_type BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getLocalMatrixHost() const {
  if (Rows() == 1 && Cols() == 1) {
    return getMatrix(0, 0)->getLocalMatrixHost();
  }
  throw Xpetra::Exceptions::RuntimeError("BlockedCrsMatrix::getLocalMatrix(): operation not supported.");
}

#ifdef HAVE_XPETRA_THYRA
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Thyra::BlockedLinearOpBase<Scalar>> BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getThyraOperator() {
  Teuchos::RCP<Thyra::LinearOpBase<Scalar>> thOp =
      Xpetra::ThyraUtils<Scalar, LO, GO, Node>::toThyra(Teuchos::rcpFromRef(*this));
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thOp));

  Teuchos::RCP<Thyra::BlockedLinearOpBase<Scalar>> thbOp =
      Teuchos::rcp_dynamic_cast<Thyra::BlockedLinearOpBase<Scalar>>(thOp);
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thbOp));
  return thbOp;
}
#endif

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
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Add(const Matrix& A, const Scalar scalarA, Matrix& B, const Scalar scalarB) const {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::Add");
  TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), Xpetra::Exceptions::RuntimeError,
                             "Matrix A is not completed");
  using Teuchos::Array;
  using Teuchos::ArrayView;

  B.scale(scalarB);

  Scalar one  = ScalarTraits<SC>::one();
  Scalar zero = ScalarTraits<SC>::zero();

  if (scalarA == zero)
    return;

  Teuchos::RCP<const Matrix> rcpA            = Teuchos::rcpFromRef(A);
  Teuchos::RCP<const CrsMatrixWrap> rcpAwrap = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(rcpA);
  TEUCHOS_TEST_FOR_EXCEPTION(rcpAwrap == Teuchos::null, Xpetra::Exceptions::BadCast,
                             "BlockedCrsMatrix::Add: matrix A must be of type CrsMatrixWrap.");
  Teuchos::RCP<const CrsMatrix> crsA = rcpAwrap->getCrsMatrix();

  size_t maxNumEntries = crsA->getLocalMaxNumRowEntries();

  size_t numEntries;
  Array<GO> inds(maxNumEntries);
  Array<SC> vals(maxNumEntries);

  RCP<const Map> rowMap = crsA->getRowMap();
  RCP<const Map> colMap = crsA->getColMap();

  ArrayView<const GO> rowGIDs = crsA->getRowMap()->getLocalElementList();
  for (size_t i = 0; i < crsA->getLocalNumRows(); i++) {
    GO row = rowGIDs[i];
    crsA->getGlobalRowCopy(row, inds(), vals(), numEntries);

    if (scalarA != one)
      for (size_t j = 0; j < numEntries; ++j)
        vals[j] *= scalarA;

    B.insertGlobalValues(row, inds(0, numEntries), vals(0, numEntries));  // insert should be ok, since blocks in BlockedCrsOpeartor do not overlap!
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateDefaultView() {
  XPETRA_MONITOR("XpetraBlockedCrsMatrix::CreateDefaultView");

  // Create default view
  this->defaultViewLabel_ = "point";
  this->CreateView(this->GetDefaultViewLabel(), getRangeMap(), getDomainMap());

  // Set current view
  this->currentViewLabel_ = this->GetDefaultViewLabel();
}

}  // namespace Xpetra

#define XPETRA_BLOCKEDCRSMATRIX_SHORT
#endif /* XPETRA_BLOCKEDCRSMATRIX_DEF_HPP */
