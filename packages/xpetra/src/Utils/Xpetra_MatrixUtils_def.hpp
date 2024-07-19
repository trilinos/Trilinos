// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PACKAGES_XPETRA_SUP_MATRIX_UTILS_DEF_HPP_
#define PACKAGES_XPETRA_SUP_MATRIX_UTILS_DEF_HPP_

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_Map.hpp"
#include "Xpetra_MapUtils.hpp"
#include "Xpetra_StridedMap.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_MapExtractor.hpp"
#include "Xpetra_MapExtractorFactory.hpp"
#include "Xpetra_Matrix.hpp"
#include "Xpetra_MatrixFactory.hpp"
#include "Xpetra_BlockedCrsMatrix.hpp"
#include "Xpetra_MatrixMatrix.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraMultiVector.hpp"
#include <Tpetra_RowMatrixTransposer.hpp>
#include <Tpetra_Details_extractBlockDiagonal.hpp>
#include <Tpetra_Details_scaleBlockDiagonal.hpp>
#endif

#include "Xpetra_MatrixUtils_decl.hpp"

namespace Xpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> MatrixUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::xpetraGidNumbering2ThyraGidNumbering(
    const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& input) {
  RCP<const Map> map   = MapUtils::shrinkMapGIDs(*(input.getMap()), *(input.getMap()));
  RCP<MultiVector> ret = MultiVectorFactory::Build(map, input.getNumVectors(), true);
  for (size_t c = 0; c < input.getNumVectors(); c++) {
    Teuchos::ArrayRCP<const Scalar> data = input.getData(c);
    for (size_t r = 0; r < input.getLocalLength(); r++) {
      ret->replaceLocalValue(Teuchos::as<LocalOrdinal>(r), c, data[r]);
    }
  }

  return ret;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> MatrixUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::findColumnSubMap(
    const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& input,
    const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& domainMap) {
  RCP<const Teuchos::Comm<int>> comm = input.getRowMap()->getComm();

  // build an overlapping version of mySpecialMap
  Teuchos::Array<GlobalOrdinal> ovlUnknownStatusGids;
  Teuchos::Array<GlobalOrdinal> ovlFoundStatusGids;

  // loop over global column map of A and find all GIDs where it is not sure, whether they are special or not
  for (size_t i = 0; i < input.getColMap()->getLocalNumElements(); i++) {
    GlobalOrdinal gcid = input.getColMap()->getGlobalElement(i);
    if (domainMap.isNodeGlobalElement(gcid) == false) {
      ovlUnknownStatusGids.push_back(gcid);
    }
  }

  // We need a locally replicated list of all DOF gids of the (non-overlapping) range map of A10
  // Communicate the number of DOFs on each processor
  std::vector<int> myUnknownDofGIDs(comm->getSize(), 0);
  std::vector<int> numUnknownDofGIDs(comm->getSize(), 0);
  myUnknownDofGIDs[comm->getRank()] = ovlUnknownStatusGids.size();
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, comm->getSize(), &myUnknownDofGIDs[0], &numUnknownDofGIDs[0]);

  // create array containing all DOF GIDs
  size_t cntUnknownDofGIDs = 0;
  for (int p = 0; p < comm->getSize(); p++) cntUnknownDofGIDs += numUnknownDofGIDs[p];
  std::vector<GlobalOrdinal> lUnknownDofGIDs(cntUnknownDofGIDs, -1);  // local version to be filled
  std::vector<GlobalOrdinal> gUnknownDofGIDs(cntUnknownDofGIDs, -1);  // global version after communication
  // calculate the offset and fill chunk of memory with local data on each processor
  size_t cntUnknownOffset = 0;
  for (int p = 0; p < comm->getRank(); p++) cntUnknownOffset += numUnknownDofGIDs[p];
  for (size_t k = 0; k < Teuchos::as<size_t>(ovlUnknownStatusGids.size()); k++) {
    lUnknownDofGIDs[k + cntUnknownOffset] = ovlUnknownStatusGids[k];
  }
  if (cntUnknownDofGIDs > 0)
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, Teuchos::as<int>(cntUnknownDofGIDs), &lUnknownDofGIDs[0], &gUnknownDofGIDs[0]);
  std::sort(gUnknownDofGIDs.begin(), gUnknownDofGIDs.end());
  gUnknownDofGIDs.erase(std::unique(gUnknownDofGIDs.begin(), gUnknownDofGIDs.end()), gUnknownDofGIDs.end());

  // loop through all GIDs with unknown status
  for (size_t k = 0; k < gUnknownDofGIDs.size(); k++) {
    GlobalOrdinal curgid = gUnknownDofGIDs[k];
    if (domainMap.isNodeGlobalElement(curgid)) {
      ovlFoundStatusGids.push_back(curgid);  // curgid is in special map (on this processor)
    }
  }

  std::vector<int> myFoundDofGIDs(comm->getSize(), 0);
  std::vector<int> numFoundDofGIDs(comm->getSize(), 0);
  myFoundDofGIDs[comm->getRank()] = ovlFoundStatusGids.size();
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, comm->getSize(), &myFoundDofGIDs[0], &numFoundDofGIDs[0]);

  // create array containing all DOF GIDs
  size_t cntFoundDofGIDs = 0;
  for (int p = 0; p < comm->getSize(); p++) cntFoundDofGIDs += numFoundDofGIDs[p];
  std::vector<GlobalOrdinal> lFoundDofGIDs(cntFoundDofGIDs, -1);  // local version to be filled
  std::vector<GlobalOrdinal> gFoundDofGIDs(cntFoundDofGIDs, -1);  // global version after communication
  // calculate the offset and fill chunk of memory with local data on each processor
  size_t cntFoundOffset = 0;
  for (int p = 0; p < comm->getRank(); p++) cntFoundOffset += numFoundDofGIDs[p];
  for (size_t k = 0; k < Teuchos::as<size_t>(ovlFoundStatusGids.size()); k++) {
    lFoundDofGIDs[k + cntFoundOffset] = ovlFoundStatusGids[k];
  }
  if (cntFoundDofGIDs > 0)
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, Teuchos::as<int>(cntFoundDofGIDs), &lFoundDofGIDs[0], &gFoundDofGIDs[0]);

  Teuchos::Array<GlobalOrdinal> ovlDomainMapArray;
  for (size_t i = 0; i < input.getColMap()->getLocalNumElements(); i++) {
    GlobalOrdinal gcid = input.getColMap()->getGlobalElement(i);
    if (domainMap.isNodeGlobalElement(gcid) == true ||
        std::find(gFoundDofGIDs.begin(), gFoundDofGIDs.end(), gcid) != gFoundDofGIDs.end()) {
      ovlDomainMapArray.push_back(gcid);
    }
  }
  RCP<Map> ovlDomainMap = MapFactory::Build(domainMap.lib(), Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(), ovlDomainMapArray(), 0, comm);
  return ovlDomainMap;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> MatrixUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SplitMatrix(
    const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& input,
    Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>> rangeMapExtractor,
    Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>> domainMapExtractor,
    Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>> columnMapExtractor,
    bool bThyraMode) {
  size_t numRows = rangeMapExtractor->NumMaps();
  size_t numCols = domainMapExtractor->NumMaps();

  TEUCHOS_TEST_FOR_EXCEPTION(rangeMapExtractor->getThyraMode() == true, Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: RangeMapExtractor must not use Thyra style numbering of GIDs. The MapExtractor must contain all GIDs of the full range map in order to define a proper splitting.")
  TEUCHOS_TEST_FOR_EXCEPTION(domainMapExtractor->getThyraMode() == true, Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: DomainMapExtractor must not use Thyra style numbering of GIDs. The MapExtractor must contain all GIDs of the full domain map in order to define a proper splitting.")

  RCP<const Map> fullRangeMap  = rangeMapExtractor->getFullMap();
  RCP<const Map> fullDomainMap = domainMapExtractor->getFullMap();

  TEUCHOS_TEST_FOR_EXCEPTION(fullRangeMap->getMaxAllGlobalIndex() != input.getRowMap()->getMaxAllGlobalIndex(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: RangeMapExtractor incompatible to row map of input matrix.")
  TEUCHOS_TEST_FOR_EXCEPTION(fullRangeMap->getGlobalNumElements() != input.getRowMap()->getGlobalNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: RangeMapExtractor incompatible to row map of input matrix.")
  TEUCHOS_TEST_FOR_EXCEPTION(fullRangeMap->getLocalNumElements() != input.getRowMap()->getLocalNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: RangeMapExtractor incompatible to row map of input matrix.")
  TEUCHOS_TEST_FOR_EXCEPTION(fullRangeMap->getMaxAllGlobalIndex() != input.getRangeMap()->getMaxAllGlobalIndex(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: RangeMapExtractor incompatible to row map of input matrix.")
  TEUCHOS_TEST_FOR_EXCEPTION(fullRangeMap->getGlobalNumElements() != input.getRangeMap()->getGlobalNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: RangeMapExtractor incompatible to row map of input matrix.")
  TEUCHOS_TEST_FOR_EXCEPTION(fullRangeMap->getLocalNumElements() != input.getRangeMap()->getLocalNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: RangeMapExtractor incompatible to row map of input matrix.")

  TEUCHOS_TEST_FOR_EXCEPTION(fullDomainMap->getMaxAllGlobalIndex() != input.getColMap()->getMaxAllGlobalIndex(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: DomainMapExtractor incompatible to domain map of input matrix. fullDomainMap->getMaxAllGlobalIndex() = " << fullDomainMap->getMaxAllGlobalIndex() << " vs. input.getColMap()->getMaxAllGlobalIndex() = " << input.getColMap()->getMaxAllGlobalIndex())
  TEUCHOS_TEST_FOR_EXCEPTION(fullDomainMap->getMaxAllGlobalIndex() != input.getDomainMap()->getMaxAllGlobalIndex(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: DomainMapExtractor incompatible to domain map of input matrix.")
  TEUCHOS_TEST_FOR_EXCEPTION(fullDomainMap->getGlobalNumElements() != input.getDomainMap()->getGlobalNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: DomainMapExtractor incompatible to domain map of input matrix.")
  TEUCHOS_TEST_FOR_EXCEPTION(fullDomainMap->getLocalNumElements() != input.getDomainMap()->getLocalNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: DomainMapExtractor incompatible to domain map of input matrix.")

  // check column map extractor
  Teuchos::RCP<const MapExtractor> myColumnMapExtractor = Teuchos::null;
  if (columnMapExtractor == Teuchos::null) {
    TEUCHOS_TEST_FOR_EXCEPTION(domainMapExtractor->getThyraMode() == true, Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: Auto generation of column map extractor not supported for Thyra style numbering.");
    // This code is always executed, since we always provide map extractors in Xpetra numbering!
    std::vector<Teuchos::RCP<const Map>> ovlxmaps(numCols, Teuchos::null);
    for (size_t c = 0; c < numCols; c++) {
      // TODO: is this routine working correctly?
      Teuchos::RCP<const Map> colMap = MatrixUtils::findColumnSubMap(input, *(domainMapExtractor->getMap(c)));
      ovlxmaps[c]                    = colMap;
    }
    RCP<const Map> fullColMap = MapUtils::concatenateMaps(ovlxmaps);
    // This MapExtractor is always in Xpetra mode!
    myColumnMapExtractor = MapExtractorFactory::Build(fullColMap, ovlxmaps);
  } else
    myColumnMapExtractor = columnMapExtractor;  // use user-provided column map extractor.

  // all above MapExtractors are always in Xpetra mode
  // build separate ones containing Thyra mode GIDs (if necessary)
  Teuchos::RCP<const MapExtractor> thyRangeMapExtractor  = Teuchos::null;
  Teuchos::RCP<const MapExtractor> thyDomainMapExtractor = Teuchos::null;
  Teuchos::RCP<const MapExtractor> thyColMapExtractor    = Teuchos::null;
  if (bThyraMode == true) {
    // build Thyra-style map extractors
    std::vector<Teuchos::RCP<const Map>> thyRgMapExtractorMaps(numRows, Teuchos::null);
    for (size_t r = 0; r < numRows; r++) {
      RCP<const Map> rMap               = rangeMapExtractor->getMap(r);
      RCP<const Map> shrinkedMap        = MapUtils::shrinkMapGIDs(*rMap, *rMap);
      RCP<const StridedMap> strRangeMap = Teuchos::rcp_dynamic_cast<const StridedMap>(rMap);
      if (strRangeMap != Teuchos::null) {
        std::vector<size_t> strInfo          = strRangeMap->getStridingData();
        GlobalOrdinal offset                 = strRangeMap->getOffset();
        LocalOrdinal stridedBlockId          = strRangeMap->getStridedBlockId();
        RCP<const StridedMap> strShrinkedMap = Teuchos::rcp(new StridedMap(shrinkedMap, strInfo, shrinkedMap->getIndexBase(), stridedBlockId, offset));
        thyRgMapExtractorMaps[r]             = strShrinkedMap;
      } else {
        thyRgMapExtractorMaps[r] = shrinkedMap;
      }
      TEUCHOS_TEST_FOR_EXCEPTION(thyRgMapExtractorMaps[r]->getLocalNumElements() != rMap->getLocalNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: Thyra-style range map extractor contains faulty data.")
    }
    RCP<const Map> fullThyRangeMap = MapUtils::concatenateMaps(thyRgMapExtractorMaps);
    thyRangeMapExtractor           = MapExtractorFactory::Build(fullThyRangeMap, thyRgMapExtractorMaps, true);
    std::vector<Teuchos::RCP<const Map>> thyDoMapExtractorMaps(numCols, Teuchos::null);
    std::vector<Teuchos::RCP<const Map>> thyColMapExtractorMaps(numCols, Teuchos::null);
    for (size_t c = 0; c < numCols; c++) {
      RCP<const Map> cMap = domainMapExtractor->getMap(c);

      RCP<const Map> shrinkedDomainMap   = MapUtils::shrinkMapGIDs(*cMap, *cMap);
      RCP<const StridedMap> strDomainMap = Teuchos::rcp_dynamic_cast<const StridedMap>(cMap);
      if (strDomainMap != Teuchos::null) {
        std::vector<size_t> strInfo                = strDomainMap->getStridingData();
        GlobalOrdinal offset                       = strDomainMap->getOffset();
        LocalOrdinal stridedBlockId                = strDomainMap->getStridedBlockId();
        RCP<const StridedMap> strShrinkedDomainMap = Teuchos::rcp(new StridedMap(shrinkedDomainMap, strInfo, shrinkedDomainMap->getIndexBase(), stridedBlockId, offset));
        thyDoMapExtractorMaps[c]                   = strShrinkedDomainMap;
      } else {
        thyDoMapExtractorMaps[c] = shrinkedDomainMap;
      }
      RCP<const Map> colMap           = myColumnMapExtractor->getMap(c);
      RCP<const Map> shrinkedColMap   = MapUtils::shrinkMapGIDs(*colMap, *cMap);
      RCP<const StridedMap> strColMap = Teuchos::rcp_dynamic_cast<const StridedMap>(colMap);
      if (strColMap != Teuchos::null) {
        std::vector<size_t> strInfo             = strColMap->getStridingData();
        GlobalOrdinal offset                    = strColMap->getOffset();
        LocalOrdinal stridedBlockId             = strColMap->getStridedBlockId();
        RCP<const StridedMap> strShrinkedColMap = Teuchos::rcp(new StridedMap(shrinkedColMap, strInfo, shrinkedColMap->getIndexBase(), stridedBlockId, offset));
        thyColMapExtractorMaps[c]               = strShrinkedColMap;
      } else {
        thyColMapExtractorMaps[c] = shrinkedColMap;
      }

      TEUCHOS_TEST_FOR_EXCEPTION(thyColMapExtractorMaps[c]->getLocalNumElements() != colMap->getLocalNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: Thyra-style column map extractor contains faulty data.")
      TEUCHOS_TEST_FOR_EXCEPTION(thyDoMapExtractorMaps[c]->getLocalNumElements() != cMap->getLocalNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: Thyra-style domain map extractor contains faulty data.")
    }
    RCP<const Map> fullThyDomainMap = MapUtils::concatenateMaps(thyDoMapExtractorMaps);
    RCP<const Map> fullThyColumnMap = MapUtils::concatenateMaps(thyColMapExtractorMaps);
    thyDomainMapExtractor           = MapExtractorFactory::Build(fullThyDomainMap, thyDoMapExtractorMaps, true);
    thyColMapExtractor              = MapExtractorFactory::Build(fullThyColumnMap, thyColMapExtractorMaps, true);
  }
  // create submatrices
  std::vector<Teuchos::RCP<Matrix>> subMatrices(numRows * numCols, Teuchos::null);
  for (size_t r = 0; r < numRows; r++) {
    for (size_t c = 0; c < numCols; c++) {
      // create empty CrsMatrix objects
      // make sure that the submatrices are defined using the right row maps (either Thyra or xpetra style)
      // Note: we're reserving a little bit too much memory for the submatrices, but should be still reasonable
      if (bThyraMode == true)
        subMatrices[r * numCols + c] = MatrixFactory::Build(thyRangeMapExtractor->getMap(r, true), input.getLocalMaxNumRowEntries());
      else
        subMatrices[r * numCols + c] = MatrixFactory::Build(rangeMapExtractor->getMap(r), input.getLocalMaxNumRowEntries());
    }
  }

  // We need a vector which lives on the column map of input and stores the block id that the column belongs to.
  // create a vector on the domain map. Loop over it and fill in the corresponding block id numbers
  // create a vector on the column map and import the data
  // Importer: source map is non-overlapping. Target map is overlapping
  // call colMap.Import(domMap,Importer,Insert)
  // do the same with "Add" to make sure only one processor is responsible for the different GIDs!
#if 0  // TAW needs to be fixed (does not compile for Scalar=complex)
    typedef Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>  VectorFactory;
    typedef Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>  Vector;
    typedef Xpetra::ImportFactory<LocalOrdinal, GlobalOrdinal, Node> ImportFactory;
    typedef Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> Import;

    RCP<Vector> doCheck = VectorFactory::Build(input.getDomainMap(), true);
    doCheck->putScalar(1.0);
    RCP<Vector> coCheck = VectorFactory::Build(input.getColMap(), true);
    RCP<Import> imp = ImportFactory::Build(input.getDomainMap(), input.getColMap());
    coCheck->doImport(*doCheck, *imp, Xpetra::ADD);
    TEUCHOS_TEST_FOR_EXCEPTION(coCheck->normInf() != Teuchos::ScalarTraits< Scalar >::magnitude(1.0), Xpetra::Exceptions::RuntimeError, "Xpetra::MatrixUtils::SplitMatrix: error when distributing data.");

    doCheck->putScalar(-1.0);
    coCheck->putScalar(-1.0);

    Teuchos::ArrayRCP< Scalar > doCheckData = doCheck->getDataNonConst(0);
    for (size_t rrr = 0; rrr < input.getDomainMap()->getLocalNumElements(); rrr++) {
      // global row id to extract data from global monolithic matrix
      GlobalOrdinal id = input.getDomainMap()->getGlobalElement(rrr); // LID -> GID (column)

      // Find the block id in range map extractor that belongs to same global id.
      size_t BlockId = domainMapExtractor->getMapIndexForGID(id);

      doCheckData[rrr] = Teuchos::as<Scalar>(BlockId);
    }

    coCheck->doImport(*doCheck, *imp, Xpetra::INSERT);

    Teuchos::ArrayRCP< Scalar > coCheckData = coCheck->getDataNonConst(0);
#endif
  // loop over all rows of input matrix
  for (size_t rr = 0; rr < input.getRowMap()->getLocalNumElements(); rr++) {
    // global row id to extract data from global monolithic matrix
    GlobalOrdinal growid = input.getRowMap()->getGlobalElement(rr);  // LID -> GID (column)

    // Find the block id in range map extractor that belongs to same global row id.
    // We assume the global ids to be unique in the map extractor
    // The MapExtractor objects may be constructed using the thyra mode. However, we use
    // global GID ids (as we have a single input matrix). The only problematic thing could
    // be that the GIDs of the input matrix are inconsistent with the GIDs in the map extractor.
    // Note, that for the Thyra mode the GIDs in the map extractor are generated using the ordering
    // of the blocks.
    size_t rowBlockId = rangeMapExtractor->getMapIndexForGID(growid);

    // global row id used for subblocks to insert information
    GlobalOrdinal subblock_growid = growid;  // for Xpetra-style numbering the global row ids are not changing
    if (bThyraMode == true) {
      // find local row id associated with growid in the corresponding subblock
      LocalOrdinal lrowid = rangeMapExtractor->getMap(rowBlockId)->getLocalElement(growid);
      // translate back local row id to global row id for the subblock
      subblock_growid = thyRangeMapExtractor->getMap(rowBlockId, true)->getGlobalElement(lrowid);
    }

    // extract matrix entries from input matrix
    // we use global ids since we have to determine the corresponding
    // block column ids using the global ids anyway
    Teuchos::ArrayView<const LocalOrdinal> indices;
    Teuchos::ArrayView<const Scalar> vals;
    input.getLocalRowView(rr, indices, vals);

    std::vector<Teuchos::Array<GlobalOrdinal>> blockColIdx(numCols, Teuchos::Array<GlobalOrdinal>());
    std::vector<Teuchos::Array<Scalar>> blockColVals(numCols, Teuchos::Array<Scalar>());

    for (size_t i = 0; i < (size_t)indices.size(); i++) {
      // gobal column id to extract data from full monolithic matrix
      GlobalOrdinal gcolid = input.getColMap()->getGlobalElement(indices[i]);

      size_t colBlockId = myColumnMapExtractor->getMapIndexForGID(gcolid);  // old buggy thing
      // size_t colBlockId = Teuchos::as<size_t>(coCheckData[indices[i]]);

      // global column id used for subblocks to insert information
      GlobalOrdinal subblock_gcolid = gcolid;  // for Xpetra-style numbering the global col ids are not changing
      if (bThyraMode == true) {
        // find local col id associated with gcolid in the corresponding subblock
        LocalOrdinal lcolid = myColumnMapExtractor->getMap(colBlockId)->getLocalElement(gcolid);
        // translate back local col id to global col id for the subblock
        subblock_gcolid = thyColMapExtractor->getMap(colBlockId, true)->getGlobalElement(lcolid);
      }
      blockColIdx[colBlockId].push_back(subblock_gcolid);
      blockColVals[colBlockId].push_back(vals[i]);
    }

    for (size_t c = 0; c < numCols; c++) {
      subMatrices[rowBlockId * numCols + c]->insertGlobalValues(subblock_growid, blockColIdx[c].view(0, blockColIdx[c].size()), blockColVals[c].view(0, blockColVals[c].size()));
    }
  }

  // call fill complete on subblocks and create BlockedCrsOperator
  RCP<BlockedCrsMatrix> bA = Teuchos::null;
  if (bThyraMode == true) {
    for (size_t r = 0; r < numRows; r++) {
      for (size_t c = 0; c < numCols; c++) {
        subMatrices[r * numCols + c]->fillComplete(thyDomainMapExtractor->getMap(c, true), thyRangeMapExtractor->getMap(r, true));
      }
    }
    bA = Teuchos::rcp(new BlockedCrsMatrix(thyRangeMapExtractor, thyDomainMapExtractor, 10 /*input.getRowMap()->getLocalMaxNumRowEntries()*/));
  } else {
    for (size_t r = 0; r < numRows; r++) {
      for (size_t c = 0; c < numCols; c++) {
        subMatrices[r * numCols + c]->fillComplete(domainMapExtractor->getMap(c), rangeMapExtractor->getMap(r));
      }
    }
    bA = Teuchos::rcp(new BlockedCrsMatrix(rangeMapExtractor, domainMapExtractor, 10 /*input.getRowMap()->getLocalMaxNumRowEntries()*/));
  }

  for (size_t r = 0; r < numRows; r++) {
    for (size_t c = 0; c < numCols; c++) {
      bA->setMatrix(r, c, subMatrices[r * numCols + c]);
    }
  }
  return bA;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MatrixUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CheckRepairMainDiagonal(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& Ac,
                                                                                     bool const& repairZeroDiagonals, Teuchos::FancyOStream& fos,
                                                                                     const typename Teuchos::ScalarTraits<Scalar>::magnitudeType threshold,
                                                                                     const Scalar replacementValue) {
  using TST = typename Teuchos::ScalarTraits<Scalar>;
  using Teuchos::rcp_dynamic_cast;

  GlobalOrdinal gZeroDiags;
  bool usedEfficientPath = false;

#ifdef HAVE_MUELU_TPETRA
  RCP<CrsMatrixWrap> crsWrapAc = rcp_dynamic_cast<CrsMatrixWrap>(Ac);
  RCP<TpetraCrsMatrix> tpCrsAc;
  if (!crsWrapAc.is_null())
    tpCrsAc = rcp_dynamic_cast<TpetraCrsMatrix>(crsWrapAc->getCrsMatrix());

  if (!tpCrsAc.is_null()) {
    auto tpCrsGraph = tpCrsAc->getTpetra_CrsMatrix()->getCrsGraph();
    size_t numRows  = Ac->getRowMap()->getLocalNumElements();
    typedef typename Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>::offset_device_view_type offset_type;
    using range_type = Kokkos::RangePolicy<LocalOrdinal, typename Node::execution_space>;
    auto offsets     = offset_type(Kokkos::ViewAllocateWithoutInitializing("offsets"), numRows);
    tpCrsGraph->getLocalDiagOffsets(offsets);

    const size_t STINV =
        Tpetra::Details::OrdinalTraits<typename offset_type::value_type>::invalid();

    if (repairZeroDiagonals) {
      // Make sure that the matrix has all its diagonal entries, so
      // we can fix them in-place.

      LO numMissingDiagonalEntries = 0;

      Kokkos::parallel_reduce(
          "countMissingDiagonalEntries",
          range_type(0, numRows),
          KOKKOS_LAMBDA(const LO i, LO& missing) {
            missing += (offsets(i) == STINV);
          },
          numMissingDiagonalEntries);

      GlobalOrdinal gNumMissingDiagonalEntries;
      Teuchos::reduceAll(*(Ac->getRowMap()->getComm()), Teuchos::REDUCE_SUM, Teuchos::as<GlobalOrdinal>(numMissingDiagonalEntries),
                         Teuchos::outArg(gNumMissingDiagonalEntries));

      if (gNumMissingDiagonalEntries == 0) {
        // Matrix has all diagonal entries, now we fix them

        auto lclA = tpCrsAc->getTpetra_CrsMatrix()->getLocalMatrixDevice();

        using ATS      = Kokkos::ArithTraits<Scalar>;
        using impl_ATS = Kokkos::ArithTraits<typename ATS::val_type>;

        LO lZeroDiags                                = 0;
        typename ATS::val_type impl_replacementValue = replacementValue;

        Kokkos::parallel_reduce(
            "fixSmallDiagonalEntries",
            range_type(0, numRows),
            KOKKOS_LAMBDA(const LO i, LO& fixed) {
              const auto offset = offsets(i);
              auto curRow       = lclA.row(i);
              if (impl_ATS::magnitude(curRow.value(offset)) <= threshold) {
                curRow.value(offset) = impl_replacementValue;
                fixed += 1;
              }
            },
            lZeroDiags);

        Teuchos::reduceAll(*(Ac->getRowMap()->getComm()), Teuchos::REDUCE_SUM, Teuchos::as<GlobalOrdinal>(lZeroDiags),
                           Teuchos::outArg(gZeroDiags));

        usedEfficientPath = true;
      }
    } else {
      // We only want to count up small diagonal entries, but not
      // fix them. So missing diagonal entries are not an issue.

      auto lclA = tpCrsAc->getTpetra_CrsMatrix()->getLocalMatrixDevice();

      using ATS      = Kokkos::ArithTraits<Scalar>;
      using impl_ATS = Kokkos::ArithTraits<typename ATS::val_type>;

      LO lZeroDiags = 0;

      Kokkos::parallel_reduce(
          "detectSmallDiagonalEntries",
          range_type(0, numRows),
          KOKKOS_LAMBDA(const LO i, LO& small) {
            const auto offset = offsets(i);
            if (offset == STINV)
              small += 1;
            else {
              auto curRow = lclA.row(i);
              if (impl_ATS::magnitude(curRow.value(offset)) <= threshold) {
                small += 1;
              }
            }
          },
          lZeroDiags);

      Teuchos::reduceAll(*(Ac->getRowMap()->getComm()), Teuchos::REDUCE_SUM, Teuchos::as<GlobalOrdinal>(lZeroDiags),
                         Teuchos::outArg(gZeroDiags));

      usedEfficientPath = true;
    }
  }
#endif

  if (!usedEfficientPath) {
    RCP<const Map> rowMap = Ac->getRowMap();
    RCP<Vector> diagVec   = VectorFactory::Build(rowMap);
    Ac->getLocalDiagCopy(*diagVec);

    LocalOrdinal lZeroDiags                 = 0;
    Teuchos::ArrayRCP<const Scalar> diagVal = diagVec->getData(0);

    for (size_t i = 0; i < rowMap->getLocalNumElements(); i++) {
      if (TST::magnitude(diagVal[i]) <= threshold) {
        lZeroDiags++;
      }
    }

    Teuchos::reduceAll(*(rowMap->getComm()), Teuchos::REDUCE_SUM, Teuchos::as<GlobalOrdinal>(lZeroDiags),
                       Teuchos::outArg(gZeroDiags));

    if (repairZeroDiagonals && gZeroDiags > 0) {
      /*
        TAW: If Ac has empty rows, put a 1 on the diagonal of Ac. Be aware that Ac might have empty rows AND
        columns.  The columns might not exist in the column map at all.  It would be nice to add the entries
        to the original matrix Ac. But then we would have to use insertLocalValues. However we cannot add
        new entries for local column indices that do not exist in the column map of Ac (at least Epetra is
        not able to do this).  With Tpetra it is also not possible to add new entries after the FillComplete
        call with a static map, since the column map already exists and the diagonal entries are completely
        missing.  Here we build a diagonal matrix with zeros on the diagonal and ones on the diagonal for
        the rows where Ac has empty rows We have to build a new matrix to be able to use insertGlobalValues.
        Then we add the original matrix Ac to our new block diagonal matrix and use the result as new
        (non-singular) matrix Ac.  This is very inefficient.

        If you know something better, please let me know.
      */
      RCP<Matrix> fixDiagMatrix = MatrixFactory::Build(rowMap, 1);
      Teuchos::Array<GlobalOrdinal> indout(1);
      Teuchos::Array<Scalar> valout(1);
      for (size_t r = 0; r < rowMap->getLocalNumElements(); r++) {
        if (TST::magnitude(diagVal[r]) <= threshold) {
          GlobalOrdinal grid = rowMap->getGlobalElement(r);
          indout[0]          = grid;
          valout[0]          = replacementValue;
          fixDiagMatrix->insertGlobalValues(grid, indout(), valout());
        }
      }
      {
        Teuchos::TimeMonitor m1(*Teuchos::TimeMonitor::getNewTimer("CheckRepairMainDiagonal: fillComplete1"));
        fixDiagMatrix->fillComplete(Ac->getDomainMap(), Ac->getRangeMap());
      }

      RCP<Matrix> newAc;
      Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TwoMatrixAdd(*Ac, false, 1.0, *fixDiagMatrix, false, 1.0, newAc, fos);
      if (Ac->IsView("stridedMaps"))
        newAc->CreateView("stridedMaps", Ac);

      Ac            = Teuchos::null;  // free singular matrix
      fixDiagMatrix = Teuchos::null;
      Ac            = newAc;  // set fixed non-singular matrix

      // call fillComplete with optimized storage option set to true
      // This is necessary for new faster Epetra MM kernels.
      Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList());
      p->set("DoOptimizeStorage", true);
      {
        Teuchos::TimeMonitor m1(*Teuchos::TimeMonitor::getNewTimer("CheckRepairMainDiagonal: fillComplete2"));
        Ac->fillComplete(p);
      }
    }  // end repair
  }

  // print some output
  fos << "CheckRepairMainDiagonal: " << (repairZeroDiagonals ? "repaired " : "found ")
      << gZeroDiags << " too small entries (threshold = " << threshold << ") on main diagonal of Ac." << std::endl;

#ifdef HAVE_XPETRA_DEBUG  // only for debugging
  {
    // check whether Ac has been repaired...
    RCP<const Map> rowMap = Ac->getRowMap();
    RCP<Vector> diagVec   = VectorFactory::Build(rowMap);
    Teuchos::ArrayRCP<const Scalar> diagVal;
    Ac->getLocalDiagCopy(*diagVec);
    diagVal = diagVec->getData(0);
    for (size_t r = 0; r < Ac->getRowMap()->getLocalNumElements(); r++) {
      if (TST::magnitude(diagVal[r]) <= threshold) {
        fos << "Error: there are too small entries left on diagonal after repair..." << std::endl;
        break;
      }
    }
  }
#endif
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MatrixUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::RelativeDiagonalBoost(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& A,
                                                                                   const Teuchos::ArrayView<const double>& relativeThreshold, Teuchos::FancyOStream& fos) {
  Teuchos::TimeMonitor m1(*Teuchos::TimeMonitor::getNewTimer("RelativeDiagonalBoost"));

  TEUCHOS_TEST_FOR_EXCEPTION(A->GetFixedBlockSize() != relativeThreshold.size() && relativeThreshold.size() != 1, Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::RelativeDiagonal Boost:  Either A->GetFixedBlockSize() != relativeThreshold.size() OR relativeThreshold.size() == 1");

  LocalOrdinal numPDEs = A->GetFixedBlockSize();
  typedef typename Teuchos::ScalarTraits<Scalar> TST;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MT;
  Scalar zero = TST::zero();
  Scalar one  = TST::one();

  // Get the diagonal
  RCP<Vector> diag = VectorFactory::Build(A->getRowMap());
  A->getLocalDiagCopy(*diag);
  Teuchos::ArrayRCP<const Scalar> dataVal = diag->getData(0);
  size_t N                                = A->getRowMap()->getLocalNumElements();

  // Compute the diagonal maxes for each PDE
  std::vector<MT> l_diagMax(numPDEs), g_diagMax(numPDEs);
  for (size_t i = 0; i < N; i++) {
    int pde = (int)(i % numPDEs);
    if ((int)i < numPDEs)
      l_diagMax[pde] = TST::magnitude(dataVal[i]);
    else
      l_diagMax[pde] = std::max(l_diagMax[pde], TST::magnitude(dataVal[i]));
  }
  Teuchos::reduceAll(*A->getRowMap()->getComm(), Teuchos::REDUCE_MAX, numPDEs, l_diagMax.data(), g_diagMax.data());

  // Apply the diagonal maxes via matrix-matrix addition
  RCP<Matrix> boostMatrix = MatrixFactory::Build(A->getRowMap(), 1);
  Teuchos::Array<GlobalOrdinal> index(1);
  Teuchos::Array<Scalar> value(1);
  for (size_t i = 0; i < N; i++) {
    GlobalOrdinal GRID = A->getRowMap()->getGlobalElement(i);
    int pde            = (int)(i % numPDEs);
    index[0]           = GRID;
    if (TST::magnitude(dataVal[i]) < relativeThreshold[pde] * g_diagMax[pde])
      value[0] = relativeThreshold[pde] * g_diagMax[pde] - TST::magnitude(dataVal[i]);
    else
      value[0] = zero;
    boostMatrix->insertGlobalValues(GRID, index(), value());
  }
  boostMatrix->fillComplete(A->getDomainMap(), A->getRangeMap());

  // FIXME: We really need an add that lets you "add into"
  RCP<Matrix> newA;
  Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TwoMatrixAdd(*A, false, one, *boostMatrix, false, one, newA, fos);
  if (A->IsView("stridedMaps"))
    newA->CreateView("stridedMaps", A);
  A = newA;
  A->fillComplete();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MatrixUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::extractBlockDiagonal(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                                                                                  Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& diagonal) {
  const UnderlyingLib lib = A.getRowMap()->lib();

  if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_XPETRA_EPETRA)
    throw(Xpetra::Exceptions::RuntimeError("Xpetra::MatrixUtils::extractBlockDiagonal not available for Epetra."));
#endif  // HAVE_XPETRA_EPETRA
  } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
    const Tpetra::CrsMatrix<SC, LO, GO, NO>& At = Xpetra::Helpers<SC, LO, GO, NO>::Op2TpetraCrs(A);
    Tpetra::MultiVector<SC, LO, GO, NO>& Dt     = Xpetra::toTpetra(diagonal);
    Tpetra::Details::extractBlockDiagonal(At, Dt);
#endif  // HAVE_XPETRA_TPETRA
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MatrixUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::inverseScaleBlockDiagonal(Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& blockDiagonal,
                                                                                       bool doTranspose,
                                                                                       Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& toBeScaled) {
  const UnderlyingLib lib = blockDiagonal.getMap()->lib();

  if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_XPETRA_EPETRA)
    throw(Xpetra::Exceptions::RuntimeError("Xpetra::MatrixUtils::inverseScaleBlockDiagonal not available for Epetra."));
#endif  // HAVE_XPETRA_EPETRA
  } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_XPETRA_TPETRA
    Tpetra::MultiVector<SC, LO, GO, NO>& Dt = Xpetra::toTpetra(blockDiagonal);
    Tpetra::MultiVector<SC, LO, GO, NO>& St = Xpetra::toTpetra(toBeScaled);
    Tpetra::Details::inverseScaleBlockDiagonal(Dt, doTranspose, St);
#endif  // HAVE_XPETRA_TPETRA
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MatrixUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::checkLocalRowMapMatchesColMap(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A) {
  RCP<const Map> rowMap = A.getRowMap();
  RCP<const Map> colMap = A.getColMap();
  bool fail             = false;
  if (rowMap->lib() == Xpetra::UseTpetra) {
    auto tpRowMap = Teuchos::rcp_dynamic_cast<const TpetraMap>(rowMap, true)->getTpetra_Map();
    auto tpColMap = Teuchos::rcp_dynamic_cast<const TpetraMap>(colMap, true)->getTpetra_Map();
    fail          = !tpColMap->isLocallyFitted(*tpRowMap);
  } else {
    RCP<const Teuchos::Comm<int>> comm = rowMap->getComm();
    LO numRows                         = Teuchos::as<LocalOrdinal>(rowMap->getLocalNumElements());

    for (LO rowLID = 0; rowLID < numRows; rowLID++) {
      GO rowGID = rowMap->getGlobalElement(rowLID);
      LO colGID = colMap->getGlobalElement(rowLID);
      if (rowGID != colGID) {
        fail = true;
        std::cerr << "On rank " << comm->getRank() << ", LID " << rowLID << " is GID " << rowGID << " in the rowmap, but GID " << colGID << " in the column map.\n";
      }
    }
  }
  TEUCHOS_TEST_FOR_EXCEPTION(fail, Exceptions::RuntimeError,
                             "Local parts of row and column map do not match!");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MatrixUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::convertMatrixToStridedMaps(
    Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> matrix,
    std::vector<size_t>& rangeStridingInfo, std::vector<size_t>& domainStridingInfo) {
  RCP<const StridedMap> stridedRowMap = StridedMapFactory::Build(matrix->getRowMap(), rangeStridingInfo, -1, 0);
  RCP<const StridedMap> stridedColMap = StridedMapFactory::Build(matrix->getColMap(), domainStridingInfo, -1, 0);

  if (matrix->IsView("stridedMaps") == true) matrix->RemoveView("stridedMaps");
  matrix->CreateView("stridedMaps", stridedRowMap, stridedColMap);
}

}  // end namespace Xpetra

#endif  // PACKAGES_XPETRA_SUP_MATRIX_UTILS_DEF_HPP_
