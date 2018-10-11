// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef PACKAGES_XPETRA_SUP_MATRIX_UTILS_HPP_
#define PACKAGES_XPETRA_SUP_MATRIX_UTILS_HPP_

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

namespace Xpetra {

/*!
  @class MatrixUtils
  @brief Xpetra utility class for common matrix-related routines

  The routines should be independent from Epetra/Tpetra and be purely implemented in Xpetra.
  Other matrix-related routines are out-sourced into other helper classes (e.g. MatrixMatrix for
  MM multiplication and addition).

*/
template <class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
class MatrixUtils {
#undef XPETRA_MATRIXUTILS_SHORT
#include "Xpetra_UseShortNames.hpp"

public:

  static Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > xpetraGidNumbering2ThyraGidNumbering(
      const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& input) {
    typedef Xpetra::MapUtils<LocalOrdinal, GlobalOrdinal, Node>  MapUtils;
    typedef Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>  MultiVector;
    typedef Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>  MultiVectorFactory;
    RCP<const Map> map = MapUtils::shrinkMapGIDs(*(input.getMap()),*(input.getMap()));
    RCP<MultiVector> ret = MultiVectorFactory::Build(map, input.getNumVectors(), true);
    for (size_t c = 0; c < input.getNumVectors(); c++) {
      Teuchos::ArrayRCP< const Scalar > data = input.getData(c);
      for (size_t r = 0; r < input.getLocalLength(); r++) {
        ret->replaceLocalValue(Teuchos::as<LocalOrdinal>(r), c, data[r]);
      }
    }

    return ret;
  }

  static Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > findColumnSubMap(
      const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& input,
      const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& domainMap ) {

    RCP< const Teuchos::Comm<int> > comm = input.getRowMap()->getComm();

    // build an overlapping version of mySpecialMap
    Teuchos::Array<GlobalOrdinal> ovlUnknownStatusGids;
    Teuchos::Array<GlobalOrdinal> ovlFoundStatusGids;

    // loop over global column map of A and find all GIDs where it is not sure, whether they are special or not
    for(size_t i = 0; i<input.getColMap()->getNodeNumElements(); i++) {
      GlobalOrdinal gcid = input.getColMap()->getGlobalElement(i);
      if(domainMap.isNodeGlobalElement(gcid) == false) {
        ovlUnknownStatusGids.push_back(gcid);
      }
    }

    // We need a locally replicated list of all DOF gids of the (non-overlapping) range map of A10
    // Communicate the number of DOFs on each processor
    std::vector<int> myUnknownDofGIDs(comm->getSize(),0);
    std::vector<int> numUnknownDofGIDs(comm->getSize(),0);
    myUnknownDofGIDs[comm->getRank()] = ovlUnknownStatusGids.size();
    Teuchos::reduceAll(*comm,Teuchos::REDUCE_MAX,comm->getSize(),&myUnknownDofGIDs[0],&numUnknownDofGIDs[0]);

    // create array containing all DOF GIDs
    size_t cntUnknownDofGIDs = 0;
    for(int p = 0; p < comm->getSize(); p++) cntUnknownDofGIDs += numUnknownDofGIDs[p];
    std::vector<GlobalOrdinal> lUnknownDofGIDs(cntUnknownDofGIDs,-1); // local version to be filled
    std::vector<GlobalOrdinal> gUnknownDofGIDs(cntUnknownDofGIDs,-1); // global version after communication
    // calculate the offset and fill chunk of memory with local data on each processor
    size_t cntUnknownOffset = 0;
    for(int p = 0; p < comm->getRank(); p++) cntUnknownOffset += numUnknownDofGIDs[p];
    for(size_t k=0; k < Teuchos::as<size_t>(ovlUnknownStatusGids.size()); k++) {
      lUnknownDofGIDs[k+cntUnknownOffset] = ovlUnknownStatusGids[k];
    }
    if(cntUnknownDofGIDs > 0)
      Teuchos::reduceAll(*comm,Teuchos::REDUCE_MAX,Teuchos::as<int>(cntUnknownDofGIDs),&lUnknownDofGIDs[0],&gUnknownDofGIDs[0]);
    std::sort(gUnknownDofGIDs.begin(), gUnknownDofGIDs.end());
    gUnknownDofGIDs.erase(std::unique(gUnknownDofGIDs.begin(), gUnknownDofGIDs.end()), gUnknownDofGIDs.end());

    // loop through all GIDs with unknown status
    for(size_t k=0; k < gUnknownDofGIDs.size(); k++) {
      GlobalOrdinal curgid = gUnknownDofGIDs[k];
      if(domainMap.isNodeGlobalElement(curgid)) {
        ovlFoundStatusGids.push_back(curgid); // curgid is in special map (on this processor)
      }
    }

    std::vector<int> myFoundDofGIDs(comm->getSize(),0);
    std::vector<int> numFoundDofGIDs(comm->getSize(),0);
    myFoundDofGIDs[comm->getRank()] = ovlFoundStatusGids.size();
    Teuchos::reduceAll(*comm,Teuchos::REDUCE_MAX,comm->getSize(),&myFoundDofGIDs[0],&numFoundDofGIDs[0]);

    // create array containing all DOF GIDs
    size_t cntFoundDofGIDs = 0;
    for(int p = 0; p < comm->getSize(); p++) cntFoundDofGIDs += numFoundDofGIDs[p];
    std::vector<GlobalOrdinal> lFoundDofGIDs(cntFoundDofGIDs,-1); // local version to be filled
    std::vector<GlobalOrdinal> gFoundDofGIDs(cntFoundDofGIDs,-1); // global version after communication
    // calculate the offset and fill chunk of memory with local data on each processor
    size_t cntFoundOffset = 0;
    for(int p = 0; p < comm->getRank(); p++) cntFoundOffset += numFoundDofGIDs[p];
    for(size_t k=0; k < Teuchos::as<size_t>(ovlFoundStatusGids.size()); k++) {
      lFoundDofGIDs[k+cntFoundOffset] = ovlFoundStatusGids[k];
    }
    if(cntFoundDofGIDs > 0)
      Teuchos::reduceAll(*comm,Teuchos::REDUCE_MAX,Teuchos::as<int>(cntFoundDofGIDs),&lFoundDofGIDs[0],&gFoundDofGIDs[0]);

    Teuchos::Array<GlobalOrdinal> ovlDomainMapArray;
    for(size_t i = 0; i<input.getColMap()->getNodeNumElements(); i++) {
      GlobalOrdinal gcid = input.getColMap()->getGlobalElement(i);
      if(domainMap.isNodeGlobalElement(gcid) == true ||
         std::find(gFoundDofGIDs.begin(), gFoundDofGIDs.end(), gcid) != gFoundDofGIDs.end()) {
        ovlDomainMapArray.push_back(gcid);
      }
    }
    RCP<Map> ovlDomainMap = MapFactory::Build (domainMap.lib(),Teuchos::OrdinalTraits<GlobalOrdinal>::invalid(),ovlDomainMapArray(),0,comm);
    return ovlDomainMap;
  }

  /** Given a matrix A split it into a nxm blocked matrix using the map extractors.

    @param input Input matrix, must already have had 'FillComplete()' called.
    @param rangeMapExtractor MapExtractor object describing the splitting of rows of the output block matrix
    @param domainMapExtractor MapExtractor object describing the splitting of columns of the output block matrix
    @param columnMapExtractor (not fully clear whether we need that. is always Teuchos::null)
    @param bThyraMode If true, build a n x n blocked operator using Thyra GIDs
  */
  static Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > SplitMatrix(
                       const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& input,
                       Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rangeMapExtractor,
                       Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > domainMapExtractor,
                       Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > columnMapExtractor = Teuchos::null,
                       bool bThyraMode = false) {
    typedef Xpetra::MapUtils<LocalOrdinal, GlobalOrdinal, Node>  MapUtils;
    typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractor;
    typedef Xpetra::MapExtractorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractorFactory;
    typedef Xpetra::MatrixUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>  MatrixUtils;

    size_t numRows  = rangeMapExtractor->NumMaps();
    size_t numCols  = domainMapExtractor->NumMaps();

    TEUCHOS_TEST_FOR_EXCEPTION(rangeMapExtractor->getThyraMode() == true, Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: RangeMapExtractor must not use Thyra style numbering of GIDs. The MapExtractor must contain all GIDs of the full range map in order to define a proper splitting.")
    TEUCHOS_TEST_FOR_EXCEPTION(domainMapExtractor->getThyraMode() == true, Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: DomainMapExtractor must not use Thyra style numbering of GIDs. The MapExtractor must contain all GIDs of the full domain map in order to define a proper splitting.")

    RCP<const Map> fullRangeMap  = rangeMapExtractor->getFullMap();
    RCP<const Map> fullDomainMap = domainMapExtractor->getFullMap();

    TEUCHOS_TEST_FOR_EXCEPTION(fullRangeMap->getMaxAllGlobalIndex() != input.getRowMap()->getMaxAllGlobalIndex(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: RangeMapExtractor incompatible to row map of input matrix.")
    TEUCHOS_TEST_FOR_EXCEPTION(fullRangeMap->getGlobalNumElements() != input.getRowMap()->getGlobalNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: RangeMapExtractor incompatible to row map of input matrix.")
    TEUCHOS_TEST_FOR_EXCEPTION(fullRangeMap->getNodeNumElements()   != input.getRowMap()->getNodeNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: RangeMapExtractor incompatible to row map of input matrix.")
    TEUCHOS_TEST_FOR_EXCEPTION(fullRangeMap->getMaxAllGlobalIndex() != input.getRangeMap()->getMaxAllGlobalIndex(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: RangeMapExtractor incompatible to row map of input matrix.")
    TEUCHOS_TEST_FOR_EXCEPTION(fullRangeMap->getGlobalNumElements() != input.getRangeMap()->getGlobalNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: RangeMapExtractor incompatible to row map of input matrix.")
    TEUCHOS_TEST_FOR_EXCEPTION(fullRangeMap->getNodeNumElements()   != input.getRangeMap()->getNodeNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: RangeMapExtractor incompatible to row map of input matrix.")

    TEUCHOS_TEST_FOR_EXCEPTION(fullDomainMap->getMaxAllGlobalIndex() != input.getColMap()->getMaxAllGlobalIndex(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: DomainMapExtractor incompatible to domain map of input matrix. fullDomainMap->getMaxAllGlobalIndex() = " << fullDomainMap->getMaxAllGlobalIndex() << " vs. input.getColMap()->getMaxAllGlobalIndex() = " << input.getColMap()->getMaxAllGlobalIndex())
    TEUCHOS_TEST_FOR_EXCEPTION(fullDomainMap->getMaxAllGlobalIndex() != input.getDomainMap()->getMaxAllGlobalIndex(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: DomainMapExtractor incompatible to domain map of input matrix.")
    TEUCHOS_TEST_FOR_EXCEPTION(fullDomainMap->getGlobalNumElements() != input.getDomainMap()->getGlobalNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: DomainMapExtractor incompatible to domain map of input matrix.")
    TEUCHOS_TEST_FOR_EXCEPTION(fullDomainMap->getNodeNumElements()   != input.getDomainMap()->getNodeNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: DomainMapExtractor incompatible to domain map of input matrix.")

    // check column map extractor
    Teuchos::RCP<const MapExtractor> myColumnMapExtractor = Teuchos::null;
    if(columnMapExtractor == Teuchos::null) {
      TEUCHOS_TEST_FOR_EXCEPTION(domainMapExtractor->getThyraMode() == true, Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: Auto generation of column map extractor not supported for Thyra style numbering.");
      // This code is always executed, since we always provide map extractors in Xpetra numbering!
      std::vector<Teuchos::RCP<const Map> > ovlxmaps(numCols, Teuchos::null);
      for(size_t c = 0; c < numCols; c++) {
        // TODO: is this routine working correctly?
        Teuchos::RCP<const Map> colMap = MatrixUtils::findColumnSubMap(input, *(domainMapExtractor->getMap(c)));
        ovlxmaps[c] = colMap;
      }
      RCP<const Map> fullColMap = MapUtils::concatenateMaps(ovlxmaps);
      // This MapExtractor is always in Xpetra mode!
      myColumnMapExtractor = MapExtractorFactory::Build(fullColMap,ovlxmaps);
    } else
      myColumnMapExtractor = columnMapExtractor; // use user-provided column map extractor.

    // all above MapExtractors are always in Xpetra mode
    // build separate ones containing Thyra mode GIDs (if necessary)
    Teuchos::RCP<const MapExtractor> thyRangeMapExtractor  = Teuchos::null;
    Teuchos::RCP<const MapExtractor> thyDomainMapExtractor = Teuchos::null;
    Teuchos::RCP<const MapExtractor> thyColMapExtractor    = Teuchos::null;
    if(bThyraMode == true) {
      // build Thyra-style map extractors
      std::vector<Teuchos::RCP<const Map> > thyRgMapExtractorMaps(numRows, Teuchos::null);
      for (size_t r = 0; r < numRows; r++) {
        RCP<const Map> rMap = rangeMapExtractor->getMap(r);
        RCP<const Map> shrinkedMap = MapUtils::shrinkMapGIDs(*rMap,*rMap);
        RCP<const StridedMap> strRangeMap = Teuchos::rcp_dynamic_cast<const StridedMap>(rMap);
        if(strRangeMap != Teuchos::null) {
          std::vector<size_t> strInfo = strRangeMap->getStridingData();
          GlobalOrdinal offset = strRangeMap->getOffset();
          LocalOrdinal stridedBlockId = strRangeMap->getStridedBlockId();
          RCP<const StridedMap> strShrinkedMap = Teuchos::rcp(new StridedMap(shrinkedMap, strInfo, shrinkedMap->getIndexBase(), stridedBlockId, offset));
          thyRgMapExtractorMaps[r] = strShrinkedMap;
        } else {
          thyRgMapExtractorMaps[r] = shrinkedMap;
        }
        TEUCHOS_TEST_FOR_EXCEPTION(thyRgMapExtractorMaps[r]->getNodeNumElements()  != rMap->getNodeNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: Thyra-style range map extractor contains faulty data.")
      }
      RCP<const Map> fullThyRangeMap = MapUtils::concatenateMaps(thyRgMapExtractorMaps);
      thyRangeMapExtractor = MapExtractorFactory::Build(fullThyRangeMap,thyRgMapExtractorMaps,true);
      std::vector<Teuchos::RCP<const Map> > thyDoMapExtractorMaps (numCols, Teuchos::null);
      std::vector<Teuchos::RCP<const Map> > thyColMapExtractorMaps(numCols, Teuchos::null);
      for (size_t c = 0; c < numCols; c++) {
        RCP<const Map> cMap   = domainMapExtractor->getMap(c);

        RCP<const Map> shrinkedDomainMap = MapUtils::shrinkMapGIDs(*cMap,*cMap);
        RCP<const StridedMap> strDomainMap = Teuchos::rcp_dynamic_cast<const StridedMap>(cMap);
        if(strDomainMap != Teuchos::null) {
          std::vector<size_t> strInfo = strDomainMap->getStridingData();
          GlobalOrdinal offset = strDomainMap->getOffset();
          LocalOrdinal stridedBlockId = strDomainMap->getStridedBlockId();
          RCP<const StridedMap> strShrinkedDomainMap = Teuchos::rcp(new StridedMap(shrinkedDomainMap, strInfo, shrinkedDomainMap->getIndexBase(), stridedBlockId, offset));
          thyDoMapExtractorMaps[c]  = strShrinkedDomainMap;
        } else {
          thyDoMapExtractorMaps[c]  = shrinkedDomainMap;
        }
        RCP<const Map> colMap = myColumnMapExtractor->getMap(c);
        RCP<const Map> shrinkedColMap = MapUtils::shrinkMapGIDs(*colMap,*cMap);
        RCP<const StridedMap> strColMap = Teuchos::rcp_dynamic_cast<const StridedMap>(colMap);
        if(strColMap != Teuchos::null) {
          std::vector<size_t> strInfo = strColMap->getStridingData();
          GlobalOrdinal offset = strColMap->getOffset();
          LocalOrdinal stridedBlockId = strColMap->getStridedBlockId();
          RCP<const StridedMap> strShrinkedColMap = Teuchos::rcp(new StridedMap(shrinkedColMap, strInfo, shrinkedColMap->getIndexBase(), stridedBlockId, offset));
          thyColMapExtractorMaps[c]  = strShrinkedColMap;
        } else {
          thyColMapExtractorMaps[c]  = shrinkedColMap;
        }

        TEUCHOS_TEST_FOR_EXCEPTION(thyColMapExtractorMaps[c]->getNodeNumElements() != colMap->getNodeNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: Thyra-style column map extractor contains faulty data.")
        TEUCHOS_TEST_FOR_EXCEPTION(thyDoMapExtractorMaps[c]->getNodeNumElements()  != cMap->getNodeNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: Thyra-style domain map extractor contains faulty data.")
      }
      RCP<const Map> fullThyDomainMap = MapUtils::concatenateMaps(thyDoMapExtractorMaps);
      RCP<const Map> fullThyColumnMap = MapUtils::concatenateMaps(thyColMapExtractorMaps);
      thyDomainMapExtractor = MapExtractorFactory::Build(fullThyDomainMap,thyDoMapExtractorMaps,true);
      thyColMapExtractor    = MapExtractorFactory::Build(fullThyColumnMap,thyColMapExtractorMaps,true);
    }
    // create submatrices
    std::vector<Teuchos::RCP<Matrix> > subMatrices(numRows*numCols, Teuchos::null);
    for (size_t r = 0; r < numRows; r++) {
      for (size_t c = 0; c < numCols; c++) {
        // create empty CrsMatrix objects
        // make sure that the submatrices are defined using the right row maps (either Thyra or xpetra style)
        // Note: we're reserving a little bit too much memory for the submatrices, but should be still reasonable
        if(bThyraMode == true)
          subMatrices[r*numCols+c] = MatrixFactory::Build (thyRangeMapExtractor->getMap(r,true),input.getNodeMaxNumRowEntries());
        else
          subMatrices[r*numCols+c] = MatrixFactory::Build (rangeMapExtractor->getMap(r),input.getNodeMaxNumRowEntries());
      }
    }

    // We need a vector which lives on the column map of input and stores the block id that the column belongs to.
    // create a vector on the domain map. Loop over it and fill in the corresponding block id numbers
    // create a vector on the column map and import the data
    // Importer: source map is non-overlapping. Target map is overlapping
    // call colMap.Import(domMap,Importer,Insert)
    // do the same with "Add" to make sure only one processor is responsible for the different GIDs!
#if 0 // TAW needs to be fixed (does not compile for Scalar=complex)
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
    for (size_t rrr = 0; rrr < input.getDomainMap()->getNodeNumElements(); rrr++) {
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
    for (size_t rr = 0; rr < input.getRowMap()->getNodeNumElements(); rr++) {

      // global row id to extract data from global monolithic matrix
      GlobalOrdinal growid = input.getRowMap()->getGlobalElement(rr); // LID -> GID (column)

      // Find the block id in range map extractor that belongs to same global row id.
      // We assume the global ids to be unique in the map extractor
      // The MapExtractor objects may be constructed using the thyra mode. However, we use
      // global GID ids (as we have a single input matrix). The only problematic thing could
      // be that the GIDs of the input matrix are inconsistent with the GIDs in the map extractor.
      // Note, that for the Thyra mode the GIDs in the map extractor are generated using the ordering
      // of the blocks.
      size_t rowBlockId = rangeMapExtractor->getMapIndexForGID(growid);

      // global row id used for subblocks to insert information
      GlobalOrdinal subblock_growid = growid; // for Xpetra-style numbering the global row ids are not changing
      if(bThyraMode == true) {
        // find local row id associated with growid in the corresponding subblock
        LocalOrdinal lrowid = rangeMapExtractor->getMap(rowBlockId)->getLocalElement(growid);
        // translate back local row id to global row id for the subblock
        subblock_growid = thyRangeMapExtractor->getMap(rowBlockId,true)->getGlobalElement(lrowid);
      }

      // extract matrix entries from input matrix
      // we use global ids since we have to determine the corresponding
      // block column ids using the global ids anyway
      Teuchos::ArrayView<const LocalOrdinal> indices;
      Teuchos::ArrayView<const Scalar> vals;
      input.getLocalRowView(rr, indices, vals);

      std::vector<Teuchos::Array<GlobalOrdinal> > blockColIdx (numCols, Teuchos::Array<GlobalOrdinal>());
      std::vector<Teuchos::Array<Scalar> >        blockColVals(numCols, Teuchos::Array<Scalar>());

      for(size_t i=0; i<(size_t)indices.size(); i++) {
        // gobal column id to extract data from full monolithic matrix
        GlobalOrdinal gcolid = input.getColMap()->getGlobalElement(indices[i]);

        size_t colBlockId = myColumnMapExtractor->getMapIndexForGID(gcolid); // old buggy thing
        //size_t colBlockId = Teuchos::as<size_t>(coCheckData[indices[i]]);

        // global column id used for subblocks to insert information
        GlobalOrdinal subblock_gcolid = gcolid; // for Xpetra-style numbering the global col ids are not changing
        if(bThyraMode == true) {
          // find local col id associated with gcolid in the corresponding subblock
          LocalOrdinal lcolid = myColumnMapExtractor->getMap(colBlockId)->getLocalElement(gcolid);
          // translate back local col id to global col id for the subblock
          subblock_gcolid = thyColMapExtractor->getMap(colBlockId,true)->getGlobalElement(lcolid);
        }
        blockColIdx [colBlockId].push_back(subblock_gcolid);
        blockColVals[colBlockId].push_back(vals[i]);
      }

      for (size_t c = 0; c < numCols; c++) {
        subMatrices[rowBlockId*numCols+c]->insertGlobalValues(subblock_growid,blockColIdx[c].view(0,blockColIdx[c].size()),blockColVals[c].view(0,blockColVals[c].size()));
      }

    }

    // call fill complete on subblocks and create BlockedCrsOperator
    RCP<BlockedCrsMatrix> bA = Teuchos::null;
    if(bThyraMode == true) {
      for (size_t r = 0; r < numRows; r++) {
        for (size_t c = 0; c < numCols; c++) {
          subMatrices[r*numCols+c]->fillComplete(thyDomainMapExtractor->getMap(c,true), thyRangeMapExtractor->getMap(r,true));
        }
      }
      bA = Teuchos::rcp(new BlockedCrsMatrix(thyRangeMapExtractor, thyDomainMapExtractor, 10 /*input.getRowMap()->getNodeMaxNumRowEntries()*/));
    } else {
      for (size_t r = 0; r < numRows; r++) {
        for (size_t c = 0; c < numCols; c++) {
          subMatrices[r*numCols+c]->fillComplete(domainMapExtractor->getMap(c), rangeMapExtractor->getMap(r));
        }
      }
      bA = Teuchos::rcp(new BlockedCrsMatrix(rangeMapExtractor, domainMapExtractor, 10 /*input.getRowMap()->getNodeMaxNumRowEntries()*/));
    }

    for (size_t r = 0; r < numRows; r++) {
      for (size_t c = 0; c < numCols; c++) {
        bA->setMatrix(r,c,subMatrices[r*numCols+c]);
      }
    }
    return bA;
  }

  /** Given a matrix A, detect too small diagonals and replace any found with ones. */

  static void CheckRepairMainDiagonal(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& Ac,
                                 bool const &repairZeroDiagonals, Teuchos::FancyOStream &fos,
                                 const typename Teuchos::ScalarTraits<Scalar>::magnitudeType threshold = Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero())
  {
    typedef typename Teuchos::ScalarTraits<Scalar> TST;
    Scalar one = TST::one();

    Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList());
    p->set("DoOptimizeStorage", true);

    RCP<const Map> rowMap = Ac->getRowMap();
    RCP<Vector> diagVec = VectorFactory::Build(rowMap);
    Ac->getLocalDiagCopy(*diagVec);

    LocalOrdinal lZeroDiags = 0;
    Teuchos::ArrayRCP< Scalar > diagVal = diagVec->getDataNonConst(0);

    for (size_t i = 0; i < rowMap->getNodeNumElements(); i++) {
      if (TST::magnitude(diagVal[i]) <= threshold) {
        lZeroDiags++;
      }
    }
    GlobalOrdinal gZeroDiags;
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
      for (size_t r = 0; r < rowMap->getNodeNumElements(); r++) {
        if (TST::magnitude(diagVal[r]) <= threshold) {
          GlobalOrdinal grid = rowMap->getGlobalElement(r);
          indout[0] = grid;
          valout[0] = one;
          fixDiagMatrix->insertGlobalValues(grid,indout(), valout());
        }
      }
      {
        Teuchos::TimeMonitor m1(*Teuchos::TimeMonitor::getNewTimer("CheckRepairMainDiagonal: fillComplete1"));
        fixDiagMatrix->fillComplete(Ac->getDomainMap(),Ac->getRangeMap());
      }

      RCP<Matrix> newAc;
      Xpetra::MatrixMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::TwoMatrixAdd(*Ac, false, 1.0, *fixDiagMatrix, false, 1.0, newAc, fos);
      if (Ac->IsView("stridedMaps"))
        newAc->CreateView("stridedMaps", Ac);

      Ac = Teuchos::null;  // free singular matrix
      fixDiagMatrix = Teuchos::null;
      Ac = newAc;          // set fixed non-singular matrix

      // call fillComplete with optimized storage option set to true
      // This is necessary for new faster Epetra MM kernels.
      {
        Teuchos::TimeMonitor m1(*Teuchos::TimeMonitor::getNewTimer("CheckRepairMainDiagonal: fillComplete2"));
        Ac->fillComplete(p);
      }
    } // end repair



    // print some output
    fos << "CheckRepairMainDiagonal: " << (repairZeroDiagonals ? "repaired " : "found ")
              << gZeroDiags << " too small entries on main diagonal of Ac." << std::endl;

#ifdef HAVE_XPETRA_DEBUG // only for debugging
    // check whether Ac has been repaired...
    Ac->getLocalDiagCopy(*diagVec);
    Teuchos::ArrayRCP< Scalar > diagVal2 = diagVec->getDataNonConst(0);
    for (size_t r = 0; r < Ac->getRowMap()->getNodeNumElements(); r++) {
      if (TST::magnitude(diagVal[r]) <= threshold) {
        fos << "Error: there are too small entries left on diagonal after repair..." << std::endl;
        break;
      }
    }
#endif
  } //CheckRepairMainDiagonal

};

} // end namespace Xpetra

#define XPETRA_MATRIXUTILS_SHORT

#endif // PACKAGES_XPETRA_SUP_MATRIX_UTILS_HPP_
