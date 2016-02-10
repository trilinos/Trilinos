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

#include "Xpetra_Matrix.hpp"
#include "Xpetra_MatrixMatrix.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"

#include "Xpetra_Map.hpp"
#include "Xpetra_StridedMap.hpp"
#include "Xpetra_StridedMapFactory.hpp"
#include "Xpetra_MapExtractor.hpp"
#include "Xpetra_MatrixFactory.hpp"
#include "Xpetra_BlockedCrsMatrix.hpp"

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

  static Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > concatenateMaps(std::vector<Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > > & subMaps) {

    // merge submaps to global map
    std::vector<GlobalOrdinal> gids;
    for(size_t tt = 0; tt<subMaps.size(); ++tt) {
      Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > subMap = subMaps[tt];
      for(LocalOrdinal l = 0; l < Teuchos::as<LocalOrdinal>(subMap->getNodeNumElements()); ++l) {
        GlobalOrdinal gid = subMap->getGlobalElement(l);
        gids.push_back(gid);
      }
    }

    const GO INVALID = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
    //std::sort(gids.begin(), gids.end());
    //gids.erase(std::unique(gids.begin(), gids.end()), gids.end());
    Teuchos::ArrayView<GO> gidsView(&gids[0], gids.size());
    Teuchos::RCP<Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > fullMap = Xpetra::MapFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(subMaps[0]->lib(), INVALID, gidsView, subMaps[0]->getIndexBase(), subMaps[0]->getComm());
    return fullMap;
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
  */
  static Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > SplitMatrix(
                       const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& input,
                       Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rangeMapExtractor,
                       Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > domainMapExtractor,
                       Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > columnMapExtractor = Teuchos::null) {
    size_t numRows  = rangeMapExtractor->NumMaps();
    size_t numCols  = domainMapExtractor->NumMaps();

    RCP<const Map> fullRangeMap  = rangeMapExtractor->getFullMap();
    RCP<const Map> fullDomainMap = domainMapExtractor->getFullMap();

    TEUCHOS_TEST_FOR_EXCEPTION(fullRangeMap->getMaxAllGlobalIndex() != input.getRowMap()->getMaxAllGlobalIndex(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: RangeMapExtractor incompatible to row map of input matrix.")
    TEUCHOS_TEST_FOR_EXCEPTION(fullRangeMap->getGlobalNumElements() != input.getRowMap()->getGlobalNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: RangeMapExtractor incompatible to row map of input matrix.")
    TEUCHOS_TEST_FOR_EXCEPTION(fullRangeMap->getNodeNumElements()   != input.getRowMap()->getNodeNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: RangeMapExtractor incompatible to row map of input matrix.")
    TEUCHOS_TEST_FOR_EXCEPTION(fullRangeMap->getMaxAllGlobalIndex() != input.getRangeMap()->getMaxAllGlobalIndex(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: RangeMapExtractor incompatible to row map of input matrix.")
    TEUCHOS_TEST_FOR_EXCEPTION(fullRangeMap->getGlobalNumElements() != input.getRangeMap()->getGlobalNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: RangeMapExtractor incompatible to row map of input matrix.")
    TEUCHOS_TEST_FOR_EXCEPTION(fullRangeMap->getNodeNumElements()   != input.getRangeMap()->getNodeNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: RangeMapExtractor incompatible to row map of input matrix.")

    TEUCHOS_TEST_FOR_EXCEPTION(fullDomainMap->getMaxAllGlobalIndex() != input.getColMap()->getMaxAllGlobalIndex(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: DomainMapExtractor incompatible to domain map of input matrix.")
    TEUCHOS_TEST_FOR_EXCEPTION(fullDomainMap->getMaxAllGlobalIndex() != input.getDomainMap()->getMaxAllGlobalIndex(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: DomainMapExtractor incompatible to domain map of input matrix.")
    TEUCHOS_TEST_FOR_EXCEPTION(fullDomainMap->getGlobalNumElements() != input.getDomainMap()->getGlobalNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: DomainMapExtractor incompatible to domain map of input matrix.")
    TEUCHOS_TEST_FOR_EXCEPTION(fullDomainMap->getNodeNumElements()   != input.getDomainMap()->getNodeNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: DomainMapExtractor incompatible to domain map of input matrix.")

    // check column map extractor
    Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > myColumnMapExtractor = Teuchos::null;
    if(columnMapExtractor == Teuchos::null) {
      TEUCHOS_TEST_FOR_EXCEPTION(domainMapExtractor->getThyraMode() == true, Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: Auto generation of column map extractor not supported for Thyra style numbering.");
      std::vector<Teuchos::RCP<const Map> > ovlxmaps(numCols, Teuchos::null);
      for(size_t c = 0; c < numCols; c++) {
        Teuchos::RCP<const Map> colMap =
          Xpetra::MatrixUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::findColumnSubMap(input, *(domainMapExtractor->getMap(c))); // TODO what about Thyra style numbering? Probably no changes necessary as we obtain a pseudo map with unique GIDs.
        ovlxmaps[c] = colMap;
      }

      Teuchos::RCP<const Map> fullColMap =
          Xpetra::MatrixUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::concatenateMaps(ovlxmaps);

      myColumnMapExtractor = Xpetra::MapExtractorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(fullColMap,ovlxmaps);
    } else
      myColumnMapExtractor = columnMapExtractor; // use user-provided column map extractor.

    // create submatrices
    std::vector<Teuchos::RCP<CrsMatrix> > subMatrices(numRows*numCols, Teuchos::null);

    for (size_t r = 0; r < numRows; r++) {
      for (size_t c = 0; c < numCols; c++) {
        // create empty CrsMatrix objects (we're reserving a little bit too much memory, but should be still reasonable)
        subMatrices[r*numCols+c] = Xpetra::CrsMatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build (rangeMapExtractor->getMap(r),input.getNodeMaxNumRowEntries());
      }
    }

    // loop over all rows of input matrix
    for (size_t rr = 0; rr < input.getRowMap()->getNodeNumElements(); rr++) {

      GlobalOrdinal growid = input.getRowMap()->getGlobalElement(rr); // LID -> GID (column)

      // Find the block id in range map extractor that belongs to same global row id.
      // We assume the global ids to be unique in the map extractor
      // The MapExtractor objects may be constructed using the thyra mode. However, we use
      // global GID ids (as we have a single input matrix). The only problematic thing could
      // be that the GIDs of the input matrix are inconsistent with the GIDs in the map extractor.
      // Note, that for the Thyra mode the GIDs in the map extractor are generated using the ordering
      // of the blocks.
      size_t rowBlockId = rangeMapExtractor->getMapIndexForGID(growid);

      // extract matrix entries from input matrix
      // we use global ids since we have to determine the corresponding
      // block column ids using the global ids anyway
      Teuchos::ArrayView<const LocalOrdinal> indices;
      Teuchos::ArrayView<const Scalar> vals;
      input.getLocalRowView(rr, indices, vals);

      std::vector<Teuchos::Array<GlobalOrdinal> > blockColIdx (numCols, Teuchos::Array<GlobalOrdinal>());
      std::vector<Teuchos::Array<Scalar> >        blockColVals(numCols, Teuchos::Array<Scalar>());

      for(size_t i=0; i<(size_t)indices.size(); i++) {
        GlobalOrdinal gcolid = input.getColMap()->getGlobalElement(indices[i]);
        size_t colBlockId = myColumnMapExtractor->getMapIndexForGID(gcolid);
        blockColIdx [colBlockId].push_back(gcolid);
        blockColVals[colBlockId].push_back( vals[i]);
      }

      for (size_t c = 0; c < numCols; c++) {
        subMatrices[rowBlockId*numCols+c]->insertGlobalValues(growid,blockColIdx[c].view(0,blockColIdx[c].size()),blockColVals[c].view(0,blockColVals[c].size()));
      }

    }

    for (size_t r = 0; r < numRows; r++) {
      for (size_t c = 0; c < numCols; c++) {
        subMatrices[r*numCols+c]->fillComplete(domainMapExtractor->getMap(c), rangeMapExtractor->getMap(r));
      }
    }

    Teuchos::RCP<BlockedCrsMatrix> bA = Teuchos::rcp(new BlockedCrsMatrix(
        rangeMapExtractor, domainMapExtractor, 10 /*input.getRowMap()->getNodeMaxNumRowEntries()*/));


    for (size_t r = 0; r < numRows; r++) {
      for (size_t c = 0; c < numCols; c++) {
        bA->setMatrix(r,c,subMatrices[r*numCols+c]);
      }
    }

    return bA;
  }

};

} // end namespace Xpetra

#define XPETRA_MATRIXUTILS_SHORT

#endif // PACKAGES_XPETRA_SUP_MATRIX_UTILS_HPP_
