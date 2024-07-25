// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_BlockingTpetra.hpp"
#include "Teko_Utilities.hpp"

#include "Tpetra_Vector.hpp"
#include "Tpetra_Map.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

namespace Teko {
namespace TpetraHelpers {
namespace Blocking {

/** Build maps to make other conversions. This function builds a map
 * using a vector of global ids local to this processor.  It also builds
 * a seperate map that (globally) starts from zero. For instance if the following
 * GIDs are passed in PID = 0, GID = [0 2 4 6 8] and PID = 1, GID = [10 12 14]
 * the the two maps created are
 *    Global Map = [(PID=0,GID=[0 2 4 6 8]), (PID=1,GID=[10 12 14])]
 *    Contiguous Map = [(PID=0,GID=[0 1 2 3 4]), (PID=1,GID=[5 6 7])]
 *
 * \param[in] gid Local global IDs to use
 * \param[in] comm Communicator to use in construction of the maps
 *
 * \returns A pair of maps: (Global Map, Contiguous Map)
 */
const MapPair buildSubMap(const std::vector<GO>& gid, const Teuchos::Comm<int>& comm) {
  using GST                                     = Tpetra::global_size_t;
  const GST invalid                             = Teuchos::OrdinalTraits<GST>::invalid();
  Teuchos::RCP<Tpetra::Map<LO, GO, NT> > gidMap = rcp(
      new Tpetra::Map<LO, GO, NT>(invalid, Teuchos::ArrayView<const GO>(gid), 0, rcpFromRef(comm)));
  Teuchos::RCP<Tpetra::Map<LO, GO, NT> > contigMap =
      rcp(new Tpetra::Map<LO, GO, NT>(invalid, gid.size(), 0, rcpFromRef(comm)));

  return std::make_pair(gidMap, contigMap);
}

/** Build the Export/Import objects that take the single vector global map and
 * build individual sub maps.
 *
 * \param[in] baseMap Single global vector map.
 * \param[in] maps Pair of maps containing the global sub map and the contiguous sub map.
 *                 These come directly from <code>buildSubMap</code>.
 *
 * \returns A pair containing pointers to the Import/Export objects.
 */
const ImExPair buildExportImport(const Tpetra::Map<LO, GO, NT>& baseMap, const MapPair& maps) {
  return std::make_pair(rcp(new Tpetra::Import<LO, GO, NT>(rcpFromRef(baseMap), maps.first)),
                        rcp(new Tpetra::Export<LO, GO, NT>(maps.first, rcpFromRef(baseMap))));
}

/** Using a list of map pairs created by <code>buildSubMap</code>, buidl the corresponding
 * multi-vector objects.
 *
 * \param[in] maps Map pairs created by <code>buildSubMap</code>
 * \param[in,out] vectors Vector objects created using the Contiguous maps
 * \param[in] count Number of multivectors to build.
 */
void buildSubVectors(const std::vector<MapPair>& maps,
                     std::vector<RCP<Tpetra::MultiVector<ST, LO, GO, NT> > >& vectors, int count) {
  std::vector<MapPair>::const_iterator mapItr;

  // loop over all maps
  for (mapItr = maps.begin(); mapItr != maps.end(); mapItr++) {
    // add new elements to vectors
    RCP<Tpetra::MultiVector<ST, LO, GO, NT> > mv =
        rcp(new Tpetra::MultiVector<ST, LO, GO, NT>((*mapItr).second, count));
    vectors.push_back(mv);
  }
}

/** Copy the contents of a global vector into many sub-vectors created by
 * <code>buildSubVectors</code>.
 *
 * \param[in,out] many Sub-vector to be filled by this operation created by
 * <code>buildSubVectors</code>. \param[in] one The source vector. \param[in] subImport A list of
 * import objects to use in copying.
 */
void one2many(std::vector<RCP<Tpetra::MultiVector<ST, LO, GO, NT> > >& many,
              const Tpetra::MultiVector<ST, LO, GO, NT>& single,
              const std::vector<RCP<Tpetra::Import<LO, GO, NT> > >& subImport) {
  // std::vector<RCP<Epetra_Vector> >::const_iterator vecItr;
  std::vector<RCP<Tpetra::MultiVector<ST, LO, GO, NT> > >::const_iterator vecItr;
  std::vector<RCP<Tpetra::Import<LO, GO, NT> > >::const_iterator impItr;

  // using Importers fill the sub vectors from the mama vector
  for (vecItr = many.begin(), impItr = subImport.begin(); vecItr != many.end();
       ++vecItr, ++impItr) {
    // for ease of access to the destination
    RCP<Tpetra::MultiVector<ST, LO, GO, NT> > destVec = *vecItr;

    // extract the map with global indicies from the current vector
    const Tpetra::Map<LO, GO, NT>& globalMap = *(*impItr)->getTargetMap();

    // build the import vector as a view on the destination
    GO lda      = destVec->getStride();
    GO destSize = destVec->getGlobalLength() * destVec->getNumVectors();
    std::vector<ST> destArray(destSize);
    Teuchos::ArrayView<ST> destVals(destArray);
    destVec->get1dCopy(destVals, lda);
    Tpetra::MultiVector<ST, LO, GO, NT> importVector(rcpFromRef(globalMap), destVals, lda,
                                                     destVec->getNumVectors());

    // perform the import
    importVector.doImport(single, **impItr, Tpetra::INSERT);

    Tpetra::Import<LO, GO, NT> importer(destVec->getMap(), destVec->getMap());
    importVector.replaceMap(destVec->getMap());
    destVec->doExport(importVector, importer, Tpetra::INSERT);
  }
}

/** Copy the contents of many sub vectors (created from a contigous sub maps) to a single global
 * vector. This should have the map used to create the Export/Import objects in the
 * <code>buildExportImport</code> function. If more then one sub vector contains values for a
 * particular GID in the single vector then the value in the final vector will be that of the last
 * sub vector (in the list <code>many</code>).
 *
 * \param[in,out] one The single vector to be filled by this operation.
 * \param[in] many Sub-vectors created by <code>buildSubVectors</code> used to fill
 * <code>one</code>. \param[in] subExport A list of export objects to use in copying.
 */
void many2one(Tpetra::MultiVector<ST, LO, GO, NT>& one,
              const std::vector<RCP<const Tpetra::MultiVector<ST, LO, GO, NT> > >& many,
              const std::vector<RCP<Tpetra::Export<LO, GO, NT> > >& subExport) {
  // std::vector<RCP<const Epetra_Vector> >::const_iterator vecItr;
  std::vector<RCP<const Tpetra::MultiVector<ST, LO, GO, NT> > >::const_iterator vecItr;
  std::vector<RCP<Tpetra::Export<LO, GO, NT> > >::const_iterator expItr;

  // using Exporters fill the empty vector from the sub-vectors
  for (vecItr = many.begin(), expItr = subExport.begin(); vecItr != many.end();
       ++vecItr, ++expItr) {
    // for ease of access to the source
    RCP<const Tpetra::MultiVector<ST, LO, GO, NT> > srcVec = *vecItr;

    // extract the map with global indicies from the current vector
    const Tpetra::Map<LO, GO, NT>& globalMap = *(*expItr)->getSourceMap();

    // build the export vector as a view of the destination
    GO lda     = srcVec->getStride();
    GO srcSize = srcVec->getGlobalLength() * srcVec->getNumVectors();
    std::vector<ST> srcArray(srcSize);
    Teuchos::ArrayView<ST> srcVals(srcArray);
    srcVec->get1dCopy(srcVals, lda);
    Tpetra::MultiVector<ST, LO, GO, NT> exportVector(rcpFromRef(globalMap), srcVals, lda,
                                                     srcVec->getNumVectors());

    // perform the export
    one.doExport(exportVector, **expItr, Tpetra::INSERT);
  }
}

/** This function will return an IntVector that is constructed with a column map.
 * The vector will be filled with -1 if there is not a corresponding entry in the
 * sub-block row map. The other columns will be filled with the contiguous row map
 * values.
 */
RCP<Tpetra::Vector<GO, LO, GO, NT> > getSubBlockColumnGIDs(
    const Tpetra::CrsMatrix<ST, LO, GO, NT>& A, const MapPair& mapPair) {
  RCP<const Tpetra::Map<LO, GO, NT> > blkGIDMap    = mapPair.first;
  RCP<const Tpetra::Map<LO, GO, NT> > blkContigMap = mapPair.second;

  // fill index vector for rows
  Tpetra::Vector<GO, LO, GO, NT> rIndex(A.getRowMap(), true);
  for (size_t i = 0; i < A.getLocalNumRows(); i++) {
    // LID is need to get to contiguous map
    LO lid = blkGIDMap->getLocalElement(
        A.getRowMap()->getGlobalElement(i));  // this LID makes me nervous
    if (lid > -1)
      rIndex.replaceLocalValue(i, blkContigMap->getGlobalElement(lid));
    else
      rIndex.replaceLocalValue(i, -1);
  }

  // get relavant column indices

  Tpetra::Import<LO, GO, NT> import(A.getRowMap(), A.getColMap());
  RCP<Tpetra::Vector<GO, LO, GO, NT> > cIndex =
      rcp(new Tpetra::Vector<GO, LO, GO, NT>(A.getColMap(), true));
  cIndex->doImport(rIndex, import, Tpetra::INSERT);

  return cIndex;
}

// build a single subblock Epetra_CrsMatrix
RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > buildSubBlock(
    int i, int j, const RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> >& A,
    const std::vector<MapPair>& subMaps) {
  // get the number of variables families
  int numVarFamily = subMaps.size();

  TEUCHOS_ASSERT(i >= 0 && i < numVarFamily);
  TEUCHOS_ASSERT(j >= 0 && j < numVarFamily);
  const RCP<const Tpetra::Map<LO, GO, NT> > gRowMap = subMaps[i].first;   // new GIDs
  const RCP<const Tpetra::Map<LO, GO, NT> > rowMap  = subMaps[i].second;  // contiguous GIDs
  const RCP<const Tpetra::Map<LO, GO, NT> > colMap  = subMaps[j].second;

  const RCP<Tpetra::Vector<GO, LO, GO, NT> > plocal2ContigGIDs =
      getSubBlockColumnGIDs(*A, subMaps[j]);
  Tpetra::Vector<GO, LO, GO, NT>& local2ContigGIDs = *plocal2ContigGIDs;

  // get entry information
  LO numMyRows     = rowMap->getLocalNumElements();
  LO maxNumEntries = A->getLocalMaxNumRowEntries();

  // for extraction
  auto indices = typename Tpetra::CrsMatrix<ST, LO, GO, NT>::nonconst_local_inds_host_view_type(
      Kokkos::ViewAllocateWithoutInitializing("rowIndices"), maxNumEntries);
  auto values = typename Tpetra::CrsMatrix<ST, LO, GO, NT>::nonconst_values_host_view_type(
      Kokkos::ViewAllocateWithoutInitializing("rowIndices"), maxNumEntries);

  // for insertion
  std::vector<GO> colIndices(maxNumEntries);
  std::vector<ST> colValues(maxNumEntries);

  // for counting
  std::vector<size_t> nEntriesPerRow(numMyRows);

  const size_t invalid = Teuchos::OrdinalTraits<size_t>::invalid();

  // insert each row into subblock
  // Count the number of entries per row in the new matrix
  for (LO localRow = 0; localRow < numMyRows; localRow++) {
    size_t numEntries = invalid;
    GO globalRow      = gRowMap->getGlobalElement(localRow);
    LO lid            = A->getRowMap()->getLocalElement(globalRow);
    TEUCHOS_ASSERT(lid > -1);

    A->getLocalRowCopy(lid, indices, values, numEntries);

    LO numOwnedCols = 0;
    for (size_t localCol = 0; localCol < numEntries; localCol++) {
      TEUCHOS_ASSERT(indices(localCol) > -1);

      // if global id is not owned by this column

      GO gid = local2ContigGIDs.getData()[indices[localCol]];
      if (gid == -1) continue;  // in contiguous row
      numOwnedCols++;
    }
    nEntriesPerRow[localRow] = numOwnedCols;
  }

  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > mat = rcp(
      new Tpetra::CrsMatrix<ST, LO, GO, NT>(rowMap, Teuchos::ArrayView<size_t>(nEntriesPerRow)));

  // insert each row into subblock
  // let FillComplete handle column distribution
  for (LO localRow = 0; localRow < numMyRows; localRow++) {
    size_t numEntries = invalid;
    GO globalRow      = gRowMap->getGlobalElement(localRow);
    LO lid            = A->getRowMap()->getLocalElement(globalRow);
    GO contigRow      = rowMap->getGlobalElement(localRow);
    TEUCHOS_ASSERT(lid > -1);

    A->getLocalRowCopy(lid, indices, values, numEntries);

    LO numOwnedCols = 0;
    for (size_t localCol = 0; localCol < numEntries; localCol++) {
      TEUCHOS_ASSERT(indices(localCol) > -1);

      // if global id is not owned by this column

      GO gid = local2ContigGIDs.getData()[indices(localCol)];
      if (gid == -1) continue;  // in contiguous row

      colIndices[numOwnedCols] = gid;
      colValues[numOwnedCols]  = values(localCol);
      numOwnedCols++;
    }

    // insert it into the new matrix
    colIndices.resize(numOwnedCols);
    colValues.resize(numOwnedCols);
    mat->insertGlobalValues(contigRow, Teuchos::ArrayView<GO>(colIndices),
                            Teuchos::ArrayView<ST>(colValues));
    colIndices.resize(maxNumEntries);
    colValues.resize(maxNumEntries);
  }

  // fill it and automagically optimize the storage
  mat->fillComplete(colMap, rowMap);

  return mat;
}

// build a single subblock Epetra_CrsMatrix
void rebuildSubBlock(int i, int j, const RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> >& A,
                     const std::vector<MapPair>& subMaps, Tpetra::CrsMatrix<ST, LO, GO, NT>& mat) {
  // get the number of variables families
  int numVarFamily = subMaps.size();

  TEUCHOS_ASSERT(i >= 0 && i < numVarFamily);
  TEUCHOS_ASSERT(j >= 0 && j < numVarFamily);

  const Tpetra::Map<LO, GO, NT>& gRowMap = *subMaps[i].first;   // new GIDs
  const Tpetra::Map<LO, GO, NT>& rowMap  = *subMaps[i].second;  // contiguous GIDs
  const Tpetra::Map<LO, GO, NT>& colMap  = *subMaps[j].second;

  const RCP<Tpetra::Vector<GO, LO, GO, NT> > plocal2ContigGIDs =
      getSubBlockColumnGIDs(*A, subMaps[j]);
  Tpetra::Vector<GO, LO, GO, NT>& local2ContigGIDs = *plocal2ContigGIDs;

  mat.resumeFill();
  mat.setAllToScalar(0.0);

  // get entry information
  LO numMyRows     = rowMap.getLocalNumElements();
  LO maxNumEntries = A->getLocalMaxNumRowEntries();

  // for extraction
  auto indices = typename Tpetra::CrsMatrix<ST, LO, GO, NT>::nonconst_local_inds_host_view_type(
      Kokkos::ViewAllocateWithoutInitializing("rowIndices"), maxNumEntries);
  auto values = typename Tpetra::CrsMatrix<ST, LO, GO, NT>::nonconst_values_host_view_type(
      Kokkos::ViewAllocateWithoutInitializing("rowIndices"), maxNumEntries);

  // for insertion
  std::vector<GO> colIndices(maxNumEntries);
  std::vector<ST> colValues(maxNumEntries);

  const size_t invalid = Teuchos::OrdinalTraits<size_t>::invalid();

  // insert each row into subblock
  // let FillComplete handle column distribution
  for (LO localRow = 0; localRow < numMyRows; localRow++) {
    size_t numEntries = invalid;
    GO globalRow      = gRowMap.getGlobalElement(localRow);
    LO lid            = A->getRowMap()->getLocalElement(globalRow);
    GO contigRow      = rowMap.getGlobalElement(localRow);
    TEUCHOS_ASSERT(lid > -1);

    A->getLocalRowCopy(lid, indices, values, numEntries);

    LO numOwnedCols = 0;
    for (size_t localCol = 0; localCol < numEntries; localCol++) {
      TEUCHOS_ASSERT(indices(localCol) > -1);

      // if global id is not owned by this column
      GO gid = local2ContigGIDs.getData()[indices(localCol)];
      if (gid == -1) continue;  // in contiguous row

      colIndices[numOwnedCols] = gid;
      colValues[numOwnedCols]  = values(localCol);
      numOwnedCols++;
    }

    // insert it into the new matrix
    colIndices.resize(numOwnedCols);
    colValues.resize(numOwnedCols);
    mat.sumIntoGlobalValues(contigRow, Teuchos::ArrayView<GO>(colIndices),
                            Teuchos::ArrayView<ST>(colValues));
    colIndices.resize(maxNumEntries);
    colValues.resize(maxNumEntries);
  }
  mat.fillComplete(rcpFromRef(colMap), rcpFromRef(rowMap));
}

}  // namespace Blocking
}  // namespace TpetraHelpers
}  // namespace Teko
