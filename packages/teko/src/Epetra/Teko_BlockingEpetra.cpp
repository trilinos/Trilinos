// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_BlockingEpetra.hpp"

#include "Epetra_IntVector.h"
#include "Epetra_LocalMap.h"

using Teuchos::RCP;
using Teuchos::rcp;

namespace Teko {
namespace Epetra {
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
const MapPair buildSubMap(const std::vector<int> &gid, const Epetra_Comm &comm) {
  Teuchos::RCP<Epetra_Map> gidMap    = rcp(new Epetra_Map(-1, gid.size(), &gid[0], 0, comm));
  Teuchos::RCP<Epetra_Map> contigMap = rcp(new Epetra_Map(-1, gid.size(), 0, comm));

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
const ImExPair buildExportImport(const Epetra_Map &baseMap, const MapPair &maps) {
  return std::make_pair(rcp(new Epetra_Import(*maps.first, baseMap)),
                        rcp(new Epetra_Export(*maps.first, baseMap)));
}

/** Using a list of map pairs created by <code>buildSubMap</code>, buidl the corresponding
 * multi-vector objects.
 *
 * \param[in] maps Map pairs created by <code>buildSubMap</code>
 * \param[in,out] vectors Vector objects created using the Contiguous maps
 * \param[in] count Number of multivectors to build.
 */
void buildSubVectors(const std::vector<MapPair> &maps,
                     std::vector<RCP<Epetra_MultiVector> > &vectors, int count) {
  std::vector<MapPair>::const_iterator mapItr;

  // loop over all maps
  for (mapItr = maps.begin(); mapItr != maps.end(); mapItr++) {
    // add new elements to vectors
    RCP<Epetra_MultiVector> mv = rcp(new Epetra_MultiVector(*(*mapItr).second, count));
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
void one2many(std::vector<RCP<Epetra_MultiVector> > &many, const Epetra_MultiVector &single,
              const std::vector<RCP<Epetra_Import> > &subImport) {
  // std::vector<RCP<Epetra_Vector> >::const_iterator vecItr;
  std::vector<RCP<Epetra_MultiVector> >::const_iterator vecItr;
  std::vector<RCP<Epetra_Import> >::const_iterator impItr;

  // using Importers fill the sub vectors from the mama vector
  for (vecItr = many.begin(), impItr = subImport.begin(); vecItr != many.end();
       ++vecItr, ++impItr) {
    // for ease of access to the destination
    RCP<Epetra_MultiVector> destVec = *vecItr;

    // extract the map with global indicies from the current vector
    const Epetra_BlockMap &globalMap = (*impItr)->TargetMap();

    // build the import vector as a view on the destination
    Epetra_MultiVector importVector(View, globalMap, destVec->Values(), destVec->Stride(),
                                    destVec->NumVectors());

    // perform the import
    importVector.Import(single, **impItr, Insert);
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
void many2one(Epetra_MultiVector &one, const std::vector<RCP<const Epetra_MultiVector> > &many,
              const std::vector<RCP<Epetra_Export> > &subExport) {
  // std::vector<RCP<const Epetra_Vector> >::const_iterator vecItr;
  std::vector<RCP<const Epetra_MultiVector> >::const_iterator vecItr;
  std::vector<RCP<Epetra_Export> >::const_iterator expItr;

  // using Exporters fill the empty vector from the sub-vectors
  for (vecItr = many.begin(), expItr = subExport.begin(); vecItr != many.end();
       ++vecItr, ++expItr) {
    // for ease of access to the source
    RCP<const Epetra_MultiVector> srcVec = *vecItr;

    // extract the map with global indicies from the current vector
    const Epetra_BlockMap &globalMap = (*expItr)->SourceMap();

    // build the export vector as a view of the destination
    Epetra_MultiVector exportVector(View, globalMap, srcVec->Values(), srcVec->Stride(),
                                    srcVec->NumVectors());
    one.Export(exportVector, **expItr, Insert);
  }
}

/** This function will return an IntVector that is constructed with a column map.
 * The vector will be filled with -1 if there is not a corresponding entry in the
 * sub-block row map. The other columns will be filled with the contiguous row map
 * values.
 */
RCP<Epetra_IntVector> getSubBlockColumnGIDs(const Epetra_CrsMatrix &A, const MapPair &mapPair) {
  RCP<const Epetra_Map> blkGIDMap    = mapPair.first;
  RCP<const Epetra_Map> blkContigMap = mapPair.second;

  // fill index vector for rows
  Epetra_IntVector rIndex(A.RowMap(), true);
  for (int i = 0; i < A.NumMyRows(); i++) {
    // LID is need to get to contiguous map
    int lid = blkGIDMap->LID(A.GRID(i));  // this LID makes me nervous
    if (lid > -1)
      rIndex[i] = blkContigMap->GID(lid);
    else
      rIndex[i] = -1;
  }

  // get relavant column indices
  Epetra_Import import(A.ColMap(), A.RowMap());
  RCP<Epetra_IntVector> cIndex = rcp(new Epetra_IntVector(A.ColMap(), true));
  cIndex->Import(rIndex, import, Insert);

  return cIndex;
}

// build a single subblock Epetra_CrsMatrix
RCP<Epetra_CrsMatrix> buildSubBlock(int i, int j, const Epetra_CrsMatrix &A,
                                    const std::vector<MapPair> &subMaps) {
  // get the number of variables families
  int numVarFamily = subMaps.size();

  TEUCHOS_ASSERT(i >= 0 && i < numVarFamily);
  TEUCHOS_ASSERT(j >= 0 && j < numVarFamily);

  const Epetra_Map &gRowMap = *subMaps[i].first;   // new GIDs
  const Epetra_Map &rowMap  = *subMaps[i].second;  // contiguous GIDs
  const Epetra_Map &colMap  = *subMaps[j].second;

  const RCP<Epetra_IntVector> plocal2ContigGIDs = getSubBlockColumnGIDs(A, subMaps[j]);
  Epetra_IntVector &local2ContigGIDs            = *plocal2ContigGIDs;

  RCP<Epetra_CrsMatrix> mat = rcp(new Epetra_CrsMatrix(Copy, rowMap, 0));

  // get entry information
  int numMyRows     = rowMap.NumMyElements();
  int maxNumEntries = A.MaxNumEntries();

  // for extraction
  std::vector<int> indices(maxNumEntries);
  std::vector<double> values(maxNumEntries);

  // for insertion
  std::vector<int> colIndices(maxNumEntries);
  std::vector<double> colValues(maxNumEntries);

  // insert each row into subblock
  // let FillComplete handle column distribution
  for (int localRow = 0; localRow < numMyRows; localRow++) {
    int numEntries = -1;
    int globalRow  = gRowMap.GID(localRow);
    int lid        = A.RowMap().LID(globalRow);
    int contigRow  = rowMap.GID(localRow);
    TEUCHOS_ASSERT(lid > -1);

    int err = A.ExtractMyRowCopy(lid, maxNumEntries, numEntries, &values[0], &indices[0]);
    TEUCHOS_ASSERT(err == 0);

    int numOwnedCols = 0;
    for (int localCol = 0; localCol < numEntries; localCol++) {
      TEUCHOS_ASSERT(indices[localCol] > -1);

      // if global id is not owned by this column
      int gid = local2ContigGIDs[indices[localCol]];
      if (gid == -1) continue;  // in contiguous row

      colIndices[numOwnedCols] = gid;
      colValues[numOwnedCols]  = values[localCol];
      numOwnedCols++;
    }

    // insert it into the new matrix
    mat->InsertGlobalValues(contigRow, numOwnedCols, &colValues[0], &colIndices[0]);
  }

  // fill it and automagically optimize the storage
  mat->FillComplete(colMap, rowMap);

  return mat;
}

// build a single subblock Epetra_CrsMatrix
void rebuildSubBlock(int i, int j, const Epetra_CrsMatrix &A, const std::vector<MapPair> &subMaps,
                     Epetra_CrsMatrix &mat) {
  // get the number of variables families
  int numVarFamily = subMaps.size();

  TEUCHOS_ASSERT(i >= 0 && i < numVarFamily);
  TEUCHOS_ASSERT(j >= 0 && j < numVarFamily);

  const Epetra_Map &gRowMap = *subMaps[i].first;   // new GIDs
  const Epetra_Map &rowMap  = *subMaps[i].second;  // contiguous GIDs

  const RCP<Epetra_IntVector> plocal2ContigGIDs = getSubBlockColumnGIDs(A, subMaps[j]);
  Epetra_IntVector &local2ContigGIDs            = *plocal2ContigGIDs;

  mat.PutScalar(0.0);

  // get entry information
  int numMyRows     = rowMap.NumMyElements();
  int maxNumEntries = A.MaxNumEntries();

  // for extraction
  std::vector<int> indices(maxNumEntries);
  std::vector<double> values(maxNumEntries);

  // for insertion
  std::vector<int> colIndices(maxNumEntries);
  std::vector<double> colValues(maxNumEntries);

  // insert each row into subblock
  // let FillComplete handle column distribution
  for (int localRow = 0; localRow < numMyRows; localRow++) {
    int numEntries = -1;
    int globalRow  = gRowMap.GID(localRow);
    int lid        = A.RowMap().LID(globalRow);
    int contigRow  = rowMap.GID(localRow);
    TEUCHOS_ASSERT(lid > -1);

    int err = A.ExtractMyRowCopy(lid, maxNumEntries, numEntries, &values[0], &indices[0]);
    TEUCHOS_ASSERT(err == 0);

    int numOwnedCols = 0;
    for (int localCol = 0; localCol < numEntries; localCol++) {
      TEUCHOS_ASSERT(indices[localCol] > -1);

      // if global id is not owned by this column
      int gid = local2ContigGIDs[indices[localCol]];
      if (gid == -1) continue;  // in contiguous row

      colIndices[numOwnedCols] = gid;
      colValues[numOwnedCols]  = values[localCol];
      numOwnedCols++;
    }

    // insert it into the new matrix
    mat.SumIntoGlobalValues(contigRow, numOwnedCols, &colValues[0], &colIndices[0]);
  }
}

}  // namespace Blocking
}  // namespace Epetra
}  // namespace Teko
