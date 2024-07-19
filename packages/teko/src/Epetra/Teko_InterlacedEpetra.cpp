// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_InterlacedEpetra.hpp"

#include <vector>

using Teuchos::RCP;
using Teuchos::rcp;

namespace Teko {
namespace Epetra {
namespace Strided {

// this assumes that there are numGlobals with numVars each interlaced
// i.e. for numVars = 2 (u,v) then the vector is
//    [u_0,v_0,u_1,v_1,u_2,v_2, ..., u_(numGlobals-1),v_(numGlobals-1)]
void buildSubMaps(int numGlobals, int numVars, const Epetra_Comm& comm,
                  std::vector<std::pair<int, RCP<Epetra_Map> > >& subMaps) {
  std::vector<int> vars;

  // build vector describing the sub maps
  for (int i = 0; i < numVars; i++) vars.push_back(1);

  // build all the submaps
  buildSubMaps(numGlobals, vars, comm, subMaps);
}

// build maps to make other conversions
void buildSubMaps(const Epetra_Map& globalMap, const std::vector<int>& vars,
                  const Epetra_Comm& comm,
                  std::vector<std::pair<int, Teuchos::RCP<Epetra_Map> > >& subMaps) {
  buildSubMaps(globalMap.NumGlobalElements(), globalMap.NumMyElements(), globalMap.MinMyGID(), vars,
               comm, subMaps);
}

// build maps to make other conversions
void buildSubMaps(int numGlobals, const std::vector<int>& vars, const Epetra_Comm& comm,
                  std::vector<std::pair<int, Teuchos::RCP<Epetra_Map> > >& subMaps) {
  std::vector<int>::const_iterator varItr;

  // compute total number of variables
  int numGlobalVars = 0;
  for (varItr = vars.begin(); varItr != vars.end(); ++varItr) numGlobalVars += *varItr;

  // must be an even number of globals
  TEUCHOS_ASSERT((numGlobals % numGlobalVars) == 0);

  Epetra_Map sampleMap(numGlobals / numGlobalVars, 0, comm);

  buildSubMaps(numGlobals, numGlobalVars * sampleMap.NumMyElements(),
               numGlobalVars * sampleMap.MinMyGID(), vars, comm, subMaps);
}

// build maps to make other conversions
void buildSubMaps(int numGlobals, int numMyElements, int minMyGID, const std::vector<int>& vars,
                  const Epetra_Comm& comm,
                  std::vector<std::pair<int, Teuchos::RCP<Epetra_Map> > >& subMaps) {
  std::vector<int>::const_iterator varItr;

  // compute total number of variables
  int numGlobalVars = 0;
  for (varItr = vars.begin(); varItr != vars.end(); ++varItr) numGlobalVars += *varItr;

  // must be an even number of globals
  TEUCHOS_ASSERT((numGlobals % numGlobalVars) == 0);
  TEUCHOS_ASSERT((numMyElements % numGlobalVars) == 0);
  TEUCHOS_ASSERT((minMyGID % numGlobalVars) == 0);

  int numBlocks  = numMyElements / numGlobalVars;
  int minBlockID = minMyGID / numGlobalVars;

  subMaps.clear();

  // index into local block in strided map
  int blockOffset = 0;
  for (varItr = vars.begin(); varItr != vars.end(); ++varItr) {
    int numLocalVars = *varItr;
    int numAllElmts  = numLocalVars * numGlobals / numGlobalVars;
    int numMyElmts   = numLocalVars * numBlocks;

    // create global arrays describing the as of yet uncreated maps
    std::vector<int> subGlobals;
    std::vector<int> contigGlobals;  // the contiguous globals

    // loop over each block of variables
    int count = 0;
    for (int blockNum = 0; blockNum < numBlocks; blockNum++) {
      // loop over each local variable in the block
      for (int local = 0; local < numLocalVars; ++local) {
        // global block number = minGID+blockNum
        // block begin global id = numGlobalVars*(minGID+blockNum)
        // global id block offset = blockOffset+local
        subGlobals.push_back((minBlockID + blockNum) * numGlobalVars + blockOffset + local);

        // also build the contiguous IDs
        contigGlobals.push_back(numLocalVars * minBlockID + count);
        count++;
      }
    }

    // sanity check
    assert((unsigned int)numMyElmts == subGlobals.size());

    // create the map with contiguous elements and the map with global elements
    RCP<Epetra_Map> subMap = rcp(new Epetra_Map(numAllElmts, numMyElmts, &subGlobals[0], 0, comm));
    RCP<Epetra_Map> contigMap =
        rcp(new Epetra_Map(numAllElmts, numMyElmts, &contigGlobals[0], 0, comm));

    Teuchos::set_extra_data(contigMap, "contigMap", Teuchos::inOutArg(subMap));
    subMaps.push_back(std::make_pair(numLocalVars, subMap));

    // update the block offset
    blockOffset += numLocalVars;
  }
}

void buildExportImport(const Epetra_Map& baseMap,
                       const std::vector<std::pair<int, RCP<Epetra_Map> > >& subMaps,
                       std::vector<RCP<Epetra_Export> >& subExport,
                       std::vector<RCP<Epetra_Import> >& subImport) {
  std::vector<std::pair<int, RCP<Epetra_Map> > >::const_iterator mapItr;

  // build importers and exporters
  for (mapItr = subMaps.begin(); mapItr != subMaps.end(); ++mapItr) {
    // exctract basic map
    const Epetra_Map& map = *(mapItr->second);

    // add new elements to vectors
    subImport.push_back(rcp(new Epetra_Import(map, baseMap)));
    subExport.push_back(rcp(new Epetra_Export(map, baseMap)));
  }
}

void buildSubVectors(const std::vector<std::pair<int, RCP<Epetra_Map> > >& subMaps,
                     std::vector<RCP<Epetra_MultiVector> >& subVectors, int count) {
  std::vector<std::pair<int, RCP<Epetra_Map> > >::const_iterator mapItr;

  // build vectors
  for (mapItr = subMaps.begin(); mapItr != subMaps.end(); ++mapItr) {
    // exctract basic map
    const Epetra_Map& map =
        *(Teuchos::get_extra_data<RCP<Epetra_Map> >(mapItr->second, "contigMap"));

    // add new elements to vectors
    RCP<Epetra_MultiVector> mv = rcp(new Epetra_MultiVector(map, count));
    Teuchos::set_extra_data(mapItr->second, "globalMap", Teuchos::inOutArg(mv));
    subVectors.push_back(mv);
  }
}

void associateSubVectors(const std::vector<std::pair<int, RCP<Epetra_Map> > >& subMaps,
                         std::vector<RCP<const Epetra_MultiVector> >& subVectors) {
  std::vector<std::pair<int, RCP<Epetra_Map> > >::const_iterator mapItr;
  std::vector<RCP<const Epetra_MultiVector> >::iterator vecItr;

  TEUCHOS_ASSERT(subMaps.size() == subVectors.size());

  // associate the sub vectors with the subMaps
  for (mapItr = subMaps.begin(), vecItr = subVectors.begin(); mapItr != subMaps.end();
       ++mapItr, ++vecItr)
    Teuchos::set_extra_data(mapItr->second, "globalMap", Teuchos::inOutArg(*vecItr));
}

// build a single subblock Epetra_CrsMatrix
RCP<Epetra_CrsMatrix> buildSubBlock(int i, int j, const Epetra_CrsMatrix& A,
                                    const std::vector<std::pair<int, RCP<Epetra_Map> > >& subMaps) {
  // get the number of variables families
  int numVarFamily = subMaps.size();

  TEUCHOS_ASSERT(i >= 0 && i < numVarFamily);
  TEUCHOS_ASSERT(j >= 0 && j < numVarFamily);

  const Epetra_Map& gRowMap = *subMaps[i].second;
  const Epetra_Map& rowMap =
      *Teuchos::get_extra_data<RCP<Epetra_Map> >(subMaps[i].second, "contigMap");
  const Epetra_Map& colMap =
      *Teuchos::get_extra_data<RCP<Epetra_Map> >(subMaps[j].second, "contigMap");
  int colFamilyCnt = subMaps[j].first;

  // compute the number of global variables
  // and the row and column block offset
  int numGlobalVars  = 0;
  int rowBlockOffset = 0;
  int colBlockOffset = 0;
  for (int k = 0; k < numVarFamily; k++) {
    numGlobalVars += subMaps[k].first;

    // compute block offsets
    if (k < i) rowBlockOffset += subMaps[k].first;
    if (k < j) colBlockOffset += subMaps[k].first;
  }

  // copy all global rows to here
  Epetra_Import import(gRowMap, A.RowMap());
  Epetra_CrsMatrix localA(Copy, gRowMap, 0);
  localA.Import(A, import, Insert);

  RCP<Epetra_CrsMatrix> mat = rcp(new Epetra_CrsMatrix(Copy, rowMap, 0));

  // get entry information
  int numMyRows     = rowMap.NumMyElements();
  int maxNumEntries = A.GlobalMaxNumEntries();

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
    int contigRow  = rowMap.GID(localRow);

    TEUCHOS_ASSERT(globalRow >= 0);
    TEUCHOS_ASSERT(contigRow >= 0);

    // extract a global row copy
    int err =
        localA.ExtractGlobalRowCopy(globalRow, maxNumEntries, numEntries, &values[0], &indices[0]);
    TEUCHOS_ASSERT(err == 0);

    int numOwnedCols = 0;
    for (int localCol = 0; localCol < numEntries; localCol++) {
      int globalCol = indices[localCol];

      // determinate which block this column ID is in
      int block = globalCol / numGlobalVars;

      bool inFamily = true;

      // test the beginning of the block
      inFamily &= (block * numGlobalVars + colBlockOffset <= globalCol);
      inFamily &= ((block * numGlobalVars + colBlockOffset + colFamilyCnt) > globalCol);

      // is this column in the variable family
      if (inFamily) {
        int familyOffset = globalCol - (block * numGlobalVars + colBlockOffset);

        colIndices[numOwnedCols] = block * colFamilyCnt + familyOffset;
        colValues[numOwnedCols]  = values[localCol];

        numOwnedCols++;
      }
    }

    // insert it into the new matrix
    mat->InsertGlobalValues(contigRow, numOwnedCols, &colValues[0], &colIndices[0]);
  }

  // fill it and automagically optimize the storage
  mat->FillComplete(colMap, rowMap);

  return mat;
}

// rebuild a single subblock Epetra_CrsMatrix
void rebuildSubBlock(int i, int j, const Epetra_CrsMatrix& A,
                     const std::vector<std::pair<int, RCP<Epetra_Map> > >& subMaps,
                     Epetra_CrsMatrix& mat) {
  // get the number of variables families
  int numVarFamily = subMaps.size();

  TEUCHOS_ASSERT(i >= 0 && i < numVarFamily);
  TEUCHOS_ASSERT(j >= 0 && j < numVarFamily);
  TEUCHOS_ASSERT(mat.Filled());

  const Epetra_Map& gRowMap = *subMaps[i].second;
  const Epetra_Map& rowMap =
      *Teuchos::get_extra_data<RCP<Epetra_Map> >(subMaps[i].second, "contigMap");
  int colFamilyCnt = subMaps[j].first;

  // compute the number of global variables
  // and the row and column block offset
  int numGlobalVars  = 0;
  int rowBlockOffset = 0;
  int colBlockOffset = 0;
  for (int k = 0; k < numVarFamily; k++) {
    numGlobalVars += subMaps[k].first;

    // compute block offsets
    if (k < i) rowBlockOffset += subMaps[k].first;
    if (k < j) colBlockOffset += subMaps[k].first;
  }

  // copy all global rows to here
  Epetra_Import import(gRowMap, A.RowMap());
  Epetra_CrsMatrix localA(Copy, gRowMap, 0);
  localA.Import(A, import, Insert);

  // clear out the old matrix
  mat.PutScalar(0.0);

  // get entry information
  int numMyRows     = rowMap.NumMyElements();
  int maxNumEntries = A.GlobalMaxNumEntries();

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
    int contigRow  = rowMap.GID(localRow);

    TEUCHOS_ASSERT(globalRow >= 0);
    TEUCHOS_ASSERT(contigRow >= 0);

    // extract a global row copy
    int err =
        localA.ExtractGlobalRowCopy(globalRow, maxNumEntries, numEntries, &values[0], &indices[0]);
    TEUCHOS_ASSERT(err == 0);

    int numOwnedCols = 0;
    for (int localCol = 0; localCol < numEntries; localCol++) {
      int globalCol = indices[localCol];

      // determinate which block this column ID is in
      int block = globalCol / numGlobalVars;

      bool inFamily = true;

      // test the beginning of the block
      inFamily &= (block * numGlobalVars + colBlockOffset <= globalCol);
      inFamily &= ((block * numGlobalVars + colBlockOffset + colFamilyCnt) > globalCol);

      // is this column in the variable family
      if (inFamily) {
        int familyOffset = globalCol - (block * numGlobalVars + colBlockOffset);

        colIndices[numOwnedCols] = block * colFamilyCnt + familyOffset;
        colValues[numOwnedCols]  = values[localCol];

        numOwnedCols++;
      }
    }

    // insert it into the new matrix
    mat.SumIntoGlobalValues(contigRow, numOwnedCols, &colValues[0], &colIndices[0]);
  }
}

// collect subvectors into a single global vector
void many2one(Epetra_MultiVector& one, const std::vector<RCP<const Epetra_MultiVector> >& many,
              const std::vector<RCP<Epetra_Export> >& subExport) {
  // std::vector<RCP<const Epetra_Vector> >::const_iterator vecItr;
  std::vector<RCP<const Epetra_MultiVector> >::const_iterator vecItr;
  std::vector<RCP<Epetra_Export> >::const_iterator expItr;

  // using Exporters fill the empty vector from the sub-vectors
  for (vecItr = many.begin(), expItr = subExport.begin(); vecItr != many.end();
       ++vecItr, ++expItr) {
    // for ease of access to the source
    RCP<const Epetra_MultiVector> srcVec = *vecItr;

    // extract the map with global indicies from the current vector
    const Epetra_Map& globalMap = *(Teuchos::get_extra_data<RCP<Epetra_Map> >(srcVec, "globalMap"));

    // build the export vector as a view of the destination
    Epetra_MultiVector exportVector(View, globalMap, srcVec->Values(), srcVec->Stride(),
                                    srcVec->NumVectors());
    one.Export(exportVector, **expItr, Insert);
  }
}

// distribute one global vector into a many subvectors
void one2many(std::vector<RCP<Epetra_MultiVector> >& many, const Epetra_MultiVector& single,
              const std::vector<RCP<Epetra_Import> >& subImport) {
  // std::vector<RCP<Epetra_Vector> >::const_iterator vecItr;
  std::vector<RCP<Epetra_MultiVector> >::const_iterator vecItr;
  std::vector<RCP<Epetra_Import> >::const_iterator impItr;

  // using Importers fill the sub vectors from the mama vector
  for (vecItr = many.begin(), impItr = subImport.begin(); vecItr != many.end();
       ++vecItr, ++impItr) {
    // for ease of access to the destination
    RCP<Epetra_MultiVector> destVec = *vecItr;

    // extract the map with global indicies from the current vector
    const Epetra_Map& globalMap =
        *(Teuchos::get_extra_data<RCP<Epetra_Map> >(destVec, "globalMap"));

    // build the import vector as a view on the destination
    Epetra_MultiVector importVector(View, globalMap, destVec->Values(), destVec->Stride(),
                                    destVec->NumVectors());

    // perform the import
    importVector.Import(single, **impItr, Insert);
  }
}

}  // namespace Strided
}  // end namespace Epetra
}  // end namespace Teko
