// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_InterlacedTpetra.hpp"
#include "Tpetra_Import.hpp"

#include <vector>

using Teuchos::RCP;
using Teuchos::rcp;

namespace Teko {
namespace TpetraHelpers {
namespace Strided {

// this assumes that there are numGlobals with numVars each interlaced
// i.e. for numVars = 2 (u,v) then the vector is
//    [u_0,v_0,u_1,v_1,u_2,v_2, ..., u_(numGlobals-1),v_(numGlobals-1)]
void buildSubMaps(GO numGlobals, int numVars, const Teuchos::Comm<int>& comm,
                  std::vector<std::pair<int, RCP<Tpetra::Map<LO, GO, NT> > > >& subMaps) {
  std::vector<int> vars;

  // build vector describing the sub maps
  for (int i = 0; i < numVars; i++) vars.push_back(1);

  // build all the submaps
  buildSubMaps(numGlobals, vars, comm, subMaps);
}

// build maps to make other conversions
void buildSubMaps(const Tpetra::Map<LO, GO, NT>& globalMap, const std::vector<int>& vars,
                  const Teuchos::Comm<int>& comm,
                  std::vector<std::pair<int, Teuchos::RCP<Tpetra::Map<LO, GO, NT> > > >& subMaps) {
  buildSubMaps(globalMap.getGlobalNumElements(), globalMap.getLocalNumElements(),
               globalMap.getMinGlobalIndex(), vars, comm, subMaps);
}

// build maps to make other conversions
void buildSubMaps(GO numGlobals, const std::vector<int>& vars, const Teuchos::Comm<int>& comm,
                  std::vector<std::pair<int, Teuchos::RCP<Tpetra::Map<LO, GO, NT> > > >& subMaps) {
  std::vector<int>::const_iterator varItr;

  // compute total number of variables
  int numGlobalVars = 0;
  for (varItr = vars.begin(); varItr != vars.end(); ++varItr) numGlobalVars += *varItr;

  // must be an even number of globals
  TEUCHOS_ASSERT((numGlobals % numGlobalVars) == 0);

  Tpetra::Map<LO, GO, NT> sampleMap(numGlobals / numGlobalVars, 0, rcpFromRef(comm));

  buildSubMaps(numGlobals, numGlobalVars * sampleMap.getLocalNumElements(),
               numGlobalVars * sampleMap.getMinGlobalIndex(), vars, comm, subMaps);
}

// build maps to make other conversions
void buildSubMaps(GO numGlobals, LO numMyElements, GO minMyGID, const std::vector<int>& vars,
                  const Teuchos::Comm<int>& comm,
                  std::vector<std::pair<int, Teuchos::RCP<Tpetra::Map<LO, GO, NT> > > >& subMaps) {
  std::vector<int>::const_iterator varItr;

  // compute total number of variables
  int numGlobalVars = 0;
  for (varItr = vars.begin(); varItr != vars.end(); ++varItr) numGlobalVars += *varItr;

  // must be an even number of globals
  TEUCHOS_ASSERT((numGlobals % numGlobalVars) == 0);
  TEUCHOS_ASSERT((numMyElements % numGlobalVars) == 0);
  TEUCHOS_ASSERT((minMyGID % numGlobalVars) == 0);

  LO numBlocks  = numMyElements / numGlobalVars;
  GO minBlockID = minMyGID / numGlobalVars;

  subMaps.clear();

  // index into local block in strided map
  GO blockOffset = 0;
  for (varItr = vars.begin(); varItr != vars.end(); ++varItr) {
    LO numLocalVars = *varItr;
    GO numAllElmts  = numLocalVars * numGlobals / numGlobalVars;
#ifndef NDEBUG
    LO numMyElmts = numLocalVars * numBlocks;
#endif

    // create global arrays describing the as of yet uncreated maps
    std::vector<GO> subGlobals;
    std::vector<GO> contigGlobals;  // the contiguous globals

    // loop over each block of variables
    LO count = 0;
    for (LO blockNum = 0; blockNum < numBlocks; blockNum++) {
      // loop over each local variable in the block
      for (LO local = 0; local < numLocalVars; ++local) {
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
    assert((size_t)numMyElmts == subGlobals.size());

    // create the map with contiguous elements and the map with global elements
    RCP<Tpetra::Map<LO, GO, NT> > subMap    = rcp(new Tpetra::Map<LO, GO, NT>(
        numAllElmts, Teuchos::ArrayView<GO>(subGlobals), 0, rcpFromRef(comm)));
    RCP<Tpetra::Map<LO, GO, NT> > contigMap = rcp(new Tpetra::Map<LO, GO, NT>(
        numAllElmts, Teuchos::ArrayView<GO>(contigGlobals), 0, rcpFromRef(comm)));

    Teuchos::set_extra_data(contigMap, "contigMap", Teuchos::inOutArg(subMap));
    subMaps.push_back(std::make_pair(numLocalVars, subMap));

    // update the block offset
    blockOffset += numLocalVars;
  }
}

void buildExportImport(const Tpetra::Map<LO, GO, NT>& baseMap,
                       const std::vector<std::pair<int, RCP<Tpetra::Map<LO, GO, NT> > > >& subMaps,
                       std::vector<RCP<Tpetra::Export<LO, GO, NT> > >& subExport,
                       std::vector<RCP<Tpetra::Import<LO, GO, NT> > >& subImport) {
  std::vector<std::pair<int, RCP<Tpetra::Map<LO, GO, NT> > > >::const_iterator mapItr;

  // build importers and exporters
  for (mapItr = subMaps.begin(); mapItr != subMaps.end(); ++mapItr) {
    // exctract basic map
    const Tpetra::Map<LO, GO, NT>& map = *(mapItr->second);

    // add new elements to vectors
    subImport.push_back(rcp(new Tpetra::Import<LO, GO, NT>(rcpFromRef(baseMap), rcpFromRef(map))));
    subExport.push_back(rcp(new Tpetra::Export<LO, GO, NT>(rcpFromRef(map), rcpFromRef(baseMap))));
  }
}

void buildSubVectors(const std::vector<std::pair<int, RCP<Tpetra::Map<LO, GO, NT> > > >& subMaps,
                     std::vector<RCP<Tpetra::MultiVector<ST, LO, GO, NT> > >& subVectors,
                     int count) {
  std::vector<std::pair<int, RCP<Tpetra::Map<LO, GO, NT> > > >::const_iterator mapItr;

  // build vectors
  for (mapItr = subMaps.begin(); mapItr != subMaps.end(); ++mapItr) {
    // exctract basic map
    const Tpetra::Map<LO, GO, NT>& map =
        *(Teuchos::get_extra_data<RCP<Tpetra::Map<LO, GO, NT> > >(mapItr->second, "contigMap"));

    // add new elements to vectors
    RCP<Tpetra::MultiVector<ST, LO, GO, NT> > mv =
        rcp(new Tpetra::MultiVector<ST, LO, GO, NT>(rcpFromRef(map), count));
    Teuchos::set_extra_data(mapItr->second, "globalMap", Teuchos::inOutArg(mv));
    subVectors.push_back(mv);
  }
}

void associateSubVectors(
    const std::vector<std::pair<int, RCP<Tpetra::Map<LO, GO, NT> > > >& subMaps,
    std::vector<RCP<const Tpetra::MultiVector<ST, LO, GO, NT> > >& subVectors) {
  std::vector<std::pair<int, RCP<Tpetra::Map<LO, GO, NT> > > >::const_iterator mapItr;
  std::vector<RCP<const Tpetra::MultiVector<ST, LO, GO, NT> > >::iterator vecItr;

  TEUCHOS_ASSERT(subMaps.size() == subVectors.size());

  // associate the sub vectors with the subMaps
  for (mapItr = subMaps.begin(), vecItr = subVectors.begin(); mapItr != subMaps.end();
       ++mapItr, ++vecItr)
    Teuchos::set_extra_data(mapItr->second, "globalMap", Teuchos::inOutArg(*vecItr),
                            Teuchos::POST_DESTROY, false);
}

// build a single subblock Epetra_CrsMatrix
RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > buildSubBlock(
    int i, int j, const RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> >& A,
    const std::vector<std::pair<int, RCP<Tpetra::Map<LO, GO, NT> > > >& subMaps) {
  // get the number of variables families
  int numVarFamily = subMaps.size();

  TEUCHOS_ASSERT(i >= 0 && i < numVarFamily);
  TEUCHOS_ASSERT(j >= 0 && j < numVarFamily);

  const Tpetra::Map<LO, GO, NT>& gRowMap = *subMaps[i].second;
  const Tpetra::Map<LO, GO, NT>& rowMap =
      *Teuchos::get_extra_data<RCP<Tpetra::Map<LO, GO, NT> > >(subMaps[i].second, "contigMap");
  const Tpetra::Map<LO, GO, NT>& colMap =
      *Teuchos::get_extra_data<RCP<Tpetra::Map<LO, GO, NT> > >(subMaps[j].second, "contigMap");
  int colFamilyCnt = subMaps[j].first;

  // compute the number of global variables
  // and the row and column block offset
  GO numGlobalVars  = 0;
  GO rowBlockOffset = 0;
  GO colBlockOffset = 0;
  for (int k = 0; k < numVarFamily; k++) {
    numGlobalVars += subMaps[k].first;

    // compute block offsets
    if (k < i) rowBlockOffset += subMaps[k].first;
    if (k < j) colBlockOffset += subMaps[k].first;
  }

  // copy all global rows to here
  Tpetra::Import<LO, GO, NT> import(A->getRowMap(), rcpFromRef(gRowMap));
  Tpetra::CrsMatrix<ST, LO, GO, NT> localA(rcpFromRef(gRowMap), 0);
  localA.doImport(*A, import, Tpetra::INSERT);

  // get entry information
  LO numMyRows     = rowMap.getLocalNumElements();
  LO maxNumEntries = A->getGlobalMaxNumRowEntries();

  // for extraction
  auto indices = typename Tpetra::CrsMatrix<ST, LO, GO, NT>::nonconst_global_inds_host_view_type(
      Kokkos::ViewAllocateWithoutInitializing("rowIndices"), maxNumEntries);
  auto values = typename Tpetra::CrsMatrix<ST, LO, GO, NT>::nonconst_values_host_view_type(
      Kokkos::ViewAllocateWithoutInitializing("rowIndices"), maxNumEntries);

  // for counting row sizes
  std::vector<size_t> numEntriesPerRow(numMyRows, 0);

  const size_t invalid = Teuchos::OrdinalTraits<size_t>::invalid();

  // Count the sizes of each row, using same logic as insertion below
  for (LO localRow = 0; localRow < numMyRows; localRow++) {
    size_t numEntries = invalid;
    GO globalRow      = gRowMap.getGlobalElement(localRow);
    GO contigRow      = rowMap.getGlobalElement(localRow);

    TEUCHOS_ASSERT(globalRow >= 0);
    TEUCHOS_ASSERT(contigRow >= 0);

    // extract a global row copy
    localA.getGlobalRowCopy(globalRow, indices, values, numEntries);
    LO numOwnedCols = 0;
    for (size_t localCol = 0; localCol < numEntries; localCol++) {
      GO globalCol = indices(localCol);

      // determinate which block this column ID is in
      int block = globalCol / numGlobalVars;

      bool inFamily = true;

      // test the beginning of the block
      inFamily &= (block * numGlobalVars + colBlockOffset <= globalCol);
      inFamily &= ((block * numGlobalVars + colBlockOffset + colFamilyCnt) > globalCol);

      // is this column in the variable family
      if (inFamily) {
        numOwnedCols++;
      }
    }
    numEntriesPerRow[localRow] += numOwnedCols;
  }

  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > mat = rcp(new Tpetra::CrsMatrix<ST, LO, GO, NT>(
      rcpFromRef(rowMap), Teuchos::ArrayView<const size_t>(numEntriesPerRow)));

  // for insertion
  std::vector<GO> colIndices(maxNumEntries);
  std::vector<ST> colValues(maxNumEntries);

  // insert each row into subblock
  // let FillComplete handle column distribution
  for (LO localRow = 0; localRow < numMyRows; localRow++) {
    size_t numEntries = invalid;
    GO globalRow      = gRowMap.getGlobalElement(localRow);
    GO contigRow      = rowMap.getGlobalElement(localRow);

    TEUCHOS_ASSERT(globalRow >= 0);
    TEUCHOS_ASSERT(contigRow >= 0);

    // extract a global row copy
    localA.getGlobalRowCopy(globalRow, indices, values, numEntries);
    LO numOwnedCols = 0;
    for (size_t localCol = 0; localCol < numEntries; localCol++) {
      GO globalCol = indices(localCol);

      // determinate which block this column ID is in
      int block = globalCol / numGlobalVars;

      bool inFamily = true;

      // test the beginning of the block
      inFamily &= (block * numGlobalVars + colBlockOffset <= globalCol);
      inFamily &= ((block * numGlobalVars + colBlockOffset + colFamilyCnt) > globalCol);

      // is this column in the variable family
      if (inFamily) {
        GO familyOffset = globalCol - (block * numGlobalVars + colBlockOffset);

        colIndices[numOwnedCols] = block * colFamilyCnt + familyOffset;
        colValues[numOwnedCols]  = values(localCol);

        numOwnedCols++;
      }
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
  mat->fillComplete(rcpFromRef(colMap), rcpFromRef(rowMap));

  return mat;
}

// rebuild a single subblock Epetra_CrsMatrix
void rebuildSubBlock(int i, int j, const RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> >& A,
                     const std::vector<std::pair<int, RCP<Tpetra::Map<LO, GO, NT> > > >& subMaps,
                     Tpetra::CrsMatrix<ST, LO, GO, NT>& mat) {
  // get the number of variables families
  int numVarFamily = subMaps.size();

  TEUCHOS_ASSERT(i >= 0 && i < numVarFamily);
  TEUCHOS_ASSERT(j >= 0 && j < numVarFamily);
  TEUCHOS_ASSERT(mat.isFillComplete());

  const Tpetra::Map<LO, GO, NT>& gRowMap = *subMaps[i].second;
  const Tpetra::Map<LO, GO, NT>& rowMap =
      *Teuchos::get_extra_data<RCP<Tpetra::Map<LO, GO, NT> > >(subMaps[i].second, "contigMap");
  const Tpetra::Map<LO, GO, NT>& colMap =
      *Teuchos::get_extra_data<RCP<Tpetra::Map<LO, GO, NT> > >(subMaps[j].second, "contigMap");
  int colFamilyCnt = subMaps[j].first;

  // compute the number of global variables
  // and the row and column block offset
  GO numGlobalVars  = 0;
  GO rowBlockOffset = 0;
  GO colBlockOffset = 0;
  for (int k = 0; k < numVarFamily; k++) {
    numGlobalVars += subMaps[k].first;

    // compute block offsets
    if (k < i) rowBlockOffset += subMaps[k].first;
    if (k < j) colBlockOffset += subMaps[k].first;
  }

  // copy all global rows to here
  Tpetra::Import<LO, GO, NT> import(A->getRowMap(), rcpFromRef(gRowMap));
  Tpetra::CrsMatrix<ST, LO, GO, NT> localA(rcpFromRef(gRowMap), 0);
  localA.doImport(*A, import, Tpetra::INSERT);

  // clear out the old matrix
  mat.resumeFill();
  mat.setAllToScalar(0.0);

  // get entry information
  LO numMyRows     = rowMap.getLocalNumElements();
  GO maxNumEntries = A->getGlobalMaxNumRowEntries();

  // for extraction
  auto indices = typename Tpetra::CrsMatrix<ST, LO, GO, NT>::nonconst_global_inds_host_view_type(
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
    GO contigRow      = rowMap.getGlobalElement(localRow);

    TEUCHOS_ASSERT(globalRow >= 0);
    TEUCHOS_ASSERT(contigRow >= 0);

    // extract a global row copy
    localA.getGlobalRowCopy(globalRow, indices, values, numEntries);

    LO numOwnedCols = 0;
    for (size_t localCol = 0; localCol < numEntries; localCol++) {
      GO globalCol = indices(localCol);

      // determinate which block this column ID is in
      int block = globalCol / numGlobalVars;

      bool inFamily = true;

      // test the beginning of the block
      inFamily &= (block * numGlobalVars + colBlockOffset <= globalCol);
      inFamily &= ((block * numGlobalVars + colBlockOffset + colFamilyCnt) > globalCol);

      // is this column in the variable family
      if (inFamily) {
        GO familyOffset = globalCol - (block * numGlobalVars + colBlockOffset);

        colIndices[numOwnedCols] = block * colFamilyCnt + familyOffset;
        colValues[numOwnedCols]  = values(localCol);

        numOwnedCols++;
      }
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

// collect subvectors into a single global vector
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
    const Tpetra::Map<LO, GO, NT>& globalMap =
        *(Teuchos::get_extra_data<RCP<Tpetra::Map<LO, GO, NT> > >(srcVec, "globalMap"));

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

// distribute one global vector into a many subvectors
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
    const Tpetra::Map<LO, GO, NT>& globalMap =
        *(Teuchos::get_extra_data<RCP<Tpetra::Map<LO, GO, NT> > >(destVec, "globalMap"));

    // build the import vector as a view on the destination
    GO destLDA  = destVec->getStride();
    GO destSize = destVec->getGlobalLength() * destVec->getNumVectors();
    std::vector<ST> destArray(destSize);
    Teuchos::ArrayView<ST> destVals(destArray);
    destVec->get1dCopy(destVals, destLDA);
    Tpetra::MultiVector<ST, LO, GO, NT> importVector(rcpFromRef(globalMap), destVals, destLDA,
                                                     destVec->getNumVectors());

    // perform the import
    importVector.doImport(single, **impItr, Tpetra::INSERT);

    Tpetra::Import<LO, GO, NT> importer(destVec->getMap(), destVec->getMap());
    importVector.replaceMap(destVec->getMap());
    destVec->doImport(importVector, importer, Tpetra::INSERT);
  }
}

}  // namespace Strided
}  // namespace TpetraHelpers
}  // end namespace Teko
