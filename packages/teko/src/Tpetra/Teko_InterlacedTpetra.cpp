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
#include "Tpetra_Details_makeColMap_decl.hpp"
#include "KokkosSparse_SortCrs.hpp"

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

// build a single subblock Tpetra::CrsMatrix
RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > buildSubBlock(
    int i, int j, const RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> >& A,
    const std::vector<std::pair<int, RCP<Tpetra::Map<LO, GO, NT> > > >& subMaps) {
  // get the number of variables families
  int numVarFamily = subMaps.size();

  TEUCHOS_ASSERT(i >= 0 && i < numVarFamily);
  TEUCHOS_ASSERT(j >= 0 && j < numVarFamily);

  const Tpetra::Map<LO, GO, NT>& gRowMap = *subMaps[i].second;
  const RCP<const Tpetra::Map<LO, GO, NT> > rowMap =
      Teuchos::get_extra_data<RCP<Tpetra::Map<LO, GO, NT> > >(subMaps[i].second, "contigMap");
  const RCP<const Tpetra::Map<LO, GO, NT> > domainMap =
      Teuchos::get_extra_data<RCP<Tpetra::Map<LO, GO, NT> > >(subMaps[j].second, "contigMap");
  const RCP<const Tpetra::Map<LO, GO, NT> > rangeMap = rowMap;
  GO colFamilyCnt                                    = subMaps[j].first;

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

  // Build the sub-block on the device, mirroring the Blocking path.
  //
  // The sub-block's rows are a subset of A's rows owned by this same process,
  // so we read A's local device matrix directly (mapping sub-block global row
  // -> A local row) instead of importing all rows into a temporary CrsMatrix.
  // The interlaced column-membership test and contiguous-column renumbering
  // are pure arithmetic, so the entire extraction runs on the device with no
  // host<->device transfers and no host-side insertGlobalValues loop.
  LO numMyRows = rowMap->getLocalNumElements();

  using local_matrix_type      = Tpetra::CrsMatrix<ST, LO, GO, NT>::local_matrix_device_type;
  using row_map_type           = local_matrix_type::row_map_type::non_const_type;
  using values_type            = local_matrix_type::values_type::non_const_type;
  using index_type             = local_matrix_type::index_type::non_const_type;
  using matrix_execution_space = typename local_matrix_type::execution_space;
  using device_type            = typename NT::device_type;

  auto A_dev        = A->getLocalMatrixDevice();
  auto gRowMap_dev  = gRowMap.getLocalMap();
  auto A_rowmap_dev = A->getRowMap()->getLocalMap();
  auto A_colmap_dev = A->getColMap()->getLocalMap();

  // Count the entries owned by this sub-block in each row and build the
  // row-pointer prefix sum in a single scan.
  auto prefixSumEntriesPerRow = row_map_type(
      Kokkos::ViewAllocateWithoutInitializing("prefixSumEntriesPerRow"), numMyRows + 1);

  LO totalNumOwnedCols = 0;
  Kokkos::parallel_scan(
      Kokkos::RangePolicy<Kokkos::Schedule<Kokkos::Dynamic>, matrix_execution_space>(0, numMyRows),
      KOKKOS_LAMBDA(const LO localRow, LO& sumNumEntries, bool finalPass) {
        GO globalRow             = gRowMap_dev.getGlobalElement(localRow);
        LO lid                   = A_rowmap_dev.getLocalElement(globalRow);
        const auto sparseRowView = A_dev.row(lid);

        LO numOwnedCols = 0;
        for (auto localCol = 0; localCol < sparseRowView.length; localCol++) {
          GO globalCol  = A_colmap_dev.getGlobalElement(sparseRowView.colidx(localCol));
          GO block      = globalCol / numGlobalVars;
          bool inFamily = (block * numGlobalVars + colBlockOffset <= globalCol) &&
                          ((block * numGlobalVars + colBlockOffset + colFamilyCnt) > globalCol);
          if (inFamily) numOwnedCols++;
        }

        if (finalPass) {
          prefixSumEntriesPerRow(localRow) = sumNumEntries;
          if (localRow == (numMyRows - 1))
            prefixSumEntriesPerRow(numMyRows) = sumNumEntries + numOwnedCols;
        }
        sumNumEntries += numOwnedCols;
      },
      totalNumOwnedCols);

  auto columnIndices = Kokkos::View<GO*, device_type>(
      Kokkos::ViewAllocateWithoutInitializing("columnIndices"), totalNumOwnedCols);
  auto values = values_type(Kokkos::ViewAllocateWithoutInitializing("values"), totalNumOwnedCols);

  // Extract the contiguous column GIDs and values for each row.
  LO maxNumEntriesSubblock = 0;
  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<Kokkos::Schedule<Kokkos::Dynamic>, matrix_execution_space>(0, numMyRows),
      KOKKOS_LAMBDA(const LO localRow, LO& maxNumEntries) {
        GO globalRow             = gRowMap_dev.getGlobalElement(localRow);
        LO lid                   = A_rowmap_dev.getLocalElement(globalRow);
        const auto sparseRowView = A_dev.row(lid);

        LO colId      = 0;
        LO colIdStart = prefixSumEntriesPerRow[localRow];
        for (auto localCol = 0; localCol < sparseRowView.length; localCol++) {
          GO globalCol  = A_colmap_dev.getGlobalElement(sparseRowView.colidx(localCol));
          GO block      = globalCol / numGlobalVars;
          bool inFamily = (block * numGlobalVars + colBlockOffset <= globalCol) &&
                          ((block * numGlobalVars + colBlockOffset + colFamilyCnt) > globalCol);
          if (!inFamily) continue;

          GO familyOffset                   = globalCol - (block * numGlobalVars + colBlockOffset);
          columnIndices(colId + colIdStart) = block * colFamilyCnt + familyOffset;
          values(colId + colIdStart)        = sparseRowView.value(localCol);
          colId++;
        }
        maxNumEntries = Kokkos::max(maxNumEntries, colId);
      },
      Kokkos::Max<LO>(maxNumEntriesSubblock));

  // Build the column map from the contiguous column GIDs, convert to local
  // column indices, and sort each row.
  Teuchos::RCP<const Tpetra::Map<LO, GO, NT> > colMap;
  Tpetra::Details::makeColMap<LO, GO, NT>(colMap, domainMap, columnIndices);
  TEUCHOS_ASSERT(colMap);

  auto colMap_dev = colMap->getLocalMap();
  auto localColumnIndices =
      index_type(Kokkos::ViewAllocateWithoutInitializing("localColumnIndices"), totalNumOwnedCols);
  Kokkos::parallel_for(
      Kokkos::RangePolicy<Kokkos::Schedule<Kokkos::Dynamic>, matrix_execution_space>(
          0, totalNumOwnedCols),
      KOKKOS_LAMBDA(const LO index) {
        localColumnIndices(index) = colMap_dev.getLocalElement(columnIndices(index));
      });

  KokkosSparse::sort_crs_matrix<matrix_execution_space, row_map_type, index_type, values_type>(
      prefixSumEntriesPerRow, localColumnIndices, values);

  auto lcl_mat = Tpetra::CrsMatrix<ST, LO, GO, NT>::local_matrix_device_type(
      "localMat", numMyRows, maxNumEntriesSubblock, totalNumOwnedCols, values,
      prefixSumEntriesPerRow, localColumnIndices);

  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > mat =
      rcp(new Tpetra::CrsMatrix<ST, LO, GO, NT>(lcl_mat, rowMap, colMap, domainMap, rangeMap));

  return mat;
}

// rebuild a single subblock Tpetra::CrsMatrix
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
  GO colFamilyCnt = subMaps[j].first;

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

  // clear out the old matrix
  mat.resumeFill();
  mat.setAllToScalar(0.0);

  // get entry information
  LO numMyRows = rowMap.getLocalNumElements();

  // Perform the re-assembly on the device, mirroring the Blocking path.
  //
  // The sub-block's rows are a subset of A's rows and are owned by this same
  // process, so we can read them directly out of A's local device matrix
  // (mapping sub-block global row -> A local row) instead of doing a redundant
  // doImport into a temporary CrsMatrix on every rebuild.  The interlaced
  // column-membership test and the contiguous-column renumbering are pure
  // arithmetic, so the whole loop runs in a single parallel_for with no
  // host<->device transfers of A's values (which is what made the old,
  // getGlobalRowCopy/sumIntoGlobalValues host loop expensive on GPU builds).
  using matrix_execution_space =
      typename Tpetra::CrsMatrix<ST, LO, GO, NT>::local_matrix_device_type::execution_space;

  auto A_dev         = A->getLocalMatrixDevice();
  auto mat_dev       = mat.getLocalMatrixDevice();
  auto gRowMap_dev   = gRowMap.getLocalMap();
  auto A_rowmap_dev  = A->getRowMap()->getLocalMap();
  auto A_colmap_dev  = A->getColMap()->getLocalMap();
  auto matColMap_dev = mat.getColMap()->getLocalMap();

  const auto invalidLO = Teuchos::OrdinalTraits<LO>::invalid();

  Kokkos::parallel_for(
      Kokkos::RangePolicy<Kokkos::Schedule<Kokkos::Dynamic>, matrix_execution_space>(0, numMyRows),
      KOKKOS_LAMBDA(const LO localRow) {
        GO globalRow             = gRowMap_dev.getGlobalElement(localRow);
        LO lid                   = A_rowmap_dev.getLocalElement(globalRow);
        const auto sparseRowView = A_dev.row(lid);

        for (auto localCol = 0; localCol < sparseRowView.length; localCol++) {
          GO globalCol = A_colmap_dev.getGlobalElement(sparseRowView.colidx(localCol));

          // determine which block this column ID is in
          GO block = globalCol / numGlobalVars;

          // is this column in the variable family
          bool inFamily = (block * numGlobalVars + colBlockOffset <= globalCol) &&
                          ((block * numGlobalVars + colBlockOffset + colFamilyCnt) > globalCol);
          if (!inFamily) continue;

          GO familyOffset = globalCol - (block * numGlobalVars + colBlockOffset);
          GO contigCol    = block * colFamilyCnt + familyOffset;

          LO lidCol = matColMap_dev.getLocalElement(contigCol);
          if (lidCol == invalidLO) continue;

          auto value = sparseRowView.value(localCol);
          mat_dev.sumIntoValues(localRow, &lidCol, 1, &value, true, false);
        }
      });

  mat.fillComplete(rcpFromRef(colMap), rcpFromRef(rowMap));
}

// collect subvectors into a single global vector
void many2one(Tpetra::MultiVector<ST, LO, GO, NT>& one,
              const std::vector<RCP<const Tpetra::MultiVector<ST, LO, GO, NT> > >& many,
              const std::vector<RCP<Tpetra::Export<LO, GO, NT> > >& subExport) {
  // std::vector<RCP<const Tpetra::Vector> >::const_iterator vecItr;
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
  // std::vector<RCP<Tpetra::Vector> >::const_iterator vecItr;
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
