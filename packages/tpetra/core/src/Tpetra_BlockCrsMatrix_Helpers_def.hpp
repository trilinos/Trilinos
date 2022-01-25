// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// ************************************************************************
// @HEADER

#ifndef TPETRA_BLOCKCRSMATRIX_HELPERS_DEF_HPP
#define TPETRA_BLOCKCRSMATRIX_HELPERS_DEF_HPP

/// \file Tpetra_BlockCrsMatrix_Helpers_def.hpp

#include "Tpetra_BlockCrsMatrix.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_HashTable.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include <ctime>
#include <fstream>

namespace Tpetra {

  template<class Scalar, class LO, class GO, class Node>
  void blockCrsMatrixWriter(BlockCrsMatrix<Scalar,LO,GO,Node> const &A, std::string const &fileName) {
    Teuchos::ParameterList pl;
    std::ofstream out;
    out.open(fileName.c_str());
    blockCrsMatrixWriter(A, out, pl);
  }

  template<class Scalar, class LO, class GO, class Node>
  void blockCrsMatrixWriter(BlockCrsMatrix<Scalar,LO,GO,Node> const &A, std::string const &fileName, Teuchos::ParameterList const &params) {
    std::ofstream out;
    out.open(fileName.c_str());
    blockCrsMatrixWriter(A, out, params);
  }

  template<class Scalar, class LO, class GO, class Node>
  void blockCrsMatrixWriter(BlockCrsMatrix<Scalar,LO,GO,Node> const &A, std::ostream &os) {
    Teuchos::ParameterList pl;
    blockCrsMatrixWriter(A, os, pl);
  }

  template<class Scalar, class LO, class GO, class Node>
  void blockCrsMatrixWriter(BlockCrsMatrix<Scalar,LO,GO,Node> const &A, std::ostream &os, Teuchos::ParameterList const &params) {

    using Teuchos::RCP;
    using Teuchos::rcp;

    typedef Teuchos::OrdinalTraits<GO>                     TOT;
    typedef BlockCrsMatrix<Scalar, LO, GO, Node>           block_crs_matrix_type;
    typedef Tpetra::Import<LO, GO, Node>                   import_type;
    typedef Tpetra::Map<LO, GO, Node>                      map_type;
    typedef Tpetra::MultiVector<GO, LO, GO, Node>          mv_type;
    typedef Tpetra::CrsGraph<LO, GO, Node>                 crs_graph_type;

    RCP<const map_type> const rowMap = A.getRowMap(); //"mesh" map
    RCP<const Teuchos::Comm<int> > comm = rowMap->getComm();
    const int myRank = comm->getRank();
    const size_t numProcs = comm->getSize();

    // If true, force use of the import strip-mining infrastructure.  This is useful for debugging on one process.
    bool alwaysUseParallelAlgorithm = false;
    if (params.isParameter("always use parallel algorithm"))
      alwaysUseParallelAlgorithm = params.get<bool>("always use parallel algorithm");
    bool printMatrixMarketHeader = true;
    if (params.isParameter("print MatrixMarket header"))
      printMatrixMarketHeader = params.get<bool>("print MatrixMarket header");

    if (printMatrixMarketHeader && myRank==0) {
      std::time_t now = std::time(NULL);

      const std::string dataTypeStr =
        Teuchos::ScalarTraits<Scalar>::isComplex ? "complex" : "real";

      // Explanation of first line of file:
      // - "%%MatrixMarket" is the tag for Matrix Market format.
      // - "matrix" is what we're printing.
      // - "coordinate" means sparse (triplet format), rather than dense.
      // - "real" / "complex" is the type (in an output format sense,
      //   not in a C++ sense) of each value of the matrix.  (The
      //   other two possibilities are "integer" (not really necessary
      //   here) and "pattern" (no values, just graph).)
      os << "%%MatrixMarket matrix coordinate " << dataTypeStr << " general" << std::endl;
      os << "% time stamp: " << ctime(&now);
      os << "% written from " << numProcs << " processes" << std::endl;
      os << "% point representation of Tpetra::BlockCrsMatrix" << std::endl;
      size_t numRows = A.getGlobalNumRows();
      size_t numCols = A.getGlobalNumCols();
      os << "% " << numRows << " block rows, " << numCols << " block columns" << std::endl;
      const LO blockSize = A.getBlockSize();
      os << "% block size " << blockSize << std::endl;
      os << numRows*blockSize << " " << numCols*blockSize << " " << A.getGlobalNumEntries()*blockSize*blockSize << std::endl;
    }

    if (numProcs==1 && !alwaysUseParallelAlgorithm) {
      writeMatrixStrip(A,os,params);
    } else {
      size_t numRows = rowMap->getNodeNumElements();

      //Create source map
      RCP<const map_type> allMeshGidsMap = rcp(new map_type(TOT::invalid(), numRows, A.getIndexBase(), comm));
      //Create and populate vector of mesh GIDs corresponding to this pid's rows.
      //This vector will be imported one pid's worth of information at a time to pid 0.
      mv_type allMeshGids(allMeshGidsMap,1);
      Teuchos::ArrayRCP<GO> allMeshGidsData = allMeshGids.getDataNonConst(0);

      for (size_t i=0; i<numRows; i++)
        allMeshGidsData[i] = rowMap->getGlobalElement(i);
      allMeshGidsData = Teuchos::null;

      // Now construct a RowMatrix on PE 0 by strip-mining the rows of the input matrix A.
      size_t stripSize = allMeshGids.getGlobalLength() / numProcs;
      size_t remainder = allMeshGids.getGlobalLength() % numProcs;
      size_t curStart = 0;
      size_t curStripSize = 0;
      Teuchos::Array<GO> importMeshGidList;
      for (size_t i=0; i<numProcs; i++) {
        if (myRank==0) { // Only PE 0 does this part
          curStripSize = stripSize;
          if (i<remainder) curStripSize++; // handle leftovers
          importMeshGidList.resize(curStripSize); // Set size of vector to max needed
          for (size_t j=0; j<curStripSize; j++) importMeshGidList[j] = j + curStart + A.getIndexBase();
          curStart += curStripSize;
        }
        // The following import map should be non-trivial only on PE 0.
        TEUCHOS_TEST_FOR_EXCEPTION(myRank>0 && curStripSize!=0,
          std::runtime_error, "Tpetra::blockCrsMatrixWriter: (pid "
          << myRank << ") map size should be zero, but is " << curStripSize);
        RCP<map_type> importMeshGidMap = rcp(new map_type(TOT::invalid(), importMeshGidList(), A.getIndexBase(), comm));
        import_type gidImporter(allMeshGidsMap, importMeshGidMap);
        mv_type importMeshGids(importMeshGidMap, 1);
        importMeshGids.doImport(allMeshGids, gidImporter, INSERT);

        // importMeshGids now has a list of GIDs for the current strip of matrix rows.
        // Use these values to build another importer that will get rows of the matrix.

        // The following import map will be non-trivial only on PE 0.
        Teuchos::ArrayRCP<const GO> importMeshGidsData = importMeshGids.getData(0);
        Teuchos::Array<GO> importMeshGidsGO;
        importMeshGidsGO.reserve(importMeshGidsData.size());
        for (typename Teuchos::ArrayRCP<const GO>::size_type j=0; j<importMeshGidsData.size(); ++j)
          importMeshGidsGO.push_back(importMeshGidsData[j]);
        RCP<const map_type> importMap = rcp(new map_type(TOT::invalid(), importMeshGidsGO(), rowMap->getIndexBase(), comm) );

        import_type importer(rowMap,importMap );
        size_t numEntriesPerRow = A.getCrsGraph().getGlobalMaxNumRowEntries();
        RCP<crs_graph_type> graph = createCrsGraph(importMap,numEntriesPerRow);
        RCP<const map_type> domainMap = A.getCrsGraph().getDomainMap();
        graph->doImport(A.getCrsGraph(), importer, INSERT);
        graph->fillComplete(domainMap, importMap);

        block_crs_matrix_type importA(*graph, A.getBlockSize());
        importA.doImport(A, importer, INSERT);

        // Finally we are ready to write this strip of the matrix
        writeMatrixStrip(importA, os, params);
      }
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  void writeMatrixStrip(BlockCrsMatrix<Scalar,LO,GO,Node> const &A, std::ostream &os, Teuchos::ParameterList const &params) {
    using Teuchos::RCP;
    using map_type = Tpetra::Map<LO, GO, Node>;
    using bcrs_type = BlockCrsMatrix<Scalar,LO,GO,Node>;
    using bcrs_local_inds_host_view_type = typename bcrs_type::local_inds_host_view_type;
    using bcrs_values_host_view_type = typename bcrs_type::values_host_view_type;
    using impl_scalar_type = typename bcrs_type::impl_scalar_type;

    size_t numRows = A.getGlobalNumRows();
    RCP<const map_type> rowMap = A.getRowMap();
    RCP<const map_type> colMap = A.getColMap();
    RCP<const Teuchos::Comm<int> > comm = rowMap->getComm();
    const int myRank = comm->getRank();

    const size_t meshRowOffset = rowMap->getIndexBase();
    const size_t meshColOffset = colMap->getIndexBase();
    TEUCHOS_TEST_FOR_EXCEPTION(meshRowOffset != meshColOffset,
      std::runtime_error, "Tpetra::writeMatrixStrip: "
      "mesh row index base != mesh column index base");

    if (myRank !=0) {

      TEUCHOS_TEST_FOR_EXCEPTION(A.getNodeNumRows() != 0,
        std::runtime_error, "Tpetra::writeMatrixStrip: pid "
        << myRank << " should have 0 rows but has " << A.getNodeNumRows());
      TEUCHOS_TEST_FOR_EXCEPTION(A.getNodeNumCols() != 0,
        std::runtime_error, "Tpetra::writeMatrixStrip: pid "
        << myRank << " should have 0 columns but has " << A.getNodeNumCols());

    } else {

      TEUCHOS_TEST_FOR_EXCEPTION(numRows != A.getNodeNumRows(),
        std::runtime_error, "Tpetra::writeMatrixStrip: "
        "number of rows on pid 0 does not match global number of rows");


      int err = 0;
      const LO blockSize = A.getBlockSize();
      const size_t numLocalRows = A.getNodeNumRows();
      bool precisionChanged=false;
      int oldPrecision = 0; // avoid "unused variable" warning
      if (params.isParameter("precision")) {
        oldPrecision = os.precision(params.get<int>("precision"));
        precisionChanged=true;
      }
      int pointOffset = 1;
      if (params.isParameter("zero-based indexing")) {
        if (params.get<bool>("zero-based indexing") == true)
          pointOffset = 0;
      }

      size_t localRowInd;
      for (localRowInd = 0; localRowInd < numLocalRows; ++localRowInd) {

        // Get a view of the current row.
        bcrs_local_inds_host_view_type localColInds;
        bcrs_values_host_view_type vals;
        LO numEntries;
        A.getLocalRowView (localRowInd, localColInds, vals); numEntries = localColInds.extent(0);
        GO globalMeshRowID = rowMap->getGlobalElement(localRowInd) - meshRowOffset;

        for (LO k = 0; k < numEntries; ++k) {
          GO globalMeshColID = colMap->getGlobalElement(localColInds[k]) - meshColOffset;
          const impl_scalar_type* curBlock = vals.data() + blockSize * blockSize * k;
          // Blocks are stored in row-major format.
          for (LO j = 0; j < blockSize; ++j) {
            GO globalPointRowID = globalMeshRowID * blockSize + j + pointOffset;
            for (LO i = 0; i < blockSize; ++i) {
              GO globalPointColID = globalMeshColID * blockSize + i + pointOffset;
              const impl_scalar_type curVal = curBlock[i + j * blockSize];

              os << globalPointRowID << " " << globalPointColID << " ";
              if (Teuchos::ScalarTraits<impl_scalar_type>::isComplex) {
                // Matrix Market format wants complex values to be
                // written as space-delimited pairs.  See Bug 6469.
                os << Teuchos::ScalarTraits<impl_scalar_type>::real (curVal) << " "
                   << Teuchos::ScalarTraits<impl_scalar_type>::imag (curVal);
              }
              else {
                os << curVal;
              }
              os << std::endl;
            }
          }
        }
      }
      if (precisionChanged)
        os.precision(oldPrecision);
      TEUCHOS_TEST_FOR_EXCEPTION(err != 0,
        std::runtime_error, "Tpetra::writeMatrixStrip: "
        "error getting view of local row " << localRowInd);

    }

  }

  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<BlockCrsMatrix<Scalar, LO, GO, Node> >
  convertToBlockCrsMatrix(const Tpetra::CrsMatrix<Scalar, LO, GO, Node>& pointMatrix, const LO &blockSize)
  {

      /*
        ASSUMPTIONS:

           1) In point matrix, all entries associated with a little block are present (even if they are zero).
           2) For given mesh DOF, point DOFs appear consecutively and in ascending order in row & column maps.
           3) Point column map and block column map are ordered consistently.
      */

      using Teuchos::Array;
      using Teuchos::ArrayView;
      using Teuchos::RCP;

      typedef Tpetra::BlockCrsMatrix<Scalar,LO,GO,Node> block_crs_matrix_type;
      typedef Tpetra::Map<LO,GO,Node>                   map_type;
      typedef Tpetra::CrsGraph<LO,GO,Node>              crs_graph_type;
      typedef Tpetra::CrsMatrix<Scalar, LO,GO,Node>     crs_matrix_type;

      const map_type &pointRowMap = *(pointMatrix.getRowMap());
      RCP<const map_type> meshRowMap = createMeshMap<LO,GO,Node>(blockSize, pointRowMap);

      const map_type &pointColMap = *(pointMatrix.getColMap());
      RCP<const map_type> meshColMap = createMeshMap<LO,GO,Node>(blockSize, pointColMap);

      const map_type &pointDomainMap = *(pointMatrix.getDomainMap());
      RCP<const map_type> meshDomainMap = createMeshMap<LO,GO,Node>(blockSize, pointDomainMap);

      const map_type &pointRangeMap = *(pointMatrix.getRangeMap());
      RCP<const map_type> meshRangeMap = createMeshMap<LO,GO,Node>(blockSize, pointRangeMap);

      // Use graph ctor that provides column map and upper bound on nonzeros per row.
      // We can use static profile because the point graph should have at least as many entries per
      // row as the mesh graph.
      RCP<crs_graph_type> meshCrsGraph = rcp(new crs_graph_type(meshRowMap, meshColMap,
                                                 pointMatrix.getGlobalMaxNumRowEntries()));
      // Fill the graph by walking through the matrix.  For each mesh row, we query the collection of point
      // rows associated with it. The point column ids are converted to mesh column ids and put into an array.
      // As each point row collection is finished, the mesh column ids are sorted, made unique, and inserted
      // into the mesh graph.
      typename crs_matrix_type::local_inds_host_view_type pointColInds;
      typename crs_matrix_type::values_host_view_type pointVals;
      Array<GO> meshColGids;
      meshColGids.reserve(pointMatrix.getGlobalMaxNumRowEntries());
      //again, I assume that point GIDs associated with a mesh GID are consecutive.
      //if they are not, this will break!!
      for (size_t i=0; i<pointMatrix.getNodeNumRows()/blockSize; i++) {
        for (int j=0; j<blockSize; ++j) {
          LO rowLid = i*blockSize+j;
          pointMatrix.getLocalRowView(rowLid,pointColInds,pointVals); //TODO optimization: Since I don't care about values,
                                                                      //TODO I should use the graph instead.
          for (size_t k=0; k<pointColInds.size(); ++k) {
            GO meshColInd = pointColMap.getGlobalElement(pointColInds[k]) / blockSize;
            meshColGids.push_back(meshColInd);
          }
        }
        //List of mesh GIDs probably contains duplicates because we looped over all point rows in the block.
        //Sort and make unique.
        std::sort(meshColGids.begin(), meshColGids.end());
        meshColGids.erase( std::unique(meshColGids.begin(), meshColGids.end()), meshColGids.end() );
        meshCrsGraph->insertGlobalIndices(meshRowMap->getGlobalElement(i), meshColGids());
        meshColGids.clear();
      }
      meshCrsGraph->fillComplete(meshDomainMap,meshRangeMap);

      //create and populate the block matrix
      RCP<block_crs_matrix_type> blockMatrix = rcp(new block_crs_matrix_type(*meshCrsGraph, blockSize));

      //preallocate the maximum number of (dense) block entries needed by any row
      int maxBlockEntries = blockMatrix->getNodeMaxNumRowEntries();
      Array<Array<Scalar>> blocks(maxBlockEntries);
      for (int i=0; i<maxBlockEntries; ++i)
        blocks[i].reserve(blockSize*blockSize);
      std::map<int,int> bcol2bentry;             //maps block column index to dense block entries
      std::map<int,int>::iterator iter;
      //Fill the block matrix.  We must do this in local index space.
      //TODO: Optimization: We assume the blocks are fully populated in the point matrix.  This means
      //TODO: on the first point row in the block row, we know that we're hitting new block col indices.
      //TODO: on other rows, we know the block col indices have all been seen before
      //int offset;
      //if (pointMatrix.getIndexBase()) offset = 0;
      //else                     offset = 1;
      for (size_t i=0; i<pointMatrix.getNodeNumRows()/blockSize; i++) {
        int blkCnt=0; //how many unique block entries encountered so far in current block row
        for (int j=0; j<blockSize; ++j) {
          LO rowLid = i*blockSize+j;
          pointMatrix.getLocalRowView(rowLid,pointColInds,pointVals);
          for (size_t k=0; k<pointColInds.size(); ++k) {
            //convert point column to block col
            LO meshColInd = pointColInds[k] / blockSize;
            iter = bcol2bentry.find(meshColInd);
            if (iter == bcol2bentry.end()) {
              //new block column
              bcol2bentry[meshColInd] = blkCnt;
              blocks[blkCnt].push_back(pointVals[k]);
              blkCnt++;
            } else {
              //block column found previously
              int littleBlock = iter->second;
              blocks[littleBlock].push_back(pointVals[k]);
            }
          }
        }
        // TODO This inserts the blocks one block entry at a time.  It is probably more efficient to
        // TODO store all the blocks in a block row contiguously so they can be inserted with a single call.
        for (iter=bcol2bentry.begin(); iter != bcol2bentry.end(); ++iter) {
          LO localBlockCol = iter->first;
          Scalar *vals = (blocks[iter->second]).getRawPtr();
          blockMatrix->replaceLocalValues(i, &localBlockCol, vals, 1);
        }

        //Done with block row.  Zero everything out.
        for (int j=0; j<maxBlockEntries; ++j)
          blocks[j].clear();
        blkCnt = 0;
        bcol2bentry.clear();
      }

      return blockMatrix;

  }

  template<class LO, class GO, class Node>
  Teuchos::RCP<const Tpetra::Map<LO,GO,Node> >
  createMeshMap (const LO& blockSize, const Tpetra::Map<LO,GO,Node>& pointMap)
  {
    typedef Teuchos::OrdinalTraits<Tpetra::global_size_t> TOT;
    typedef Tpetra::Map<LO,GO,Node> map_type;

    //calculate mesh GIDs
    Teuchos::ArrayView<const GO> pointGids = pointMap.getNodeElementList();
    Teuchos::Array<GO> meshGids;
    GO indexBase = pointMap.getIndexBase();

    // Use hash table to track whether we've encountered this GID previously.  This will happen
    // when striding through the point DOFs in a block.  It should not happen otherwise.
    // I don't use sort/make unique because I don't want to change the ordering.
    meshGids.reserve(pointGids.size());
    Tpetra::Details::HashTable<GO,int> hashTable(pointGids.size());
    for (int i=0; i<pointGids.size(); ++i) {
      GO meshGid = (pointGids[i]-indexBase) / blockSize + indexBase;
      if (hashTable.get(meshGid) == -1) {
        hashTable.add(meshGid,1);   //(key,value)
        meshGids.push_back(meshGid);
      }
    }

    Teuchos::RCP<const map_type> meshMap = Teuchos::rcp( new map_type(TOT::invalid(), meshGids(), 0, pointMap.getComm()) );
    return meshMap;

  }


  template<class LO, class GO, class Node>
  Teuchos::RCP<const Tpetra::Map<LO,GO,Node> >
  createPointMap (const LO& blockSize, const Tpetra::Map<LO,GO,Node>& blockMap)
  {
    typedef Teuchos::OrdinalTraits<Tpetra::global_size_t> TOT;
    typedef Tpetra::Map<LO,GO,Node> map_type;

    //calculate mesh GIDs
    Teuchos::ArrayView<const GO> blockGids = blockMap.getNodeElementList();
    Teuchos::Array<GO> pointGids(blockGids.size() * blockSize);
    GO indexBase = blockMap.getIndexBase();

    for(LO i=0, ct=0; i<(LO)blockGids.size(); i++) {
      GO base = (blockGids[i] - indexBase)* blockSize + indexBase;
      for(LO j=0; j<blockSize; j++, ct++) 
	pointGids[i*blockSize+j] = base+j;
    }

    Teuchos::RCP<const map_type> pointMap = Teuchos::rcp( new map_type(TOT::invalid(), pointGids(), 0, blockMap.getComm()) );
    return pointMap;

  }


  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<CrsMatrix<Scalar, LO, GO, Node> >
  convertToCrsMatrix(const Tpetra::BlockCrsMatrix<Scalar, LO, GO, Node>& blockMatrix) {

    using Teuchos::Array;
    using Teuchos::ArrayView;
    using Teuchos::RCP;
    
    typedef Tpetra::BlockCrsMatrix<Scalar,LO,GO,Node> block_crs_matrix_type;
    typedef Tpetra::Map<LO,GO,Node>                   map_type;
    typedef Tpetra::CrsMatrix<Scalar, LO,GO,Node>     crs_matrix_type;    

    using crs_graph_type           = typename block_crs_matrix_type::crs_graph_type;
    using local_graph_device_type  = typename crs_matrix_type::local_graph_device_type;
    using local_matrix_device_type = typename crs_matrix_type::local_matrix_device_type;
    using row_map_type             = typename local_graph_device_type::row_map_type::non_const_type;
    using entries_type             = typename local_graph_device_type::entries_type::non_const_type;
    using values_type              = typename local_matrix_device_type::values_type::non_const_type;

    using row_map_type_const       = typename local_graph_device_type::row_map_type;
    using entries_type_const       = typename local_graph_device_type::entries_type;

    using little_block_type        = typename block_crs_matrix_type::const_little_block_type;
    using offset_type              = typename row_map_type::non_const_value_type;
    using column_type              = typename entries_type::non_const_value_type;

    using execution_space = typename Node::execution_space;
    using range_type = Kokkos::RangePolicy<execution_space, LO>;


    LO blocksize = blockMatrix.getBlockSize();
    const offset_type bs2 = blocksize * blocksize; 
    size_t block_nnz = blockMatrix.getNodeNumEntries();
    size_t point_nnz = block_nnz * bs2;

    // We can get these from the blockMatrix directly
    RCP<const map_type> pointDomainMap = blockMatrix.getDomainMap();
    RCP<const map_type> pointRangeMap  = blockMatrix.getRangeMap();
  
    // We need to generate the row/col Map ourselves.
    RCP<const map_type> blockRowMap = blockMatrix.getRowMap();
    RCP<const map_type> pointRowMap = createPointMap<LO,GO,Node>(blocksize, *blockRowMap);

    RCP<const map_type> blockColMap = blockMatrix.getColMap();
    RCP<const map_type> pointColMap = createPointMap<LO,GO,Node>(blocksize, *blockColMap);

    // Get the last few things

    const crs_graph_type & blockGraph = blockMatrix.getCrsGraph();
    LO point_rows = (LO) pointRowMap->getNodeNumElements();
    LO block_rows = (LO) blockRowMap->getNodeNumElements();
    auto blockValues = blockMatrix.getValuesDevice();
    auto blockLocalGraph = blockGraph.getLocalGraphDevice();
    row_map_type_const blockRowptr = blockLocalGraph.row_map;
    entries_type_const blockColind = blockLocalGraph.entries;

    // Generate the point matrix rowptr / colind / values
    row_map_type rowptr("row_map", point_rows+1);
    entries_type colind("entries", point_nnz);
    values_type values("values",  point_nnz);

    // Pre-fill the rowptr
    Kokkos::parallel_for("fillRowPtr",range_type(0,block_rows),KOKKOS_LAMBDA(const LO i) {
	offset_type base = blockRowptr[i];
	offset_type diff = blockRowptr[i+1] - base;
	for(LO j=0; j<blocksize; j++) {
	  rowptr[i*blocksize +j] = base*bs2 + j*diff*blocksize;
	}

	// Fill the last guy, if we're on the final entry
	if(i==block_rows-1) {
	  rowptr[block_rows*blocksize] = blockRowptr[i+1] * bs2;
	}	  
      });


    Kokkos::parallel_for("copyEntriesAndValues",range_type(0,block_rows),KOKKOS_LAMBDA(const LO i) {
	const offset_type blkBeg    = blockRowptr[i];
	const offset_type numBlocks = blockRowptr[i+1] - blkBeg;

	// For each block in the row...
	for (offset_type block=0; block < numBlocks; block++) {
	  column_type point_col_base = blockColind[blkBeg + block] * blocksize;	  
	  little_block_type my_block(blockValues.data () + (blkBeg+block) * bs2, blocksize, blocksize);
                                   
	  // For each entry in the block...
	  for(LO little_row=0; little_row<blocksize; little_row++) {
	    offset_type point_row_offset = rowptr[i*blocksize + little_row];
	    for(LO little_col=0; little_col<blocksize; little_col++) {
	      colind[point_row_offset + block*blocksize + little_col] = point_col_base + little_col;
	      values[point_row_offset + block*blocksize + little_col] = my_block(little_row,little_col);
	    }
	  }

	}
      });

    // Generate the matrix
    RCP<crs_matrix_type> pointCrsMatrix = rcp(new crs_matrix_type(pointRowMap, pointColMap, 0));
    pointCrsMatrix->setAllValues(rowptr,colind,values);

    // FUTURE OPTIMIZATION: Directly compute import/export, rather than letting ESFC do it
    pointCrsMatrix->expertStaticFillComplete(pointDomainMap,pointRangeMap);
    
    return pointCrsMatrix;
  }


} // namespace Tpetra

//
// Explicit instantiation macro for blockCrsMatrixWriter (various
// overloads), writeMatrixStrip, and convertToBlockCrsMatrix.
//
// Must be expanded from within the Tpetra namespace!
//
#define TPETRA_BLOCKCRSMATRIX_HELPERS_INSTANT(S,LO,GO,NODE) \
  template void blockCrsMatrixWriter(BlockCrsMatrix<S,LO,GO,NODE> const &A, std::string const &fileName); \
  template void blockCrsMatrixWriter(BlockCrsMatrix<S,LO,GO,NODE> const &A, std::string const &fileName, Teuchos::ParameterList const &params); \
  template void blockCrsMatrixWriter(BlockCrsMatrix<S,LO,GO,NODE> const &A, std::ostream &os); \
  template void blockCrsMatrixWriter(BlockCrsMatrix<S,LO,GO,NODE> const &A, std::ostream &os, Teuchos::ParameterList const &params); \
  template void writeMatrixStrip(BlockCrsMatrix<S,LO,GO,NODE> const &A, std::ostream &os, Teuchos::ParameterList const &params); \
  template Teuchos::RCP<BlockCrsMatrix<S, LO, GO, NODE> > convertToBlockCrsMatrix(const CrsMatrix<S, LO, GO, NODE>& pointMatrix, const LO &blockSize); \
  template Teuchos::RCP<CrsMatrix<S, LO, GO, NODE> > convertToCrsMatrix(const BlockCrsMatrix<S, LO, GO, NODE>& blockMatrix);

//
// Explicit instantiation macro for createMeshMap / createPointMap
//
#define TPETRA_CREATEMESHMAP_INSTANT(LO,GO,NODE) \
  template Teuchos::RCP<const Map<LO,GO,NODE> > createMeshMap  (const LO& blockSize, const Map<LO,GO,NODE>& pointMap); \
  template Teuchos::RCP<const Map<LO,GO,NODE> > createPointMap (const LO& blockSize, const Map<LO,GO,NODE>& blockMap);

#endif // TPETRA_BLOCKCRSMATRIX_HELPERS_DEF_HPP
