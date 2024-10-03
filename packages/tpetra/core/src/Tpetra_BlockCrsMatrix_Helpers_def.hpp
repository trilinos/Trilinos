// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
      size_t numRows = rowMap->getLocalNumElements();

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

      TEUCHOS_TEST_FOR_EXCEPTION(A.getLocalNumRows() != 0,
        std::runtime_error, "Tpetra::writeMatrixStrip: pid "
        << myRank << " should have 0 rows but has " << A.getLocalNumRows());
      TEUCHOS_TEST_FOR_EXCEPTION(A.getLocalNumCols() != 0,
        std::runtime_error, "Tpetra::writeMatrixStrip: pid "
        << myRank << " should have 0 columns but has " << A.getLocalNumCols());

    } else {

      TEUCHOS_TEST_FOR_EXCEPTION(numRows != A.getLocalNumRows(),
        std::runtime_error, "Tpetra::writeMatrixStrip: "
        "number of rows on pid 0 does not match global number of rows");


      int err = 0;
      const LO blockSize = A.getBlockSize();
      const size_t numLocalRows = A.getLocalNumRows();
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
  Teuchos::RCP<Tpetra::CrsGraph<LO, GO, Node> >
  getBlockCrsGraph(const Tpetra::CrsMatrix<Scalar, LO, GO, Node>& pointMatrix, const LO &blockSize, bool use_LID)
  {
      TEUCHOS_FUNC_TIME_MONITOR_DIFF("Tpetra::getBlockCrsGraph", getBlockCrsGraph0);
      /*
        ASSUMPTIONS:

           1) In point matrix, all entries associated with a little block are present (even if they are zero).
           2) For given mesh DOF, point DOFs appear consecutively and in ascending order in row & column maps.
           3) Point column map and block column map are ordered consistently.
      */

      using Teuchos::Array;
      using Teuchos::ArrayView;
      using Teuchos::RCP;

      typedef Tpetra::Map<LO,GO,Node>                   map_type;
      typedef Tpetra::CrsGraph<LO,GO,Node>              crs_graph_type;
      typedef Tpetra::CrsMatrix<Scalar, LO,GO,Node>     crs_matrix_type;

      using local_graph_device_type  = typename crs_matrix_type::local_graph_device_type;
      using row_map_type             = typename local_graph_device_type::row_map_type::non_const_type;
      using entries_type             = typename local_graph_device_type::entries_type::non_const_type;

      using offset_type              = typename row_map_type::non_const_value_type;

      using execution_space = typename Node::execution_space;
      using range_type = Kokkos::RangePolicy<execution_space, LO>;

      const map_type &pointRowMap = *(pointMatrix.getRowMap());
      RCP<const map_type> meshRowMap, meshColMap, meshDomainMap, meshRangeMap;

      const map_type &pointColMap = *(pointMatrix.getColMap());
      const map_type &pointDomainMap = *(pointMatrix.getDomainMap());
      const map_type &pointRangeMap = *(pointMatrix.getRangeMap());

      {
        TEUCHOS_FUNC_TIME_MONITOR_DIFF("Tpetra::getBlockCrsGraph::createMeshMaps", getBlockCrsGraph1);
        meshRowMap = createMeshMap<LO,GO,Node>(blockSize, pointRowMap, use_LID);
        meshColMap = createMeshMap<LO,GO,Node>(blockSize, pointColMap, use_LID);
        meshDomainMap = createMeshMap<LO,GO,Node>(blockSize, pointDomainMap, use_LID);
        meshRangeMap = createMeshMap<LO,GO,Node>(blockSize, pointRangeMap, use_LID);
        Kokkos::DefaultExecutionSpace().fence();
      }

      if(meshColMap.is_null()) throw std::runtime_error("ERROR: Cannot create mesh colmap");

      auto localMeshColMap = meshColMap->getLocalMap();
      auto localPointColMap = pointColMap.getLocalMap();
      auto localPointRowMap = pointRowMap.getLocalMap();

      RCP<crs_graph_type> meshCrsGraph;

      const offset_type bs2 = blockSize * blockSize;

      if (use_LID) {
        TEUCHOS_FUNC_TIME_MONITOR_DIFF("Tpetra::getBlockCrsGraph::LID", getBlockCrsGraph2);
        auto pointLocalGraph = pointMatrix.getCrsGraph()->getLocalGraphDevice();
        auto pointRowptr = pointLocalGraph.row_map;
        auto pointColind = pointLocalGraph.entries;

        TEUCHOS_TEST_FOR_EXCEPTION(pointColind.extent(0) % bs2 != 0,
          std::runtime_error, "Tpetra::getBlockCrsGraph: "
          "local number of non zero entries is not a multiple of blockSize^2 ");

        LO block_rows = pointRowptr.extent(0) == 0 ? 0 : (pointRowptr.extent(0)-1)/blockSize;
        row_map_type blockRowptr("blockRowptr", block_rows+1);
        entries_type blockColind("blockColind", pointColind.extent(0)/(bs2));

        Kokkos::parallel_for("fillMesh",range_type(0,block_rows), KOKKOS_LAMBDA(const LO i) {

          const LO offset_b = pointRowptr(i*blockSize)/bs2;
          const LO offset_b_max = pointRowptr((i+1)*blockSize)/bs2;

          if (i==block_rows-1)
            blockRowptr(i+1) = offset_b_max;
          blockRowptr(i) = offset_b;

          const LO offset_p = pointRowptr(i*blockSize);

          for (LO k=0; k<offset_b_max-offset_b; ++k) {
            blockColind(offset_b + k) = pointColind(offset_p + k * blockSize)/blockSize;
          }
        });

        meshCrsGraph = rcp(new crs_graph_type(meshRowMap, meshColMap, blockRowptr, blockColind));
        meshCrsGraph->fillComplete(meshDomainMap,meshRangeMap);
        Kokkos::DefaultExecutionSpace().fence();
      }
      else {
        TEUCHOS_FUNC_TIME_MONITOR_DIFF("Tpetra::getBlockCrsGraph::GID", getBlockCrsGraph3);
        auto pointLocalGraph = pointMatrix.getCrsGraph()->getLocalGraphDevice();
        auto pointRowptr = pointLocalGraph.row_map;
        auto pointColind = pointLocalGraph.entries;

        TEUCHOS_TEST_FOR_EXCEPTION(pointColind.extent(0) % bs2 != 0,
          std::runtime_error, "Tpetra::getBlockCrsGraph: "
          "local number of non zero entries is not a multiple of blockSize^2 ");

        LO block_rows = pointRowptr.extent(0) == 0 ? 0 : (pointRowptr.extent(0)-1)/blockSize;
        row_map_type blockRowptr("blockRowptr", block_rows+1);
        entries_type blockColind("blockColind", pointColind.extent(0)/(bs2));

        Kokkos::parallel_for("fillMesh",range_type(0,block_rows), KOKKOS_LAMBDA(const LO i) {

          const LO offset_b = pointRowptr(i*blockSize)/bs2;
          const LO offset_b_max = pointRowptr((i+1)*blockSize)/bs2;

          if (i==block_rows-1)
            blockRowptr(i+1) = offset_b_max;
          blockRowptr(i) = offset_b;

          const LO offset_p = pointRowptr(i*blockSize);
          const LO offset_p_max = pointRowptr((i+1)*blockSize);

          LO filled_block = 0;
          for (LO p_i=0; p_i<offset_p_max-offset_p; ++p_i) {
            auto bcol_GID = localPointColMap.getGlobalElement(pointColind(offset_p + p_i))/blockSize;
            auto bcol_LID = localMeshColMap.getLocalElement(bcol_GID);

            bool visited = false;
            for (LO k=0; k<filled_block; ++k) {
              if (blockColind(offset_b + k) == bcol_LID)
                visited = true;
            }
            if (!visited) {
              blockColind(offset_b + filled_block) = bcol_LID;
              ++filled_block;
            }
          }
        });

        meshCrsGraph = rcp(new crs_graph_type(meshRowMap, meshColMap, blockRowptr, blockColind));
        meshCrsGraph->fillComplete(meshDomainMap,meshRangeMap);
        Kokkos::DefaultExecutionSpace().fence();
      }

      return meshCrsGraph;
  }

  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<BlockCrsMatrix<Scalar, LO, GO, Node> >
  convertToBlockCrsMatrix(const Tpetra::CrsMatrix<Scalar, LO, GO, Node>& pointMatrix, const LO &blockSize, bool use_LID)
  {
      TEUCHOS_FUNC_TIME_MONITOR_DIFF("Tpetra::convertToBlockCrsMatrix", convertToBlockCrsMatrix0);
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
      typedef Tpetra::CrsMatrix<Scalar, LO,GO,Node>     crs_matrix_type;

      using local_graph_device_type  = typename crs_matrix_type::local_graph_device_type;
      using local_matrix_device_type = typename crs_matrix_type::local_matrix_device_type;
      using row_map_type             = typename local_graph_device_type::row_map_type::non_const_type;
      using values_type              = typename local_matrix_device_type::values_type::non_const_type;

      using offset_type              = typename row_map_type::non_const_value_type;

      using execution_space = typename Node::execution_space;
      using range_type = Kokkos::RangePolicy<execution_space, LO>;

      RCP<block_crs_matrix_type> blockMatrix;

      const offset_type bs2 = blockSize * blockSize;

      auto meshCrsGraph = getBlockCrsGraph(pointMatrix, blockSize, use_LID);

      if (use_LID) {
        TEUCHOS_FUNC_TIME_MONITOR_DIFF("Tpetra::convertToBlockCrsMatrix::LID", convertToBlockCrsMatrix1);
        auto pointLocalGraph = pointMatrix.getCrsGraph()->getLocalGraphDevice();
        auto pointRowptr = pointLocalGraph.row_map;
        auto pointColind = pointLocalGraph.entries;

        offset_type block_rows = pointRowptr.extent(0) == 0 ? 0 : (pointRowptr.extent(0)-1)/blockSize;
        values_type blockValues("values",  meshCrsGraph->getLocalNumEntries()*bs2);
        auto pointValues = pointMatrix.getLocalValuesDevice (Access::ReadOnly);
        auto blockRowptr = meshCrsGraph->getLocalGraphDevice().row_map;

        Kokkos::parallel_for("copyblockValues",range_type(0,block_rows),KOKKOS_LAMBDA(const LO i) {
          const offset_type blkBeg    = blockRowptr[i];
          const offset_type numBlocks = blockRowptr[i+1] - blkBeg;

          // For each block in the row...
          for (offset_type block=0; block < numBlocks; block++) {

            // For each entry in the block...
            for(LO little_row=0; little_row<blockSize; little_row++) {
              offset_type point_row_offset = pointRowptr[i*blockSize + little_row];
              for(LO little_col=0; little_col<blockSize; little_col++) {
                blockValues((blkBeg+block) * bs2 + little_row * blockSize + little_col) = 
                  pointValues[point_row_offset + block*blockSize + little_col];
              }
            }

          }
          });
        blockMatrix = rcp(new block_crs_matrix_type(*meshCrsGraph, blockValues, blockSize));
        Kokkos::DefaultExecutionSpace().fence();
      }
      else {
        TEUCHOS_FUNC_TIME_MONITOR_DIFF("Tpetra::convertToBlockCrsMatrix::GID", convertToBlockCrsMatrix2);
        auto localMeshColMap = meshCrsGraph->getColMap()->getLocalMap();
        auto localPointColMap = pointMatrix.getColMap()->getLocalMap();

        values_type blockValues("values",  meshCrsGraph->getLocalNumEntries()*bs2);
        {
          auto pointLocalGraph = pointMatrix.getCrsGraph()->getLocalGraphDevice();
          auto pointRowptr = pointLocalGraph.row_map;
          auto pointColind = pointLocalGraph.entries;

          offset_type block_rows = pointRowptr.extent(0) == 0 ? 0 : (pointRowptr.extent(0)-1)/blockSize;
          auto pointValues = pointMatrix.getLocalValuesDevice (Access::ReadOnly);
          auto blockRowptr = meshCrsGraph->getLocalGraphDevice().row_map;
          auto blockColind = meshCrsGraph->getLocalGraphDevice().entries;

          row_map_type pointGColind("pointGColind", pointColind.extent(0));

          Kokkos::parallel_for("computePointGColind",range_type(0,pointColind.extent(0)),KOKKOS_LAMBDA(const LO i) {
            pointGColind(i) = localPointColMap.getGlobalElement(pointColind(i));
          });

          row_map_type blockGColind("blockGColind", blockColind.extent(0));

          Kokkos::parallel_for("computeBlockGColind",range_type(0,blockGColind.extent(0)),KOKKOS_LAMBDA(const LO i) {
            blockGColind(i) = localMeshColMap.getGlobalElement(blockColind(i));
          });

          Kokkos::parallel_for("copyblockValues",range_type(0,block_rows),KOKKOS_LAMBDA(const LO i) {
            const offset_type blkBeg    = blockRowptr[i];
            const offset_type numBlocks = blockRowptr[i+1] - blkBeg;

            for (offset_type point_i=0; point_i < pointRowptr[i*blockSize + 1] - pointRowptr[i*blockSize]; point_i++) {

              offset_type block_inv=static_cast<offset_type>(-1);
              offset_type little_col_inv=static_cast<offset_type>(-1);
              for (offset_type block_2=0; block_2 < numBlocks; block_2++) {
                for (LO little_col_2=0; little_col_2 < blockSize; little_col_2++) {
                  if (blockGColind(blkBeg+block_2)*blockSize + little_col_2 == pointGColind(pointRowptr[i*blockSize] + point_i)) {
                    block_inv = block_2;
                    little_col_inv = little_col_2;
                    break;
                  }
                }
                if (block_inv!=static_cast<offset_type>(-1))
                  break;
              }

              for(LO little_row=0; little_row<blockSize; little_row++) {
                blockValues((blkBeg+block_inv) * bs2 + little_row * blockSize + little_col_inv) = pointValues[pointRowptr[i*blockSize+little_row] + point_i];
              }
            }
            });
        }
        blockMatrix = rcp(new block_crs_matrix_type(*meshCrsGraph, blockValues, blockSize));
        Kokkos::DefaultExecutionSpace().fence();
      }

      return blockMatrix;

  }

  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LO, GO, Node>>
  fillLogicalBlocks(const Tpetra::CrsMatrix<Scalar, LO, GO, Node>& pointMatrix, const LO &blockSize)
  {
    using crs_t           = Tpetra::CrsMatrix<Scalar, LO, GO, Node>;
    using execution_space = typename Node::execution_space;
    using dev_row_view_t  = typename crs_t::local_graph_device_type::row_map_type::non_const_type;
    using dev_col_view_t  = typename crs_t::local_graph_device_type::entries_type::non_const_type;
    using dev_val_view_t  = typename crs_t::local_matrix_device_type::values_type::non_const_type;
    using range_type      = Kokkos::RangePolicy<execution_space, size_t>;
    using team_policy     = Kokkos::TeamPolicy<execution_space>;
    using member_type     = typename team_policy::member_type;
    using scratch_view    = Kokkos::View<bool*, typename execution_space::scratch_memory_space>;
    using Ordinal         = typename dev_row_view_t::non_const_value_type;

    // Get structure / values
    auto local = pointMatrix.getLocalMatrixDevice();
    auto row_ptrs = local.graph.row_map;
    auto col_inds = local.graph.entries;
    auto values = local.values;
    const auto nrows = pointMatrix.getLocalNumRows();

    //
    // Populate all active blocks, they must be fully populated with entries so fill with zeroes
    //

    // Make row workspace views
    dev_row_view_t new_rowmap("new_rowmap", nrows+1);
    const auto blocks_per_row = nrows / blockSize; // assumes square matrix
    dev_row_view_t active_block_row_map("active_block_row_map", blocks_per_row + 1);
    const int max_threads = execution_space::concurrency();
    assert(blockSize > 1);
    assert(nrows % blockSize == 0);
    const int mem_level = 1;
    const int bytes = scratch_view::shmem_size(blocks_per_row);

    if (max_threads >= blockSize) {
      // Prefer 1 team per block since this will require a lot less scratch memory
      team_policy tp(blocks_per_row, blockSize);

      // Count active blocks
      Kokkos::parallel_for("countActiveBlocks", tp.set_scratch_size(mem_level, Kokkos::PerTeam(bytes)), KOKKOS_LAMBDA(const member_type& team) {
        Ordinal block_row = team.league_rank();

        scratch_view row_block_active(team.team_scratch(mem_level), blocks_per_row);
        Kokkos::single(
          Kokkos::PerTeam(team), [&] () {
            for (size_t row_block_idx = 0; row_block_idx < blocks_per_row; ++row_block_idx) {
              row_block_active(row_block_idx) = false;
            }
        });
        team.team_barrier();

        // All threads in a team scan a blocks-worth of rows to see which
        // blocks are active
        Kokkos::parallel_for(
          Kokkos::TeamThreadRange(team, blockSize), [&] (Ordinal block_offset) {

            Ordinal row     = block_row*blockSize + block_offset;
            Ordinal row_itr = row_ptrs(row);
            Ordinal row_end = row_ptrs(row+1);

            for (size_t row_block_idx = 0; row_block_idx < blocks_per_row; ++row_block_idx) {
              const Ordinal first_possible_col_in_block      = row_block_idx * blockSize;
              const Ordinal first_possible_col_in_next_block = (row_block_idx+1) * blockSize;
              Ordinal curr_nnz_col = col_inds(row_itr);
              while (curr_nnz_col >= first_possible_col_in_block && curr_nnz_col < first_possible_col_in_next_block && row_itr < row_end) {
                // This block has at least one nnz in this row
                row_block_active(row_block_idx) = true;
                ++row_itr;
                if (row_itr == row_end) break;
                curr_nnz_col = col_inds(row_itr);
              }
            }
        });

        team.team_barrier();

        Kokkos::single(
          Kokkos::PerTeam(team), [&] () {
            Ordinal count = 0;
            for (size_t row_block_idx = 0; row_block_idx < blocks_per_row; ++row_block_idx) {
              if (row_block_active(row_block_idx)) {
                ++count;
              }
            }
            active_block_row_map(block_row) = count;
        });
      });
    }
    else {
      // We don't have enough parallelism to make a thread team handle a block, so just
      // have 1 thread per block
      team_policy tp(blocks_per_row, 1);

      // Count active blocks
      Kokkos::parallel_for("countActiveBlocks", tp.set_scratch_size(mem_level, Kokkos::PerTeam(bytes)), KOKKOS_LAMBDA(const member_type& team) {
        Ordinal block_row = team.league_rank();

        scratch_view row_block_active(team.team_scratch(mem_level), blocks_per_row);
        for (size_t row_block_idx = 0; row_block_idx < blocks_per_row; ++row_block_idx) {
          row_block_active(row_block_idx) = false;
        }

        // One thread scans a blocks-worth of rows to see which blocks are active
        for (int block_offset=0; block_offset < blockSize; ++block_offset) {
          Ordinal row     = block_row*blockSize + block_offset;
          Ordinal row_itr = row_ptrs(row);
          Ordinal row_end = row_ptrs(row+1);

          for (size_t row_block_idx = 0; row_block_idx < blocks_per_row; ++row_block_idx) {
            const Ordinal first_possible_col_in_block      = row_block_idx * blockSize;
            const Ordinal first_possible_col_in_next_block = (row_block_idx+1) * blockSize;
            Ordinal curr_nnz_col = col_inds(row_itr);
            while (curr_nnz_col >= first_possible_col_in_block && curr_nnz_col < first_possible_col_in_next_block && row_itr < row_end) {
              // This block has at least one nnz in this row
              row_block_active(row_block_idx) = true;
              ++row_itr;
              if (row_itr == row_end) break;
              curr_nnz_col = col_inds(row_itr);
            }
          }
        }

        Ordinal count = 0;
        for (size_t row_block_idx = 0; row_block_idx < blocks_per_row; ++row_block_idx) {
          if (row_block_active(row_block_idx)) {
            ++count;
          }
        }
        active_block_row_map(block_row) = count;
      });
    }

    Ordinal nnz_block_count = 0;
#if KOKKOSKERNELS_VERSION >= 40199
    KokkosKernels::Impl::kk_exclusive_parallel_prefix_sum<
      execution_space>(active_block_row_map.extent(0), active_block_row_map, nnz_block_count);
#else
    KokkosKernels::Impl::kk_exclusive_parallel_prefix_sum<
      dev_row_view_t, execution_space>(active_block_row_map.extent(0), active_block_row_map, nnz_block_count);
#endif
    dev_col_view_t block_col_ids("block_col_ids", nnz_block_count);

    // Find active blocks
    if (max_threads >= blockSize) {
      // Prefer 1 team per block since this will require a lot less scratch memory
      team_policy tp(blocks_per_row, blockSize);

      Kokkos::parallel_for("findActiveBlocks", tp.set_scratch_size(mem_level, Kokkos::PerTeam(bytes)), KOKKOS_LAMBDA(const member_type& team) {
        Ordinal block_row = team.league_rank();

        scratch_view row_block_active(team.team_scratch(mem_level), blocks_per_row);
        Kokkos::single(
          Kokkos::PerTeam(team), [&] () {
            for (size_t row_block_idx = 0; row_block_idx < blocks_per_row; ++row_block_idx) {
              row_block_active(row_block_idx) = false;
            }
        });
        team.team_barrier();

        // All threads in a team scan a blocks-worth of rows to see which
        // blocks are active
        Kokkos::parallel_for(
          Kokkos::TeamThreadRange(team, blockSize), [&] (Ordinal block_offset) {

            Ordinal row     = block_row*blockSize + block_offset;
            Ordinal row_itr = row_ptrs(row);
            Ordinal row_end = row_ptrs(row+1);

            for (size_t row_block_idx = 0; row_block_idx < blocks_per_row; ++row_block_idx) {
              const Ordinal first_possible_col_in_block      = row_block_idx * blockSize;
              const Ordinal first_possible_col_in_next_block = (row_block_idx+1) * blockSize;
              Ordinal curr_nnz_col = col_inds(row_itr);
              while (curr_nnz_col >= first_possible_col_in_block && curr_nnz_col < first_possible_col_in_next_block && row_itr < row_end) {
                // This block has at least one nnz in this row
                row_block_active(row_block_idx) = true;
                ++row_itr;
                if (row_itr == row_end) break;
                curr_nnz_col = col_inds(row_itr);
              }
            }
        });

        team.team_barrier();

        Kokkos::single(
          Kokkos::PerTeam(team), [&] () {
            Ordinal offset = active_block_row_map[block_row];
            for (size_t row_block_idx = 0; row_block_idx < blocks_per_row; ++row_block_idx) {
              if (row_block_active(row_block_idx)) {
                block_col_ids(offset) = row_block_idx;
                ++offset;
              }
            }
        });
      });
    }
    else {
      team_policy tp(blocks_per_row, 1);

      Kokkos::parallel_for("findActiveBlocks", tp.set_scratch_size(mem_level, Kokkos::PerTeam(bytes)), KOKKOS_LAMBDA(const member_type& team) {
        Ordinal block_row = team.league_rank();

        scratch_view row_block_active(team.team_scratch(mem_level), blocks_per_row);
        for (size_t row_block_idx = 0; row_block_idx < blocks_per_row; ++row_block_idx) {
          row_block_active(row_block_idx) = false;
        }

        // One thread scans a blocks-worth of rows to see which blocks are active
        for (int block_offset=0; block_offset < blockSize; ++block_offset) {
          Ordinal row     = block_row*blockSize + block_offset;
          Ordinal row_itr = row_ptrs(row);
          Ordinal row_end = row_ptrs(row+1);

          for (size_t row_block_idx = 0; row_block_idx < blocks_per_row; ++row_block_idx) {
            const Ordinal first_possible_col_in_block      = row_block_idx * blockSize;
            const Ordinal first_possible_col_in_next_block = (row_block_idx+1) * blockSize;
            Ordinal curr_nnz_col = col_inds(row_itr);
            while (curr_nnz_col >= first_possible_col_in_block && curr_nnz_col < first_possible_col_in_next_block && row_itr < row_end) {
              // This block has at least one nnz in this row
              row_block_active(row_block_idx) = true;
              ++row_itr;
              if (row_itr == row_end) break;
              curr_nnz_col = col_inds(row_itr);
            }
          }
        }

        Ordinal offset = active_block_row_map[block_row];
        for (size_t row_block_idx = 0; row_block_idx < blocks_per_row; ++row_block_idx) {
          if (row_block_active(row_block_idx)) {
            block_col_ids(offset) = row_block_idx;
            ++offset;
          }
        }
      });
    }

    // Sizing
    Kokkos::parallel_for("sizing", range_type(0, nrows), KOKKOS_LAMBDA(const size_t row) {
      const auto block_row = row / blockSize;
      const Ordinal block_row_begin = active_block_row_map(block_row);
      const Ordinal block_row_end   = active_block_row_map(block_row+1);
      const Ordinal row_nnz         = (block_row_end - block_row_begin) * blockSize;
      new_rowmap(row) = row_nnz;
    });

    Ordinal new_nnz_count = 0;
#if KOKKOSKERNELS_VERSION >= 40199
    KokkosKernels::Impl::kk_exclusive_parallel_prefix_sum<
      execution_space>(new_rowmap.extent(0), new_rowmap, new_nnz_count);
#else
    KokkosKernels::Impl::kk_exclusive_parallel_prefix_sum<
      dev_row_view_t, execution_space>(new_rowmap.extent(0), new_rowmap, new_nnz_count);
#endif
    // Now populate cols and vals
    dev_col_view_t new_col_ids("new_col_ids", new_nnz_count);
    dev_val_view_t new_vals("new_vals",       new_nnz_count);
    Kokkos::parallel_for("entries", range_type(0, nrows), KOKKOS_LAMBDA(const size_t row) {
      Ordinal row_itr = row_ptrs(row);
      Ordinal row_end = row_ptrs(row+1);
      Ordinal row_itr_new = new_rowmap(row);

      Ordinal block_row       = row / blockSize;
      Ordinal block_row_begin = active_block_row_map(block_row);
      Ordinal block_row_end   = active_block_row_map(block_row+1);

      for (Ordinal row_block_idx = block_row_begin; row_block_idx < block_row_end; ++row_block_idx) {
        const Ordinal block_col                        = block_col_ids(row_block_idx);
        const Ordinal first_possible_col_in_block      = block_col * blockSize;
        const Ordinal first_possible_col_in_next_block = (block_col+1) * blockSize;
        for (Ordinal possible_col = first_possible_col_in_block; possible_col < first_possible_col_in_next_block; ++possible_col, ++row_itr_new) {
          new_col_ids(row_itr_new) = possible_col;
          Ordinal curr_nnz_col = col_inds(row_itr);
          if (possible_col == curr_nnz_col && row_itr < row_end) {
            // Already a non-zero entry
            new_vals(row_itr_new) = values(row_itr);
            ++row_itr;
          }
        }
      }
    });

    // Create new, filled CRS
    auto crs_row_map = pointMatrix.getRowMap();
    auto crs_col_map = pointMatrix.getColMap();
    Teuchos::RCP<crs_t> result = Teuchos::rcp(new crs_t(crs_row_map, crs_col_map, new_rowmap, new_col_ids, new_vals));
    result->fillComplete();
    return result;
  }

  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LO, GO, Node>>
  unfillFormerBlockCrs(const Tpetra::CrsMatrix<Scalar, LO, GO, Node>& pointMatrix)
  {
    using crs_t           = Tpetra::CrsMatrix<Scalar, LO, GO, Node>;
    using dev_row_view_t  = typename crs_t::local_graph_device_type::row_map_type::non_const_type;
    using dev_col_view_t  = typename crs_t::local_graph_device_type::entries_type::non_const_type;
    using dev_val_view_t  = typename crs_t::local_matrix_device_type::values_type::non_const_type;
    using impl_scalar_t   = typename Kokkos::ArithTraits<Scalar>::val_type;
    using STS             = Kokkos::ArithTraits<impl_scalar_t>;
    using Ordinal         = typename dev_row_view_t::non_const_value_type;
    using execution_space = typename Node::execution_space;
    using range_type      = Kokkos::RangePolicy<execution_space, size_t>;

    // Get structure / values
    auto local = pointMatrix.getLocalMatrixHost();
    auto row_ptrs = local.graph.row_map;
    auto col_inds = local.graph.entries;
    auto values = local.values;
    const auto nrows = pointMatrix.getLocalNumRows();

    // Sizing and rows
    dev_row_view_t new_rowmap("new_rowmap", nrows+1);
    const impl_scalar_t zero = STS::zero();
    Kokkos::parallel_for("sizing", range_type(0, nrows), KOKKOS_LAMBDA(const size_t row) {
      const Ordinal row_nnz_begin = row_ptrs(row);
      Ordinal row_nnzs = 0;
      for (Ordinal row_nnz = row_nnz_begin; row_nnz < row_ptrs(row+1); ++row_nnz) {
        const impl_scalar_t value = values(row_nnz);
        if (value != zero) {
          ++row_nnzs;
        }
      }
      new_rowmap(row) = row_nnzs;
    });

    Ordinal real_nnzs = 0;
#if KOKKOSKERNELS_VERSION >= 40199
    KokkosKernels::Impl::kk_exclusive_parallel_prefix_sum<
      execution_space>(new_rowmap.extent(0), new_rowmap, real_nnzs);
#else
    KokkosKernels::Impl::kk_exclusive_parallel_prefix_sum<
      dev_row_view_t, execution_space>(new_rowmap.extent(0), new_rowmap, real_nnzs);
#endif
    // Now populate cols and vals
    dev_col_view_t new_col_ids("new_col_ids", real_nnzs);
    dev_val_view_t new_vals("new_vals",       real_nnzs);
    Kokkos::parallel_for("entries", range_type(0, nrows), KOKKOS_LAMBDA(const size_t row) {
      Ordinal row_nnzs = 0;
      const Ordinal old_row_nnz_begin = row_ptrs(row);
      const Ordinal new_row_nnz_begin = new_rowmap(row);
      for (Ordinal old_row_nnz = old_row_nnz_begin; old_row_nnz < row_ptrs(row+1); ++old_row_nnz) {
        const impl_scalar_t value = values(old_row_nnz);
        if (value != zero) {
          new_col_ids(new_row_nnz_begin + row_nnzs) = col_inds(old_row_nnz);
          new_vals(new_row_nnz_begin + row_nnzs) = value;
          ++row_nnzs;
        }
      }
    });

    // Create new, unfilled CRS
    auto crs_row_map = pointMatrix.getRowMap();
    auto crs_col_map = pointMatrix.getColMap();
    Teuchos::RCP<crs_t> result = Teuchos::rcp(new crs_t(crs_row_map, crs_col_map, new_rowmap, new_col_ids, new_vals));
    result->fillComplete();
    return result;
  }

  template<class LO, class GO, class Node>
  Teuchos::RCP<const Tpetra::Map<LO,GO,Node> >
  createMeshMap (const LO& blockSize, const Tpetra::Map<LO,GO,Node>& pointMap, bool use_LID)
  {
    typedef Teuchos::OrdinalTraits<Tpetra::global_size_t> TOT;
    typedef Tpetra::Map<LO,GO,Node> map_type;

    if(use_LID) {

      using execution_space = typename Node::execution_space;
      using range_type      = Kokkos::RangePolicy<execution_space, size_t>;

      auto pointGlobalID = pointMap.getMyGlobalIndicesDevice();
      LO block_rows = pointGlobalID.extent(0)/blockSize;
      Kokkos::View<GO*, typename map_type::device_type> meshGlobalID("meshGlobalID", block_rows);

      Kokkos::parallel_for("fillMeshMap",range_type(0,block_rows), KOKKOS_LAMBDA(const LO i) {
        meshGlobalID(i) = pointGlobalID(i*blockSize)/blockSize;
      });

      Teuchos::RCP<const map_type> meshMap = Teuchos::rcp( new map_type(TOT::invalid(), meshGlobalID, 0, pointMap.getComm()) );
      return meshMap;
    }
    else {
      //calculate mesh GIDs
      Teuchos::ArrayView<const GO> pointGids = pointMap.getLocalElementList();
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
  }


  template<class LO, class GO, class Node>
  Teuchos::RCP<const Tpetra::Map<LO,GO,Node> >
  createPointMap (const LO& blockSize, const Tpetra::Map<LO,GO,Node>& blockMap)
  {
    typedef Teuchos::OrdinalTraits<Tpetra::global_size_t> TOT;
    typedef Tpetra::Map<LO,GO,Node> map_type;

    //calculate mesh GIDs
    Teuchos::ArrayView<const GO> blockGids = blockMap.getLocalElementList();
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
    size_t block_nnz = blockMatrix.getLocalNumEntries();
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
    LO point_rows = (LO) pointRowMap->getLocalNumElements();
    LO block_rows = (LO) blockRowMap->getLocalNumElements();
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
  template Teuchos::RCP<BlockCrsMatrix<S, LO, GO, NODE> > convertToBlockCrsMatrix(const CrsMatrix<S, LO, GO, NODE>& pointMatrix, const LO &blockSize, bool use_LID); \
  template Teuchos::RCP<CrsMatrix<S, LO, GO, NODE> > fillLogicalBlocks(const CrsMatrix<S, LO, GO, NODE>& pointMatrix, const LO &blockSize); \
  template Teuchos::RCP<CrsMatrix<S, LO, GO, NODE> > unfillFormerBlockCrs(const CrsMatrix<S, LO, GO, NODE>& pointMatrix); \
  template Teuchos::RCP<CrsMatrix<S, LO, GO, NODE> > convertToCrsMatrix(const BlockCrsMatrix<S, LO, GO, NODE>& blockMatrix); \
  template Teuchos::RCP<CrsGraph<LO, GO, NODE>> getBlockCrsGraph(const Tpetra::CrsMatrix<S, LO, GO, NODE>& pointMatrix, const LO &blockSize, bool use_local_ID);

//
// Explicit instantiation macro for createMeshMap / createPointMap
//
#define TPETRA_CREATEMESHMAP_INSTANT(LO,GO,NODE) \
  template Teuchos::RCP<const Map<LO,GO,NODE> > createMeshMap  (const LO& blockSize, const Map<LO,GO,NODE>& pointMap, bool use_local_ID); \
  template Teuchos::RCP<const Map<LO,GO,NODE> > createPointMap (const LO& blockSize, const Map<LO,GO,NODE>& blockMap);

#endif // TPETRA_BLOCKCRSMATRIX_HELPERS_DEF_HPP
