// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <vector>
#include <numeric>
#include <iostream>
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>

// Create and return a simple example CrsMatrix, with row distribution
// over the given Map.
int main(int argc, char *argv[]) {

  using Teuchos::rcp;
  using Teuchos::RCP;
  using map_type = Tpetra::Map<>;
  using LO = typename map_type::local_ordinal_type;
  using GO = typename map_type::global_ordinal_type;
  using crs_matrix_type = Tpetra::CrsMatrix<>;

  Teuchos::oblackholestream blackHole;
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);
  {
    Teuchos::RCP<const Teuchos::Comm<int>> comm = Tpetra::getDefaultComm();
    const int myRank = comm->getRank();

    TEUCHOS_TEST_FOR_EXCEPTION(comm->getSize() != 2, std::logic_error,
                               "Example requires 2 process ranks");

    const Tpetra::global_size_t numGblRows = 9;
    const Tpetra::global_size_t numGblCols = 14;
    const GO indexBase = 0;

    // Create the row and column index lists on each processor
    const LO numRow = (myRank == 0 ? 5 : 4);
    const LO numCol = (myRank == 0 ? 8 : 6);
    std::vector<GO> rowIndices(numRow);
    std::vector<GO> colIndices(numCol);
    {
      GO rowStart = (myRank == 0) ? 0 : 5;
      std::iota(rowIndices.begin(), rowIndices.end(), rowStart);
      GO colStart = (myRank == 0) ? 0 : 2;
      std::iota(colIndices.begin(), colIndices.end(), colStart);
    }

    // Create the row map
    RCP<const map_type> rowMap = rcp(
        new map_type(numGblRows, rowIndices.data(), numRow, indexBase, comm));

    // Create the column map
    RCP<const map_type> colMap = rcp(
        new map_type(numGblCols, colIndices.data(), numCol, indexBase, comm));

    // Create a Tpetra sparse matrix whose rows have distribution
    // given by the row Map and column Map.
    RCP<crs_matrix_type> A(new crs_matrix_type(rowMap, colMap, 0));

    // This tells the Trilinos test framework that the test passed.
    if (myRank == 0) {
      std::cout << "End Result: TEST PASSED" << std::endl;
    }
    // ScopeGuard's destructor calls MPI_Finalize, if its constructor
    // called MPI_Init.  Likewise, it calls Kokkos::finalize, if its
    // constructor called Kokkos::initialize.
  }
  return 0;
}
