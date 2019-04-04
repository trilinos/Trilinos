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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
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
