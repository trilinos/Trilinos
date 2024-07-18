// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <iostream>

// Create and return a simple example CrsMatrix, with row distribution
// over the given Map.
//
// CrsMatrixType: The type of the Tpetra::CrsMatrix specialization to use.
template<class CrsMatrixType>
Teuchos::RCP<const CrsMatrixType>
createMatrix(const Teuchos::RCP<const typename CrsMatrixType::map_type>& map)
{
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::tuple;
  typedef Tpetra::global_size_t GST;

  // Fetch typedefs from the Tpetra::CrsMatrix.
  using scalar_type = typename CrsMatrixType::scalar_type;
  using GO = typename CrsMatrixType::global_ordinal_type;

  // Create a Tpetra::Matrix using the Map, with dynamic allocation.
  RCP<CrsMatrixType> A(new CrsMatrixType(map, 3));

  // Add rows one at a time.  Off diagonal values will always be -1.
  const scalar_type two    = static_cast<scalar_type>( 2.0);
  const scalar_type negOne = static_cast<scalar_type>(-1.0);
  const GST numGlobalIndices = map->getGlobalNumElements();

  // The list of global elements owned by this MPI process.
  auto myGlobalIndices = map->getMyGlobalIndices();
  for(GO i=0; i<static_cast<GO>(myGlobalIndices.size()); i++) {
    const GO i_global = myGlobalIndices(i);

    // Can't insert local indices without a column map, so we insert
    // global indices here.
    if(i_global == 0) {
      A->insertGlobalValues(i_global,
                            tuple(i_global, i_global+1),
                            tuple(two, negOne));
    } else if(static_cast<GST>(i_global) == numGlobalIndices - 1) {
      A->insertGlobalValues(i_global,
                            tuple(i_global-1, i_global),
                            tuple(negOne, two));
    } else {
      A->insertGlobalValues(i_global,
                            tuple(i_global-1, i_global, i_global+1),
                            tuple(negOne, two, negOne));
    }
  }
  // Finish up the matrix.
  A->fillComplete();
  return A;
}

int
main(int argc, char *argv[])
{
  using Teuchos::rcp;
  using Teuchos::RCP;
  using GST = Tpetra::global_size_t;
  using scalar_type = Tpetra::CrsMatrix<>::scalar_type;
  using crs_matrix_type = Tpetra::CrsMatrix<scalar_type>;
  using map_type = Tpetra::Map<>;
  using global_ordinal_type = Tpetra::Map<>::global_ordinal_type;

  Teuchos::oblackholestream blackHole;
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);
  {
    // Get a communicator corresponding to MPI_COMM_WORLD
    Teuchos::RCP<const Teuchos::Comm<int>> comm = Tpetra::getDefaultComm();
    const int myRank = comm->getRank();

    // The global number of rows in the matrix A to create.  We scale
    // this relative to the number of (MPI) processes, so that no matter
    // how many MPI processes you run, every process will have 10 rows.
    const GST numGlobalIndices = 10 * comm->getSize();
    const global_ordinal_type indexBase = 0;

    // Construct a Map that puts approximately the same number of
    // equations on each processor.
    RCP<const map_type> globalMap =
      rcp(new map_type(numGlobalIndices, indexBase, comm, Tpetra::GloballyDistributed));

    // Create a sparse matrix using globalMap
    RCP<const crs_matrix_type> A = createMatrix<crs_matrix_type>(globalMap);

    // This tells the Trilinos test framework that the test passed.
    if(myRank == 0) {
      std::cout << "End Result: TEST PASSED" << std::endl;
    }
    // ScopeGuard's destructor calls MPI_Finalize, if its constructor
    // called MPI_Init.  Likewise, it calls Kokkos::finalize, if its
    // constructor called Kokkos::initialize.
  }
  return 0;
}
