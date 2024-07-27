// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_Map.hpp>
#include <Tpetra_Core.hpp>

int
main (int argc, char *argv[])
{
  using std::endl;
  using Teuchos::RCP;

  typedef Tpetra::Map<> map_type;
  typedef typename map_type::global_ordinal_type global_ordinal_type;

  Teuchos::oblackholestream blackHole;
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);
  {

    // Get a communicator corresponding to MPI_COMM_WORLD
    Teuchos::RCP<const Teuchos::Comm<int>> comm = Tpetra::getDefaultComm();
    const int myRank = comm->getRank();

    // Create sizes for map
    const size_t numLocalEntries = 5;
    const Tpetra::global_size_t numGlobalEntries = comm->getSize() * numLocalEntries;
    const global_ordinal_type indexBase = 0;

    // Create the Tpetra Map object
    RCP<const map_type> contigMap =
      rcp(new map_type(numGlobalEntries, indexBase, comm));

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! contigMap->isContiguous(), std::logic_error,
      "The supposedly contiguous Map isn't contiguous.");

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
