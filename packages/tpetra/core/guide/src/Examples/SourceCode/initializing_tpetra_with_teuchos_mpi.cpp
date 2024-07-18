// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_Core.hpp>
#include <Tpetra_Version.hpp>

  int
main (int argc, char *argv[])
{

  Tpetra::ScopeGuard tpetraScope(&argc, &argv);
  {
    // Get a communicator corresponding to MPI_COMM_WORLD
    Teuchos::RCP<const Teuchos::Comm<int>> comm = Tpetra::getDefaultComm();

    // Get my process' rank, and the total number of processes.
    // Equivalent to MPI_Comm_rank resp. MPI_Comm_size.
    const int myRank = comm->getRank();
    const int numProcs = comm->getSize();

    if (myRank == 0) {
      std::cout << "Total number of processes: " << numProcs << std::endl;
    }

    if (myRank == 0) {
      // On (MPI) Process 0, print out the Tpetra software version.
      std::cout << Tpetra::version() << std::endl << std::endl;
    }

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
