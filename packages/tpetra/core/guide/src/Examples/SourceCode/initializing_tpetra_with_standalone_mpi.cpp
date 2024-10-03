// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Your code is an existing MPI code, so it presumably
// includes mpi.h directly.
#include <mpi.h>
#include <Teuchos_DefaultMpiComm.hpp> // wrapper for MPI_Comm
#include <Tpetra_Version.hpp> // Tpetra version string

int
main (int argc, char *argv[])
{
  using Teuchos::Comm;
  using Teuchos::MpiComm;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // It's bad form to ignore the error codes returned by MPI functions,
  // but we do so here for brevity.
  (void) MPI_Init(&argc, &argv);

  // This code takes the place of whatever you do to get an MPI_Comm.
  MPI_Comm yourComm = MPI_COMM_WORLD;

  // Wrap the MPI_Comm.
  RCP<const Comm<int> > comm(new MpiComm<int>(yourComm));

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

  // Since you called MPI_Init, you are responsible for calling
  // MPI_Finalize.
  (void) MPI_Finalize();

  // This tells the Trilinos test framework that the test passed.
  if (myRank == 0) {
    std::cout << "End Result: TEST PASSED" << std::endl;
  }
  return 0;
}

