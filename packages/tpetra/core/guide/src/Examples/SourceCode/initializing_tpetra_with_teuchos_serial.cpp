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
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::rcp;

    // Tpetra's default communicator will use a 1-process comm.
#ifdef HAVE_TPETRA_MPI
  Tpetra::ScopeGuard tpetraScope(&argc, &argv, MPI_COMM_SELF);
#else
  // Not building with MPI, so default comm won't use MPI.
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);
#endif // HAVE_TPETRA_MPI

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();

  // With a "serial" communicator, the rank is always 0,
  // and the number of processes is always 1.
  const int myRank = comm->getRank();
  const int numProcs = comm->getSize();

  if (myRank == 0) {
    std::cout << "Total number of processes: " << numProcs << std::endl;
  }

  if (comm->getRank () == 0) {
    // On Process 0, print out the Tpetra software version.
    std::cout << Tpetra::version() << std::endl << std::endl;
  }

  // This tells the Trilinos test framework that the test passed.
  if (myRank == 0) {
    std::cout << "End Result: TEST PASSED" << std::endl;
  }
  return 0;
}
