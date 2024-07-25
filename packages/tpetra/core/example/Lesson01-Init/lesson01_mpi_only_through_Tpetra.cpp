// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*!
\example lesson01_mpi_only_through_Tpetra.cpp
\brief Initialization example for a code that only uses MPI through Tpetra.

\ref Tpetra_Lesson01 gives a full description of this example.
*/

//
// This example includes MPI initialization, getting a Teuchos::Comm
// communicator, and printing out Tpetra version information.
//

#include <Tpetra_Core.hpp>
#include <Tpetra_Version.hpp>

// Do something with the given communicator.  In this case, we just
// print Tpetra's version to stdout on Process 0.
void
exampleRoutine (const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
{
  if (comm->getRank () == 0) {
    // On (MPI) Process 0, print out the Tpetra software version.
    std::cout << Tpetra::version () << std::endl << std::endl;
  }
}

int
main (int argc, char *argv[])
{
  // These "using" declarations make the code more concise, in that
  // you don't have to write the namespace along with the class or
  // object name.  This is especially helpful with commonly used
  // things like std::endl.
  using std::cout;
  using std::endl;

  // Start up MPI, if using MPI.  Trilinos doesn't have to be built
  // with MPI; it's called a "serial" build if you build without MPI.
  // Tpetra::ScopeGuard hides this implementation detail.
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);

  {
    // Never let Tpetra objects persist after either MPI_Finalize or
    // Kokkos::finalize has been called.  This is because the objects'
    // destructors may need to call MPI or Kokkos functions.  In
    // particular, never create Tpetra objects at main scope.

    // Get a pointer to the communicator object representing
    // MPI_COMM_WORLD.  The function knows whether or not we built with
    // MPI support.  If we didn't build with MPI, we'll get a
    // "communicator" with size 1, whose only process has rank 0.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();

    // Get my process' rank, and the total number of processes.
    // Equivalent to MPI_Comm_rank resp. MPI_Comm_size.
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();

    if (myRank == 0) {
      cout << "Total number of processes: " << numProcs << endl;
    }

    // Do something with the new communicator.
    exampleRoutine (comm);

    // This tells the Trilinos test framework that the test passed.
    if (myRank == 0) {
      cout << "End Result: TEST PASSED" << endl;
    }
    // ScopeGuard's destructor calls MPI_Finalize, if its constructor
    // called MPI_Init.  Likewise, it calls Kokkos::finalize, if its
    // constructor called Kokkos::initialize.
  }

  // You don't have to do anything here!  Just return from main().
  // Isn't that helpful?
  return 0;
}
