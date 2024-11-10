// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*!
\example lesson01_no_mpi.cpp
\brief Initialization example for integrating Tpetra into an existing non-MPI application.

\ref Tpetra_Lesson01 gives a full description of this example.
*/

// ... Your other include files go here ...
#include <Tpetra_Core.hpp>
#include <Tpetra_Version.hpp>

// Do something with the given communicator.  In this case, we just
// print Tpetra's version to stdout on Process 0.
void
exampleRoutine (const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
{
  if (comm->getRank () == 0) {
    // On Process 0, print out the Tpetra software version.
    std::cout << Tpetra::version () << std::endl << std::endl;
  }
}

int
main (int argc, char *argv[])
{
  using std::cout;
  using std::endl;
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Tpetra's default communicator will use a 1-process comm.
#ifdef HAVE_TPETRA_MPI
  Tpetra::ScopeGuard tpetraScope (&argc, &argv, MPI_COMM_SELF);
#else
  // Not building with MPI, so default comm won't use MPI.
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
#endif // HAVE_TPETRA_MPI
  {
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();

    // With a single-process communicator, the rank is always 0, and
    // the number of processes is always 1.
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();

    if (myRank == 0) {
      cout << "Total number of processes: " << numProcs << endl;
    }

    // Test the two assertions in the previous comment.
    // TEUCHOS_TEST_FOR_EXCEPTION is a macro defined in the Teuchos
    // package that takes three arguments: a bool expression, an
    // exception to throw if the expression evaluates to true, and a
    // message (interpreted as if it follows a "<<" after an
    // std::ostream) to include in the exception.  The macro includes
    // useful line number and file information in the exception
    // message, as well as a place where you can set a breakpoint in a
    // debugger right before the exception is thrown.

    TEUCHOS_TEST_FOR_EXCEPTION
      (myRank != 0, std::logic_error,
       "This is a single-MPI-process example, but the calling process' "
       "rank is " << myRank << " != 0.  Please report this bug.");

    TEUCHOS_TEST_FOR_EXCEPTION
      (numProcs != 1, std::logic_error,
       "This is a single-MPI-process example, but the number of "
       "processes in the Teuchos::Comm is " << numProcs << " != 1.  "
       "Please report this bug.");

    // Do something with the new communicator.
    exampleRoutine (comm);

    // This tells the Trilinos test framework that the test passed.
    if (myRank == 0) {
      cout << "End Result: TEST PASSED" << endl;
    }
  }
  return 0;
}
