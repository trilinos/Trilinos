// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*!
\example lesson01_mpi_on_its_own.cpp
\brief Initialization example for an existing MPI application that wants to use Tpetra.

\ref Tpetra_Lesson01 gives a full description of this example.
*/

//
// This example shows how to wrap the MPI_Comm (MPI communicator) that
// you are using, so that Tpetra can use it as well.  it includes MPI
// initialization, getting a Teuchos::Comm communicator, and printing
// out Tpetra version information.
//

// Your code is an existing MPI code, so it presumably includes mpi.h directly.
#include <mpi.h>

#include <Teuchos_DefaultMpiComm.hpp> // wrapper for MPI_Comm
#include <Tpetra_Version.hpp> // Tpetra version string

//
// ... Your other include files go here ...
//

// Do something with the given communicator.  In this case, we just
// print Tpetra's version to stdout on Process 0 in the given
// communicator.
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
  // things like std::endl or Teuchos::RCP.
  using std::cout;
  using std::endl;
  using Teuchos::Comm;
  using Teuchos::MpiComm;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // We assume that your code calls MPI_Init.  It's bad form
  // to ignore the error codes returned by MPI functions, but
  // we do so here for brevity.
  (void) MPI_Init (&argc, &argv);

  // This code takes the place of whatever you do to get an MPI_Comm.
  MPI_Comm yourComm = MPI_COMM_WORLD;

  {
    // Never create Tpetra objects at main scope.  Their destructors
    // must be called before MPI_Finalize and Kokkos::finalize are
    // called.

    // If your code plans to use MPI on its own, as well as through
    // Trilinos, consider giving Trilinos a copy of your MPI_Comm
    // (created via MPI_Comm_dup) rather than your MPI_Comm directly.
    // Trilinos may in the future duplicate the MPI_Comm
    // automatically, but it does not currently do this.  Duplicating
    // the MPI_Comm is not necessary, but may make it easier for you
    // to overlap asynchronous communication operations performed by
    // Trilinos with those performed by your code.

    // Wrap the MPI_Comm.  If you wrap it in this way, you are telling
    // Trilinos that you are responsible for calling MPI_Comm_free on
    // your MPI_Comm after use, if necessary.  (It's not necessary for
    // MPI_COMM_WORLD.)  There is a way to tell Trilinos to call
    // MPI_Comm_free itself; we don't show it here.  (It involves
    // passing the result of Teuchos::opaqueWrapper to MpiComm's
    // constructor.)

    RCP<const Comm<int> > comm (new MpiComm<int> (yourComm));

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
  }

  // If you need to call MPI_Comm_free on your MPI_Comm, now would be
  // the time to do so, before calling MPI_Finalize.  You may also
  // automate this process; ask the tutorial presenter for more
  // information.

  // Since you called MPI_Init, you are responsible for calling
  // MPI_Finalize.
  (void) MPI_Finalize ();
  return 0;
}
