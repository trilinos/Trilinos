/*!
\example lesson01_mpi_on_its_own.cpp
\brief Initialization example for an existing MPI application that wants to use Epetra.

\ref Epetra_Lesson01 gives a full description of this example.
*/

//
// This example shows how to wrap the MPI_Comm (MPI communicator) that
// you are using, so that Epetra can use it as well.  it includes MPI
// initialization, wrapping your MPI_Comm in an Epetra communicator
// wrapper, and printing out Epetra version information.
//

// This defines useful macros like HAVE_MPI, which is defined if and
// only if Epetra was built with MPI enabled.
#include <Epetra_config.h>

#ifdef HAVE_MPI
// Your code is an existing MPI code, so it presumably includes mpi.h directly.
#  include <mpi.h>
// Epetra's wrapper for MPI_Comm.  This header file only exists if
// Epetra was built with MPI enabled.
#  include <Epetra_MpiComm.h>
#else
#  error "This example requires MPI in order to build."
#endif // HAVE_MPI

#include <Epetra_Version.h>

//
// ... Your other include files go here ...
//


// Do something with the given communicator.  In this case, we just
// print Epetra's version to the given output stream, on Process 0.
void
exampleRoutine (const Epetra_Comm& comm,
                std::ostream& out)
{
  if (comm.MyPID () == 0) {
    // On (MPI) Process 0, print out the Epetra software version.
    out << Epetra_Version () << std::endl << std::endl;
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

  // We assume that your code calls MPI_Init.  It's bad form
  // to ignore the error codes returned by MPI functions, but
  // we do so here for brevity.
  (void) MPI_Init (&argc, &argv);

  // This code takes the place of whatever you do to get an MPI_Comm.
  MPI_Comm yourComm = MPI_COMM_WORLD;

  // If your code plans to use MPI on its own, as well as through
  // Trilinos, you should strongly consider giving Trilinos a copy
  // of your MPI_Comm (created via MPI_Comm_dup).  Trilinos may in
  // the future duplicate the MPI_Comm automatically, but it does
  // not currently do this.

  // Wrap the MPI_Comm.  You are responsible for calling MPI_Comm_free
  // on your MPI_Comm after use, if necessary.  (It's not necessary or
  // legal to do this for built-in communicators like MPI_COMM_WORLD
  // or MPI_COMM_SELF.)
  Epetra_MpiComm comm (yourComm);

  // Epetra_Comm has methods that wrap basic MPI functionality.
  // MyPID() is equivalent to MPI_Comm_rank; it returns my process'
  // rank.  NumProc() is equivalent to MPI_Comm_size; it returns the
  // total number of processes in the communicator.
  const int myRank = comm.MyPID ();
  const int numProcs = comm.NumProc ();

  if (myRank == 0) {
    cout << "Total number of processes: " << numProcs << endl;
  }

  // Do something with the new Epetra communicator.
  exampleRoutine (comm, cout);

  // This tells the Trilinos test framework that the test passed.
  if (myRank == 0) {
    cout << "End Result: TEST PASSED" << endl;
  }

  // If you need to call MPI_Comm_free on your MPI_Comm, now would be
  // the time to do so, before calling MPI_Finalize.

  // Since you called MPI_Init, you are responsible for calling
  // MPI_Finalize after you are done using MPI.
  (void) MPI_Finalize ();
  return 0;
}
