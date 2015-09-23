/*!
\example lesson01_mpi_only_through_Epetra.cpp
\brief Initialization example for a code that only uses MPI through Epetra.

\ref Epetra_Lesson01 gives a full description of this example.
*/

//
// This example includes conditional MPI initialization, getting an
// Epetra communicator wrapper, and printing out Epetra version
// information.
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
#  include <Epetra_SerialComm.h>
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

#ifdef HAVE_MPI
  // Start up MPI, if using MPI.  Trilinos doesn't have to be built
  // with MPI; it's called a "serial" build if you build without MPI.
  //
  // It's bad form to ignore the error codes returned by MPI
  // functions, but we do so here for brevity.
  (void) MPI_Init (&argc, &argv);

  // Wrap MPI_COMM_WORLD in an Epetra communicator wrapper.
  // Epetra_MpiComm is a subclass of Epetra_Comm, so you may use it
  // wherever an Epetra_Comm is required.
  Epetra_MpiComm comm (MPI_COMM_WORLD);
#else
  // Make a "serial" (non-MPI) communicator.  It doesn't actually
  // "communicate," because it only has one process, whose rank is
  // always 0.  Epetra_SerialComm is a subclass of Epetra_Comm, so you
  // may use it wherever an Epetra_Comm is required.
  Epetra_SerialComm comm;
#endif

  // Epetra_Comm has methods that wrap basic MPI functionality.
  // MyPID() is equivalent to MPI_Comm_rank, and NumProc() to
  // MPI_Comm_size.
  //
  // With a "serial" communicator, the rank is always 0, and the
  // number of processes is always 1.
  const int myRank = comm.MyPID ();
  const int numProcs = comm.NumProc ();

  if (myRank == 0) {
    cout << "Total number of processes: " << numProcs << endl;
  }

  // Do something with the new Epetra communicator.
  exampleRoutine (comm, cout);

  // This tells the Trilinos test framework that the test passed.
  if (comm.MyPID () == 0) {
    cout << "End Result: TEST PASSED" << endl;
  }

#ifdef HAVE_MPI
  // Since you called MPI_Init, you are responsible for calling
  // MPI_Finalize after you are done using MPI.
  (void) MPI_Finalize ();
#endif // HAVE_MPI

  return 0;
}
