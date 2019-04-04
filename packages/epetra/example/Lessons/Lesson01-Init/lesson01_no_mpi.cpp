/*!
\example lesson01_no_mpi.cpp
\brief Initialization example for integrating Epetra into an existing non-MPI application.

\ref Epetra_Lesson01 gives a full description of this example.
*/

#include <Epetra_config.h>
// Wrapper for a "communicator" containing only one process.  This
// header file always exists, whether or not Epetra was built with MPI
// enabled.
#include <Epetra_SerialComm.h>
#include <Epetra_Version.h>

#include <cstdlib>
#include <sstream>
#include <stdexcept>

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
main (int /* argc */, char * /* argv */[])
{
  // These "using" declarations make the code more concise, in that
  // you don't have to write the namespace along with the class or
  // object name.  This is especially helpful with commonly used
  // things like std::endl.
  using std::cout;
  using std::endl;

  // Make a "serial" (non-MPI) communicator.  It doesn't actually
  // "communicate," because it only has one process, whose rank is
  // always 0.  Epetra_SerialComm is a subclass of Epetra_Comm, so you
  // may use it wherever an Epetra_Comm is required.
  Epetra_SerialComm comm;

  // Epetra_Comm has methods that wrap basic MPI functionality.
  // MyPID() is equivalent to MPI_Comm_rank, and NumProc() to
  // MPI_Comm_size.
  //
  // With a "serial" communicator, the rank is always 0, and the
  // number of processes is always 1.
  const int myRank = comm.MyPID ();
  const int numProcs = comm.NumProc ();

  // Test the two assertions in the previous comment.
  if (numProcs != 1) {
    std::ostringstream err;
    err << "This is a serial (non-MPI) example, but the number of processes "
        << "in the Epetra_Comm is " << numProcs << " != 1.  Please report "
        << "this bug.";
    throw std::logic_error (err.str ());
  }
  if (myRank != 0) {
    std::ostringstream err;
    err << "This is a serial (non-MPI) example, but the rank of the calling "
      "process in the Epetra_Comm is " << myRank << " != 0.  Please report "
      "this bug.";
    throw std::logic_error (err.str ());
  }

  // Do something with the new communicator.
  exampleRoutine (comm, cout);

  // This tells the Trilinos test framework that the test passed.
  if (myRank == 0) {
    cout << "End Result: TEST PASSED" << endl;
  }

  return 0;
}
