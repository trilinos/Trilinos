/*!
\example lesson01_mpi_only_through_Tpetra.cpp
\brief Initialization example for a code that only uses MPI through Tpetra.

\ref Tpetra_Lesson01 gives a full description of this example.
*/

//
// This example includes MPI initialization, getting a Teuchos::Comm
// communicator, and printing out Tpetra version information.
//

#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Version.hpp>
#include <Teuchos_GlobalMPISession.hpp>

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
  // GlobalMPISession hides this implementation detail.
  //
  // Note the third argument.  If you pass GlobalMPISession the
  // address of an std::ostream, it will print a one-line status
  // message with the rank on each MPI process.  This may be
  // undesirable if running with a large number of MPI processes.
  // You can avoid printing anything here by passing in either
  // NULL or the address of a Teuchos::oblackholestream.
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, NULL);

  // Get a pointer to the communicator object representing
  // MPI_COMM_WORLD.  getDefaultPlatform.getComm() doesn't create a
  // new object every time you call it; it just returns the same
  // communicator each time.  Thus, you can call it anywhere and get
  // the same communicator.  (This is handy if you don't want to pass
  // a communicator around everywhere, though it's always better to
  // parameterize your algorithms on the communicator.)
  //
  // "Tpetra::DefaultPlatform" knows whether or not we built with MPI
  // support.  If we didn't build with MPI, we'll get a "communicator"
  // with size 1, whose only process has rank 0.
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();

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

  // GlobalMPISession calls MPI_Finalize() in its destructor, if
  // appropriate.  You don't have to do anything here!  Just return
  // from main().  Isn't that helpful?
  return 0;
}
