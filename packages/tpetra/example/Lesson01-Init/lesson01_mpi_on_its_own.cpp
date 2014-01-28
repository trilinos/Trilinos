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
// ... Your other include files go here ...
#include <Teuchos_DefaultMpiComm.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Version.hpp>
#include <Teuchos_oblackholestream.hpp>

void
exampleRoutine (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                std::ostream& out)
{
  // Print out the Tpetra software version information.
  out << Tpetra::version() << std::endl << std::endl;
}

int
main (int argc, char *argv[])
{
  // These "using" declarations make the code more concise, in that
  // you don't have to write the namespace along with the class or
  // object name.  This is especially helpful with commonly used
  // things like std::endl or Teuchos::RCP.
  using std::endl;
  using Teuchos::Comm;
  using Teuchos::MpiComm;
  using Teuchos::opaqueWrapper;
  using Teuchos::RCP;
  using Teuchos::rcp;

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


  // Wrap the MPI_Comm.  If you wrap it in this way, you are telling
  // Trilinos that you are responsible for calling MPI_Comm_free on
  // your MPI_Comm after use, if necessary.  (It's not necessary for
  // MPI_COMM_WORLD.)  There is a way to tell Trilinos to call
  // MPI_Comm_free itself; we don't show it here.

  RCP<const Comm<int> > comm = rcp (new MpiComm<int> (opaqueWrapper (yourComm)));

  // Get my process' rank, and the total number of processes.
  // Equivalent to MPI_Comm_rank resp. MPI_Comm_size.  We don't
  // actually use numProcs in this code, so I've commented it out to
  // avoid a compiler warning for an unused variable.
  const int myRank = comm->getRank();
  //const int numProcs = comm->getSize();

  // The stream to which to write output.  Only MPI Process 0
  // gets to write to stdout; the other MPI processes get a
  // "black hole stream" which discards output (see above).
  Teuchos::oblackholestream blackHole;
  std::ostream& out = (myRank == 0) ? std::cout : blackHole;

  // We have a communicator and an output stream.
  // Let's do something with them!
  exampleRoutine (comm, out);

  // This tells the Trilinos test framework that the test passed.
  if (myRank == 0) {
    std::cout << "End Result: TEST PASSED" << std::endl;
  }

  // If you need to call MPI_Comm_free on your MPI_Comm, now would
  // be the time to do so, before calling MPI_Finalize.  You may also
  // automate this process; ask the tutorial presenter for more information.

  // Since you called MPI_Init, you are responsible for calling MPI_Finalize.
  (void) MPI_Finalize ();
  return 0;
}
