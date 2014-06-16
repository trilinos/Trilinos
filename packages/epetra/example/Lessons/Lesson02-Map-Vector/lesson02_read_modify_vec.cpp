/*!
\example lesson02_read_modify_vec.cpp
\brief Read and modify the entries of a vector (Epetra_Vector),
  using local indices.

\ref Epetra_Lesson02 explains this example in detail.
*/

// This defines useful macros like HAVE_MPI, which is defined if and
// only if Epetra was built with MPI enabled.
#include <Epetra_config.h>

#ifdef HAVE_MPI
#  include <mpi.h>
// Epetra's wrapper for MPI_Comm.  This header file only exists if
// Epetra was built with MPI enabled.
#  include <Epetra_MpiComm.h>
#else
#  include <Epetra_SerialComm.h>
#endif // HAVE_MPI

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_Version.h>

void
exampleRoutine (const Epetra_Comm& comm,
                std::ostream& out)
{
  using std::endl;

  // Print out the Epetra software version.
  if (comm.MyPID () == 0) {
    out << Epetra_Version () << endl << endl;
  }

  // The type of global indices.  You could just set this to int,
  // but we want the example to work for Epetra64 as well.
#ifdef EPETRA_NO_32BIT_GLOBAL_INDICES
  // Epetra was compiled only with 64-bit global index support, so use
  // 64-bit global indices.
  typedef long long global_ordinal_type;
#else
  // Epetra was compiled with 32-bit global index support.  If
  // EPETRA_NO_64BIT_GLOBAL_INDICES is defined, it does not also
  // support 64-bit indices.
  typedef int global_ordinal_type;
#endif // EPETRA_NO_32BIT_GLOBAL_INDICES

  //////////////////////////////////////////////////////////////////////
  // Create an Epetra_Map
  //////////////////////////////////////////////////////////////////////

  // The total (global, i.e., over all MPI processes) number of
  // entries in the Map.
  //
  // For this example, we scale the global number of entries in the
  // Map with the number of MPI processes.  That way, you can run this
  // example with any number of MPI processes and every process will
  // still have a positive number of entries.
  const global_ordinal_type numGlobalEntries = comm.NumProc () * 5;

  // Index base of the Map.  We choose zero-based (C-style) indexing.
  const global_ordinal_type indexBase = 0;

  // Construct a Map that puts the same number of equations on each
  // MPI process.
  Epetra_Map contigMap (numGlobalEntries, indexBase, comm);

  //////////////////////////////////////////////////////////////////////
  // Create an Epetra_Vector
  //////////////////////////////////////////////////////////////////////

  // Create a Vector with the Map we created above.
  // This version of the constructor will fill in the vector with zeros.
  Epetra_Vector x (contigMap);

  //////////////////////////////////////////////////////////////////////
  // Fill the Vector with a single number, or with random numbers
  //////////////////////////////////////////////////////////////////////

  // Set all entries of x to 42.0.
  (void) x.PutScalar (42.0);

  // Print the norm of x.
  double theNorm = 0.0;
  (void) x.Norm2 (&theNorm);
  out << "Norm of x (all entries are 42.0): " << theNorm << endl;

  // Set the entries of x to (pseudo)random numbers.  Please don't
  // consider this a good parallel pseudorandom number generator.
  (void) x.Random ();

  // Print the norm of x.
  (void) x.Norm2 (&theNorm);
  out << "Norm of x (random numbers): " << theNorm << endl;

  //////////////////////////////////////////////////////////////////////
  // Read the entries of the Vector
  //////////////////////////////////////////////////////////////////////

  {
    const int localLength = x.MyLength ();

    // Count the local number of entries less than 0.5.
    // Use local indices to access the entries of x_data.
    int localCount = 0;
    for (int localIndex = 0; localIndex < localLength; ++localIndex) {
      if (x[localIndex] < 0.5) {
        ++localCount;
      }
    }

    int globalCount = 0;
    (void) comm.SumAll (&localCount, &globalCount, 1);

    // Find the total number of entries less than 0.5,
    // over all processes in the Vector's communicator.
    out << "x has " << globalCount << " entr"
        << (globalCount != 1 ? "ies" : "y")
        << " less than 0.5." << endl;
  }

  //////////////////////////////////////////////////////////////////////
  // Modify the entries of the Vector
  //////////////////////////////////////////////////////////////////////

  {
    // Use local indices to access the entries of x_data.
    const int localLength = x.MyLength ();
    for (int localIndex = 0; localIndex < localLength; ++localIndex) {
      // Add the value of the local index to every entry of x.
      x[localIndex] += static_cast<double> (localIndex);
    }

    // Print the norm of x.
    theNorm = 0.0;
    (void) x.Norm2 (&theNorm);
    out << "Norm of x (modified random numbers): " << theNorm << endl;
  }
}

//
// The same main() driver routine as in the previous Epetra lesson.
//
int
main (int argc, char *argv[])
{
  using std::cout;
  using std::endl;

#ifdef HAVE_MPI
  MPI_Init (&argc, &argv);
  Epetra_MpiComm comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif // HAVE_MPI

  if (comm.MyPID () == 0) {
    cout << "Total number of processes: " << comm.NumProc () << endl;
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
