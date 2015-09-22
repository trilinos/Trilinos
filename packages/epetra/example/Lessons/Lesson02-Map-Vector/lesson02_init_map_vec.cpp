/*!
\example lesson02_init_map_vec.cpp
\brief Create data distributions (Epetra_Map), and create vectors
  (Epetra_Vector) distributed according to those distributions.

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

#include <stdexcept>

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
  // Create some Epetra_Map objects
  //////////////////////////////////////////////////////////////////////

  //
  // Epetra has local and global Maps.  Local maps describe objects
  // that are replicated over all participating MPI processes.  Global
  // maps describe distributed objects.  You can do imports and
  // exports between local and global maps; this is how you would turn
  // locally replicated objects into distributed objects and vice
  // versa.
  //

  // The total (global, i.e., over all MPI processes) number of
  // entries in the Map.  This has the same type as that of global
  // indices, so it can represent very large values if Epetra was
  // built with 64-bit global index support.
  //
  // For this example, we scale the global number of entries in the
  // Map with the number of MPI processes.  That way, you can run this
  // example with any number of MPI processes and every process will
  // still have a positive number of entries.
  const global_ordinal_type numGlobalEntries = comm.NumProc () * 5;

  // Tpetra can index the entries of a Map starting with 0 (C style),
  // 1 (Fortran style), or any base you want.  1-based indexing is
  // handy when interfacing with Fortran.  We choose 0-based indexing
  // here.  This also has the same type as that of global indices.
  const global_ordinal_type indexBase = 0;

  // Construct a Map that puts the same number of equations on each
  // (MPI) process.  The Epetra_Comm is passed in by value, but that's
  // OK, because Epetra_Comm has shallow copy semantics.  (Its copy
  // constructor and assignment operator do not call MPI_Comm_dup;
  // they just pass along the MPI_Comm.)
  Epetra_Map contigMap (numGlobalEntries, indexBase, comm);

  // contigMap is contiguous by construction.
  if (! contigMap.LinearMap ()) {
    throw std::logic_error ("The supposedly contiguous Map isn't contiguous.");
  }

  // Let's create a second Map.  It will have the same number of
  // global entries per process, but will distribute them differently,
  // in round-robin (1-D cyclic) fashion instead of contiguously.

  // We'll use the version of the Map constructor that takes, on each
  // MPI process, a list of the global indices in the Map belonging to
  // that process.  You can use this constructor to construct an
  // overlapping (also called "not 1-to-1") Map, in which one or more
  // entries are owned by multiple processes.  We don't do that here;
  // we make a nonoverlapping (also called "1-to-1") Map.
  const int numGblIndsPerProc = 5;
  global_ordinal_type* gblIndList = new global_ordinal_type [numGblIndsPerProc];

  const int numProcs = comm.NumProc ();
  const int myRank = comm.MyPID ();
  for (int k = 0; k < numGblIndsPerProc; ++k) {
    gblIndList[k] = myRank + k*numProcs;
  }

  Epetra_Map cyclicMap (numGlobalEntries, numGblIndsPerProc,
                        gblIndList, indexBase, comm);
  // The above constructor makes a deep copy of the input index list,
  // so it's safe to deallocate that list after this constructor
  // completes.
  if (gblIndList != NULL) {
    delete [] gblIndList;
    gblIndList = NULL;
  }

  // If there's more than one MPI process in the communicator,
  // then cyclicMap is definitely NOT contiguous.
  if (comm.NumProc () > 1 && cyclicMap.LinearMap ()) {
    throw std::logic_error ("The cyclic Map claims to be contiguous.");
  }

  // contigMap and cyclicMap should always be compatible.  However, if
  // the communicator contains more than 1 process, then contigMap and
  // cyclicMap are NOT the same.
  // if (! contigMap.isCompatible (*cyclicMap)) {
  //   throw std::logic_error ("contigMap should be compatible with cyclicMap, "
  //                           "but it's not.");
  // }
  if (comm.NumProc () > 1 && contigMap.SameAs (cyclicMap)) {
    throw std::logic_error ("contigMap should not be the same as cyclicMap.");
  }

  //////////////////////////////////////////////////////////////////////
  // We have maps now, so we can create vectors.
  //////////////////////////////////////////////////////////////////////

  // Create an Epetra_Vector with the contiguous Map we created above.
  // This version of the constructor will fill the vector with zeros.
  // The Vector constructor takes a Map by value, but that's OK,
  // because Epetra_Map has shallow copy semantics.  It uses reference
  // counting internally to avoid copying data unnecessarily.
  Epetra_Vector x (contigMap);

  // The copy constructor performs a deep copy.
  // x and y have the same Map.
  Epetra_Vector y (x);

  // Create a Vector with the 1-D cyclic Map.  Calling the constructor
  // with false for the second argument leaves the data uninitialized,
  // so that you can fill it later without paying the cost of
  // initially filling it with zeros.
  Epetra_Vector z (cyclicMap, false);

  // Set the entries of z to (pseudo)random numbers.  Please don't
  // consider this a good parallel pseudorandom number generator.
  (void) z.Random ();

  // Set the entries of x to all ones.
  (void) x.PutScalar (1.0);

  // Define some constants for use below.
  const double alpha = 3.14159;
  const double beta = 2.71828;
  const double gamma = -10.0;

  // x = beta*x + alpha*z
  //
  // This is a legal operation!  Even though the Maps of x and z are
  // not the same, their Maps are compatible.  Whether it makes sense
  // or not depends on your application.
  (void) x.Update (alpha, z, beta);

  (void) y.PutScalar (42.0); // Set all entries of y to 42.0
  // y = gamma*y + alpha*x + beta*z
  y.Update (alpha, x, beta, z, gamma);

  // Compute the 2-norm of y.
  //
  // The norm may have a different type than scalar_type.
  // For example, if scalar_type is complex, then the norm is real.
  // The ScalarTraits "traits class" gives us the type of the norm.
  double theNorm = 0.0;
  (void) y.Norm2 (&theNorm);

  // Print the norm of y on Proc 0.
  out << "Norm of y: " << theNorm << endl;
}

//
// The same main() driver routine as in the first Epetra lesson.
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
