/*!
\example lesson05_redistribution.cpp
\brief Epetra parallel data redistribution, using Map and Export.

\ref Epetra_Lesson05 explains this example in detail.
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

#include <Epetra_CrsMatrix.h>
#include <Epetra_Export.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_Version.h>

#include <sstream>
#include <stdexcept>


// The type of global indices.  You could just set this to int,
// but we want the example to work for Epetra64 as well.
#ifdef EPETRA_NO_32BIT_GLOBAL_INDICES
// Epetra was compiled only with 64-bit global index support,
// so use 64-bit global indices.
#    define EXAMPLE_USES_64BIT_GLOBAL_INDICES 1
typedef long long global_ordinal_type;
#else
#  ifdef EPETRA_NO_64BIT_GLOBAL_INDICES
// Epetra was compiled only with 32-bit global index support,
// so use 32-bit global indices.
typedef int global_ordinal_type;
#  else
#    define EXAMPLE_USES_64BIT_GLOBAL_INDICES 1
// Epetra was compiled with both 64-bit and 32-bit global index
// support.  Use 64-bit global indices for maximum generality.
typedef long long global_ordinal_type;
#  endif // EPETRA_NO_64BIT_GLOBAL_INDICES
#endif // EPETRA_NO_32BIT_GLOBAL_INDICES


// Create and return a pointer to an example CrsMatrix, with row
// distribution over the given Map.  The caller is responsible for
// freeing the result.
Epetra_CrsMatrix*
createCrsMatrix (const Epetra_Map& map)
{
  const Epetra_Comm& comm = map.Comm ();

  // Create an Epetra_CrsMatrix using the Map, with dynamic allocation.
  Epetra_CrsMatrix* A = new Epetra_CrsMatrix (Copy, map, 3);

  // The list of global indices owned by this MPI process.

  const global_ordinal_type* myGblElts = NULL;
  global_ordinal_type numGblElts = 0;
#ifdef EXAMPLE_USES_64BIT_GLOBAL_INDICES
  myGblElts = map.MyGlobalElements64 ();
  numGblElts = map.NumGlobalElements64 ();
#else
  myGblElts = map.MyGlobalElements ();
  numGblElts = map.NumGlobalElements ();
#endif // EXAMPLE_USES_64BIT_GLOBAL_INDICES

  // The number of global indices owned by this MPI process.
  const int numMyElts = map.NumMyElements ();

  // In general, tests like this really should synchronize across all
  // processes.  However, the likely cause for this case is a
  // misconfiguration of Epetra, so we expect it to happen on all
  // processes, if it happens at all.
  if (numMyElts > 0 && myGblElts == NULL) {
    throw std::logic_error ("Failed to get the list of global indices");
  }

  // Local error code for use below.
  int lclerr = 0;

  // Fill the sparse matrix, one row at a time.
  double tempVals[3];
  global_ordinal_type tempGblInds[3];
  for (int i = 0; i < numMyElts; ++i) {
    // A(0, 0:1) = [2, -1]
    if (myGblElts[i] == 0) {
      tempVals[0] = 2.0;
      tempVals[1] = -1.0;
      tempGblInds[0] = myGblElts[i];
      tempGblInds[1] = myGblElts[i] + 1;
      if (lclerr == 0) {
        lclerr = A->InsertGlobalValues (myGblElts[i], 2, tempVals, tempGblInds);
      }
      if (lclerr != 0) {
        break;
      }
    }
    // A(N-1, N-2:N-1) = [-1, 2]
    else if (myGblElts[i] == numGblElts - 1) {
      tempVals[0] = -1.0;
      tempVals[1] = 2.0;
      tempGblInds[0] = myGblElts[i] - 1;
      tempGblInds[1] = myGblElts[i];
      if (lclerr == 0) {
        lclerr = A->InsertGlobalValues (myGblElts[i], 2, tempVals, tempGblInds);
      }
      if (lclerr != 0) {
        break;
      }
    }
    // A(i, i-1:i+1) = [-1, 2, -1]
    else {
      tempVals[0] = -1.0;
      tempVals[1] = 2.0;
      tempVals[2] = -1.0;
      tempGblInds[0] = myGblElts[i] - 1;
      tempGblInds[1] = myGblElts[i];
      tempGblInds[2] = myGblElts[i] + 1;
      if (lclerr == 0) {
        lclerr = A->InsertGlobalValues (myGblElts[i], 3, tempVals, tempGblInds);
      }
      if (lclerr != 0) {
        break;
      }
    }
  }

  // If any process failed to insert at least one entry, throw.
  int gblerr = 0;
  (void) comm.MaxAll (&lclerr, &gblerr, 1);
  if (gblerr != 0) {
    if (A != NULL) {
      delete A;
    }
    throw std::runtime_error ("Some process failed to insert an entry.");
  }

  // Tell the sparse matrix that we are done adding entries to it.
  gblerr = A->FillComplete ();
  if (gblerr != 0) {
    if (A != NULL) {
      delete A;
    }
    std::ostringstream os;
    os << "A->FillComplete() failed with error code " << gblerr << ".";
    throw std::runtime_error (os.str ());
  }

  return A;
}


void
example (const Epetra_Comm& comm)
{
  // The global number of rows in the matrix A to create.  We scale
  // this relative to the number of (MPI) processes, so that no matter
  // how many MPI processes you run, every process will have 10 rows.
  const global_ordinal_type numGblElts = 10 * comm.NumProc ();
  // The global min global index in all the Maps here.
  const global_ordinal_type indexBase = 0;

  // Local error code for use below.
  //
  // In the ideal case, we would use this to emulate behavior like
  // that of Haskell's Maybe in the context of MPI.  That is, if one
  // process experiences an error, we don't want to abort early and
  // cause the other processes to deadlock on MPI communication
  // operators.  Rather, we want to chain along the local error state,
  // until we reach a point where it's natural to pass along that
  // state with other processes.  For example, if one is doing an
  // MPI_Allreduce anyway, it makes sense to pass along one more bit
  // of information: whether the calling process is in a local error
  // state.  Epetra's interface doesn't let one chain the local error
  // state in this way, so we use extra collectives below to propagate
  // that state.  The code below uses very conservative error checks;
  // typical user code would not need to be so conservative and could
  // therefore avoid all the all-reduces.
  int lclerr = 0;

  // Construct a Map that is global (not locally replicated), but puts
  // all the equations on MPI Proc 0.
  const int procZeroMapNumLclElts = (comm.MyPID () == 0) ?
    numGblElts :
    static_cast<global_ordinal_type> (0);
  Epetra_Map procZeroMap (numGblElts, procZeroMapNumLclElts, indexBase, comm);

  // Construct a Map that puts approximately the same number of
  // equations on each processor.
  Epetra_Map globalMap (numGblElts, indexBase, comm);

  // Create a sparse matrix using procZeroMap.
  Epetra_CrsMatrix* A = createCrsMatrix (procZeroMap);
  if (A == NULL) {
    lclerr = 1;
  }

  // Make sure that sparse matrix creation succeeded.  Normally you
  // don't have to check this; we are being extra conservative because
  // this example is also a test.  Even though the matrix's rows live
  // entirely on Process 0, the matrix is nonnull on all processes in
  // its Map's communicator.
  int gblerr = 0;
  (void) comm.MaxAll (&lclerr, &gblerr, 1);
  if (gblerr != 0) {
    throw std::runtime_error ("createCrsMatrix returned NULL on at least one "
                              "process.");
  }

  //
  // We've created a sparse matrix whose rows live entirely on MPI
  // Process 0.  Now we want to distribute it over all the processes.
  //

  // Redistribute the matrix.  Since both the source and target Maps
  // are one-to-one, we could use either an Import or an Export.  If
  // only the source Map were one-to-one, we would have to use an
  // Import; if only the target Map were one-to-one, we would have to
  // use an Export.  We do not allow redistribution using Import or
  // Export if neither source nor target Map is one-to-one.

  // Make an export object with procZeroMap as the source Map, and
  // globalMap as the target Map.  The Export type has the same
  // template parameters as a Map.  Note that Export does not depend
  // on the Scalar template parameter of the objects it
  // redistributes.  You can reuse the same Export for different
  // Tpetra object types, or for Tpetra objects of the same type but
  // different Scalar template parameters (e.g., Scalar=float or
  // Scalar=double).
  Epetra_Export exporter (procZeroMap, globalMap);

  // Make a new sparse matrix whose row map is the global Map.
  Epetra_CrsMatrix B (Copy, globalMap, 0);

  // Redistribute the data, NOT in place, from matrix A (which lives
  // entirely on Proc 0) to matrix B (which is distributed evenly over
  // the processes).
  //
  // Export() has collective semantics, so we must always call it on
  // all processes collectively.  This is why we don't select on
  // lclerr, as we do for the local operations above.
  lclerr = B.Export (*A, exporter, Insert);

  // Make sure that the Export succeeded.  Normally you don't have to
  // check this; we are being extra conservative because this example
  // example is also a test.  We test both min and max, since lclerr
  // may be negative, zero, or positive.
  gblerr = 0;
  (void) comm.MinAll (&lclerr, &gblerr, 1);
  if (gblerr != 0) {
    throw std::runtime_error ("Export() failed on at least one process.");
  }
  (void) comm.MaxAll (&lclerr, &gblerr, 1);
  if (gblerr != 0) {
    throw std::runtime_error ("Export() failed on at least one process.");
  }

  // FillComplete has collective semantics, so we must always call it
  // on all processes collectively.  This is why we don't select on
  // lclerr, as we do for the local operations above.
  lclerr = B.FillComplete ();

  // Make sure that FillComplete succeeded.  Normally you don't have
  // to check this; we are being extra conservative because this
  // example is also a test.  We test both min and max, since lclerr
  // may be negative, zero, or positive.
  gblerr = 0;
  (void) comm.MinAll (&lclerr, &gblerr, 1);
  if (gblerr != 0) {
    throw std::runtime_error ("B.FillComplete() failed on at least one process.");
  }
  (void) comm.MaxAll (&lclerr, &gblerr, 1);
  if (gblerr != 0) {
    throw std::runtime_error ("B.FillComplete() failed on at least one process.");
  }

  if (A != NULL) {
    delete A;
  }
}


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

  const int myRank = comm.MyPID ();
  const int numProcs = comm.NumProc ();

  if (myRank == 0) {
    // Print out the Epetra software version.
    cout << Epetra_Version () << endl << endl
         << "Total number of processes: " << numProcs << endl;
  }

  example (comm); // Run the whole example.

  // This tells the Trilinos test framework that the test passed.
  if (myRank == 0) {
    cout << "End Result: TEST PASSED" << endl;
  }

#ifdef HAVE_MPI
  (void) MPI_Finalize ();
#endif // HAVE_MPI

  return 0;
}
