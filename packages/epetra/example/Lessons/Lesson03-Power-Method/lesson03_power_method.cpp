/*!
\example lesson03_power_method.cpp
\brief Example showing how to fill and compute with a Epetra sparse matrix

\ref Epetra_Lesson03 explains this example in detail.
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
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_Version.h>

#include <sstream>
#include <stdexcept>

// Implementation of the power method for estimating the eigenvalue of
// maximum magnitude of a matrix.  This function returns the
// eigenvalue estimate.
//
// A: The sparse matrix or operator, as an Epetra_Operator.
// niters: Maximum number of iterations of the power method.
// tolerance: If the 2-norm of the residual A*x-lambda*x (for the
//   current eigenvalue estimate lambda) is less than this, stop
//   iterating.
// out: output stream to which to print the current status of the
//   power method.
double
powerMethod (const Epetra_Operator& A,
             const int niters,
             const double tolerance);

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

  // The number of rows and columns in the matrix.
  const global_ordinal_type numGlobalElements = 50;

  // Construct a Map that puts approximately the same number of
  // equations on each processor.
  const global_ordinal_type indexBase = 0;
  Epetra_Map map (numGlobalElements, indexBase, comm);

  // Get the list of global indices that this process owns.  In this
  // example, this is unnecessary, because we know that we created a
  // contiguous Map (see above).  (Thus, we really only need the min
  // and max global index on this process.)  However, in general, we
  // don't know what global indices the Map owns, so if we plan to add
  // entries into the sparse matrix using global indices, we have to
  // get the list of global indices this process owns.
  const int numMyElements = map.NumMyElements ();

  global_ordinal_type* myGlobalElements = NULL;

#ifdef EPETRA_NO_32BIT_GLOBAL_INDICES
  myGlobalElements = map.MyGlobalElements64 ();
#else
  myGlobalElements = map.MyGlobalElements ();
#endif // EPETRA_NO_32BIT_GLOBAL_INDICES

  // In general, tests like this really should synchronize across all
  // processes.  However, the likely cause for this case is a
  // misconfiguration of Epetra, so we expect it to happen on all
  // processes, if it happens at all.
  if (numMyElements > 0 && myGlobalElements == NULL) {
    throw std::logic_error ("Failed to get the list of global indices");
  }

  if (myRank == 0) {
    cout << endl << "Creating the sparse matrix" << endl;
  }

  // Create a Epetra sparse matrix whose rows have distribution given
  // by the Map.  The max number of entries per row is 3.
  Epetra_CrsMatrix A (Copy, map, 3);

  // Local error code for use below.
  int lclerr = 0;

  // Fill the sparse matrix, one row at a time.  InsertGlobalValues
  // adds entries to the sparse matrix, using global column indices.
  // It changes both the graph structure and the values.
  double tempVals[3];
  global_ordinal_type tempGblInds[3];
  for (int i = 0; i < numMyElements; ++i) {
    // A(0, 0:1) = [2, -1]
    if (myGlobalElements[i] == 0) {
      tempVals[0] = 2.0;
      tempVals[1] = -1.0;
      tempGblInds[0] = myGlobalElements[i];
      tempGblInds[1] = myGlobalElements[i] + 1;
      if (lclerr == 0) {
        lclerr = A.InsertGlobalValues (myGlobalElements[i], 2, tempVals, tempGblInds);
      }
      if (lclerr != 0) {
        break;
      }
    }
    // A(N-1, N-2:N-1) = [-1, 2]
    else if (myGlobalElements[i] == numGlobalElements - 1) {
      tempVals[0] = -1.0;
      tempVals[1] = 2.0;
      tempGblInds[0] = myGlobalElements[i] - 1;
      tempGblInds[1] = myGlobalElements[i];
      if (lclerr == 0) {
        lclerr = A.InsertGlobalValues (myGlobalElements[i], 2, tempVals, tempGblInds);
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
      tempGblInds[0] = myGlobalElements[i] - 1;
      tempGblInds[1] = myGlobalElements[i];
      tempGblInds[2] = myGlobalElements[i] + 1;
      if (lclerr == 0) {
        lclerr = A.InsertGlobalValues (myGlobalElements[i], 3, tempVals, tempGblInds);
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
    throw std::runtime_error ("Some process failed to insert an entry.");
  }

  // Tell the sparse matrix that we are done adding entries to it.
  gblerr = A.FillComplete ();
  if (gblerr != 0) {
    std::ostringstream os;
    os << "A.FillComplete() failed with error code " << gblerr << ".";
    throw std::runtime_error (os.str ());
  }

  // Number of iterations
  const int niters = 500;
  // Desired (absolute) residual tolerance
  const double tolerance = 1.0e-2;

  // Run the power method and report the result.
  double lambda = powerMethod (A, niters, tolerance);
  if (myRank == 0) {
    cout << endl << "Estimated max eigenvalue: " << lambda << endl;
  }

  //
  // Now we're going to change values in the sparse matrix and run the
  // power method again.
  //

  //
  // Increase diagonal dominance
  //
  if (myRank == 0) {
    cout << endl << "Increasing magnitude of A(0,0), solving again" << endl;
  }

  if (A.RowMap ().MyGID (0)) {
    // Get a copy of the row with with global index 0.  Modify the
    // diagonal entry of that row.  Submit the modified values to the
    // matrix.
    const global_ordinal_type gidOfFirstRow = 0;
    // Since 0 is a GID in the row Map on the calling process,
    // lidOfFirstRow is a valid LID of that GID in the row Map.
    const int lidOfFirstRow = A.RowMap ().LID (gidOfFirstRow);
    int numEntriesInRow = A.NumMyEntries (lidOfFirstRow);
    double* rowvals = new double [numEntriesInRow];
    global_ordinal_type* rowinds = new global_ordinal_type [numEntriesInRow];

    // Get a copy of the entries and column indices of global row 0.
    // Get global column indices, so that we can figure out which
    // entry corresponds to the diagonal entry in this row.  (The row
    // Map and column Map of the matrix may differ, so their local
    // indices may not be the same.)
    //
    // Note that it's legal (though we don't exercise it in this
    // example) for the row Map of the sparse matrix not to be one to
    // one.  This means that more than one process might own entries
    // in the first row.  In general, multiple processes might own the
    // (0,0) entry, so that the global A(0,0) value is really the sum
    // of all processes' values for that entry.  However, scaling the
    // entry by a constant factor distributes across that sum, so it's
    // OK to do so.
    if (lclerr == 0) {
      lclerr = A.ExtractGlobalRowCopy (gidOfFirstRow,
                                       numEntriesInRow, numEntriesInRow,
                                       rowvals, rowinds);
    }
    if (lclerr == 0) { // no error
      for (int i = 0; i < numEntriesInRow; ++i) {
        if (rowinds[i] == gidOfFirstRow) {
          // We have found the diagonal entry; modify it.
          rowvals[i] *= 10.0;
        }
      }
      // "Replace global values" means modify the values, but not the
      // structure of the sparse matrix.  If the specified columns
      // aren't already populated in this row on this process, then this
      // method fails and returns nonzero.  Since we have already called
      // FillComplete() on this matrix, we may not modify its graph
      // structure any more.
      if (lclerr == 0) {
        lclerr = A.ReplaceGlobalValues (gidOfFirstRow, numEntriesInRow,
                                        rowvals, rowinds);
      }
    }

    if (rowvals != NULL) {
      delete [] rowvals;
    }
    if (rowinds != NULL) {
      delete [] rowinds;
    }
  }

  // If the owning process(es) of global row 0 failed to replace the
  // one entry, throw.
  gblerr = 0;
  (void) comm.MaxAll (&lclerr, &gblerr, 1);
  if (gblerr != 0) {
    throw std::runtime_error ("One of the owning process(es) of global "
                              "row 0 failed to replace an entry.");
  }

  // Run the power method again.
  lambda = powerMethod (A, niters, tolerance);
  if (myRank == 0) {
    cout << endl << "Estimated max eigenvalue: " << lambda << endl;
  }

  // This tells the Trilinos test framework that the test passed.
  if (myRank == 0) {
    cout << "End Result: TEST PASSED" << endl;
  }

#ifdef HAVE_MPI
  (void) MPI_Finalize ();
#endif // HAVE_MPI

  return 0;
}


double
powerMethod (const Epetra_Operator& A,
             const int niters,
             const double tolerance)
{
  using std::cout;
  using std::endl;

  // An Operator doesn't have a Comm, but its domain Map does.
  const Epetra_Comm& comm = A.OperatorDomainMap ().Comm ();
  const int myRank = comm.MyPID ();

  // Create three vectors for iterating the power method.  Since the
  // power method computes z = A*q, q should be in the domain of A and
  // z should be in the range.  (Obviously the power method requires
  // that the domain and the range are equal, but it's a good idea to
  // get into the habit of thinking whether a particular vector
  // "belongs" in the domain or range of the matrix.)  The residual
  // vector "resid" is of course in the range of A.
  Epetra_Vector q (A.OperatorDomainMap ());
  Epetra_Vector z (A.OperatorRangeMap ());
  Epetra_Vector resid (A.OperatorRangeMap ());

  // Local error code for use below.
  int lclerr = 0;

  // Fill the iteration vector z with random numbers to start.  Don't
  // have grand expectations about the quality of our pseudorandom
  // number generator; this is usually good enough for eigensolvers.
  //
  // If this call fails, just let the code below finish before trying
  // to catch the error globally.  Ditto for other calls below.
  lclerr = z.Random ();

  // lambda: the current approximation of the eigenvalue of maximum magnitude.
  // normz: the 2-norm of the current iteration vector z.
  // residual: the 2-norm of the current residual vector "resid"
  double lambda = 0.0;
  double normz = 0.0;
  double residual = 0.0;

  const double zero = 0.0;
  const double one = 1.0;

  // How often to report progress in the power method.  Reporting
  // progress requires computing a residual which can be expensive.
  // However, if you don't compute the residual often enough, you
  // might keep iterating even after you've converged.
  const int reportFrequency = 10;

  // Do the power method, until the method has converged or the
  // maximum iteration count has been reached.
  for (int iter = 0; iter < niters; ++iter) {
    // If you feel confident that your code is correct, you may omit
    // everything having to do with lclerr, and just write the following:
    //
    // z.Norm2 (&normz);         // Compute the 2-norm of z
    // q.Scale (one / normz, z); // q := z / normz
    // A.Apply (q, z);           // z := A * q
    // q.Dot (z, &lambda);       // Approx. max eigenvalue

    lclerr = (lclerr == 0) ? z.Norm2 (&normz) : lclerr;
    lclerr = (lclerr == 0) ? q.Scale (one / normz, z) : lclerr;
    lclerr = (lclerr == 0) ? A.Apply (q, z) : lclerr;
    lclerr = (lclerr == 0) ? q.Dot (z, &lambda) : lclerr;

    // Compute and report the residual norm every reportFrequency
    // iterations, or if we've reached the maximum iteration count.
    if (iter % reportFrequency == 0 || iter + 1 == niters) {
      // If you feel confident that your code is correct, you may omit
      // everything having to do with lclerr, and just write the
      // following:
      //
      // resid.Update (one, z, -lambda, q, zero); // z := A*q - lambda*q
      // resid.Norm2 (&residual); // 2-norm of the residual vector

      lclerr = (lclerr == 0) ? resid.Update (one, z, -lambda, q, zero) : lclerr;
      lclerr = (lclerr == 0) ? resid.Norm2 (&residual) : lclerr;

      if (myRank == 0) {
        cout << "Iteration " << iter << ":" << endl
             << "- lambda = " << lambda << endl
             << "- ||A*q - lambda*q||_2 = " << residual << endl;
      }
    }
    if (residual < tolerance) {
      if (myRank == 0) {
        cout << "Converged after " << iter << " iterations" << endl;
      }
      break;
    } else if (iter + 1 == niters) {
      if (myRank == 0) {
        cout << "Failed to converge after " << niters << " iterations" << endl;
      }
      break;
    }
  }

  // If any process failed to insert at least one entry, throw.
  int gblerr = 0;
  (void) comm.MaxAll (&lclerr, &gblerr, 1);
  if (gblerr != 0) {
    throw std::runtime_error ("The power method failed in some way.");
  }

  return lambda;
}
