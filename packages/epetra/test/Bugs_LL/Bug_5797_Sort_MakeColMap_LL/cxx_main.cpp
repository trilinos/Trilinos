// Program for testing Epetra64 implementation.
// Builds a Laplacian matrx using either Epetra or Epetra64.
// Multiplies it by the vector of all ones and checks norms.
//
// To run:
//    mpirun -np # test64.exe [n]
// where n is an optional argument specifying the number of matrix rows.
// Default n == 25.
//
// Two macro definitions below control behavior:
//    ITYPE:  int --> use Epetra
//            long long --> use Epetra64
//    OFFSET_EPETRA64:  Add this value to each row/column index.  Resulting
//                      matrix rows/columns are indexed from 
//                      OFFSET_EPETRA64 to OFFSET_EPETRA64+n-1.

#include <limits.h>
#define ITYPE long long
//#define OFFSET_EPETRA64 ((long long)(2)*(long long)INT_MAX)
//#define OFFSET_EPETRA64 (INT_MAX-5)
#define OFFSET_EPETRA64 0

#include <stdio.h>

#include "Epetra_ConfigDefs.h"

#ifdef EPETRA_MPI
#include <mpi.h>
#include "Epetra_MpiComm.h"
#define FINALIZE MPI_Finalize()
#else
#include "Epetra_SerialComm.h"
#define FINALIZE
#endif
 
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"

////////////////////////////////////////////////////////////////
// Build the following Laplacian matrix using Epetra or Epetra64
//
//    | 2 -1  ...    -1...           |
//    |-1  3 -1 ...     -1...        |
//    |        ...                   |
//    |-1       -1  4 -1...   -1...  |
//    |        ...                   |
//    |         -1            -1 3 -1|
//    |            -1           -1  2|
//
////////////////////////////////////////////////////////////////

#define MIN(a,b) ((a) < (b) ? (a) : (b))

int main(int narg, char *arg[])
{
  using std::cout;

#ifdef EPETRA_MPI  
  // Initialize MPI  
  MPI_Init(&narg,&arg);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  int me = comm.MyPID();
  int np = comm.NumProc();

  ITYPE nGlobalRows = 25;
  if (narg > 1) 
    nGlobalRows = (ITYPE) atol(arg[1]);

  bool verbose = (nGlobalRows < 26);

  // Linear map similar to Trilinos default, 
  // but want to allow adding OFFSET_EPETRA64 to the indices.
  int nMyRows = (int) (nGlobalRows / np + (nGlobalRows % np > me));
  ITYPE myFirstRow = (ITYPE)(me * (nGlobalRows / np) + MIN(nGlobalRows%np, me));
  ITYPE *myGlobalRows = new ITYPE[nMyRows];
  for (int i = 0; i < nMyRows; i++)
    myGlobalRows[i] = (ITYPE)i + myFirstRow + OFFSET_EPETRA64;
  Epetra_Map *rowMap = new Epetra_Map(-1, nMyRows, &myGlobalRows[0], 0, comm);
  if (verbose) rowMap->Print(cout);

  // Create an integer vector nnzPerRow that is used to build the Epetra Matrix.
  // nnzPerRow[i] is the number of entries for the ith local equation
  std::vector<int> nnzPerRow(nMyRows+1, 0);

  // Also create lists of the nonzeros to be assigned to processors.
  // To save programming time and complexity, these vectors are allocated 
  // bigger than they may actually be needed.
  std::vector<ITYPE> iv(5*nMyRows+1);
  std::vector<ITYPE> jv(5*nMyRows+1);
  std::vector<double> vv(5*nMyRows+1);

  // Generate the nonzeros for the Laplacian matrix.
  ITYPE nMyNonzeros = 0;
  for (ITYPE i = 0, myrowcnt = 0; i < nGlobalRows; i++) {
    if (rowMap->MyGID(i+OFFSET_EPETRA64)) { 
      int idegree = 0;
      // This processor owns this row; add nonzeros.
      if (i > 0) {
        iv[nMyNonzeros] = i + OFFSET_EPETRA64;
        jv[nMyNonzeros] = i-1 + OFFSET_EPETRA64;
        vv[nMyNonzeros] = -1;
        if (verbose)
          cout << "(" << iv[nMyNonzeros] << "," << jv[nMyNonzeros] << ")="
               << vv[nMyNonzeros] << " on processor " << me
               << " in " << myrowcnt << endl;
        nMyNonzeros++;
        nnzPerRow[myrowcnt]++;
        idegree++;
      }

      if (i < nGlobalRows - 1) {
        iv[nMyNonzeros] = i + OFFSET_EPETRA64;
        jv[nMyNonzeros] = i+1 + OFFSET_EPETRA64;
        vv[nMyNonzeros] = -1;
        if (verbose) 
          cout << "(" << iv[nMyNonzeros] << "," << jv[nMyNonzeros] << ")="
               << vv[nMyNonzeros] << " on processor " << me
               << " in " << myrowcnt << endl;
        nMyNonzeros++;
        nnzPerRow[myrowcnt]++;
        idegree++;
      }

      if (i+5 < nGlobalRows) {
        iv[nMyNonzeros] = i + OFFSET_EPETRA64;
        jv[nMyNonzeros] = i+5 + OFFSET_EPETRA64;
        vv[nMyNonzeros] = -1;
        if (verbose)
          cout << "(" << iv[nMyNonzeros] << "," << jv[nMyNonzeros] << ")="
               << vv[nMyNonzeros] << " on processor " << me
               << " in " << myrowcnt << endl;
        nMyNonzeros++;
        nnzPerRow[myrowcnt]++;
        idegree++;
      }

      if (i-5 >= 0) {
        iv[nMyNonzeros] = i + OFFSET_EPETRA64;
        jv[nMyNonzeros] = i-5 + OFFSET_EPETRA64;
        vv[nMyNonzeros] = -1;
        if (verbose)
          cout << "(" << iv[nMyNonzeros] << "," << jv[nMyNonzeros] << ")="
               << vv[nMyNonzeros] << " on processor " << me
               << " in " << myrowcnt << endl;
        nMyNonzeros++;
        nnzPerRow[myrowcnt]++;
        idegree++;
      }

      iv[nMyNonzeros] = i + OFFSET_EPETRA64;
      jv[nMyNonzeros] = i + OFFSET_EPETRA64;
      vv[nMyNonzeros] = idegree;
      if (verbose) 
        cout << "(" << iv[nMyNonzeros] << "," << jv[nMyNonzeros] << ")="
             << vv[nMyNonzeros] << " on processor " << me
             << " in " << myrowcnt << endl;
      nMyNonzeros++;
      nnzPerRow[myrowcnt]++;

      myrowcnt++;
    }
  }

  // Create an Epetra_Matrix
  // TODO:  Should be able to use StaticProfile == true, but there appears to
  // TODO:  be a bug in Epetra64 (or in my usage of it).  I filed bug 5791
  // TODO:  and will use StaticProfile == false for now.  When the bug is fixed,
  // TODO:  we should switch back to StaticProfile == true, as it will be more
  // TODO:  efficient.
  // #warning Once bug 5791 is fixed, change StaticProfile back to true.
  // -- Done.  Chetan.
  Epetra_CrsMatrix *A = new Epetra_CrsMatrix(Copy, *rowMap, &nnzPerRow[0], true);

  // Insert the nonzeros.
  int info;
  ITYPE sum = 0;
  for (int i=0; i < nMyRows; i++) {
    if (nnzPerRow[i]) {
      if (verbose) {
        cout << "InsertGlobalValus row " << iv[sum]
             << " count " << nnzPerRow[i] 
             << " cols " << jv[sum] << " " << jv[sum+1] << " ";
        if (nnzPerRow[i] == 3) cout << jv[sum+2];
        cout << endl;
      }
      info = A->InsertGlobalValues(iv[sum],nnzPerRow[i],&vv[sum],&jv[sum]);
      assert(info==0);
      sum += nnzPerRow[i];
    }
  }

  // Finish up
  info = A->FillComplete();
  assert(info==0);
  if (verbose) A->Print(cout);

  // Sanity test:  Product of matrix and vector of ones should have norm == 0
  // and max/min/mean values of 0
  Epetra_Vector sanity(A->RangeMap());
  Epetra_Vector sanityres(A->DomainMap());
  sanity.PutScalar(1.);
  A->Multiply(false, sanity, sanityres);

  double jjone, jjtwo, jjmax;
  sanityres.Norm1(&jjone);
  sanityres.Norm2(&jjtwo);
  sanityres.NormInf(&jjmax);
  if (me == 0)
    cout << "SanityTest norms 1/2/inf: " << jjone << " "
                                         << jjtwo << " " << jjmax << endl;

  bool test_failed = (jjone != 0) || (jjtwo != 0) || (jjmax != 0);

  sanityres.MinValue(&jjone);
  sanityres.MeanValue(&jjtwo);
  sanityres.MaxValue(&jjmax);
  if (me == 0)
    cout << "SanityTest values min/max/avg: " << jjone << " "
                                              << jjmax << " " << jjtwo << endl;

  test_failed = test_failed || (jjone != 0) || (jjtwo != 0) || (jjmax != 0);

  if (me == 0) {
    if(test_failed)
      cout << "Bug_5797_Sort_MakeColMap_LL tests FAILED" << endl;
  }

  delete A;
  delete rowMap;
  delete [] myGlobalRows;

  FINALIZE;
}

