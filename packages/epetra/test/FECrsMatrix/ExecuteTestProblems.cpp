#include "Epetra_BLAS.h"
#include "ExecuteTestProblems.h"
#include "Epetra_Comm.h"
#include "Epetra_Vector.h"
  int MatrixTests(const Epetra_Map & Map, const Epetra_LocalMap & LocalMap, int NumVectors,
		      bool verbose)
  {
    const Epetra_Comm & Comm = Map.Comm();
    int ierr = 0, i;
    int IndexBase = 0;

    /* get ID of this processor */

    int MyPID   = Comm.MyPID();


    Epetra_MultiVector A(LocalMap, NumVectors);

    return(ierr);
  }

int MultiVectorTests(const Epetra_Map & Map, int NumVectors, bool verbose)
{
  const Epetra_Comm & Comm = Map.Comm();
  int ierr = 0, i;
  
  /* get number of processors and the name of this processor */
  
  int numprocs = Comm.NumProc();
  int MyPID   = Comm.MyPID();
  
  // Construct FECrsMatrix
  
  if (verbose) cout << "constructing Epetra_FECrsMatrix" << endl;

  //
  //we'll set up a tri-diagonal matrix.
  //

  int numGlobalRows = Map.NumGlobalElements();
  int numLocalRows = Map.NumMyElements();
  int minLocalRow = Map.MinMyGID();
  int rowLengths = 3;

  Epetra_FECrsMatrix A(Copy, Map, rowLengths);
 
  if (verbose) {
    cout << "calling A.SumIntoGlobalValues with 1-D data array"<<endl;
  }

  int numCols = 3;
  int* ptIndices = new int[numCols];
  for(int k=0; k<numCols; ++k) {
    ptIndices[k] = minLocalRow+k;
  }

  double* values_1d = new double[numCols*numCols];
  for(int j=0; j<numCols*numCols; ++j) {
    values_1d[j] = 3.0;
  }

  EPETRA_TEST_ERR( A.SumIntoGlobalValues(numCols, ptIndices,
                                         numCols, ptIndices,
                                         values_1d,
                                         Epetra_FECrsMatrix::ROW_MAJOR),ierr);

  //For an extreme test, we'll have all processors sum into all rows.

  int minGID = Map.MinAllGID();

  //For now we're going to assume that there's just one point associated with
  //each GID (element).

  double* ptCoefs = new double[3];

  {for(int i=0; i<numGlobalRows; ++i) {
    if (i>0 && i<numGlobalRows-1) {
      ptIndices[0] = minGID+i-1;
      ptIndices[1] = minGID+i;
      ptIndices[2] = minGID+i+1;
      ptCoefs[0] = -1.0;
      ptCoefs[1] = 2.0;
      ptCoefs[2] = -1.0;
      numCols = 3;
    }
    else if (i == 0) {
      ptIndices[0] = minGID+i;
      ptIndices[1] = minGID+i+1;
      ptCoefs[0] = 2.0;
      ptCoefs[1] = -1.0;
      numCols = 2;
    }
    else {
      ptIndices[0] = minGID+i-1;
      ptIndices[1] = minGID+i;
      ptCoefs[0] = -1.0;
      ptCoefs[1] = 2.0;
      numCols = 2;
    }

    int row = minGID+i;

    if (i%2==0) {
      if (verbose) {
        cout << "calling A.SumIntoGlobalValues with "<<numCols<<" values"<<endl;
      }
      EPETRA_TEST_ERR( A.SumIntoGlobalValues(1, &row,
                                             numCols, ptIndices,
                                             &ptCoefs,
                                             Epetra_FECrsMatrix::ROW_MAJOR), ierr );
    }
    else {
      if (verbose) {
        cout << "calling A.ReplaceGlobalValues with "<<numCols<<" values"<<endl;
      }
      EPETRA_TEST_ERR( A.ReplaceGlobalValues(1, &row,
                                             numCols, ptIndices,
                                             &ptCoefs,
                                             Epetra_FECrsMatrix::ROW_MAJOR), ierr );
    }
  }}

  if (verbose) {
    cout << "calling A.GlobalAssemble()" << endl;
  }

  EPETRA_TEST_ERR( A.GlobalAssemble(), ierr );

  if (verbose) {
  cout << "after globalAssemble"<<endl;
  A.Print(cout);
  }

  delete [] values_1d;
  delete [] ptIndices;
  delete [] ptCoefs;

  return(ierr);
}

