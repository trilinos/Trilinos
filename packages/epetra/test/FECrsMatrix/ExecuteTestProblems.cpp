#include "Epetra_BLAS.h"
#include "ExecuteTestProblems.h"
#include "Epetra_Comm.h"
#include "Epetra_Vector.h"

int MatrixTests(const Epetra_Map & Map, const Epetra_LocalMap & LocalMap,
		int NumVectors, bool verbose)
  {
    const Epetra_Comm & Comm = Map.Comm();
    int ierr = 0;

    /* get ID of this processor */

    int MyPID   = Comm.MyPID();


    Epetra_MultiVector A(LocalMap, NumVectors);

    return(ierr);
  }

int Drumm1(const Epetra_Map& map, bool verbose)
{
  //Simple 2-element problem (element as in "finite-element") from
  //Clif Drumm. Two triangular elements, one per processor, as shown
  //here:
  //
  //   *----*
  //  3|\  2|
  //   | \  |
  //   | 0\1|
  //   |   \|
  //   *----*
  //  0    1
  //
  //Element 0 on processor 0, element 1 on processor 1.
  //Processor 0 will own nodes 0,1 and processor 1 will own nodes 2,3.
  //Each processor will pass a 3x3 element-matrix to Epetra_FECrsMatrix.
  //After GlobalAssemble(), the matrix should be as follows:
  //
  //         row 0: 2  1  0  1
  //proc 0   row 1: 1  4  1  2
  //----------------------------------
  //         row 2: 0  1  2  1
  //proc 1   row 3: 1  2  1  4
  //

  int numProcs = map.Comm().NumProc();
  int localProc = map.Comm().MyPID();

  if (numProcs != 2) return(0);

  int indexBase = 0, ierr = 0;

  if (localProc == 0) {
    int numMyNodes = 2;
    int* myNodes = new int[numMyNodes];
    myNodes[0] = 0;
    myNodes[1] = 1;

    double* values = new double[9];
    values[0] = 2.0;
    values[1] = 1.0;
    values[2] = 1.0;
    values[3] = 1.0;
    values[4] = 2.0;
    values[5] = 1.0;
    values[6] = 1.0;
    values[7] = 1.0;
    values[8] = 2.0;

    Epetra_Map Map(-1, numMyNodes, myNodes, indexBase, map.Comm());

    delete [] myNodes;
    numMyNodes = 3;
    myNodes = new int[numMyNodes];
    myNodes[0] = 0;
    myNodes[1] = 1;
    myNodes[2] = 3;

    int rowLengths = 3;
    Epetra_FECrsMatrix A(Copy, Map, rowLengths);

    EPETRA_TEST_ERR( A.SumIntoGlobalValues(numMyNodes, myNodes,
			  numMyNodes, myNodes,
			  values, Epetra_FECrsMatrix::ROW_MAJOR),ierr);

    EPETRA_TEST_ERR( A.GlobalAssemble(), ierr );

    A.Print(cout);

    delete [] myNodes;
    delete [] values;
  }
  else {
    int numMyNodes = 2;
    int* myNodes = new int[numMyNodes];
    myNodes[0] = 2;
    myNodes[1] = 3;

    double* values = new double[9];
    values[0] = 2.0;
    values[1] = 1.0;
    values[2] = 1.0;
    values[3] = 1.0;
    values[4] = 2.0;
    values[5] = 1.0;
    values[6] = 1.0;
    values[7] = 1.0;
    values[8] = 2.0;

    Epetra_Map Map(-1, numMyNodes, myNodes, indexBase, map.Comm());

    int rowLengths = 3;
    Epetra_FECrsMatrix A(Copy, Map, rowLengths);

    delete [] myNodes;
    numMyNodes = 3;
    myNodes = new int[numMyNodes];
    myNodes[0] = 1;
    myNodes[1] = 2;
    myNodes[2] = 3;

    EPETRA_TEST_ERR( A.SumIntoGlobalValues(numMyNodes, myNodes,
			  numMyNodes, myNodes,
			  values, Epetra_FECrsMatrix::ROW_MAJOR),ierr);

    EPETRA_TEST_ERR( A.GlobalAssemble(), ierr );

    A.Print(cout);

    delete [] myNodes;
    delete [] values;
  }

  return(0);
}

int Drumm2(const Epetra_Map& map, bool verbose)
{
  //Simple 2-element problem (element as in "finite-element") from
  //Clif Drumm. Two triangular elements, one per processor, as shown
  //here:
  //
  //   *----*
  //  3|\  2|
  //   | \  |
  //   | 0\1|
  //   |   \|
  //   *----*
  //  0    1
  //
  //Element 0 on processor 0, element 1 on processor 1.
  //Processor 0 will own nodes 0,1,3 and processor 1 will own node 2.
  //Each processor will pass a 3x3 element-matrix to Epetra_FECrsMatrix.
  //After GlobalAssemble(), the matrix should be as follows:
  //
  //         row 0: 2  1  0  1
  //proc 0   row 1: 1  4  1  2
  //         row 2: 0  1  2  1
  //----------------------------------
  //proc 1   row 3: 1  2  1  4
  //

  int numProcs = map.Comm().NumProc();
  int localProc = map.Comm().MyPID();

  if (numProcs != 2) return(0);

  int indexBase = 0, ierr = 0;

  if (localProc == 0) {
    int numMyNodes = 3;
    int* myNodes = new int[numMyNodes];
    myNodes[0] = 0;
    myNodes[1] = 1;
    myNodes[2] = 3;

    double* values = new double[9];
    values[0] = 2.0;
    values[1] = 1.0;
    values[2] = 1.0;
    values[3] = 1.0;
    values[4] = 2.0;
    values[5] = 1.0;
    values[6] = 1.0;
    values[7] = 1.0;
    values[8] = 2.0;

    Epetra_Map Map(-1, numMyNodes, myNodes, indexBase, map.Comm());

    int rowLengths = 3;
    Epetra_FECrsMatrix A(Copy, Map, rowLengths);

    EPETRA_TEST_ERR( A.SumIntoGlobalValues(numMyNodes, myNodes,
			  numMyNodes, myNodes,
			  values, Epetra_FECrsMatrix::ROW_MAJOR),ierr);

    EPETRA_TEST_ERR( A.GlobalAssemble(), ierr );

    A.Print(cout);

    //now let's make sure we can do a matvec with this matrix.
    Epetra_Vector x(Map), y(Map);
    x.PutScalar(1.0);
    EPETRA_TEST_ERR( A.Multiply(false, x, y), ierr);

    if (verbose) {
      cout << "y = A*x, x=1.0's"<<endl;
    }

    y.Print(cout);

    delete [] myNodes;
    delete [] values;
  }
  else {
    int numMyNodes = 1;
    int* myNodes = new int[numMyNodes];
    myNodes[0] = 2;

    double* values = new double[9];
    values[0] = 2.0;
    values[1] = 1.0;
    values[2] = 1.0;
    values[3] = 1.0;
    values[4] = 2.0;
    values[5] = 1.0;
    values[6] = 1.0;
    values[7] = 1.0;
    values[8] = 2.0;

    Epetra_Map Map(-1, numMyNodes, myNodes, indexBase, map.Comm());

    int rowLengths = 3;
    Epetra_FECrsMatrix A(Copy, Map, rowLengths);

    delete [] myNodes;
    numMyNodes = 3;
    myNodes = new int[numMyNodes];
    myNodes[0] = 1;
    myNodes[1] = 2;
    myNodes[2] = 3;

    EPETRA_TEST_ERR( A.SumIntoGlobalValues(numMyNodes, myNodes,
			  numMyNodes, myNodes,
			  values, Epetra_FECrsMatrix::ROW_MAJOR),ierr);

    EPETRA_TEST_ERR( A.GlobalAssemble(), ierr );

    A.Print(cout);

    //now let's make sure we can do a matvec with this matrix.
    Epetra_Vector x(Map), y(Map);
    x.PutScalar(1.0);
    EPETRA_TEST_ERR( A.Multiply(false, x, y), ierr);

    y.Print(cout);

    delete [] myNodes;
    delete [] values;
  }

  return(0);
}

int MultiVectorTests(const Epetra_Map & Map, int NumVectors, bool verbose)
{
  const Epetra_Comm & Comm = Map.Comm();
  int ierr = 0;
  
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
  }
  A.Print(cout);

  delete [] values_1d;
  delete [] ptIndices;
  delete [] ptCoefs;

  return(ierr);
}

