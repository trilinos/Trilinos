#include "Epetra_BLAS.h"
#include "ExecuteTestProblems.h"
#include "Epetra_Comm.h"
#include "Epetra_Vector.h"
#include "Epetra_FEVector.h"

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

    if (verbose) {
    A.Print(cout);
    }

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

    if (verbose) {
    A.Print(cout);
    }

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

    if (verbose) {
    A.Print(cout);
    }

    //now let's make sure we can do a matvec with this matrix.
    Epetra_Vector x(Map), y(Map);
    x.PutScalar(1.0);
    EPETRA_TEST_ERR( A.Multiply(false, x, y), ierr);

    if (verbose&&localProc==0) {
      cout << "y = A*x, x=1.0's"<<endl;
    }

    if (verbose) {
    y.Print(cout);
    }

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

    if (verbose) {
    A.Print(cout);
    }

    //now let's make sure we can do a matvec with this matrix.
    Epetra_Vector x(Map), y(Map);
    x.PutScalar(1.0);
    EPETRA_TEST_ERR( A.Multiply(false, x, y), ierr);

    if (verbose) {
    y.Print(cout);
    }

    delete [] myNodes;
    delete [] values;
  }

  return(0);
}

int Drumm3(const Epetra_Map& map, bool verbose)
{
  const Epetra_Comm & Comm = map.Comm();
  int ierr = 0, i;

  /* get number of processors and the name of this processor */

  int Numprocs = Comm.NumProc();
  int MyPID   = Comm.MyPID();

  if (Numprocs != 2) return(0);

  int NumGlobalRows = 4;
  int IndexBase = 0;
  Epetra_Map Map(NumGlobalRows, IndexBase, Comm);

  // Construct FECrsMatrix

  int NumLocalRows = Map.NumMyElements();
  int MinLocalRow = Map.MinMyGID();
  int NumEntriesPerRow = 3;

  Epetra_FECrsMatrix A(Copy, Map, NumEntriesPerRow);

  double ElementArea = 0.5;
  
  int NumCols = 3;
  int* Indices = new int[NumCols];

  if(MyPID==0)  // indices corresponding to element 0 on processor 0
  {
    Indices[0] = 0;
    Indices[1] = 1;
    Indices[2] = 3;
  }
  else if(MyPID==1)  // indices corresponding to element 1 on processor 1
  {
    Indices[0] = 1;
    Indices[1] = 2;
    Indices[2] = 3;
  }

  double* Values = new double[NumCols*NumCols];

// removal term
  Values[0] = 2*ElementArea/12.;
  Values[1] = 1*ElementArea/12.;
  Values[2] = 1*ElementArea/12.;
  Values[3] = 1*ElementArea/12.;
  Values[4] = 2*ElementArea/12.;
  Values[5] = 1*ElementArea/12.;
  Values[6] = 1*ElementArea/12.;
  Values[7] = 1*ElementArea/12.;
  Values[8] = 2*ElementArea/12.;

  A.SumIntoGlobalValues(NumCols, Indices,
                        Values,
                        Epetra_FECrsMatrix::ROW_MAJOR);

  A.GlobalAssemble();

//  A.Print(cout);

// Create vectors for CG algorithm

  Epetra_FEVector* bptr = new Epetra_FEVector(A.RowMap());
  Epetra_FEVector* x0ptr = new Epetra_FEVector(A.RowMap());

  Epetra_FEVector& b = *bptr;
  Epetra_FEVector& x0 = *x0ptr;

  // source terms
  NumCols = 2;

  if(MyPID==0)  // indices corresponding to element 0 on processor 0
  {
    Indices[0] = 0;
    Indices[1] = 3;

    Values[0] = 1./2.;
    Values[1] = 1./2.;

   }
   else
   {
    Indices[0] = 1;
    Indices[1] = 2;

    Values[0] = 0;
    Values[1] = 0;
   }

  b.SumIntoGlobalValues(NumCols, Indices, Values);

  b.GlobalAssemble();

  if (verbose&&MyPID==0) cout << "b:" << endl;
  if (verbose) {
  b.Print(cout);
  }

  x0 = b;

  if (verbose&&MyPID==0) {
  cout << "x:"<<endl;
  }

  if (verbose) {
  x0.Print(cout);
  }

  delete [] Values;
  delete [] Indices;

  delete bptr;
  delete x0ptr;

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
  
  if (verbose&&MyPID==0) cout << "constructing Epetra_FECrsMatrix" << endl;

  //
  //we'll set up a tri-diagonal matrix.
  //

  int numGlobalRows = Map.NumGlobalElements();
  int numLocalRows = Map.NumMyElements();
  int minLocalRow = Map.MinMyGID();
  int rowLengths = 3;

  Epetra_FECrsMatrix A(Copy, Map, rowLengths);
 
  if (verbose&&MyPID==0) {
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
      if (verbose&&MyPID==0) {
        cout << "calling A.SumIntoGlobalValues with "<<numCols<<" values"<<endl;
      }
      EPETRA_TEST_ERR( A.SumIntoGlobalValues(1, &row,
                                             numCols, ptIndices,
                                             &ptCoefs,
                                             Epetra_FECrsMatrix::ROW_MAJOR), ierr );
    }
    else {
      if (verbose&&MyPID==0) {
        cout << "calling A.ReplaceGlobalValues with "<<numCols<<" values"<<endl;
      }
      EPETRA_TEST_ERR( A.ReplaceGlobalValues(1, &row,
                                             numCols, ptIndices,
                                             &ptCoefs,
                                             Epetra_FECrsMatrix::ROW_MAJOR), ierr );
    }
  }}

  if (verbose&&MyPID==0) {
    cout << "calling A.GlobalAssemble()" << endl;
  }

  EPETRA_TEST_ERR( A.GlobalAssemble(), ierr );

  if (verbose&&MyPID==0) {
  cout << "after globalAssemble"<<endl;
  }
  if (verbose) {
  A.Print(cout);
  }

  delete [] values_1d;
  delete [] ptIndices;
  delete [] ptCoefs;

  return(ierr);
}

