#include "Epetra_BLAS.h"
#include "ExecuteTestProblems.h"
#include "Epetra_Comm.h"
#include "Epetra_Vector.h"
  int MatrixTests(const Epetra_BlockMap & Map, const Epetra_LocalMap & LocalMap, int NumVectors,
		      bool verbose)
  {
    const Epetra_Comm & Comm = Map.Comm();
    int ierr = 0;
    double *residual = new double[NumVectors];

    /* get ID of this processor */

    int MyPID   = Comm.MyPID();


    Epetra_MultiVector A(LocalMap, NumVectors);

    delete [] residual;
    
    return(ierr);
  }

int MultiVectorTests(const Epetra_BlockMap & Map, int NumVectors, bool verbose)
{
  const Epetra_Comm & Comm = Map.Comm();
  int ierr = 0;
  
  /* get number of processors and the name of this processor */
  
  // int NumProc = Comm.getNumProc();
  int MyPID   = Comm.MyPID();
  
  // Construct FEVector
  
  if (verbose&&MyPID==0) cout << "constructing Epetra_FEVector" << endl;

  Epetra_FEVector A(Map);
 
  //For an extreme test, we'll have each processor sum-in a 1.0 for All
  //global ids.

  int minGID = Map.MinAllGID();
  int maxGID = Map.MaxAllGID();
  int numGlobalIDs = Map.NumGlobalElements();

  //For now we're going to have just one point associated with
  //each GID (element).

  int* ptIndices = new int[numGlobalIDs];
  double* ptCoefs = new double[numGlobalIDs];

  {for(int i=0; i<numGlobalIDs; ++i) {
    ptIndices[i] = minGID+i;
    ptCoefs[i] = 1.0;
  }}

  if (verbose&&MyPID==0) {
    cout << "calling A.SumIntoGlobalValues with " << numGlobalIDs << " values"<<endl;
  }
  EPETRA_TEST_ERR( A.SumIntoGlobalValues(numGlobalIDs, ptIndices, ptCoefs), ierr);

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

  delete [] ptIndices;
  delete [] ptCoefs;

  return(ierr);
}

int fevec1(Epetra_Comm& Comm, bool verbose)
{
  int ierr = 0;
  int NumGlobalRows = 4;
  int indexBase = 0;
  Epetra_Map Map(NumGlobalRows, indexBase, Comm);

  int Numprocs = Comm.NumProc();
  int MyPID   = Comm.MyPID();

  if (Numprocs != 2) return(0);

  int NumLocalRows = Map.NumMyElements();
  int MinLocalRow = Map.MinMyGID();

  int NumCols = 3;
  int* Indices = new int[NumCols];
 
  double* Values = new double[NumCols];
 
// Create vectors

  Epetra_FEVector b(Map);
  Epetra_FEVector x0(Map);
 
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

  EPETRA_TEST_ERR( b.SumIntoGlobalValues(NumCols, Indices, Values),
		   ierr);

  EPETRA_TEST_ERR( b.GlobalAssemble(), ierr);

  if (verbose&&MyPID==0) {
    cout << "b:"<<endl;
  }

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

  return(0);
}

int fevec2(Epetra_Comm& Comm, bool verbose)
{
  int ierr = 0;
  int NumGlobalRows = 12;
  int elemSize = 3;
  int indexBase = 0;
  Epetra_BlockMap Map(NumGlobalRows, elemSize, indexBase, Comm);

  int Numprocs = Comm.NumProc();
  int MyPID   = Comm.MyPID();

  if (Numprocs != 2) return(0);

  int NumLocalRows = Map.NumMyElements();
  int MinLocalRow = Map.MinMyGID();

  int NumCols = 3;
  int* Indices = new int[NumCols];
  int* numValuesPerID = new int[NumCols];
  for(int i=0; i<NumCols; ++i) {
    numValuesPerID[i] = elemSize;
  }
 
  double* Values = new double[NumCols*elemSize];
 
// Create vectors

  Epetra_FEVector b(Map);
  Epetra_FEVector x0(Map);
 
// source terms
  NumCols = 2;
 
  if(MyPID==0)  // indices corresponding to element 0 on processor 0
  {
    Indices[0] = 0;
    Indices[1] = 3;
 
    Values[0] = 1./2.;
    Values[1] = 1./2.;
    Values[2] = 1./2.;
    Values[3] = 1./2.;
    Values[4] = 1./2.;
    Values[5] = 1./2.;

   }
   else
   {
    Indices[0] = 1;
    Indices[1] = 2;
 
    Values[0] = 0;
    Values[1] = 0;
    Values[2] = 0;
    Values[3] = 0;
    Values[4] = 0;
    Values[5] = 0;
   }

  EPETRA_TEST_ERR( b.SumIntoGlobalValues(NumCols, Indices,
					 numValuesPerID, Values),
		   ierr);

  EPETRA_TEST_ERR( b.GlobalAssemble(), ierr);

  if (verbose&&MyPID==0) {
    cout << "b:"<<endl;
  }

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
  delete [] numValuesPerID;

  return(0);
}
