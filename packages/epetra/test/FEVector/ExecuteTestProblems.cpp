#include "Epetra_BLAS.h"
#include "ExecuteTestProblems.h"
#include "Epetra_Comm.h"
#include "Epetra_Vector.h"
  int MatrixTests(const Epetra_BlockMap & Map, const Epetra_LocalMap & LocalMap, int NumVectors,
		      bool verbose)
  {
    const Epetra_Comm & Comm = Map.Comm();
    int ierr = 0, i;
    int IndexBase = 0;
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
  int ierr = 0, i;
  
  /* get number of processors and the name of this processor */
  
  // int NumProc = Comm.getNumProc();
  int MyPID   = Comm.MyPID();
  
  // Construct FEVector
  
  if (verbose) cout << "constructing Epetra_FEVector" << endl;

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

  if (verbose) {
    cout << "calling A.SumIntoGlobalValues with " << numGlobalIDs << " values"<<endl;
  }
  EPETRA_TEST_ERR( A.SumIntoGlobalValues(numGlobalIDs, ptIndices, ptCoefs), ierr);

  if (verbose) {
    cout << "calling A.GlobalAssemble()" << endl;
  }

  EPETRA_TEST_ERR( A.GlobalAssemble(), ierr );

  if (verbose) {
  cout << "after globalAssemble"<<endl;
  }
  A.Print(cout);

  return(ierr);
}

int fevec1(Epetra_Comm& Comm, bool verbose)
{
  int ierr;
  int NumGlobalRows = 4;
  int indexBase = 0;
  Epetra_Map Map(NumGlobalRows, indexBase, Comm);

  int Numprocs = Comm.NumProc();
  int MyPID   = Comm.MyPID();

  if (Numprocs != 2) return(0);

  int NumLocalRows = Map.NumMyElements();
  int MinLocalRow = Map.MinMyGID();

  double ElementArea = 0.5;

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

  if (verbose) {
    cout << "b:"<<endl;
  }

  b.Print(cout);

  x0 = b;

  if (verbose) {
    cout << "x:"<<endl;
  }

  x0.Print(cout);

  delete [] Values;
  delete [] Indices;

  return(0);
}
