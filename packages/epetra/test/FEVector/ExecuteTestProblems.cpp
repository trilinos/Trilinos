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

  cout << "after globalAssemble"<<endl;
  A.Print(cout);

  return(ierr);
}

