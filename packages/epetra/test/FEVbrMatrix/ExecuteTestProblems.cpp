#include "Epetra_BLAS.h"
#include "ExecuteTestProblems.h"
#include "Epetra_Comm.h"
#include "Epetra_Vector.h"
#include <Epetra_FEVector.h>
  int MatrixTests(const Epetra_Map & Map, const Epetra_LocalMap & LocalMap, int NumVectors,
		      bool verbose)
  {
    const Epetra_Comm & Comm = Map.Comm();
    int ierr = 0;

    /* get ID of this processor */

    int MyPID   = Comm.MyPID();


    Epetra_MultiVector A(LocalMap, NumVectors);

    return(ierr);
  }

int quad1(const Epetra_Map& map, bool verbose)
{
  const Epetra_Comm & Comm = map.Comm();
  int numProcs = Comm.NumProc();
  int localProc = Comm.MyPID();

  //Set up a simple finite-element mesh containing 2-D quad elements, 1 per proc.
  //
  //   *-----*-----*-----*
  //  0|    2|    4|    6|
  //   | 0   | 1   | np-1|
  //   |     |     |     |
  //   *-----*-----*-----*
  //  1     3     5     7
  //
  //In the above drawing, 'np' means num-procs. node-numbers are to the
  //lower-left of each node (*).
  //
  //Each processor owns element 'localProc', and each processor owns
  //nodes 'localProc*2' and 'localProc*2+1' except for the last processor,
  //which also owns the last two nodes.
  //
  //There will be 3 degrees-of-freedom per node, so each element-matrix is
  //of size 12x12.
  //

  int myFirstNode = localProc*2;
  int myLastNode = localProc*2+1;
  if (localProc == numProcs-1) {
    myLastNode += 2;
  }

  int numMyNodes = myLastNode - myFirstNode + 1;
  int* myNodes = new int[numMyNodes];
  int i, j, ierr;
  for(i=0; i<numMyNodes; ++i) {
    myNodes[i] = myFirstNode + i;
  }

  int dofPerNode = 3; //degrees-of-freedom per node
  int indexBase = 0;

  Epetra_BlockMap blkMap(-1, numMyNodes, myNodes, dofPerNode,
			 indexBase, Comm);

  int rowLengths = 4; //each element-matrix will have 4 block-columns.
                      //the rows of the assembled matrix will be longer than
                      //this, but we don't need to worry about that because the
                      //VbrMatrix will add memory as needed. For a real
                      //application where efficiency is a concern, better
                      //performance would be obtained by giving a more accurate
                      //row-length here.

  Epetra_FEVbrMatrix A(Copy, blkMap, rowLengths);

  int nodesPerElem = 4;
  int* elemNodes = new int[nodesPerElem];
  for(i=0; i<nodesPerElem; ++i) elemNodes[i] = myFirstNode+i;

  int elemMatrixDim = nodesPerElem*dofPerNode;
  int len = elemMatrixDim*elemMatrixDim;
  double* elemMatrix = new double[len];

  //In an actual finite-element problem, we would calculate and fill 
  //meaningful element stiffness matrices. But for this simple matrix assembly
  //test, we're just going to fill our element matrix with 1.0's. This will
  //make it easy to see whether the matrix is correct after it's assembled.

  for(i=0; i<len; ++i) elemMatrix[i] = 1.0;

  //For filling in the matrix block-entries, we would ordinarily have to
  //carefully copy, or set pointers to, appropriate sections of the
  //element-matrix. But for this simple case we know that the element-matrix
  //is all 1's, so we'll just set our block-entry pointer to point to the
  //beginning of the element-matrix and leave it at that.
  //Note that the matrix class will refer to dofPerNode X dofPerNode (==9)
  //positions in the memory pointed to by 'blockEntry'.

  double* blockEntry = elemMatrix;

  //The element-matrix is a 4x4 (nodesPerElem X nodesPerElem) matrix of
  //3x3 block-entries. We'll now load our element-matrix into the global
  //matrix by looping over it and loading block-entries individually.

  for(i=0; i<nodesPerElem; ++i) {
    int blkrow = myFirstNode+i;
    EPETRA_TEST_ERR( A.BeginInsertGlobalValues(blkrow, nodesPerElem, elemNodes),
		    ierr);

    for(j=0; j<nodesPerElem; ++j) {
      EPETRA_TEST_ERR( A.SubmitBlockEntry( blockEntry, dofPerNode,
					  dofPerNode, dofPerNode), ierr);
    }

    EPETRA_TEST_ERR( A.EndSubmitEntries(), ierr);
  }

  EPETRA_TEST_ERR( A.GlobalAssemble(), ierr);

  if (verbose) {
    cout << "after globalAssemble"<<endl;
  }
  A.Print(cout);

  delete [] elemMatrix;
  delete [] myNodes;
  delete [] elemNodes;

  return(0);
}

int quad2(const Epetra_Map& map, bool verbose)
{
  const Epetra_Comm & Comm = map.Comm();
  int numProcs = Comm.NumProc();
  int localProc = Comm.MyPID();

  //Set up a simple finite-element mesh containing 2-D quad elements, 2 per proc.
  //(This test is similar to quad1() above, except for having 2 elements-per-proc
  // rather than 1.)
  //
  //   *-----*-----*-----*-------*
  //  0|    2|    4|    6|      8|
  //   | 0   | 1   | 2   | 2*np-1|
  //   |     |     |     |       |
  //   *-----*-----*-----*-------*
  //  1     3     5     7       9
  //
  //In the above drawing, 'np' means num-procs. node-numbers are to the
  //lower-left of each node (*).
  //
  //Each processor owns element 'localProc' and 'localProc+1', and each processor
  //owns nodes 'localProc*4' through 'localProc*4+3' except for the last
  //processor, which also owns the last two nodes in the mesh.
  //
  //There will be 3 degrees-of-freedom per node, so each element-matrix is
  //of size 12x12.
  //

  int myFirstNode = localProc*4;
  int myLastNode = localProc*4+3;
  if (localProc == numProcs-1) {
    myLastNode += 2;
  }

  int numMyElems = 2;

  int numMyNodes = myLastNode - myFirstNode + 1;
  int* myNodes = new int[numMyNodes];
  int i, j, ierr;
  for(i=0; i<numMyNodes; ++i) {
    myNodes[i] = myFirstNode + i;
  }

  int dofPerNode = 3; //degrees-of-freedom per node
  int indexBase = 0;

  Epetra_BlockMap blkMap(-1, numMyNodes, myNodes, dofPerNode,
			 indexBase, Comm);

  int rowLengths = 4; //each element-matrix will have 4 block-columns.
                      //the rows of the assembled matrix will be longer than
                      //this, but we don't need to worry about that because the
                      //VbrMatrix will add memory as needed. For a real
                      //application where efficiency is a concern, better
                      //performance would be obtained by giving a more accurate
                      //row-length here.

  Epetra_FEVbrMatrix A(Copy, blkMap, rowLengths);

  int nodesPerElem = 4;
  int* elemNodes = new int[nodesPerElem];

  int elemMatrixDim = nodesPerElem*dofPerNode;
  int len = elemMatrixDim*elemMatrixDim;
  double* elemMatrix = new double[len];

  //In an actual finite-element problem, we would calculate and fill 
  //meaningful element stiffness matrices. But for this simple matrix assembly
  //test, we're just going to fill our element matrix with 1.0's. This will
  //make it easy to see whether the matrix is correct after it's assembled.

  for(i=0; i<len; ++i) elemMatrix[i] = 1.0;

  //For filling in the matrix block-entries, we would ordinarily have to
  //carefully copy, or set pointers to, appropriate sections of the
  //element-matrix. But for this simple case we know that the element-matrix
  //is all 1's, so we'll just set our block-entry pointer to point to the
  //beginning of the element-matrix and leave it at that.
  //Note that the matrix class will refer to dofPerNode X dofPerNode (==9)
  //positions in the memory pointed to by 'blockEntry'.

  double* blockEntry = elemMatrix;

  //Each element-matrix is a 4x4 (nodesPerElem X nodesPerElem) matrix of
  //3x3 block-entries. We'll now load our element-matrices into the global
  //matrix by looping over them and loading block-entries individually.

  int firstNode = myFirstNode;
  for(int el=0; el<numMyElems; ++el) {
    for(i=0; i<nodesPerElem; ++i) {
   
      for(int n=0; n<nodesPerElem; ++n) elemNodes[n] = firstNode+n;

      int blkrow = firstNode+i;
      EPETRA_TEST_ERR( A.BeginInsertGlobalValues(blkrow, nodesPerElem, elemNodes),
		       ierr);

      for(j=0; j<nodesPerElem; ++j) {
	EPETRA_TEST_ERR( A.SubmitBlockEntry( blockEntry, dofPerNode,
					     dofPerNode, dofPerNode), ierr);
      }

      int this_err = A.EndSubmitEntries();
      if (this_err < 0) {
	cerr << "error in quad2, A.EndSubmitEntries(): " << this_err << endl;
	return(this_err);
      }
    }

    firstNode += 2;
  }

  EPETRA_TEST_ERR( A.GlobalAssemble(), ierr);

  if (verbose) {
    cout << "after globalAssemble"<<endl;
  }
  A.Print(cout);

  //now let's make sure that we can perform a matvec...
  Epetra_FEVector x(blkMap), y(blkMap);
  x.PutScalar(1.0);

  EPETRA_TEST_ERR( A.Multiply(false, x, y), ierr);

  if (verbose) {
    cout << "quad2, y:"<<endl;
  }
  y.Print(cout);

  delete [] elemMatrix;
  delete [] myNodes;
  delete [] elemNodes;

  return(0);
}

int MultiVectorTests(const Epetra_Map & Map, int NumVectors, bool verbose)
{
  const Epetra_Comm & Comm = Map.Comm();
  int ierr = 0, i, j;
  
  /* get number of processors and the name of this processor */
  
  int numprocs = Comm.NumProc();
  int MyPID   = Comm.MyPID();
  
  // Construct FEVbrMatrix
  
  if (verbose) cout << "constructing Epetra_FEVbrMatrix" << endl;

  //
  //we'll set up a tri-diagonal matrix.
  //

  int numGlobalRows = Map.NumGlobalElements();
  int numLocalRows = Map.NumMyElements();
  int minLocalRow = Map.MinMyGID();
  int rowLengths = 3;

  Epetra_FEVbrMatrix A(Copy, Map, rowLengths);
 
  if (verbose) {
    cout << "calling A.InsertGlobalValues with 1-D data array"<<endl;
  }

  int numCols = 3;
  int* ptIndices = new int[numCols];
  for(int k=0; k<numCols; ++k) {
    ptIndices[k] = minLocalRow+k;
  }

  double* values_1d = new double[numCols*numCols];
  for(j=0; j<numCols*numCols; ++j) {
    values_1d[j] = 3.0;
  }

  //For an extreme test, we'll have all processors sum into all rows.

  int minGID = Map.MinAllGID();

  //For now we're going to assume that there's just one point associated with
  //each GID (element).

  double* ptCoefs = new double[3];

  {for(i=0; i<numGlobalRows; ++i) {
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
      ptIndices[2] = minGID+i+2;
      ptCoefs[0] = 2.0;
      ptCoefs[1] = -1.0;
      ptCoefs[2] = -1.0;
      numCols = 3;
    }
    else {
      ptIndices[0] = minGID+i-2;
      ptIndices[1] = minGID+i-1;
      ptIndices[2] = minGID+i;
      ptCoefs[0] = -1.0;
      ptCoefs[1] = -1.0;
      ptCoefs[2] = 2.0;
      numCols = 3;
    }

    int row = minGID+i;

    EPETRA_TEST_ERR( A.BeginInsertGlobalValues(row, rowLengths, ptIndices), ierr);

    for(j=0; j<rowLengths; ++j) {
      EPETRA_TEST_ERR( A.SubmitBlockEntry(&(ptCoefs[j]), 1, 1, 1), ierr);
    }

    EPETRA_TEST_ERR( A.EndSubmitEntries(), ierr);

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
