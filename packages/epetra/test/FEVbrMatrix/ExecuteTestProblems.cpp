//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER


#include "Epetra_BLAS.h"
#include "ExecuteTestProblems.h"
#include "Epetra_Comm.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_FEVbrMatrix.h"
#include <Epetra_FEVector.h>

#include "../epetra_test_err.h"


int quad1(const Epetra_Map& map, bool verbose)
{
  const Epetra_Comm & Comm = map.Comm();
  int numProcs = Comm.NumProc();
  int localProc = Comm.MyPID();

  Comm.Barrier();
  if (verbose && localProc == 0) {
    cout << "====================== quad1 =============================="<<endl;
  }

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
  //of size 12x12. (4 nodes per element, X 3 dof per node)
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

  int rowLengths = 3; //each element-matrix will have 4 block-columns.
                      //the rows of the assembled matrix will be longer than
                      //this, but we don't need to worry about that because the
                      //VbrMatrix will add memory as needed. For a real
                      //application where efficiency is a concern, better
                      //performance would be obtained by giving more accurate
                      //row-lengths here.

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
      for(int ii=0; ii<dofPerNode*dofPerNode; ii++) {
        blockEntry[ii] = blkrow+elemNodes[j];
      }
      EPETRA_TEST_ERR( A.SubmitBlockEntry( blockEntry, dofPerNode,
					  dofPerNode, dofPerNode), ierr);
    }

    int err = A.EndSubmitEntries();
    if (err < 0) {
      cout << "quad1: error in A.EndSubmitEntries: "<<err<<endl;
      return(-1);
    }
  }

  EPETRA_TEST_ERR( A.GlobalAssemble(), ierr);

  if (verbose && localProc==0) {
    cout << "after globalAssemble"<<endl;
  }
  if (verbose) {
    A.Print(cout);
  }

  int numMyRows = A.NumMyRows();
  int correct_numMyRows = dofPerNode*numMyNodes;

  if (numMyRows != correct_numMyRows) {
    cout << "proc " << localProc << ", numMyRows("<<numMyRows<<") doesn't match"
	 << " correct_numMyRows("<<correct_numMyRows<<")."<<endl;
    return(-1);
  }

  int numMyNonzeros = A.NumMyNonzeros();
  int correct_numMyNonzeros = nodesPerElem*nodesPerElem*dofPerNode*dofPerNode;

  if (numProcs > 1) {
    if (localProc == numProcs-1) {
      correct_numMyNonzeros += dofPerNode*dofPerNode*4;
    }
    else if (localProc > 0) {
      correct_numMyNonzeros -= dofPerNode*dofPerNode*4;
    }
    else {
      //localProc==0 && numProcs > 1
      correct_numMyNonzeros -= dofPerNode*dofPerNode*8;
    }
  }

  if (numMyNonzeros != correct_numMyNonzeros) {
    cout << "proc " << localProc << ", numMyNonzeros(" << numMyNonzeros
	 <<") != correct_numMyNonzeros("<<correct_numMyNonzeros<<")"<<endl;
    return(-1);
  }


  delete [] elemMatrix;
  delete [] myNodes;
  delete [] elemNodes;

  Comm.Barrier();

  return(0);
}

int quad2(const Epetra_Map& map, bool verbose)
{
  const Epetra_Comm & Comm = map.Comm();
  int numProcs = Comm.NumProc();
  int localProc = Comm.MyPID();
  Comm.Barrier();
  if (verbose && localProc == 0) {
    cout << "====================== quad2 =============================="<<endl;
  }

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

  if (verbose && localProc==0) {
    cout << "after globalAssemble"<<endl;
  }

  if (verbose) {
    A.Print(cout);
  }

  //now let's make sure that we can perform a matvec...
  Epetra_FEVector x(blkMap, 1), y(blkMap, 1);
  x.PutScalar(1.0);

  EPETRA_TEST_ERR( A.Multiply(false, x, y), ierr);

  if (verbose && localProc==0) {
    cout << "quad2, y:"<<endl;
  }
  if (verbose) {
    y.Print(cout);
  }

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
  
  int MyPID   = Comm.MyPID();
  
  // Construct FEVbrMatrix
  
  if (verbose && MyPID==0) cout << "constructing Epetra_FEVbrMatrix" << endl;

  //
  //we'll set up a tri-diagonal matrix.
  //

  int numGlobalRows = Map.NumGlobalElements();
  int minLocalRow = Map.MinMyGID();
  int rowLengths = 3;

  Epetra_FEVbrMatrix A(Copy, Map, rowLengths);
 
  if (verbose && MyPID==0) {
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

int four_quads(const Epetra_Comm& Comm, bool preconstruct_graph, bool verbose)
{   
  if (verbose) {
    cout << "******************* four_quads ***********************"<<endl;
  }

  //This function assembles a matrix representing a finite-element mesh
  //of four 2-D quad elements. There are 9 nodes in the problem. The
  //same problem is assembled no matter how many processors are being used
  //(within reason). It may not work if more than 9 processors are used.
  //
  //  *------*------*
  // 6|     7|     8|
  //  | E2   | E3   |
  //  *------*------*
  // 3|     4|     5|
  //  | E0   | E1   |
  //  *------*------*
  // 0      1      2
  //
  //Nodes are denoted by * with node-numbers below and left of each node.
  //E0, E1 and so on are element-numbers.
  //
  //Each processor will contribute a sub-matrix of size 4x4, filled with 1's,
  //for each element. Thus, the coefficient value at position 0,0 should end up
  //being 1.0*numProcs, the value at position 4,4 should be 1.0*4*numProcs, etc.
  //
  //Depending on the number of processors being used, the locations of the
  //specific matrix positions (in terms of which processor owns them) will vary.
  //
  
  int numProcs = Comm.NumProc();
  
  int numNodes = 9;
  int numElems = 4;
  int numNodesPerElem = 4;

  int blockSize = 1;
  int indexBase = 0;

  //Create a map using epetra-defined linear distribution.
  Epetra_BlockMap map(numNodes, blockSize, indexBase, Comm);

  Epetra_CrsGraph* graph = NULL;

  int* nodes = new int[numNodesPerElem];
  int i, j, k, err = 0;

  if (preconstruct_graph) {
    graph = new Epetra_CrsGraph(Copy, map, 1);

    //we're going to fill the graph with indices, but remember it will only
    //accept indices in rows for which map.MyGID(row) is true.

    for(i=0; i<numElems; ++i) {
      switch(i) {
      case 0:
        nodes[0] = 0; nodes[1] = 1; nodes[2] = 4; nodes[3] = 3;
        break;
      case 1:
        nodes[0] = 1; nodes[1] = 2; nodes[2] = 5; nodes[3] = 4;
        break;
      case 2:
        nodes[0] = 3; nodes[1] = 4; nodes[2] = 7; nodes[3] = 6;
        break;
      case 3:
        nodes[0] = 4; nodes[1] = 5; nodes[2] = 8; nodes[3] = 7;
        break;
      }

      for(j=0; j<numNodesPerElem; ++j) {
        if (map.MyGID(nodes[j])) {
          err = graph->InsertGlobalIndices(nodes[j], numNodesPerElem,
                                           nodes);
          if (err<0) return(err);
        }
      }
    }

    EPETRA_CHK_ERR( graph->FillComplete() );
  }

  Epetra_FEVbrMatrix* A = NULL;

  if (preconstruct_graph) {
    A = new Epetra_FEVbrMatrix(Copy, *graph);
  }
  else {
    A = new Epetra_FEVbrMatrix(Copy, map, 1);
  }

  //EPETRA_CHK_ERR( A->PutScalar(0.0) );

  double* values_1d = new double[numNodesPerElem*numNodesPerElem];
  double** values_2d = new double*[numNodesPerElem];

  for(i=0; i<numNodesPerElem*numNodesPerElem; ++i) values_1d[i] = 1.0;

  int offset = 0;
  for(i=0; i<numNodesPerElem; ++i) {
    values_2d[i] = &(values_1d[offset]);
    offset += numNodesPerElem;
  }

  for(i=0; i<numElems; ++i) {
    switch(i) {
    case 0:
      nodes[0] = 0; nodes[1] = 1; nodes[2] = 4; nodes[3] = 3;
      break;

    case 1:
      nodes[0] = 1; nodes[1] = 2; nodes[2] = 5; nodes[3] = 4;
      break;

    case 2:
      nodes[0] = 3; nodes[1] = 4; nodes[2] = 7; nodes[3] = 6;
      break;

     case 3:
      nodes[0] = 4; nodes[1] = 5; nodes[2] = 8; nodes[3] = 7;
      break;
    }

    for(j=0; j<numNodesPerElem; ++j) {
      if (preconstruct_graph) {
	err = A->BeginSumIntoGlobalValues(nodes[j], numNodesPerElem, nodes);
	if (err<0) return(err);
      }
      else {
	err = A->BeginInsertGlobalValues(nodes[j], numNodesPerElem, nodes);
	if (err<0) return(err);
      }
    
      for(k=0; k<numNodesPerElem; ++k) {
	err = A->SubmitBlockEntry(values_1d, blockSize, blockSize, blockSize);
	if (err<0) return(err);
      }

      err = A->EndSubmitEntries();
      if (err<0) return(err);
    }
  }

  EPETRA_CHK_ERR( A->GlobalAssemble() );

  Epetra_FEVbrMatrix* Acopy = new Epetra_FEVbrMatrix(*A);

  if (verbose) {
    cout << "A:"<<*A << endl;
    cout << "Acopy:"<<*Acopy<<endl;
  }

  Epetra_Vector x(A->RowMap()), y(A->RowMap());

  x.PutScalar(1.0); y.PutScalar(0.0);

  Epetra_Vector x2(Acopy->RowMap()), y2(Acopy->RowMap());

  x2.PutScalar(1.0); y2.PutScalar(0.0);

  A->Multiply(false, x, y);

  Acopy->Multiply(false, x2, y2);

  double ynorm2, y2norm2;

  y.Norm2(&ynorm2);
  y2.Norm2(&y2norm2);
  if (ynorm2 != y2norm2) {
    cerr << "norm2(A*ones) != norm2(*Acopy*ones)"<<endl;
    return(-99);
  }

  Epetra_FEVbrMatrix* Acopy2 =
    new Epetra_FEVbrMatrix(Copy, A->RowMap(), A->ColMap(), 1);

  *Acopy2 = *Acopy;

  Epetra_Vector x3(Acopy->RowMap()), y3(Acopy->RowMap());

  x3.PutScalar(1.0); y3.PutScalar(0.0);

  Acopy2->Multiply(false, x3, y3);

  double y3norm2;
  y3.Norm2(&y3norm2);

  if (y3norm2 != y2norm2) {
    cerr << "norm2(Acopy*ones) != norm2(Acopy2*ones)"<<endl;
    return(-999);
  }

  int len = 20;
  int* indices = new int[len];
  double* values = new double[len];
  int numIndices;

  if (map.MyGID(0)) {
    int lid = map.LID(0);
    EPETRA_CHK_ERR( A->ExtractMyRowCopy(lid, len, numIndices,
					values, indices) );
    if (numIndices != 4) {
      return(-1);
    }
    if (indices[0] != lid) {
      return(-2);
    }

    if (values[0] != 1.0*numProcs) {
      cout << "ERROR: values[0] ("<<values[0]<<") should be "<<numProcs<<endl;
      return(-3);
    }
  }

  if (map.MyGID(4)) {
    int lid = map.LID(4);
    EPETRA_CHK_ERR( A->ExtractMyRowCopy(lid, len, numIndices,
					values, indices) );

    if (numIndices != 9) {
      return(-4);
    }
    int lcid = A->LCID(4);
//     if (indices[lcid] != 4) {
//       cout << "ERROR: indices[4] ("<<indices[4]<<") should be "
// 	   <<A->LCID(4)<<endl;
//       return(-5);
//     }
    if (values[lcid] != 4.0*numProcs) {
      cout << "ERROR: values["<<lcid<<"] ("<<values[lcid]<<") should be "
	   <<4*numProcs<<endl;
      return(-6);
    }
  }

  delete [] values_2d;
  delete [] values_1d;
  delete [] nodes;
  delete [] indices;
  delete [] values;

  delete A;
  delete Acopy2;
  delete Acopy;
  delete graph;

  return(0);
}
