//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER

#include "Epetra_BLAS.h"
#include "ExecuteTestProblems.h"
#include "Epetra_Comm.h"
#include "Epetra_Vector.h"
#include "Epetra_FEVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"


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

    EPETRA_TEST_ERR( A.InsertGlobalValues(numMyNodes, myNodes,
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

    EPETRA_TEST_ERR( A.InsertGlobalValues(numMyNodes, myNodes,
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

    EPETRA_TEST_ERR( A.InsertGlobalValues(numMyNodes, myNodes,
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

    EPETRA_TEST_ERR( A.InsertGlobalValues(numMyNodes, myNodes,
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

  A.InsertGlobalValues(NumCols, Indices,
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

  int myPID = Comm.MyPID();
  int numProcs = Comm.NumProc();

  int numNodes = 9;
  int numElems = 4;
  int numNodesPerElem = 4;

  int indexBase = 0;

  //Create a map using epetra-defined linear distribution.
  Epetra_Map map(numNodes, indexBase, Comm);

  Epetra_CrsGraph* graph = NULL;

  int* nodes = new int[numNodesPerElem];
  int i, j, err = 0;

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

  Epetra_FECrsMatrix* A = NULL;

  if (preconstruct_graph) {
    A = new Epetra_FECrsMatrix(Copy, *graph);
  }
  else {
    A = new Epetra_FECrsMatrix(Copy, map, 1);
  }

  EPETRA_CHK_ERR( A->PutScalar(0.0) );

  double* values_1d = new double[numNodesPerElem*numNodesPerElem];
  double** values_2d = new double*[numNodesPerElem];

  for(i=0; i<numNodesPerElem*numNodesPerElem; ++i) values_1d[i] = 1.0;

  int offset = 0;
  for(i=0; i<numNodesPerElem; ++i) {
    values_2d[i] = &(values_1d[offset]);
    offset += numNodesPerElem;
  }

  int format = Epetra_FECrsMatrix::ROW_MAJOR;
  Epetra_IntSerialDenseVector epetra_nodes(View, nodes, numNodesPerElem);
  Epetra_SerialDenseMatrix epetra_values(View, values_1d, numNodesPerElem,
					 numNodesPerElem, numNodesPerElem);

  for(i=0; i<numElems; ++i) {
    switch(i) {
    case 0:
      nodes[0] = 0; nodes[1] = 1; nodes[2] = 4; nodes[3] = 3;
      if (preconstruct_graph) {
	err = A->SumIntoGlobalValues(epetra_nodes,
				     epetra_values, format);
	if (err<0) return(err);
      }
      else {
	err = A->InsertGlobalValues(epetra_nodes,
				    epetra_values, format);
	if (err<0) return(err);
      }
      break;

    case 1:
      nodes[0] = 1; nodes[1] = 2; nodes[2] = 5; nodes[3] = 4;
      if (preconstruct_graph) {
	err = A->SumIntoGlobalValues(numNodesPerElem, nodes,
				     values_2d, format);
	if (err<0) return(err);
      }
      else {
	err = A->InsertGlobalValues(numNodesPerElem, nodes,
				    values_2d, format);
	if (err<0) return(err);
      }
      break;

    case 2:
      nodes[0] = 3; nodes[1] = 4; nodes[2] = 7; nodes[3] = 6;
      if (preconstruct_graph) {
	err = A->SumIntoGlobalValues(numNodesPerElem, nodes,
				     numNodesPerElem, nodes,
				     values_1d, format);
	if (err<0) return(err);
      }
      else {
	err = A->InsertGlobalValues(numNodesPerElem, nodes,
				    numNodesPerElem, nodes,
				    values_1d, format);
	if (err<0) return(err);
      }
      break;

     case 3:
      nodes[0] = 4; nodes[1] = 5; nodes[2] = 8; nodes[3] = 7;
      if (preconstruct_graph) {
	err = A->SumIntoGlobalValues(numNodesPerElem, nodes,
				     numNodesPerElem, nodes,
				     values_2d, format);
	if (err<0) return(err);
      }
      else {
	err = A->InsertGlobalValues(numNodesPerElem, nodes,
				    numNodesPerElem, nodes,
				    values_2d, format);
	if (err<0) return(err);
      }
      break;
    }
  }

  EPETRA_CHK_ERR( A->GlobalAssemble() );

  if (verbose) {
    cout << *A << endl;
  }

  int len = 20;
  int* indices = new int[len];
  double* values = new double[len];
  int numIndices;

  if (map.MyGID(0)) {
    EPETRA_CHK_ERR( A->ExtractGlobalRowCopy(0, len, numIndices,
					    values, indices) );
    if (numIndices != 4) {
      return(-1);
    }
    if (indices[0] != 0) {
      return(-2);
    }

    if (values[0] != 1.0*numProcs) {
      cout << "ERROR: values[0] ("<<values[0]<<") should be "<<numProcs<<endl;
      return(-3);
    }
  }

  if (map.MyGID(4)) {
    EPETRA_CHK_ERR( A->ExtractGlobalRowCopy(4, len, numIndices,
					    values, indices) );

    if (numIndices != 9) {
      return(-4);
    }
    int lcid = A->LCID(4);
    if (lcid<0) {
      return(-5);
    }
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
  delete graph;

  return(0);
}
