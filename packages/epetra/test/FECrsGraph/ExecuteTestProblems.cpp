//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include "Epetra_BLAS.h"
#include "ExecuteTestProblems.h"
#include "Epetra_Comm.h"
#include "Epetra_Vector.h"
#include "Epetra_FEVector.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
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
  //Each processor will pass a 3x3 element-connectivity-matrix to
  //Epetra_FECrsGraph.
  //After GlobalAssemble(), the graph should be as follows:
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

  int numMyNodes = 2;
  int* myNodes = new int[numMyNodes];

  if (localProc == 0) {
    myNodes[0] = 0;
    myNodes[1] = 1;
  }
  else {
    myNodes[0] = 2;
    myNodes[1] = 3;
  }

  Epetra_Map Map(-1, numMyNodes, myNodes, indexBase, map.Comm());

  delete [] myNodes;
  numMyNodes = 3;
  myNodes = new int[numMyNodes];

  if (localProc == 0) {
    myNodes[0] = 0;
    myNodes[1] = 1;
    myNodes[2] = 3;
  }
  else {
    myNodes[0] = 1;
    myNodes[1] = 2;
    myNodes[2] = 3;
  }

  int rowLengths = 3;
  Epetra_FECrsGraph A(Copy, Map, rowLengths);

  EPETRA_TEST_ERR( A.InsertGlobalIndices(numMyNodes, myNodes,
					 numMyNodes, myNodes),ierr);

  EPETRA_TEST_ERR( A.GlobalAssemble(), ierr );

  if (verbose) {
    A.Print(cout);
  }

  delete [] myNodes;

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
  //Each processor will pass a 3x3 element-connectivity-matrix to
  //Epetra_FECrsGraph.
  //After GlobalAssemble(), the graph should be as follows:
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
  int numMyNodes = 3;
  int* myNodes = new int[numMyNodes];

  if (localProc == 0) {
    myNodes[0] = 0;
    myNodes[1] = 1;
    myNodes[2] = 3;
  }
  else {
    numMyNodes = 1;
    myNodes[0] = 2;
  }

  Epetra_Map Map(-1, numMyNodes, myNodes, indexBase, map.Comm());

  int rowLengths = 3;
  Epetra_FECrsGraph A(Copy, Map, rowLengths);

  if (localProc != 0) {
    numMyNodes = 3;
    myNodes[0] = 1;
    myNodes[1] = 2;
    myNodes[2] = 3;
  }

  EPETRA_TEST_ERR( A.InsertGlobalIndices(numMyNodes, myNodes,
					 numMyNodes, myNodes),ierr);

  EPETRA_TEST_ERR( A.GlobalAssemble(), ierr );

  if (verbose) {
    A.Print(cout);
  }

  delete [] myNodes;

  return(0);
}

int four_quads(const Epetra_Comm& Comm, bool preconstruct_graph, bool verbose)
{
  if (verbose) {
    cout << "******************* four_quads ***********************"<<endl;
  }

  //This function assembles a matrix representing a finite-element
  //mesh of four 2-D quad elements. There are 9 nodes in the problem. The
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

  Epetra_FECrsGraph* graph = NULL;

  int* nodes = new int[numNodesPerElem];
  int i, j, err = 0;

  if (preconstruct_graph) {
    graph = new Epetra_FECrsGraph(Copy, map, 1);

    //we're going to fill the graph with indices, by passing our
    //connectivity lists.
    //FECrsGraph should accept indices in all rows, regardless of
    //whether map.MyGID(row) is true.

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

      err = graph->InsertGlobalIndices(numNodesPerElem, nodes,
                                       numNodesPerElem, nodes);
      if (err < 0) {
        cerr << "ERROR, FECrsGraph error in InsertGlobalIndices, err="
          << err << endl;
        return(-1);
      }
    }

    EPETRA_CHK_ERR( graph->GlobalAssemble() );
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

  err = A->GlobalAssemble();
  if (err < 0) {
    return(err);
  }

  Epetra_Vector x(A->RowMap()), y(A->RowMap());

  x.PutScalar(1.0); y.PutScalar(0.0);

  Epetra_FECrsMatrix Acopy(*A);

  Epetra_Vector x2(Acopy.RowMap()), y2(Acopy.RowMap());

  x2.PutScalar(1.0); y2.PutScalar(0.0);

  A->Multiply(false, x, y);

  Acopy.Multiply(false, x2, y2);

  double ynorm2, y2norm2;

  y.Norm2(&ynorm2);
  y2.Norm2(&y2norm2);
  if (ynorm2 != y2norm2) {
    cerr << "norm2(A*ones) != norm2(Acopy*ones)"<<endl;
    return(-99);
  }

  err = Acopy.GlobalAssemble();
  if (err < 0) {
    return(err);
  }

  if (verbose) {
    cout << "A:"<<endl<<*A << endl;
    cout << "Acopy:"<<endl<<Acopy << endl;
  }

  Epetra_FECrsMatrix Acopy2(Copy, A->RowMap(), A->ColMap(), 1);

  Acopy2 = Acopy;

  Epetra_Vector x3(Acopy.RowMap()), y3(Acopy.RowMap());

  x3.PutScalar(1.0); y3.PutScalar(0.0);

  Acopy2.Multiply(false, x3, y3);

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

// now let's do the checks for Acopy...

  if (map.MyGID(0)) {
    EPETRA_CHK_ERR( Acopy.ExtractGlobalRowCopy(0, len, numIndices,
                                            values, indices) );
    if (numIndices != 4) {
      return(-1);
    }
    if (indices[0] != 0) {
      return(-2);
    }

    if (values[0] != 1.0*numProcs) {
      cout << "ERROR: Acopy.values[0] ("<<values[0]<<") should be "<<numProcs<<endl;
      return(-3);
    }
  }

  if (map.MyGID(4)) {
    EPETRA_CHK_ERR( Acopy.ExtractGlobalRowCopy(4, len, numIndices,
                                            values, indices) );

    if (numIndices != 9) {
      return(-4);
    }
    int lcid = A->LCID(4);
    if (lcid<0) {
      return(-5);
    }
    if (values[lcid] != 4.0*numProcs) {
      cout << "ERROR: Acopy.values["<<lcid<<"] ("<<values[lcid]<<") should be "
           <<4*numProcs<<endl;
      return(-6);
    }
  }

// now let's do the checks for Acopy2...

  if (map.MyGID(0)) {
    EPETRA_CHK_ERR( Acopy2.ExtractGlobalRowCopy(0, len, numIndices,
                                            values, indices) );
    if (numIndices != 4) {
      return(-1);
    }
    if (indices[0] != 0) {
      return(-2);
    }

    if (values[0] != 1.0*numProcs) {
      cout << "ERROR: Acopy2.values[0] ("<<values[0]<<") should be "<<numProcs<<endl;
      return(-3);
    }
  }

  if (map.MyGID(4)) {
    EPETRA_CHK_ERR( Acopy2.ExtractGlobalRowCopy(4, len, numIndices,
                                            values, indices) );

    if (numIndices != 9) {
      return(-4);
    }
    int lcid = A->LCID(4);
    if (lcid<0) {
      return(-5);
    }
    if (values[lcid] != 4.0*numProcs) {
      cout << "ERROR: Acopy2.values["<<lcid<<"] ("<<values[lcid]<<") should be "
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

