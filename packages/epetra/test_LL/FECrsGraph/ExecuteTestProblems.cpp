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
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
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

  Epetra_Map Map((long long) -1, numMyNodes, myNodes, indexBase, map.Comm());

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
    A.Print(std::cout);
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

  Epetra_Map Map((long long) -1, numMyNodes, myNodes, indexBase, map.Comm());

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
    A.Print(std::cout);
  }

  delete [] myNodes;

  return(0);
}

int four_quads(const Epetra_Comm& Comm, bool preconstruct_graph, bool verbose)
{
  if (verbose) {
    std::cout << "******************* four_quads ***********************"<<std::endl;
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

  int numProcs = Comm.NumProc();

  long long numNodes = 9;
  int numElems = 4;
  int numNodesPerElem = 4;

  int indexBase = 0;

  //Create a map using epetra-defined linear distribution.
  Epetra_Map map(numNodes, indexBase, Comm);

  Epetra_FECrsGraph* graph = NULL;

  long long* nodes = new long long[numNodesPerElem];
  int i, err = 0;

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
        std::cerr << "ERROR, FECrsGraph error in InsertGlobalIndices, err="
          << err << std::endl;
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
  Epetra_LongLongSerialDenseVector epetra_nodes(View, nodes, numNodesPerElem);
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
    std::cerr << "norm2(A*ones) != norm2(Acopy*ones)"<<std::endl;
    return(-99);
  }

  err = Acopy.GlobalAssemble();
  if (err < 0) {
    return(err);
  }

  if (verbose) {
    std::cout << "A:"<<std::endl<<*A << std::endl;
    std::cout << "Acopy:"<<std::endl<<Acopy << std::endl;
  }

  Epetra_FECrsMatrix Acopy2(Copy, A->RowMap(), A->ColMap(), 1);

  Acopy2 = Acopy;

  Epetra_Vector x3(Acopy.RowMap()), y3(Acopy.RowMap());

  x3.PutScalar(1.0); y3.PutScalar(0.0);

  Acopy2.Multiply(false, x3, y3);

  double y3norm2;
  y3.Norm2(&y3norm2);

  if (y3norm2 != y2norm2) {
    std::cerr << "norm2(Acopy*ones) != norm2(Acopy2*ones)"<<std::endl;
    return(-999);
  }

  int len = 20;
  long long* indices = new long long[len];
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
      std::cout << "ERROR: values[0] ("<<values[0]<<") should be "<<numProcs<<std::endl;
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
      std::cout << "ERROR: values["<<lcid<<"] ("<<values[lcid]<<") should be "
	   <<4*numProcs<<std::endl;
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
      std::cout << "ERROR: Acopy.values[0] ("<<values[0]<<") should be "<<numProcs<<std::endl;
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
      std::cout << "ERROR: Acopy.values["<<lcid<<"] ("<<values[lcid]<<") should be "
           <<4*numProcs<<std::endl;
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
      std::cout << "ERROR: Acopy2.values[0] ("<<values[0]<<") should be "<<numProcs<<std::endl;
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
      std::cout << "ERROR: Acopy2.values["<<lcid<<"] ("<<values[lcid]<<") should be "
           <<4*numProcs<<std::endl;
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

int Young1(const Epetra_Comm& Comm, bool verbose)
{
  //This is a test case submitted by Joe Young with bug 2421. It runs
  //only on 2 processors.
  if (Comm.NumProc() != 2) {
    return(0);
  }

  // Give rows 0-2 to proc 0 and 3-5 to proc 1
  int             RowIndices[3];
  if (Comm.MyPID() == 0) {
    RowIndices[0] = 0;
    RowIndices[1] = 1;
    RowIndices[2] = 2;
  } else {
    RowIndices[0] = 3;
    RowIndices[1] = 4;
    RowIndices[2] = 5;
  }
  Epetra_Map      RangeMap((long long) -1, 3, RowIndices, 0, Comm);
  Epetra_Map & RowMap = RangeMap;

  // Define a second map that gives col 0 to proc 0 and col 1 to proc 1 
  int             ColIndices[1];
  if (Comm.MyPID() == 0) {
    ColIndices[0] = 0;
  }
  else {
    ColIndices[0] = 1;
  }
  Epetra_Map      DomainMap((long long) -1, 1, ColIndices, 0, Comm);

  // Construct a graph where both processors only insert into local
  // elements
  Epetra_FECrsGraph BrokenGraph(Copy, RowMap, 2);
  for (int i = 0; i < RangeMap.NumMyElements(); i++) {
    int             ig = RowIndices[i];
    int             jgs[2] = { 0, 1 };
    BrokenGraph.InsertGlobalIndices(1, &ig, 2, jgs);
  }
  BrokenGraph.GlobalAssemble(DomainMap, RangeMap);

  // Check the size of the matrix that would be created from the graph 
  long long numCols1 = BrokenGraph.NumGlobalCols();
  if (verbose) {
    std::cout << "Number of global rows in the graph where only "
        "local elements were inserted: " << BrokenGraph.NumGlobalRows()
      << std::endl;
    std::cout << "Number of global cols in the graph where only "
        "local elements were inserted: " << BrokenGraph.NumGlobalCols()
      << std::endl;
  }
  // Construct a graph where both processors insert into global elements
  Epetra_FECrsGraph Graph(Copy, RowMap, 2);
  for (int i = 0; i < 6; i++) {
    int             ig = i;
    int             jgs[2] = { 0, 1 };
    Graph.InsertGlobalIndices(1, &ig, 2, jgs);
  }
  Graph.GlobalAssemble(DomainMap, RangeMap);

  // Check the size of the matrix that would be created from the graph 
  long long numCols2 = Graph.NumGlobalCols();
  if (verbose) {
    std::cout << "Number of global rows in the graph where "
        "global elements were inserted: " << Graph.NumGlobalRows()
       << std::endl;
    std::cout << "Number of global cols in the graph where "
        "global elements were inserted: " << Graph.NumGlobalCols()
      << std::endl;
  }

  if (numCols1 != numCols2) return(-1);
  return(0);
}

int rectangular(const Epetra_Comm& Comm, bool verbose)
{
  int mypid = Comm.MyPID();
  int numlocalrows = 3;
  Epetra_Map rowmap((long long) -1, numlocalrows, 0, Comm);

  int numglobalrows = numlocalrows*Comm.NumProc();

  int numcols = 2*numglobalrows;

  Epetra_FECrsGraph fegraph(Copy, rowmap, numcols);

  int* cols = new int[numcols];
  for(int j=0; j<numcols; ++j) cols[j] = j;

  Epetra_Map domainmap((long long) -1, numcols, 0, Comm);

  int firstlocalrow = numlocalrows*mypid;
  int lastlocalrow = numlocalrows*(mypid+1)-1;

  for(int i=0; i<numglobalrows; ++i) {
    //if i is a local row, then skip it. We want each processor to only
    //load rows that belong on other processors.
    if (i >= firstlocalrow && i <= lastlocalrow) continue;

    EPETRA_CHK_ERR( fegraph.InsertGlobalIndices(1, &i, numcols, &(cols[0])) );
  }

  EPETRA_CHK_ERR( fegraph.GlobalAssemble(domainmap, rowmap) );

  if (verbose) {
    std::cout << "********************** fegraph **********************" << std::endl;
    std::cout << fegraph << std::endl;
  }

  delete [] cols;

  return(0);
}

