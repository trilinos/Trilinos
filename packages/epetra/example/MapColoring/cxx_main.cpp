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

#include <assert.h>
#include <stdio.h>
#include <set>
#include <algorithm>
#include <vector>
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MapColoring.h"
#include "Epetra_CrsMatrix.h"
#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

using namespace std;

int main(int argc, char * argv[]) {

// Set up the Epetra communicator
#ifdef EPETRA_MPI
  MPI_Init(&argc, &argv);
  int size;      // Number of MPI processes
  int rank;      // My process ID
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  int size = 1;  // Serial case (not using MPI)
  int rank = 0;  // Processor ID
  Epetra_SerialComm Comm;
#endif
  //  cout << Comm << endl;

  int MyPID    = Comm.MyPID();
  int NumProc  = Comm.NumProc();
  bool verbose = (MyPID == 0);
  cout << "MyPID   = " << MyPID   << endl
       << "NumProc = " << NumProc << endl;

  // Get the problem size from the command line argument
  if (argc < 2 || argc > 3) {
    if (verbose)
      cout << "Usage: " << argv[0] << " nx [ny]" << endl;
    exit(1);
  } // end if
  int nx = atoi(argv[1]);            // Get the dimensions for a 1D or 2D
  int ny = 1;                        // central difference problem
  if (argc == 3) ny = atoi(argv[2]);
  int NumGlobalElements = nx * ny;
  if (NumGlobalElements < NumProc) {
    if (verbose)
      cout << "numGlobalElements = " << NumGlobalElements 
	   << " cannot be < number of processors = " << NumProc << endl;
    exit(1);
  } // end if

  // Epetra distribution map
  int IndexBase = 0;       // Zero-based indices
  Epetra_Map Map(NumGlobalElements, IndexBase, Comm);
  //  if (verbose) cout << Map << endl;

  // Extract the global indices of the elements local to this processor
  int NumMyElements = Map.NumMyElements();
  int * MyGlobalElements = new int[NumMyElements];
  Map.MyGlobalElements(MyGlobalElements);
  for (int p = 0; p < NumProc; p++) 
    if (p == MyPID) {
      cout << endl << "Processor " << MyPID << ": Global elements = ";
      for (int i = 0; i < NumMyElements; i++)
	cout << MyGlobalElements[i] << " ";
      cout << endl;
      Comm.Barrier();
    } // end if

  // Create the number of non-zeros for a tridiagonal (1D problem) or banded
  // (2D problem) matrix
  int * NumNz = new int[NumMyElements];
  int global_i;
  int global_j;
  for (int i = 0; i < NumMyElements; i++) {
    NumNz[i] = 5;
    global_j = MyGlobalElements[i] / nx;
    global_i = MyGlobalElements[i] - global_j * nx;
    if (global_i == 0)    NumNz[i] -= 1;  // By having separate statements,
    if (global_i == nx-1) NumNz[i] -= 1;  // this works for 2D as well as 1D
    if (global_j == 0)    NumNz[i] -= 1;  // systems (i.e. nx x 1 or 1 x ny)
    if (global_j == ny-1) NumNz[i] -= 1;  // or even a 1 x 1 system
  }
  if (verbose) {
    cout << endl << "NumNz: ";
    for (int i = 0; i < NumMyElements; i++) cout << NumNz[i] << " ";
    cout << endl;
  } // end if

  // Create the Epetra Compressed Row Sparse matrix
  // Note: the actual values for the matrix entries are meaningless for
  // this exercise, but I assign them anyway.
  Epetra_CrsMatrix A(Copy, Map, NumNz);

  double * Values = new double[4];
  for (int i = 0; i < 4; i++) Values[i] = -1.0;
  int * Indices = new int[4];
  double diag = 2.0;
  if (ny > 1) diag = 4.0;
  int NumEntries;

  for (int i = 0; i < NumMyElements; i++) {
    global_j = MyGlobalElements[i] / nx;
    global_i = MyGlobalElements[i] - global_j * nx;
    NumEntries = 0;
    if (global_j > 0 && ny > 1)
      Indices[NumEntries++] = global_i   + (global_j-1)*nx;
    if (global_i > 0)
      Indices[NumEntries++] = global_i-1 +  global_j   *nx;
    if (global_i < nx-1)
      Indices[NumEntries++] = global_i+1 +  global_j   *nx;
    if (global_j < ny-1 && ny > 1)
      Indices[NumEntries++] = global_i   + (global_j+1)*nx;

    // Put in the off-diagonal terms
    assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values,
				Indices) == 0);
    // Put in the diagonal entry
    assert(A.InsertGlobalValues(MyGlobalElements[i], 1, &diag,
				MyGlobalElements+i) == 0);
  } // end i loop

  // Finish up matrix construction
  delete [] Values;
  delete [] Indices;
  assert(A.FillComplete() == 0);
  // cout << endl << A << endl;

  // Create the local distance-1 adjancency graph
  // This is essentially a transpose of the Epetra_CrsGraph, where off-
  // processor couplings are ignored and global indexes are converted to
  // local.  We use the C++ standard libraries vector and set, since we
  // don't know how many nonzeroes we will end up with for each column.

  vector< set<int> > adj1(NumMyElements);
  for (int lr = 0; lr < adj1.size(); lr++) {
    int lrid;   // Local row ID
    double * Values = new double[NumNz[lr]];
    int * Indices = new int[NumNz[lr]];
    assert(A.ExtractMyRowCopy(lr, NumNz[lr], NumNz[lr], Values, Indices) == 0);
    for (int i = 0; i < NumNz[lr]; i++) {
      lrid = A.LRID(Indices[i]);
      if (lrid >= 0) adj1[lrid].insert(lr);
    } // end i loop
    delete [] Values;
    delete [] Indices;
  } // end lr loop

  if (verbose) {
    cout << endl;
    for (int lr = 0; lr < NumMyElements; lr++) {
      cout << "adj1[" << lr << "] = { ";
      for (set<int>::const_iterator p = adj1[lr].begin(); p != adj1[lr].end();
	   p++) cout << *p << " ";
      cout << "}" << endl;
    } // end lr loop
  } // end if

  // Create the local distance-2 adjancency graph
  // This is built from the distance-1 adjancency graph.  We use the C++
  // standard libraries vector and set, since we don't know how many
  // nonzeroes we will end up with for each column.

  vector< set<int> > adj2(NumMyElements);
  for (int lc = 0; lc < NumMyElements; lc++) {
    for (set<int>::const_iterator p = adj1[lc].begin(); p != adj1[lc].end();
	 p++) {
      int lrid;    // Local row ID
      double * Values = new double[NumNz[*p]];
      int * Indices = new int[NumNz[*p]];
      assert(A.ExtractMyRowCopy(*p, NumNz[*p], NumNz[*p], Values, Indices)
	     == 0);
      for (int i = 0; i < NumNz[*p]; i++) {
	lrid = A.LRID(Indices[i]);
	if (lrid >= 0) adj2[lc].insert(lrid);
      } // end i loop
      delete [] Values;
      delete [] Indices;
    } // end p loop
  } // end lc loop

  cout << endl;
  for (int lc = 0; lc < NumMyElements; lc++) {
    cout << "adj2[" << lc << "] = { ";
    for (set<int>::const_iterator p = adj2[lc].begin(); p != adj2[lc].end();
	 p++) cout << *p << " ";
    cout << "}" << endl;
  } // end lc loop

  // Now that we have the local distance-2 adjacency graph, we can compute a
  // color map using a greedy algorithm.  The first step is to compute Delta,
  // the maximum size (degree) of adj1.
  size_t Delta = 0;
  for (int i = 0; i < NumMyElements; i++)
    Delta = max(Delta, adj1[i].size());
  cout << endl << "Delta = " << Delta << endl << endl;
  
  // Now create a color map and initialize all values to 0, which
  // indicates that none of the columns have yet been colored.
  int * color_map = new int[NumMyElements];
  for (int i = 0; i < NumMyElements; i++) color_map[i] = 0;

  // Apply the distance-2 greedy coloring algorithm
  for (int column = 0; column < NumMyElements; column++) {
    set<int> allowedColors;    // Create the set of allowed colors
    for (int i = 1; i < Delta*Delta+1; i++) allowedColors.insert(i);
    for (set<int>::const_iterator p = adj2[column].begin();
	 p != adj2[column].end(); p++) if (color_map[*p] > 0)
	   allowedColors.erase(color_map[*p]);
    color_map[column] = *(allowedColors.begin());
    cout << "color_map[" << column << "] = " << color_map[column] << endl;
  } // end col loop

  // New section to Epetra_MapColoring
  Epetra_MapColoring C1(Map, color_map);

  cout << C1;
  
  // Clean up

  delete [] MyGlobalElements;
  delete [] NumNz;
  delete [] color_map;

  cout << endl << argv[0] << " done." << endl;

} // end main

