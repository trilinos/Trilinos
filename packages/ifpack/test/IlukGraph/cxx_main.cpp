//@HEADER
// ***********************************************************************
// 
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
// ***********************************************************************
//@HEADER

#include "Ifpack_ConfigDefs.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "Epetra_Map.h"
#include "Ifpack_Version.h"
#include "Epetra_CrsGraph.h"
#include "Ifpack_IlukGraph.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

// Prototype

int check(Epetra_CrsGraph& L, Epetra_CrsGraph& U, Ifpack_IlukGraph& LU, 
	  int NumGlobalRows1, int NumMyRows1, int LevelFill1, bool verbose);

 int main(int argc, char *argv[])
{
  int ierr = 0, i, j;
  int NumIndices;
  int * Indices;
  int nx, ny;

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );

#else

  int size = 1; // Serial case (not using MPI)
  int rank = 0;

  Epetra_SerialComm Comm;
#endif

  bool verbose = false;

  int nextarg = 1;
  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') {
    verbose = true;
    nextarg++;
  }

  // char tmp;
  // if (rank==0) cout << "Press any key to continue..."<< endl;
  // if (rank==0) cin >> tmp;
  // Comm.Barrier();

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  if (verbose && MyPID==0)
    cout << Ifpack_Version() << endl << endl;

  if (verbose) cout << Comm <<endl;

  int sqrtNumProc = (int) ceil(sqrt((double) NumProc));

  bool verbose1 = verbose;
  verbose = verbose && (MyPID==0);

  if (verbose1 && argc != 4) {
    nx = 10;
    ny = 12*NumProc;
    cout << "Setting nx = " << nx << ", ny = " << ny << endl;
  }
    else if (!verbose1 && argc != 3) {
    nx = 10;
    ny = 12*NumProc;
    }
  else {
    nx = atoi(argv[nextarg++]);
    if (nx<3) {cout << "nx = " << nx << ": Must be greater than 2 for meaningful graph." << endl; exit(1);}
    ny = atoi(argv[nextarg++]);
    if (ny<3) {cout << "ny = " << ny << ": Must be greater than 2 for meaningful graph." << endl; exit(1);}
  }
  
  int NumGlobalPoints = nx*ny;
  int IndexBase = 0;

  if (verbose)
    cout << "\n\n*****Building 5 point matrix, Level 1 and 2 filled matrices for" << endl
	 << "  nx = " << nx << ",  ny = " << ny << endl<< endl;


  // Create a 5 point stencil graph, level 1 fill of it and level 2 fill of it

  Epetra_Map Map(NumGlobalPoints, IndexBase, Comm);

  int NumMyPoints = Map.NumMyPoints();

  Epetra_CrsGraph A(Copy, Map, 5);
  Epetra_CrsGraph L0(Copy, Map, Map, 2);
  Epetra_CrsGraph U0(Copy, Map, Map, 2);
  Epetra_CrsGraph L1(Copy, Map, Map, 3);
  Epetra_CrsGraph U1(Copy, Map, Map, 3);
  Epetra_CrsGraph L2(Copy, Map, Map, 4);
  Epetra_CrsGraph U2(Copy, Map, Map, 4);
  
  // Add  rows one-at-a-time

  Indices = new int[4]; // Work space
  
  for (j=0; j<ny; j++) {
    for (i=0; i<nx; i++) {
      int Row = i+j*nx;
      if (Map.MyGID(Row)) { // Only work on rows I own
	
	//**** Work on lower triangle of all three matrices ****
	
	// Define entries (i-1,j), (i,j-1)
	
	int k = 0;
	if (i>0)    Indices[k++] = i-1 + j   *nx;
	if (j>0)    Indices[k++] = i   +(j-1)*nx;
	
	// Define lower triangular terms of original matrix and L(0)
	assert(A.InsertGlobalIndices(Row, k, Indices)>=0);
	assert(L0.InsertGlobalIndices(Row, k, Indices)>=0);
	
	// Define entry (i+1,j-1)
	if ((i<nx-1) && (j>0   )) Indices[k++] = i+1 +(j-1)*nx;
	
	
	// Define lower triangle of level(1) fill matrix
	assert(L1.InsertGlobalIndices(Row, k, Indices)>=0);
	
	// Define entry (i+2, j-1)
	
	if ((i<nx-2) && (j>0   )) Indices[k++] = i+2 +(j-1)*nx;
	
	// Define lower triangle of level(2) fill matrix
	assert(L2.InsertGlobalIndices(Row, k, Indices)>=0);
	
	// Define main diagonal of original matrix
	assert(A.InsertGlobalIndices(Row, 1, &Row)>=0);
	
	k = 0; // Reset index counter
	
	//**** Work on upper triangle of all three matrices ****
	
	// Define entries (i+1,j), ( i,j+1)
	
	if (i<nx-1) Indices[k++] = i+1 + j   *nx;
	if (j<ny-1) Indices[k++] = i   +(j+1)*nx;
	
	// Define upper  triangular terms of original matrix and L(0)
	assert(A.InsertGlobalIndices(Row, k, Indices)>=0);
	assert(U0.InsertGlobalIndices(Row, k, Indices)>=0);
	
	// Define entry (i-1,j+1)
	
	if ((i>0   ) && (j<ny-1)) Indices[k++] = i-1 +(j+1)*nx;
	
	// Define upper triangle of level(1) fill matrix
	assert(U1.InsertGlobalIndices(Row, k, Indices)>=0);
	
	// Define entry (i-2, j+1)
	
	if ((i>1   ) && (j<ny-1)) Indices[k++] = i-2 +(j+1)*nx;
	
	// Define upper triangle of level(2) fill matrix
	assert(U2.InsertGlobalIndices(Row, k, Indices)>=0);
      }
    }
  }

  delete [] Indices;

  // Finish up
  if (verbose) cout << "\n\nCompleting A" << endl<< endl;
  assert(A.TransformToLocal()==0);
  if (verbose) cout << "\n\nCompleting L0" << endl<< endl;
  assert(L0.TransformToLocal()==0);
  if (verbose) cout << "\n\nCompleting U0" << endl<< endl;
  assert(U0.TransformToLocal()==0);
  if (verbose) cout << "\n\nCompleting L1" << endl<< endl;
  assert(L1.TransformToLocal()==0);
  if (verbose) cout << "\n\nCompleting U1" << endl<< endl;
  assert(U1.TransformToLocal()==0);
  if (verbose) cout << "\n\nCompleting L2" << endl<< endl;
  assert(L2.TransformToLocal()==0);
  if (verbose) cout << "\n\nCompleting U2" << endl<< endl;
  assert(U2.TransformToLocal()==0);

  if (verbose) cout << "\n\n*****Testing ILU(0) constructor on A" << endl<< endl;

  Ifpack_IlukGraph ILU0(A, 0, 0);
  assert(ILU0.ConstructFilledGraph()==0);

  assert(check(L0, U0, ILU0, NumGlobalPoints, NumMyPoints, 0, verbose)==0);

  if (verbose) cout << "\n\n*****Testing ILU(1) constructor on A" << endl<< endl;

  Ifpack_IlukGraph ILU1(A, 1, 0);
  assert(ILU1.ConstructFilledGraph()==0);

  assert(check(L1, U1, ILU1, NumGlobalPoints, NumMyPoints, 1, verbose)==0);

  if (verbose) cout << "\n\n*****Testing ILU(2) constructor on A" << endl<< endl;

  Ifpack_IlukGraph ILU2(A, 2, 0);
  assert(ILU2.ConstructFilledGraph()==0);

  assert(check(L2, U2, ILU2, NumGlobalPoints, NumMyPoints, 2, verbose)==0);

  if (verbose) cout << "\n\n*****Testing copy constructor" << endl<< endl;

  Ifpack_IlukGraph ILUC(ILU2);
  
  assert(check(L2, U2, ILUC, NumGlobalPoints, NumMyPoints, 2, verbose)==0);

  if (verbose) cout << "\n\n*****Testing copy constructor" << endl<< endl;

  Ifpack_IlukGraph * OverlapGraph;
  for (int overlap = 1; overlap < 4; overlap++) {
    if (verbose) cout << "\n\n*********************************************" << endl;
    if (verbose) cout << "\n\nConstruct Level 1 fill with Overlap = " << overlap << ".\n\n" << endl;
    
    OverlapGraph = new Ifpack_IlukGraph(A, 1, overlap);
    assert(OverlapGraph->ConstructFilledGraph()==0);
    
    
    if (verbose) {
      cout << "Number of Global Rows     = " << OverlapGraph->NumGlobalRows() << endl;
      cout << "Number of Global Nonzeros = " << OverlapGraph->NumGlobalNonzeros() << endl;
      cout << "Number of Local Rows     = " << OverlapGraph->NumMyRows() << endl;
      cout << "Number of Local Nonzeros = " << OverlapGraph->NumMyNonzeros() << endl;
    }
    delete OverlapGraph;
  }
    
  if (verbose1) {
    // Test ostream << operator (if verbose1)
    // Construct a Map that puts 6 equations on each PE
    
    int NumElements1 = 6;
    int NumPoints1 = NumElements1;

    // Create an integer vector NumNz that is used to build the Petra Matrix.
    // NumNz[i] is the Number of terms for the ith global equation on this processor
    
    int * NumNz1 = new int[NumPoints1];
    
    // We are building a tridiagonal matrix where each row has (-1 2 -1)
    // So we need 2 off-diagonal terms (except for the first and last equation)
    
    for (i=0; i<NumPoints1; i++)
      if (i==0 || i == NumPoints1-1)
	NumNz1[i] = 2;
      else
	NumNz1[i] = 3;
    
    // Create a Epetra_Matrix
    
    Epetra_Map Map1(NumPoints1, NumPoints1, 1, Comm);
    Epetra_CrsGraph A1(Copy, Map1, NumNz1);
    
    // Add  rows one-at-a-time
    // Need some vectors to help
    // Off diagonal Values will always be -1
    
    
    int *Indices1 = new int[2];
    int NumEntries1;
    
    for (i=0; i<NumPoints1; i++)
      {
	if (i==0)
	  {
	    Indices1[0] = 2;
	    NumEntries1 = 1;
	  }
	else if (i == NumPoints1-1)
	  {
	    Indices1[0] = NumPoints1-1;
	    NumEntries1 = 1;
	  }
	else
	  {
	    Indices1[0] = i;
	    Indices1[1] = i+2;
	    NumEntries1 = 2;
	  }
	assert(A1.InsertGlobalIndices(i+1, NumEntries1, Indices1)==0);
	int ip1 = i+1;
	assert(A1.InsertGlobalIndices(ip1, 1, &ip1)==0); // Put in the diagonal entry
      }
    
    // Finish up
    assert(A1.TransformToLocal()==0);
    
    if (verbose) cout << "\n\nPrint out tridiagonal matrix with IndexBase = 1.\n\n" << endl;
    cout << A1 << endl;

    if (verbose) cout << "\n\nConstruct Level 1 fill with IndexBase = 1.\n\n" << endl;

    Ifpack_IlukGraph ILU11(A1, 1, 0);
    assert(ILU11.ConstructFilledGraph()==0);

    if (verbose) cout << "\n\nPrint out Level 1 ILU graph of tridiagonal matrix with IndexBase = 1.\n\n" << endl;
    if (verbose1) cout << ILU11 << endl;

    
  // Release all objects
  delete [] NumNz1;
  delete [] Indices1;

  }

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif


/* end main
*/
return ierr ;
}

int check(Epetra_CrsGraph& L, Epetra_CrsGraph& U, Ifpack_IlukGraph& LU, 
	  int NumGlobalRows1, int NumMyRows1, int LevelFill1, bool verbose) {  

  int i, j;
  int NumIndices, * Indices;
  int NumIndices1, * Indices1;

  bool debug = true;

  Epetra_CrsGraph& L1 = LU.L_Graph();
  Epetra_CrsGraph& U1 = LU.U_Graph();

  // Test entries and count nonzeros

  int Nout = 0;

  for (i=0; i<LU.NumMyRows(); i++) {

    assert(L.ExtractMyRowView(i, NumIndices, Indices)==0);
    assert(L1.ExtractMyRowView(i, NumIndices1, Indices1)==0);
    assert(NumIndices==NumIndices1);
    for (j=0; j<NumIndices1; j++) {
      if (debug &&(Indices[j]!=Indices1[j])) {
	int MyPID = L.RowMap().Comm().MyPID();
	cout << "Proc " << MyPID
	     << " Local Row = " << i
	     << "  L.Indices["<< j <<"]  = " << Indices[j]
	     << " L1.Indices["<< j <<"] = " << Indices1[j] << endl;
      }
      assert(Indices[j]==Indices1[j]);
    }
    Nout += (NumIndices-NumIndices1);

    assert(U.ExtractMyRowView(i, NumIndices, Indices)==0);
    assert(U1.ExtractMyRowView(i, NumIndices1, Indices1)==0);
    assert(NumIndices==NumIndices1);
    for (j=0; j<NumIndices1; j++)  {
      if (debug &&(Indices[j]!=Indices1[j])) {
	int MyPID = L.RowMap().Comm().MyPID();
	cout << "Proc " << MyPID
	     << " Local Row = " << i
	     << "  U.Indices["<< j <<"]  = " << Indices[j]
	     << " U1.Indices["<< j <<"] = " << Indices1[j] << endl;
      }
      assert(Indices[j]==Indices1[j]);
    }
    Nout += (NumIndices-NumIndices1);
  }

  // Test query functions

  int NumGlobalRows = LU.NumGlobalRows();
  if (verbose) cout << "\n\nNumber of Global Rows = " << NumGlobalRows << endl<< endl;

  assert(NumGlobalRows==NumGlobalRows1);

  int NumGlobalNonzeros = LU.NumGlobalNonzeros();
  if (verbose) cout << "\n\nNumber of Global Nonzero entries = " 
		    << NumGlobalNonzeros << endl<< endl;

  int NoutG = 0;

  L.RowMap().Comm().SumAll(&Nout, &NoutG, 1);

  assert(NumGlobalNonzeros==L.NumGlobalNonzeros()+U.NumGlobalNonzeros()-NoutG);

  int NumMyRows = LU.NumMyRows();
  if (verbose) cout << "\n\nNumber of Rows = " << NumMyRows << endl<< endl;

  assert(NumMyRows==NumMyRows1);

  int NumMyNonzeros = LU.NumMyNonzeros();
  if (verbose) cout << "\n\nNumber of Nonzero entries = " << NumMyNonzeros << endl<< endl;

  assert(NumMyNonzeros==L.NumMyNonzeros()+U.NumMyNonzeros()-Nout);

  if (verbose) cout << "\n\nLU check OK" << endl<< endl;

  return(0);
}
