//@HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
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
// ***********************************************************************
//@HEADER

// CrsMatrix_MapColoring Test routine
#include <Epetra_ConfigDefs.h>
#include "EpetraExt_Version.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#endif

#include "Epetra_SerialComm.h"
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_IntVector.h"
#include "Epetra_MapColoring.h"
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"
#include "../epetra_test_err.h"

#include <vector>

// prototype
void printColoring (const Epetra_MapColoring & ColorMap, Epetra_CrsGraph *A, bool verbose);

int main(int argc, char *argv[]) {


#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init( &argc, &argv );
  //int size, rank; // Number of MPI processes, My process ID

  //MPI_Comm_size(MPI_COMM_WORLD, &size);
  //MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#else

  //int size = 1; // Serial case (not using MPI)
  //int rank = 0;

#endif

  bool verbose = false;

  int nx = 5;
  int ny = 5;

  if( argc > 1 )
  {
    if( argc > 4 )
    {
      cout << "Usage: " << argv[0] << " [-v [nx [ny]]]" << endl;
      exit(1);
    }

    int loc = 1;
    // Check if we should print results to standard out
    if(argv[loc][0]=='-' && argv[loc][1]=='v')
    { verbose = true; ++loc; }

    if (loc < argc) nx = atoi( argv[loc++] );
    if( loc < argc) ny = atoi( argv[loc] );
  }

#ifdef EPETRA_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  bool verbose1 = false;
  if(verbose) verbose1 = (MyPID==0);

  if(verbose1)
    cout << EpetraExt::EpetraExt_Version() << endl << endl;

  Comm.Barrier();

  if(verbose) cout << Comm << endl << flush;
  Comm.Barrier();

  int NumGlobalElements = nx * ny;
  if( NumGlobalElements < NumProc )
  {
    cout << "NumGlobalElements = " << NumGlobalElements <<
            " cannot be < number of processors = " << NumProc;
    exit(1);
  } 
	
  int IndexBase = 0;
  Epetra_Map Map( NumGlobalElements, IndexBase, Comm );

  // Extract the global indices of the elements local to this processor
  int NumMyElements = Map.NumMyElements();
  std::vector<int> MyGlobalElements( NumMyElements );
  Map.MyGlobalElements( &MyGlobalElements[0] );
  if( verbose ) cout << Map;

  // Create the number of non-zeros for a tridiagonal (1D problem) or banded
  // (2D problem) matrix
  std::vector<int> NumNz( NumMyElements, 5 );
  int global_i;
  int global_j;
  for (int i = 0; i < NumMyElements; ++i)
  {
    global_j = MyGlobalElements[i] / nx;
    global_i = MyGlobalElements[i] - global_j * nx;
    if (global_i == 0)    NumNz[i] -= 1;  // By having separate statements,
    if (global_i == nx-1) NumNz[i] -= 1;  // this works for 2D as well as 1D
    if (global_j == 0)    NumNz[i] -= 1;  // systems (i.e. nx x 1 or 1 x ny)
    if (global_j == ny-1) NumNz[i] -= 1;  // or even a 1 x 1 system
  }
  if(verbose)
  { 
    cout << endl << "NumNz: ";
    for (int i = 0; i < NumMyElements; i++) cout << NumNz[i] << " ";
    cout << endl;
  } // end if
  
  // Create the Epetra Compressed Row Sparse Graph
  Epetra_CrsGraph A( Copy, Map, &NumNz[0] );
  
  std::vector<int> Indices(5);
  int NumEntries;
  
  for (int i = 0; i < NumMyElements; ++i )
  {
    global_j = MyGlobalElements[i] / nx;
    global_i = MyGlobalElements[i] - global_j * nx;
    NumEntries = 0;
    // (i,j-1) entry
    if (global_j > 0 && ny > 1)
      Indices[NumEntries++] = global_i   + (global_j-1)*nx;
    // (i-1,j) entry
    if (global_i > 0)
      Indices[NumEntries++] = global_i-1 +  global_j   *nx;
    // (i,j) entry
    Indices[NumEntries++] = MyGlobalElements[i];
    // (i+1,j) entry
    if (global_i < nx-1)
      Indices[NumEntries++] = global_i+1 +  global_j   *nx;
    // (i,j+1) entry
    if (global_j < ny-1 && ny > 1)
      Indices[NumEntries++] = global_i   + (global_j+1)*nx;

    // Insert the global indices
    assert( A.InsertGlobalIndices( MyGlobalElements[i], NumEntries, &Indices[0] )  == 0 );
  } // end i loop

  // Finish up graph construction
  assert(A.FillComplete() == 0);

  EpetraExt::CrsGraph_MapColoring
    Greedy0MapColoringTransform( EpetraExt::CrsGraph_MapColoring::GREEDY,
		                 0, false, verbose );
  Epetra_MapColoring & Greedy0ColorMap = Greedy0MapColoringTransform( A );
  printColoring(Greedy0ColorMap, &A,verbose);

  EpetraExt::CrsGraph_MapColoring
    Greedy1MapColoringTransform( EpetraExt::CrsGraph_MapColoring::GREEDY,
		                 1, false, verbose );
  Epetra_MapColoring & Greedy1ColorMap = Greedy1MapColoringTransform( A );
  printColoring(Greedy1ColorMap, &A,verbose);

  EpetraExt::CrsGraph_MapColoring
    Greedy2MapColoringTransform( EpetraExt::CrsGraph_MapColoring::GREEDY,
		                 2, false, verbose );
  Epetra_MapColoring & Greedy2ColorMap = Greedy2MapColoringTransform( A );
  printColoring(Greedy2ColorMap, &A,verbose);

  EpetraExt::CrsGraph_MapColoring
    Lubi0MapColoringTransform( EpetraExt::CrsGraph_MapColoring::LUBY,
		               0, false, verbose );
  Epetra_MapColoring & Lubi0ColorMap = Lubi0MapColoringTransform( A );
  printColoring(Lubi0ColorMap, &A,verbose);

  EpetraExt::CrsGraph_MapColoring
    Lubi1MapColoringTransform( EpetraExt::CrsGraph_MapColoring::LUBY,
		               1, false, verbose );
  Epetra_MapColoring & Lubi1ColorMap = Lubi1MapColoringTransform( A );
  printColoring(Lubi1ColorMap, &A,verbose);

  EpetraExt::CrsGraph_MapColoring
    Lubi2MapColoringTransform( EpetraExt::CrsGraph_MapColoring::LUBY,
		               2, false, verbose );
  Epetra_MapColoring & Lubi2ColorMap = Lubi2MapColoringTransform( A );
  printColoring(Lubi2ColorMap, &A,verbose);

#ifdef EPETRA_MPI
  if( verbose ) cout << "Parallel Map Coloring 1!\n";
  EpetraExt::CrsGraph_MapColoring
    Parallel1MapColoringTransform( EpetraExt::CrsGraph_MapColoring::PSEUDO_PARALLEL,
		                   0, false, verbose );
  Epetra_MapColoring & Parallel1ColorMap = Parallel1MapColoringTransform( A );
  printColoring(Parallel1ColorMap, &A,verbose);

  if( verbose ) cout << "Parallel Map Coloring 2!\n";
  EpetraExt::CrsGraph_MapColoring
    Parallel2MapColoringTransform( EpetraExt::CrsGraph_MapColoring::JONES_PLASSMAN,
		                   0, false, verbose );
  Epetra_MapColoring & Parallel2ColorMap = Parallel2MapColoringTransform( A );
  printColoring(Parallel2ColorMap, &A,verbose);
#endif


#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return 0;
}

void printColoring (const Epetra_MapColoring & ColorMap, Epetra_CrsGraph * A, bool verbose) {

  int NumColors = ColorMap.NumColors();
  int * ListOfColors = ColorMap.ListOfColors();

  EpetraExt::CrsGraph_MapColoringIndex MapColoringIndexTransform( ColorMap );
  vector<Epetra_IntVector> & ColIndices = MapColoringIndexTransform( *A );

  if( verbose )
  {

    cout << endl;
    cout << "***************************************\n";
    cout << "Column Indexing by Color:\n";
    cout << "***************************************\n";
    cout << endl;
    for( int i = 0; i < NumColors; ++i )
    {
      cout << "COLOR: " << ListOfColors[i] << endl;
      cout << ColIndices[i];
    }
    cout << endl;
  }

  return;
}
