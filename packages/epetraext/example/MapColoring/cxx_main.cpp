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
#include "Epetra_CrsMatrix.h"
#include "Epetra_IntVector.h"
#include "Epetra_MapColoring.h"
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"
#include "EpetraExt_BlockMapIn.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "../../test/epetra_test_err.h"

#include <vector>

int main(int argc, char *argv[]) {

  int returnierr=0;

#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init( &argc, &argv );
#endif

  bool verbose = false;

  if( argc < 2 || argc > 3 )
  {
    cout << "Usage: " << argv[0] << " [-v] base_name \n"
         << "\t{ Files: base_name_matrix contains matrix-market data and \n"
         << "\t  base_name_map contains the needed map.\n" << endl;
    exit(1);
  }

  int loc = 1;
  // Check if we should print results to standard out
  if(argv[loc][0]=='-' && argv[loc][1]=='v')
  { verbose = true; ++loc; }

#ifdef EPETRA_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();

  bool verbose1 = false;
  if(verbose) verbose1 = (MyPID==0);

  if(verbose1)
    cout << EpetraExt::EpetraExt_Version() << endl << endl;

  Comm.Barrier();

  string filename1(argv[loc++]);
  string filename2 = filename1 + "_map";
  filename1 += "_matrix";

  if(verbose1) cout << "Reading Epetra_BlockMap from file : " << filename2;
  Epetra_BlockMap* BlkMap(0);
  EpetraExt::MatrixMarketFileToBlockMap(filename2.c_str(), Comm, BlkMap);
  if(verbose1) cout << " Done." << endl;

  if(verbose1) cout << "Converting Epetra_BlockMap to Epetra_Map ... ";
  Epetra_Map* Map = new Epetra_Map( BlkMap->NumGlobalElements(),
                                    BlkMap->NumMyElements(),
                                    BlkMap->MyGlobalElements(),
                                    BlkMap->IndexBase(),
                                    BlkMap->Comm());
  if(verbose1) cout << " Done." << endl;

  if(verbose1) cout << "Reading Epetra_CrsMatrix from file : " << filename1;
  Epetra_CrsMatrix* A_Matrix(0);
  EpetraExt::MatrixMarketFileToCrsMatrix( filename1.c_str(), *Map, A_Matrix);
  if(verbose1) cout << " Done." << endl;

  A_Matrix->TransformToLocal();

  if(verbose) cout << Comm << endl << flush;
  Comm.Barrier();

  // Obtain the Epetra Compressed Row Sparse Graph
  const Epetra_CrsGraph & constA = A_Matrix->Graph();
  Epetra_CrsGraph & A = const_cast<Epetra_CrsGraph&>(constA);

  // Finish up graph construction
  assert(A.TransformToLocal() == 0);

  EpetraExt::CrsGraph_MapColoring
    Greedy0MapColoringTransform( EpetraExt::CrsGraph_MapColoring::ALGO_GREEDY,
		                verbose, 0 );
  Epetra_MapColoring & Greedy0ColorMap = Greedy0MapColoringTransform( A );

  EpetraExt::CrsGraph_MapColoring
    Greedy1MapColoringTransform( EpetraExt::CrsGraph_MapColoring::ALGO_GREEDY,
		                verbose, 1 );
  Epetra_MapColoring & Greedy1ColorMap = Greedy1MapColoringTransform( A );

  EpetraExt::CrsGraph_MapColoring
    Greedy2MapColoringTransform( EpetraExt::CrsGraph_MapColoring::ALGO_GREEDY,
		                verbose, 2 );
  Epetra_MapColoring & Greedy2ColorMap = Greedy2MapColoringTransform( A );

  EpetraExt::CrsGraph_MapColoring
    Lubi0MapColoringTransform( EpetraExt::CrsGraph_MapColoring::ALGO_LUBI,
		               verbose, 0);
  Epetra_MapColoring & Lubi0ColorMap = Lubi0MapColoringTransform( A );

  EpetraExt::CrsGraph_MapColoring
    Lubi1MapColoringTransform( EpetraExt::CrsGraph_MapColoring::ALGO_LUBI,
		               verbose, 1);
  Epetra_MapColoring & Lubi1ColorMap = Lubi1MapColoringTransform( A );

  EpetraExt::CrsGraph_MapColoring
    Lubi2MapColoringTransform( EpetraExt::CrsGraph_MapColoring::ALGO_LUBI,
		               verbose, 2);
  Epetra_MapColoring & Lubi2ColorMap = Lubi2MapColoringTransform( A );

  if( verbose ) cout << "Parallel Map Coloring 1!\n";
  EpetraExt::CrsGraph_MapColoring
    Parallel1MapColoringTransform( EpetraExt::CrsGraph_MapColoring::ALGO_GREEDY,
		               verbose, 0, true);
  Epetra_MapColoring & Parallel1ColorMap = Parallel1MapColoringTransform( A );

  if( verbose ) cout << "Parallel Map Coloring 2!\n";
  EpetraExt::CrsGraph_MapColoring
    Parallel2MapColoringTransform( EpetraExt::CrsGraph_MapColoring::ALGO_GREEDY,
		               verbose, 0, true, false);
  Epetra_MapColoring & Parallel2ColorMap = Parallel2MapColoringTransform( A );

  int NumColors = Greedy0ColorMap.NumColors();
  int * ListOfColors = Greedy0ColorMap.ListOfColors();

  EpetraExt::CrsGraph_MapColoringIndex MapColoringIndexTransform( Greedy0ColorMap );
  vector<Epetra_IntVector> & ColIndices = MapColoringIndexTransform( A );

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

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return returnierr;
}

