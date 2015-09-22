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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

// CrsGraph_BTF Test routine

#include <EpetraExt_ConfigDefs.h>
#include "EpetraExt_Version.h"

#include "Epetra_Time.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#endif

#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "EpetraExt_AmesosBTF_CrsMatrix.h"

int main(int argc, char *argv[]) {

  int returnierr=0;

  using std::cout;
  using std::endl;
  using std::flush;

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (size > 1) {
    cout << "This example cannot be run on more than one processor!" << endl;
    MPI_Finalize();
    returnierr = -1;
    return returnierr;
  }

#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;


#ifdef EPETRA_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  if (!verbose) Comm.SetTracebackMode(0); // This should shut down any error traceback reporting

  if (verbose) {
    cout << EpetraExt::EpetraExt_Version() << endl << endl;
    cout << Comm << endl << flush;
  }

  Comm.Barrier();

  int NumMyElements = 3;
 
  Epetra_Map Map( NumMyElements, 0, Comm );
  
  Epetra_CrsGraph Graph( Copy, Map, 1 );

  int index[2];
  index[0] = 2;
  Graph.InsertGlobalIndices( 0, 1, &index[0] );
  index[0] = 0;
  index[1] = 2;
  Graph.InsertGlobalIndices( 1, 2, &index[0] );
  index[0] = 1;
  Graph.InsertGlobalIndices( 2, 1, &index[0] );

  Graph.FillComplete();

  // Create an Epetra::CrsMatrix
  Epetra_CrsMatrix Matrix( Copy, Graph );
  double value[2];
  index[0] = 2; value[0] = 3.0;
  Matrix.ReplaceMyValues( 0, 1, &value[0], &index[0] );
  index[0] = 0; index[1] = 2;
  value[0] = 2.0; value[1] = 2.5;
  Matrix.ReplaceMyValues( 1, 2, &value[0], &index[0] );
  index[0] = 1; value[0] = 1.0;
  Matrix.ReplaceMyValues( 2, 1, &value[0], &index[0] );
  Matrix.FillComplete();

  EpetraExt::AmesosBTF_CrsMatrix BTFTrans( 0.0, true, verbose );
  Epetra_CrsMatrix & NewMatrix = BTFTrans( Matrix );

  if (verbose) {
    cout << "*************** PERFORMING BTF TRANSFORM ON CRS_MATRIX **************" <<endl<<endl;
    cout << "CrsMatrix *before* BTF transform: " << endl << endl;
    cout << Matrix << endl;
  }

  BTFTrans.fwd();

  if (verbose) {
    cout << "CrsMatrix *after* BTF transform: " << endl << endl;
    cout << NewMatrix << endl;
  }

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif
  
  return returnierr;
}
  
