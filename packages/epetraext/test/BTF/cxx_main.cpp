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
#include "EpetraExt_BTF_CrsGraph.h"
#include "../epetra_test_err.h"


int main(int argc, char *argv[]) {

  int i, ierr=0, returnierr=0;

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#else

  int size = 1; // Serial case (not using MPI)
  int rank = 0;

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

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  if (verbose) {
    cout << EpetraExt::EpetraExt_Version() << endl << endl;
    cout << Comm << endl << flush;
  }

  Comm.Barrier();
  bool verbose1 = verbose;

  if (verbose) verbose = (MyPID==0);

  int NumMyElements = 3;
  int NumGlobalElements = NumMyElements;
  int IndexBase = 0;
 
  cout << "MyPID: " << MyPID << ", NumMyElements: " << NumMyElements << endl;

  Epetra_Map Map( NumMyElements, 0, Comm );
  cout << Map << endl;
  
  Epetra_CrsGraph Graph( Copy, Map, 1 );

  int index = 2;
  Graph.InsertGlobalIndices( 0, 1, &index );
  index = 0;
  Graph.InsertGlobalIndices( 1, 1, &index );
  index = 1;
  Graph.InsertGlobalIndices( 2, 1, &index );

  Graph.TransformToLocal();
  cout << Graph << endl;

  EpetraExt::CrsGraph_BTF BTFTransform;
  Epetra_CrsGraph & NewGraph = BTFTransform( Graph );

  cout << NewGraph << endl;

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return returnierr;
}

