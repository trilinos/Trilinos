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

// Epetra_FECrsGraph Test routine

#include "Epetra_Time.h"
#include "Epetra_Map.h"
#include "ExecuteTestProblems.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "../epetra_test_err.h"
#include "Epetra_Version.h"

int main(int argc, char *argv[]) {

  int ierr = 0;

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

#else

  int size = 1; // Serial case (not using MPI)
  int rank = 0;
  Epetra_SerialComm Comm;

#endif

//  Comm.SetTracebackMode(0); // This should shut down any error tracing
  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

#ifdef EPETRA_MPI
  int localverbose = verbose ? 1 : 0;
  int globalverbose=0;
  MPI_Allreduce(&localverbose, &globalverbose, 1, MPI_INT, MPI_SUM,
		MPI_COMM_WORLD);
  verbose = (globalverbose>0);
#endif

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc(); 

  if (verbose && MyPID==0)
    cout << Epetra_Version() << endl << endl;

  if (verbose) cout << Comm <<endl;

  // Redefine verbose to only print on PE 0
  //if (verbose && rank!=0) verbose = false;

  int NumVectors = 1;
  int NumMyElements = 4;
  int NumGlobalElements = NumMyElements*NumProc;
  int IndexBase = 0;
  
  Epetra_Map Map(NumGlobalElements, NumMyElements, IndexBase, Comm);

  EPETRA_TEST_ERR( Drumm1(Map, verbose),ierr);

  EPETRA_TEST_ERR( Drumm2(Map, verbose),ierr);


  bool preconstruct_graph = false;

  EPETRA_TEST_ERR( four_quads(Comm, preconstruct_graph, verbose), ierr);

  preconstruct_graph = true;

  EPETRA_TEST_ERR( four_quads(Comm, preconstruct_graph, verbose), ierr);

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return ierr;
}

