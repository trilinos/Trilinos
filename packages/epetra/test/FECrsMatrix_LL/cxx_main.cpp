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


// Epetra_FECrsMatrix Test routine

#include "Epetra_Time.h"
#include "Epetra_Map.h"
#include "Epetra_FECrsMatrix.h"
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
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

#else

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

  int NumMyElements = 4;
  long long NumGlobalElements = NumMyElements*NumProc;
  int IndexBase = 0;
  
  Epetra_Map Map(NumGlobalElements, NumMyElements, IndexBase, Comm);

  EPETRA_TEST_ERR( Drumm1(Map, verbose),ierr);

  EPETRA_TEST_ERR( Drumm2(Map, verbose),ierr);

  EPETRA_TEST_ERR( Drumm3(Map, verbose),ierr);


  bool preconstruct_graph = false;

  EPETRA_TEST_ERR( four_quads(Comm, preconstruct_graph, verbose), ierr);

  preconstruct_graph = true;

  EPETRA_TEST_ERR( four_quads(Comm, preconstruct_graph, verbose), ierr);

  EPETRA_TEST_ERR( submatrix_formats(Comm, verbose), ierr);

  EPETRA_TEST_ERR( rectangular(Comm, verbose), ierr);

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return ierr;
}

