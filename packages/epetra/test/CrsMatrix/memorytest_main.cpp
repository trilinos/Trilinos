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

// This program tests the memory management system of the class CrsMatrix (memory leak, invalid free).
// It should be run using valgrind.

#include <vector>

#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
//#include "../epetra_test_err.h"
#include "Epetra_Version.h"

int main(int argc, char *argv[])
{
  int ierr = 0;

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc, &argv);
  int rank; // My process ID

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );

#else

  int rank = 0;
  Epetra_SerialComm Comm;

#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  int verbose_int = verbose ? 1 : 0;
  Comm.Broadcast(&verbose_int, 1, 0);
  verbose = verbose_int==1 ? true : false;

  Comm.SetTracebackMode(0); // This should shut down any error traceback reporting
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  if(verbose && MyPID==0)
    cout << Epetra_Version() << std::endl << std::endl;

  if (verbose) cout << "Processor "<<MyPID<<" of "<< NumProc
		    << " is alive."<<endl;

  // unused: bool verbose1 = verbose;

  // Redefine verbose to only print on PE 0
  if(verbose && rank!=0) 
    verbose = false;

  if (verbose) cout << "Test the memory management system of the class CrsMatrix (memory leak, invalid free)" << std::endl;

  //
  // Test 1
  //
  
  if(Comm.NumProc() == 1) { // this is a sequential test

    if (verbose) cout << "* Using Copy, ColMap, Variable number of indices per row and Static profile (cf. bug #5499)" << std::endl;
    
    // Row Map
    Epetra_Map RowMap(2, 0, Comm);
    
    // ColMap  
    std::vector<int> colids(2);
    colids[0]=0;
    colids[1]=1;
    Epetra_Map ColMap(-1, 2, &colids[0], 0, Comm);

    // NumEntriesPerRow
    std::vector<int> NumEntriesPerRow(2);
    NumEntriesPerRow[0]=2;
    NumEntriesPerRow[1]=2;
    
    // Test
    Epetra_CrsMatrix A(Copy, RowMap, ColMap, &NumEntriesPerRow[0], true);
    A.FillComplete();
    
  }

  /*
    if (bool) {
    if (verbose) cout << endl << "tests FAILED" << endl << endl;
    }
    else {*/
  if (verbose) cout << endl << "tests PASSED" << endl << endl;
  /*    } */

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return ierr;
}
