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


// Epetra_FEVector Test routine

#include "Epetra_Time.h"
#include "Epetra_BlockMap.h"
#include "Epetra_FEVector.h"
#include "ExecuteTestProblems.h"
#include "../src/Epetra_test_functions.h"
#include "Epetra_Comm.h"
#include "../epetra_test_err.h"
#include "Epetra_Version.h"

int main(int argc, char *argv[]) {

  int ierr = 0;

  //epetra_test::create_comm initializes mpi if EPETRA_MPI is defined
  Epetra_Comm* epetra_comm = epetra_test::create_comm(argc, argv);
  Epetra_Comm& Comm = *epetra_comm;

  int NumProc = Comm.NumProc(); 

  bool verbose =
    epetra_test::global_check_for_flag_on_proc_0("-v",argc,argv, Comm);

  if (verbose) cout << Epetra_Version() << endl << endl;

  int NumVectors = 1;
  int NumMyElements = 4;
  int NumGlobalElements = NumMyElements*NumProc;
  int IndexBase = 0;
  int ElementSize = 1;

  Epetra_BlockMap BlockMap(NumGlobalElements, NumMyElements,
                           ElementSize, IndexBase, Comm);
  BlockMap.SetTracebackMode(0); // This should shut down any error tracing

  EPETRA_TEST_ERR(MultiVectorTests(BlockMap, NumVectors, verbose),ierr);

  Comm.Barrier();
  if (verbose)cout << endl;
  Comm.Barrier();

  EPETRA_TEST_ERR( fevec0(Comm, verbose), ierr);

  Comm.Barrier();
  if (verbose)cout << endl;
  Comm.Barrier();

  EPETRA_TEST_ERR( fevec1(Comm, verbose), ierr);

  Comm.Barrier();
  if (verbose)cout << endl;
  Comm.Barrier();

  EPETRA_TEST_ERR( fevec2(Comm, verbose), ierr);

  Comm.Barrier();
  if (verbose)cout << endl;
  Comm.Barrier();

  EPETRA_TEST_ERR( fevec3(Comm, verbose), ierr);

  Comm.Barrier();
  if (verbose)cout << endl;
  Comm.Barrier();

  EPETRA_TEST_ERR( fevec4(Comm, verbose), ierr);

  Comm.Barrier();
  if (verbose)cout << endl;
  Comm.Barrier();

  EPETRA_TEST_ERR( fevec5(Comm, verbose), ierr);

  Comm.Barrier();
  if (verbose)cout << endl;
  Comm.Barrier();

  EPETRA_TEST_ERR( fevec6(Comm, verbose), ierr);
  delete epetra_comm;

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif
 
  return ierr;
}

