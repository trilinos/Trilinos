/*
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
*/

#include "Epetra_ConfigDefs.h"
#include "Epetra_IntMultiVector.h"
#include "Epetra_BlockMap.h"
#include "Epetra_CombineMode.h"
#include "Epetra_Import.h"
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
  int global_ierr = 0;
  int numErr = 0;

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int rank; // My process ID

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );

#else

  int rank = 0;
  Epetra_SerialComm Comm;

#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if(argc > 1) {
    if(argv[1][0]=='-' && argv[1][1]=='v') {
      verbose = true;
    }
  }

  Comm.SetTracebackMode(0); // This should shut down any error traceback reporting
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  if(verbose && MyPID==0)
    cout << Epetra_Version() << endl << endl;

  if(verbose) cout << "Processor " << MyPID << " of " << NumProc << " is alive." << endl;

  // Redefine verbose to only print on PE 0
  if(verbose && rank != 0)
		verbose = false;

  int NumMyEquations = 6;
  int NumGlobalEquations = NumMyEquations*NumProc+EPETRA_MIN(NumProc,3);
  if(MyPID < 3)
    NumMyEquations++;

  // Construct a Map that puts approximately the same Number of equations on each processor

  Epetra_BlockMap Map(NumGlobalEquations, NumMyEquations, 1, 0, Comm);

  // Get update list and number of local equations from newly created Map
  std::vector<int> MyGlobalElements(Map.NumMyElements());
  Map.MyGlobalElements(MyGlobalElements.data());

  // Construct the IntMultiVector
  const int numVecs = 3;
  Epetra_IntMultiVector myIMV(Map, numVecs);

  // Extract a view to populate the IntMultiVector
  int myLDA;
  int* myVals;
  ierr = myIMV.ExtractView(&myVals, &myLDA);
  Comm.MaxAll(&ierr, &global_ierr, 1);
  EPETRA_TEST_ERR(global_ierr, numErr);
  for(int i = 0; i < myIMV.MyLength(); ++i) {
    for(int j = 0; j < numVecs; ++j) {
      myVals[j*myLDA + i] = numVecs*i + j;
    }
  }

  // Test: replace local values
  myIMV.ReplaceMyValue(3, 2, 42);
  ierr = (myIMV[2][3] == 42 ? 0 : 1);
  Comm.MaxAll(&ierr, &global_ierr, 1);
  EPETRA_TEST_ERR(global_ierr, numErr);

  // Test: sum into local values
  myIMV.SumIntoMyValue(3, 1, 42);
  ierr = (myIMV[1][3] == 3*numVecs + 1 + 42 ? 0 : 1);
  Comm.MaxAll(&ierr, &global_ierr, 1);
  EPETRA_TEST_ERR(global_ierr, numErr);

  // Test: replace and sum global values on rank 2
  if(MyPID == 2) {
    myIMV.ReplaceGlobalValue(Map.GID(3), 2, 48);
    ierr = (myIMV[2][3] == 48 ? 0 : 1);
    EPETRA_TEST_ERR(ierr, numErr);
    myIMV.SumIntoGlobalValue(Map.GID(3), 1, 48);
    ierr = (myIMV[1][3] == 3*numVecs + 1 + 42 + 48 ? 0 : 1);
    EPETRA_TEST_ERR(ierr, numErr);
  }
  Comm.MaxAll(&ierr, &global_ierr, 1);
  EPETRA_TEST_ERR(global_ierr, numErr);

  // Construct a Map that puts all the elements on proc 0
  Epetra_BlockMap Map0(NumGlobalEquations, (MyPID == 0 ? NumGlobalEquations : 0), 1, 0, Comm);
  Epetra_IntMultiVector myIMV0(Map0, numVecs);

  Epetra_Import myImporter(Map0, Map);
  ierr = myIMV0.Import(myIMV, myImporter, Epetra_CombineMode::Insert);
  Comm.MaxAll(&ierr, &global_ierr, 1);
  EPETRA_TEST_ERR(global_ierr, numErr);

  if(MyPID == 0) {
    ierr = (myIMV0[1][ 3] == 3*numVecs + 1 + 42      ? 0 : 1);
    EPETRA_TEST_ERR(ierr, numErr);
    ierr = (myIMV0[1][10] == 3*numVecs + 1 + 42      ? 0 : 1);
    EPETRA_TEST_ERR(ierr, numErr);
    ierr = (myIMV0[1][17] == 3*numVecs + 1 + 42 + 48 ? 0 : 1);
    EPETRA_TEST_ERR(ierr, numErr);
    ierr = (myIMV0[1][24] == 3*numVecs + 1 + 42      ? 0 : 1);
    EPETRA_TEST_ERR(ierr, numErr);
  }


#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif
/* end main */
	return(ierr);
}
