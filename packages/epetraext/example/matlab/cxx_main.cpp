/*
//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
// ***********************************************************************
//@HEADER
*/

#include "Epetra_ConfigDefs.h"

#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_DataAccess.h"
#include "Epetra_CrsMatrix.h"

#include "EpetraExt_MatlabEngine.h" // contains the EpetraExt_MatlabEngine class

int main(int argc, char *argv[]) {

// standard Epetra MPI/Serial Comm startup	
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  
  int MyPID = comm.MyPID();
  int ierr = 0;
  bool verbose = (0 == MyPID);
  bool reportErrors = (0 == MyPID);
  // setup MatlabEngine
  if (verbose) cout << "going to startup a matlab process...\n";
  EpetraExt::EpetraExt_MatlabEngine engine (comm);
  if (verbose) cout << "matlab started\n";
  
  // setup an array of doubles to be used for the examples
  int M = 20;
  int numGlobalElements = M * comm.NumProc();
  int N = 3;
  int numMyEntries = M * N;
  double* A = new double[numMyEntries];
  double* Aptr = A;
  int startValue = numMyEntries * MyPID;

  for(int col=0; col < N; col++) {
	for(int row=0; row < M; row++) {
          *Aptr++ = startValue++;
      }
  }

  // setup an array of ints to be used for the examples
  int* intA = new int[numMyEntries];
  int* intAptr = intA;
  int intStartValue = numMyEntries * MyPID;
  for(int i=0; i < M*N; i++) {
      *intAptr++ = intStartValue++;
  }
  
  // construct a map to be used by distributed objects
  Epetra_Map map (numGlobalElements, 0, comm);
  
  // CrsMatrix example
  // constructs a globally distributed CrsMatrix and then puts it into Matlab
  if (verbose) cout << " constructing CrsMatrix...\n";
  Epetra_CrsMatrix crsMatrix (Copy, map, N);
  int* indices = new int[N];
  for (int col=0; col < N; col++) {
    indices[col] = col;	  
  }
  
  double value = startValue;
  double* values = new double[numMyEntries];
  int minMyGID = map.MinMyGID();
  for (int row=0; row < M; row++) {
    for (int col=0; col < N; col++) {
      values[col] = value++;
    }
      
    crsMatrix.InsertGlobalValues(minMyGID + row, N, values, indices);
  }
  
  crsMatrix.FillComplete();
  if (verbose) cout << " CrsMatrix constructed\n";
  if (verbose) cout << " putting CrsMatrix into Matlab as CRSM\n";
  ierr = engine.PutRowMatrix(crsMatrix, "CRSM", false);
  if (ierr) {
    if (reportErrors) cout << "There was an error in engine.PutRowMatrix(crsMatrix, \"CRSM\", false): " << ierr << endl;
  }
  
  // BlockMap example
  // puts a map into Matlab
  if (verbose) cout << " putting Map into Matlab as MAP\n";
  ierr = engine.PutBlockMap(map, "MAP", false);
  if (ierr) {
    if (reportErrors) cout << "There was an error in engine.PutBlockMap(map, \"MAP\", false);: " << ierr << endl;
  }
  
  // MultiVector example
  // constructs a globally distributed MultiVector and then puts it into Matlab
  if (verbose) cout << " constructing MultiVector...\n";
  Epetra_MultiVector multiVector (Copy, map, A, M, N);
  if (verbose) cout << " MultiVector constructed\n";
  if (verbose) cout << " putting MultiVector into Matlab as MV\n";
  ierr = engine.PutMultiVector(multiVector, "MV");
  if (ierr) {
    if (reportErrors) cout << "There was an error in engine.PutMultiVector(multiVector, \"MV\"): " << ierr << endl;
  }
  
  // SerialDenseMatrix example
  // constructs a SerialDenseMatrix on every PE
  if (verbose) cout << " constructing a SerialDenseMatrix...\n";
  Epetra_SerialDenseMatrix sdMatrix (Copy, A, M, M, N);
  if (verbose) cout << " SerialDenseMatrix constructed\n";
  if (verbose) cout << " putting SerialDenseMatrix from PE0 into Matlab as SDM_PE0\n";
  // since the third parameter is left out, the SerialDenseMatrix from PE0 is used by default
  ierr = engine.PutSerialDenseMatrix(sdMatrix, "SDM_PE0");
  if (ierr) {
    if (reportErrors) cout << "There was an error in engine.PutSerialDenseMatrix(sdMatrix, \"SDM_PE0\"): " << ierr << endl;
  }
  if (comm.NumProc() > 1) {
    if (verbose) cout << " putting SerialDenseMatrix from PE1 into Matlab as SDM_PE1\n";
    // specifying 1 as the third parameter will put the SerialDenseMatrix from PE1 into Matlab
    ierr = engine.PutSerialDenseMatrix(sdMatrix, "SDM_PE1", 1);
    if (ierr) {
      if (reportErrors) cout << "There was an error in engine.PutSerialDenseMatrix(sdMatrix, \"SDM_PE1\", 1): " << ierr << endl;
    }
  }
  

  // SerialDenseVector example
  // constructs a SerialDenseVector on every PE
  if (verbose) cout << " constructing a SerialDenseVector...\n";
  Epetra_SerialDenseVector sdVector (Copy, A, M);
  if (verbose) cout << " SerialDenseVector constructed\n";
  // since the third parameter is left out, the SerialDenseMatrix from PE0 is used by default
  if (verbose) cout << " putting SerialDenseVector from PE0 into Matlab as SDV_PE0\n";
  ierr = engine.PutSerialDenseMatrix(sdVector, "SDV_PE0");
  if (ierr) {
    if (reportErrors) cout << "There was an error in engine.PutSerialDenseMatrix(sdVector, \"SDV_PE0\"): " << ierr << endl;
  }
  if (comm.NumProc() > 1) {
    if (verbose) cout << " putting SerialDenseVector from PE1 into Matlab as SDV_PE1\n";
    // specifying 1 as the third parameter will put the SerialDenseVector from PE1 into Matlab
    ierr = engine.PutSerialDenseMatrix(sdVector, "SDV_PE1", 1);
    if (ierr) {
      if (reportErrors) cout << "There was an error in engine.PutSerialDenseMatrix(sdMatrix, \"SDV_PE1\", 1): " << ierr << endl;
    }
  }

  // IntSerialDenseMatrix example
  // constructs a IntSerialDenseMatrix on every PE
  if (verbose) cout << " constructing a IntSerialDenseMatrix...\n";
  Epetra_IntSerialDenseMatrix isdMatrix (Copy, intA, M, M, N);
  if (verbose) cout << " IntSerialDenseMatrix constructed\n";
  // since the third parameter is left out, the IntSerialDenseMatrix from PE0 is used by default
  if (verbose) cout << " putting IntSerialDenseMatrix from PE0 into Matlab as ISDM_PE0\n";
  ierr = engine.PutIntSerialDenseMatrix(isdMatrix, "ISDM_PE0");
  if (ierr) {
    if (reportErrors) cout << "There was an error in engine.PutIntSerialDenseMatrix(isdMatrix, \"ISDM_PE0\"): " << ierr << endl;
  }
  if (comm.NumProc() > 1) {
    if (verbose) cout << " putting IntSerialDenseMatrix from PE1 into Matlab as ISDM_PE1\n";
    // specifying 1 as the third parameter will put the IntSerialDenseMatrix from PE1 into Matlab
    ierr = engine.PutIntSerialDenseMatrix(isdMatrix, "ISDM_PE1", 1);
    if (ierr) {
      if (reportErrors) cout << "There was an error in engine.PutSerialDenseMatrix(isdMatrix, \"ISDM_PE1\", 1): " << ierr << endl;
    }
  }


  // IntSerialDenseVector example
  // constructs a IntSerialDenseVector on every PE
  if (verbose) cout << " constructing a IntSerialDenseVector...\n";
  Epetra_IntSerialDenseVector isdVector (Copy, intA, M);
  if (verbose) cout << " IntSerialDenseVector constructed\n";
  // since the third parameter is left out, the IntSerialDenseVector from PE0 is used by default
  if (verbose) cout << " putting IntSerialDenseVector from PE0 into Matlab as ISDV_PE0\n";
  ierr = engine.PutIntSerialDenseMatrix(isdVector, "ISDV_PE0");
  if (ierr) {
    if (reportErrors) cout << "There was an error in engine.PutIntSerialDenseMatrix(isdVector, \"ISDV_PE0\"): " << ierr << endl;
  }
  if (comm.NumProc() > 1) {
    if (verbose) cout << " putting IntSerialDenseVector from PE1 into Matlab as ISDV_PE1\n";
    // specifying 1 as the third parameter will put the IntSerialDenseVector from PE1 into Matlab
    ierr = engine.PutIntSerialDenseMatrix(isdVector, "ISDV_PE1", 1);
    if (ierr) {
      if (reportErrors) cout << "There was an error in engine.PutSerialDenseMatrix(isdVector, \"ISDV_PE1\", 1): " << ierr << endl;
    }
  }
  
  // entering a while loop on PE0 will keep the Matlab workspace alive
  /*
  if (MyPID == 0)
  while(1) {
    // do nothing
  }
  */

  const int bufSize = 200;
  char s [bufSize];
  const int matlabBufferSize = 1024 * 16;
  char matlabBuffer [matlabBufferSize];
  
  // send some commands to Matlab and output the result to stdout
  engine.EvalString("whos", matlabBuffer, matlabBufferSize);
  if (verbose) cout << matlabBuffer << endl;
  engine.EvalString("SDV_PE0", matlabBuffer, matlabBufferSize);
  if (verbose) cout << matlabBuffer << endl;
  if (comm.NumProc() > 1) {
    engine.EvalString("SDV_PE1", matlabBuffer, matlabBufferSize);
    if (verbose) cout << matlabBuffer << endl;
  }
  
  // the following allows user interaction with Matlab
  if (MyPID == 0)
  while(1) {
      // Prompt the user and get a string
      printf(">> ");
      if (fgets(s, bufSize, stdin) == NULL) {
          printf("Bye\n");
          break ;
      }
      printf ("command :%s:\n", s) ;
      
      // send the command to MATLAB
      // output goes to stdout
      ierr = engine.EvalString(s, matlabBuffer, matlabBufferSize);
      if (ierr != 0) {
          printf("there was an error: %d", ierr);
		  ierr = 0;
      }
      else {
      	  printf("Matlab Output:\n%s", matlabBuffer);
      }
  }
  
  if (verbose) cout << endl << " all done\n";

// standard finalizer for Epetra MPI Comms
#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  // we need to delete engine because the MatlabEngine finalizer shuts down the Matlab process associated with this example
  // if we don't delete the Matlab engine, then this example application will not shut down properly
  delete &engine;
  return(0);
}
