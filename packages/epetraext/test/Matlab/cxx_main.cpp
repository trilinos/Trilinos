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
//#include "EpetraExt_Version.h"
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

#include "EpetraExt_MatlabEngine.h"

// the following deal with matlab provided headers:
#include "engine.h"
#include "mex.h"
#undef printf  // matlab has its own printf that we don't want to use

#define BUFSIZE 200
#define MATLABBUF 1024 * 16

int main(int argc, char *argv[]) {
cout << "going to setup MPI...\n";

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  cout << "mpit setup complete\n";

  int MyPID = comm.MyPID();

  char s [BUFSIZE] ;
  char matlabBuffer [MATLABBUF];
  cout << "going to init matlab\n";
  EpetraExt::EpetraExt_MatlabEngine * enginePtr = new EpetraExt::EpetraExt_MatlabEngine(comm);
  EpetraExt::EpetraExt_MatlabEngine & engine = *enginePtr;
  cout << "matlab started\n";
  
  /* GetCrsMatrix test
  engine.EvalString("CRSM=sparse(eye(8,10))", matlabBuffer, MATLABBUF);
  cout << matlabBuffer << endl;
  int myM=4;
  int M = myM * comm.NumProc();
  int N = 10;  
  Epetra_Map getMap (M, 0, comm);
  Epetra_Map colMap(N, N, 0, comm);
  Epetra_CrsMatrix getCRSM (Copy, getMap, colMap, N);
  double colValue = 0;
  for(int row=myM*MyPID; row < myM*(MyPID+1); row++) {
    getCRSM.InsertGlobalValues(row, 1, &colValue, &row);
  }
  getCRSM.FillComplete(colMap, getMap);
  //getCRSM.FillComplete();
  int ierr = engine.GetCrsMatrix("CRSM", getCRSM, false);
  if (ierr) {
    cout << "engine.GetCrsMatrix(\"CRSM\", getCRSM, false) failed" << endl;
  }
  
  cout << getCRSM << endl;
  
  engine.EvalString("whos", matlabBuffer, MATLABBUF);
  cout << matlabBuffer << endl;
  */
  
  /* GetIntSerialDenseMatrix test
  engine.EvalString("ISDM=rand(8,2)*100", matlabBuffer, MATLABBUF);
  cout << matlabBuffer << endl;
  int procToGet = 1;
  int M = 8;
  int N = 2;  
  int* A = new int[M*N];
  Epetra_IntSerialDenseMatrix getISDM (View, A, M, M, N);
  int ierr = engine.GetIntSerialDenseMatrix("ISDM", getISDM, procToGet);
  if (ierr) {
    cout << "engine.GetIntSerialDenseMatrix(\"ISDM\", getISDM, procToGet) failed" << endl;
  }
  
  if (MyPID == 1) cout << getISDM << endl;
  */
  
  /* GetSerialDenseMatrix test
  engine.EvalString("SDM=rand(8,2)", matlabBuffer, MATLABBUF);
  cout << matlabBuffer << endl;
  int procToGet = 1;
  int M = 8;
  int N = 2;  
  double* A = new double[M*N];
  Epetra_SerialDenseMatrix getSDM (View, A, M, M, N);
  int ierr = engine.GetSerialDenseMatrix("SDM", getSDM, procToGet);
  if (ierr) {
    cout << "engine.GetSerialDenseMatrix(\"SDM\", getSDM, procToGet) failed" << endl;
  }
  
  if (MyPID == 1) cout << getSDM << endl;
  */
  
  /* GetMultiVector test
  if (comm.NumProc() != 2) {
    if (MyPID == 0) cout << "Error: this test must be run with exactly two PE." << endl;
    delete &engine;
    #ifdef EPETRA_MPI
    MPI_Finalize();
    #endif
    return(-1);
  }
  engine.EvalString("MV=rand(8,2)", matlabBuffer, MATLABBUF);
  cout << matlabBuffer << endl;
  int myM = 4;
  int M = myM * comm.NumProc();
  int N = 2;
  Epetra_Map getMap (M, 0, comm);
  double* A = new double[myM*N];
  Epetra_MultiVector getMV (View, getMap, A, myM, N);
  cout << "MultiVector created" << endl;
  int ierr = engine.GetMultiVector("MV", getMV);
  if (ierr) {
    cout << "engine.GetMultiVector(\"MV\", getMV) failed" << endl;
  }
  
  cout << getMV << endl;
  */
  
  /* CrsMatrix test
  int numGlobalElements = 8;
  int numMyElements = 8/comm.NumProc();
  int M=numGlobalElements/comm.NumProc();
  int N=10;
  int* myGlobalElements = new int[M];
  
  int minGID = 0;
  int startIndex = minGID + M*comm.MyPID();
  for(int i=0; i < M; i++) {
      myGlobalElements[i] = startIndex++;
  }
  

  int* colMapGIDs = new int[N];
  for (int i=0; i < N; i++) {
    colMapGIDs[i] = i;
  }
  //Epetra_Map map (numGlobalElements, numMyElements, myGlobalElements, minGID, comm);
  Epetra_Map map (numGlobalElements, minGID, comm);
  Epetra_Map colMap (N, N, colMapGIDs, 0, comm);		   

  Epetra_CrsMatrix crsMatrix (Copy, map, colMap, N);
  
  cout << "crs matrix created\n";
  
  //int indices[8] = {-4,-3,-2,-1,0,1,2,3};
  int indices[10] = {0,1,2,3,4,5,6,7,8,9};
  double* values = new double[N];
  double value = M * N * comm.MyPID();
  for (int i=0; i < M; i++) {
    if (i % 2 == 0) {
      for (int j=0; j < N; j++) {
      	values[j] = value++;
      }
      
      crsMatrix.InsertGlobalValues(myGlobalElements[i], N, values, indices);
    }
  }
  
  cout << "done putting values\n";
  crsMatrix.FillComplete(colMap, map);
  cout << "done filling crsMatrix and calling crsMatrix.FillComplete()\n";
  
  cout << crsMatrix;
  
  cout << "done printing crsMatrix\n";
  
  //cout << map;
  //cout << "done printing map\n";
  
  int ierr = engine.PutRowMatrix(crsMatrix, "TEST", true);
  //int ierr = engine.PutBlockMap(map, "TEST", true);
  //cout << "done calling engine.PutRowMatrix(crsMatrix, \"TEST\", false)\n";
  if (ierr != 0) {
    cout << "engine.PutRowMatrix(crsMatrix, \"TEST\") returned nonzero result: " << ierr << "\n";
    return(-1);
  }
  
  */

  /* MultiVector test
  cout << MyPID << " going to do multivector test...\n";
  int numGlobalElements = 100;
  int M = numGlobalElements/comm.NumProc();
  int N = 3;
  int numMyElements = M * N;
  double* A = new double[numMyElements];
  double* Aptr = A;
  int startValue = 0;

  cout << MyPID << " allocated space for A, now filling A\n";
  for(int col=0; col < N; col++) {
	startValue = (col * numGlobalElements) + (M * MyPID);
	for(int row=0; row < M; row++) {
          *Aptr++ = row+startValue;
      }
  }
  cout << MyPID << " A filled\n";

  Epetra_Map map (numGlobalElements, 0, comm);
  Epetra_MultiVector multiVector (Copy, map, A, M, N);
  //cout << multiVector;
  engine.PutMultiVector(multiVector, "TEST");
  */
  
  /*SerialDenseMatrix test
  cout << MyPID << " going to do SerialDenseMatrix test...\n";
  double* A = new double[30];
  cout << MyPID << " allocated space for A, now filling A\n";
  double* Aptr = A;
  int M = 5;
  int N = 6;
  int startValue = M*N*comm.MyPID();
  for(int i=0; i < M*N; i++) {
      *Aptr++ = i + startValue;
  }
  cout << MyPID << " A filled\n";

  Epetra_SerialDenseMatrix sdMatrix (View, A, M, M, N);
  engine.PutSerialDenseMatrix(sdMatrix, "TEST", 0);
  cout << sdMatrix;
  */

  /* SerialDenseVector test
  double* A = new double[30];
  double* Aptr = A;
  int length = 30;
  for(int i=0; i < length; i++) {
      *Aptr++ = i;
  }

  Epetra_SerialDenseVector sdVector (Copy, A, length);
  engine.PutSerialDenseMatrix(sdVector, "SDVECTOR");
  cout << sdVector;
  */

  /*IntSerialDenseMatrix test
  cout << MyPID << " going to do IntSerialDenseMatrix test...\n";
  int* A = new int[30];
  cout << MyPID << " allocated space for A, now filling A\n";
  int* Aptr = A;
  int M = 5;
  int N = 6;
  int startValue = M*N*comm.MyPID();
  for(int i=0; i < M*N; i++) {
      *Aptr++ = i + startValue;
  }
  cout << MyPID << " A filled\n";
  Epetra_IntSerialDenseMatrix isdMatrix (Copy, A, M, M, N);
  cout << isdMatrix;
  engine.PutIntSerialDenseMatrix(isdMatrix, "TEST", 0);
  */


  /* SerialDenseVector test
  int* A = new int[30];
  int* Aptr = A;
  int length = 30;
  for(int i=0; i < length; i++) {
      *Aptr++ = i;
  }

  Epetra_IntSerialDenseVector isdVector (Copy, A, length);
  engine.PutIntSerialDenseMatrix(isdVector, "ISDVECTOR");
  cout << isdVector;
  */

  /*while(1) {

	// do nothing
	}*/

  /*if (comm.NumProc() == 1) {
  int err;
  while(1) {
      // Prompt the user and get a string
      printf(">> ");
      if (fgets(s, BUFSIZE, stdin) == NULL) {
          printf("Bye\n");
          break ;
      }
      printf ("command :%s:\n", s) ;
      
      // Send the command to MATLAB
      // output goes to stdout
      err = engine.EvalString(s, matlabBuffer, MATLABBUF);
      if (err != 0) {
          printf("there was an error: %d", err);
		  err = 0;
      }
      else {
      	  printf("Matlab Output:\n%s", matlabBuffer);
      }
  }
  }*/

  //delete engine ;

  /*
  engine.EvalString("size(TEST)", matlabBuffer, MATLABBUF);
  cout << matlabBuffer << "\n";
  engine.EvalString("TEST", matlabBuffer, MATLABBUF);
  cout << matlabBuffer << "\n";
  */
  
  cout << "\n" << comm.MyPID() << " all done\n";
  
#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  delete enginePtr;

  return(0);
 
}
