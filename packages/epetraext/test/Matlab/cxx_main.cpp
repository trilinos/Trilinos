/*@HEADER
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
*/

#include "Epetra_ConfigDefs.h"
#include "EpetraExt_Version.h"
#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Trilinos_Util.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_DataAccess.h"

#include "EpetraExt_MatlabEngine.h"

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
  EpetraExt::EpetraExt_MatlabEngine engine (comm);
  cout << "matlab started";

  ///* CrsMatrix test
  
  //*/

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

  
  engine.EvalString("size(TEST)", matlabBuffer, MATLABBUF);
  cout << matlabBuffer << "\n";
  engine.EvalString("TEST", matlabBuffer, MATLABBUF);
  cout << matlabBuffer << "\n";

  cout << "\n" << comm.MyPID() << " all done\n";
  return(0);
  
}
