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


#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Flops.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "../epetra_test_err.h"
#include "Epetra_Version.h"

// prototypes

int main(int argc, char *argv[])
{
  int ierr = 0;

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
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  int verbose_int = verbose ? 1 : 0;
  Comm.Broadcast(&verbose_int, 1, 0);
  verbose = verbose_int==1 ? true : false;


  //  char tmp;
  //  if (rank==0) cout << "Press any key to continue..."<< std::endl;
  //  if (rank==0) cin >> tmp;
  //  Comm.Barrier();

  Comm.SetTracebackMode(0); // This should shut down any error traceback reporting
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  if(verbose && MyPID==0)
    cout << Epetra_Version() << std::endl << std::endl;

  if (verbose) cout << "Processor "<<MyPID<<" of "<< NumProc
		    << " is alive."<<endl;

  bool verbose1 = verbose;

  // Redefine verbose to only print on PE 0
  if(verbose && rank!=0) verbose = false;
  
  int NumMyEquations = 1;
  long long NumGlobalEquations = NumProc; 
  
  // Get update list and number of local equations from newly created Map
  long long* MyGlobalElementsLL = new long long[NumMyEquations];
  
  MyGlobalElementsLL[0] = 2000000000+MyPID;

  // Construct a Map that puts approximately the same Number of equations on each processor

  Epetra_Map MapLL(NumGlobalEquations, NumMyEquations, MyGlobalElementsLL, 0, Comm);

  EPETRA_TEST_ERR(MapLL.GlobalIndicesInt(),ierr);
  EPETRA_TEST_ERR(!(MapLL.GlobalIndicesLongLong()),ierr);

  // Create an integer vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation on this processor

  int* NumNzLL = new int[NumMyEquations];
  NumNzLL[0] = 0;  

  // Create int types meant to add to long long matrix for test of failure
  int* MyIntGlobalElementsLL = new int[NumMyEquations]; 
  MyIntGlobalElementsLL[0] = 20000+MyPID;

  // Create a long long Epetra_Matrix
   Epetra_CrsMatrix A_LL(Copy, MapLL, NumNzLL);
  EPETRA_TEST_ERR(A_LL.IndicesAreGlobal(),ierr);
  EPETRA_TEST_ERR(A_LL.IndicesAreLocal(),ierr);

  // Insert values
  double one = 1.0;
  // Try to add ints which should fail and be caught as an int
  try {
    A_LL.InsertGlobalValues(MyIntGlobalElementsLL[0], 1, &one, MyIntGlobalElementsLL+0);
  } catch(int i) {
    EPETRA_TEST_ERR(!(i==-1),ierr);
  }
  // Add long longs which should succeed
  EPETRA_TEST_ERR(!(A_LL.InsertGlobalValues(MyGlobalElementsLL[0], 1, &one, MyGlobalElementsLL+0)==0),ierr);
  EPETRA_TEST_ERR(!(A_LL.IndicesAreGlobal()),ierr);
  EPETRA_TEST_ERR(!(A_LL.FillComplete(false)==0),ierr);
  EPETRA_TEST_ERR(!(A_LL.IndicesAreLocal()),ierr);
  

  // Get update list and number of local equations from newly created Map
  int* MyGlobalElementsInt = new int[NumMyEquations];
  
  MyGlobalElementsInt[0] = 2000+MyPID;

  // Construct a Map that puts approximately the same Number of equations on each processor

  Epetra_Map MapInt(NumGlobalEquations, NumMyEquations, MyGlobalElementsInt, 0, Comm);

  EPETRA_TEST_ERR(!(MapInt.GlobalIndicesInt()),ierr);
  EPETRA_TEST_ERR(MapInt.GlobalIndicesLongLong(),ierr);

  // Create an integer vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation on this processor

  int* NumNzInt = new int[NumMyEquations];
  NumNzInt[0] = 0; 

  // Create int types meant to add to long long matrix for test of failure
  long long* MyLLGlobalElementsInt = new long long[NumMyEquations]; 
  MyLLGlobalElementsInt[0] = 2000000000+MyPID;

  // Create a int Epetra_Matrix
  Epetra_CrsMatrix A_Int(Copy, MapInt, NumNzInt);
  EPETRA_TEST_ERR(A_Int.IndicesAreGlobal(),ierr);
  EPETRA_TEST_ERR(A_Int.IndicesAreLocal(),ierr);

  // Insert values
  try {
    A_Int.InsertGlobalValues(MyLLGlobalElementsInt[0], 1, &one, MyLLGlobalElementsInt+0);
  } catch(int i) {
    EPETRA_TEST_ERR(!(i==-1),ierr);
  }
  // Add long longs which should succeed
  EPETRA_TEST_ERR(!(A_Int.InsertGlobalValues(MyGlobalElementsInt[0], 1, &one, MyGlobalElementsInt+0)==0),ierr);
  EPETRA_TEST_ERR(!(A_Int.IndicesAreGlobal()),ierr);
  EPETRA_TEST_ERR(!(A_Int.FillComplete(false)==0),ierr);
  EPETRA_TEST_ERR(!(A_Int.IndicesAreLocal()),ierr);
  
  delete [] MyGlobalElementsLL;
  delete [] NumNzLL;
  delete [] MyIntGlobalElementsLL;
  delete [] MyGlobalElementsInt; 
  delete [] NumNzInt;
  delete [] MyLLGlobalElementsInt;

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return ierr ;
}
