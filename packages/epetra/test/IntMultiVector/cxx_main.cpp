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


// Epetra_BlockMap Test routine

#include "Epetra_Time.h"
#include "Epetra_BlockMap.h"
#include "Epetra_IntMultiVector.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "ExecuteTestProblems.h"
#include "../epetra_test_err.h"
#include "Epetra_Version.h"

// Restored MultiVector tests
int main(int argc, char *argv[]) {

  int ierr = 0, i, j;

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int rank; // My process ID

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

#else

  int rank = 0;
  Epetra_SerialComm Comm;

#endif

  Comm.SetTracebackMode(0); // This should shut down any error tracing
  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  //  char tmp;
  //  if (rank==0) cout << "Press any key to continue..."<< endl;
  //  if (rank==0) cin >> tmp;
  //  Comm.Barrier();

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  if (verbose && MyPID==0)
    cout << Epetra_Version() << endl << endl;

  if (verbose) cout << Comm <<endl;

  bool verbose1 = verbose;

  // Redefine verbose to only print on PE 0
  if (verbose && rank!=0) verbose = false;

  int NumMyElements = 10000;
  int NumMyElements1 = NumMyElements; // Needed for localmap
  int NumGlobalElements = NumMyElements*NumProc+EPETRA_MIN(NumProc,3);
  if (MyPID < 3) NumMyElements++;
  int IndexBase = 0;
  int ElementSize = 7;
  int NumVectors = 4;

  // Test LocalMap constructor
  // and Epetra-defined uniform linear distribution constructor

  if (verbose) cout << "\n*********************************************************" << endl;
  if (verbose) cout << "Checking Epetra_LocalMap(NumMyElements1, IndexBase, Comm)" << endl;
  if (verbose) cout << "     and Epetra_BlockMap(NumGlobalElements, ElementSize, IndexBase, Comm)" << endl;
  if (verbose) cout << "*********************************************************" << endl;

  Epetra_LocalMap *LocalMap = new Epetra_LocalMap(NumMyElements1, IndexBase,
                              Comm);
  Epetra_BlockMap * BlockMap = new Epetra_BlockMap(NumGlobalElements, ElementSize, IndexBase, Comm);
  EPETRA_TEST_ERR(IntMultiVectorTests(*BlockMap, NumVectors, verbose),ierr);

  delete BlockMap;

  // Test User-defined linear distribution constructor

  if (verbose) cout << "\n*********************************************************" << endl;
  if (verbose) cout << "Checking Epetra_BlockMap(NumGlobalElements, NumMyElements, ElementSize, IndexBase, Comm)" << endl;
  if (verbose) cout << "*********************************************************" << endl;

  BlockMap = new Epetra_BlockMap(NumGlobalElements, NumMyElements, ElementSize, IndexBase, Comm);

  EPETRA_TEST_ERR(IntMultiVectorTests(*BlockMap, NumVectors, verbose),ierr);

  delete BlockMap;

  // Test User-defined arbitrary distribution constructor
  // Generate Global Element List.  Do in reverse for fun!

  int * MyGlobalElements = new int[NumMyElements];
  int MaxMyGID = (Comm.MyPID()+1)*NumMyElements-1+IndexBase;
  if (Comm.MyPID()>2) MaxMyGID+=3;
  for (i = 0; i<NumMyElements; i++) MyGlobalElements[i] = MaxMyGID-i;

  if (verbose) cout << "\n*********************************************************" << endl;
  if (verbose) cout << "Checking Epetra_BlockMap(NumGlobalElements, NumMyElements, MyGlobalElements,  ElementSize, IndexBase, Comm)" << endl;
  if (verbose) cout << "*********************************************************" << endl;

  BlockMap = new Epetra_BlockMap(NumGlobalElements, NumMyElements, MyGlobalElements, ElementSize,
		      IndexBase, Comm);
  EPETRA_TEST_ERR(IntMultiVectorTests(*BlockMap, NumVectors, verbose),ierr);

  delete BlockMap;

  int * ElementSizeList = new int[NumMyElements];
  int NumMyEquations = 0;
  int NumGlobalEquations = 0;
  for (i = 0; i<NumMyElements; i++)
    {
      ElementSizeList[i] = i%6+2; // blocksizes go from 2 to 7
      NumMyEquations += ElementSizeList[i];
    }
  ElementSize = 7; // Set to maximum for use in checkmap
  NumGlobalEquations = Comm.NumProc()*NumMyEquations;

  // Adjust NumGlobalEquations based on processor ID
  if (Comm.NumProc() > 3)
    {
      if (Comm.MyPID()>2)
	NumGlobalEquations += 3*((NumMyElements)%6+2);
      else
	NumGlobalEquations -= (Comm.NumProc()-3)*((NumMyElements-1)%6+2);
    }

  if (verbose) cout << "\n*********************************************************" << endl;
  if (verbose) cout << "Checking Epetra_BlockMap(NumGlobalElements, NumMyElements, MyGlobalElements,  ElementSizeList, IndexBase, Comm)" << endl;
  if (verbose) cout << "*********************************************************" << endl;

  BlockMap = new Epetra_BlockMap(NumGlobalElements, NumMyElements, MyGlobalElements, ElementSizeList,
		      IndexBase, Comm);
  EPETRA_TEST_ERR(IntMultiVectorTests(*BlockMap, NumVectors, verbose),ierr);

  // Test Copy constructor

  if (verbose) cout << "\n*********************************************************" << endl;
  if (verbose) cout << "Checking Epetra_BlockMap(*BlockMap)" << endl;
  if (verbose) cout << "*********************************************************" << endl;

  Epetra_BlockMap * BlockMap1 = new Epetra_BlockMap(*BlockMap);

  EPETRA_TEST_ERR(IntMultiVectorTests(*BlockMap, NumVectors, verbose),ierr);

  delete [] ElementSizeList;
  delete [] MyGlobalElements;
  delete BlockMap;
  delete BlockMap1;


  // Test Petra-defined uniform linear distribution constructor

  if (verbose) cout << "\n*********************************************************" << endl;
  if (verbose) cout << "Checking Epetra_Map(NumGlobalElements, IndexBase, Comm)" << endl;
  if (verbose) cout << "*********************************************************" << endl;

  Epetra_Map * Map = new Epetra_Map(NumGlobalElements, IndexBase, Comm);
  EPETRA_TEST_ERR(IntMultiVectorTests(*Map, NumVectors, verbose),ierr);

  delete Map;

  // Test User-defined linear distribution constructor

  if (verbose) cout << "\n*********************************************************" << endl;
  if (verbose) cout << "Checking Epetra_Map(NumGlobalElements, NumMyElements, IndexBase, Comm)" << endl;
  if (verbose) cout << "*********************************************************" << endl;

  Map = new Epetra_Map(NumGlobalElements, NumMyElements, IndexBase, Comm);

  EPETRA_TEST_ERR(IntMultiVectorTests(*Map, NumVectors, verbose),ierr);

  delete Map;

  // Test User-defined arbitrary distribution constructor
  // Generate Global Element List.  Do in reverse for fun!

  MyGlobalElements = new int[NumMyElements];
  MaxMyGID = (Comm.MyPID()+1)*NumMyElements-1+IndexBase;
  if (Comm.MyPID()>2) MaxMyGID+=3;
  for (i = 0; i<NumMyElements; i++) MyGlobalElements[i] = MaxMyGID-i;

  if (verbose) cout << "\n*********************************************************" << endl;
  if (verbose) cout << "Checking Epetra_Map(NumGlobalElements, NumMyElements, MyGlobalElements,  IndexBase, Comm)" << endl;
  if (verbose) cout << "*********************************************************" << endl;

  Map = new Epetra_Map(NumGlobalElements, NumMyElements, MyGlobalElements,
		      IndexBase, Comm);
  EPETRA_TEST_ERR(IntMultiVectorTests(*Map, NumVectors, verbose),ierr);

  // Test Copy constructor

  if (verbose) cout << "\n*********************************************************" << endl;
  if (verbose) cout << "Checking Epetra_Map(*Map)" << endl;
  if (verbose) cout << "*********************************************************" << endl;

  Epetra_Map Map1(*Map);

  EPETRA_TEST_ERR(IntMultiVectorTests(*Map, NumVectors, verbose),ierr);

  delete [] MyGlobalElements;
  delete Map;
  delete LocalMap;

  if (verbose1)
    {
      // Test MultiVector ostream operator with Petra-defined uniform linear distribution constructor
      // and a small vector

      Epetra_Map Map4(100, IndexBase, Comm);
      int * Ip = new int[200];
      for (j=0; j<2; j++)
	for (i=0; i<100; i++)
	  Ip[i+j*100] = i+j*100;
      Epetra_IntMultiVector I(View, Map4,Ip, 100, 2);

      if (verbose) cout << "\n\nTesting ostream operator:  Multivector  should be 100-by-2 and print i,j indices"
	   << endl << endl;
      cout << I << endl;

      Epetra_BlockMap Map5(-1, 25, 4, IndexBase, Comm);
      Epetra_IntMultiVector I1(View, Map5,Ip, 100, 2);
      if (verbose) cout << "\n\nTesting ostream operator:  Same Multivector as before except using BlockMap of 25x4"
	   << endl << endl;
      cout << I1 << endl;

      if (verbose) cout << "Traceback Mode value = " << I.GetTracebackMode() << endl;
      delete [] Ip;
    }

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return ierr;
}

