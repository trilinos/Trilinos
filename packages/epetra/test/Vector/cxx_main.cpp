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

// Epetra_Vector Test routine

#include "Epetra_Time.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "BuildTestProblems.h"
#include "ExecuteTestProblems.h"
#include "../epetra_test_err.h"
#include "Epetra_Version.h"

int main(int argc, char *argv[]) {

  int ierr = 0, i;

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

  Comm.SetTracebackMode(0); // This should shut down any error tracing
  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  //  char tmp;
  //  if (rank==0) cout << "Press any key to continue..."<< endl;
  //  if (rank==0) cin >> tmp;
  //  Comm.Barrier();

  Comm.SetTracebackMode(0); // This should shut down any error traceback reporting
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
  
  // Test LocalMap constructor
  // and Petra-defined uniform linear distribution constructor

  if (verbose) cout << "\n*********************************************************" << endl;
  if (verbose) cout << "Checking Epetra_LocalMap(NumMyElements1, IndexBase, Comm)" << endl;
  if (verbose) cout << "     and Epetra_BlockMap(NumGlobalElements, ElementSize, IndexBase, Comm)" << endl;
  if (verbose) cout << "*********************************************************" << endl;

  Epetra_LocalMap *LocalMap = new Epetra_LocalMap(NumMyElements1, IndexBase,
                              Comm);
  Epetra_BlockMap * BlockMap = new Epetra_BlockMap(NumGlobalElements, ElementSize, IndexBase, Comm);
  EPETRA_TEST_ERR(VectorTests(*BlockMap, verbose),ierr);

  EPETRA_TEST_ERR(MatrixTests(*BlockMap, *LocalMap, verbose),ierr);

  delete BlockMap;

  // Test User-defined linear distribution constructor

  if (verbose) cout << "\n*********************************************************" << endl;
  if (verbose) cout << "Checking Epetra_BlockMap(NumGlobalElements, NumMyElements, ElementSize, IndexBase, Comm)" << endl;
  if (verbose) cout << "*********************************************************" << endl;

  BlockMap = new Epetra_BlockMap(NumGlobalElements, NumMyElements, ElementSize, IndexBase, Comm);

  EPETRA_TEST_ERR(VectorTests(*BlockMap, verbose),ierr);

  EPETRA_TEST_ERR(MatrixTests(*BlockMap, *LocalMap, verbose),ierr);

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
  EPETRA_TEST_ERR(VectorTests(*BlockMap, verbose),ierr);

  EPETRA_TEST_ERR(MatrixTests(*BlockMap, *LocalMap, verbose),ierr);

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
  EPETRA_TEST_ERR(VectorTests(*BlockMap, verbose),ierr);

  EPETRA_TEST_ERR(MatrixTests(*BlockMap, *LocalMap, verbose),ierr);

  // Test Copy constructor

  if (verbose) cout << "\n*********************************************************" << endl;
  if (verbose) cout << "Checking Epetra_BlockMap(*BlockMap)" << endl;
  if (verbose) cout << "*********************************************************" << endl;

  Epetra_BlockMap * BlockMap1 = new Epetra_BlockMap(*BlockMap);

  EPETRA_TEST_ERR(VectorTests(*BlockMap, verbose),ierr);

  EPETRA_TEST_ERR(MatrixTests(*BlockMap, *LocalMap, verbose),ierr);

  delete [] ElementSizeList;
  delete [] MyGlobalElements;
  delete BlockMap;
  delete BlockMap1;


  // Test Petra-defined uniform linear distribution constructor

  if (verbose) cout << "\n*********************************************************" << endl;
  if (verbose) cout << "Checking Epetra_Map(NumGlobalElements, IndexBase, Comm)" << endl;
  if (verbose) cout << "*********************************************************" << endl;

  Epetra_Map * Map = new Epetra_Map(NumGlobalElements, IndexBase, Comm);
  EPETRA_TEST_ERR(VectorTests(*Map, verbose),ierr);

  EPETRA_TEST_ERR(MatrixTests(*Map, *LocalMap, verbose),ierr);

  delete Map;

  // Test User-defined linear distribution constructor

  if (verbose) cout << "\n*********************************************************" << endl;
  if (verbose) cout << "Checking Epetra_Map(NumGlobalElements, NumMyElements, IndexBase, Comm)" << endl;
  if (verbose) cout << "*********************************************************" << endl;

  Map = new Epetra_Map(NumGlobalElements, NumMyElements, IndexBase, Comm);

  EPETRA_TEST_ERR(VectorTests(*Map, verbose),ierr);

  EPETRA_TEST_ERR(MatrixTests(*Map, *LocalMap, verbose),ierr);

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
  EPETRA_TEST_ERR(VectorTests(*Map, verbose),ierr);

  EPETRA_TEST_ERR(MatrixTests(*Map, *LocalMap, verbose),ierr);

  // Test Copy constructor

  if (verbose) cout << "\n*********************************************************" << endl;
  if (verbose) cout << "Checking Epetra_Map(*Map)" << endl;
  if (verbose) cout << "*********************************************************" << endl;
 
  Epetra_Map Map1(*Map);

  EPETRA_TEST_ERR(VectorTests(*Map, verbose),ierr);

  EPETRA_TEST_ERR(MatrixTests(*Map, *LocalMap, verbose),ierr);

  delete [] MyGlobalElements;
  delete Map;

  if (verbose1)
    {
      // Test Vector MFLOPS for 2D Dot Product
      int M = 1;
      int K = 1000000;
      Epetra_Map Map2(-1, K, IndexBase, Comm);
      Epetra_LocalMap Map3(M, IndexBase, Comm);
      
      Epetra_Vector A(Map2);A.Random();
      Epetra_Vector B(Map2);B.Random();
      Epetra_Vector C(Map3);C.Random();

      // Test Epetra_Vector label
      char* VecLabel = A.Label();
      const char* VecLabel1 = "Epetra::Vector";
      if (verbose) cout << endl << endl <<"This should say " << VecLabel1 << ": " << VecLabel << endl << endl << endl;
      EPETRA_TEST_ERR(strcmp(VecLabel1,VecLabel),ierr);
      if (verbose) cout << "Testing Assignment operator" << endl;

      double tmp1 = 1.00001* (double) (MyPID+1);
      double tmp2 = tmp1;
      A[1] = tmp1;
      tmp2 = A[1];
      cout << "On PE "<< MyPID << "  A[1] should equal = " << tmp1;
      if (tmp1==tmp2) cout << " and it does!" << endl;
      else cout << " but it equals " << tmp2;
 
      Comm.Barrier();
	  
      if (verbose) cout << endl << endl << "Testing MFLOPs" << endl;
      Epetra_Flops counter;
      C.SetFlopCounter(counter);
      Epetra_Time mytimer(Comm);
      C.Multiply('T', 'N', 0.5, A, B, 0.0);
      double Multiply_time = mytimer.ElapsedTime();
      double Multiply_flops = C.Flops();
      if (verbose) cout << "\n\nTotal FLOPs = " << Multiply_flops << endl;
      if (verbose) cout << "Total Time  = " << Multiply_time << endl;
      if (verbose) cout << "MFLOPs      = " << Multiply_flops/Multiply_time/1000000.0 << endl;

      Comm.Barrier(); 
	  
      // Test Vector ostream operator with Petra-defined uniform linear distribution constructor
      // and a small vector
      
      Epetra_Map Map4(100, IndexBase, Comm);
      double * Dp = new double[100]; 
      for (i=0; i<100; i++)
	Dp[i] = i;
      Epetra_Vector D(View, Map4,Dp);
	  
      if (verbose) cout << "\n\nTesting ostream operator:  Multivector  should be 100-by-2 and print i,j indices" 
	   << endl << endl;
      cout << D << endl;

      if (verbose) cout << "Traceback Mode value = " << D.GetTracebackMode() << endl;
      delete [] Dp;
    }

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return ierr;
}

