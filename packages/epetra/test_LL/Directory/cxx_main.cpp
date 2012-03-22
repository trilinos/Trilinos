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
#include "Epetra_Map.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#endif
#include "Epetra_SerialComm.h"
#include "Epetra_Util.h"
#include "../epetra_test_err.h"
#include "Epetra_Version.h"
#include "Epetra_Directory.h"

int directory_test_1(Epetra_Comm& Comm);
int directory_test_2(Epetra_Comm& Comm);
int directory_test_3(Epetra_Comm& Comm);
int directory_test_4(Epetra_Comm& Comm);
int directory_test_5(Epetra_Comm& Comm);

int main(int argc, char *argv[]) {
  bool verbose = false;
  // Check if we should print results to standard out
  if (argc > 1) {
    if ((argv[1][0] == '-') && (argv[1][1] == 'v')) {
      verbose = true;
    }
  }

  int returnierr = 0;

#ifdef EPETRA_MPI

  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  if (!verbose) {
    Comm.SetTracebackMode(0); // This should shut down any error traceback reporting
  }
  int MyPID = Comm.MyPID();

  int verbose_int = verbose ? 1 : 0;
  Comm.Broadcast(&verbose_int, 1, 0);
  verbose = verbose_int==1 ? true : false;

  if (verbose && MyPID==0)
    cout << Epetra_Version() << endl << endl;

  EPETRA_TEST_ERR( directory_test_1(Comm), returnierr );

  EPETRA_TEST_ERR( directory_test_2(Comm), returnierr );

  EPETRA_TEST_ERR( directory_test_3(Comm), returnierr );

  EPETRA_TEST_ERR( directory_test_4(Comm), returnierr );

  EPETRA_TEST_ERR( directory_test_5(Comm), returnierr );

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  if (MyPID == 0) {
    if (returnierr == 0) {
      cout << "Epetra_Directory tests passed."<<endl;
    }
    else {
      cout << "Epetra_Directory tests failed."<<endl;
    }
  }

  return returnierr;
}

int directory_test_1(Epetra_Comm& Comm)
{
  //set up a map with arbitrary distribution of IDs, but with unique
  //processor ID ownership (i.e., each ID only appears on 1 processor)

  int myPID = Comm.MyPID();
  int numProcs = Comm.NumProc();

  if (numProcs < 2) return(0);

  int myFirstID = (myPID+1)*(myPID+1);
  int myNumIDs = 3+myPID;

  int* myIDs = new int[myNumIDs];
  int i;
  for(i=0; i<myNumIDs; ++i) {
    myIDs[i] = myFirstID+i;
  }

  Epetra_BlockMap blkmap(-1, myNumIDs, myIDs, 1, 0, Comm);

  Epetra_Directory* directory = Comm.CreateDirectory(blkmap);

  int proc = myPID+1;
  if (proc >= numProcs) proc = 0;

  int procNumIDs = 3+proc;
  int procFirstID = (proc+1)*(proc+1);
  int procLastID = procFirstID+procNumIDs - 1;

  int queryProc1 = -1;
  int queryProc2 = -1;

  int err = directory->GetDirectoryEntries(blkmap, 1, &procFirstID,
					   &queryProc1, NULL, NULL);
  err += directory->GetDirectoryEntries(blkmap, 1, &procLastID,
					&queryProc2, NULL, NULL);
  delete directory;
  delete [] myIDs;

  if (queryProc1 != proc || queryProc2 != proc) {
    return(-1);
  }

  return(0);
}

int directory_test_2(Epetra_Comm& Comm)
{
  //set up a Epetra_BlockMap with arbitrary distribution of IDs, but with unique
  //processor ID ownership (i.e., each ID only appears on 1 processor)
  //
  //the thing that makes this Epetra_BlockMap nasty is that higher-numbered
  //processors own lower IDs.

  int myPID = Comm.MyPID();
  int numProcs = Comm.NumProc();

  if (numProcs < 2) return(0);

  int myFirstID = (numProcs-myPID)*(numProcs-myPID);
  int myNumIDs = 3;

  int* myIDs = new int[myNumIDs];
  int i;
  for(i=0; i<myNumIDs; ++i) {
    myIDs[i] = myFirstID+i;
  }

  Epetra_BlockMap blkmap(-1, myNumIDs, myIDs, 1, 0, Comm);

  Epetra_Directory* directory = Comm.CreateDirectory(blkmap);

  int proc = myPID+1;
  if (proc >= numProcs) proc = 0;

  int procNumIDs = 3;
  int procFirstID = (numProcs-proc)*(numProcs-proc);
  int procLastID = procFirstID+procNumIDs - 1;

  int queryProc1 = -1;
  int queryProc2 = -1;

  int err = directory->GetDirectoryEntries(blkmap, 1, &procFirstID,
					   &queryProc1, NULL, NULL);
  err += directory->GetDirectoryEntries(blkmap, 1, &procLastID,
					&queryProc2, NULL, NULL);
  delete directory;
  delete [] myIDs;

  if (queryProc1 != proc || queryProc2 != proc) {
    return(-1);
  }

  return(0);
}

int directory_test_3(Epetra_Comm& Comm)
{
  //set up a map with arbitrary distribution of IDs, including non-unique
  //processor ID ownership (i.e., some IDs appear on more than 1 processor)

  int myPID = Comm.MyPID();
  int numProcs = Comm.NumProc();

  if (numProcs < 2) return(0);

  int myFirstID = (myPID+1)*(myPID+1);
  int myNumIDs = 4;

  int* myIDs = new int[myNumIDs];
  int i;
  for(i=0; i<myNumIDs-1; ++i) {
    myIDs[i] = myFirstID+i;
  }

  int nextProc = myPID+1;
  if (nextProc >= numProcs) nextProc = 0;

  int nextProcFirstID = (nextProc+1)*(nextProc+1);
  myIDs[myNumIDs-1] = nextProcFirstID;

  Epetra_BlockMap blkmap(-1, myNumIDs, myIDs, 1, 0, Comm);

  Epetra_Directory* directory = Comm.CreateDirectory(blkmap);

  bool uniqueGIDs = directory->GIDsAllUniquelyOwned();

  delete directory;
  delete [] myIDs;

  if (uniqueGIDs) {
    return(-1);
  }

  return(0);
}

int directory_test_4(Epetra_Comm& Comm)
{
  int myPID = Comm.MyPID();
  int numProcs = Comm.NumProc();

  if (numProcs < 2) return(0);

  //Set up a map with overlapping ranges of GIDs.
  int num = 5;
  int numMyGIDs = 2*num;
  int myFirstGID = myPID*num;

  int* myGIDs = new int[numMyGIDs];

  for(int i=0; i<numMyGIDs; ++i) {
    myGIDs[i] = myFirstGID+i;
  }

  Epetra_Map overlappingmap((long long)-1, numMyGIDs, myGIDs, 0, Comm);

  delete [] myGIDs;

  long long numGlobal0 = overlappingmap.NumGlobalElements();

  Epetra_Map uniquemap1 =
    Epetra_Util::Create_OneToOne_Map(overlappingmap);

  bool use_high_sharing_proc = true;

  Epetra_Map uniquemap2 =
    Epetra_Util::Create_OneToOne_Map(overlappingmap, use_high_sharing_proc);

  long long numGlobal1 = uniquemap1.NumGlobalElements();
  long long numGlobal2 = uniquemap2.NumGlobalElements();

  //The two one-to-one maps should have the same number of global elems.
  if (numGlobal1 != numGlobal2) {
    return(-1);
  }

  //The number of global elems should be greater in the original map
  //than in the one-to-one map.
  if (numGlobal0 <= numGlobal1) {
    return(-2);
  }

  int numLocal1 = uniquemap1.NumMyElements();
  int numLocal2 = uniquemap2.NumMyElements();

  //If this is proc 0 or proc numProcs-1, then the number of
  //local elements should be different in the two one-to-one maps.
  if ((myPID==0 || myPID==numProcs-1) && numLocal1 == numLocal2) {
    return(-3);
  }

  return(0);
}

int directory_test_5(Epetra_Comm& Comm)
{
  int myPID = Comm.MyPID();
  int numProcs = Comm.NumProc();

  if (numProcs < 2) return(0);

  //Set up a map with overlapping ranges of GIDs.
  int num = 5;
  int numMyGIDs = 2*num;
  int myFirstGID = myPID*num;

  int* myGIDs = new int[numMyGIDs];
  int* sizes = new int[numMyGIDs];

  for(int i=0; i<numMyGIDs; ++i) {
    myGIDs[i] = myFirstGID+i;
    sizes[i] = myFirstGID+i+1;
  }

  Epetra_BlockMap overlappingmap(-1, numMyGIDs, myGIDs, sizes, 0, Comm);

  delete [] myGIDs;
  delete [] sizes;

  long long numGlobal0 = overlappingmap.NumGlobalElements();

  Epetra_BlockMap uniquemap1 =
    Epetra_Util::Create_OneToOne_BlockMap(overlappingmap);

  bool use_high_sharing_proc = true;

  Epetra_BlockMap uniquemap2 =
    Epetra_Util::Create_OneToOne_BlockMap(overlappingmap, use_high_sharing_proc);

  long long numGlobal1 = uniquemap1.NumGlobalElements();
  long long numGlobal2 = uniquemap2.NumGlobalElements();

  //The two one-to-one maps should have the same number of global elems.
  if (numGlobal1 != numGlobal2) {
    return(-1);
  }

  //The number of global elems should be greater in the original map
  //than in the one-to-one map.
  if (numGlobal0 <= numGlobal1) {
    return(-2);
  }

  int numLocal1 = uniquemap1.NumMyElements();
  int numLocal2 = uniquemap2.NumMyElements();

  //If this is proc 0 or proc numProcs-1, then the number of
  //local elements should be different in the two one-to-one maps.
  if ((myPID==0 || myPID==numProcs-1) && numLocal1 == numLocal2) {
    return(-3);
  }

  return(0);
}
