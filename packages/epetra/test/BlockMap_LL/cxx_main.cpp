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
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#endif
#include "Epetra_SerialComm.h"
#include "checkmap.h"
#include "../epetra_test_err.h"
#include "Epetra_Version.h"

int main(int argc, char *argv[]) {
  bool verbose = false;
  // Check if we should print results to standard out
  if (argc > 1) 
		if ((argv[1][0] == '-') && (argv[1][1] == 'v')) 
			verbose = true;

  int i;
	int ierr = 0;
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
  int NumProc = Comm.NumProc();

  int verbose_int = verbose ? 1 : 0;
  Comm.Broadcast(&verbose_int, 1, 0);
  verbose = verbose_int==1 ? true : false;

  if (verbose && MyPID==0)
    cout << Epetra_Version() << endl << endl;

  if (verbose) cout << Comm << endl << flush;
  Comm.Barrier();
  bool verbose1 = verbose;
  if (verbose) verbose = (MyPID==0);

  int NumMyElements = 10000;
  long long NumGlobalElements = NumMyElements*NumProc+EPETRA_MIN(NumProc,3);
  if (MyPID < 3) NumMyElements++;
  int IndexBase = 0;
  int ElementSize = 7;
  bool DistributedGlobal = (NumGlobalElements>NumMyElements);
  bool IsOneToOne = true;
  Epetra_BlockMap * Map;

  // Test exceptions

  if (verbose) 
    cout << "*******************************************************************************************" << endl
	 << "        Testing Exceptions (Expect error messages if EPETRA_NO_ERROR_REPORTS is not defined" << endl
	 << "*******************************************************************************************" << endl
	 << endl << endl;

  try {
    if (verbose) cout << "Checking Epetra_BlockMap(-2, ElementSize, IndexBase, Comm)" << endl;
    Epetra_BlockMap TestMap(-2LL, ElementSize, IndexBase, Comm);
  }
  catch (int Error) {
    if (Error != -1) {
      if (Error != 0) {
				EPETRA_TEST_ERR(Error,returnierr);
				if (verbose) cout << "Error code should be -1" << endl;
      }
      else { // Error == 0
				cout << "Error code = " << Error << "Should be -1" << endl;
				returnierr += 1;
      }
    }
    else if (verbose) cout << "Checked OK\n\n" << endl;
  }

  try {
    if (verbose) cout << "Checking Epetra_BlockMap(2, 3, ElementSize, IndexBase, Comm)" << endl;
    Epetra_BlockMap TestMap(2LL, 3, ElementSize, IndexBase, Comm);
  }
  catch (int Error) {
    if (Error != -4) {
      if (Error != 0) {
				EPETRA_TEST_ERR(Error,returnierr);
				if (verbose) cout << "Error code should be -4" << endl;
      }
      else { // Error == 0
				cout << "Error code = " << Error << "Should be -4" << endl;
				returnierr += 1;
      }
    }
    else if (verbose) cout << "Checked OK\n\n" << endl;
  }

  if (verbose) cerr << flush;
  if (verbose) cout << flush;
  Comm.Barrier();
  if (verbose) 
    cout << endl << endl
				 << "*******************************************************************************************" << endl
				 << "        Testing valid constructor now......................................................" << endl
				 << "*******************************************************************************************" << endl
				 << endl << endl;
  // Test Epetra-defined uniform linear distribution constructor
  Map = new Epetra_BlockMap(NumGlobalElements, ElementSize, IndexBase, Comm);
  if (verbose) cout << "Checking Epetra_BlockMap(NumGlobalElements, ElementSize, IndexBase, Comm)" << endl;
  ierr = checkmap(*Map, NumGlobalElements, NumMyElements, 0, ElementSize, 0,
		  NumGlobalElements*ElementSize, NumMyElements*ElementSize,
		  IndexBase, Comm, DistributedGlobal,IsOneToOne);

  EPETRA_TEST_ERR(ierr,returnierr);
  if (verbose && ierr==0) cout << "Checked OK\n\n" <<endl;

  delete Map;

  // Test User-defined linear distribution constructor
  Map = new Epetra_BlockMap(NumGlobalElements, NumMyElements, ElementSize, IndexBase, Comm);

  if (verbose) cout << "Checking Epetra_BlockMap(NumGlobalElements, NumMyElements, ElementSize, IndexBase, Comm)" << endl;
  ierr = checkmap(*Map, NumGlobalElements, NumMyElements, 0, ElementSize, 0,
                   NumGlobalElements*ElementSize, NumMyElements*ElementSize,
		  IndexBase, Comm, DistributedGlobal,IsOneToOne);

  EPETRA_TEST_ERR(ierr,returnierr);
  if (verbose && ierr==0) cout << "Checked OK\n\n" <<endl;

  delete Map;

  // Test User-defined arbitrary distribution constructor and fill MyGlobalElements
  // such that the map is not one-to-one.
  int NumMyElems = 5;
  long long NumGlobalElems = (Comm.NumProc()+1)*NumMyElems;
  long long myFirstElem = Comm.MyPID()*NumMyElems;
  if (Comm.MyPID() == 0) NumMyElems *= 2;

  long long* myElems = new long long[NumMyElems];
  for(int ii=0; ii<NumMyElems; ++ii) {
    myElems[ii] = myFirstElem + ii;
  }

  Map = new Epetra_BlockMap(NumGlobalElems, NumMyElems, myElems, 1, 0, Comm);

  if (verbose) cout << "Checking non-oneToOne Epetra_BlockMap(...)"<<endl;
  ierr = Map->IsOneToOne() == false ? 0 : -1;

  //this Map is 1-to-1 if we're running on 1 processor, otherwise it
  //should not be 1-to-1.
  if (Comm.NumProc() > 1) {
    EPETRA_TEST_ERR(ierr,returnierr);
  }
  if (verbose && ierr==0) cout << "Checked OK\n\n" <<endl;

  delete [] myElems;
  delete Map;

  // Test User-defined arbitrary distribution constructor
  // Generate Global Element List.  Do in reverse for fun!

  long long * MyGlobalElements = new long long[NumMyElements];
  int MaxMyGID = (Comm.MyPID()+1)*NumMyElements-1+IndexBase;
  if (Comm.MyPID()>2) 
		MaxMyGID+=3;
  for (i = 0; i<NumMyElements; i++) 
		MyGlobalElements[i] = MaxMyGID-i;

  Map = new Epetra_BlockMap(NumGlobalElements, NumMyElements, MyGlobalElements, ElementSize,
                            IndexBase, Comm);

  if (verbose) cout << "Checking Epetra_BlockMap(NumGlobalElements, NumMyElements, MyGlobalElements,  ElementSize, IndexBase, Comm)" << endl;
  ierr = checkmap(*Map, NumGlobalElements, NumMyElements, MyGlobalElements, ElementSize, 0,
                  NumGlobalElements*ElementSize, NumMyElements*ElementSize,
                  IndexBase, Comm, DistributedGlobal,IsOneToOne);

  EPETRA_TEST_ERR(ierr,returnierr);
  if (verbose && ierr==0) cout << "Checked OK\n\n" <<endl;

  Epetra_BlockMap * Map3 = new Epetra_BlockMap(*Map);// A map to test the SameAs method later

  delete Map;

  int * ElementSizeList = new int[NumMyElements];
  int NumMyEquations = 0;
  int NumGlobalEquations = 0;
  for (i = 0; i<NumMyElements; i++) {
		ElementSizeList[i] = i%6 + 2; // elementsizes go from 2 to 7
		NumMyEquations += ElementSizeList[i];
	}
  ElementSize = 7; // Set to maximum for use in checkmap
  NumGlobalEquations = Comm.NumProc()*NumMyEquations;
	
  // Adjust NumGlobalEquations based on processor ID
  if (Comm.NumProc() > 3) {
		if (Comm.MyPID()>2)
			NumGlobalEquations += 3*((NumMyElements)%6+2);
		else 
			NumGlobalEquations -= (Comm.NumProc()-3)*((NumMyElements-1)%6+2);
	}
  Map = new Epetra_BlockMap(NumGlobalElements, NumMyElements, MyGlobalElements, ElementSizeList,
														IndexBase, Comm);
  if (verbose) cout << "Checking Epetra_BlockMap(NumGlobalElements, NumMyElements, MyGlobalElements,  ElementSizeList, IndexBase, Comm)" << endl;
  ierr = checkmap(*Map, NumGlobalElements, NumMyElements, MyGlobalElements, ElementSize, ElementSizeList,
		  NumGlobalEquations, NumMyEquations,
		  IndexBase, Comm, DistributedGlobal,IsOneToOne);
	
  EPETRA_TEST_ERR(ierr,returnierr);
  if (verbose && ierr==0) cout << "Checked OK\n\n" <<endl;

  // Test Copy constructor
  Epetra_BlockMap * Map1 = new Epetra_BlockMap(*Map);

  // Test SameAs() method
  bool same = Map1->SameAs(*Map);
	EPETRA_TEST_ERR(!(same==true),returnierr);// should return true since Map1 is a copy of Map

  Epetra_BlockMap * Map2 = new Epetra_BlockMap(NumGlobalElements,NumMyElements,MyGlobalElements,ElementSizeList,IndexBase,Comm);
  same = Map2->SameAs(*Map);
	EPETRA_TEST_ERR(!(same==true),returnierr); // Map and Map2 were created with the same sets of parameters
  delete Map2;

  // now test SameAs() on some maps that are different

  Map2 = new Epetra_BlockMap(NumGlobalElements,NumMyElements,MyGlobalElements,ElementSizeList,IndexBase-1,Comm);
  same = Map2->SameAs(*Map);
	EPETRA_TEST_ERR(!(same==false),returnierr); // IndexBases are different
  delete Map2;
	
  int *ElementSizeList1 = new int[NumMyElements];
  for (i=0; i<NumMyElements; i++) 
		ElementSizeList1[i] = i%5 + 2; // element sizes go from 2 to 6
  Map2 = new Epetra_BlockMap(NumGlobalElements,NumMyElements,MyGlobalElements,ElementSizeList1,IndexBase,Comm);
  same = Map2->SameAs(*Map);
	EPETRA_TEST_ERR(!(same==false),returnierr); // ElementSizes are different
  delete [] ElementSizeList1;
  delete Map2;

  same = Map3->SameAs(*Map);
	EPETRA_TEST_ERR(!(same==false),returnierr); // Map3 saved from an earlier test
  delete Map3;

  // Back to testing copy constructor
  if (verbose) cout << "Checking Epetra_BlockMap(*Map)" << endl;
  ierr = checkmap(*Map1, NumGlobalElements, NumMyElements, MyGlobalElements, ElementSize, ElementSizeList,
		  NumGlobalEquations, NumMyEquations,
		  IndexBase, Comm, DistributedGlobal,IsOneToOne);

  EPETRA_TEST_ERR(ierr,returnierr);
  if (verbose && ierr==0) cout << "Checked OK\n\n" <<endl;

  if (verbose1) {
    if (verbose) cout << "Test ostream << operator" << endl << flush;
  }
    // Build a small map for test cout.  Use 10 elements from current map
    long long * MyEls = Map->MyGlobalElements64();
    int * MySz  = Map->ElementSizeList();
    int IndBase = Map->IndexBase();
    int MyLen = EPETRA_MIN(10+Comm.MyPID(),Map->NumMyElements());
    Epetra_BlockMap * SmallMap = new Epetra_BlockMap((long long)-1, MyLen, MyEls, MySz, IndBase, Comm);
    if (verbose1) {
    cout << *SmallMap;
    }
    delete SmallMap;

  delete Map;
  delete Map1;
	

  //create a map where proc 1 has no local elements, then check to make sure that
  //if NumMyElements == 0, then MaxMyGID < MinMyGID.

  if (MyPID == 1) {
    Map1 = new Epetra_BlockMap(-1, 0, (long long*)0, (int*)0, IndexBase, Comm);
  }
  else {
    Map1 = new Epetra_BlockMap(-1, NumMyElements, MyGlobalElements, ElementSizeList, IndexBase, Comm);
  }

  int numMyElems = Map1->NumMyElements();
  if (MyPID == 1) {
    EPETRA_TEST_ERR(!(numMyElems == 0), returnierr);
    long long maxgid = Map1->MaxMyGID64();
    long long mingid = Map1->MinMyGID64();
    EPETRA_TEST_ERR( !(maxgid<mingid), returnierr);
  }

  delete[] ElementSizeList;
  delete[] MyGlobalElements;
  delete Map1;

  // test reference counting
  ierr = 0;

  if (verbose) 
    cout << endl << endl
	 << "*******************************************************************************************" << endl
	 << "        Testing reference counting now....................................................." << endl
	 << "*******************************************************************************************" << endl << endl;

  Epetra_BlockMap b1(NumGlobalElements, NumMyElements, ElementSize, IndexBase, Comm);
  int b1count = b1.ReferenceCount();
  const Epetra_BlockMapData* b1addr = b1.DataPtr();
  EPETRA_TEST_ERR(!(b1count==1),ierr); // count should be 1
  if(verbose) cout << "Default constructor. \nb1= " << b1count << "  " << b1addr << endl;
	
  Epetra_BlockMap* b2 = new Epetra_BlockMap(b1);
  int b2count = b2->ReferenceCount();
  const Epetra_BlockMapData* b2addr = b2->DataPtr();
  int b1countold = b1count;
  b1count = b1.ReferenceCount();
  EPETRA_TEST_ERR(!(b2count==b1count && b1count==(b1countold+1)),ierr); // both counts should be 2
  EPETRA_TEST_ERR(!(b1addr==b2addr),ierr); // addresses should be same
  if(verbose) cout << "Copy constructor. \nb1= " << b1count << "  " << b1addr << "\nb2= " << b2count << "  " << b2addr << endl;

  delete b2;
  b1countold = b1count;
  b1count = b1.ReferenceCount();
  EPETRA_TEST_ERR(!(b1count==b1countold-1), ierr); // count should have decremented (to 1)
  EPETRA_TEST_ERR(!(b1addr==b1.DataPtr()), ierr); // b1addr should be unchanged
  if(verbose) cout << "b2 destroyed. \nb1= " << b1count << "  " << b1addr << endl;

  { // inside of braces to test stack deallocation.
    if(verbose) cout << "Assignment operator, post construction" << endl;
    Epetra_BlockMap b3(NumGlobalElements, NumMyElements, ElementSize, IndexBase-1, Comm);
    int b3count = b3.ReferenceCount();
    const Epetra_BlockMapData* b3addr = b3.DataPtr();
    EPETRA_TEST_ERR(!(b3count==1),ierr); // b3count should be 1 initially
    EPETRA_TEST_ERR(!(b1addr!=b3addr),ierr); // b1 and b3 should have different ptr addresses
    if(verbose) cout << "Prior to assignment: \nb1= " << b1count << "  " << b1addr << "\nb3= " << b3count << "  " << b3addr << endl;
    b3 = b1;
    b3count = b3.ReferenceCount();
    b3addr = b3.DataPtr();
    b1countold = b1count;
    b1count = b1.ReferenceCount();
    EPETRA_TEST_ERR(!(b3count==b1count && b1count==b1countold+1),ierr); // both counts should be 2
    EPETRA_TEST_ERR(!(b1addr==b3addr),ierr); // addresses should be same
    if(verbose) cout << "After assignment: \nb1= " << b1count << "  " << b1addr << "\nb3= " << b3count << "  " << b3addr << endl;
  }
  b1countold = b1count;
  b1count = b1.ReferenceCount();
  EPETRA_TEST_ERR(!(b1count==b1countold-1), ierr); // count should have decremented (to 1)
  EPETRA_TEST_ERR(!(b1addr==b1.DataPtr()), ierr); // b1addr should be unchanged
  if (verbose) cout << "b3 destroyed. \nb1= " << b1count << "  " << b1addr << endl;

  EPETRA_TEST_ERR(ierr,returnierr);
  if (verbose && (ierr == 0)) cout << "Checked OK\n\n" <<endl;
  // done with reference counting testing

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return returnierr;
}
