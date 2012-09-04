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


// Epetra_Map (and Epetra_LocalMap) Test routine
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "checkmap.h"
#include "../epetra_test_err.h"
#include "Epetra_Version.h"

int checkMapDataClass(Epetra_Comm& Comm, int verbose);
int checkLocalMapDataClass(Epetra_Comm& Comm, int verbose);

int main(int argc, char *argv[]) {

  int ierr=0, returnierr=0;

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;


  if (!verbose) {
    Comm.SetTracebackMode(0); // This should shut down any error traceback reporting
  }
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  if (verbose && MyPID==0)
    cout << Epetra_Version() << endl << endl;

  if (verbose) cout << Comm << endl;

  bool verbose1 = verbose;
  if (verbose) verbose = (MyPID==0);

  int NumMyElements = 10000;
  int NumMyElements1 = NumMyElements; // Used for local map
  long long NumGlobalElements = NumMyElements*NumProc+EPETRA_MIN(NumProc,3);
  if (MyPID < 3) NumMyElements++;
  int IndexBase = 0;
  bool DistributedGlobal = (NumGlobalElements>NumMyElements);
  
  Epetra_Map* Map;

  // Test exceptions

  if (verbose)
    cout << "*******************************************************************************************" << endl
	 << "        Testing Exceptions (Expect error messages if EPETRA_NO_ERROR_REPORTS is not defined" << endl
	 << "*******************************************************************************************" << endl
	 << endl << endl;

  try {
    if (verbose) cout << "Checking Epetra_Map(-2, IndexBase, Comm)" << endl;
    Map = new Epetra_Map((long long)-2, IndexBase, Comm);
  }
  catch (int Error) {
    if (Error!=-1) {
      if (Error!=0) {
	EPETRA_TEST_ERR(Error,returnierr);
	if (verbose) cout << "Error code should be -1" << endl;
      }
      else {
	cout << "Error code = " << Error << "Should be -1" << endl;
	returnierr+=1;
      }
    }
    else if (verbose) cout << "Checked OK\n\n" << endl;
  }

  try {
    if (verbose) cout << "Checking Epetra_Map(2, 3, IndexBase, Comm)" << endl;
    Map = new Epetra_Map((long long)2, 3, IndexBase, Comm);
  }
  catch (int Error) {
    if (Error!=-4) {
      if (Error!=0) {
	EPETRA_TEST_ERR(Error,returnierr);
	if (verbose) cout << "Error code should be -4" << endl;
      }
      else {
	cout << "Error code = " << Error << "Should be -4" << endl;
	returnierr+=1;
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
  Map = new Epetra_Map(NumGlobalElements, IndexBase, Comm);
  if (verbose) cout << "Checking Epetra_Map(NumGlobalElements, IndexBase, Comm)" << endl;
  ierr = checkmap(*Map, NumGlobalElements, NumMyElements, 0, 
		  IndexBase, Comm, DistributedGlobal);

  EPETRA_TEST_ERR(ierr,returnierr);
  if (verbose && ierr==0) cout << "Checked OK\n\n" <<endl;

  delete Map;

  // Test User-defined linear distribution constructor
  Map = new Epetra_Map(NumGlobalElements, NumMyElements, IndexBase, Comm);

  if (verbose) cout << "Checking Epetra_Map(NumGlobalElements, NumMyElements, IndexBase, Comm)" << endl;
  ierr = checkmap(*Map, NumGlobalElements, NumMyElements, 0, 
		  IndexBase, Comm, DistributedGlobal);

  EPETRA_TEST_ERR(ierr,returnierr);
  if (verbose && ierr==0) cout << "Checked OK\n\n" <<endl;
  delete Map;

  // Test User-defined arbitrary distribution constructor
  // Generate Global Element List.  Do in reverse for fun!

  long long * MyGlobalElements = new long long[NumMyElements];
  int MaxMyGID = (Comm.MyPID()+1)*NumMyElements-1+IndexBase;
  if (Comm.MyPID()>2) MaxMyGID+=3;
  for (int i = 0; i<NumMyElements; i++) MyGlobalElements[i] = MaxMyGID-i;

  Map = new Epetra_Map(NumGlobalElements, NumMyElements, MyGlobalElements, 
											 IndexBase, Comm);
  if (verbose) cout << "Checking Epetra_Map(NumGlobalElements, NumMyElements, MyGlobalElements,  IndexBase, Comm)" << endl;
  ierr = checkmap(*Map, NumGlobalElements, NumMyElements, MyGlobalElements, 
									IndexBase, Comm, DistributedGlobal);

  EPETRA_TEST_ERR(ierr,returnierr);
  if (verbose && ierr==0) cout << "Checked OK\n\n" <<endl;
  // Test Copy constructor
  Epetra_Map* Map1 = new Epetra_Map(*Map);

  // Test SameAs() method
  bool same = Map1->SameAs(*Map);
  EPETRA_TEST_ERR(!(same==true),ierr);// should return true since Map1 is a copy of Map

  Epetra_BlockMap* Map2 = new Epetra_Map(NumGlobalElements, NumMyElements, MyGlobalElements,  IndexBase, Comm);
  same = Map2->SameAs(*Map);
  EPETRA_TEST_ERR(!(same==true),ierr); // Map and Map2 were created with the same sets of parameters
  delete Map2;

  // now test SameAs() on a map that is different

  Map2 =  new Epetra_Map(NumGlobalElements, NumMyElements, MyGlobalElements, IndexBase-1, Comm);
  same = Map2->SameAs(*Map);
  EPETRA_TEST_ERR(!(same==false),ierr); // IndexBases are different
  delete Map2;

  // Back to testing copy constructor
  if (verbose) cout << "Checking Epetra_Map(*Map)" << endl;
  ierr = checkmap(*Map1, NumGlobalElements, NumMyElements, MyGlobalElements, 
		  IndexBase, Comm, DistributedGlobal);

  EPETRA_TEST_ERR(ierr,returnierr);
  if (verbose && ierr==0) cout << "Checked OK\n\n" <<endl;
  Epetra_Map* SmallMap = 0;
  if (verbose1) {
    // Build a small map for test cout.  Use 10 elements from current map
    long long* MyEls = Map->MyGlobalElements64();
    int IndBase = Map->IndexBase();
    int MyLen = EPETRA_MIN(10+Comm.MyPID(),Map->NumMyElements());
    SmallMap = new Epetra_Map((long long)-1, MyLen, MyEls, IndBase, Comm);
  }

  delete [] MyGlobalElements;
  delete Map;
  delete Map1;

	// Test reference-counting in Epetra_Map
	if (verbose) cout << "Checking Epetra_Map reference counting" << endl;
	ierr = checkMapDataClass(Comm, verbose);
	EPETRA_TEST_ERR(ierr,returnierr);
  if (verbose && ierr==0) cout << "Checked OK\n\n" <<endl;

  // Test LocalMap constructor
  Epetra_LocalMap* LocalMap = new Epetra_LocalMap((long long)NumMyElements1, IndexBase, Comm);
  if (verbose) cout << "Checking Epetra_LocalMap(NumMyElements1, IndexBase, Comm)" << endl;
  ierr = checkmap(*LocalMap, NumMyElements1, NumMyElements1, 0, IndexBase, Comm, false);

  EPETRA_TEST_ERR(ierr,returnierr);
  if (verbose && ierr==0) cout << "Checked OK\n\n" <<endl;
  // Test Copy constructor
  Epetra_LocalMap* LocalMap1 = new Epetra_LocalMap(*LocalMap);
  if (verbose) cout << "Checking Epetra_LocalMap(*LocalMap)" << endl;
  ierr = checkmap(*LocalMap1, NumMyElements1, NumMyElements1, 0, IndexBase, Comm, false);

  EPETRA_TEST_ERR(ierr,returnierr);
  if (verbose && ierr==0) cout << "Checked OK\n\n" <<endl;
  delete LocalMap1;
  delete LocalMap;

	// Test reference-counting in Epetra_LocalMap
	if (verbose) cout << "Checking Epetra_LocalMap reference counting" << endl;
	ierr = checkLocalMapDataClass(Comm, verbose);
	EPETRA_TEST_ERR(ierr,returnierr);
  if (verbose && ierr==0) cout << "Checked OK\n\n" <<endl;

	// Test output
  if (verbose1) {
    if (verbose) cout << "Test ostream << operator" << endl << flush;
    cout << *SmallMap;
    delete SmallMap;
  }

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return returnierr;
}

int checkMapDataClass(Epetra_Comm& Comm, int verbose) {
	int returnierr = 0;
	long long NumGlobalElements = 1000;
	int IndexBase = 0;

	Epetra_Map m1(NumGlobalElements, IndexBase, Comm);
	int m1count = m1.ReferenceCount();
	const Epetra_BlockMapData* m1addr = m1.DataPtr();
	EPETRA_TEST_ERR(!(m1count==1),returnierr); // count should be 1
	if(verbose) cout << "Default constructor. \nm1= " << m1count << "  " << m1addr << endl;
	
	Epetra_Map* m2 = new Epetra_Map(m1);
	int m2count = m2->ReferenceCount();
	const Epetra_BlockMapData* m2addr = m2->DataPtr();
	int m1countold = m1count;
	m1count = m1.ReferenceCount();
	EPETRA_TEST_ERR(!(m2count==m1count && m1count==(m1countold+1)),returnierr); // both counts should be 2
	EPETRA_TEST_ERR(!(m1addr==m2addr),returnierr); // addresses should be same
	if(verbose) cout << "Copy constructor. \nm1= " << m1count << "  " << m1addr 
									 << "\nm2= " << m2count << "  " << m2addr << endl;

	delete m2;
	m1countold = m1count;
	m1count = m1.ReferenceCount();
	EPETRA_TEST_ERR(!(m1count == m1countold-1), returnierr); // count should have decremented (to 1)
	EPETRA_TEST_ERR(!(m1addr == m1.DataPtr()),returnierr); // m1addr should be unchanged
	if(verbose) cout << "m2 destroyed. \nm1= " << m1count << "  " << m1addr << endl;

	{ // inside of braces to test stack deallocation.
		if(verbose) cout << "Assignment operator, post construction" << endl;
		Epetra_Map m3(NumGlobalElements, IndexBase+1, Comm);
		int m3count = m3.ReferenceCount();
		const Epetra_BlockMapData* m3addr = m3.DataPtr();
		EPETRA_TEST_ERR(!(m3count==1),returnierr); // m3count should be 1 initially
		EPETRA_TEST_ERR(!(m1addr!=m3addr),returnierr); // m1 and m3 should have different ptr addresses
		if(verbose) cout << "Prior to assignment: \nm1= " << m1count << "  " << m1addr 
										 << "\nm3= " << m3count << "  " << m3addr << endl;
		m3 = m1;
		m3count = m3.ReferenceCount();
		m3addr = m3.DataPtr();
		m1countold = m1count;
		m1count = m1.ReferenceCount();
		EPETRA_TEST_ERR(!(m3count==m1count && m1count==m1countold+1),returnierr); // both counts should be 2
		EPETRA_TEST_ERR(!(m1addr==m3addr),returnierr); // addresses should be same
		if(verbose) cout << "After assignment: \nm1= " << m1count << "  " << m1addr 
										 << "\nm3= " << m3count << "  " << m3addr << endl;
	}
	m1countold = m1count;
	m1count = m1.ReferenceCount();
	EPETRA_TEST_ERR(!(m1count==m1countold-1), returnierr); // count should have decremented (to 1)
	EPETRA_TEST_ERR(!(m1addr== m1.DataPtr()),returnierr); // m1addr should be unchanged
	if(verbose) cout << "m3 destroyed. \nm1= " << m1count << "  " << m1addr << endl;

	return(returnierr);
}

int checkLocalMapDataClass(Epetra_Comm& Comm, int verbose) {
	int returnierr = 0;
	int NumMyElements = 100;
	int IndexBase = 0;

	Epetra_LocalMap m1(NumMyElements, IndexBase, Comm);
	int m1count = m1.ReferenceCount();
	const Epetra_BlockMapData* m1addr = m1.DataPtr();
	EPETRA_TEST_ERR(!(m1count==1),returnierr); // count should be 1
	if(verbose) cout << "Default constructor. \nm1= " << m1count << "  " << m1addr << endl;
	
	Epetra_LocalMap* m2 = new Epetra_LocalMap(m1);
	int m2count = m2->ReferenceCount();
	const Epetra_BlockMapData* m2addr = m2->DataPtr();
	int m1countold = m1count;
	m1count = m1.ReferenceCount();
	EPETRA_TEST_ERR(!(m2count==m1count && m1count==(m1countold+1)),returnierr); // both counts should be 2
	EPETRA_TEST_ERR(!(m1addr==m2addr),returnierr); // addresses should be same
	if(verbose) cout << "Copy constructor. \nm1= " << m1count << "  " << m1addr 
									 << "\nm2= " << m2count << "  " << m2addr << endl;

	delete m2;
	m1countold = m1count;
	m1count = m1.ReferenceCount();
	EPETRA_TEST_ERR(!(m1count == m1countold-1), returnierr); // count should have decremented (to 1)
	EPETRA_TEST_ERR(!(m1addr == m1.DataPtr()),returnierr); // m1addr should be unchanged
	if(verbose) cout << "m2 destroyed. \nm1= " << m1count << "  " << m1addr << endl;

	{ // inside of braces to test stack deallocation.
		if(verbose) cout << "Assignment operator, post construction" << endl;
		Epetra_LocalMap m3(NumMyElements, IndexBase+1, Comm);
		int m3count = m3.ReferenceCount();
		const Epetra_BlockMapData* m3addr = m3.DataPtr();
		EPETRA_TEST_ERR(!(m3count==1),returnierr); // m3count should be 1 initially
		EPETRA_TEST_ERR(!(m1addr!=m3addr),returnierr); // m1 and m3 should have different ptr addresses
		if(verbose) cout << "Prior to assignment: \nm1= " << m1count << "  " << m1addr 
										 << "\nm3= " << m3count << "  " << m3addr << endl;
		m3 = m1;
		m3count = m3.ReferenceCount();
		m3addr = m3.DataPtr(); // cast int* to int
		m1countold = m1count;
		m1count = m1.ReferenceCount();
		EPETRA_TEST_ERR(!(m3count==m1count && m1count==m1countold+1),returnierr); // both counts should be 2
		EPETRA_TEST_ERR(!(m1addr==m3addr),returnierr); // addresses should be same
		if(verbose) cout << "After assignment: \nm1= " << m1count << "  " << m1addr 
										 << "\nm3= " << m3count << "  " << m3addr << endl;
	}
	m1countold = m1count;
	m1count = m1.ReferenceCount();
	EPETRA_TEST_ERR(!(m1count==m1countold-1), returnierr); // count should have decremented (to 1)
	EPETRA_TEST_ERR(!(m1addr== m1.DataPtr()),returnierr); // m1addr should be unchanged
	if(verbose) cout << "m3 destroyed. \nm1= " << m1count << "  " << m1addr << endl;

	return(returnierr);
}
