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

#include "Epetra_CrsGraph.h"
#include "Epetra_Map.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "../epetra_test_err.h"
#include "Epetra_Version.h"

// Prototype
int check(Epetra_CrsGraph& A, int NumMyRows1, int NumGlobalRows1, int NumMyNonzeros1,
					int NumGlobalNonzeros1, int* MyGlobalElements, bool verbose);

int checkSharedOwnership(Epetra_Comm& Comm, bool verbose);
int checkCopyAndAssignment(Epetra_Comm& Comm, bool verbose);

int main(int argc, char *argv[]) {
  int ierr = 0;
  int i;
  int forierr = 0;
  int* Indices;
  bool debug = true;

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );

#else

  int size = 1; // Serial case (not using MPI)
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

  //char tmp;
  //if (rank==0) cout << "Press any key to continue..."<< endl;
  //if (rank==0) cin >> tmp;
  //Comm.Barrier();

  Comm.SetTracebackMode(0); // This should shut down any error traceback reporting
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  if(verbose && MyPID==0)
    cout << Epetra_Version() << endl << endl;

  if(verbose) cout << "Processor "<<MyPID<<" of "<< NumProc << " is alive." << endl;

  bool verbose1 = verbose;

  // Redefine verbose to only print on PE 0
  if(verbose && rank != 0) 
		verbose = false;

  int NumMyEquations = 5;
  int NumGlobalEquations = NumMyEquations*NumProc+EPETRA_MIN(NumProc,3);
  if(MyPID < 3) 
    NumMyEquations++;

  // Construct a Map that puts approximately the same Number of equations on each processor

  Epetra_Map& Map = *new Epetra_Map(NumGlobalEquations, NumMyEquations, 0, Comm);
  
  // Get update list and number of local equations from newly created Map
  int* MyGlobalElements = new int[Map.NumMyElements()];
  Map.MyGlobalElements(MyGlobalElements);

  // Create an integer vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation on this processor

  int* NumNz = new int[NumMyEquations];

  // We are building a tridiagonal matrix where each row has (-1 2 -1)
  // So we need 2 off-diagonal terms (except for the first and last equation)

  for(i=0; i<NumMyEquations; i++)
    if(MyGlobalElements[i]==0 || MyGlobalElements[i] == NumGlobalEquations-1)
      NumNz[i] = 1;
    else
      NumNz[i] = 2;

  // Create a Epetra_CrsGraph

  Epetra_CrsGraph& A = *new Epetra_CrsGraph(Copy, Map, NumNz);
  EPETRA_TEST_ERR(A.IndicesAreGlobal(),ierr);
  EPETRA_TEST_ERR(A.IndicesAreLocal(),ierr);
  
  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1

  Indices = new int[2];
  int NumEntries;

  forierr = 0;
  for(i = 0; i < NumMyEquations; i++) {
    if(MyGlobalElements[i] == 0) {
			Indices[0] = 1;
			NumEntries = 1;
		}
    else if(MyGlobalElements[i] == NumGlobalEquations-1) {
			Indices[0] = NumGlobalEquations-2;
			NumEntries = 1;
		}
    else {
			Indices[0] = MyGlobalElements[i]-1;
			Indices[1] = MyGlobalElements[i]+1;
			NumEntries = 2;
		}
		forierr += !(A.InsertGlobalIndices(MyGlobalElements[i], NumEntries, Indices)==0);
		forierr += !(A.InsertGlobalIndices(MyGlobalElements[i], 1, MyGlobalElements+i)>0); // Put in the diagonal entry (should cause realloc)
  }
  EPETRA_TEST_ERR(forierr,ierr);
	//A.PrintGraphData(cout);
  delete[] Indices;
  
  // Finish up
  EPETRA_TEST_ERR(!(A.IndicesAreGlobal()),ierr);
  EPETRA_TEST_ERR(!(A.FillComplete()==0),ierr);
  EPETRA_TEST_ERR(!(A.IndicesAreLocal()),ierr);
  EPETRA_TEST_ERR(A.StorageOptimized(),ierr);

	A.OptimizeStorage();

  EPETRA_TEST_ERR(!(A.StorageOptimized()),ierr);
  EPETRA_TEST_ERR(A.UpperTriangular(),ierr);
  EPETRA_TEST_ERR(A.LowerTriangular(),ierr);

  if(verbose) cout << "\n*****Testing variable entry constructor\n" << endl;

  int NumMyNonzeros = 3 * NumMyEquations;
  if(A.LRID(0) >= 0) 
		NumMyNonzeros--; // If I own first global row, then there is one less nonzero
  if(A.LRID(NumGlobalEquations-1) >= 0) 
		NumMyNonzeros--; // If I own last global row, then there is one less nonzero

  EPETRA_TEST_ERR(check(A, NumMyEquations, NumGlobalEquations, NumMyNonzeros, 3*NumGlobalEquations-2, 
												MyGlobalElements, verbose),ierr);
  forierr = 0;
  for(i = 0; i < NumMyEquations; i++) 
		forierr += !(A.NumGlobalIndices(MyGlobalElements[i])==NumNz[i]+1);
  EPETRA_TEST_ERR(forierr,ierr);
  for(i = 0; i < NumMyEquations; i++) 
		forierr += !(A.NumMyIndices(i)==NumNz[i]+1);
  EPETRA_TEST_ERR(forierr,ierr);

  if(verbose) cout << "NumIndices function check OK" << endl;

  delete &A;

  if(debug) Comm.Barrier();

  if(verbose) cout << "\n*****Testing constant entry constructor\n" << endl;

  Epetra_CrsGraph& AA = *new Epetra_CrsGraph(Copy, Map, 5);
  
  if(debug) Comm.Barrier();

  for(i = 0; i < NumMyEquations; i++) 
		AA.InsertGlobalIndices(MyGlobalElements[i], 1, MyGlobalElements+i);

  // Note:  All processors will call the following Insert routines, but only the processor
  //        that owns it will actually do anything

  int One = 1;
  if(AA.MyGlobalRow(0)) {
    EPETRA_TEST_ERR(!(AA.InsertGlobalIndices(0, 0, &One)==0),ierr);
  }
  else 
		EPETRA_TEST_ERR(!(AA.InsertGlobalIndices(0, 1, &One)==-1),ierr);
  EPETRA_TEST_ERR(!(AA.FillComplete()==0),ierr);
  EPETRA_TEST_ERR(AA.StorageOptimized(),ierr);
  EPETRA_TEST_ERR(!(AA.UpperTriangular()),ierr);
  EPETRA_TEST_ERR(!(AA.LowerTriangular()),ierr);

  if(debug) Comm.Barrier();
  EPETRA_TEST_ERR(check(AA, NumMyEquations, NumGlobalEquations, NumMyEquations, NumGlobalEquations, 
												MyGlobalElements, verbose),ierr);

  if(debug) Comm.Barrier();

  forierr = 0;
  for(i = 0; i < NumMyEquations; i++) 
		forierr += !(AA.NumGlobalIndices(MyGlobalElements[i])==1);
  EPETRA_TEST_ERR(forierr,ierr);

  if(verbose) cout << "NumIndices function check OK" << endl;

  if(debug) Comm.Barrier();

  if(verbose) cout << "\n*****Testing copy constructor\n" << endl;

  Epetra_CrsGraph& B = *new Epetra_CrsGraph(AA);
  delete &AA;

  EPETRA_TEST_ERR(check(B, NumMyEquations, NumGlobalEquations, NumMyEquations, NumGlobalEquations, 
												MyGlobalElements, verbose),ierr);

  forierr = 0;
  for(i = 0; i < NumMyEquations; i++) 
		forierr += !(B.NumGlobalIndices(MyGlobalElements[i])==1);
  EPETRA_TEST_ERR(forierr,ierr);

  if(verbose) cout << "NumIndices function check OK" << endl;

  if(debug) Comm.Barrier();

  if(verbose) cout << "\n*****Testing post construction modifications\n" << endl;

  EPETRA_TEST_ERR(!(B.InsertGlobalIndices(0, 1, &One)==-2),ierr);
  delete &B;

  // Release all objects
  delete[] MyGlobalElements;
  delete[] NumNz;
  delete &Map;
			

  if (verbose1) {
    // Test ostream << operator (if verbose1)
    // Construct a Map that puts 2 equations on each PE
    
    int NumMyElements1 = 4;
    int NumMyEquations1 = NumMyElements1;
    int NumGlobalEquations1 = NumMyEquations1*NumProc;

    Epetra_Map& Map1 = *new Epetra_Map(NumGlobalEquations1, NumMyElements1, 1, Comm);
    
    // Get update list and number of local equations from newly created Map
    int* MyGlobalElements1 = new int[Map1.NumMyElements()];
    Map1.MyGlobalElements(MyGlobalElements1);
    
    // Create an integer vector NumNz that is used to build the Petra Matrix.
    // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation on this processor
    
    int* NumNz1 = new int[NumMyEquations1];
    
    // We are building a tridiagonal matrix where each row has (-1 2 -1)
    // So we need 2 off-diagonal terms (except for the first and last equation)
    
    for(i = 0; i < NumMyEquations1; i++)
      if(MyGlobalElements1[i]==1 || MyGlobalElements1[i] == NumGlobalEquations1)
				NumNz1[i] = 1;
      else
				NumNz1[i] = 2;
    
    // Create a Epetra_Graph using 1-based arithmetic
    
    Epetra_CrsGraph& A1 = *new Epetra_CrsGraph(Copy, Map1, NumNz1);
    
    // Add  rows one-at-a-time
    // Need some vectors to help
    // Off diagonal Values will always be -1
    
    
    int* Indices1 = new int[2];
    int NumEntries1;

    forierr = 0;
    for(i = 0; i < NumMyEquations1; i++) {
			if(MyGlobalElements1[i]==1) {
				Indices1[0] = 2;
				NumEntries1 = 1;
			}
			else if(MyGlobalElements1[i] == NumGlobalEquations1) {
	    Indices1[0] = NumGlobalEquations1-1;
	    NumEntries1 = 1;
			}
			else {
				Indices1[0] = MyGlobalElements1[i]-1;
				Indices1[1] = MyGlobalElements1[i]+1;
				NumEntries1 = 2;
			}
			forierr += !(A1.InsertGlobalIndices(MyGlobalElements1[i], NumEntries1, Indices1)==0);
			forierr += !(A1.InsertGlobalIndices(MyGlobalElements1[i], 1, MyGlobalElements1+i)>0); // Put in the diagonal entry
    }
    EPETRA_TEST_ERR(forierr,ierr);
		
    // Finish up
    EPETRA_TEST_ERR(!(A1.FillComplete()==0),ierr);
    
    if(verbose) cout << "Print out tridiagonal matrix, each part on each processor. Index base is one.\n" << endl;
    cout << A1 << endl;
    
  // Release all objects
  delete[] NumNz1;
  delete[] Indices1;
  delete[] MyGlobalElements1;

  delete &A1;
  delete &Map1;
  }

	// Test copy constructor, op=, and reference-counting
	int tempierr = 0;
	if(verbose) cout << "\n*****Checking cpy ctr, op=, and reference counting." << endl;
	tempierr = checkCopyAndAssignment(Comm, verbose);
	EPETRA_TEST_ERR(tempierr, ierr);
	if(verbose && (tempierr == 0)) cout << "Checked OK." << endl;

	// Test shared-ownership code (not implemented yet)
	tempierr = 0;
	if(verbose) cout << "\n*****Checking shared-ownership tests." << endl;
	tempierr = checkSharedOwnership(Comm, verbose);
	EPETRA_TEST_ERR(tempierr, ierr);
	if(verbose && (tempierr == 0)) cout << "Checked OK." << endl;

			
#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif
/* end main */
	return(ierr);
}

//==============================================================================
int checkSharedOwnership(Epetra_Comm& Comm, bool verbose) {
	// check to make sure each function returns 1 when it should
	// check to make sure each function doesn't return 1 when it shouldn't
	int ierr = 0;

	// initialize Map
	const int NumGlobalElements = 10;
	const int IndexBase = 0;
	Epetra_Map Map1(NumGlobalElements, IndexBase, Comm);
	// initialize Graphs
	const int NumIndicesPerRow = 5;
	Epetra_CrsGraph SoleOwner(Copy, Map1, NumIndicesPerRow);
	Epetra_CrsGraph SharedOrig(Copy, Map1, NumIndicesPerRow);
	Epetra_CrsGraph SharedOwner(SharedOrig);
	// arrays used by Insert & Remove
	Epetra_IntSerialDenseVector array1(2);
	array1[0] = NumIndicesPerRow / 2;
	array1[1] = array1[0] + 1;
	Epetra_IntSerialDenseVector array2(NumIndicesPerRow);
	for(int i = 0; i < NumIndicesPerRow; i++)
		array2[i] = i;
	// output variables (declaring them here lets us comment out indiv. tests)
	int soleOutput, sharedOutput;

	// InsertMyIndices
	if(verbose) cout << "InsertMyIndices..." << endl;
	soleOutput = SoleOwner.InsertMyIndices(0, 2, array1.Values());
	sharedOutput = SharedOwner.InsertMyIndices(0, 2, array1.Values());
	EPETRA_TEST_ERR(!(soleOutput == 0), ierr);
	EPETRA_TEST_ERR(!(sharedOutput == 1), ierr);
	if(verbose && ierr > 0) cout << "soleOutput = " << soleOutput << " sharedOutput = " << sharedOutput << endl;

	// SortIndices
	if(verbose) cout << "SortIndices..." << endl;
	soleOutput = SoleOwner.SortIndices();
	sharedOutput = SharedOwner.SortIndices();
	EPETRA_TEST_ERR(!(soleOutput == 0), ierr);
	EPETRA_TEST_ERR(!(sharedOutput == 1), ierr);
	if(verbose && ierr > 0) cout << "soleOutput = " << soleOutput << " sharedOutput = " << sharedOutput << endl;

	// RemoveRedundantIndices
	if(verbose) cout << "RemoveRedundantIndices..." << endl;
	SoleOwner.InsertGlobalIndices(0, 1, array1.Values());
	SharedOwner.InsertGlobalIndices(0, 1, array1.Values());
	soleOutput = SoleOwner.RemoveRedundantIndices();
	sharedOutput = SharedOwner.RemoveRedundantIndices();
	EPETRA_TEST_ERR(!(soleOutput == 0), ierr);
	EPETRA_TEST_ERR(!(sharedOutput == 1), ierr);
	if(verbose && ierr > 0) cout << "soleOutput = " << soleOutput << " sharedOutput = " << sharedOutput << endl;

	// FillComplete (#1)
	if(verbose) cout << "FillComplete..." << endl;
	soleOutput = SoleOwner.FillComplete();
	sharedOutput = SharedOwner.FillComplete();
	EPETRA_TEST_ERR(!(soleOutput == 0), ierr);
	EPETRA_TEST_ERR(!(sharedOutput == 1), ierr);
	if(verbose && ierr > 0) cout << "soleOutput = " << soleOutput << " sharedOutput = " << sharedOutput << endl;

	// OptimizeStorage
	if(verbose) cout << "OptimizeStorage..." << endl;
	soleOutput = SoleOwner.OptimizeStorage();
	sharedOutput = SharedOwner.OptimizeStorage();
	EPETRA_TEST_ERR(!(soleOutput == 0), ierr);
	EPETRA_TEST_ERR(!(sharedOutput == 1), ierr);
	if(verbose && ierr > 0) cout << "soleOutput = " << soleOutput << " sharedOutput = " << sharedOutput << endl;

	// RemoveMyIndices (#1)
	if(verbose) cout << "RemoveMyIndices..." << endl;
	soleOutput = SoleOwner.RemoveMyIndices(0, 1, &array1[1]);
	sharedOutput = SharedOwner.RemoveMyIndices(0, 1, &array1[1]);
	EPETRA_TEST_ERR(!(soleOutput == 0), ierr);
	EPETRA_TEST_ERR(!(sharedOutput == 1), ierr);
	if(verbose && ierr > 0) cout << "soleOutput = " << soleOutput << " sharedOutput = " << sharedOutput << endl;

	// RemoveMyIndices (#2)
	if(verbose) cout << "RemoveMyIndices(#2)..." << endl;
	soleOutput = SoleOwner.RemoveMyIndices(0);
	sharedOutput = SharedOwner.RemoveMyIndices(0);
	EPETRA_TEST_ERR(!(soleOutput == 0), ierr);
	EPETRA_TEST_ERR(!(sharedOutput == 1), ierr);
	if(verbose && ierr > 0) cout << "soleOutput = " << soleOutput << " sharedOutput = " << sharedOutput << endl;

	// FillComplete (#2)
	if(verbose) cout << "FillComplete(#2)..." << endl;
	soleOutput = SoleOwner.FillComplete(SoleOwner.DomainMap(), SoleOwner.RangeMap());
	sharedOutput = SharedOwner.FillComplete(SharedOwner.DomainMap(), SharedOwner.RangeMap());
	EPETRA_TEST_ERR(!(soleOutput == 0), ierr);
	EPETRA_TEST_ERR(!(sharedOutput == 1), ierr);
	if(verbose && ierr > 0) cout << "soleOutput = " << soleOutput << " sharedOutput = " << sharedOutput << endl;

	{
		// make new Graphs so that we can insert Global instead of Local
		// inside of new scope so that we can use same names
		Epetra_CrsGraph SoleOwner(Copy, Map1, NumIndicesPerRow);
		Epetra_CrsGraph SharedOrig(Copy, Map1, NumIndicesPerRow);
		Epetra_CrsGraph SharedOwner(SharedOrig);
		
		int GlobalRow = SoleOwner.GRID(0);

		// InsertGlobalIndices
		if(verbose) cout << "InsertGlobalIndices..." << endl;
		soleOutput = SoleOwner.InsertGlobalIndices(GlobalRow, 2, array2.Values());
		sharedOutput = SharedOwner.InsertGlobalIndices(GlobalRow, 2, array2.Values());
		EPETRA_TEST_ERR(!(soleOutput == 0), ierr);
		EPETRA_TEST_ERR(!(sharedOutput == 1), ierr);
		if(verbose && ierr > 0) cout << "soleOutput = " << soleOutput << " sharedOutput = " << sharedOutput << endl;
		
		// RemoveGlobalIndices (#1)
		if(verbose) cout << "RemoveGlobalIndices..." << endl;
		soleOutput = SoleOwner.RemoveGlobalIndices(GlobalRow, 1, &array2[1]);
		sharedOutput = SharedOwner.RemoveGlobalIndices(GlobalRow, 1, &array2[1]);
		EPETRA_TEST_ERR(!(soleOutput == 0), ierr);
		EPETRA_TEST_ERR(!(sharedOutput == 1), ierr);
		if(verbose && ierr > 0) cout << "soleOutput = " << soleOutput << " sharedOutput = " << sharedOutput << endl;
		
		// RemoveGlobalIndices (#2)
		if(verbose) cout << "RemoveGlobalIndices(#2)..." << endl;
		soleOutput = SoleOwner.RemoveGlobalIndices(GlobalRow);
		sharedOutput = SharedOwner.RemoveGlobalIndices(GlobalRow);
		EPETRA_TEST_ERR(!(soleOutput == 0), ierr);
		EPETRA_TEST_ERR(!(sharedOutput == 1), ierr);
		if(verbose && ierr > 0) cout << "soleOutput = " << soleOutput << " sharedOutput = " << sharedOutput << endl;
		}


	// *PROT* InsertIndices
	// *PROT* MakeIndicesLocal

	return(ierr);

}

//==============================================================================
int checkCopyAndAssignment(Epetra_Comm& Comm, bool verbose) {
	int ierr = 0;
	// check to make sure that copy ctr and op= are working correctly
	// (making level 1 deep copies)

	// create initial Map and Graph
	const int NumIndicesPerRow = 10;
	const int NumGlobalElements = 50;
	const int IndexBase = 0;
	Epetra_Map Map1(NumGlobalElements, IndexBase, Comm);
	Epetra_CrsGraph Graph1(Copy, Map1, NumIndicesPerRow);
	int g1count = Graph1.ReferenceCount();
	const Epetra_CrsGraphData* g1addr = Graph1.DataPtr();
	EPETRA_TEST_ERR(!(g1count == 1), ierr);
	if(verbose) cout << "Graph1 created (def ctr). data addr = " << g1addr << " ref. count = " << g1count << endl;

	//test copy constructor
	{
		Epetra_CrsGraph Graph2(Graph1);
		int g2count = Graph2.ReferenceCount();
		const Epetra_CrsGraphData* g2addr = Graph2.DataPtr();
		EPETRA_TEST_ERR(!(g2count == g1count+1), ierr);
		EPETRA_TEST_ERR(!(g2addr == g1addr), ierr);
		if(verbose) cout << "Graph2 created (cpy ctr). data addr = " << g2addr << " ref. count = " << g2count << endl;
	}
	int g1newcount = Graph1.ReferenceCount();
	const Epetra_CrsGraphData* g1newaddr = Graph1.DataPtr();
	EPETRA_TEST_ERR(!(g1newcount == g1count), ierr);
	EPETRA_TEST_ERR(!(g1newaddr = g1addr), ierr);
	if(verbose) cout << "Graph2 destroyed. Graph1 data addr = " << g1newaddr << " ref. count = " << g1newcount << endl;

	//test op=
	Epetra_CrsGraph Graph3(Copy, Map1, NumIndicesPerRow+1);
	int g3count = Graph3.ReferenceCount();
	const Epetra_CrsGraphData* g3addr = Graph3.DataPtr();
	EPETRA_TEST_ERR(!(g3addr != g1addr), ierr);
	if(verbose) cout << "Graph3 created (op= before). data addr = " << g3addr << " ref. count = " << g3count << endl;
	Graph3 = Graph1;
	g3count = Graph3.ReferenceCount();
	g3addr = Graph3.DataPtr();
	EPETRA_TEST_ERR(!(g3count == g1count+1), ierr);
	EPETRA_TEST_ERR(!(g3addr == g1addr), ierr);
	if(verbose) cout << "Graph3 set equal to Graph1 (op= after). data addr = " << g3addr << " ref. count = " << g3count << endl;

	return(ierr);
}

//==============================================================================
int check(Epetra_CrsGraph& A, int NumMyRows1, int NumGlobalRows1, int NumMyNonzeros1,
	  int NumGlobalNonzeros1, int* MyGlobalElements, bool verbose)
{  

  int ierr = 0;
	int i;
	int j;
	int forierr = 0;
  int NumGlobalIndices;
  int NumMyIndices;
	int* MyViewIndices;
  int MaxNumIndices = A.MaxNumIndices();
  int* MyCopyIndices = new int[MaxNumIndices];
  int* GlobalCopyIndices = new int[MaxNumIndices];

  // Test query functions

  int NumMyRows = A.NumMyRows();
  if(verbose) cout << "Number of local Rows = " << NumMyRows << endl;

  EPETRA_TEST_ERR(!(NumMyRows==NumMyRows1),ierr);

  int NumMyNonzeros = A.NumMyNonzeros();
  if(verbose) cout << "Number of local Nonzero entries = " << NumMyNonzeros << endl;

  EPETRA_TEST_ERR(!(NumMyNonzeros==NumMyNonzeros1),ierr);

  int NumGlobalRows = A.NumGlobalRows();
  if(verbose) cout << "Number of global Rows = " << NumGlobalRows << endl;

  EPETRA_TEST_ERR(!(NumGlobalRows==NumGlobalRows1),ierr);

  int NumGlobalNonzeros = A.NumGlobalNonzeros();
  if(verbose) cout << "Number of global Nonzero entries = " << NumGlobalNonzeros << endl;

  EPETRA_TEST_ERR(!(NumGlobalNonzeros==NumGlobalNonzeros1),ierr);

  // GlobalRowView should be illegal (since we have local indices)

  EPETRA_TEST_ERR(!(A.ExtractGlobalRowView(A.RowMap().MaxMyGID(), NumGlobalIndices, GlobalCopyIndices)==-2),ierr);

  // Other binary tests

  EPETRA_TEST_ERR(A.NoDiagonal(),ierr);
  EPETRA_TEST_ERR(!(A.Filled()),ierr);
  EPETRA_TEST_ERR(!(A.Sorted()),ierr);
  EPETRA_TEST_ERR(!(A.MyGRID(A.RowMap().MaxMyGID())),ierr);
  EPETRA_TEST_ERR(!(A.MyGRID(A.RowMap().MinMyGID())),ierr);
  EPETRA_TEST_ERR(A.MyGRID(1+A.RowMap().MaxMyGID()),ierr);
  EPETRA_TEST_ERR(A.MyGRID(-1+A.RowMap().MinMyGID()),ierr);
  EPETRA_TEST_ERR(!(A.MyLRID(0)),ierr);
  EPETRA_TEST_ERR(!(A.MyLRID(NumMyRows-1)),ierr);
  EPETRA_TEST_ERR(A.MyLRID(-1),ierr);
  EPETRA_TEST_ERR(A.MyLRID(NumMyRows),ierr);
    
  forierr = 0;
  for(i = 0; i < NumMyRows; i++) {
    int Row = A.GRID(i);
    A.ExtractGlobalRowCopy(Row, MaxNumIndices, NumGlobalIndices, GlobalCopyIndices);
    A.ExtractMyRowView(i, NumMyIndices, MyViewIndices);
    forierr += !(NumGlobalIndices==NumMyIndices);
    for(j = 1; j < NumMyIndices; j++) EPETRA_TEST_ERR(!(MyViewIndices[j-1]<MyViewIndices[j]),ierr);
    for(j = 0; j < NumGlobalIndices; j++) {
			forierr += !(GlobalCopyIndices[j]==A.GCID(MyViewIndices[j]));
			forierr += !(A.LCID(GlobalCopyIndices[j])==MyViewIndices[j]);
    }
  }
  EPETRA_TEST_ERR(forierr,ierr);
  forierr = 0;
  for(i = 0; i < NumMyRows; i++) {
    int Row = A.GRID(i);
    A.ExtractGlobalRowCopy(Row, MaxNumIndices, NumGlobalIndices, GlobalCopyIndices);
    A.ExtractMyRowCopy(i, MaxNumIndices, NumMyIndices, MyCopyIndices);
    forierr += !(NumGlobalIndices==NumMyIndices);
    for(j = 1; j < NumMyIndices; j++) 
			EPETRA_TEST_ERR(!(MyCopyIndices[j-1]<MyCopyIndices[j]),ierr);
    for(j = 0; j < NumGlobalIndices; j++) {
			forierr += !(GlobalCopyIndices[j]==A.GCID(MyCopyIndices[j]));
			forierr += !(A.LCID(GlobalCopyIndices[j])==MyCopyIndices[j]);
    }

  }
  EPETRA_TEST_ERR(forierr,ierr);

  delete[] MyCopyIndices;
  delete[] GlobalCopyIndices;

  if(verbose) cout << "Rows sorted check OK" << endl;

  return(ierr);
}
