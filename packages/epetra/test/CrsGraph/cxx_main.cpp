#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsGraph.h"
#include "Epetra_Map.h"

// Prototype

int check(Epetra_CrsGraph& A, int NumMyRows1, int NumGlobalRows1, int NumMyNonzeros1,
	  int NumGlobalNonzeros1, int * MyGlobalElements, bool verbose);

 int main(int argc, char *argv[])
{
  int ierr = 0, i, j;
  int NumIndices;
  int * Indices;
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
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  //  char tmp;
  //  if (rank==0) cout << "Press any key to continue..."<< endl;
  //  if (rank==0) cin >> tmp;
  //  Comm.Barrier();

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();
  if (verbose) cout << "Processor "<<MyPID<<" of "<< NumProc
		    << " is alive."<<endl;

  bool verbose1 = verbose;

  // Redefine verbose to only print on PE 0
  if (verbose && rank!=0) verbose = false;

  int NumMyEquations = 10000;
  int NumGlobalEquations = NumMyEquations*NumProc+EPETRA_MIN(NumProc,3);
  if (MyPID < 3) NumMyEquations++;
  int IndexBase = 0;
  int ElementSize = 7;
  bool DistributedGlobal = (NumGlobalEquations>NumMyEquations);

  // Construct a Map that puts approximately the same Number of equations on each processor

  Epetra_Map& Map = *new Epetra_Map(NumGlobalEquations, NumMyEquations, 0, Comm);
  
  // Get update list and number of local equations from newly created Map
  int * MyGlobalElements = new int[Map.NumMyElements()];
  Map.MyGlobalElements(MyGlobalElements);

  // Create an integer vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation on this processor

  int * NumNz = new int[NumMyEquations];

  // We are building a tridiagonal matrix where each row has (-1 2 -1)
  // So we need 2 off-diagonal terms (except for the first and last equation)

  for (i=0; i<NumMyEquations; i++)
    if (MyGlobalElements[i]==0 || MyGlobalElements[i] == NumGlobalEquations-1)
      NumNz[i] = 1;
    else
      NumNz[i] = 2;

  // Create a Epetra_CrsGraph

  Epetra_CrsGraph& A = *new Epetra_CrsGraph(Copy, Map, NumNz);
  assert(!A.IndicesAreGlobal());
  assert(!A.IndicesAreLocal());
  
  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1


  Indices = new int[2];
  int NumEntries;
  
  for (i=0; i<NumMyEquations; i++)
    {
    if (MyGlobalElements[i]==0)
      {
	Indices[0] = 1;
	NumEntries = 1;
      }
    else if (MyGlobalElements[i] == NumGlobalEquations-1)
      {
	Indices[0] = NumGlobalEquations-2;
	NumEntries = 1;
      }
    else
      {
	Indices[0] = MyGlobalElements[i]-1;
	Indices[1] = MyGlobalElements[i]+1;
	NumEntries = 2;
      }
     assert(A.InsertGlobalIndices(MyGlobalElements[i], NumEntries, Indices)==0);
     assert(A.InsertGlobalIndices(MyGlobalElements[i], 1, MyGlobalElements+i)>0); // Put in the diagonal entry (should cause realloc)
    }

  delete [] Indices;
  
  // Finish up
  assert(A.IndicesAreGlobal());
  assert(A.TransformToLocal()==0);
  assert(A.IndicesAreLocal());
  assert(!A.StorageOptimized());
  A.OptimizeStorage();
  assert(A.StorageOptimized());
  assert(!A.UpperTriangular());
  assert(!A.LowerTriangular());

  if (verbose) cout << "\n\n*****Testing variable entry constructor" << endl<< endl;

  int NumMyNonzeros = 3*NumMyEquations;
  if (A.LRID(0)>=0) NumMyNonzeros--; // If I own first global row, then there is one less nonzero
  if (A.LRID(NumGlobalEquations-1)>=0) NumMyNonzeros--; // If I own last global row, then there is one less nonzero

  assert(check(A, NumMyEquations, NumGlobalEquations, NumMyNonzeros, 3*NumGlobalEquations-2, 
	       MyGlobalElements, verbose)==0);

  for (i=0; i<NumMyEquations; i++) assert(A.NumGlobalIndices(MyGlobalElements[i])==NumNz[i]+1);
  for (i=0; i<NumMyEquations; i++) assert(A.NumMyIndices(i)==NumNz[i]+1);

  if (verbose) cout << "\n\nNumIndices function check OK" << endl<< endl;

  delete &A;

  if (debug) Comm.Barrier();

  if (verbose) cout << "\n\n*****Testing constant entry constructor" << endl<< endl;

  Epetra_CrsGraph & AA = *new Epetra_CrsGraph(Copy, Map, 5);
  
  if (debug) Comm.Barrier();

  for (i=0; i< NumMyEquations; i++) AA.InsertGlobalIndices(MyGlobalElements[i], 1, MyGlobalElements+i);

  // Note:  All processors will call the following Insert routines, but only the processor
  //        that owns it will actually do anything

  int One = 1;
  if (AA.MyGlobalRow(0)) assert(AA.InsertGlobalIndices(0, 0, &One)==0);
  else assert(AA.InsertGlobalIndices(0, 1, &One)==-1);
  assert(AA.TransformToLocal()==0);
  assert(!AA.StorageOptimized());
  assert(AA.UpperTriangular());
  assert(AA.LowerTriangular());
  
  if (debug) Comm.Barrier();
  assert(check(AA, NumMyEquations, NumGlobalEquations, NumMyEquations, NumGlobalEquations, 
	       MyGlobalElements, verbose)==0);

  if (debug) Comm.Barrier();

  for (i=0; i<NumMyEquations; i++) assert(AA.NumGlobalIndices(MyGlobalElements[i])==1);

  if (verbose) cout << "\n\nNumIndices function check OK" << endl<< endl;

  if (debug) Comm.Barrier();

  if (verbose) cout << "\n\n*****Testing copy constructor" << endl<< endl;

  Epetra_CrsGraph & B = *new Epetra_CrsGraph(AA);
  delete &AA;

  assert(check(B, NumMyEquations, NumGlobalEquations, NumMyEquations, NumGlobalEquations, 
	       MyGlobalElements, verbose)==0);

  for (i=0; i<NumMyEquations; i++) assert(B.NumGlobalIndices(MyGlobalElements[i])==1);

  if (verbose) cout << "\n\nNumIndices function check OK" << endl<< endl;

  if (debug) Comm.Barrier();

  if (verbose) cout << "\n\n*****Testing post construction modifications" << endl<< endl;

  assert(B.InsertGlobalIndices(0, 1, &One)==-2);
  delete &B;

  // Release all objects
  delete [] MyGlobalElements;
  delete [] NumNz;
  delete &Map;
			

  if (verbose1) {
    // Test ostream << operator (if verbose1)
    // Construct a Map that puts 2 equations on each PE
    
    int NumMyElements1 = 4;
    int NumMyEquations1 = NumMyElements1;
    int NumGlobalEquations1 = NumMyEquations1*NumProc;

    Epetra_Map& Map1 = *new Epetra_Map(NumGlobalEquations1, NumMyElements1, 1, Comm);
    
    // Get update list and number of local equations from newly created Map
    int * MyGlobalElements1 = new int[Map1.NumMyElements()];
    Map1.MyGlobalElements(MyGlobalElements1);
    
    // Create an integer vector NumNz that is used to build the Petra Matrix.
    // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation on this processor
    
    int * NumNz1 = new int[NumMyEquations1];
    
    // We are building a tridiagonal matrix where each row has (-1 2 -1)
    // So we need 2 off-diagonal terms (except for the first and last equation)
    
    for (i=0; i<NumMyEquations1; i++)
      if (MyGlobalElements1[i]==1 || MyGlobalElements1[i] == NumGlobalEquations1)
	NumNz1[i] = 1;
      else
	NumNz1[i] = 2;
    
    // Create a Epetra_Graph using 1-based arithmetic
    
    Epetra_CrsGraph& A1 = *new Epetra_CrsGraph(Copy, Map1, NumNz1);
    
    // Add  rows one-at-a-time
    // Need some vectors to help
    // Off diagonal Values will always be -1
    
    
    int *Indices1 = new int[2];
    int NumEntries1;
    
    for (i=0; i<NumMyEquations1; i++)
      {
	if (MyGlobalElements1[i]==1)
	  {
	    Indices1[0] = 2;
	    NumEntries1 = 1;
	  }
	else if (MyGlobalElements1[i] == NumGlobalEquations1)
	  {
	    Indices1[0] = NumGlobalEquations1-1;
	    NumEntries1 = 1;
	  }
	else
	  {
	    Indices1[0] = MyGlobalElements1[i]-1;
	    Indices1[1] = MyGlobalElements1[i]+1;
	    NumEntries1 = 2;
	  }
	assert(A1.InsertGlobalIndices(MyGlobalElements1[i], NumEntries1, Indices1)==0);
	assert(A1.InsertGlobalIndices(MyGlobalElements1[i], 1, MyGlobalElements1+i)>0); // Put in the diagonal entry
      }
    
    // Finish up
    assert(A1.TransformToLocal()==0);
    
    if (verbose) cout << "\n\nPrint out tridiagonal matrix, each part on each processor. Index base is one.\n\n" << endl;
    cout << A1 << endl;
    
  // Release all objects
  delete [] NumNz1;
  delete [] Indices1;
  delete [] MyGlobalElements1;

  delete &A1;
  delete &Map1;
  }
			
#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif
/* end main
*/
return ierr ;
}

int check(Epetra_CrsGraph& A, int NumMyRows1, int NumGlobalRows1, int NumMyNonzeros1,
	  int NumGlobalNonzeros1, int * MyGlobalElements, bool verbose) {  

  int i, j;
  int NumGlobalIndices;
  int NumMyIndices, * MyViewIndices, *GlobalViewIndices;
  int MaxNumIndices = A.MaxNumIndices();
  int * MyCopyIndices = new int[MaxNumIndices];
  int * GlobalCopyIndices = new int[MaxNumIndices];

  // Test query functions

  int NumMyRows = A.NumMyRows();
  if (verbose) cout << "\n\nNumber of local Rows = " << NumMyRows << endl<< endl;

  assert(NumMyRows==NumMyRows1);

  int NumMyNonzeros = A.NumMyNonzeros();
  if (verbose) cout << "\n\nNumber of local Nonzero entries = " << NumMyNonzeros << endl<< endl;

  assert(NumMyNonzeros==NumMyNonzeros1);

  int NumGlobalRows = A.NumGlobalRows();
  if (verbose) cout << "\n\nNumber of global Rows = " << NumGlobalRows << endl<< endl;

  assert(NumGlobalRows==NumGlobalRows1);

  int NumGlobalNonzeros = A.NumGlobalNonzeros();
  if (verbose) cout << "\n\nNumber of global Nonzero entries = " << NumGlobalNonzeros << endl<< endl;

  assert(NumGlobalNonzeros==NumGlobalNonzeros1);

  // GlobalRowView should be illegal (since we have local indices)

  assert(A.ExtractGlobalRowView(A.RowMap().MaxMyGID(), NumGlobalIndices, GlobalCopyIndices)==-2);

  // Other binary tests

  assert(!A.NoDiagonal());
  assert(A.Filled());
  assert(A.Sorted());
  assert(A.MyGRID(A.RowMap().MaxMyGID()));
  assert(A.MyGRID(A.RowMap().MinMyGID()));
  assert(!A.MyGRID(1+A.RowMap().MaxMyGID()));
  assert(!A.MyGRID(-1+A.RowMap().MinMyGID()));
  assert(A.MyLRID(0));
  assert(A.MyLRID(NumMyRows-1));
  assert(!A.MyLRID(-1));
  assert(!A.MyLRID(NumMyRows));
    

  for (i=0; i<NumMyRows; i++) {
    int Row = A.GRID(i);
    A.ExtractGlobalRowCopy(Row, MaxNumIndices, NumGlobalIndices, GlobalCopyIndices);
    A.ExtractMyRowView(i, NumMyIndices, MyViewIndices);
    assert(NumGlobalIndices==NumMyIndices);
    for (j=1; j<NumMyIndices; j++) assert(MyViewIndices[j-1]<MyViewIndices[j]);
    for (j=0; j<NumGlobalIndices; j++) {
	assert(GlobalCopyIndices[j]==A.GCID(MyViewIndices[j]));
	assert(A.LCID(GlobalCopyIndices[j])==MyViewIndices[j]);
    }
  }

  for (i=0; i<NumMyRows; i++) {
    int Row = A.GRID(i);
    A.ExtractGlobalRowCopy(Row, MaxNumIndices, NumGlobalIndices, GlobalCopyIndices);
    A.ExtractMyRowCopy(i, MaxNumIndices, NumMyIndices, MyCopyIndices);
    assert(NumGlobalIndices==NumMyIndices);
    for (j=1; j<NumMyIndices; j++) assert(MyCopyIndices[j-1]<MyCopyIndices[j]);
    for (j=0; j<NumGlobalIndices; j++) {
	assert(GlobalCopyIndices[j]==A.GCID(MyCopyIndices[j]));
	assert(A.LCID(GlobalCopyIndices[j])==MyCopyIndices[j]);
    }

  }

  delete [] MyCopyIndices;
  delete [] GlobalCopyIndices;

  if (verbose) cout << "\n\nRows sorted check OK" << endl<< endl;

  return(0);
}
