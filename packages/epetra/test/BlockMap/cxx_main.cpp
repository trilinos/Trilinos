// Epetra_BlockMap Test routine

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#endif
#include "Epetra_SerialComm.h"
#include "Epetra_Time.h"
#include "Epetra_BlockMap.h"
#include "checkmap.h"

int main(int argc, char *argv[]) {

  int i;

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#else

  int size = 1; // Serial case (not using MPI)
  int rank = 0;

#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;


#ifdef EPETRA_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  int temp;
  if (Comm.MyPID()==0)
    cin >> temp;
  Comm.Barrier();
  Comm.SetTracebackMode(0); // This should shut down any error traceback reporting
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();
  cout << Comm << endl << flush;
  Comm.Barrier();
  bool verbose1 = verbose;
  if (verbose) verbose = (MyPID==0);

  int NumMyElements = 10000;
  int NumGlobalElements = NumMyElements*NumProc+EPETRA_MIN(NumProc,3);
  if (MyPID < 3) NumMyElements++;
  int IndexBase = 0;
  int ElementSize = 7;
  bool DistributedGlobal = (NumGlobalElements>NumMyElements);
  
  Epetra_BlockMap * Map;

  // Test exceptions

  if (verbose) 
    cout << "*******************************************************************************************" << endl
	 << "        Testing Exceptions (Expect error messages if EPETRA_NO_ERROR_REPORTS is not defined" << endl
	 << "*******************************************************************************************" << endl
	 << endl << endl;

  try {
    if (verbose) cout << "Checking Epetra_BlockMap(-2, ElementSize, IndexBase, Comm)" << endl;
    Map = new Epetra_BlockMap(-2, ElementSize, IndexBase, Comm);
  }
  catch (int Error) {
    if (Error!=-1) {
      if (verbose) cout << "Error code = " << Error << "Should be -1" << endl;
      return 1;
    }
    if (verbose) cout << "Checked OK\n\n" << endl;
  }

  try {
    if (verbose) cout << "Checking Epetra_BlockMap(2, 3, ElementSize, IndexBase, Comm)" << endl;
    Map = new Epetra_BlockMap(2, 3, ElementSize, IndexBase, Comm);
  }
  catch (int Error) {
    if (Error!=-4) {
      if (verbose) cout << "Error code = " << Error << "Should be -4" << endl;
      return 1;
    }
    if (verbose) cout << "Checked OK\n\n" << endl;
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
  int ierr = checkmap(*Map, NumGlobalElements, NumMyElements, 0, ElementSize, 0,
		      NumGlobalElements*ElementSize, NumMyElements*ElementSize,
		      IndexBase, Comm, DistributedGlobal);

  if (verbose)
    if (ierr==0) cout << "Checked OK\n\n" <<endl;
    else cout << "Error code: "<< ierr << endl;

  assert(ierr==0);

  delete Map;

  // Test User-defined linear distribution constructor
  Map = new Epetra_BlockMap(NumGlobalElements, NumMyElements, ElementSize, IndexBase, Comm);

  if (verbose) cout << "Checking Epetra_BlockMap(NumGlobalElements, NumMyElements, ElementSize, IndexBase, Comm)" << endl;
  ierr = checkmap(*Map, NumGlobalElements, NumMyElements, 0, ElementSize, 0,
		  NumGlobalElements*ElementSize, NumMyElements*ElementSize,
		  IndexBase, Comm, DistributedGlobal);

  if (verbose)
    if (ierr==0) cout << "Checked OK\n\n" <<endl;
    else cout << "Error code: "<< ierr << endl;

  assert(ierr==0);

  delete Map;

  // Test User-defined arbitrary distribution constructor
  // Generate Global Element List.  Do in reverse for fun!

  int * MyGlobalElements = new int[NumMyElements];
  int MaxMyGID = (Comm.MyPID()+1)*NumMyElements-1+IndexBase;
  if (Comm.MyPID()>2) MaxMyGID+=3;
  for (i = 0; i<NumMyElements; i++) MyGlobalElements[i] = MaxMyGID-i;

  Map = new Epetra_BlockMap(NumGlobalElements, NumMyElements, MyGlobalElements, ElementSize,
		      IndexBase, Comm);
  if (verbose) cout << "Checking Epetra_BlockMap(NumGlobalElements, NumMyElements, MyGlobalElements,  ElementSize, IndexBase, Comm)" << endl;
  ierr = checkmap(*Map, NumGlobalElements, NumMyElements, MyGlobalElements, ElementSize, 0,
		  NumGlobalElements*ElementSize, NumMyElements*ElementSize,
		  IndexBase, Comm, DistributedGlobal);

  if (verbose)
    if (ierr==0) cout << "Checked OK\n\n" <<endl;
    else cout << "Error code: "<< ierr << endl;

  assert(ierr==0);
  Epetra_BlockMap * Map3 = new Epetra_BlockMap(*Map);// A map to test the SameAs method later

  delete Map;

  int * ElementSizeList = new int[NumMyElements];
  int NumMyEquations = 0;
  int NumGlobalEquations = 0;
  for (i = 0; i<NumMyElements; i++) 
    {
      ElementSizeList[i] = i%6 + 2; // elementsizes go from 2 to 7
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
  Map = new Epetra_BlockMap(NumGlobalElements, NumMyElements, MyGlobalElements, ElementSizeList,
		      IndexBase, Comm);
  if (verbose) cout << "Checking Epetra_BlockMap(NumGlobalElements, NumMyElements, MyGlobalElements,  ElementSizeList, IndexBase, Comm)" << endl;
  ierr = checkmap(*Map, NumGlobalElements, NumMyElements, MyGlobalElements, ElementSize, ElementSizeList,
		  NumGlobalEquations, NumMyEquations,
		  IndexBase, Comm, DistributedGlobal);

  if (verbose)
    if (ierr==0) cout << "Checked OK\n\n" <<endl;
    else cout << "Error code: "<< ierr << endl;

  assert(ierr==0);

  // Test Copy constructor
  Epetra_BlockMap * Map1 = new Epetra_BlockMap(*Map);

  // Test SameAs() method
  bool same = Map1->SameAs(*Map);
  assert(same==true);// should return true since Map1 is a copy of Map

  Epetra_BlockMap * Map2 = new Epetra_BlockMap(NumGlobalElements,NumMyElements,MyGlobalElements,ElementSizeList,IndexBase,Comm);
  same = Map2->SameAs(*Map);
  assert(same==true); // Map and Map2 were created with the same sets of parameters
  delete Map2;

  // now test SameAs() on some maps that are different

  Map2 = new Epetra_BlockMap(NumGlobalElements,NumMyElements,MyGlobalElements,ElementSizeList,IndexBase-1,Comm);
  same = Map2->SameAs(*Map);
  assert(same==false); // IndexBases are different
  delete Map2;

  int *ElementSizeList1 = new int[NumMyElements];
  for (i=0; i<NumMyElements; i++) ElementSizeList1[i] = i%5 + 2; // element sizes go from 2 to 6
  Map2 = new Epetra_BlockMap(NumGlobalElements,NumMyElements,MyGlobalElements,ElementSizeList1,IndexBase,Comm);
  same = Map2->SameAs(*Map);
  assert(same==false); // ElementSizes are different
  delete [] ElementSizeList1;
  delete Map2;

  same = Map3->SameAs(*Map);
  assert(same==false); // Map3 saved from an earlier test
  delete Map3;


  if (verbose) cout << "Checking Epetra_BlockMap(*Map)" << endl;
  ierr = checkmap(*Map1, NumGlobalElements, NumMyElements, MyGlobalElements, ElementSize, ElementSizeList,
		  NumGlobalEquations, NumMyEquations,
		  IndexBase, Comm, DistributedGlobal);

  if (verbose)
    if (ierr==0) cout << "Checked OK\n\n" <<endl;
    else cout << "Error code: "<< ierr << endl;

  assert(ierr==0);

  if (verbose1) {
    if (verbose) cout << "Test ostream << operator" << endl << flush;
    // Build a small map for test cout.  Use 10 elements from current map
    int * MyEls = Map->MyGlobalElements();
    int * MySz  = Map->ElementSizeList();
    int IndBase = Map->IndexBase();
    int MyLen = EPETRA_MIN(10+Comm.MyPID(),Map->NumMyElements());
    Epetra_BlockMap * SmallMap = new Epetra_BlockMap(-1, MyLen, MyEls, MySz, IndBase, Comm);
    cout << *SmallMap;
    delete SmallMap;
  }
	

  delete [] ElementSizeList;
  delete [] MyGlobalElements;
  delete Map;
  delete Map1;


#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return 0;
}

