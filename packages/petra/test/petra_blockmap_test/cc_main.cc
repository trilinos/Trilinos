// Petra_BlockMap Test routine

#include <math.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <assert.h>
#ifdef PETRA_MPI
#include <mpi.h>
#endif
#include "Petra_Comm.h"
#include "Petra_Time.h"
#include "Petra_BlockMap.h"
#include "checkmap.h"

int main(int argc, char *argv[]) {

  int i;

#ifdef PETRA_MPI

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


#ifdef PETRA_MPI
  Petra_Comm & Comm = *new Petra_Comm( MPI_COMM_WORLD );
#else
  Petra_Comm & Comm = *new Petra_Comm();
#endif

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();
  if (verbose) cout << "Processor "<<MyPID<<" of "<< NumProc
		    << " is alive."<<endl;

  bool verbose1 = verbose;
  verbose = (MyPID==0);

  int NumMyElements = 10000;
  int NumGlobalElements = NumMyElements*NumProc+minfn(NumProc,3);
  if (MyPID < 3) NumMyElements++;
  int IndexBase = 0;
  int ElementSize = 7;
  bool DistributedGlobal = (NumGlobalElements>NumMyElements);
  

  // Test Petra-defined uniform linear distribution constructor
  Petra_BlockMap * Map = new Petra_BlockMap(NumGlobalElements, ElementSize, IndexBase, Comm);
  if (verbose) cout << "Checking Petra_BlockMap(NumGlobalElements, ElementSize, IndexBase, Comm)" << endl;
  int ierr = checkmap(*Map, NumGlobalElements, NumMyElements, 0, ElementSize, 0,
		      NumGlobalElements*ElementSize, NumMyElements*ElementSize,
		      IndexBase, Comm, DistributedGlobal);

  if (verbose)
    if (ierr==0) cout << "Checked OK\n\n" <<endl;
    else cout << "Error code: "<< ierr << endl;

  assert(ierr==0);

  delete Map;

  // Test User-defined linear distribution constructor
  Map = new Petra_BlockMap(NumGlobalElements, NumMyElements, ElementSize, IndexBase, Comm);

  if (verbose) cout << "Checking Petra_BlockMap(NumGlobalElements, NumMyElements, ElementSize, IndexBase, Comm)" << endl;
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

  Map = new Petra_BlockMap(NumGlobalElements, NumMyElements, MyGlobalElements, ElementSize,
		      IndexBase, Comm);
  if (verbose) cout << "Checking Petra_BlockMap(NumGlobalElements, NumMyElements, MyGlobalElements,  ElementSize, IndexBase, Comm)" << endl;
  ierr = checkmap(*Map, NumGlobalElements, NumMyElements, MyGlobalElements, ElementSize, 0,
		  NumGlobalElements*ElementSize, NumMyElements*ElementSize,
		  IndexBase, Comm, DistributedGlobal);

  if (verbose)
    if (ierr==0) cout << "Checked OK\n\n" <<endl;
    else cout << "Error code: "<< ierr << endl;

  assert(ierr==0);

  delete Map;

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
  Map = new Petra_BlockMap(NumGlobalElements, NumMyElements, MyGlobalElements, ElementSizeList,
		      IndexBase, Comm);
  if (verbose) cout << "Checking Petra_BlockMap(NumGlobalElements, NumMyElements, MyGlobalElements,  ElementSizeList, IndexBase, Comm)" << endl;
  ierr = checkmap(*Map, NumGlobalElements, NumMyElements, MyGlobalElements, ElementSize, ElementSizeList,
		  NumGlobalEquations, NumMyEquations,
		  IndexBase, Comm, DistributedGlobal);

  if (verbose)
    if (ierr==0) cout << "Checked OK\n\n" <<endl;
    else cout << "Error code: "<< ierr << endl;

  assert(ierr==0);

  // Test Copy constructor
  Petra_BlockMap * Map1 = new Petra_BlockMap(*Map);

  if (verbose) cout << "Checking Petra_BlockMap(*Map)" << endl;
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
    int MyLen = minfn(10+Comm.MyPID(),Map->NumMyElements());
    Petra_BlockMap * SmallMap = new Petra_BlockMap(-1, MyLen, MyEls, MySz, IndBase, Comm);
    cout << *SmallMap;
    delete SmallMap;
  }
	

  delete [] ElementSizeList;
  delete [] MyGlobalElements;
  delete Map;
  delete Map1;


  delete &Comm;

#ifdef PETRA_MPI
  MPI_Finalize();
#endif

  return 0;
}

