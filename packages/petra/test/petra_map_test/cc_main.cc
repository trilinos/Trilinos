// Petra_Comm Test routine

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
#include "Petra_Map.h"
#include "Petra_LocalMap.h"
#include "checkmap.h"

int main(int argc, char *argv[]) {


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
  int NumMyElements1 = NumMyElements; // Used for local map
  int NumGlobalElements = NumMyElements*NumProc+minfn(NumProc,3);
  if (MyPID < 3) NumMyElements++;
  int IndexBase = 0;
  bool DistributedGlobal = (NumGlobalElements>NumMyElements);
  

  // Test Petra-defined uniform linear distribution constructor
  Petra_Map * Map = new Petra_Map(NumGlobalElements, IndexBase, Comm);
  if (verbose) cout << "Checking Petra_Map(NumGlobalElements, IndexBase, Comm)" << endl;
  int ierr = checkmap(*Map, NumGlobalElements, NumMyElements, 0, 
		  IndexBase, Comm, DistributedGlobal);

  if (verbose)
    if (ierr==0) cout << "Checked OK\n\n" <<endl;
    else cout << "Error code: "<< ierr << endl;

  assert(ierr==0);

  delete Map;

  // Test User-defined linear distribution constructor
  Map = new Petra_Map(NumGlobalElements, NumMyElements, IndexBase, Comm);

  if (verbose) cout << "Checking Petra_Map(NumGlobalElements, NumMyElements, IndexBase, Comm)" << endl;
  ierr = checkmap(*Map, NumGlobalElements, NumMyElements, 0, 
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
  for (int i = 0; i<NumMyElements; i++) MyGlobalElements[i] = MaxMyGID-i;

  Map = new Petra_Map(NumGlobalElements, NumMyElements, MyGlobalElements, 
		      IndexBase, Comm);
  if (verbose) cout << "Checking Petra_Map(NumGlobalElements, NumMyElements, MyGlobalElements,  IndexBase, Comm)" << endl;
  ierr = checkmap(*Map, NumGlobalElements, NumMyElements, MyGlobalElements, 
		  IndexBase, Comm, DistributedGlobal);

  if (verbose)
    if (ierr==0) cout << "Checked OK\n\n" <<endl;
    else cout << "Error code: "<< ierr << endl;

  assert(ierr==0);

  // Test Copy constructor
  Petra_Map * Map1 = new Petra_Map(*Map);

  if (verbose) cout << "Checking Petra_Map(*Map)" << endl;
  ierr = checkmap(*Map1, NumGlobalElements, NumMyElements, MyGlobalElements, 
		  IndexBase, Comm, DistributedGlobal);

  if (verbose)
    if (ierr==0) cout << "Checked OK\n\n" <<endl;
    else cout << "Error code: "<< ierr << endl;

  assert(ierr==0);

  Petra_Map * SmallMap = 0;
  if (verbose1) {
    // Build a small map for test cout.  Use 10 elements from current map
    int * MyEls = Map->MyGlobalElements();
    int IndBase = Map->IndexBase();
    int MyLen = minfn(10+Comm.MyPID(),Map->NumMyElements());
    SmallMap = new Petra_Map(-1, MyLen, MyEls, IndBase, Comm);
  }

  delete [] MyGlobalElements;
  delete Map;
  delete Map1;


  // Test LocalMap constructor
  Petra_LocalMap *LocalMap = new Petra_LocalMap(NumMyElements1, IndexBase, 
						Comm);
  if (verbose) cout << "Checking Petra_LocalMap(NumMyElements1, IndexBase, Comm)" << endl;
  ierr = checkmap(*LocalMap, NumMyElements1, NumMyElements1, 0, 
		  IndexBase, Comm, false);

  if (verbose)
    if (ierr==0) cout << "Checked OK\n\n" <<endl;
    else cout << "Error code: "<< ierr << endl;

  assert(ierr==0);

  // Test Copy constructor
  Petra_LocalMap * LocalMap1 = new Petra_LocalMap(*LocalMap);
  if (verbose) cout << "Checking Petra_LocalMap(*LocalMap)" << endl;
  ierr = checkmap(*LocalMap1, NumMyElements1, NumMyElements1, 
		  0, IndexBase, Comm, false);

  if (verbose)
    if (ierr==0) cout << "Checked OK\n\n" <<endl;
    else cout << "Error code: "<< ierr << endl;

  assert(ierr==0);

  delete LocalMap1;
  delete LocalMap;
  if (verbose1) {
    if (verbose) cout << "Test ostream << operator" << endl << flush;
    cout << *SmallMap;
    delete SmallMap;
  }




  delete &Comm;

#ifdef PETRA_MPI
  MPI_Finalize();
#endif

  return 0;
}

