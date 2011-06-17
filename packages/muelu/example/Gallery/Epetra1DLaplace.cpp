#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <string>
#include <cmath>
#include <vector>

#include <Teuchos_CommandLineProcessor.hpp>

#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Version.h"

#include <string>
#include <sstream>
#include <fstream>
std::string PrintMemoryUsage() {
  std::ostringstream mem;
  std::ifstream proc("/proc/self/status");
  std::string s;
  while(getline(proc, s), !proc.fail()) {
    if(s.substr(0, 6) == "VmSize") {
      mem << s;
      return mem.str();
    }
  }
  return mem.str();
}

int main(int argc, char *argv[])
{
  int ierr = 0, i;

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  int myRank = Comm.MyPID();

  int numGlobalElements = 10000000;
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("numGlobalElements",&numGlobalElements,"Global problem size.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  
  Epetra_Map Map(numGlobalElements, 0, Comm);
  
  int NumMyElements = Map.NumMyElements();
  
  std::vector<int> MyGlobalElements(NumMyElements);
  Map.MyGlobalElements(&MyGlobalElements[0]);
  
  std::vector<int> NumNz(NumMyElements);
  
  for (i=0; i<NumMyElements; i++)
    //     if (MyGlobalElements[i]==0 || MyGlobalElements[i] == numGlobalElements-1)
    //       NumNz[i] = 2;
    //     else
    NumNz[i] = 3;
  
  std::cout << myRank << ": " << "Inital memory usage: " << PrintMemoryUsage() << std::endl;

  Epetra_CrsMatrix A(Copy, Map, &NumNz[0]);
  
  std::cout << myRank << ": " << "Memory after CrsMatrix constructor: " << PrintMemoryUsage() << std::endl;

  std::vector<double> Values(2);
  Values[0] = -1.0; Values[1] = -1.0;
  std::vector<int> Indices(2);
  double two = 2.0;
  int NumEntries;
  
  for (i=0; i<NumMyElements; i++)
    {
    if (MyGlobalElements[i]==0)
      {
	Indices[0] = 1;
	NumEntries = 1;
      }
    else if (MyGlobalElements[i] == numGlobalElements-1)
      {
	Indices[0] = numGlobalElements-2;
	NumEntries = 1;
      }
    else
      {
	Indices[0] = MyGlobalElements[i]-1;
	Indices[1] = MyGlobalElements[i]+1;
	NumEntries = 2;
      }
     ierr = A.InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
     assert(ierr==0);
     // Put in the diagonal entry
     ierr = A.InsertGlobalValues(MyGlobalElements[i], 1, &two, &MyGlobalElements[i]);
     assert(ierr==0);
    }
   
  std::cout << myRank << ": " << "Memory after InsertGlobalValues(): " << PrintMemoryUsage() << std::endl;

  ierr = A.FillComplete();
  assert(ierr==0);

  std::cout << myRank << ": " << "Memory after FillComplete(): " << PrintMemoryUsage() << std::endl;

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

  return ierr;
}
