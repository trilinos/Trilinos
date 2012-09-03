/*--------------------------------------------------------------------*/
#include "Epetra_FEVector.h"

#ifdef HAVE_MPI

#include <time.h>
#include "mpi.h"
#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"

int main(int argCount, char **argValue)
{
  int ierr;
  MPI_Init(&argCount,&argValue);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  const int rank = Comm.MyPID();

  // Construct a Map 
  int nGlobalElements = 1000000;
  Epetra_Map Map(nGlobalElements, 0, Comm);

  // Create a vector
  Epetra_FEVector b(Map, 1);

  time_t startTime = 0;
  if (rank == 0) {
    startTime = time(0);
  }

  // Fill matrix on the master process
  if (rank == 0) {
    double values[1];
    int    indices[1];

    for (int globalRowIdx=0; globalRowIdx<nGlobalElements; ++globalRowIdx) {
      indices[0] = globalRowIdx;
      values[0] = 3.2 + globalRowIdx*0.01;

      if (globalRowIdx % 10000 == 0) {
  cerr << "About to insert row " << globalRowIdx << "\n";
      }

      ierr = b.ReplaceGlobalValues(1, (const int *)&indices[0], 
           (const double *)&values[0]);
      assert(ierr==0);
    }
  }

  double insertionTime = 0;
  if (rank == 0) {
    time_t endTime = time(0);
    insertionTime = difftime(endTime, startTime);
  }

  // Finish up
  ierr = b.GlobalAssemble();
  assert(ierr==0);

  if (rank == 0) {
    cerr << "insertion time = " << insertionTime << " (seconds)\n";
  }


  MPI_Finalize();

  return 0;
}
#else
int main(int,char**)
{
  return 0;
}
#endif
/*--------------------------------------------------------------------*/
