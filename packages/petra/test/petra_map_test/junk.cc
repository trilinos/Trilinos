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

  Petra_Map * map = new Petra_Map(100, 0, MPI_COMM_WORLD);

  cout << *map << endl;

#ifdef PETRA_MPI
  MPI_Finalize();
#endif

  return 0;
}

