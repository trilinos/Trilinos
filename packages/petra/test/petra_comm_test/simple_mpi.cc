// Petra_Comm Test routine

#include <math.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <iostream.h>
#include <mpi.h>

int main(int argc, char *argv[]) {

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  cout << "Process "<<rank<<" of "<<size<<" is alive." <<endl;

  MPI_Finalize();
  return 0;
}
