/*TEST
CCFLAGS=""
OUTPUT="""siMPI running on rank 0
"""
STATUS=0
TEST*/

#include "mpi.h"

int main(int argc, char** argv) {
  int rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  printf("siMPI running on rank %d\n",rank);
  MPI_Finalize();
  return 0;
}
