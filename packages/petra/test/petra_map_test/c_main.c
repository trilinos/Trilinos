/*
  file main.c
*/
#include <stdio.h>
#include "mpi.h"
#include "Petra_c_wrappers.h"

int main(int argc, char *argv[]) {
  int mypid;
  PETRA_COMM petra_comm;
  /*
     Start executable statements.
  */
  MPI_Init(&argc,&argv);
  petra_comm = petra_comm_create( MPI_COMM_WORLD );
  /*
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  */
  mypid = petra_comm_getmypid(petra_comm);
  printf("MyPID = %d\n",mypid);
  MPI_Finalize();

  return 0;
}

/*
  end of file main.c
*/
