
/*TEST
SKIP=1
PATH='tests/testReduce1.c'
CCFLAGS=""
INPUT=""
OUTPUT=''
STATUS=0
TEST*/

/******************************************************************/
/* FILE  ***********      testReduce1.c        ********************/
/******************************************************************/
/* Author : Lisa ALano July 16 2002                               */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/

#if 0
CCFLAGS = TEST FLAGS 
ARGS = None
INPUT = EOF 
OUTPUT = myNumber = 1010 and total = 222
  ERROR: -5001 MPI_REDUCE : Invalid buffer pointer: Arguments must specify different buffers (no aliasing)
  ERROR: -5001 MPI aborting...
STATUS = 1 
#endif

#include <stdio.h>
#include <string.h>
#include "mpi.h"

int main(int argc, char**argv) 
{
  int my_rank;
  int p;
  int myNumber, total;
  int retval;

  myNumber = 1010;
  total = 222;
  printf("myNumber = %d and total = %d\n", myNumber, total);

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  retval = MPI_Reduce (&myNumber, &myNumber, 1, MPI_INT, MPI_SUM, my_rank, MPI_COMM_WORLD);
  printf("myNumber = %d and total = %d\n", myNumber, total);
  printf("retval = %d\n", retval);

  MPI_Finalize();
  return 0;
}
