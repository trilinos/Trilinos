
/*TEST
PATH='tests/testReduce2.c'
CCFLAGS=""
INPUT=""
OUTPUT='myNumber = 1010 and total = 222\nmyNumber = 1010 and total = 1010\nretval = 0\n'
STATUS=0
TEST*/

/******************************************************************/
/* FILE  ***********      testReduce3.c        ********************/
/******************************************************************/
/* Author : Lisa ALano July 16 2002                               */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/

#if 0
CCFLAGS = TEST FLAGS 
ARGS = None
INPUT = EOF 
OUTPUT =  myNumber = 1010 and total = 222
  myNumber = 1010 and total = 1010
  retval = 0 
STATUS = 0 
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

  retval = MPI_Reduce (&myNumber, &total, 1, MPI_INT, MPI_MAX, my_rank, MPI_COMM_WORLD);
  printf("myNumber = %d and total = %d\n", myNumber, total);
  printf("retval = %d\n", retval);

  MPI_Finalize();
  return 0;
}
