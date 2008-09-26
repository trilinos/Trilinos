
/*TEST
PATH='tests/testRankSize1.c'
CCFLAGS=""
INPUT=""
OUTPUT='My_rank is 0.\nThe number of processes is 1.\n'
STATUS=0
TEST*/

/******************************************************************/
/* FILE  ***********     testRankSize1.c       ********************/
/******************************************************************/
/* Author : Lisa Alano June 5 2002                                */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/
/******************************************************************/

#if 0
CCFLAGS = None 
ARGS = None
INPUT = EOF 
OUTPUT = Error with my_rank (double instead of int).
STATUS = 1
#endif


#include <stdio.h>
#include <string.h>
#include "mpi.h"

int main(int argc, char**argv) 
{
  int my_rank;
  int p;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  printf("My_rank is %d.\n",my_rank);
  printf("The number of processes is %d.\n",p);

  MPI_Finalize();
  return 0;
}
