
/*TEST
SKIP=1
PATH='tests/testRankSize7.c'
CCFLAGS=""
INPUT=""
OUTPUT=''
STATUS=0
TEST*/

/******************************************************************/
/* FILE  ***********    testRankSize7.c	       ********************/
/******************************************************************/
/* Author : Lisa Alano June 5 2002                                */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/
/******************************************************************/

#if 0
CCFLAGS = None 
ARGS = None 
INPUT = EOF 
OUTPUT = No MPI_Init call. 
STATUS = 1 
#endif

#include "mpi.h"

int main(int argc, char**argv) 
{
  int my_rank;
  int p;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Finalize();
  return 0;
}
