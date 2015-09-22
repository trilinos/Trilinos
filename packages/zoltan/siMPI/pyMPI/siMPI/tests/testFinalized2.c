
/*TEST
PATH='tests/testFinalized2.c'
CCFLAGS=""
INPUT=""
OUTPUT="ERROR: -5016 MPI_FINALIZED: Invalid pointer.\nERROR: -5015 MPI aborting...\n"
STATUS=1
TEST*/

/******************************************************************/
/* FILE  ***********     testFinalized2.c      ********************/
/******************************************************************/
/* Author : Lisa Alano July 25 2002                               */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/
/******************************************************************/

#include <stdio.h>
#include "mpi.h"

int main(int argc, char**argv) 
{
  MPI_Init(&argc, &argv);
  MPI_Finalized(NULL);
  MPI_Finalize();
  return 0;
}
