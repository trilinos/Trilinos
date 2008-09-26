/*TEST
PATH='tests/testFinalized2.c'
CCFLAGS=""
INPUT=""
OUTPUT="ERROR: -5016 MPI_FINALIZE: MPI has not been initialized.\nERROR: -5016 MPI aborting...\n"
STATUS=1
TEST*/


/******************************************************************/
/* FILE  ***********    testInitFinalize2.c    ********************/
/******************************************************************/
/* Author : Lisa ALano June 7 2002                                */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/
/******************************************************************/

#include "mpi.h"

int main(int argc, char**argv) 
{
  MPI_Finalize();
  return 0;
}
