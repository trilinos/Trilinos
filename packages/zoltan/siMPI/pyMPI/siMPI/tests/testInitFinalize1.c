
/*TEST
PATH='tests/testInitFinalize1.c'
CCFLAGS=""
INPUT=""
OUTPUT="ERROR: -5016 MPI was already finalized\nERROR: -5016 Error with Init status.\nERROR: -5016 MPI aborting..."
STATUS=1
TEST*/

/******************************************************************/
/* FILE  ***********     testInitFinalize1.c   ********************/
/******************************************************************/
/* Author : Lisa ALano June 7 2002                                */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/
/******************************************************************/

#include "mpi.h"

int main(int argc, char**argv) 
{
  MPI_Init(&argc, &argv);
  MPI_Finalize();
  MPI_Init(&argc, &argv);
  return 0;
}
