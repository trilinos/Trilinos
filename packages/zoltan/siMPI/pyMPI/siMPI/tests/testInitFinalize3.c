
/*TEST
PATH='tests/testInitFinalize3.c'
CCFLAGS=""
INPUT=""
OUTPUT=''
STATUS=0
TEST*/

/******************************************************************/
/* FILE  ***********    testInitFinalize3.c    ********************/
/******************************************************************/
/* Author : Lisa ALano June 7 2002                                */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/
/******************************************************************/

#include "mpi.h"

int main(int argc, char**argv) 
{
  MPI_Init(&argc, &argv);
  return 0;
}
