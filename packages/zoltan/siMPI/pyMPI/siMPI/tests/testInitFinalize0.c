
/*TEST
PATH='tests/testInitFinalize0.c'
CCFLAGS=""
INPUT=""
OUTPUT=''
STATUS=0
TEST*/

/******************************************************************/
/* FILE  ***********    testInitFinalize0.c    ********************/
/******************************************************************/
/* Author : Lisa ALano June 5 2002                                */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/
/******************************************************************/

#if 0
CCFLAGS = None 
ARGS = None
INPUT = EOF 
OUTPUT = None
STATUS = 0 
#endif


#include "mpi.h"

int main(int argc, char**argv) 
{
  MPI_Init(&argc, &argv);
  MPI_Finalize();
  return 0;
}
