
/*TEST
PATH='tests/testSendRecv4.c'
CCFLAGS=""
INPUT=""
OUTPUT=''
STATUS=0
TEST*/

/******************************************************************/
/* FILE  ***********      testSendRecv4.c      ********************/
/******************************************************************/
/* Author : Lisa Alano June 5 2002                                */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/
/******************************************************************/

#if 0
CCFLAGS = None 
ARGS = None 
INPUT = EOF 
OUTPUT = None 
       There is no matching receive for the send
STATUS = 0 
#endif

#include <stdio.h>
#include <string.h>
#include "mpi.h"

int main(int argc, char**argv) 
{
  int my_rank, p;
  char message1[50];
  int source, dest, tag; 

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);


  source = tag = dest = 0;
  sprintf(message1, "Hello there");
  MPI_Send(message1, strlen(message1)+1, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
  
  MPI_Finalize();
  return 0;
}
