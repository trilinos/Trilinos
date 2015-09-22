
/*TEST
PATH='tests/testSendRecv0.c'
CCFLAGS=""
INPUT=""
OUTPUT='Hello there\n'
STATUS=0
TEST*/

/******************************************************************/
/* FILE  ***********      testSendRecv0.c      ********************/
/******************************************************************/
/* Author : Lisa Alano June 5 2002                                */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/
/******************************************************************/

#if 0
CCFLAGS = None 
ARGS = None
INPUT = EOF 
OUTPUT = Hello there
STATUS = 0 
#endif

#include <stdio.h>
#include <string.h>
#include "mpi.h"

int main(int argc, char**argv) 
{
  int my_rank;
  int p;
  char message1[50];
  char message2[50];
  int source, dest, tag; 
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  source = tag = dest = 0;
  sprintf(message1, "Hello there");
  MPI_Send(message1, strlen(message1)+1, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
  MPI_Recv(message2, 50, MPI_CHAR, source, tag, MPI_COMM_WORLD, &status);
  
  printf("%s\n", message2);

  MPI_Finalize();
  return 0;
}
