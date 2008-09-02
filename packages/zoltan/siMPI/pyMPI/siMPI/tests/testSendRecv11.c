
/*TEST
PATH='tests/testSendRecv11.c'
CCFLAGS=""
INPUT=""
OUTPUT='0 1 2 3 4 5 6 7 8 9 \n'
STATUS=0
TEST*/

/******************************************************************/
/* FILE  ***********     testSendRecv11.c      ********************/
/******************************************************************/
/* Author : Lisa Alano June 5 2002                                */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/
/******************************************************************/

#if 0
CCFLAGS = None 
ARGS = None
INPUT = EOF 
OUTPUT = 0 1 2 3 4 5 6 7 8 9
STATUS = 0 
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

int main(int argc, char**argv) 
{
  int my_rank;
  int p, i;
  int message1[15];
  int message2[15];
  int source, dest, tag; 
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  
  for (i=0; i<10; i++)
    message1[i] = i;

  source = tag = dest = 0;
  MPI_Send(message1, 15, MPI_INT, dest, tag, MPI_COMM_WORLD);
  MPI_Recv(message2, 15, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
 
  for (i=0; i<10; i++) 
    printf("%d ", message2[i]);
  printf("\n");

  MPI_Finalize();
  return 0;
}
