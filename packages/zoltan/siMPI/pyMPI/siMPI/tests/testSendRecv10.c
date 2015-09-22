
/*TEST
PATH='tests/testSendRecv10.c'
CCFLAGS=""
INPUT=""
OUTPUT='0 1 2 3 4 5 6 7 8 9 \n--\n0 1 2 3 4 5 6 7 8 9\n'
STATUS=0
TEST*/

/******************************************************************/
/* FILE  ***********    testSendRecv10.c       ********************/
/******************************************************************/
/* Author : Lisa Alano June 5 2002                                */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/
/******************************************************************/

#include <stdio.h>
#include <string.h>
#include "mpi.h"

int main(int argc, char**argv) 
{
  int my_rank;
  int p, i;
  char message1[15];
  char message2[15];
  int source, dest, tag; 
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  for(i = 0 ; i < 10 ; i++) {
      message1[i] = i;
  }
  for(i=0; i<10; i++) {
     printf("%d ", message1[i]);
  }
  printf("\n--\n");
  source = tag = dest = 0;
  MPI_Send(message1, 15, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
  MPI_Recv(message2, 15, MPI_CHAR, source, tag, MPI_COMM_WORLD, &status);

  for(i=0; i<10; i++)
  {
    printf("%d ", message2[i]);
  }
  MPI_Finalize();
  return 0;
}
