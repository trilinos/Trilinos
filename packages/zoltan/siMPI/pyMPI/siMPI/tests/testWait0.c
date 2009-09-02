
/*TEST
PATH='tests/testWait0.c'
CCFLAGS=""
INPUT=""
OUTPUT='Hello there\n'
STATUS=0
TEST*/

/******************************************************************/
/* FILE  ***********      testWait0.c         ********************/
/******************************************************************/
/* Author : Lisa Alano July 29 2002                               */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/
/******************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mpi.h"

int main(int argc, char**argv) 
{
  int my_rank;
  int p;
  char message1[50];
  char message2[50];
  int source, dest, tag; 
  MPI_Request request1;
  MPI_Request request2;
  MPI_Status* status;
  status = (MPI_Status *) malloc (sizeof(MPI_Status));

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  source = tag = dest = 0;
  sprintf(message1, "Hello there");
  MPI_Isend(message1, strlen(message1)+1, MPI_CHAR, dest, tag, MPI_COMM_WORLD, &request1);
  MPI_Irecv(message2, 50, MPI_CHAR, source, tag, MPI_COMM_WORLD, &request2);
  MPI_Wait (&request1, status); 
  MPI_Wait (&request2, status); 
  printf("%s\n", message2);

  MPI_Finalize();
  return 0;
}
