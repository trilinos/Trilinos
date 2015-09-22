
/*TEST
PATH='tests/testIsend1.c'
CCFLAGS=""
INPUT=""
OUTPUT='First-> <no value yet>\nSecond-> Hello there\n'
STATUS=0
TEST*/

/******************************************************************/
/* FILE  ***********      testIsend1.c         ********************/
/******************************************************************/
/* Author : Lisa Alano July 29 2002                               */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/
/******************************************************************/

#if 0
CCFLAGS = None 
ARGS = None
INPUT = EOF 
STATUS = 0 
OUTPUT =
First-> <no value yet>
Second-> Hello there 
#endif

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
  MPI_Request request;
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  source = tag = dest = 0;
  sprintf(message1, "Hello there");
  sprintf(message2, "<no value yet>");
  MPI_Irecv(message2, 50, MPI_CHAR, source, tag, MPI_COMM_WORLD, &request);
  printf("First-> %s\n", message2);
  MPI_Isend(message1, strlen(message1)+1, MPI_CHAR, dest, tag, MPI_COMM_WORLD, &request);
  MPI_Irecv(message2, 50, MPI_CHAR, source, tag, MPI_COMM_WORLD, &request);
  MPI_Wait(&request,&status);
  printf("Second-> %s\n", message2);


  MPI_Finalize();
  return 0;
}
