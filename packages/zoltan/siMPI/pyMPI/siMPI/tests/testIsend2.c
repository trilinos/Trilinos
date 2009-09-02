
/*TEST
PATH='tests/testIsend2.c'
CCFLAGS=""
INPUT=""
OUTPUT='request->valid = 500\n'
STATUS=0
TEST*/

/******************************************************************/
/* FILE  ***********      testIsend2.c         ********************/
/******************************************************************/
/* Author : Lisa Alano July 29 2002                               */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/
/******************************************************************/

#if 0
CCFLAGS = None 
ARGS = None
INPUT = EOF 
OUTPUT = request-> = 0 
  
STATUS = 0 
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

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  source = tag = dest = 0;
  sprintf(message1, "Hello there");
  MPI_Irecv(message2, 50, MPI_CHAR, source, tag, MPI_COMM_WORLD, &request);
  printf("request->valid = %d\n",request->valid);
   
  /*printf("%s\n", message2);*/ /*This will be garbage*/

  MPI_Finalize();
  return 0;
}
