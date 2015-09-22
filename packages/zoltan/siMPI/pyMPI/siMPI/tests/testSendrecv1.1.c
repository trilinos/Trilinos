
/*TEST
SKIP=1
PATH='tests/testSendrecv1.1.c'
CCFLAGS=""
INPUT=""
OUTPUT=''
STATUS=0
TEST*/

/******************************************************************/
/* FILE  ***********      testSendrecv1.c      ********************/
/******************************************************************/
/* Author : Lisa Alano July 29 2002                               */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/
/******************************************************************/

#if 0
CCFLAGS = None 
ARGS = None
INPUT = EOF 
OUTPUT = ERROR: -5001 MPI_RECV argument error
  ERROR: -5001 MPI aborting... 
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
  MPI_Sendrecv(message1, strlen(message1)+1, MPI_CHAR, dest, tag, NULL, 50, MPI_CHAR, source, tag, MPI_COMM_WORLD, &status);
  
  printf("%s\n", message2);

  MPI_Finalize();
  return 0;
}
