
/*TEST
SKIP=1
PATH='tests/testTypeContig0.c'
CCFLAGS=""
INPUT=""
OUTPUT=''
STATUS=0
TEST*/

#include "mpi.h"
#include <stdio.h>
#define SIZE 4

int main(argc,argv)
int argc;
char *argv[];  {
int numtasks, rank, source=0, dest, tag=1, i;
float a[SIZE] =
  {1.0, 2.0, 3.0, 4.0};
float b[SIZE];

MPI_Status stat;
MPI_Datatype rowtype;

MPI_Init(&argc,&argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

MPI_Type_contiguous(SIZE, MPI_FLOAT, &rowtype);
MPI_Type_commit(&rowtype);

if (rank == 0) {
  for (i=0; i<numtasks; i++)
    MPI_Send(&a[i], 1, rowtype, i, tag, MPI_COMM_WORLD);
}

MPI_Recv(b, SIZE, MPI_FLOAT, source, tag, MPI_COMM_WORLD, &stat);
printf("rank= %d  b= %3.1f %3.1f %3.1f %3.1f\n",
         rank,b[0],b[1],b[2],b[3]);

MPI_Finalize();
}

