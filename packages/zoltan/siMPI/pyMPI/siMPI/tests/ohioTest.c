
/*TEST
PATH='tests/ohioTest.c'
CCFLAGS=""
INPUT=""
OUTPUT='PE:0 result is 1.000000 + 2.000000i\n'
STATUS=0
SKIP=1
TEST*/

#include <stdio.h>
#include "mpi.h"

/* Ohio Supercomputer Center */ 

/* Modified to single rank -- pjm */


typedef struct {
   double real,imag;
} complex;

void cprod(complex *in, complex *inout, int *len, MPI_Datatype *dptr) {
   int i;
   complex c;

   for (i=0; i<*len; ++i) {
      c.real=(*in).real * (*inout).real - (*in).imag * (*inout).imag;
      c.imag=(*in).real * (*inout).imag + (*in).imag * (*inout).real;
      *inout=c;
      in++;
      inout++;
   }
}

int main (int argc, char *argv[]) {
   int rank;
   int root;
   complex source,result;

   MPI_Op myop;
   MPI_Datatype ctype;
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   MPI_Type_contiguous(2,MPI_DOUBLE,&ctype);
   MPI_Type_commit(&ctype);
   MPI_Op_create(cprod,1,&myop);
   root=0;
   source.real=rank+1;
   source.imag=rank+2;
   result.real = 0;
   result.imag = 0;
   MPI_Reduce(&source,&result,1,ctype,myop,root,MPI_COMM_WORLD);
   if(rank==root) printf ("PE:%d result is %lf + %lfi\n",rank, result.real, result.imag);
   MPI_Finalize();
}

/* PE:0 result is 1.000000 + 2.000000i */
