#include <stdlib.h>
#include <stdio.h>
#include "spblas.h"

void cblas_dusmm(int m, int nrhs, int k, double alpha, SPBLASMAT *A,
		 double *x, int xstride, double beta, double *b, int bstride)
 
/*  Compute sparse matrix time dense matrix multiply. Only works for VBR now.
*/

{
  int ncpu, mchunk, icpu, istrt, istop;

  ncpu = omp_get_num_procs();
  /*ncpu = 2;*/
  /*printf("Number of CPUS = %d\n",ncpu);*/
  if (ncpu > 2 || m < 2 ) abort();  /* Can't handle these cases yet */

  mchunk = MAX(1,(m+1)/ncpu);

  for (icpu = 0; icpu<ncpu; icpu++)
     {
	istrt = icpu * mchunk;
     istop = MIN(m, istrt+mchunk);

     cblas_dusmm1(m, nrhs, k, alpha, A, x, xstride, beta, b, bstride,
                  istrt, istop);
     }
}

