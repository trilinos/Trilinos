#include <stdlib.h>
#include <stdio.h>
#include "spblas.h"

void cblas_dusmm1(int m, int nrhs, int k, double alpha, SPBLASMAT *A,
		 double *x, int xstride, double beta, double *b, int bstride, int istrt, int istop)
 
/*  Compute sparse matrix time dense matrix multiply. Only works for VBR now.
*/

{
  int i, j, irhs;
  int *indx,  *bindx, * rpntr,  *cpntr, *bpntrb, *bpntre;
  int nrow, ncol;
  double *val, *xptr, *Aptr, *bptr, d_one = 1.0;

  val = A->val;
  indx = A->indx;
  bindx = A->bindx;
  rpntr = A->rpntr;
  cpntr = A->cpntr;
  bpntrb = A->bpntrb;
  bpntre = A->bpntre;

  printf(" ISTRT = %d, ISTOP = %d\n",istrt,istop);

  /* Compute Matrix multiply block entry by block entry using DGEMM */

  for (i=istrt; i<istop; i++)
    {
      printf("Computing row = %d\n",i);
      nrow = rpntr[i+1] - rpntr[i];
      bptr = b+rpntr[i];

      if (beta == 0.0)
	for (irhs=0; irhs<nrhs; irhs++)
	  for (j=0; j<nrow; j++)
	    bptr[j+irhs*bstride] = 0.0;
      else
	for (irhs=0; irhs<nrhs; irhs++)
	  for (j=0; j<nrow; j++)
	    bptr[j+irhs*bstride] *= beta;

      for (j=bpntrb[i]; j<bpntre[i]; j++)
	{
	  Aptr = val+indx[j];
	  xptr = x+cpntr[bindx[j]];
	  ncol = cpntr[bindx[j]+1] - cpntr[bindx[j]];
	  F77NAME(dgemm)
	    ("N", "N", &nrow, &nrhs, &ncol, &alpha, Aptr, 
	     &nrow, xptr, &xstride,
	     &d_one, bptr, &bstride);
	}
    }
	  




/* end cblas_duscr_vbr */
}
