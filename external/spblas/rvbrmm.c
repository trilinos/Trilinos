#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "spblas.h"
#include "tvbrmm.h"
/*#define DEBUG*/

int main(int argc, char *argv[])
{
  int    j, n, k, neqns, nvars, nnz, nnzmx, bnnz, nrhs, neqns_nrhs, ierr;
  int    nrhs_start, nrhs_stop;
  int    xstride, bstride;
  int    *indx, *bindx, *bpntrb, *bpntre, *cpntr, *rpntr, nb;
  int    izero = 0, ione = 1, i, min_calls;
  char title[80];
  SPBLASMAT  A, A1;

  char *read_file_name;
  double *val, *x, *b_reference, *b_computed;
  double factor = 2.0, dneg_one = - 1.0, done = 1.0, diff;
  double tstart, tstop, total_ops, mflops_ref, mflops;
  double alpha, beta;
  double    min_ops = 25000000; /* This value roughly determines the total
				   number of flops executed between calls
				   to the timer.  Fine grain timers may not
				   need a large min_ops value. */

  cblas_dusmm_read(argv[1], &n, &nrhs_stop, &k, &alpha, &A1, &x, &xstride, &beta,
		   &b_reference, &bstride); 

  b_computed  = (double *) calloc(bstride*nrhs_stop,sizeof(double));
   
  /* Create Sparse Blas Matrix Handle for VBR matrix */
  
  val = A1.val;
  indx = A1.indx;
  bindx = A1.bindx;
  rpntr = A1.rpntr;
  cpntr = A1.cpntr;
  bpntrb = A1.bpntrb;
  bpntre = A1.bpntre;

  cblas_duscr_vbr(n, val, indx, bindx, rpntr, cpntr, bpntrb, bpntre, &A);

  for (nrhs=1; nrhs <= nrhs_stop; nrhs++)
    {
      total_ops = (double)nrhs_stop;
      total_ops = total_ops*A.nops_per_rhs;
      
      /*      min_calls = MAX((int)(min_ops/total_ops),1); */

      min_calls = 1;

      total_ops *= (double)min_calls;
      
      /* Call reference version of MM routine */
      
      _TIMER(&tstart);
      for (i=0; i<min_calls; i++)
	cblas_dusmm_ref(n, nrhs, n, alpha, &A, x, xstride, beta, 
			b_reference, bstride);
      _TIMER(&tstop);
      mflops_ref = (total_ops * 1.0E-6)/(tstop-tstart);
      /*printf("Time for reference = %8.3g\n\n",tstop-tstart);*/
      
      
      /* Call optimized  version of MM routine */
      
      _TIMER(&tstart);
      for (i=0; i<min_calls; i++)
	cblas_dusmm(n, nrhs, n, alpha, &A, x, xstride, beta, b_computed, bstride);
      _TIMER(&tstop);
      mflops = (total_ops * 1.0E-6)/(tstop-tstart);
      /*printf("Time for optimized = %8.3g\n\n",tstop-tstart);*/
      
      /* Compare results */
      neqns = bstride;
      neqns_nrhs = neqns*nrhs;
      F77NAME(daxpy) 
	(&neqns_nrhs, &dneg_one, b_computed, &ione, b_reference, &ione);
      
      diff = F77NAME(dnrm2) (&neqns_nrhs, b_reference, &ione);

      nvars = xstride;
      nnz = indx[bpntre[n-1]];
      bnnz = bpntre[n-1];

      printf("Number of equations       = %d\n",neqns);
      printf("Number of variables       = %d\n",nvars);
      printf("Number of block equations = %d\n",n);
      printf("Number of nonzeros        = %d\n",nnz);
      printf("Number of block nonzeros  = %d\n",bnnz);
      printf("Reference MFLOPS          = %8.3g\n",mflops_ref);
      printf("Computed MFLOPS           = %8.3g\n",mflops_ref);
      printf("Difference in solutions   = %4.1e\n",diff);
#ifdef DEBUG   
      /* Print out computed RHS (for debugging only) */
      for (i=0; i<neqns_nrhs; i++)
	printf("b[%d] = %e\n",i+1,b_computed[i]); 
#endif

    } /* End of nrhs loop */
  return 0 ;
}
