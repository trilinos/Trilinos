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
  int    j, n, neqns, nnz, nnzmx, bnnz, nrhs, neqns_nrhs, ierr;
  int    nrhs_start, nrhs_stop;
  int    *indx, *bindx, *bpntr, *cpntr, nb, nx, ny, nz;
  int    izero = 0, ione = 1, i, min_calls;
  char title[80];
  SPBLASMAT  A;

  double *val, *x, *b_reference, *b_computed;
  double factor = 2.0, dneg_one = - 1.0, done = 1.0, diff;
  double tstart, tstop, total_ops, mflops_ref, mflops;
  double    min_ops = 25000000; /* This value roughly determines the total
				   number of flops executed between calls
				   to the timer.  Fine grain timers may not
				   need a large min_ops value. */

  if (argc < 2 || (argc < 3 && argv[1][0] == '-' && argv[1][1] == 'h'))
    { 
      printf("\n%s nb nrhs nx [ny [nz]]\n\n",argv[0]);
      printf("     -H - Print header for output (optional)\n");
      printf("     nb - Block entry size\n");
      printf("     nrhs - Number of RHS.  If nrhs = A, loop from nrhs=1-15\n");
      printf("     nx,ny,nz - Nodes in X,Y,Z\n\n");
      printf("If ny and nz undefined, nz=ny=nx is assumed\n");
      printf("If nz undefined, nz=ny is assumed\n\n");
      printf(
	     "Example: Prints header, then computes and prints results for two cases. \n\n");
      printf("  %s -H \n",argv[0]);
      printf("  %s 4 4 4 \n",argv[0]);
      printf("  %s 4 8 4 \n",argv[0]);

      exit(1);
    }
  if (argv[1][0] == '-')
    {
      printf(
	     "Block | # of | Dim. of   |# of  | # of   |Total # | Total #  | MFLOPS  | MFLOPS   |        | Total\n");
      printf(
	     "Size  | RHS  | XYZ cube  |Bl Eq |BlEntry |eqns    | entries  |Reference|Optimized | Error  | Op/call\n");
      printf(
	     "------|------|-----------|------|--------|--------|----------|---------|----------|--------|--------\n");
      exit(0);
    }
  nb = atoi(argv[1]); /* Convert arg to integer */
  if (nb < 0)
    {
      printf("nb must be greater than zero\n");
      exit(1);
    }

  if (argv[2][0] == 'A') /* If nrhs = A, then do a range of 1 <= nrhs <= 15 */
    {
      nrhs_start = 1;
      nrhs_stop  = 15;
      nrhs = nrhs_stop;
    }
  else
    {
      nrhs = atoi(argv[2]); /* Convert arg to integer */
      if (nrhs < 0)
	{
	  printf("nrhs must be greater than zero\n");
	  exit(1);
	}
      nrhs_start = nrhs;
      nrhs_stop  = nrhs;
    }
  
  if (argc < 5)
    {
      nx = atoi(argv[3]); /* Convert arg to integer */
      if (nx < 0)
	{
	  printf("nx must be greater than zero\n");
	  exit(1);
	}
      else
	{
	  nz = nx; ny = nx;
	}
    }
  else if (argc < 6)
    {
      nx = atoi(argv[3]); /* Convert arg to integer */
      ny = atoi(argv[4]); /* Convert arg to integer */
      if (nx < 0 || ny < 0)
	{
	  printf("nx and ny must be greater than zero\n");
	  exit(1);
	}
      else
	{
	  nz = ny;
	}
    }
  else
    {
      nx = atoi(argv[3]); /* Convert arg to integer */
      ny = atoi(argv[4]); /* Convert arg to integer */
      nz = atoi(argv[5]); /* Convert arg to integer */
    }

  /* Allocate space for matrix */
  n = nx * ny * nz;
  neqns = n * nb;
  nnzmx = 27 * n;
  bnnz = nb * nb * nnzmx;
  neqns_nrhs = neqns*nrhs;
   
  bpntr       = (int *) calloc(n+1,sizeof(int));
  cpntr       = (int *) calloc(n+1,sizeof(int));
  x           = (double *) calloc(neqns_nrhs,sizeof(double));
  b_reference = (double *) calloc(neqns_nrhs,sizeof(double));
  b_computed  = (double *) calloc(neqns_nrhs,sizeof(double));
   
  indx        = (int *) calloc(nnzmx,sizeof(int));
  bindx       = (int *) calloc(nnzmx,sizeof(int));
   
  val         = (double *) calloc(bnnz,sizeof(double));
  if (val == NULL) perror("Error: Not enough space to create matrix");
   
  F77NAME(d27ptgen)
    (&nnzmx, &nx, &ny, &nz, &nb, &factor, &izero, title, &n, 
     val, bpntr, bindx, &nnz, &ierr);
  if (ierr != 0)
    {
      printf("Error in generating matrix. Ran out of memory.\n");
      exit(1);
    }
#ifdef DEBUG   
  /* Print out val (for debugging only) */
    for (i=0; i<nb*nb*nnz; i++)
      printf("val[%d] = %e\n",i+1,val[i]);
#endif

  /* Generate partition vector */
  cpntr[0] = 0;
  for (i=0; i<n; i++) cpntr[i+1] = cpntr[i]+nb;

  /* Set pointer and index vectors to zero base */
  for (i=0; i<n+1; i++) bpntr[i]--;
  for (i=0; i<nnz; i++) bindx[i]--;

#ifdef DEBUG   
  /* set x to increasing integer  numbers (use for debugging) */

    for (i=0; i<neqns; i++)
      x[i] = i+1; 
    for (j=1; j<nrhs; j++)
      for (i=0;i<neqns;i++)
	x[i+j*neqns] = x[i];

  /* Print out LHS (for debugging only) */
    for (i=0; i<neqns_nrhs; i++)
      printf("x[%d] = %e\n",i+1,x[i]);

#else

  /* set x to random numbers */

  for (i=0; i<neqns_nrhs; i++)
    x[i] = drand48();

#endif

  /* Set values in indx */
  indx[0] = bpntr[0];
  for (i=0; i< n; i++)
    for (j=bpntr[i]; j<bpntr[i+1]; j++)
      indx[j+1] = indx[j]+(cpntr[i+1]-cpntr[i])*(cpntr[bindx[j]+1]-cpntr[bindx[j]]);

  /* Create Sparse Blas Matrix Handle for VBR matrix */
  
  cblas_duscr_vbr(n, val, indx, bindx, cpntr, cpntr, bpntr, &(bpntr[1]), &A);

  for (nrhs=nrhs_start; nrhs <= nrhs_stop; nrhs++)
    {
      total_ops = (double)nrhs;
      total_ops = total_ops*A.nops_per_rhs;
      
      min_calls = MAX((int)(min_ops/total_ops),1);
      total_ops *= (double)min_calls;
      
      /* Call reference version of MM routine */
      
      _TIMER(&tstart);
      for (i=0; i<min_calls; i++)
	cblas_dusmm_ref(n, nrhs, n, 1.0, &A, x, neqns, 0.0, b_reference, neqns);
      _TIMER(&tstop);
      mflops_ref = (total_ops * 1.0E-6)/(tstop-tstart);
      /*printf("Time for reference = %8.3g\n\n",tstop-tstart);*/
      
      
      /* Call optimized  version of MM routine */
      
      _TIMER(&tstart);
      for (i=0; i<min_calls; i++)
	cblas_dusmm(n, nrhs, n, 1.0, &A, x, neqns, 0.0, b_computed, neqns);
      _TIMER(&tstop);
      mflops = (total_ops * 1.0E-6)/(tstop-tstart);
      /*printf("Time for optimized = %8.3g\n\n",tstop-tstart);*/
      
      /* Compare results */
      F77NAME(daxpy) 
	(&neqns_nrhs, &dneg_one, b_computed, &ione, b_reference, &ione);
      
      diff = F77NAME(dnrm2) (&neqns_nrhs, b_reference, &ione);
      printf(
   " %4d | %4d |%3d,%3d,%3d|%5d |%7d |%7d |%9d |%8.3g | %8.3g | %4.1e | %4.2e\n",
   nb, nrhs, nx,ny,nz, n, nnz, neqns, nb*nnz, mflops_ref, mflops, diff, total_ops);

#ifdef DEBUG   
      /* Print out computed RHS (for debugging only) */
      for (i=0; i<neqns_nrhs; i++)
	printf("b[%d] = %e\n",i+1,b_computed[i]); 
#endif

    } /* End of nrhs loop */
  return 0 ;
}
