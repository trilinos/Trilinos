#include <stdlib.h>
#include <stdio.h>
#include "spblas.h"

void cblas_dusmm_read(char *read_file_name, 
		      int *m, int *nrhs, int *k, double *alpha, SPBLASMAT *A,
		      double **x, int *xstride, double *beta, double **b, int *bstride)
 
/*  Compute sparse matrix time dense matrix multiply. Only works for VBR now.
*/

{
  FILE *read_file ;
  int i, j, jj, irhs;
  int *indx,  *bindx, *rpntr,  *cpntr, *bpntrb, *bpntre;
  int nrow, ncol, n_blk_col;
  int itmp;
  double dtmp;
  double *val;

  /* Open read file */

  read_file = fopen(read_file_name, "r");
  

  /* Read out parameters */

  fscanf(read_file,"%d %d %d %d %d %d",m, &n_blk_col, nrhs, k, xstride, bstride);
  fscanf(read_file,"%lf %lf",alpha, beta);

  /* Read in each array */
  rpntr = calloc(*m+1,sizeof(int));
  cpntr = calloc(n_blk_col+1,sizeof(int));
  bpntrb = calloc(*m+1,sizeof(int));
  bpntre = calloc(*m+1,sizeof(int));

  for (i=0; i<*m+1; i++)
    {
      fscanf(read_file,"%d",&itmp);
      rpntr[i] = itmp;
    }

  for (i=0; i<n_blk_col+1; i++)
    {
      fscanf(read_file,"%d",&itmp);
      cpntr[i] = itmp;
    }

  for (i=0; i<*m+1; i++)
    {
      fscanf(read_file,"%d",&itmp);
      bpntrb[i] = itmp;
    }

  for (i=0; i<*m+1; i++)
    {
      fscanf(read_file,"%d",&itmp);
      bpntre[i] = itmp;
    }

  bindx = calloc(bpntre[*m-1]+1,sizeof(int));
  indx  = calloc(bpntre[*m-1]+1,sizeof(int));

  for (i=0; i<*m; i++)
      for (j=bpntrb[i]; j<bpntre[i]; j++)
	{
	  fscanf(read_file,"%d",&itmp);
	  bindx[j] = itmp;
	}

  for (i=0; i<*m; i++)
      for (j=bpntrb[i]; j<bpntre[i]; j++)
	{
	  fscanf(read_file,"%d",&itmp);
	  indx[j] = itmp;
	}

  /* indx is one element longer, read it too */

  fscanf(read_file,"%d",&itmp);
  indx[bpntre[*m-1]] = itmp;

  val  = calloc(itmp+1,sizeof(double));

  for (i=0; i<*m; i++)
      for (j=bpntrb[i]; j<bpntre[i]; j++)
	for (jj = indx[j]; jj<indx[j+1]; jj++)
	  {
	    fscanf(read_file,"%lf",&dtmp);
	    val[jj] = dtmp;
	  }

  *x  = calloc(*nrhs*(*xstride),sizeof(double));
  *b  = calloc(*nrhs*(*bstride),sizeof(double));

  for (irhs=0; irhs<*nrhs; irhs++)
    for (j=irhs*(*xstride); j<*xstride+irhs*(*xstride); j++)
      {
	fscanf(read_file,"%lf",&dtmp);
	(*x)[j] = dtmp;
      }

  for (irhs=0; irhs<*nrhs; irhs++)
    for (j=irhs*(*bstride); j<*bstride+irhs*(*bstride); j++)
      {
	fscanf(read_file,"%lf",&dtmp);
	(*b)[j] = dtmp;
      }

  fclose(read_file);

  A->val = val ;
  A->indx = indx;
  A->bindx = bindx;
  A->rpntr = rpntr;
  A->cpntr = cpntr;
  A->bpntrb = bpntrb;
  A->bpntre = bpntre;

/* end cblas_duscr_read */
}
