#include <stdlib.h>
#include <stdio.h>
#include "spblas.h"

void cblas_dusmm_dump(char *dump_file_name, 
		      int m, int nrhs, int k, double alpha, SPBLASMAT *A,
		      double *x, int xstride, double beta, double *b, int bstride)
 
/*  Compute sparse matrix time dense matrix multiply. Only works for VBR now.
*/

{
  FILE *dump_file ;
  int i, j, jj, irhs;
  int *indx,  *bindx, * rpntr,  *cpntr, *bpntrb, *bpntre;
  int nrow, ncol, n_blk_col;
  double *val, *xptr, *Aptr, *bptr, d_one = 1.0;

  val = A->val;
  indx = A->indx;
  bindx = A->bindx;
  rpntr = A->rpntr;
  cpntr = A->cpntr;
  bpntrb = A->bpntrb;
  bpntre = A->bpntre;

  /* Open dump file */

  dump_file = fopen(dump_file_name, "w");
  
  /*Find length of column partition vector */

  n_blk_col = 0;
  for (i=0; i<m; i++)
      for (j=bpntrb[i]; j<bpntre[i]; j++)
	{
	  n_blk_col = MAX(n_blk_col,bindx[j]);
	  n_blk_col = MAX(n_blk_col,bindx[j]+1);
	}

  /* Dump out parameters */

  fprintf(dump_file,"%d %d %d %d %d %d\n",m, n_blk_col, nrhs, k, xstride, bstride);
  fprintf(dump_file,"%22.16e %22.16e\n",alpha, beta);

  /* Dump out each array */

  for (i=0; i<m+1; i++)
    fprintf(dump_file,"%d\n",rpntr[i]);

  for (i=0; i<n_blk_col+1; i++)
    fprintf(dump_file,"%d\n",cpntr[i]);

  for (i=0; i<m+1; i++)
    fprintf(dump_file,"%d\n",bpntrb[i]);

  for (i=0; i<m+1; i++)
    fprintf(dump_file,"%d\n",bpntre[i]);

  for (i=0; i<m; i++)
      for (j=bpntrb[i]; j<bpntre[i]; j++)
	fprintf(dump_file,"%d\n",bindx[j]);

  for (i=0; i<m; i++)
      for (j=bpntrb[i]; j<bpntre[i]; j++)
	fprintf(dump_file,"%d\n",indx[j]);

  /* indx is one element longer, dump it too */

  fprintf(dump_file,"%d\n",indx[bpntre[m-1]]);


  for (i=0; i<m; i++)
      for (j=bpntrb[i]; j<bpntre[i]; j++)
	for (jj = indx[j]; jj<indx[j+1]; jj++)
	  fprintf(dump_file,"%22.16e\n",val[jj]);

  for (irhs=0; irhs<nrhs; irhs++)
    for (j=irhs*xstride; j<xstride+irhs*xstride; j++)
      fprintf(dump_file,"%22.16e\n",x[j]);

  for (irhs=0; irhs<nrhs; irhs++)
    for (j=irhs*bstride; j<bstride+irhs*bstride; j++)
      fprintf(dump_file,"%22.16e\n",b[j]);

  fclose(dump_file);

/* end cblas_duscr_dump */
}
