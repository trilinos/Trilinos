void  scscmv (int isym, int indexbase, int m, int n, 
	      double *val, int *indx, int *pntr, 
	      double *x, double *y)
{
    int i, j, jbgn, jend;


/*     Performs the matrix-vector operation

                      y = A*x 

       where x and y are vectors and A is a sparse matrix stored
       in (Harwell-Boeing) compress column format. */

/*     -------------------------- 
       First executable statement 
       -------------------------- */

/* .....initialize soln */

    for (i = 0; i < m; i++)
	y[i] = 0.0;

/* .....do a series of SPAXPYs (sparse saxpys) */

    for (j = 0; j <n ; j++) 
      {
	jbgn = pntr[j]-indexbase;
	jend = pntr[j + 1]-indexbase;

	for (i = jbgn; i < jend; i++)
	  {
	    y[indx[i]-indexbase] += val[i] * x[j];
	    if (indx[i]-indexbase != j && isym) y[j] += val[i]*x[indx[i]-indexbase];
	  }
      }

    return;
} /* scscmv */

