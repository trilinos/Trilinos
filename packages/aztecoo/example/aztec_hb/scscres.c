#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "az_aztec.h"
double scscres (int isym, int m, int n, 
	      double *val, int *indx, int *pntr,
	      double *x, double *b)
{
    int i, j, ibgn, iend, ione = 1;
    double norm_tmp = 0.0, norm_b = 0.0;
    double scaled_res_norm, res_norm, *tmp, max_norm = 0.0;


/*     Computes the residual

                      res = || b - A*x ||

       where x and b are vectors and A is a sparse matrix stored
       in MSR format. */

/*     -------------------------- 
       First executable statement 
       -------------------------- */

    /* Create tmp workspace */
    tmp = (double *) calloc(m,sizeof(double));

/* .....initialize soln */

    for (i = 0; i < m; i++)
	tmp[i] = b[i];

/* .....do a series of SPAXPYs (sparse saxpys) */

    for (j = 0; j < n ; j++) 
      {
	ibgn = pntr[j];
	iend = pntr[j + 1];
	
	for (i = ibgn; i < iend; i++)
	  {
	    tmp[indx[i]] -= val[i] * x[j];
 	    if (indx[i] != j && isym) tmp[j] -= val[i]*x[indx[i]];
	  }
     }
    for (i = 0; i < m; i++)
      {
	max_norm = AZ_MAX(fabs(tmp[i]),max_norm);
	norm_tmp += tmp[i]*tmp[i];
	norm_b += b[i]*b[i];
      }
   
    res_norm = sqrt(norm_tmp);
    printf("\n\nMax norm of residual        = %12.4g\n",max_norm);
    printf(    "Two norm of residual        = %12.4g\n",res_norm);
    if (norm_b > 1.0E-7) 
      {
	   scaled_res_norm = res_norm/sqrt(norm_b);
	   printf(    "Scaled two norm of residual = %12.4g\n",scaled_res_norm);
      }

    free((void *) tmp);

    return(scaled_res_norm);

} /* scscres */

