#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "az_aztec.h"
double smsrres (int m, int n, 
	      double *val, int *indx, 
	      double *xlocal, double *x, double *b)
{
    int i, j, jbgn, jend, ione = 1;
    double sum, norm_tmp = 0.0, norm_b = 0.0;
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
	tmp[i] = b[i] - val[i] * xlocal[i];

/* .....do a series of SPDOTs (sparse dot products) */

    for (i = 0; i <m ; i++) 
      {
	jbgn = indx[i];
	jend = indx[i + 1];
	sum = 0.0;
	
	for (j = jbgn; j < jend; j++)
	  sum += val[j] * x[indx[j]];

	tmp[i] -= sum;
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

} /* smsrres */

