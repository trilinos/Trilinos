/*@HEADER
// ***********************************************************************
// 
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "spblas.h"
#define max(x,y) (( x > y ) ? x : y)     /* max function  */

double svbrres (int m, int n, int m_blk,
		double *val, int *indx, int *bindx, int *rpntr,
		int *cpntr, int *bpntrb, int *bpntre,
		double *x, double *b)
{
    int i, j, jbgn, jend, ione = 1;
    double sum, norm_tmp = 0.0, norm_b = 0.0;
    double scaled_res_norm, res_norm, *tmp, max_norm = 0.0;
    SPBLASMAT  A;



/*     Computes the residual

                      res = || b - A*x ||

       where x and b are vectors and A is a sparse matrix stored
       in MSR format. */

/*     -------------------------- 
       First executable statement 
       -------------------------- */
    /* Create sparse matrix handle */
    cblas_duscr_vbr(m_blk, val, indx, bindx, rpntr, cpntr, bpntrb, bpntre, &A);
    

    /* Create tmp workspace, set to b */

    tmp = (double *) calloc(m,sizeof(double));

    for (i = 0; i < m; i++) tmp[i] = b[i];

    /* Call DUSMM to compute residual (in tmp) */

    cblas_dusmm(m_blk, 1, n, -1.0, &A, x, m, 1.0, tmp, m);

    for (i = 0; i <m ; i++) 
      {
	max_norm = max(fabs(tmp[i]),max_norm);
	norm_tmp += tmp[i]*tmp[i];
	norm_b += b[i]*b[i];
      }
   
    res_norm = sqrt(norm_tmp);
    scaled_res_norm = res_norm/sqrt(norm_b);
    printf("\n\nMax norm of residual        = %12.4g\n",max_norm);
    printf(    "Two norm of residual        = %12.4g\n",res_norm);
    if (norm_b > 1.0E-7) 
      {
	scaled_res_norm = res_norm/sqrt(norm_b);
	printf(    "Scaled two norm of residual = %12.4g\n",scaled_res_norm);
      }
    /* Compute residual statistics */
    /*      if (res_norm > 0.2 )
    cblas_dusmm_dump("/u1/mheroux/dump_file",
		     m_blk, 1, n, -1.0, &A, x, n, 1.0, b, m);
   for (i=0; i<m_blk; i++)
      {
	printf("***** Row %d *******\n",i);
	printf("bpntrb[%d] = %d\n",i,bpntrb[i]);
	printf("bpntre[%d] = %d\n",i,bpntre[i]);
	printf("rpntr[%d] = %d\n",i,rpntr[i]);
	for (j=bpntrb[i]; j<bpntre[i]; j++)
	  {
	    printf("bindx[%d] = %d\n",j,bindx[j]);
	    printf("indx[%d] = %d\n",j,indx[j]);
	  }
	
	  
      }
	printf("rpntr[%d] = %d\n",m_blk,rpntr[m_blk]);
	j = bpntre[m_blk-1];
	printf("bindx[%d] = %d\n",j,bindx[j]);
	printf("indx[%d] = %d\n",j,indx[j]);
    printf("val[indx[bpntrb[m_blk-1]]  ] = %12.4g\n",val[indx[bpntrb[m_blk-1]]  ]);
    printf("val[indx[bpntrb[m_blk-1]]+1] = %12.4g\n",val[indx[bpntrb[m_blk-1]]+1]);
    printf("val[indx[bpntrb[m_blk-1]]+2] = %12.4g\n",val[indx[bpntrb[m_blk-1]]+2]);

    for (i = 0; i <m ; i++) 
      {
	printf("tmp[%d] = %12.4g\n",i,tmp[i]);
	printf("  x[%d] = %12.4g\n",i,  x[i]);
	printf("  b[%d] = %12.4g\n",i,  b[i]);
      }
    */
    free((void *) tmp);

    return(res_norm);

} /* svbrres */

