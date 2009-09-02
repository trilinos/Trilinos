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
#include "paz_aztec.h"
#include "prototypes.h"

void distrib_msr_matrix(int *proc_config, int *N_global,
	      int *n_nonzeros, int *N_update, int **update,
	      double **val, int **bindx,
	      double **x, double **b, double **bt, double **xexact)
#undef DEBUG 

{
  int i, n_entries, N_columns, n_global_nonzeros;
  int ii, j, row, have_xexact = 0 ;
  int kk = 0;
  int max_ii = 0, max_jj = 0;
  int ione = 1;
  double value;
  double *cnt;
  int *pntr, *bindx1, *pntr1;
  double *val1, *b1, *bt1, *x1, *xexact1;

  printf("Processor %d of %d entering distrib_matrix.\n",
	 proc_config[PAZ_node],proc_config[PAZ_N_procs]) ;

  /*************** Distribute global matrix to all processors ************/

  if(proc_config[PAZ_node] == 0)
    {
      if ((*xexact) != NULL) have_xexact = 1;
      printf("Broadcasting exact solution\n");
    }

  if(proc_config[PAZ_N_procs]  > 1) { 

      PAZ_broadcast((char *) N_global,        sizeof(int), proc_config, PAZ_PACK);
      PAZ_broadcast((char *) n_nonzeros, sizeof(int), proc_config, PAZ_PACK);
      PAZ_broadcast((char *) &have_xexact, sizeof(int), proc_config, PAZ_PACK);
      PAZ_broadcast(NULL, 0, proc_config, PAZ_SEND);

      if(proc_config[PAZ_node] != 0)
	{
	  (*bindx) = (int   *) calloc(*n_nonzeros+1,sizeof(int)) ;
	  (*val) = (double *) calloc(*n_nonzeros+1,sizeof(double)) ;
	}

      PAZ_broadcast((char *) (*bindx), sizeof(int)   *(*n_nonzeros+1), 
		   proc_config, PAZ_PACK);
      PAZ_broadcast(NULL, 0, proc_config, PAZ_SEND);
      PAZ_broadcast((char *) (*val),  sizeof(double)*(*n_nonzeros+1), 
		   proc_config, PAZ_PACK);
      PAZ_broadcast(NULL, 0, proc_config, PAZ_SEND);

      printf("Processor %d of %d done with matrix broadcast.\n",
	     proc_config[PAZ_node],proc_config[PAZ_N_procs]) ;
 
      /* Set rhs and initialize guess */
      if(proc_config[PAZ_node] != 0)
	{
	  (*b) = (double *) calloc(*N_global,sizeof(double)) ;
	  (*bt) = (double *) calloc(*N_global,sizeof(double)) ;
	  (*x) = (double *) calloc(*N_global,sizeof(double)) ;
	  if (have_xexact)
	  (*xexact) =   (double *) calloc(*N_global,sizeof(double)) ;
	}

      PAZ_broadcast((char *) (*x), sizeof(double)*(*N_global), proc_config, PAZ_PACK);
      PAZ_broadcast((char *) (*b), sizeof(double)*(*N_global), proc_config, PAZ_PACK);
      PAZ_broadcast((char *) (*bt), sizeof(double)*(*N_global), proc_config, PAZ_PACK);
      if (have_xexact)
	PAZ_broadcast((char *) 
		     (*xexact), sizeof(double)*(*N_global), proc_config, PAZ_PACK);
      PAZ_broadcast(NULL, 0, proc_config, PAZ_SEND);
      printf("Processor %d of %d done with rhs/guess broadcast.\n",
	     proc_config[PAZ_node],proc_config[PAZ_N_procs]) ;

    }

  /********************** Generate update map  *************************/

  PAZ_read_update(N_update, update, proc_config, (*N_global),
           1, PAZ_linear) ;
  
  printf("Processor %d of %d has %d rows of %d total rows.\n",
	 proc_config[PAZ_node],proc_config[PAZ_N_procs],*N_update,(*N_global)) ;

  /*************** Construct local matrix from global matrix ************/

  /* The local matrix is a copy of the rows assigned to this processor.  
     It is stored in MSR format and still has global indices (PAZ_transform
     will complete conversion to local indices.
  */

  if(proc_config[PAZ_N_procs]  > 1) { 
      n_global_nonzeros = *n_nonzeros;

      *n_nonzeros = *N_update;
      
      for (i=0; i<*N_update; i++)
	*n_nonzeros += (*bindx)[(*update)[i]+1] - (*bindx)[(*update)[i]];

      printf("Processor %d of %d has %d nonzeros of %d total nonzeros.\n",
	     proc_config[PAZ_node],proc_config[PAZ_N_procs],
	     *n_nonzeros,n_global_nonzeros) ;

#ifdef DEBUG
      { double sum1 = 0.0;
      for (i=0;i<(*N_global); i++) sum1 += (*b)[i];

      printf("Processor %d of %d has sum of b = %12.4g.\n",
	     proc_config[PAZ_node],proc_config[PAZ_N_procs],sum1) ;
      }
#endif /* DEBUG */

      /* Allocate memory for local matrix */

      bindx1 = (int   *) calloc(*n_nonzeros+1,sizeof(int)) ;
      val1 = (double *) calloc(*n_nonzeros+1,sizeof(double)) ;
      b1 =   (double *) calloc(*N_update,sizeof(double)) ;
      bt1 =   (double *) calloc(*N_update,sizeof(double)) ;
      x1 =   (double *) calloc(*N_update,sizeof(double)) ;
      if (have_xexact)
      xexact1 =   (double *) calloc(*N_update,sizeof(double)) ;
     
      bindx1[0] = *N_update+1;
      
      for (i=0; i<*N_update; i++)
	{
	  row = (*update)[i];
	  b1[i] = (*b)[row];
	  bt1[i] = (*bt)[row];
	  x1[i] = (*x)[row];
	  if (have_xexact) xexact1[i] = (*xexact)[row];
	  val1[i] = (*val)[row];
	  bindx1[i+1] = bindx1[i];

#ifdef DEBUG	  
	  printf("Proc %d of %d: Global row = %d: Local row = %d: 
                  b = %12.4g: x = %12.4g: bindx = %d: val = %12.4g \n",
		 proc_config[PAZ_node],proc_config[PAZ_N_procs], 
		 row, i, b1[i], x1[i], bindx1[i], val1[i]) ;
#endif

	  for (j = (*bindx)[row]; j < (*bindx)[row+1]; j++)
	    {
	      val1[  bindx1 [i+1] ] = (*val)[j];
	      bindx1[bindx1 [i+1] ] = (*bindx)[j];
	      bindx1[i+1] ++;
	    }
	}

      printf("Processor %d of %d done with extracting local operators.\n",
	     proc_config[PAZ_node],proc_config[PAZ_N_procs]) ;

      if (have_xexact)
	{
	  printf(
     "The residual using MSR format and exact solution on processor %d is %12.4g\n",
	      proc_config[PAZ_node],
	      smsrres (*N_update, (*N_global), val1, bindx1, xexact1, (*xexact), b1));
	}
  
      /* Release memory for global matrix, rhs and solution */
      
      free ((void *) (*val));
      free ((void *) (*bindx));
      free ((void *) (*b));
      free ((void *) (*bt));
      free ((void *) (*x));
      if (have_xexact) free((void *) *xexact);

      /* Return local matrix through same pointers. */
      
      *val = val1;
      *bindx = bindx1;
      *b = b1;
      *bt = bt1;
      *x = x1;
      if (have_xexact) *xexact = xexact1;

  }
  if (have_xexact && proc_config[PAZ_N_procs]  == 1)
    {
      printf(
	     "The residual using MSR format and exact solution on processor %d is %12.4g\n",
	     proc_config[PAZ_node],
	     smsrres (*N_update, (*N_global), (*val), (*bindx), 
		      (*xexact), (*xexact), (*b)));
    }
  
  
  printf("Processor %d of %d leaving distrib_matrix.\n",
	 proc_config[PAZ_node],proc_config[PAZ_N_procs]) ;
  
  /* end distrib_matrix */
}
