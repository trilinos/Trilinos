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

void distrib_vbr_matrix(int *proc_config,
	      int *N_global, int *N_blk_global,
	      int *n_nonzeros, int *n_blk_nonzeros, 
	      int *N_update, int **update,
	      double **val, int **indx, int **rpntr, int **cpntr,
	      int **bpntr, int **bindx,
	      double **x, double **b, double **bt, double **xexact)
#undef DEBUG 

{
  int i, n_entries, N_columns, n_global_nonzeros, n_global_blk_nonzeros;
  int N_local;
  int ii, j, row, have_xexact = 0 ;
  int kk = 0;
  int max_ii = 0, max_jj = 0;
  int ione = 1;
  double value;
  double *cnt;
  int *rpntr1, *bindx1, *bpntr1, *indx1;
  double *val1, *b1, *bt1, *x1, *xexact1;

  printf("Processor %d of %d entering distrib_matrix.\n",
	 proc_config[PAZ_node],proc_config[PAZ_N_procs]) ;

  /*************** Distribute global matrix to all processors ************/

  if(proc_config[PAZ_node] == 0)
    {
      if ((*xexact) != NULL) have_xexact = 1;
      printf("Broadcasting exact solution\n");
    }

  if(proc_config[PAZ_N_procs]  > 1)
    { 

      PAZ_broadcast((char *) N_global,      sizeof(int), proc_config, PAZ_PACK);
      PAZ_broadcast((char *) N_blk_global,  sizeof(int), proc_config, PAZ_PACK);
      PAZ_broadcast((char *) n_nonzeros,     sizeof(int), proc_config, PAZ_PACK);
      PAZ_broadcast((char *) n_blk_nonzeros, sizeof(int), proc_config, PAZ_PACK);
      PAZ_broadcast((char *) &have_xexact,   sizeof(int), proc_config, PAZ_PACK);
      PAZ_broadcast(NULL, 0, proc_config, PAZ_SEND);

      printf("Processor %d of %d done with global parameter  broadcast.\n",
	     proc_config[PAZ_node],proc_config[PAZ_N_procs]) ;

      if(proc_config[PAZ_node] != 0)
	{
      *bpntr = (int   *) calloc(*N_blk_global+1,sizeof(int)) ;
      *rpntr = (int   *) calloc(*N_blk_global+1,sizeof(int)) ;
      *bindx = (int   *) calloc(*n_blk_nonzeros+1,sizeof(int)) ;
      *indx  = (int   *) calloc(*n_blk_nonzeros+1,sizeof(int)) ;
      *val = (double *) calloc(*n_nonzeros+1,sizeof(double)) ;
      printf("Processor %d of %d done with global calloc.\n",
	     proc_config[PAZ_node],proc_config[PAZ_N_procs]) ;
}

      PAZ_broadcast((char *) (*bpntr), sizeof(int)   *(*N_blk_global+1), 
		   proc_config, PAZ_PACK);
      PAZ_broadcast((char *) (*rpntr), sizeof(int)   *(*N_blk_global+1), 
		   proc_config, PAZ_PACK);
      PAZ_broadcast((char *) (*bindx), sizeof(int)   *(*n_blk_nonzeros+1), 
		   proc_config, PAZ_PACK);
      PAZ_broadcast((char *) (*indx),  sizeof(int)   *(*n_blk_nonzeros+1), 
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

  PAZ_read_update(N_update, update, proc_config, *N_blk_global,
           1, PAZ_linear) ;

  printf("Processor %d of %d has %d rows of %d total block rows.\n",
	 proc_config[PAZ_node],proc_config[PAZ_N_procs],*N_update,*N_blk_global) ;

  /*************** Construct local matrix from global matrix ************/

  /* The local matrix is a copy of the rows assigned to this processor.  
     It is stored in MSR format and still has global indices (PAZ_transform
     will complete conversion to local indices.
  */

  if(proc_config[PAZ_N_procs]  > 1)
    { 
      n_global_nonzeros = *n_nonzeros;
      n_global_blk_nonzeros = *n_blk_nonzeros;

      *n_nonzeros = 0;
      *n_blk_nonzeros = 0;
      N_local = 0;
      
      for (i=0; i<*N_update; i++)
	{
	  row = (*update)[i];
	  *n_nonzeros     += (*indx)[(*bpntr)[row+1]] - (*indx)[(*bpntr)[row]];
	  *n_blk_nonzeros += (*bpntr)[row+1] - (*bpntr)[row];
	  N_local         += (*rpntr)[row+1] - (*rpntr)[row];
	  
	}

      printf("Processor %d of %d has %d nonzeros of %d total nonzeros.\n",
	     proc_config[PAZ_node],proc_config[PAZ_N_procs],
	     *n_nonzeros,n_global_nonzeros) ;

   printf("Processor %d of %d has %d block nonzeros of %d total block nonzeros.\n",
	     proc_config[PAZ_node],proc_config[PAZ_N_procs],
	     *n_blk_nonzeros,n_global_blk_nonzeros) ;

   printf("Processor %d of %d has %d equations of %d total equations.\n",
	     proc_config[PAZ_node],proc_config[PAZ_N_procs],
	     N_local,*N_global) ;

#ifdef DEBUG
      { double sum1 = 0.0;
      for (i=0;i<*N_global; i++) sum1 += (*b)[i];

      printf("Processor %d of %d has sum of b = %12.4g.\n",
	     proc_config[PAZ_node],proc_config[PAZ_N_procs],sum1) ;
      }
#endif /* DEBUG */

      /* Allocate memory for local matrix */

      bpntr1 = (int   *) calloc(*N_update+1,sizeof(int)) ;
      rpntr1 = (int   *) calloc(*N_update+1,sizeof(int)) ;
      bindx1 = (int   *) calloc(*n_blk_nonzeros+1,sizeof(int)) ;
      indx1  = (int   *) calloc(*n_blk_nonzeros+1,sizeof(int)) ;
      val1 = (double *) calloc(*n_nonzeros+1,sizeof(double)) ;
      b1 =   (double *) calloc(N_local,sizeof(double)) ;
      bt1 =   (double *) calloc(N_local,sizeof(double)) ;
      x1 =   (double *) calloc(N_local,sizeof(double)) ;
      if (have_xexact)
      xexact1 =   (double *) calloc(N_local,sizeof(double)) ;

      {     
	int cur_blk_size, indx_offset, len_val, row_offset, row_offset1;
	double *val_ptr, *val1_ptr;

	bpntr1[0] = 0;
	indx1[0] = 0;
	rpntr1[0] = 0;
	for (i=0; i<*N_update; i++)
	  {
	    row = (*update)[i];
	    cur_blk_size = (*rpntr)[row+1] - (*rpntr)[row];
	    rpntr1[i+1] = rpntr1[i] + cur_blk_size;
	    row_offset = (*rpntr)[row];
	    row_offset1 = rpntr1[i];
	    for (j = 0; j<cur_blk_size; j++)
	      {
		b1[row_offset1+j] = (*b)[row_offset+j];
		x1[row_offset1+j] = (*x)[row_offset+j];
		if (have_xexact) xexact1[row_offset1+j] = (*xexact)[row_offset+j];
	      }
	    bpntr1[i+1] = bpntr1[i];
	    
#ifdef DEBUG	  
	    printf("Proc %d of %d: Global row = %d: Local row = %d: 
                    b = %12.4g: x = %12.4g: bindx = %d: val = %12.4g \n",
		    proc_config[PAZ_node],proc_config[PAZ_N_procs], 
		    row, i, b1[i], x1[i], bindx1[i], val1[i]) ;
#endif
	    indx_offset = (*indx)[(*bpntr)[row]] - indx1[bpntr1[i]];
	    for (j = (*bpntr)[row]; j < (*bpntr)[row+1]; j++)
	      {
		indx1[bpntr1 [i+1] + 1] = (*indx)[j+1] - indx_offset;
		bindx1[bpntr1 [i+1] ] = (*bindx)[j];
		bpntr1[i+1] ++;
	      }
	    len_val = indx1[bpntr1[i+1]] - indx1[bpntr1[i]];
	    val_ptr = (*val)+(*indx)[(*bpntr)[row]];
	    val1_ptr = val1+indx1[bpntr1[i]];
	    for (j = 0; j<len_val; j++)
	      { 
		*val1_ptr = *val_ptr;
		val_ptr++; val1_ptr++;
	      }
	  }
      }
      printf("Processor %d of %d done with extracting local operators.\n",
	     proc_config[PAZ_node],proc_config[PAZ_N_procs]) ;

      if (have_xexact)
	{
	  printf(
     "The residual using VBR format and exact solution on processor %d is %12.4g\n",
	      proc_config[PAZ_node],
	      svbrres (N_local, *N_global, *N_update, val1, indx1, bindx1, 
		       rpntr1, (*rpntr), bpntr1, bpntr1+1,
		       (*xexact), b1));
	}
  
      /* Release memory for global matrix, rhs and solution */
      
      free ((void *) (*val));
      free ((void *) (*indx));
      free ((void *) (*bindx));
      free ((void *) (*bpntr));
      free ((void *) (*rpntr));
      free ((void *) (*b));
      free ((void *) (*bt));
      free ((void *) (*x));
      if (have_xexact) free((void *) *xexact);

      /* Return local matrix through same pointers. */
      
      *val = val1;
      *indx = indx1;
      *bindx = bindx1;
      *bpntr = bpntr1;
      *rpntr = rpntr1;
      *b = b1;
      *bt = bt1;
      *x = x1;
      if (have_xexact) *xexact = xexact1;

    }
      if (have_xexact && proc_config[PAZ_N_procs]  == 1)
	{
	  printf(
     "The residual using VBR format and exact solution on processor %d is %12.4g\n",
	      proc_config[PAZ_node],
	      svbrres (*N_global, *N_global, *N_update, (*val), (*indx), (*bindx), 
		       (*rpntr), (*rpntr), (*bpntr), (*bpntr)+1,
		       (*xexact), (*b)));
	}

  
  printf("Processor %d of %d leaving distrib_matrix.\n",
	 proc_config[PAZ_node],proc_config[PAZ_N_procs]) ;
  
  /* end distrib_matrix */
}
