// @HEADER
// ***********************************************************************
// 
//                 TriUtils: Trilinos Utilities Package
//                 Copyright (2011) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Trilinos_Util.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"

// CJ TODO FIXME: Trilinos_Util_distrib_vbr_matrix available only if 32 bit GIDs available.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES

void Trilinos_Util_distrib_vbr_matrix(const Epetra_Comm & Comm,
	      int *N_global, int *N_blk_global,
	      int *n_nonzeros, int *n_blk_nonzeros, 
	      int *N_update, int **update,
	      double **val, int **indx, int **rpntr, int **cpntr,
	      int **bpntr, int **bindx,
	      double **x, double **b, double **xexact)
#undef DEBUG 

{
  int i, n_global_nonzeros, n_global_blk_nonzeros;
  int N_local;
  int j, row, have_xexact = 0 ;
  int *rpntr1, *bindx1, *bpntr1, *indx1;
  double *val1, *b1, *x1, *xexact1=0;
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  printf("Processor %d of %d entering distrib_matrix.\n",
	 MyPID,NumProc) ;

  /*************** Distribute global matrix to all processors ************/

  if(MyPID == 0)
    {
      if ((*xexact) != NULL) have_xexact = 1;
      printf("%s", "Broadcasting exact solution\n");
    }

  if(NumProc  > 1)
    { 

      Comm.Broadcast( N_global,      1, 0);
      Comm.Broadcast( N_blk_global,  1, 0);
      Comm.Broadcast( n_nonzeros,     1, 0);
      Comm.Broadcast( n_blk_nonzeros, 1, 0);
      Comm.Broadcast( &have_xexact,   1, 0);

      printf("Processor %d of %d done with global parameter  broadcast.\n",
	     MyPID,NumProc) ;

      if(MyPID != 0)
	{
      *bpntr = (int   *) calloc(*N_blk_global+1,sizeof(int)) ;
      *rpntr = (int   *) calloc(*N_blk_global+1,sizeof(int)) ;
      *bindx = (int   *) calloc(*n_blk_nonzeros+1,sizeof(int)) ;
      *indx  = (int   *) calloc(*n_blk_nonzeros+1,sizeof(int)) ;
      *val = (double *) calloc(*n_nonzeros+1,sizeof(double)) ;
      printf("Processor %d of %d done with global calloc.\n",
	     MyPID,NumProc) ;
}

      Comm.Broadcast( (*bpntr), (*N_blk_global+1), 0);
      Comm.Broadcast( (*rpntr), (*N_blk_global+1), 0);
      Comm.Broadcast( (*bindx), (*n_blk_nonzeros+1), 0);
      Comm.Broadcast( (*indx),  (*n_blk_nonzeros+1), 0);
      Comm.Broadcast( (*val),  (*n_nonzeros+1), 0);

      printf("Processor %d of %d done with matrix broadcast.\n",
	     MyPID,NumProc) ;
 
      /* Set rhs and initialize guess */
      if(MyPID != 0)
	{
	  (*b) = (double *) calloc(*N_global,sizeof(double)) ;
	  (*x) = (double *) calloc(*N_global,sizeof(double)) ;
	  if (have_xexact)
	  (*xexact) =   (double *) calloc(*N_global,sizeof(double)) ;
	}

      Comm.Broadcast( (*x), (*N_global), 0);
      Comm.Broadcast( (*b), (*N_global), 0);
      if (have_xexact)
	Comm.Broadcast((*xexact), (*N_global), 0);
      printf("Processor %d of %d done with rhs/guess broadcast.\n",
	     MyPID,NumProc) ;

    }

  /********************** Generate update map  *************************/

  //read_update(N_update, update, proc_config, *N_blk_global, 1, linear) ;

  Epetra_Map map(*N_blk_global, 0, Comm);
  *N_update = map.NumMyElements();
  (*update) = (int *) calloc(*N_update,sizeof(int)) ;
  map.MyGlobalElements(*update);

  printf("Processor %d of %d has %d rows of %d total block rows.\n",
	 MyPID,NumProc,*N_update,*N_blk_global) ;

  /*************** Construct local matrix from global matrix ************/

  /* The local matrix is a copy of the rows assigned to this processor.  
     It is stored in MSR format and still has global indices 
  */

  if(NumProc  > 1)
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
	     MyPID,NumProc,
	     *n_nonzeros,n_global_nonzeros) ;

   printf("Processor %d of %d has %d block nonzeros of %d total block nonzeros.\n",
	     MyPID,NumProc,
	     *n_blk_nonzeros,n_global_blk_nonzeros) ;

   printf("Processor %d of %d has %d equations of %d total equations.\n",
	     MyPID,NumProc,
	     N_local,*N_global) ;

#ifdef DEBUG
      { double sum1 = 0.0;
      for (i=0;i<*N_global; i++) sum1 += (*b)[i];

      printf("Processor %d of %d has sum of b = %12.4g.\n",
	     MyPID,NumProc,sum1) ;
      }
#endif /* DEBUG */

      /* Allocate memory for local matrix */

      bpntr1 = (int   *) calloc(*N_update+1,sizeof(int)) ;
      rpntr1 = (int   *) calloc(*N_update+1,sizeof(int)) ;
      bindx1 = (int   *) calloc(*n_blk_nonzeros+1,sizeof(int)) ;
      indx1  = (int   *) calloc(*n_blk_nonzeros+1,sizeof(int)) ;
      val1 = (double *) calloc(*n_nonzeros+1,sizeof(double)) ;
      b1 =   (double *) calloc(N_local,sizeof(double)) ;
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
	    printf("Proc %d of %d: Global row = %d: Local row = %d: b = %12.4g: x = %12.4g: bindx = %d: val = %12.4g \n",
		    MyPID,NumProc, 
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
	     MyPID,NumProc) ;

      if (have_xexact)
	{
	  printf(
     "The residual using VBR format and exact solution on processor %d is %12.4g\n",
	      MyPID,
	      Trilinos_Util_svbrres (N_local, *N_global, *N_update, val1, indx1, bindx1, 
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
      free ((void *) (*x));
      if (have_xexact) free((void *) *xexact);

      /* Return local matrix through same pointers. */
      
      *val = val1;
      *indx = indx1;
      *bindx = bindx1;
      *bpntr = bpntr1;
      *rpntr = rpntr1;
      *b = b1;
      *x = x1;
      if (have_xexact) *xexact = xexact1;

    }
      if (have_xexact && NumProc  == 1)
	{
	  printf(
     "The residual using VBR format and exact solution on processor %d is %12.4g\n",
	      MyPID,
	      Trilinos_Util_svbrres (*N_global, *N_global, *N_update, (*val), (*indx), (*bindx), 
		       (*rpntr), (*rpntr), (*bpntr), (*bpntr)+1,
		       (*xexact), (*b)));
	}

  
  printf("Processor %d of %d leaving distrib_matrix.\n",
	 MyPID,NumProc) ;
  
  /* end distrib_matrix */
}

#endif
