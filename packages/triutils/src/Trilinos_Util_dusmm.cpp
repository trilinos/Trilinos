// @HEADER
// ***********************************************************************
// 
//                 TriUtils: Trilinos Utilities Package
//                 Copyright (2001) Sandia Corporation
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
// @HEADER

#include "Trilinos_Util.h"
#include "Epetra_BLAS.h"

void Trilinos_Util_dusmm(int m, int nrhs, int k, double alpha, SPBLASMAT *A,
		 double *x, int xstride, double beta, double *b, int bstride)
 
/*  Compute sparse matrix time dense matrix multiply. Only works for VBR now.
*/

{
  int i, irhs, irhs1, irhs_begin, irhs_end, j, ix, xlen;
  int minblocksize, maxblocksize;
  int *indx,  *bindx, * rpntr,  *cpntr, *bpntrb, *bpntre, *ncolvec;
  int nrow, ncol, nrhs1, j1;
  double *val,  *buffer, *bufferptr, *bufferptr1;
  double *xptr, *xptr1, *Aptr, *bptr;
  double sum0, sum1, xtmp0, xtmp1;

  Epetra_BLAS blas;
  val = A->val;
  indx = A->indx;
  bindx = A->bindx;
  rpntr = A->rpntr;
  cpntr = A->cpntr;
  bpntrb = A->bpntrb;
  bpntre = A->bpntre;
  buffer = A->buffer;
  ncolvec = A->ncolvec;
  maxblocksize = A->maxblocksize;
  minblocksize = A->minblocksize;

  /* Special case for blocksize of 1 */
  if (maxblocksize ==1)
    {
      for (i=0; i<m; i++)
	{
	  bptr = b;
	  xptr = x;
	  for (irhs=0; irhs< nrhs; irhs++)
	    {
	      if (beta == 0.0)
		bptr[i] = 0.0;
	      else
		bptr[i] *= beta;
	      
	      if (alpha == 1.0)
		for (j=bpntrb[i]; j<bpntre[i]; j++)
		  bptr[i] += xptr[bindx[j]]*val[j];
	      else if (alpha == -1.0)
		for (j=bpntrb[i]; j<bpntre[i]; j++)
		  bptr[i] -= xptr[bindx[j]]*val[j];
	      else
		for (j=bpntrb[i]; j<bpntre[i]; j++)
		  bptr[i] += alpha*xptr[bindx[j]]*val[j];
	      bptr += bstride;
	      xptr += xstride;
	    }
	}
      return;
    }	    
      
  /* Special case for blocksize of 2 */
  if (minblocksize == 2 && maxblocksize == 2)
    {
      for (i=0; i<m; i++)
	{
	  bptr = b;
	  xptr = x;
	  for (irhs=0; irhs< nrhs; irhs++)
	    {
	      sum0 = 0.0;
	      sum1 = 0.0;
	      for (j=bpntrb[i]; j<bpntre[i]; j++)
		{
		  j1=indx[j];
		  xtmp0 = xptr[2*bindx[j]  ];
		  xtmp1 = xptr[2*bindx[j]+1];
		  sum0 += xtmp0*val[j1]   + xtmp1*val[j1+2];
		  sum1 += xtmp0*val[j1+1] + xtmp1*val[j1+3];
		}
	      if (beta == 0.0)
		{
		bptr[2*i]   = alpha*sum0;
		bptr[2*i+1] = alpha*sum1;
		}
	      else
		{
		bptr[2*i]   = beta * bptr[2*i]    + alpha*sum0;
		bptr[2*i+1] = beta * bptr[2*i+1]  + alpha*sum1;
		}
	      bptr += bstride;
	      xptr += xstride;
	    }
	}
      return;
    }	    
      

  /* Compute Matrix multiply block row by block row using DGEMM */

  for (i=0; i<m; i++)
    {
      nrow = rpntr[i+1] - rpntr[i];
      ncol = ncolvec[i];
      Aptr = val+indx[bpntrb[i]];
      /* Copy RHS to buffer, Stripmine by MAXNRHS (defined in spblas.h) */
      for (irhs=0; irhs<nrhs; irhs+=MAXNRHS)
	{
	  irhs_begin = irhs;
	  irhs_end = Trilinos_Util_min(irhs+MAXNRHS, nrhs);
	  nrhs1 = irhs_end - irhs_begin;
	  bufferptr = buffer;
	  xptr = x+irhs_begin*xstride;
	  bptr = b+rpntr[i]+irhs_begin*bstride;
	  for (irhs1=irhs_begin; irhs1<irhs_end; irhs1++)
	    {
	      bufferptr1 = bufferptr;
	      for (j=bpntrb[i]; j<bpntre[i]; j++)
		{
		  xptr1 = xptr+cpntr[bindx[j]];
		  xlen = cpntr[bindx[j]+1] - cpntr[bindx[j]];
		  for (ix=0; ix<xlen; ix++)
		  {
		    *bufferptr1 = *xptr1;
		    bufferptr1++; xptr1++;
		  }
		}
	      bufferptr += ncol;
	      xptr += xstride;
	    }
	  blas.GEMM('N', 'N', nrow, nrhs1, ncol, alpha, Aptr, nrow, buffer, ncol,
	     beta, bptr, bstride);
	}
    }
	  




/* end Trilinos_Util_dusmm_vbr */
}
