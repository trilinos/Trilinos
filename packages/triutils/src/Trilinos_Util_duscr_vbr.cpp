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

void Trilinos_Util_duscr_vbr(int n, double *val, int *indx, 
                             int *bindx, int *rpntr, 
		                   int *cpntr, int *bpntrb, int *bpntre, SPBLASMAT *A)
 
/*  Create and handle for a VBR matrix.  Build any auxiliary data structures
    that might be helpful for performance.
*/

{
  int i, j, buffersize, maxbuffersize, *ncolvec, minblocksize, maxblocksize;
  int blocksize;
  double *buffer, nops_per_rhs;

  A->val = val;
  A->indx = indx;
  A->bindx = bindx;
  A->rpntr = rpntr;
  A->cpntr = cpntr;
  A->bpntrb = bpntrb;
  A->bpntre = bpntre;
  A->n = n;

  /* Create n length int vector to contain compressed column lengths */
  
  ncolvec = (int *) calloc(n, sizeof(int));
  
  /* Compute size of and build work buffer to contain RHS blocks */

  maxbuffersize = 0;
  nops_per_rhs = 0.0;
  maxblocksize = 0;
  minblocksize = n;
  for (i=0; i<n; i++) 
    {
    buffersize = 0;
    for (j = bpntrb[i];j<bpntre[i]; j++)
      {
	blocksize = cpntr[bindx[j]+1] - cpntr[bindx[j]];
	buffersize += blocksize;
	minblocksize = Trilinos_Util_min(minblocksize,blocksize);
	maxblocksize = Trilinos_Util_max(maxblocksize,blocksize);
      }
    ncolvec[i] = buffersize;
    maxbuffersize = Trilinos_Util_max(maxbuffersize, buffersize);
    minblocksize = Trilinos_Util_min(minblocksize,rpntr[i+1]-rpntr[i]);
    maxblocksize = Trilinos_Util_max(maxblocksize,rpntr[i+1]-rpntr[i]);
    nops_per_rhs += (double) 2*(rpntr[i+1]-rpntr[i])*buffersize;
  }
  /* Get buffer space */
  buffer = (double *) calloc(maxbuffersize*MAXNRHS,sizeof(double));


  A->buffersize = maxbuffersize*MAXNRHS;
  A->bufferstride = maxbuffersize;
  A->buffer = buffer;
  A->ncolvec = ncolvec;
  A->nops_per_rhs = nops_per_rhs;
  A->minblocksize = minblocksize;
  A->maxblocksize = maxblocksize;
 
/* end Trilinos_Util_duscr_vbr */
}
