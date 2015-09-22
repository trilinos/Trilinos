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
