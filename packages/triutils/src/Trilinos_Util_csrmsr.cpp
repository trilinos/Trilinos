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

int Trilinos_Util_csrmsr( int n, double *a, int *ja, int *ia, double *ao,
             int *jao, double *wk, int *iwk)
{
    int iptr, i, j, k, ii, icount;

/*
 -----------------------------------------------------------------------
 Compressed Sparse Row   to      Modified - Sparse Row 
                                 Sparse row with separate main diagonal 
 -----------------------------------------------------------------------
 converts a general sparse matrix a, ja, ia into 
 a compressed matrix using a separated diagonal (referred to as 
 the bell-labs format as it is used by bell labs semi conductor 
 group. We refer to it here as the modified sparse row format. 
 Note: this has been coded in such a way that one can overwrite 
 the output matrix onto the input matrix if desired by a call of 
 the form 

     call csrmsr (n, a, ja, ia, a, ja, wk,iwk) 

 In case ao, jao, are different from a, ja, then one can 
 use ao, jao as the work arrays in the calling sequence: 

     call csrmsr (n, a, ja, ia, ao, jao, ao,jao) 

 -----------------------------------------------------------------------
 

 on entry : 
 --------- 
 a, ja, ia = matrix in csr format. note that the 
 	     algorithm is in place: ao, jao can be the same 
            as a, ja, in which case it will be overwritten on it 
            upon return. 

 on return : 
 ----------- 

 ao, jao  = sparse matrix in modified sparse row storage format: 
 	   +  ao(1:n) contains the diagonal of the matrix. 
 	   +  ao(n+2:nnz) contains the nondiagonal elements of the 
             matrix, stored rowwise. 
 	   +  jao(n+2:nnz) : their column indices 
 	   +  jao(1:n+1) contains the pointer array for the nondiagonal 
             elements in ao(n+1:nnz) and jao(n+2:nnz). 
             i.e., for i .le. n+1 jao(i) points to beginning of row i 
 	      in arrays ao, jao. 
 	       here nnz = number of nonzero elements+1 
 work arrays: 
 ------------ 
 wk	= real work array of length n 
 iwk   = integer work array of length n+1 

 notes: 
 ------- 
        Algorithm is in place.  i.e. both: 

          call csrmsr (n, a, ja, ia, ao, jao, ao,jao) 
          (in which  ao, jao, are different from a, ja) 
           and 
          call csrmsr (n, a, ja, ia, a, ja, wk,iwk) 
          (in which  wk, jwk, are different from a, ja) 
        are OK. 
 -------- 
 -----------------------------------------------------------------------
*/ 
    icount = 0;

 // store away diagonal elements and count nonzero diagonal elements. 

    for (i = 0; i < n; i++) {
	wk[i] = 0.;
	iwk[i+1] = ia[i+1] - ia[i];
	for (k = ia[i]; k <ia[i+1]; k++) {
	    if (ja[k]==i) {
		wk[i] = a[k];
		++icount;
		--iwk[i+1];
	    }
	}
    }

 // compute total length 

    iptr = n + ia[n] - icount;

     // copy backwards (to avoid collisions) 

    for (ii = n-1; ii >=0; ii--) {
	for (k = ia[ii+1]-1; k >= ia[ii]; k--) {
	    j = ja[k];
	    if (j != ii) {
		ao[iptr] = a[k];
		jao[iptr] = j;
		--iptr;
	    }
	}
    }

  //  compute pointer values and copy wk(*) 

    jao[0] = n + 1;
    for (i = 0; i <n; i++) {
	ao[i] = wk[i];
	jao[i+1] = jao[i] + iwk[i+1];
    }
    return 0;
}

