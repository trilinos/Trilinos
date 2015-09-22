/*
//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
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
//@HEADER
*/

#include <stdio.h>
#include <stdlib.h>
#include "az_aztec.h"

/**************************************************************/
/**************************************************************/
/**************************************************************/
void AZ_lower_icc(int bindx[],double val[],int N, double rhs[])
{
/*
 * Back solve the lower triangular system given by bindx
 * and val.
 *
 */
 
   int i,j,col;

   for (i = 0 ; i < N ; i++ ) {
      for (j = bindx[i]; j < bindx[i+1]; j++ ) {
         col = bindx[j];
         rhs[col] -= (val[j]*rhs[i]);
      }
   }
   for (i = 0 ; i < N ; i++ ) rhs[i] = rhs[i]*val[i];
}

/**************************************************************/
/**************************************************************/
/**************************************************************/
void AZ_upper_icc(int bindx[],double val[],int N, double rhs[])
{
/*
 * Forward solve the upper triangular system given by bindx
 * and val.
 *
 */
   int i,j,col;

   for (i = N-1 ; i >= 0 ; i-- ) {
      for (j = bindx[i]; j < bindx[i+1]; j++ ) {
         col = bindx[j];
         rhs[i] -= (val[j]*rhs[col]);
      }
   }
}

/**************************************************************/
/**************************************************************/
/**************************************************************/

void AZ_fact_chol(int bindx[], double val[], int N, double rthresh,
                  double athresh)
{
/*
 * Compute the Cholesky factorization of the matrix
 * stored in bindx and val.
 *
 * Note: On input, it is required that the redundant
 * entries in the lower triangular part of the matrix
 * be given.
 *
 */
    int i,j,k,kk,tk;
    double *sum, temp;
    int    *mark, *diag;
    int    row, col, first, last, next_nz;

    diag  = (int    *) AZ_allocate(N*sizeof(int));
    mark  = (int    *) AZ_allocate(N*sizeof(int));
    sum   = (double *) AZ_allocate(N*sizeof(double));

    if (sum == NULL) {
       printf("Not enough memory to perform ICC factorization\n");
       exit(1);
    }

    for (i = 0 ; i < N; i++) sum[i]   = 0.0;
    for (i = 0 ; i < N; i++) mark[i]  = 0;

    /* find the diagonal entry in each row */

    first = bindx[0];
    for (i = 0 ; i < N ; i++) {
       last = bindx[i+1];
       for (j = first; j < last ; j++) 
          if (bindx[j] > i) break;
       diag[i] = j;
       first = last;
     }

    /* Shift diagonal */

      if (rthresh ==0.0) rthresh = 1.0; /* default is zero, means no rthresh */
      if (rthresh != 1.0 || athresh != 0.0)
         for (i = 0 ; i < N ; i++) {
             if (val[i] >=0)
               val[i] = val[i] * rthresh + athresh;
             else
               val[i] = val[i] * rthresh - athresh;
         }

    /* do the factorization */

    for (i = 0 ; i < N ; i++ ) {
       val[i] -= sum[i];
       for (kk = diag[i]; kk < bindx[i+1]; kk++) 
          mark[bindx[kk]] = kk+1;
       for (kk = bindx[i]; kk < diag[i]; kk++) {
          row = bindx[kk];
          for (k = diag[row]; k < bindx[row+1]; k++ )
             if (bindx[k] == i) break;
          if (bindx[k] != i) {
             printf("The matrix is not symmetric. Can not use ICC\n");
             exit(1);
          }
          temp = val[k];
          for (tk = k+1; tk < bindx[row+1]; tk++) {
             col = bindx[tk];
             if (mark[col]) val[mark[col]-1] -= temp*val[tk]*val[row];
          }
       }
       for (kk = diag[i]; kk < bindx[i+1]; kk++) {
          col = bindx[kk];
          mark[col] = 0;
          val[kk] = val[kk]/val[i];
          sum[col] += val[kk]*val[kk]*val[i];
       }
   }

   /* compress out the lower triangular part */

   next_nz = N+1;
   for (i = 0 ; i < N ; i++ ) {
      for (j = diag[i] ; j < bindx[i+1]; j++) {
         bindx[next_nz] = bindx[j];
         val[next_nz++] = val[j];
      }
   }
   for (i = 1; i <= N; i++) 
      bindx[i] = bindx[i-1] + bindx[i] - diag[i-1];
   

   for (i = 0 ; i < N ; i++) val[i] = 1./val[i];

   AZ_free(sum);
   AZ_free(mark);
   AZ_free(diag);
}
