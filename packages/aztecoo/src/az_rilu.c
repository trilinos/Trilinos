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
#include "az_aztec.h"

void AZ_fact_rilu(int N, int *nz_used, int *iu, int *iw, 
                  AZ_MATRIX *A_overlapped, double omega,
                  double rthresh, double athresh) 
{
   int i,j,k,ptr,start,end,kk,jw,M;
   int *msr_bindx;
   double s, *msr_val;
   

   msr_val   = A_overlapped->val;
   msr_bindx = A_overlapped->bindx;
   M         = msr_bindx[N];
   *nz_used  = M;

   /* Scale diagonal if necessary */
   if (rthresh==0.0) rthresh = 1.0; /* Default means no rthresh */
   if (rthresh != 1.0 || athresh != 0.0)
     for (i = 0 ; i < N ; i++ ) {
        if (msr_val[i] >=0)
           msr_val[i] = msr_val[i] * rthresh + athresh;
        else
           msr_val[i] = msr_val[i] * rthresh - athresh;
     }

   /* RILU factorization */

    for (i = 0 ; i < N ; i++ ) {
       start = msr_bindx[i];
       end   = msr_bindx[i+1];
       AZ_sort( &(msr_bindx[start]), end - start , NULL, &(msr_val[start]));
    }


      ptr = N+1;
      iu[0] = N;


      /* 
       Create the iu pointer to point on the last "l" element in msr_val 
       and create a working array iw to store the " aij" values for each 
       row i 
      */

      
      for (j = 0; j < N+1; j++) {
          iw[j] = 0;
      }
      for (i = 0; i < N; i++) {
          iw[i] = i;
          iu[i+1] = msr_bindx[i+1]-1;
          for ( j = msr_bindx[i]; j < msr_bindx[i+1]; j++ ) {
             iw [ msr_bindx[j] ] = ptr;
             if ( msr_bindx[j] < i ) {
                iu[i] = ptr;
             }
             ptr++;
          }

      /* do factorization */

          s = 0.0;
          for (k = msr_bindx[i]; k < iu[i]+1; k++) {

             /* special for the l part */

             kk = msr_bindx[k];
             msr_val[k] = msr_val[k] / msr_val [kk];
            
             for ( j = iu[kk]+1; j < msr_bindx[kk+1];j++ ) {

                jw = iw[msr_bindx[j]];
                if ( jw != 0 ) {
                   msr_val[jw] = msr_val[jw] - msr_val[k] * msr_val[j];
                }
                else{
                   s = s - msr_val[k] * msr_val[j];
                }

             }

          }

          /* compute the u_ii term */

          msr_val[i] = msr_val[i] + omega*s;

          for ( j = msr_bindx[i]; j < msr_bindx[i+1]; j++ ) {
              iw [ msr_bindx[j] ] = 0;
          }
          iw[i]=0;
       }

       for (i = 0; i < N; i++) msr_val[i] = 1/msr_val[i] ;

       /* set the pointers to fortran format */
       /* iu point to the first u element */

       for (i = 0; i < N+1; i++) {
           iu[i] = iu[i] + 2 ;
       }
       for (i = 0; i < M; i++) {
            msr_bindx[i] = msr_bindx[i] + 1 ;
       }
}
