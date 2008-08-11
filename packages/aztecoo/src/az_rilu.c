/*@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
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
