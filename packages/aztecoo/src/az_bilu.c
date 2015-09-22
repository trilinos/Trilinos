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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "az_aztec.h"
#include "az_blas_wrappers.h"
#include "az_lapack_wrappers.h"

/**************************************************************/
/**************************************************************/
/**************************************************************/

void AZ_fact_bilu(int Nrows, AZ_MATRIX *matrix, int diag_block[],
	     int pivot[])
{
/*
 * Compute an incomplete block LU factorization. The resulting
 * LU factors overwrite the matrix.
 *
 * Paramaters
 * ==========
 * 
 * Nrows              On input, number of block rows in the system.
 * 
 * matrix             On input, VBR matrix to be factored.
 *                    On output, LU factors. Note: Instead of
 *                    storing U_ii, the LU factors of U_ii are
 *                    stored. Instead of storing U_ij (j>i),
 *                    the matrix inv(U_ii) U_ij is stored.
 *                    NOTE: IT IS ASSUMED THAT THE BLOCK COLUMNS
 *                    WITHIN EACH ROW ARE GIVEN IN ORDER.
 *
 * diag_block         On input, diag_block[i] indicates the location
 *                    in indx[] and bindx[] of row i's diagonal block.
 *
 * pivot              Vector on length cpntr[Nrows]. On output, pivot
 *                    contains pivoting information that is used
 *                    during the backsolves.
 *
 */
   int i, j, k, tj, tk, si, sj, sk;
   int *ipattern;
   int largest;
   double *Sub_block, *temp;
   int info;
   int pi,pj,pp;
   char N[2], T[2];
   double alpha, beta, *val;
   int    *bindx, *indx, *cpntr, *bpntr;

   bindx = matrix->bindx;
   indx  = matrix->indx;
   cpntr = matrix->cpntr;
   bpntr = matrix->bpntr;
   val   = matrix->val;
   N[0] = 'N'; N[1] = '\0';
   T[0] = 'T'; T[1] = '\0';

   /* figure out the largest possible block */
 
   largest = 0;
   for (i = 0 ; i < Nrows; i++) {
      si = cpntr[i+1] - cpntr[i];
      if (si > largest) largest = si;
   }
   largest = largest*largest;

   ipattern = (int *) AZ_allocate((Nrows+1)*sizeof(int));
   if (ipattern == NULL) {
      printf("Not enough space in bilu\n");
      exit(1);
   }
   for (i  = 0; i < Nrows ; i++) ipattern[i] = -1;

   Sub_block = (double *) AZ_allocate(largest*sizeof(double));
   if (Sub_block == NULL) {
      printf("Not enough space in bilu\n");
      exit(1);
   }


   for (i=0;  i < Nrows ; i++ ) {
      si = cpntr[i+1] - cpntr[i];

      /* store the nonzero pattern for this row */

      for (tj = bpntr[i]; tj < bpntr[i+1]; tj++ ) 
         ipattern[bindx[tj]] = indx[tj];


      for (tj = bpntr[i]; tj < bpntr[i+1]; tj++ ) {
         j = bindx[tj];
         if (i > j) {
           sj = cpntr[j+1] - cpntr[j];
           alpha = -1.0;
           beta = 1.;

           for (tk = bpntr[j]; tk < bpntr[j+1]; tk++) {
               k = bindx[tk];
               if ( (ipattern[k] != -1) && (k > j) ) {
                  sk = cpntr[k+1] - cpntr[k];
                  DGEMM_F77(CHAR_MACRO(N[0]), CHAR_MACRO(N[0]), &si, &sk, &sj, &alpha, 
                         &(val[indx[tj]]), &si, &(val[indx[tk]]), &sj, &beta, 
                         &(val[ipattern[k]]), &si);

               }
            }



            /* take transpose */

           temp = &(val[indx[tj]]);
           pp = 0;
           for (pi = 0 ; pi < si; pi++) {
              for (pj = 0 ; pj < sj; pj++) {
                 Sub_block[pp++] = temp[ pj*si + pi];
              }
           }
           DGETRS_F77(CHAR_MACRO(T[0]), &sj , &si, &(val[indx[diag_block[j]]]),
                   &sj, &(pivot[cpntr[j]]), Sub_block,
                   &sj, &info);

            /* transpose back */

           temp = &(val[indx[tj]]);
           pp = 0;
           for (pi = 0 ; pi < sj; pi++) {
             for (pj = 0 ; pj < si; pj++) {
                temp[pp++] = Sub_block[ pj*sj + pi];
             }
           }

         }
      }

      DGETRF_F77(&si, &si, &(val[indx[diag_block[i]]]), &si, &(pivot[cpntr[i]]), 
              &info);
      if (info > 0) {
          printf("Incomplete factorization yields singular subblock\n");
          printf("Can not use this factorization.\n");
          exit(1);
      }

      for ( tk = bpntr[i]; tk < bpntr[i+1]; tk++ ) {
         k = bindx[tk];
         if (k > i) {
            sk = cpntr[k+1] - cpntr[k];

            DGETRS_F77(CHAR_MACRO(N[0]), &si , &sk, &(val[indx[diag_block[i]]]),
                    &si, &(pivot[cpntr[i]]), &(val[indx[tk]]),
                    &si, &info);
         }
      }

      /* clear the nonzero pattern for this row */

      for (tj =bpntr[i]; tj < bpntr[i+1]; tj++ ) 
         ipattern[bindx[tj]] = -1;
   }
   AZ_free(Sub_block);
   AZ_free(ipattern);
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_lower_triang_vbr_solve(int Nrows, int cpntr[], int bpntr[], 
         int indx[], int bindx[], double val[], double b[]) {
/*******************************************************************************

  Lower triangular solver for Ly = b with L stored in VBR matrix format.

  Note: In this version the diagonal blocks of L are assumed to be identity.
  Return code:     void
  ============

  Parameter list:
  ===============

  Nrows            Number of (block) rows in the matrix and L.

  val, indx        VBR MATRIX representing L. 
  bindx, rpntr,
  cpntr, bpntr
                  
  diag_block:      On input, array of size Nrows points on each diagonal block.

  b:               On input, right hand side.
                   On output, contains the result vector.

*******************************************************************************/



   int i, j, i1, si, sj, tj, ione = 1;
   double minus_one = -1., one = 1.;
   char N[2];
           
   N[0] = 'N'; N[1] = '\0';
   for ( i = 0; i < Nrows ; i++ ) {
      si = cpntr[i+1] - cpntr[i];
      i1 = cpntr[i];

      for ( tj = bpntr[i]; tj < bpntr[i+1]; tj++ ) {
         j = bindx[tj];
         sj = cpntr[j+1] - cpntr[j];
         if (j < i) {
            DGEMV_F77(CHAR_MACRO(N[0]), &si, &sj, &minus_one, &(val[indx[tj]]), &si, 
                   &(b[cpntr[j]]), &ione, &one, &(b[i1]), &ione);
         }
      }
   }
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_upper_triang_vbr_solve(int Nrows, int cpntr[], int bpntr[], int indx[],
        int bindx[], double val[], double b[], int pivot[], int diag_block[]) 
{
/*
  Upper triangular solver for Uy = b with U stored in VBR matrix format.

  NOTE: Instead of storing U_ii, the LU factors of U_ii are stored. Instead 
        of storing U_ij (j>i), the matrix inv(U_ii) U_ij is stored.

  Return code:     void
  ============

  Parameter list:
  ===============

  Nrows            On input, number of (block) rows in the matrix and L.

  val, indx        On input, VBR matrix representing U where the LU factors
  bindx, rpntr,    of U_ii are stored instead of U_ii and the matrix 
  cpntr, bpntr     inv(U_ii) U_ij is stored instead of U_ij (j>i).

  diag_block:      On input, array of size Nrows points on each diagonal block.

  b:               On input, right hand side of linear system.
                   On output, contains the result vector, y.
*/
   int i, j, i1, si, sj, tj, ione = 1, info;
   double minus_one = -1., one = 1.;
   char N[2];
           
   N[0] = 'N'; N[1] = '\0';

   for (i = Nrows-1; i >= 0 ; i--) { 
      si = cpntr[i+1]  - cpntr[i];
      i1 = cpntr[i];

      DGETRS_F77(CHAR_MACRO(N[0]), &si , &ione, &(val[indx[diag_block[i]]]),
              &si, &(pivot[cpntr[i]]), &(b[i1]),
              &si, &info);

      for (tj = bpntr[i]; tj < bpntr[i+1]; tj++ ) {
         j  = bindx[tj];
         sj = cpntr[j+1] - cpntr[j];
         if (j > i) {
            DGEMV_F77(CHAR_MACRO(N[0]), &si, &sj, &minus_one, &(val[indx[tj]]), &si, 
                   &(b[cpntr[j]]), &ione, &one, &(b[i1]), &ione);
         }
      }
   }
}
extern void AZ_funswill(int *);

void AZ_funswill(int *trash)

{
   /* Just some garbage which does nothing. This will fool the */
   /* lint compiler so we don't have lint warnings.            */
   /* The purpose of this routine was to avoid some intel compiler problems */

   *trash += 1;    *trash -= 1;
}

