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

/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 * $Name$
 *====================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "az_aztec.h"

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
void AZ_space_for_padded_matrix(int overlap, int N_nonzeros, int N, 
    int *extra_rows, int *extra_nonzeros, int N_external, int *largest)
 
{
/****************************************************************************
  Estimate the number of additional rows and nonzeros due to overlapping.
  Currently, this estimate is based on the number of external variables
  and the number of nonzeros per row.

  Author:          Ray Tuminaro, SNL, 9222
 
  Return code:     void
  ============
 
  Parameter list:
  ===============
 
  overlap:         On input, 

                      == AZ_none: nonoverlapping domain decomposition
                      == AZ_diag: use rows corresponding to external variables
                                  but only keep the diagonal for these rows.
                      == k      : Obtain rows that are a distance k away from
                                  rows owned by this processor.
  
  N_nonzeros:      On input, number of nonzeros in the unpadded matrix.

  N:               On input, number of rows in the unpadded matrix.

  extra_rows:      On output, estimate of the number of additional rows 
                   needed for padding the matrix corresponding to 'overlap'.

  extra_nonzeros:  On output, estimate of the number of additional nonzeros
                   needed for padding the matrix corresponding to 'overlap'.
 
  N_external:      On input, number of external variables corresponding to
                   the unpadded matrix.

  largest:         On output, estimate of the maximum number of nonzeros
                   in any row due to overlapping.

****************************************************************************/
    int i;
    double new_exts, d_externs;
int temp;

    if ( (overlap == 0) || (overlap == AZ_diag) ) {
       *extra_rows     = N_external;
       *extra_nonzeros = N_external;
       *largest = 1;
    }
    else if (overlap >= 1) {

       /* we must estimate in this case */

       *extra_rows     = (int) 5.5*((double) (N_external*overlap));
d_externs = (double) N_external;
new_exts =  d_externs;
for (i = 2; i <= overlap; i++ ) {
   new_exts = new_exts + 4.*(sqrt(3.14159*new_exts)+3.14159);
                       /* This formula is based on the growth in */
                       /* the surface area of a sphere.          */
                       /*   S_0 = 4 pi r^2                       */
                       /*   S_1 = 4 pi (r+1)^2                   */
                       /*       = S_0 + 4 pi + 8 pi r            */
                       /*   substitute sqrt(S_0/(4 pi)) for r    */
                       /*   S_1 = S_0 + 4 ( pi + sqrt(S_0 pi))   */
   d_externs += new_exts;
}
*extra_rows = (int) d_externs;
*extra_rows = (*extra_rows)*2 + 30;

       if (N != 0) {
          temp = N_nonzeros/N;
          *extra_nonzeros = N + (*extra_rows)*temp;
          *largest        = 3.5*N_nonzeros/N;
          *extra_rows     += 25;
          *extra_nonzeros += 25;
          *largest        += 25;
       }
       else {
          *extra_rows     = 0;
          *extra_nonzeros = 0;
          *largest        = 0;
       }
    }
    else AZ_perror("Inproper level of overlap\n");
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
int AZ_adjust_N_nz_to_fit_memory(int N,int N_int_arrays, int N_dbl_arrays)
{
/****************************************************************************
  Find (and return) the largest value of k <= N such that we can 
  successfully allocate  N_int_arrays integer arrays of size k and 
  N_dbl_arrays double arrays of size k.

  Author:          Ray Tuminaro, SNL, 9222
 
  Return code:     int
  ============
 
  Parameter list:
  ===============
 
  N:               On input, the maximum number of integers and doubles
                   that we wish to try and allocate.
 */

   double **dptr;
   int    **iptr;
   int    i;

   iptr = (int **) AZ_allocate(N_int_arrays*sizeof(int *));
   dptr = (double **) AZ_allocate(N_dbl_arrays*sizeof(double *));

   if ( (dptr == 0) || (iptr == 0) ) 
      AZ_perror("ERROR: not enough memory for preconditioner.\n");

   for (i = 0 ; i < N_int_arrays ; i++ ) 
      iptr[i] = (int    *) AZ_allocate((N+20)*sizeof(int));
   for (i = 0 ; i < N_dbl_arrays ; i++ ) 
      dptr[i] = (double *) AZ_allocate((N+20)*sizeof(double));
                                   /* add a little extra */
                                   /* for manage memory  */
 
   /* Decrease memory until the problem fits */
 
   while ( (dptr[N_dbl_arrays-1] == NULL) || 
           (iptr[N_int_arrays-1] == NULL) ) {

      for (i = N_dbl_arrays-1 ; i >= 0; i-- ) 
         if (dptr[i] != NULL) AZ_free(dptr[i]);
      for (i = N_int_arrays-1 ; i >= 0; i-- ) 
         if (iptr[i] != NULL) AZ_free(iptr[i]);

      N = (int) ( ((double) N)*.91);
      if (N == 0) AZ_perror("ERROR: not enough memory for preconditioner.\n");
 
      for (i = 0 ; i < N_int_arrays ; i++ ) 
         iptr[i] = (int    *) AZ_allocate((N+20)*sizeof(int));
      for (i = 0 ; i < N_dbl_arrays ; i++ ) 
         dptr[i] = (double *) AZ_allocate((N+20)*sizeof(double));
   }
   for (i = N_dbl_arrays-1 ; i >= 0; i-- ) AZ_free(dptr[i]);
   for (i = N_int_arrays-1 ; i >= 0; i-- ) AZ_free(iptr[i]);
   AZ_free(dptr);
   AZ_free(iptr);
   return(N);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void AZ_combine_overlapped_values(int sym_flag,int data_org[],int options[],
	double x[], int map[], double ext_vals[], int name, 
	int proc_config[])
{
  /* Add the values that are redundant. That is, add the external values 
   * to the border values that correspond to them. This will make the    
   * operator symmetric if the incomplete factorization used above was   
   * symmetric.                                                          */

  int type, total, i, j, count, st, from, N_unpadded, N;
  MPI_AZRequest request[AZ_MAX_NEIGHBORS];  /* Message handle */
  double *little;
  double scale = .5;

  N_unpadded = data_org[AZ_N_internal] + data_org[AZ_N_border];
  N          = N_unpadded + data_org[AZ_N_external];

  if (sym_flag == AZ_symmetric) scale = 1.;
  else return;

  if (options[AZ_overlap] == 0) return;

  /* unshuffle the data */

  if (options[AZ_overlap] >= 1) {
     for (i = 0 ; i < N-N_unpadded ; i++ ) ext_vals[i] = x[i + N_unpadded];
     for (i = 0 ; i < N-N_unpadded ; i++ ) 
        x[i+N_unpadded] = ext_vals[map[i]-N_unpadded];
  }


  /* first send the external points to the neighbors */
 
  type            = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + 
                     AZ_MSG_TYPE;

  /* figure out longest message to be received and allocate space for it. */

  total = 0;
  for ( i = 0 ; i < data_org[AZ_N_neigh] ; i++ )
     total += data_org[AZ_send_length+i];

  little = (double *) AZ_manage_memory(total*sizeof(double), AZ_ALLOC, name,
                                             "temp in combine", &i);


  /* post receives */

  count = 0;
  for ( i = 0 ; i < data_org[AZ_N_neigh] ; i++ ) {
     from = data_org[AZ_neighbors+i];
     (void) mdwrap_iread((void *) &(little[count]),
                  sizeof(double)*data_org[AZ_send_length+i],
                  &from, &type, request+i);
     count += data_org[AZ_send_length+i];
  }

  /* send messages */

  count = data_org[AZ_N_internal] + data_org[AZ_N_border];
  for ( i = 0 ; i < data_org[AZ_N_neigh] ; i++ ) {
     (void) mdwrap_write((void *) &(x[count]), data_org[AZ_rec_length+i]*
                     sizeof(double), data_org[AZ_neighbors+i], type, &st);
     count += data_org[AZ_rec_length+i];
  }
 
  /* receive messages and add recvd values to the send list */
 
  count = 0;
  for ( i = 0 ; i < data_org[AZ_N_neigh] ; i++ ) {
     from = data_org[AZ_neighbors+i];
     (void) mdwrap_wait((void *) &(little[count]),
                  sizeof(double)*data_org[AZ_send_length+i],
                  &from, &type, &st,request+i);
     count += data_org[AZ_send_length+i];
  }
  for ( j = 0 ; j < total; j++ ) {
     x[ data_org[AZ_send_list+j] ] += little[j];
     x[ data_org[AZ_send_list+j] ] *= scale;
  }

}
