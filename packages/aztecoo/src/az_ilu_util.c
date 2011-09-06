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
#include "az_aztec.h"

  /* externals */

extern void sort_blk_col_indx(int num_blks_row, int *bindx_start_row,
                       int *ordered_index);

extern void sort2(int, int *, int*);

extern void order_parallel(int M, double *val_old, double *val_new, 
                   int *bindx_old, int *bindx_new, int *indx_old, int *indx_new,
                   int *bpntr_old, int *bpntr_new, int *diag_block);
extern void get_diag(int M, int *bindx, int *bpntr, int *diag_block);

/******************************************************************************/

void sort_blk_col_indx(int num_blks_row, int *bindx_start_row,
                       int *ordered_index)

/*******************************************************************************

  Routine to sort the block entires for a given block row into increasing order.

  Author:
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  num_blks_row:    Number of blocks in in this row.

  bindx_start_row: On Input:  starting address for the array of block column
                              indices for this row.
                   On Output: ordered block column indicies for this row.
                              (size >= num_blks_row).

  ordered_index:   On Input:  integer array.
                   On Output: an array of indices which describes the
                              reordering for the array "bindx_start_row"
                              (size >= num_blks_row).

*******************************************************************************/

{

  /* local variables */

  int i;


  /**************************** execution begins ******************************/

  /* Initialize the ordering index vector */

  for (i = 0; i < num_blks_row; i++) ordered_index[i] = i;

  /* Sort block column index array and produce the ordering index */

  sort2(num_blks_row, bindx_start_row-1, ordered_index-1);

} /* sort_blk_col_indx */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void sort2(int n, int ra[], int rb[])

/*******************************************************************************

  Numerical Recipes C source code modified to have first argument an integer
  array.

  Sorts the array ra[1,..,n] in ascending numerical order using heapsort
  algorithm, while making the corresponding rearrangement of the array
  rb[1,..,n].

  NOTE: The arrays start at 1 instead of 0, therefore you must pass call from C
  for a zero based array as:

                   sort(n, ra-1, rb-1);


  Author:          Modified by John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  n:               Length of arrays ra and rb.

  ra:              Array to be sorted.

  rb               Second array order acording to the sorted ra array.

*******************************************************************************/

{

  /* local variables */

  int l, j, ir, i;
  int rra;
  int rrb;

  /**************************** execution begins ******************************/

  l  = (n >> 1) + 1;
  ir = n;

  for (;;) {
    if (l > 1) {
      rra = ra[--l];
      rrb = rb[l];
    }
    else {
      rra    = ra[ir];
      rrb    = rb[ir];
      ra[ir] = ra[1];
      rb[ir] = rb[1];

      if (--ir <= 1) {
        ra[1] = rra;
        rb[1] = rrb;
        return;
      }
    }

    i = l;
    j = l << 1;

    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) ++j;
      if (rra < ra[j]) {
        ra[i] = ra[j];
        rb[i] = rb[j];
        j    += (i = j);
      }
      else j = ir + 1;
    }

    ra[i] = rra;
    rb[i] = rrb;
  }

} /* sort2 */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_order(int M, double *val_old, double *val_new, int *bindx,
              int *indx_old, int *indx_new, int *bpntr, int *diag_block)

/*******************************************************************************

  For each row, reorders the blocks of the matrix (indices and values) in
  increasing order and constructs the array of pointers: diag_block to the
  diagonal blocks.

  Returns diag_block[i] = -1 if no diagonal block has been found at the i'th
  row.

  Author:          Lydie Prevost, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  M:               The number of (block) rows in the matrix.

  val_old:         Array containing the entries of the matrix before reordering.
                   The matrix is stored block-row-by-block-row. Each block entry
                   is dense and stored by columns (VBR).

  val_new:         On output, array containing the entries of the matrix after
                   reordering. The matrix is stored block-row-by-block-row. Each
                   block entry is dense and stored by columns (VBR).

  bindx:           On input, contains the block column indices of the non-zero
                   block entries of the matrix before reordering.
                   On output, contains the block column indices of the non-zero
                   block entries of the matrix after reordering.

  indx_old:        The ith element of indx_old points to the location in val_old
                   of the (0,0) entry of the ith block entry. The last element
                   is the number of nonzero entries of matrix plus one.

  indx_new:        The ith element of indx_new points to the location in val_new
                   of the (0,0) entry of the ith block entry. The last element
                   is the number of nonzero entries of matrix plus one.

  bpntr:           The ith element of bpntr points to the first block entry of
                   the ith row in bindx. The last element is the number of
                   nonzero blocks of matrix plus one.

  diag_block:      On output, array of size M points on each diagonal block.

*******************************************************************************/

{

  /* local variables */

  int     i, kk, j, ii;
  int     num_blks_row;
  int    *sort, old_blk_index, counter;
  int     new_blk;
  int    *temp_ind, size_temp_ind = 10, size_temp_val = 40;
  double *temp_val;
  int     total_vals;
  int     start, end;

  /**************************** execution begins ******************************/

  temp_ind = (int    *) AZ_allocate(size_temp_ind*sizeof(int));
  temp_val = (double *) AZ_allocate(size_temp_val*sizeof(double));
  sort     = (int    *) AZ_allocate(sizeof(int) * (M));

  if ( (temp_val == NULL) || (sort == NULL)) 
     AZ_perror("Out of space inside AZ_sort()\n");

  for (i = 0; i < M; i++) diag_block[i] = -1;

  for (i = 0; i < M; i++) {     /* loop over the rows */

    /* constructs the new array of block column indices in this row */

    num_blks_row = bpntr[i+1] - bpntr[i];

    if (num_blks_row+1 > size_temp_ind) {
      size_temp_ind = num_blks_row + 1;
      AZ_free(temp_ind);
      temp_ind = (int *) AZ_allocate(size_temp_ind * sizeof(int));
    }

    for (ii = bpntr[i]; ii <= bpntr[i+1]; ii++)
      temp_ind[ii - bpntr[i]] = indx_old[ii];

    total_vals = indx_old[bpntr[i+1]] - indx_old[bpntr[i]];

    sort_blk_col_indx(num_blks_row, bindx+bpntr[i], sort);

    /* for each block of this row computes the new indices and constructs the
       pointers on the diagonal block i */

    indx_new[0] = indx_old[0];
    for (kk = 0; kk < num_blks_row; kk++) {
      new_blk = kk + bpntr[i];

      /* index into old VBR matrix and get the block size */

      indx_new[new_blk+1] = indx_new[new_blk] +
        (temp_ind[sort[kk]+1] - temp_ind[sort[kk]]);

      if (bindx[new_blk] == i) diag_block[i] = new_blk;
    }

    /* constructs the new array containing the entries of the matrix */

    if (total_vals > size_temp_val) {
      size_temp_val = total_vals;
      AZ_free(temp_val);
      temp_val = (double *) AZ_allocate(size_temp_val * sizeof(double));
    }

    start   = indx_old[bpntr[i]];
    end     = indx_old[bpntr[i+1]];
    counter = 0;

    for (ii = start; ii < end; ii++)
      temp_val[counter++] = val_old[ii];

    for (kk = 0; kk < num_blks_row; kk++) {

      /*
       * Get old block index into the coefficient array and new block location.
       */

      old_blk_index = temp_ind[sort[kk]] - temp_ind[0];
      new_blk       = kk + bpntr[i];

      counter = 0;
      for (j = indx_new[new_blk]; j < indx_new[new_blk+1]; j++) {
        val_new[j] = temp_val[old_blk_index + counter++];
      }
    }
  }

  AZ_free((void *) sort);
  AZ_free((void *) temp_ind);
  AZ_free((void *) temp_val);

} /* order */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void order_parallel(int M, double *val_old, double *val_new, int *bindx_old,
                    int *bindx_new, int *indx_old, int *indx_new,
                    int *bpntr_old, int *bpntr_new, int *diag_block)

/*******************************************************************************

  For each row, reorders the blocks of the matrix (indices and values) in
  increasing order and constructs the array of pointers: diag_block to the
  diagonal blocks.

  Returns diag_block[i]=-1  if no diagonal block has been found at the i'th row.

  Author:          Lydie Prevost, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  M:               The number of (block) rows in the matrix.

  val_old:         Array containing the entries of the matrix before reordering.
                   The matrix is stored block-row-by-block-row. Each block entry
                   is dense and stored by columns (VBR).

  val_new:         On output, array containing the entries of the matrix after
                   reordering. The matrix is stored block-row-by-block-row. Each
                   block entry is dense and stored by columns (VBR).

  bindx_old:       Contains the block column indices of the non-zero block
                   entries of the matrix before reordering.

  bindx_new:       On output, contains the block column indices of the reordered
                   non-zero block entries of the matrix before reordering.

  indx_old:        The ith element of indx_old points to the location in val_old
                   of the (0,0) entry of the ith block entry. The last element
                   is the number of nonzero entries of matrix plus one.

  indx_new:        The ith element of indx_new points to the location in val_new
                   of the (0,0) entry of the ith block entry. The last element
                   is the number of nonzero entries of matrix plus one.

  bpntr_old:       The i'th element of bpntr points to the first block entry of
                   the i'th row in bindx_old. The last element is the number of
                   nonzero blocks of matrix plus one.

  bpntr_new:       On, output, the i'th element of bpntr points to the first
                   block entry of the i'th row in bindx_new. The last element
                   is the number of nonzero blocks of matrix plus one.

  diag_block:      On output, array of size M points on each diagonal block.

*******************************************************************************/

{

  /* local variables */

  int  i, kk, j;
  int  num_blks_row_old, num_blks_row_new;
  int *sort, ptr, compt, temp;

  /**************************** execution begins ******************************/

  /* Allocate work space */

  sort = (int *) AZ_allocate(sizeof(int) * (M));
  if (sort == NULL) {
    (void) AZ_printf_err("Error: not enough memory inside order_parallel\n"
                   "       must run a smaller problem\n");
    exit(-1);
  }

  /* Initialize */

  for (i = 0; i < M; i++) diag_block[i] = -1;
  bpntr_new[0] = bindx_new[0] = 0;

  for (i = 0; i < M; i++) {     /* loop over the rows */

    /* constructs the new array of block column indices in this row */

    num_blks_row_old = bpntr_old[i+1] - bpntr_old[i];

    /* copy old block column index array and then sort it this defines the new
       block column index array */

    for (j = 0; j < num_blks_row_old; j++)
      bindx_new[bpntr_new[i] + j] = bindx_old[bpntr_old[i] + j];

    sort_blk_col_indx(num_blks_row_old, &bindx_new[bpntr_new[i]], sort);

    /* Count the blocks that multiply internal and border unknowns */

    num_blks_row_new = 0;
    for (j = 0; j < num_blks_row_old; j++) {
      if (bindx_new[bpntr_new[i] + j] >= M) break;

      num_blks_row_new++;
    }

    bpntr_new[i+1] = bpntr_new[i] + num_blks_row_new;

    /* for each block of this row compute the new index vector and construct the
       pointers to the diagonal block i */

    for (kk = bpntr_new[i]; kk < bpntr_new[i+1]; kk++) {

      /* Define new indx vector */

      if (kk - bpntr_new[i] == 0) {
        indx_new[0] = indx_old[0];
      }
      else {
        temp         = sort[kk-1-bpntr_old[i]] + bpntr_old[i];
        indx_new[kk] = indx_new[kk-1] + (indx_old[temp+1] - indx_old[temp]);
      }

      /* Get diagonal block pointers */

      if (bindx_new[kk] == i) diag_block[i] = kk;
    }

    /* constructs the new array containing the entries of the matrix */

    for (kk = bpntr_new[i]; kk < bpntr_new[i+1]; kk++) {
      ptr   = indx_old[sort[kk-bpntr_old[i]] + bpntr_old[i]];
      compt = -1;

      for (j = indx_new[kk]; j < indx_new[kk+1]; j++) {
        compt++;
        val_new[j] = val_old[ptr + compt];
      }
    }
  }

  AZ_free((void *) sort);

} /* order_parallel */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void get_diag(int M, int *bindx, int *bpntr, int *diag_block)

/*******************************************************************************

  Constructs the array of pointers: diag_block to the diagonal blocks.  Returns
  diag_block[i] = -1 if no diagonal block has been found at the ith row.

  Author:          Lydie Prevost, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  M:               Number of (block) rows in the matrix and L.

  bindx:           Contains the block column indices of the non-zero block
                   entries.

  bpntr:           The i'th element of bpntr points to the first block entry of
                   the i'th row in bindx. The last element is the number of
                   nonzero blocks of matrix plus one.

  diag_block:      On output, array of size M points on each diagonal block.

*******************************************************************************/

{
  int i, kk;

  for (i = 0; i < M; i++) diag_block[i] = -1;

  for (i = 0; i < M; i++) {
    for (kk = bpntr[i]; kk < bpntr[i+1]; kk++) {
      if (bindx[kk] == i)
        diag_block[i] = kk;
    }
  }

} /* get_diag */




#define add_scaled_row(row)                             \
    accum_col[Ncols++] = row;                           \
    for (kk = bindx[row] ; kk <= last[row]; kk++ ) {    \
       col = bindx[kk];                                 \
       if (col >= 0) accum_col[Ncols++] = col;          \
    }

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
void AZ_MSR_mult_patterns(int bindx[], int N, int last[], int bindx_length,
	int *accum_col)
{
/*
 * Multiply two MSR matrices (actually only the sparsity patterns)
 * and store the resulting sparse matrix. That is,
 *                 C = A * B
 *
 * Note: A, B and C are all stored in the same array (bindx). Specifically,
 *       all the elements in bindx[] correspond to A. However, some columns
 *       are encoded as negative numbers. To obtain the real column number
 *       one must do the following:  real_column = -2 - bindx[k]. The matrix,
 *       B, corresponds only to the positive column indices in bindx[].
 *
 * Parameters
 * ======
 *    bindx[]       On input, bindx[] holds the two input matrices as 
 *                  described above. On output, bindx[] holds the matrix
 *                  corresponding to the product of the two input matrices.
 *                  However, any matrix entry which was not in the matrix B
 *                  is encoded as a negative column number (see above).
 *
 *    N             On input, size of the matrices to multiply.
 *
 *    last          On input, uninitialized workspace of size N.
 *
 *    bindx_length  On input, the size of the array allocated for bindx.
 *                  In general, the number of elements in bindx[] will grow and
 *                  so bindx[] should contain additional room for this growth.
 *
 *    accum_col     On input, uninitialized workspace of size 2*N.
 */

   int    i, k, kk, next_nz, Ncols;
   int    row, col, *signs;
   int    start_row, end_row, first_one, orig_col;
   int largest_col, smallest_col;


   /* move off-diagonals to the back of the array */
   /* and  initialize start/end ptrs              */

   kk = bindx_length-1;
   end_row   = bindx[N]-1; 
   for (i = N-1 ; i >= 0 ; i--) {
      start_row = bindx[i];
      last[i] = kk;
      for (k = end_row; k >= start_row; k-- ) bindx[kk--] = bindx[k];
      end_row   = start_row - 1; 
      bindx[i] = kk+1;
   }

   /* initialize the arrays */

   for (i = 0 ; i < 2*N; i++) accum_col[i] = 0;
   signs = &(accum_col[N]);

   next_nz   = N+1;

   largest_col = 0;
   for (i = 0 ; i < N ; i++) {
      if (largest_col < i) largest_col = i;
      Ncols = 0;
      add_scaled_row(i);

      for (k = bindx[i] ; k <= last[i]; k++ ) {
         if (Ncols >= N) {
            AZ_sort(accum_col, Ncols, NULL, NULL);
            AZ_rm_duplicates(accum_col, &Ncols);
         }
         row        = bindx[k];
         if (row < 0) row = -2 - row;
         add_scaled_row(row);
      }
      AZ_sort(accum_col, Ncols, NULL, NULL);
      AZ_rm_duplicates(accum_col, &Ncols);

      for (k = 0 ; k < Ncols; k++) signs[accum_col[k]] = -1;

      /* compute the smallest and largest column */
      /* that could appear in a factorization    */

      smallest_col = i;
      first_one    = next_nz;
      if (bindx[i] <= last[i]) {
         orig_col = bindx[bindx[i]];
         if (orig_col < 0)  orig_col = -2 - orig_col; 
         if (smallest_col > orig_col) smallest_col = orig_col;

         orig_col = bindx[last[i]];
         if (orig_col < 0)  orig_col = -2 - orig_col; 
         if (largest_col < orig_col) largest_col = orig_col;
      }

      /* record sign of column to be stored  */

      for (k = bindx[i]; k <= last[i]; k++) {
         orig_col = bindx[k];
         if (orig_col >= 0)  signs[orig_col] = 1;    
      }

      if (next_nz+Ncols-2 > last[i]) {
         AZ_printf_err("Not enough room for the larger sparsity pattern\n");
         exit(1);
      }
      for (k = 0 ; k < Ncols ; k++) {
         col = accum_col[k];
         if (col != i) {
            if (signs[col] == -1) col = -2 - col;

            if ((accum_col[k] <= largest_col) && 
                (accum_col[k] >= smallest_col ) ) 
                   bindx[next_nz++] = col;
         }
      }
      bindx[i] = first_one;
      last[i]  = next_nz - 1;
   }
   bindx[N] = last[N-1]+1;

}
       
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
void AZ_rm_duplicates(int array[], int *N)
{
/*
 * remove any duplicates that might appear in the SORTED
 * array 'array'.
 *
 */
  int k, kk;

  kk = 0;
  for (k = 1; k < *N; k++) {
    if (array[kk] != array[k]) {
      kk++;
      array[kk] = array[k];
    }
  }
  if (*N != 0) kk++;

  *N= kk;
}

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

int AZ_fill_sparsity_pattern(struct context *context, int ifill, int bindx[], 
			     double val[], int N)
{
/*
 * Expand the MSR matrix (bindx,val) by adding zeros such that the new
 * matrix has the same sparsity pattern as when ILU(ifill) is performed
 * on the original matrix.
 *
 */
int   length, last_one, flag, i;
int   *work1, *work2;
double temp;

   length = context->N_nz_allocated;
   last_one = bindx[N]-1;

   /* allocate the work space arrays. If there is enough space in */
   /* val[] use that for one of the two work arrays               */

   if ( (length - last_one - 2)*sizeof(double) > (N+1)*sizeof(int) ) {
      flag = 0;
      work1 = (int *) &(val[last_one+1]);
   }
   else {
      work1 = (int *) AZ_allocate((N+1)*sizeof(int));
      flag = 1;
   }
   work2 = (int *) AZ_allocate(2*(N+1)*sizeof(int));
   if (work2 == NULL) AZ_perror("Out of space in ilu.\n");

   /* Take the power of the matrix. That is, A = A^ifill. */

   for (i = 0 ; i < ifill; i++)
      AZ_MSR_mult_patterns(bindx, N, work1, length,work2);

   AZ_free(work2);
   if (flag) AZ_free(work1);

   /* Move the nonzero values into their proper location in the  */
   /* new expanded matrix. Also, decode any columns which appear */
   /* as negative numbers (see AZ_MSR_mult_patterns()).          */

   for (i = bindx[N]-1; i >=  bindx[0]; i--) {
      if (bindx[i] >= 0) {
         temp            = val[last_one];
         val[last_one--] = 0.0;
         val[i]          = temp;
      }
      else {
         bindx[i] = -2 - bindx[i];
         val[i] = 0.0;
      }
   }

   return(bindx[N]);
}

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

void AZ_sort_msr(int bindx[], double val[], int N)
{
/*
 * Sort the column numbers within each MSR row.
 */

   int i, start, last;

   for (i = 0 ; i < N ; i++ ) {
      start = bindx[i];
      last  = bindx[i+1];
      AZ_sort( &(bindx[start]), last - start , NULL, &(val[start]));
   }
}
