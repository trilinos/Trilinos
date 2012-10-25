/* ******************************************************************** */
static int rows_that_fit, rows_length, *rows, NBrows, end, start = 0;
static int subB_Nnz, Next_est, total_cols = 0;
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ml_struct.h"
#include "ml_rap.h"
#include "ml_memory.h"
#include "ml_aztec_utils.h"
#include "ml_mat_formats.h"
#include <limits.h>

/* ******************************************************************** */
/* matrix matrix multiplication                                         */ 
/* Cmatrix = Amatrix * Bmatrix                                          */
/*                                                                      */
/* Note: Matrices can be stored in chunks. See ml_rap.h for a           */
/*       description of the matrix structure ML_matrix.                 */
/* -------------------------------------------------------------------- */

#ifdef HAVE_ML_AZTECOO
/* ******************************************************************** */
void ML_blkmatmat_mult(ML_Operator *Amatrix, ML_Operator *Bmatrix,
		       ML_Operator **Cmatrix)
/* -------------------------------------------------------------------- */
{
  int    i,k, jj, next_nz, Ncols, N, Nnz_estimate, sub_i, accum_size, row;
  int    *Cbpntr, *Cbindx, *A_i_cols, rowi_N, *accum_col, row2_N;
  int    *Cindx, *Aindx;
  double *accum_val;
  double dtemp;
  double *StartAblock, *StartBblock;
  ML_Operator *current, *previous_matrix;
  int    max_nz_row_new = 0, total_nz = 0, index_length = 0;
  double A_avg_nz_per_row, B_avg_nz_per_row, estimated_nz_per_row;
  int    A_i_allocated;
  void   (*Agetrow)(ML_Operator *,int,int *,int *,int **,int **,int *,int);
  void   (*Bgetrow)(ML_Operator *,int,int *,int *,int **,int **,int *,int);
  int    *B_indx;
  int    *Bptr, *Bcols;
  double *Bvals;
  int *col_inds, B_total_Nnz, itemp,itemp2, B_allocated, B_allocated_int, hash_val, *accum_index, lots_of_space;
  int rows_that_fit, rows_length, *rows, NBrows, end, start = 0;
  int subB_Nnz, Next_est, total_cols = 0;
  /*
    double t1, t6;
  */
  int tcols, hash_used, j, *tptr;
  int *acc_col_ptr, *Bcol_ptr; 
  double  *acc_val_ptr;
  int     *Bval_ptr;
  int allzeros;
  double *Avalues;
  double *Bvalues;
  int save_ints[6];
  struct ML_vbrdata *A_VBR, *B_VBR;
  int NrowsPerBlock, NcolsPerBlock, Nghost = 0, InnerDim,iii,jjj,kkk;
  int LargestRowsPerBlock = 1, NnzPerBlock, RowOffset, RowOffsetBlocks, next_value;
  double sum, *Cvalues;
  double *tmp_val_ptr;
  int    *Ccpntr, *Crpntr, oldstart = 0;
  struct ML_vbrdata     *Cvbr_mat;
  int hashTableIsPowerOfTwo = 0;
  int nearbyIndex;
  int blocks;


/*  printf("This is an experimental routine. It basically works but ...\n");
  printf("there are some things hardwired and some thing would need\n");
  printf("to be added for both post and pre processing\n");
  printf("Here is a list:\n");
  printf("    a) This routine allows the left matrix to have variable\n");
  printf("       block heights and widths. It also allows the right matrix\n");
  printf("       to have varying heights (but the block widths need to be\n");
  printf("       constant). This is good for smoothed aggregation because\n");
  printf("       normally Ptent has constant block widths. The main bad\n");
  printf("       thing is that ML has no way to specify all of this\n");
  printf("       We simply use A->num_PDEs for all of them\n");
  printf("    b) We would need a routine to convert crs matrices to \n");
  printf("       vbr matrices. This routine would probably also have \n");
  printf("       to work with submatrices created by ML_exch_row().\n");
  printf("       Most likely this conversion would occur after ML_exch_row\n");
  printf("       as we don't want to alter that routine.\n");
  printf("    c) We would need some form of back_to_local() that works\n");
  printf("       with VBR matrices and could map global ids back to local\n");
  printf("    d) Right now the destroy is hardwired in things like \n");
  printf("       ML_2matmult() for this MSR_CSR_recursive? I'm not sure\n");
  printf("       why this is hardwired and not simply kicked off by the\n");
  printf("       operators destroy function. Anyway, we would need an\n");
  printf("       appropriate destroy for VBR with recursion. \n");
*/

  save_ints[0] = Amatrix->getrow->Nrows;
  save_ints[1] = Amatrix->outvec_leng;  
  save_ints[2] = Amatrix->invec_leng;   
  save_ints[3] = Bmatrix->getrow->Nrows;
  save_ints[4] = Bmatrix->outvec_leng;  
  save_ints[5] = Bmatrix->invec_leng;   

  /**************************************************************************/
  /* figure out if a cheaper getrow function (either ML_get_matrow_CSR or   */
  /* ML_get_row_CSR_norow_map) can be used avoiding some function overhead. */
  /* ---------------------------------------------------------------------- */

  Agetrow  =  ML_get_matrow_VBR;
  Bgetrow  =  ML_get_matrow_VBR;
  A_VBR    = (struct ML_vbrdata    *) Amatrix->data;
  Avalues  = A_VBR->val;
  B_VBR    =  (struct ML_vbrdata *)    Bmatrix->data;
  Bvalues  = B_VBR->val;

  /* Nghost is not right as it does not consider ghosts        */
  /* associated with ML_exch_row() ... but better than nothing */

  Nghost   = ML_CommInfoOP_Compute_TotalRcvLength(Bmatrix->getrow->pre_comm);

  /* Reset various lengths to reflect block quantities   */
  /*                                                     */
  /* Actually, it would be better if we really had these */
  /* quantities available as this routine does not need  */
  /* to force constant block sizes for the left matrix   */
  /* nor a constant height for the right matrix.         */


  N = Amatrix->getrow->N_block_rows;

  if (Bmatrix->invec_leng+Nghost > 0)
    NcolsPerBlock = B_VBR->cpntr[1] - B_VBR->cpntr[0];

  /********************************************************/
  /* put here to alleviate a hidden error - C. Tong       */
  if ( Amatrix->max_nz_per_row > Amatrix->getrow->N_block_rows )
    Amatrix->max_nz_per_row = Amatrix->getrow->N_block_rows;
  if ( Bmatrix->max_nz_per_row > Bmatrix->getrow->N_block_rows )
    Bmatrix->max_nz_per_row = Bmatrix->getrow->N_block_rows;
  /********************************************************/

  /*********************************************************/
  /* Count the number of external variables. Unfortunately */
  /* ML_exch_row() does not properly put the correct       */
  /* communication information ... and so this is just an  */
  /* estimate of how many externals there are.             */
  /* ------------------------------------------------------*/

  current = Bmatrix;
  Next_est = 0;
  while(current != NULL) {
    if (current->getrow->pre_comm != NULL) {
      if (current->getrow->pre_comm->total_rcv_length <= 0) {
	ML_CommInfoOP_Compute_TotalRcvLength(current->getrow->pre_comm);
      }
      Next_est += current->getrow->pre_comm->total_rcv_length;
    }
    current = current->sub_matrix;
  }

  /* Estimate the total number of columns in Bmatrix. Actually,  */
  /* we need a crude estimate. We will use twice this number for */
  /* hashing the columns. Once we are finished hashing columns   */
  /* we will count the number of columns and revise the estimate */
  /* for the number of externals and increase the size of the    */
  /* accumulator if needed.                                      */
  /* ----------------------------------------------------------- */

  accum_size = Bmatrix->invec_leng + 75;
  accum_size += 3*Next_est;
  if (Bmatrix->N_total_cols_est + 75 > accum_size ) {
    accum_size = Bmatrix->N_total_cols_est + 75;
  }

  /* the '3' above is a total kludge. the problem is that we   */
  /* need to know the number of neighbors in the external rows */
  /* and this is annoying to compute */

  /**************************************************************************/
  /* allocate space to hold a single row of Amatrix and the accumulator     */
  /* which contains the current row in the resulting matrix.                */
  /*------------------------------------------------------------------------*/

  A_i_allocated = Amatrix->max_nz_per_row + 1;
  A_i_cols  = (int    *) ML_allocate(A_i_allocated * sizeof(int) );
  Aindx     = (int *) ML_allocate(A_i_allocated * sizeof(int));
  accum_col = (int    *) ML_allocate( accum_size * sizeof(int) );
  accum_val = (double *) ML_allocate(accum_size*NcolsPerBlock*sizeof(double));
  if ( (Aindx == NULL) || (accum_val == NULL)) {
    printf("Not enough space in ML_matmatmult().\n");
    printf("trying to allocate %d %d elements \n",A_i_allocated,accum_size);
    printf("Left  matrix has %d rows \n", Amatrix->getrow->N_block_rows);
    printf("Right matrix has %d rows \n", Bmatrix->getrow->N_block_rows);
    printf("Left  matrix has %d nz per row\n", Amatrix->max_nz_per_row);
    printf("Right matrix has %d nz per row\n", Bmatrix->max_nz_per_row);
    exit(1);
  }

  /**************************************************************************/
  /* Make conservative estimates of the size needed to hold the resulting   */
  /* matrix. NOTE: These arrays can be increased later in the computation   */
  /*------------------------------------------------------------------------*/

  if ( Amatrix->getrow->N_block_rows > 0 ) {
    row = 0;
    Agetrow(Amatrix,1, &row, &A_i_allocated , &A_i_cols, &Aindx, &i,0);
    row = (Amatrix->getrow->N_block_rows-1)/2;
    Agetrow(Amatrix,1, &row, &A_i_allocated , &A_i_cols, &Aindx, &k,0);
    row = Amatrix->getrow->N_block_rows-1;
    Agetrow(Amatrix,1, &row, &A_i_allocated , &A_i_cols, &Aindx, &jj,0);
    A_avg_nz_per_row = ((double) (i+k+jj))/3.0;
  } else A_avg_nz_per_row = 100;

  if ( Bmatrix->getrow->N_block_rows > 0 ) {
    row = 0;
    Bgetrow(Bmatrix,1,&row, &accum_size, &accum_col, (int **) &accum_val, &i,0);
    row = (Bmatrix->getrow->N_block_rows-1)/2;
    Bgetrow(Bmatrix,1,&row, &accum_size, &accum_col, (int **) &accum_val, &k,0);
    row = Bmatrix->getrow->N_block_rows-1;
    Bgetrow(Bmatrix,1,&row, &accum_size, &accum_col, (int **) &accum_val,&jj,0);
    B_avg_nz_per_row = ((double) (i+k+jj))/3.0;
  } else B_avg_nz_per_row = i = k = jj = 100;

  if (i  > Bmatrix->max_nz_per_row) Bmatrix->max_nz_per_row =  i+10;
  if (k  > Bmatrix->max_nz_per_row) Bmatrix->max_nz_per_row =  k+10;
  if (jj > Bmatrix->max_nz_per_row) Bmatrix->max_nz_per_row = jj+10;

  estimated_nz_per_row = sqrt(A_avg_nz_per_row) + sqrt(B_avg_nz_per_row) - 1.;
  estimated_nz_per_row *= estimated_nz_per_row;
  Nnz_estimate = (int) (((double) Amatrix->getrow->N_block_rows) *
			estimated_nz_per_row*.75) + 100 ;
  if (Nnz_estimate <= Bmatrix->max_nz_per_row)
    Nnz_estimate = Bmatrix->max_nz_per_row + 1;

  /**************************************************************************/
  /* Allocate space for the matrix. Since 'Nnz_estimate' is just an estimate*/
  /* of the space needed, we will reduce 'Nnz_estimate' if we are not       */
  /* successful allocating space.                                           */
  /*------------------------------------------------------------------------*/

  Cbpntr     = (int    *) ML_allocate((N+1)* sizeof(int) );
  Cindx = NULL; Cbindx = NULL;
  while ( ((Cindx == NULL)||(Cbindx==NULL) || (Cvalues == NULL)) && 
	  (Nnz_estimate > Bmatrix->max_nz_per_row) ) {
    if (Cbindx  != NULL) ML_free(Cbindx);
    if (Cindx != NULL) ML_free(Cindx);
    Cbindx   = (int    *) ML_allocate( Nnz_estimate* sizeof(int) );
    Cindx  = (int    *) ML_allocate( Nnz_estimate* sizeof(int));
    Cvalues= (double *) ML_allocate( Nnz_estimate*
				     NcolsPerBlock*NcolsPerBlock* sizeof(double));
    if (Cvalues == NULL) Nnz_estimate = Nnz_estimate/2;
  }
  if ( ((Cindx == NULL)||(Cbindx==NULL)||(Cvalues==NULL)) && (N != 0)) {
    printf("Not enough space for new matrix in ML_matmatmult().\n");
    printf("trying to allocate %d elements \n",Nnz_estimate);
    printf("Left  matrix has %d rows \n", Amatrix->getrow->N_block_rows);
    printf("Right matrix has %d rows \n", Bmatrix->getrow->N_block_rows);
    printf("Left  matrix has %d nz per row\n", Amatrix->max_nz_per_row);
    printf("Right matrix has %d nz per row\n", Bmatrix->max_nz_per_row);
    exit(1);
  }
  next_value = 0;


  /**************************************************************************/
  /* Make a copy of Bmatrix and use this copy to speed up the computation.  */
  /*   1) First count the total number of nonzeros                          */
  /*   2) allocate space for Bmatrix copy.                                  */
  /*      Note: accum_index is used for hashing. We make it twice the size  */
  /*            of the accumulator for good performance.                    */
  /*   3) grab the matrix. Put the nonzero values in Bcols & B_indx.         */
  /*   4) Using a hash table, renumber the column indices so that they are  */
  /*      in [0:index_length-1]. Store the original column indices in       */
  /*      col_inds[].                                                       */
  /*------------------------------------------------------------------------*/

  if (Bmatrix->N_nonzeros <= 0) {
    B_total_Nnz = 0;
    for (i = 0; i < Bmatrix->getrow->N_block_rows; i++ ) {
      Bgetrow(Bmatrix,1,&i,&accum_size,&accum_col,(int **) &accum_val,&row2_N,0);
      B_total_Nnz += row2_N;
    }
  }
  else B_total_Nnz = Bmatrix->N_nonzeros;
  B_total_Nnz++; /* avoid possible division by zero below */

  /* ******************************************************************
     Try to make the hash table size a power of two.  Doing so makes
     hash table lookup more efficient by using the bitwise AND (&) operator
     and avoiding the modulus (%) function.

     The new lookup is ML_fast_hash, the original is ML_hash_it.  All hashing
     is located in ml/src/Utils/ml_utils.[ch].

     First check for and handle potential overflow:

           Case 1:  If accum_size is larger than 2^29, but not by much,
                    use a hash table size of 2^30.
           Case 2:  If accum_size is larger than 2^29 but less than 2^30,
                    use a hash table size of 2*accum_size.
           Case 3:  If accum_size is larger than 2^30 but not by much,
                    use a hash table size of 2^30-1 (INT_MAX).
           Case 4:  If accum_size is much larger than 2^30,
                    issue an error and abort.

     Otherwise, choose the hash table to be 2^k, where 2^k > 2*accum_size,
     k is a positive integer, and k < 31.

     Note:
     This potentially wastes some space if accum_size is close to but not
     exactly a power of 2.

     Some useful constants:

       2^31-1 = 2147483647 = INT_MAX
       2^30   = 1073741824
       2^29   =  536870912

  ****************************************************************** */

  if ( accum_size > 536870912 ) {
    if ( 1.8*accum_size < 1073741824) {
      index_length = 1073741824;
      hashTableIsPowerOfTwo = 1;
    } else if (accum_size < 1073741824) {
      index_length = 2*accum_size;
      hashTableIsPowerOfTwo = 0;
    } else if (1.8*accum_size < INT_MAX) {
      index_length = INT_MAX;
      hashTableIsPowerOfTwo = 0;
    } else {
      printf("Cannot create hash table for accumulator of size %d\n",
             accum_size);
      exit(1);
    }
  } else {
    index_length = 2;
    while (index_length < 2*accum_size) {
      index_length *= 2;
    }
    hashTableIsPowerOfTwo = 1;
  }

  accum_index  = (int *) ML_allocate(index_length*sizeof(int));
  col_inds     = (int *) ML_allocate(index_length*sizeof(int));

  Bptr      = (int    *) ML_allocate( (Bmatrix->getrow->N_block_rows+1)*sizeof(int));
  B_allocated = B_total_Nnz * 2;
  B_allocated_int = Bmatrix->blocks + 1;
  B_allocated_int = B_allocated_int*2;
  B_allocated_int =  B_total_Nnz;  /*this line needs to go sometime and be replaced by those above it*/
  lots_of_space = 0;

  Bcols     = NULL; B_indx = NULL; Bvals = NULL;
  while (Bvals == NULL || B_indx == NULL) {
    lots_of_space++;
    if (Bcols != NULL) ML_free(Bcols);
    if (B_indx != NULL) ML_free(B_indx);
    Bcols     = (int    *) ML_allocate( B_allocated_int * sizeof(int));
    B_indx    = (int    *) ML_allocate( B_allocated_int * sizeof(int));
    Bvals    = (double *) ML_allocate( B_allocated * sizeof(double));
    B_allocated /= 2;
    B_allocated_int /=2; 
  }
  ML_free(B_indx); ML_free(Bcols);

  Bcols     = (int    *) ML_allocate( B_allocated_int * sizeof(int));
  B_indx    = (int    *) ML_allocate( B_allocated_int * sizeof(int));
  Bvals    = (double *) ML_allocate( B_allocated * sizeof(double));
  if (lots_of_space != 1) lots_of_space = 0;

  if (Bmatrix->getrow->N_block_rows != 0) {
    dtemp = ((double) B_allocated)/ ( (double) B_total_Nnz);
    dtemp *= .9;
    rows_that_fit = (int) (( (double) Bmatrix->getrow->N_block_rows) * dtemp);
    if (rows_that_fit == 0) rows_that_fit++;
  }
  else rows_that_fit = 1;

  *Cmatrix        = NULL;
  previous_matrix = NULL;
  sub_i           = 0;
  next_nz         = 0;
  Cbpntr[0]        = next_nz;

  end = N;
  start = 0;

  hash_used = 0;
  RowOffset = 0;
  RowOffsetBlocks = 0;
  while (start < N) {

    itemp = 0;
    itemp2 = 0;
    Bptr[0] = 0;
    if (lots_of_space) {
      for (i = 0; i < Bmatrix->getrow->N_block_rows; i++ ) {
	    VBR_block_getrow(Bmatrix, i, &B_allocated_int, &B_allocated, &blocks, &Bcols, &B_indx, &Bvals, &row2_N, itemp, itemp2);
	    itemp += blocks;
        itemp2 += row2_N;
	    Bptr[i+1] = itemp;
      }
      subB_Nnz = itemp;
    }

    if (!lots_of_space) {
      NBrows = 0;
      rows_length = rows_that_fit*2 + 3;
      rows = (int *) ML_allocate(rows_length*sizeof(int));
      for (i=0; i < rows_length; i++) rows[i] = -1;
      ML_determine_Bblkrows(start, &end, Amatrix,
			 &rows, &rows_length, &NBrows,
			 &rows_that_fit, Agetrow);

      jj = 0;
      itemp = 0;
      ML_az_sort(rows, NBrows, NULL, NULL);
      for (i = 0; i < Bmatrix->getrow->N_block_rows; i++ ) {
	    if (rows[jj] == i) {
	      jj++;
	      k = B_allocated;
	      VBR_block_getrow(Bmatrix, i, &B_allocated_int, &B_allocated, &blocks, &Bcols, &B_indx, &Bvals, &row2_N, itemp, itemp2);
    	  itemp += blocks;
          itemp2 += row2_N;
    	}
    	Bptr[i+1] = itemp;
      }
      subB_Nnz = itemp;
    }

    for (i = 0; i < index_length; i++) accum_index[i] = -1;
    for (i = 0; i < index_length; i++) col_inds[i] = -1;

    i = 0;
    if (Bmatrix->getrow->pre_comm != NULL) 
      i = Bmatrix->getrow->pre_comm->total_rcv_length;

    /* Record a "nearby" index for use as a column index in any empty rows. */
    /* This avoids creating a communication pattern in which a single       */
    /* processor must communicate with all others in                        */
    /* ML_CommInfoOP_GenUsingGIDExternals().                                */

    if (subB_Nnz > 0) nearbyIndex = Bcols[0];
    else              nearbyIndex = -1;

    tcols = 0;
    hash_used = 0;
    if (hashTableIsPowerOfTwo)
    {
      int ilm1 = index_length-1;
      for (i = 0; i < subB_Nnz; i++) {
        if (hash_used >= ((int) (.5 * index_length)) ) {
          ML_free(accum_index);
          if (index_length >= 1073741824) /* 2^30 */
            pr_error("Exceeded largest possible hash table size\n");
          accum_index = (int *) ML_allocate(sizeof(int)*2*index_length);
          tptr = (int *) ML_allocate(sizeof(int)*2*index_length);
          if (tptr == NULL) pr_error("ML_matmat_mult: out of tptr space\n");
          for (j = 0; j < 2*index_length; j++) accum_index[j] = -1;
          for (j = 0; j < 2*index_length; j++) tptr[j] = -1;
          for (j = 0; j < i; j++)
            Bcols[j] = col_inds[Bcols[j]];
          ML_free(col_inds);  col_inds = tptr;
          tcols = 0;
          hash_used = 0;
          index_length *= 2;
          ilm1 = index_length-1;
          for (j = 0; j < i; j++) {
            ML_fast_hash(Bcols[j], col_inds,ilm1,&hash_used, &hash_val);
            if (col_inds[hash_val] == -1) tcols++;
            col_inds[hash_val] = Bcols[j];
            Bcols[j] = hash_val;
          }
        }
         
        ML_fast_hash(Bcols[i], col_inds, ilm1, &hash_used, &hash_val);
        if (col_inds[hash_val] == -1) tcols++;
        col_inds[hash_val] = Bcols[i];
        Bcols[i] = hash_val;
      } /*for (i = 0; i < subB_Nnz; i++)*/
    } 
    else {
      /* hash table size is not 2^k */
      for (i = 0; i < subB_Nnz; i++) {
        if (hash_used >= ((int) (.5 * index_length)) ) {
          ML_free(accum_index);
          accum_index = (int *) ML_allocate(sizeof(int)*2*index_length);
          tptr = (int *) ML_allocate(sizeof(int)*2*index_length);
          if (tptr == NULL) pr_error("ML_matmat_mult: out of tptr space\n");
          for (j = 0; j < 2*index_length; j++) accum_index[j] = -1;
          for (j = 0; j < 2*index_length; j++) tptr[j] = -1;
          for (j = 0; j < i; j++)
            Bcols[j] = col_inds[Bcols[j]];
          ML_free(col_inds);  col_inds = tptr;
          tcols = 0;
          hash_used = 0;
          index_length *= 2;
          for (j = 0; j < i; j++) {
            ML_hash_it(Bcols[j], col_inds, index_length, &hash_used, &hash_val);
            if (col_inds[hash_val] == -1) tcols++;
            col_inds[hash_val] = Bcols[j];
            Bcols[j] = hash_val;
          }
        }
         
        ML_hash_it(Bcols[i], col_inds, index_length, &hash_used, &hash_val);
        if (col_inds[hash_val] == -1) tcols++;
        col_inds[hash_val] = Bcols[i];
        Bcols[i] = hash_val;
      } /*for (i = 0; i < subB_Nnz; i++) */
    } /* if (hashTableIsPowerOfTwo) */

    /* Count how many columns there are, record the number of   */
    /* columns and change the size of the accumulator if needed */

    jj = 0;
    for (i = 0; i < index_length; i++) { 
      if (col_inds[i] == -1) jj++;
    }
    total_cols += (index_length - jj);
    if (accum_size < total_cols) {
      ML_free(accum_col);
      ML_free(accum_val);
      accum_size = total_cols + 1;
      accum_col = (int *) ML_allocate(accum_size*sizeof(int));
      accum_val= (double *) ML_allocate(accum_size*NcolsPerBlock*sizeof(double));
                                 /* does this need to be this large? */
      if (accum_val == NULL) 
	    pr_error("ML_matmat_mult: no room for accumulator\n");
    }



    /**************************************************************************/
    /* Perform the matrix-matrix multiply operation by computing one new row  */
    /* at a time.                                                             */
    /*------------------------------------------------------------------------*/

    Crpntr = (int *) ML_allocate(sizeof(int)*(end-start+2));
    Crpntr[0] = 0;
    NrowsPerBlock = -1;
    /*loop over matrix block rows*/
    for (i = start; i < end ; i++) {
      NrowsPerBlock = A_VBR->rpntr[i+1] - A_VBR->rpntr[i];
      Crpntr[i-start+1] = Crpntr[i-start] + NrowsPerBlock;
      NnzPerBlock = NrowsPerBlock*NcolsPerBlock;
      if (NrowsPerBlock > LargestRowsPerBlock) {
	    LargestRowsPerBlock = NrowsPerBlock;
	    ML_free(accum_val);
	    accum_val = (double *) ML_allocate(accum_size*LargestRowsPerBlock*
					   NcolsPerBlock*sizeof(double));
      }
      Ncols = 0;
      acc_col_ptr = accum_col;
      acc_val_ptr = accum_val;
      Agetrow(Amatrix,1, &i, &A_i_allocated, &A_i_cols, &Aindx, &rowi_N, 0); /*Agetrow is a pointer to a functioni that gets the block row*/

      /***********************************************************************/
      /* Take each column entry in the ith row of A and find corresponding   */
      /* row in B. Muliply each B row and sum it into the accumulator.       */
      /*---------------------------------------------------------------------*/

      /*iterate on the blocks of the blockrow*/
      for (k = 0; k < rowi_N; k++) {
	    StartAblock = &(Avalues[(int) Aindx[k]]);
    	InnerDim = (A_VBR->cpntr)[A_i_cols[k]+1] - (A_VBR->cpntr)[A_i_cols[k]]; /*this may not work in parallel due to block orderings*/
  	    jj = Bptr[A_i_cols[k]];
	    Bval_ptr = &(B_indx[jj]);
	    Bcol_ptr = &(Bcols[jj]);
        /*multiply a block of a row of A times a block row of B*/
	    while (jj++ < Bptr[A_i_cols[k]+1]) {
	      StartBblock = &(Bvals[(int) *Bval_ptr]); 
          /*we have a new block*/
	      if (accum_index[*Bcol_ptr] < 0) {
	        *acc_col_ptr = *Bcol_ptr; acc_col_ptr++;
	        for (iii = 0; iii < NcolsPerBlock; iii++) { /*this should be a GEMM call or at least tested as one*/
	          for (jjj = 0; jjj < NrowsPerBlock; jjj++) {
		        sum = 0.;
		        for (kkk = 0; kkk < InnerDim; kkk++) {
		          sum += StartAblock[jjj+NrowsPerBlock*kkk]*
		          StartBblock[kkk+InnerDim*iii];
		        }
		        *acc_val_ptr = sum; acc_val_ptr++;
	          }
            }
	        accum_index[*Bcol_ptr] = Ncols++;
	      }
	      else {
	      /* take accum_index and simply multiply it by the block size */
	        tmp_val_ptr =  &(accum_val[accum_index[*Bcol_ptr]*NnzPerBlock]);
	        for (iii = 0; iii < NcolsPerBlock; iii++) {
	          for (jjj = 0; jjj < NrowsPerBlock; jjj++) {
		        sum = 0.;
  	    	    for (kkk = 0; kkk < InnerDim; kkk++) {
		          sum += StartAblock[jjj+NrowsPerBlock*kkk]*
		          StartBblock[kkk+InnerDim*iii];
		        }
		        *tmp_val_ptr += sum; tmp_val_ptr++;
	          }
            }
	      }
	      Bcol_ptr++; Bval_ptr++;
	    }
      } /***********************************************************************/
        /* Convert back to the original column indices.                        */
        /*---------------------------------------------------------------------*/

      acc_col_ptr = accum_col;
      for (jj = 0; jj < Ncols; jj++ ) {
	    accum_index[*acc_col_ptr] = -1;
	    *acc_col_ptr = col_inds[*acc_col_ptr];
	    acc_col_ptr++;
      }

      /* check if we have enough space to store the new matrix row */
      /* If not, place the current matrix into a submatrix and     */
      /* allocate new vectors (hopefully large enough to hold the  */
      /* rest of the matrix) and a 'parent' matrix to hold the     */
      /* remaining rows.                                           */

      if (next_nz+Ncols > Nnz_estimate) {
   	  /* create sub_matrix object */

	    total_nz += next_nz;

	    Cvbr_mat = (struct ML_vbrdata *) ML_allocate(sizeof(struct ML_vbrdata));
    	Cvbr_mat->bindx       = Cbindx;
	    Cvbr_mat->val         = Cvalues;
    	Cvbr_mat->bpntr       = Cbpntr;
	    Cvbr_mat->indx        = Cindx; 
     
        Ccpntr = (int *) ML_allocate(sizeof(int)*hash_used);
        Ccpntr[0] = 0;
        for (kkk = 1; kkk <= hash_used; kkk++) 
           Ccpntr[kkk] = Ccpntr[kkk-1] + NcolsPerBlock;
    	Cvbr_mat->cpntr       = Ccpntr;
    	Cvbr_mat->rpntr       = Crpntr;
    	current = ML_Operator_Create(Amatrix->comm);
	    ML_Operator_Set_1Levels(current, Bmatrix->from, Amatrix->to);

                                                                                                                 
        if(current->ML_id != ML_ID_OP ) {
          printf("ML_Operator_Set_ApplyFunc error : wrong object.\n");
          exit(-1);
        }
        /* newly added : 8/17/00 */
        if (current->data != NULL && current->data_destroy != NULL )
        {
          current->data_destroy(current->data);
          current->data = NULL;
        }
        current->invec_leng = save_ints[5];
        current->outvec_leng = save_ints[1];
        current->data = (void*)Cvbr_mat;
        current->getrow->func_ptr = VBR_getrows;
                                                                                                                 
        current->matvec->ML_id = ML_NONEMPTY;
        if(previous_matrix != NULL)
        {
          current->matvec->Nrows = RowOffset +previous_matrix->getrow->Nrows;
          current->getrow->Nrows = RowOffset + previous_matrix->getrow->Nrows;
          current->getrow->N_block_rows = RowOffsetBlocks + previous_matrix->getrow->N_block_rows;
        }
        else
        {
          current->matvec->Nrows = RowOffset;
          current->getrow->Nrows = RowOffset;
          current->getrow->N_block_rows = RowOffsetBlocks;
        }

	    /*ML_Operator_Set_Getrow(current, RowOffset, 
                               az_vbrgetrow_wrapper);*/


        /* current->data_destroy = AZ_ML_FullClean; */

        Cbpntr[sub_i+1] = next_nz;
    	Cindx[next_nz] = next_value;
	    current->max_nz_per_row = max_nz_row_new; 
    	current->N_nonzeros     = total_nz; 
    	current->sub_matrix   = previous_matrix;
        current->getrow->columns_loc_glob = ML_GLOBAL_INDICES;

	    /* allocate space for new matrix */

    	if (i != 0) {
	      dtemp = ((double) N-i)/ ((double) i);
	      dtemp *= (double) total_nz;
    	  Nnz_estimate = (int)(1.1*dtemp);
    	  Nnz_estimate += Ncols;
    	}
    	else Nnz_estimate = Nnz_estimate*N + Ncols;

	    Cbpntr = (int    *) ML_allocate( (N-i+1)* sizeof(int) );
    	Cbindx = (int    *) ML_allocate( Nnz_estimate* sizeof(int) );
    	Cindx  = (int    *) ML_allocate( Nnz_estimate* sizeof(int));
    	Cvalues= (double *) ML_allocate( Nnz_estimate*
   					 NcolsPerBlock*NcolsPerBlock* sizeof(double));

    	if ((Cindx == NULL) || (Cbindx == NULL) || (Cvalues == NULL) ) {
    	  printf("Not enough space for matrix\n");
    	  exit(1);
    	}
        RowOffset = 0;
        RowOffsetBlocks = 0;
    	next_nz   = 0;
        next_value = 0;
    	Cbpntr[0]  = next_nz;

	    previous_matrix = current;
    	sub_i = 0;
      }

      /* store matrix row */

      allzeros = 1;
      for (k = 0; k < Ncols; k++) {
  	    Cbindx[next_nz] = accum_col[k];
	    Cindx[next_nz]  = next_value;
	    next_nz++;
	    allzeros = 0;
	    kkk=0;
	    for (iii = 0; iii < NcolsPerBlock; iii++) {
	      for (jjj = 0; jjj < NrowsPerBlock; jjj++) {
	        Cvalues[next_value++] = accum_val[k*NnzPerBlock+kkk];
	        kkk++;
	      }
	    }
      }
      /* if entire row is zero, store one entry to avoid empty row */

      if ((Ncols == 0) && (nearbyIndex != -1)) {
	    Cbindx[next_nz] = nearbyIndex;
	    Cindx[next_nz] = next_value;
	    next_nz++;
	    for (iii = 0; iii < NcolsPerBlock; iii++) {
	      for (jjj = 0; jjj < NrowsPerBlock; jjj++) {
	        Cvalues[next_value++] = accum_val[k*NnzPerBlock+kkk];
	        kkk++;
	      }
	    }
      }
      Cbpntr[sub_i+1] = next_nz;
      sub_i++;
      if (Ncols > max_nz_row_new) max_nz_row_new = Ncols;
      RowOffset += NrowsPerBlock;
      RowOffsetBlocks++;
    }
    oldstart = start;
    start = end;
  }
  Cindx[next_nz] = next_value;

  if (B_indx != NULL) ML_free(B_indx);
  if (Bcols != NULL) ML_free(Bcols);
  if (Bptr != NULL) ML_free(Bptr);
  if (accum_index != NULL) ML_free(accum_index);
  if (col_inds != NULL) ML_free(col_inds);
  ML_free(accum_col);
  ML_free(accum_val);
  ML_free(Aindx);
  ML_free(A_i_cols);

  /*test output*/
  /*for (i = 0; i <2; i++)
  {
    for (j = 0; j < 32; j++)
      printf("%lf\n", Cvalues[j+i*16]); fflush(stdout);
    printf("indx %d\n", Cindx[i*3]);fflush(stdout);
    for(j = 0; j < 4; j++)
    {
      printf("bindx %d\n", Cbindx[j+i*3]);fflush(stdout);
      printf("indx %d\n", Cindx[j+1+i*3]);fflush(stdout);
    }
  }*/


  /* create 'parent' object corresponding to the resulting matrix */

  total_nz += next_nz;
  Cvbr_mat = (struct ML_vbrdata *) ML_allocate(sizeof(struct ML_vbrdata));
  Cvbr_mat->bindx       = Cbindx;
  Cvbr_mat->val         = Cvalues;
  Cvbr_mat->bpntr       = Cbpntr;
  Cvbr_mat->indx        =  Cindx; 
  


  Ccpntr = (int *) ML_allocate(sizeof(int)*hash_used);
  Ccpntr[0] = 0;
  for (kkk = 1; kkk <= hash_used; kkk++) 
     Ccpntr[kkk] = Ccpntr[kkk-1] + NcolsPerBlock;
  Cvbr_mat->cpntr       = Ccpntr;
  Cvbr_mat->rpntr       = Crpntr;
  *Cmatrix = ML_Operator_Create(Amatrix->comm);
  ML_Operator_Set_1Levels((*Cmatrix), Bmatrix->from, Amatrix->to);

   if((*Cmatrix)->ML_id != ML_ID_OP ) {
      printf("ML_Operator_Set_ApplyFunc error : wrong object.\n");
      exit(-1);
   }
/* newly added : 8/17/00 */
   if ((*Cmatrix)->data != NULL && (*Cmatrix)->data_destroy != NULL )
   {
      (*Cmatrix)->data_destroy((*Cmatrix)->data);
      (*Cmatrix)->data = NULL;
   }
   (*Cmatrix)->invec_leng = save_ints[5];
   (*Cmatrix)->outvec_leng = save_ints[1];
   (*Cmatrix)->data = (void*)Cvbr_mat;
   (*Cmatrix)->getrow->func_ptr = VBR_getrows;
   (*Cmatrix)->getrow->columns_loc_glob = ML_GLOBAL_INDICES;
                                                                                                                 
   (*Cmatrix)->matvec->ML_id = ML_NONEMPTY;
   if(previous_matrix != NULL)
   {
     (*Cmatrix)->matvec->Nrows = RowOffset+previous_matrix->getrow->Nrows;
     (*Cmatrix)->getrow->Nrows = RowOffset+previous_matrix->getrow->Nrows;
     (*Cmatrix)->getrow->N_block_rows = RowOffsetBlocks+previous_matrix->getrow->N_block_rows;
   }
   else
   {
     (*Cmatrix)->matvec->Nrows = RowOffset;
     (*Cmatrix)->getrow->Nrows = RowOffset;
     (*Cmatrix)->getrow->N_block_rows = RowOffsetBlocks;
   }

  /*ML_Operator_Set_Getrow((*Cmatrix), RowOffset+NrowsPerBlock, 
			 az_vbrgetrow_wrapper);*/


  /* (*Cmatrix)->data_destroy = AZ_ML_FullClean; */

        
  Cindx[next_nz] = next_value;
  (*Cmatrix)->max_nz_per_row = max_nz_row_new; 
  (*Cmatrix)->N_nonzeros     = total_nz; 
  if(previous_matrix == NULL)
    printf("NULL\n");

  (*Cmatrix)->getrow->Nrows = N;
  (*Cmatrix)->sub_matrix     = previous_matrix;

  if (A_i_allocated-1 > Amatrix->max_nz_per_row) 
    Amatrix->max_nz_per_row = A_i_allocated;

  (*Cmatrix)->N_total_cols_est = hash_used;
  Bmatrix->N_total_cols_est = hash_used;


  if (Bmatrix->getrow->pre_comm != NULL) {
    ML_CommInfoOP_Clone(&((*Cmatrix)->getrow->pre_comm),
			Bmatrix->getrow->pre_comm);
  }

  Amatrix->getrow->Nrows = save_ints[0]; 
  Amatrix->outvec_leng   = save_ints[1];  
  Amatrix->invec_leng    = save_ints[2];   
  Bmatrix->getrow->Nrows = save_ints[3];
  Bmatrix->outvec_leng   = save_ints[4];  
  Bmatrix->invec_leng    = save_ints[5];   
  (*Cmatrix)->getrow->Nrows = save_ints[0]; 
  (*Cmatrix)->outvec_leng   = save_ints[1]; 
  (*Cmatrix)->invec_leng    = save_ints[5];   

} /* ML_blkmatmat_mult() */
#endif

/* ******************************************************************** */
void ML_matmat_mult(ML_Operator *Amatrix, ML_Operator *Bmatrix,
                    ML_Operator **Cmatrix)
/* -------------------------------------------------------------------- */
{
   int    i,k, jj, next_nz, Ncols, N, Nnz_estimate, sub_i, accum_size, row;
   int    *C_ptr, *Ccol, *A_i_cols, rowi_N, *accum_col, row2_N;
   double *Cval, *A_i_vals, multiplier, *accum_val, dtemp;
   ML_Operator *current, *previous_matrix, *next;
   struct ML_CSR_MSRdata *temp;
   int    max_nz_row_new = 0, total_nz = 0, index_length = 0;
   int    min_nz_row_new = 1e6;
   int    avg_nz_row_new = 0;
   double A_avg_nz_per_row, B_avg_nz_per_row, estimated_nz_per_row;
   int    A_i_allocated;
   int    flag;
   void   (*Agetrow)(ML_Operator *,int,int *,int *,int **,double **,int *,int);
   void   (*Bgetrow)(ML_Operator *,int,int *,int *,int **,double **,int *,int);
   double *Bvals;
   int    *Bptr, *Bcols;
   int *col_inds, B_total_Nnz, itemp, B_allocated, hash_val, *accum_index, lots_of_space;
   int rows_that_fit, rows_length, *rows, NBrows, end, start = 0;
   int subB_Nnz, Next_est, total_cols = 0;
   int tcols, hash_used, j, *tptr;
   int *acc_col_ptr, *Bcol_ptr; double *acc_val_ptr, *Bval_ptr;
   int allzeros;
   int hashTableIsPowerOfTwo = 0;
   int nearbyIndex;
   N  = Amatrix->getrow->Nrows;

   /**************************************************************************/
   /* figure out if a cheaper getrow function (either ML_get_matrow_CSR or   */
   /* ML_get_row_CSR_norow_map) can be used avoiding some function overhead. */
   /* ---------------------------------------------------------------------- */

   Bgetrow = ML_get_matrix_row;
   flag    = 1;
   next    = Bmatrix;
   while (next != NULL) {
      if (next->getrow->func_ptr != CSR_getrow) flag = 0;
      next = next->sub_matrix;
   }
   if (flag == 1) { 
      Bgetrow = ML_get_matrow_CSR;

      next = Bmatrix;
      while (next != NULL) {
         if (next->getrow->row_map!= NULL) flag = 0;
         next = next->sub_matrix;
      }
      if (flag == 1) Bgetrow = ML_get_row_CSR_norow_map;
   }

   Agetrow = ML_get_matrix_row;
   flag    = 1;
   next    = Amatrix;
   while (next != NULL) {
      if (next->getrow->func_ptr != CSR_getrow) flag = 0;
      next = next->sub_matrix;
   }
   if (flag == 1) Agetrow = ML_get_matrow_CSR;

   /********************************************************/
   /* put here to alleviate a hidden error - C. Tong       */
   if ( Amatrix->max_nz_per_row > Amatrix->getrow->Nrows )
      Amatrix->max_nz_per_row = Amatrix->getrow->Nrows;
   if ( Bmatrix->max_nz_per_row > Bmatrix->getrow->Nrows )
      Bmatrix->max_nz_per_row = Bmatrix->getrow->Nrows;
   /********************************************************/

   /*********************************************************/
   /* Count the number of external variables. Unfortunately */
   /* ML_exch_row() does not properly put the correct       */
   /* communication information ... and so this is just an  */
   /* estimate of how many externals there are.             */
   /* ------------------------------------------------------*/

   current = Bmatrix;
   Next_est = 0;
   while(current != NULL) {
      if (current->getrow->pre_comm != NULL) {
    if (current->getrow->pre_comm->total_rcv_length <= 0) {
#ifdef charles
      if (Amatrix->comm->ML_mypid == 0)
         printf("*** ML_matmat_mult: %s x %s\n",Amatrix->label,Bmatrix->label);
      printf("*** %d: recomputing rcv length (%u %u)\n",Bmatrix->comm->ML_mypid,Bmatrix,current); fflush(stdout);
#endif
      ML_CommInfoOP_Compute_TotalRcvLength(current->getrow->pre_comm);
    }

         Next_est += current->getrow->pre_comm->total_rcv_length;
#ifdef charles
      printf("*** %d: Nghost = %d  %d (%d x %d)\n",Bmatrix->comm->ML_mypid,Next_est,current->getrow->pre_comm->total_rcv_length,Bmatrix->outvec_leng,Bmatrix->invec_leng); fflush(stdout);
#endif
      }
#ifdef charles
      else {printf("%d: pre_comm is null?\n",Bmatrix->comm->ML_mypid); fflush(stdout);}
#endif
      current = current->sub_matrix;
   }
#ifdef charles
   printf("\n\n%d: Next_est = %d\n\n",Bmatrix->comm->ML_mypid,Next_est);
#endif
   /* Estimate the total number of columns in Bmatrix. Actually,  */
   /* we need a crude estimate. We will use twice this number for */
   /* hashing the columns. Once we are finished hashing columns   */
   /* we will count the number of columns and revise the estimate */
   /* for the number of externals and increase the size of the    */
   /* accumulator if needed.                                      */
   /* ----------------------------------------------------------- */

   accum_size = Bmatrix->invec_leng + 75;
   accum_size += 3*Next_est;
   if (Bmatrix->N_total_cols_est + 75 > accum_size ) {
     accum_size = Bmatrix->N_total_cols_est + 75;
   }

           /* the '3' above is a total kludge. the problem is that we   */
           /* need to know the number of neighbors in the external rows */
           /* and this is annoying to compute */

   /**************************************************************************/
   /* allocate space to hold a single row of Amatrix and the accumulator     */
   /* which contains the current row in the resulting matrix.                */
   /*------------------------------------------------------------------------*/

   A_i_allocated = Amatrix->max_nz_per_row + 1;
   A_i_cols  = (int    *) ML_allocate(A_i_allocated * sizeof(int) );
   A_i_vals  = (double *) ML_allocate(A_i_allocated * sizeof(double));
   accum_col = (int    *) ML_allocate( accum_size * sizeof(int) );
   accum_val = (double *) ML_allocate( accum_size * sizeof(double) );
   if ( (A_i_vals == NULL) || (accum_val == NULL)) {
      printf("Not enough space in ML_matmatmult().\n");
      printf("trying to allocate %d %d elements \n",A_i_allocated,accum_size);
      printf("Left  matrix has %d rows \n", Amatrix->getrow->Nrows);
      printf("Right matrix has %d rows \n", Bmatrix->getrow->Nrows);
      printf("Left  matrix has %d nz per row\n", Amatrix->max_nz_per_row);
      printf("Right matrix has %d nz per row\n", Bmatrix->max_nz_per_row);
      exit(1);
   }

   /**************************************************************************/
   /* Make conservative estimates of the size needed to hold the resulting   */
   /* matrix. NOTE: These arrays can be increased later in the computation   */
   /*                                                                        */
   /* We try to use the actual data from the rowptrs if we can get at 'em,   */
   /* otherwise we probe a few rows to generate a guess.                     */
   /*------------------------------------------------------------------------*/
   ML_estimate_avg_nz_per_row(Amatrix,&A_avg_nz_per_row);
   if (!A_avg_nz_per_row) A_avg_nz_per_row=100;
   ML_estimate_avg_nz_per_row(Bmatrix,&B_avg_nz_per_row);
   if (!B_avg_nz_per_row) B_avg_nz_per_row=100;

   estimated_nz_per_row = sqrt(A_avg_nz_per_row) + sqrt(B_avg_nz_per_row) - 1.;
   estimated_nz_per_row *= estimated_nz_per_row;
   Nnz_estimate = (int) (((double) Amatrix->getrow->Nrows) *
                           estimated_nz_per_row*.75) + 100;
   if (Nnz_estimate <= Bmatrix->max_nz_per_row)
      Nnz_estimate = Bmatrix->max_nz_per_row + 1;

   /**************************************************************************/
   /* Allocate space for the matrix. Since 'Nnz_estimate' is just an estimate*/
   /* of the space needed, we will reduce 'Nnz_estimate' if we are not       */
   /* successful allocating space.                                           */
   /*------------------------------------------------------------------------*/

   C_ptr     = (int    *) ML_allocate((N+1)* sizeof(int) );
   Cval = NULL; Ccol = NULL;
   while ( ((Cval == NULL)||(Ccol==NULL)) && (Nnz_estimate > Bmatrix->max_nz_per_row) ) {
      if (Ccol != NULL) ML_free(Ccol);
      Ccol  = (int    *) ML_allocate( Nnz_estimate* sizeof(int) );
      Cval  = (double *) ML_allocate( Nnz_estimate* sizeof(double));
      if (Cval == NULL) Nnz_estimate = Nnz_estimate/2;
   }
   if ( ((Cval == NULL)||(Ccol==NULL)) && (N != 0)) {
      printf("Not enough space for new matrix in ML_matmatmult().\n");
      printf("trying to allocate %d elements \n",Nnz_estimate);
      printf("Left  matrix has %d rows \n", Amatrix->getrow->Nrows);
      printf("Right matrix has %d rows \n", Bmatrix->getrow->Nrows);
      printf("Left  matrix has %d nz per row\n", Amatrix->max_nz_per_row);
      printf("Right matrix has %d nz per row\n", Bmatrix->max_nz_per_row);
      exit(1);
   }

   /**************************************************************************/
   /* Make a copy of Bmatrix and use this copy to speed up the computation.  */
   /*   1) First count the total number of nonzeros                          */
   /*   2) allocate space for Bmatrix copy.                                  */
   /*      Note: accum_index is used for hashing. We make it twice the size  */
   /*            of the accumulator for good performance.                    */
   /*   3) grab the matrix. Put the nonzero values in Bcols & Bvals.         */
   /*   4) Using a hash table, renumber the column indices so that they are  */
   /*      in [0:index_length-1]. Store the original column indices in       */
   /*      col_inds[].                                                       */
   /*------------------------------------------------------------------------*/

#ifdef TO_BE_CHECKED
   if (Bmatrix->N_nonzeros <= 0) {
      B_total_Nnz = 0;
      for (i = 0; i < Bmatrix->getrow->Nrows; i++ ) {
         Bgetrow(Bmatrix,1,&i,&accum_size, &accum_col, &accum_val, &row2_N, 0);
         B_total_Nnz += row2_N;
      }
   }
   else B_total_Nnz = Bmatrix->N_nonzeros;
#endif
   B_total_Nnz = 0;
   for (i = 0; i < Bmatrix->getrow->Nrows; i++ ) {
     Bgetrow(Bmatrix,1,&i,&accum_size, &accum_col, &accum_val, &row2_N, 0);
     B_total_Nnz += row2_N;
   }
   B_total_Nnz++; /* avoid possible division by zero below */

   /* ******************************************************************
      Try to make the hash table size a power of two.  Doing so makes
      hash table lookup more efficient by using the bitwise AND (&) operator
      and avoiding the modulus (%) function.

      The new lookup is ML_fast_hash, the original is ML_hash_it.  All hashing
      is located in ml/src/Utils/ml_utils.[ch].

      First check for and handle potential overflow:

            Case 1:  If accum_size is larger than 2^29, but not by much,
                     use a hash table size of 2^30.
            Case 2:  If accum_size is larger than 2^29 but less than 2^30,
                     use a hash table size of 2*accum_size.
            Case 3:  If accum_size is larger than 2^30 but not by much,
                     use a hash table size of 2^30-1 (INT_MAX).
            Case 4:  If accum_size is much larger than 2^30,
                     issue an error and abort.

      Otherwise, choose the hash table to be 2^k, where 2^k > 2*accum_size,
      k is a positive integer, and k < 31.

      Note:
      This potentially wastes some space if accum_size is close to but not
      exactly a power of 2.

      Some useful constants:

        2^31-1 = 2147483647 = INT_MAX
        2^30   = 1073741824
        2^29   =  536870912

   ****************************************************************** */
 
   if ( accum_size > 536870912 ) {
     if ( 1.8*accum_size < 1073741824) {
       index_length = 1073741824;
       hashTableIsPowerOfTwo = 1;
     } else if (accum_size < 1073741824) {
       index_length = 2*accum_size;
       hashTableIsPowerOfTwo = 0;
     } else if (1.8*accum_size < INT_MAX) {
       index_length = INT_MAX;
       hashTableIsPowerOfTwo = 0;
     } else {
       printf("Cannot create hash table for accumulator of size %d\n",
              accum_size);
       exit(1);
     }
   } else {
     index_length = 2;
     while (index_length < 2*accum_size) {
       index_length *= 2;
     }
     hashTableIsPowerOfTwo = 1;
   }

   accum_index  = (int *) ML_allocate(index_length*sizeof(int));
   col_inds     = (int *) ML_allocate(index_length*sizeof(int));

   Bptr      = (int    *) ML_allocate( (Bmatrix->getrow->Nrows+1)*sizeof(int));
   B_allocated = B_total_Nnz * 2;
   lots_of_space = 0;

   Bcols     = NULL; Bvals = NULL;
   while (Bvals == NULL) {
      lots_of_space++;
      if (Bcols != NULL) ML_free(Bcols);
      Bcols     = (int    *) ML_allocate( B_allocated * sizeof(int));
      Bvals     = (double *) ML_allocate( B_allocated * sizeof(double));
      B_allocated /= 2;
   }
   ML_free(Bvals); ML_free(Bcols);

   Bcols     = (int    *) ML_allocate( B_allocated * sizeof(int));
   Bvals     = (double *) ML_allocate( B_allocated * sizeof(double));
   if (lots_of_space != 1) lots_of_space = 0;

   if (Bmatrix->outvec_leng != 0) {
     dtemp = ((double) B_allocated)/ ( (double) B_total_Nnz);
     dtemp *= .9;
     rows_that_fit = (int) (( (double) Bmatrix->outvec_leng) * dtemp);
     if (rows_that_fit == 0) rows_that_fit++;
   }
   else rows_that_fit = 1;


   *Cmatrix        = NULL;
   previous_matrix = NULL;
   sub_i           = 0;
   next_nz         = 0;
   C_ptr[0]        = next_nz;

   end = N;
   start = 0;

   hash_used = 0;
   while (start < N) {

   itemp = 0;

   Bptr[0] = 0;
   if (lots_of_space) {
     for (i = 0; i < Bmatrix->getrow->Nrows; i++ ) {
       Bgetrow(Bmatrix, 1, &i, &B_allocated, &Bcols, &Bvals, &row2_N, itemp);
       itemp += row2_N;
       Bptr[i+1] = itemp;
     }
     subB_Nnz = itemp;
   }

   if (!lots_of_space) {
      NBrows = 0;
      rows_length = rows_that_fit*2 + 3;
      rows = (int *) ML_allocate(rows_length*sizeof(int));
      for (i=0; i < rows_length; i++) rows[i] = -1;
      ML_determine_Brows(start, &end, Amatrix,
               &rows, &rows_length, &NBrows,
             &rows_that_fit, Agetrow);

      jj = 0;
      itemp = 0;
      ML_az_sort(rows, NBrows, NULL, NULL);
      for (i = 0; i < Bmatrix->getrow->Nrows; i++ ) {
    if (rows[jj] == i) {
      jj++;
      k = B_allocated;
      Bgetrow(Bmatrix,1, &i, &B_allocated, &Bcols, &Bvals, &row2_N, itemp);
      itemp += row2_N;
    }
    Bptr[i+1] = itemp;
      }
      subB_Nnz = itemp;
   }

   for (i = 0; i < index_length; i++) accum_index[i] = -1;
   for (i = 0; i < index_length; i++) col_inds[i] = -1;

   i = 0;
   if (Bmatrix->getrow->pre_comm != NULL) 
      i = Bmatrix->getrow->pre_comm->total_rcv_length;
#ifdef charles
   printf("%d: before hashing, hash table size = %d invec and Nghost = (%d %d)\n",

      Bmatrix->comm->ML_mypid, index_length, Bmatrix->invec_leng, Next_est);
   fflush(stdout);
#endif

   tcols = 0;
   hash_used = 0;

   /* Record a "nearby" index for use as a column index in any empty rows. */
   /* This avoids creating a communication pattern in which a single       */
   /* processor must communicate with all others in                        */
   /* ML_CommInfoOP_GenUsingGIDExternals().                                */

   if (subB_Nnz > 0) nearbyIndex = Bcols[0];
   else              nearbyIndex = -1;

   if (hashTableIsPowerOfTwo)
   {
     int ilm1 = index_length-1;
     for (i = 0; i < subB_Nnz; i++) {
       if (hash_used >= ((int) (.5 * index_length)) ) {
#ifdef charles
         printf("%d: running out of hashing space: row = %d, tcols=%d"
                " hash_length=%d, hash_used = %d\n",
                Bmatrix->comm->ML_mypid,i,tcols,index_length,hash_used);
         fflush(stdout);
#endif
         ML_free(accum_index);
         if (index_length >= 1073741824) /* 2^30 */
           pr_error("Exceeded largest possible hash table size\n");
         accum_index = (int *) ML_allocate(sizeof(int)*2*index_length);
         tptr = (int *) ML_allocate(sizeof(int)*2*index_length);
         if (tptr == NULL) pr_error("ML_matmat_mult: out of tptr space\n");
         for (j = 0; j < 2*index_length; j++) accum_index[j] = -1;
         for (j = 0; j < 2*index_length; j++) tptr[j] = -1;
         for (j = 0; j < i; j++)
           Bcols[j] = col_inds[Bcols[j]];
         ML_free(col_inds);  col_inds = tptr;
         tcols = 0;
         hash_used = 0;
         index_length *= 2;
         ilm1 = index_length-1;
         for (j = 0; j < i; j++) {
           ML_fast_hash(Bcols[j], col_inds,ilm1,&hash_used, &hash_val);
           if (col_inds[hash_val] == -1) tcols++;
           col_inds[hash_val] = Bcols[j];
           Bcols[j] = hash_val;
         }
       }
       ML_fast_hash(Bcols[i], col_inds, ilm1, &hash_used, &hash_val);
       if (col_inds[hash_val] == -1) tcols++;
       col_inds[hash_val] = Bcols[i];
       Bcols[i] = hash_val;
     } /*for (i = 0; i < subB_Nnz; i++)*/

   } else {

     for (i = 0; i < subB_Nnz; i++) {
       if (hash_used >= ((int) (.75 * index_length)) ) {
#ifdef charles
         printf("%d: running out of hashing space: row = %d, tcols=%d"
                " hash_length=%d, hash_used = %d\n",
                Bmatrix->comm->ML_mypid,i,tcols,index_length,hash_used);
         fflush(stdout);
#endif
         ML_free(accum_index);
         if (index_length >= 1073741824) /* 2^30 */
           pr_error("Exceeded largest possible hash table size\n");
         accum_index = (int *) ML_allocate(sizeof(int)*2*index_length);
         tptr = (int *) ML_allocate(sizeof(int)*2*index_length);
         if (tptr == NULL) pr_error("ML_matmat_mult: out of tptr space\n");
         for (j = 0; j < 2*index_length; j++) accum_index[j] = -1;
         for (j = 0; j < 2*index_length; j++) tptr[j] = -1;
         for (j = 0; j < i; j++)
           Bcols[j] = col_inds[Bcols[j]];
         ML_free(col_inds);  col_inds = tptr;
         tcols = 0;
         hash_used = 0;
         index_length *= 2;
         for (j = 0; j < i; j++) {
           ML_hash_it(Bcols[j], col_inds, index_length, &hash_used, &hash_val);
           if (col_inds[hash_val] == -1) tcols++;
           col_inds[hash_val] = Bcols[j];
           Bcols[j] = hash_val;
         }
       }
         
       ML_hash_it(Bcols[i], col_inds, index_length, &hash_used, &hash_val);
       if (col_inds[hash_val] == -1) tcols++;
       col_inds[hash_val] = Bcols[i];
       Bcols[i] = hash_val;
     } /*for (i = 0; i < subB_Nnz; i++)*/

   } /*if (hashTableIsPowerOfTwo)*/

#ifdef charles
   printf("%d: after hashing %d %d\n",Bmatrix->comm->ML_mypid,tcols,hash_used);
   fflush(stdout);
#endif

   /* Count how many columns there are, record the number of   */
   /* columns and change the size of the accumulator if needed */

   jj = 0;
   for (i = 0; i < index_length; i++) { 
      if (col_inds[i] == -1) jj++;
   }
   total_cols += (index_length - jj);
   if (accum_size < total_cols) {
      ML_free(accum_col);
      ML_free(accum_val);
      accum_size = total_cols + 1;
      accum_col = (int *) ML_allocate(accum_size*sizeof(int));
      accum_val = (double *) ML_allocate(accum_size*sizeof(double));
      if (accum_val == NULL) 
    pr_error("ML_matmat_mult: no room for accumulator\n");
   }

   /**************************************************************************/
   /* Perform the matrix-matrix multiply operation by computing one new row  */
   /* at a time.                                                             */
   /*------------------------------------------------------------------------*/

   for (i = start; i < end ; i++) {
      Ncols = 0;
      acc_col_ptr = accum_col;
      acc_val_ptr = accum_val;
      Agetrow(Amatrix,1, &i, &A_i_allocated, &A_i_cols, &A_i_vals, &rowi_N, 0);

      /***********************************************************************/
      /* Take each column entry in the ith row of A and find corresponding   */
      /* row in B. Muliply each B row and sum it into the accumulator.       */
      /*---------------------------------------------------------------------*/

      if (rowi_N != 0) {
        multiplier = A_i_vals[0];
        if ( multiplier != 0.0 ) {
          jj = Bptr[*A_i_cols];
          Bval_ptr = &(Bvals[jj]);
          Bcol_ptr = &(Bcols[jj]);
          while (jj++ < Bptr[*A_i_cols+1]) {
            *acc_col_ptr =  *Bcol_ptr; acc_col_ptr++;
            *acc_val_ptr = multiplier*(*Bval_ptr); acc_val_ptr++; Bval_ptr++;
            accum_index[*Bcol_ptr] = Ncols++; Bcol_ptr++;
          }
        }
      }
      for (k = 1; k < rowi_N; k++) {
        multiplier = A_i_vals[k];
        if ( multiplier != 0.0 ) {
          jj = Bptr[A_i_cols[k]];
          Bval_ptr = &(Bvals[jj]);
          Bcol_ptr = &(Bcols[jj]);
          while (jj++ < Bptr[A_i_cols[k]+1]) {
            if (accum_index[*Bcol_ptr] < 0) {
              *acc_col_ptr = *Bcol_ptr; acc_col_ptr++;
              *acc_val_ptr = multiplier*(*Bval_ptr); acc_val_ptr++; Bval_ptr++;
              accum_index[*Bcol_ptr] = Ncols++; Bcol_ptr++;
            }
            else {
              accum_val[accum_index[*Bcol_ptr]] += multiplier*(*Bval_ptr);
              Bcol_ptr++; Bval_ptr++;
            }
          }
        }
      }
      /***********************************************************************/
      /* Convert back to the original column indices.                        */
      /*---------------------------------------------------------------------*/

      acc_col_ptr = accum_col;
      for (jj = 0; jj < Ncols; jj++ ) {
         accum_index[*acc_col_ptr] = -1;
         *acc_col_ptr = col_inds[*acc_col_ptr];
         acc_col_ptr++;
      }

      /* empty row. Let's put a zero in a nearby column. */

      if ((Ncols == 0) && (nearbyIndex != -1)) {
         accum_col[Ncols] = nearbyIndex;
         accum_val[Ncols++] = 0.0;
      }

      /* check if we have enough space to store the new matrix row */
      /* If not, place the current matrix into a submatrix and     */
      /* allocate new vectors (hopefully large enough to hold the  */
      /* rest of the matrix) and a 'parent' matrix to hold the     */
      /* remaining rows.                                           */

      if (next_nz+Ncols > Nnz_estimate) {

         /* create sub_matrix object */

         total_nz += next_nz;

         temp = (struct ML_CSR_MSRdata *) 
                ML_allocate(sizeof(struct ML_CSR_MSRdata));
         temp->columns         = Ccol;
         temp->values          = Cval;
         temp->rowptr          = C_ptr;
         current = ML_Operator_Create(Amatrix->comm);
         ML_Operator_Set_1Levels(current, Bmatrix->from, Amatrix->to);
         ML_Operator_Set_ApplyFuncData(current, Bmatrix->invec_leng, 
                       Amatrix->outvec_leng, temp,
                       i,NULL,0);

         ML_Operator_Set_Getrow(current, i, CSR_getrow);

         current->max_nz_per_row = max_nz_row_new; 
         current->min_nz_per_row = min_nz_row_new; 
         current->N_nonzeros     = total_nz; 
         current->sub_matrix   = previous_matrix;

         /* allocate space for new matrix */

         if (i != 0) {
            dtemp = ((double) N-i)/ ((double) i);
            dtemp *= (double) total_nz;
            Nnz_estimate = (int)(1.1*dtemp);
            Nnz_estimate += Ncols;
         }
         else Nnz_estimate = Nnz_estimate*N + Ncols;

         C_ptr = (int    *) ML_allocate( (N-i+1)* sizeof(int) );
         Ccol  = (int    *) ML_allocate( Nnz_estimate* sizeof(int) );
         Cval  = (double *) ML_allocate( Nnz_estimate* sizeof(double));
         if ((Cval == NULL) || (Ccol == NULL)) {
            printf("Not enough space for matrix\n");
            exit(1);
         }
         next_nz   = 0;
         C_ptr[0]  = next_nz;

         previous_matrix = current;
         sub_i = 0;
      }

      /* store matrix row */

      if (ML_Use_LowMemory() != ML_TRUE)
      {
        memcpy(&(Ccol[next_nz]),accum_col, sizeof(int)*Ncols);
        memcpy(&(Cval[next_nz]),accum_val, sizeof(double)*Ncols);
        next_nz += Ncols;
      }
      else
      {
        /* above code might be a bit faster??? */
        allzeros = 1;
        for (k = 0; k < Ncols; k++) {
          /* This 'if' might break some applications somewhere */
          /* but I can't remember who and where or why?        */
          /* For now, I want to reduce memory in alegra so I am*/
          /* putting it in. If we later want to take this out  */
          /* we should use the memcpy code above.              */
          if (accum_val[k] != 0.0) {
            allzeros = 0;
            Ccol[next_nz] = accum_col[k];
            Cval[next_nz++] = accum_val[k];
          }
        }
        /* if entire row is zero, store one entry to avoid empty row */
        if (allzeros) {
          Ccol[next_nz] = accum_col[0];
          Cval[next_nz++] = accum_val[0];
        }
      }
      C_ptr[sub_i+1] = next_nz;
      sub_i++;
      if (Ncols > max_nz_row_new) max_nz_row_new = Ncols;
      if (Ncols > 0 && Ncols < min_nz_row_new) min_nz_row_new = Ncols;
   } /*for (i = start; i < end ; i++) */
   start = end;
   } /*while (start < N) */
   if (Bvals != NULL) ML_free(Bvals);
   if (Bcols != NULL) ML_free(Bcols);
   if (Bptr != NULL) ML_free(Bptr);
   if (accum_index != NULL) ML_free(accum_index);
   if (col_inds != NULL) ML_free(col_inds);
   ML_free(accum_col);
   ML_free(accum_val);
   ML_free(A_i_vals);
   ML_free(A_i_cols);

   /* create 'parent' object corresponding to the resulting matrix */

   total_nz += next_nz;
   temp = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata));
   temp->columns          = Ccol;
   temp->values           = Cval;
   temp->rowptr           = C_ptr;
   *Cmatrix = ML_Operator_Create(Amatrix->comm);
   ML_Operator_Set_1Levels(*Cmatrix, Bmatrix->from, Amatrix->to);
   ML_Operator_Set_ApplyFuncData(*Cmatrix,Bmatrix->invec_leng, 
                    Amatrix->outvec_leng,temp,N,NULL,0);
   ML_Operator_Set_Getrow(*Cmatrix, N, CSR_getrow);
   (*Cmatrix)->getrow->Nrows = N;
   (*Cmatrix)->max_nz_per_row = max_nz_row_new;
   (*Cmatrix)->min_nz_per_row = min_nz_row_new;
   (*Cmatrix)->N_nonzeros     = total_nz;
   (*Cmatrix)->sub_matrix     = previous_matrix;
   if (A_i_allocated-1 > Amatrix->max_nz_per_row) 
      Amatrix->max_nz_per_row = A_i_allocated;

   (*Cmatrix)->N_total_cols_est = hash_used;
   Bmatrix->N_total_cols_est = hash_used;


   if (Bmatrix->getrow->pre_comm != NULL) {
     ML_CommInfoOP_Clone(&((*Cmatrix)->getrow->pre_comm),
             Bmatrix->getrow->pre_comm);
   }
} /* ML_matmat_mult() */

/*********************************************************************************/
/*Recives an ML_Operator in_matrix in a sparse data structure(csr, msr or vbr and returns  */
/*the same ML_Operator stored as a vbr matrix. 

  When submatrix != 0

  On input row_block_size and col_block_size are positive integers if either the column
  or row block sizes are fixed.  They are zero if there is not a constant block size.

  When row_block_size is zero then rpntr contains the row blocking information.
  When col_block_size is zero then cpntr contains the column blocking information.

  If either row_block_size or col_block_size is negative there is a problem and the function should
  quit with an error.

  This function assumes no empty columns.

  When submatrix = 1

  On input the only parameter that matters is the col_block_size which must be greater
  than 0.  All other parameters are ignored as there are not important to the final result.
  However rpntr and cpntr are allocated and used in the new VBR structure.

  The matrix passed in should not be a submatrix and its sub-matrix will be untouched.
                                  */
/*********************************************************************************/

/*Where data is being exchanged things can become more efficient with an isend though this has risks.  Only works on matrices with one submatrix*/
void ML_convert2vbr(ML_Operator *in_matrix, int row_block_size, int rpntr[], int col_block_size, int cpntr[], int submatrix)
{
   int i, j, k, kk, ii, ll, m;
   int jj = 0;
   int blockrows = 0;
   int blockcolumns = 0;
   int row_length; /*the length of the current point row*/
   int blockend; /*the end point of a block in terms of size*/
   int blockstart; /*the start point of a block in terms of size*/
   int blockfound;
   int A_i_allocated;
   int spaceneeded;
   int val_size;
   int temp_alc; 
   int spacetight = 0;
   int *A_i_cols; /*columns for the single point row of values*/
   double *accum_val; /*storage of a single point row of values*/
   double *vals;  /*block row storage of matrix*/
   int *bindx, *indx, *bpntr; /*quick links to needed structures
                                                also makes the code look prettier*/
   int iplus1 = 0;
   int iminus1 = 0;
   double x;
   double *temp_double;
   struct ML_vbrdata *out_data;
   int nnz = 0;
   int Nghost_nodes, total_cols = 0;
   int **all_rpntr; /*stores all rpntr's on recieve*/
   USR_REQ *recv_requests;
   int rows_sent;
   int *send_buffer;
   ML_Operator *cur_matrix;
   struct ML_CSR_MSRdata *cur_data;
   int src, mid;
   int *recv_lens;
   int neighbors;
   int mats = 0;
   ML_CommInfoOP *pre_comm;
   USR_COMM USR_comm;     

   if(submatrix == 1)
   {
     /*We need to free rpntr if it is allocated since we are using it later*/
     if(rpntr != NULL)
       ML_free(rpntr);
     /*we only need to worry about block column informtion as this is all we will use*/  
     if (col_block_size < 0)
     {
       pr_error("In function convert2vbr col_block_size is negative which is undefined.  Please see function header and pass in the appropriate value.\n");
     }
     if (col_block_size == 0)
     {
         pr_error("In function convert2vbr col_block_size is 0 which is undefined for submatrices.  Please see function header and pass in the appropriate value\n");
     }
     else
     {
       cur_matrix = in_matrix;
       /*how many submatrices we have*/
       while(cur_matrix->sub_matrix != NULL)
       {
         mats++;
         cur_matrix = cur_matrix->sub_matrix;
       }
       USR_comm = cur_matrix->comm->USR_comm;
       pre_comm = in_matrix->getrow->pre_comm;
       neighbors = pre_comm->N_neighbors;
       out_data = (struct ML_vbrdata *)cur_matrix->data;
       all_rpntr = (int**)ML_allocate(neighbors*sizeof (int*));
       recv_lens = (int*)ML_allocate(neighbors*sizeof(int));
       if((recv_lens == NULL) || (all_rpntr == NULL))
       {
          pr_error("Not enough space to allocate %d int pointers and ints in convert2vbr.\n", neighbors);
       }
       recv_requests = (USR_REQ*)ML_allocate(neighbors*sizeof(USR_REQ));
       if(recv_requests == NULL)
       {
          pr_error("Not enough space to allocate %d requests in convert2vbr.\n", neighbors);
       }

       for(i = 0; i < neighbors; i++)
       {
         if(pre_comm->neighbors[i].N_rcv > 0)
         {
           /*This will get how much data is coming*/
           ML_Comm_Irecv(&recv_lens[i], sizeof(int), &(pre_comm->neighbors[i].ML_id), &(pre_comm->neighbors[i].ML_id), USR_comm, &recv_requests[i]);
         }
       }
       /*this could be done later as the max actually sent but probably is not worth the few bytes savings*/
       send_buffer = (int*)ML_allocate((cur_matrix->getrow->N_block_rows)*sizeof(int));
       if(send_buffer == NULL)
       {
          pr_error("Not enough space to allocate %d ints in convert2vbr.\n", cur_matrix->getrow->N_block_rows);
       }

       for(i = 0; i < neighbors; i++)
       {
         if(pre_comm->neighbors[i].N_send > 0)
         {
           k = 0;
           kk = 0;
           rows_sent = 0;
           /*figure out what needs to be sent*/
           while(k < pre_comm->neighbors[i].N_send)
           {
             kk = 0;
             rows_sent++;
             while(out_data->rpntr[kk] <= pre_comm->neighbors[i].send_list[k])
             {
               kk++;
             }
             k += out_data->rpntr[kk] - out_data->rpntr[kk - 1];
           }
           k =ML_Comm_Send(&rows_sent, sizeof(int), pre_comm->neighbors[i].ML_id, cur_matrix->comm->ML_mypid, USR_comm);
         }
       }
     }
     for(i = 0; i < neighbors; i++)
     {
       if(pre_comm->neighbors[i].N_rcv > 0)
       {
         ML_Comm_CheapWait(send_buffer, j, &src, &mid, USR_comm, &recv_requests[i]);
         all_rpntr[i] = (int*)ML_allocate((recv_lens[i])*sizeof(int));
         if(send_buffer == NULL)
         {
           pr_error("Not enough space to allocate %d ints in convert2vbr.\n", recv_lens[i]);
         }
         ML_Comm_Irecv(all_rpntr[i], recv_lens[i]*sizeof(int), &(pre_comm->neighbors[i].ML_id), &(pre_comm->neighbors[i].ML_id), USR_comm, &recv_requests[i]);
       }
     }
     for(i = 0; i < neighbors; i++)
     {
       j = 0;
       kk = 0;
       k = 0;
       rows_sent = 0;
       if(pre_comm->neighbors[i].N_send > 0)
       {
         while(j < pre_comm->neighbors[i].N_send)
         {
           kk = 0;
           rows_sent++;
           while(out_data->rpntr[kk] <= pre_comm->neighbors[i].send_list[j])
           {
             kk++;
           }
             send_buffer[rows_sent-1] = out_data->rpntr[kk] - out_data->rpntr[kk - 1];
             j += out_data->rpntr[kk] - out_data->rpntr[kk - 1];
         }
         k++;
         k =ML_Comm_Send(send_buffer, (rows_sent)*sizeof(int), pre_comm->neighbors[i].ML_id, cur_matrix->comm->ML_mypid, USR_comm);
       }
     }
      for(i = 0; i < neighbors; i++)
     {
       if(pre_comm->neighbors[i].N_rcv > 0)
       {
         ML_Comm_CheapWait(send_buffer, j, &src, &mid, USR_comm, &recv_requests[i]);
       }
     }
     ML_free(recv_requests);
     ML_free(send_buffer);
     /*loop over submatricesi this will only happen once but was left in from when functionality for more than one sub_matrix was attempted but there should never be more than one submatrix*/
     m = mats;
     while(m > 0)
     {
       m--;
       cur_matrix = in_matrix;
       for(j = 0; j < m; j++)
         cur_matrix = cur_matrix->sub_matrix;
       cur_data = (struct ML_CSR_MSRdata *)cur_matrix->data;
       i = 0; /*which blockrow*/
       k = 0; /*which processor*/
       cpntr = (int*)ML_allocate((((cur_matrix->invec_leng+cur_matrix->getrow->pre_comm->total_rcv_length)/col_block_size)+1)*sizeof(int));
       rpntr = (int*)ML_allocate((cur_matrix->outvec_leng- cur_matrix->sub_matrix->outvec_leng+1)*sizeof(int));
       bpntr = (int*)ML_allocate((cur_matrix->outvec_leng- cur_matrix->sub_matrix->outvec_leng+1)*sizeof(int));
       jj = ((cur_matrix->invec_leng+cur_matrix->getrow->pre_comm->total_rcv_length)/col_block_size)*cur_matrix->outvec_leng;
       if (jj > cur_data->rowptr[cur_matrix->outvec_leng - cur_matrix->sub_matrix->outvec_leng])
         jj = cur_data->rowptr[cur_matrix->outvec_leng - cur_matrix->sub_matrix->outvec_leng];
       bindx = (int*)ML_allocate(jj*sizeof(int));
       indx = (int*)ML_allocate((jj+1)*sizeof(int));
       vals = (double*)ML_allocate(cur_data->rowptr[cur_matrix->outvec_leng - cur_matrix->sub_matrix->outvec_leng]*sizeof(double));
       if(vals == NULL)
       {
          pr_error("Not enough space to allocate %d doubles and some integers in convert2vbr.\n", cur_data->rowptr[cur_matrix->outvec_leng - cur_matrix->sub_matrix->outvec_leng]);
       }
       cpntr[0] = 0;
       /*set up column pointer for submatrix*/
       for(kk = 1; kk <= (cur_matrix->invec_leng+cur_matrix->getrow->pre_comm->total_rcv_length)/col_block_size; kk++)
       {
         cpntr[kk] = cpntr[kk-1] + col_block_size;
       }
       rpntr[0] = bpntr[0] = indx[0] = jj = 0;
       /*loop a processor worth of matrix rows*/
       while(rpntr[i] < cur_matrix->getrow->Nrows - cur_matrix->sub_matrix->getrow->Nrows)
       {
         /*loop over one processors rows*/
         for(j = 0; j < recv_lens[k]; j++)
         {
           i++; /*blockrow*/
           iminus1 = i - 1;
           rpntr[i] = rpntr[iminus1] + all_rpntr[k][j]; 
           bpntr[i] = bpntr[iminus1] + (cur_data->rowptr[rpntr[iminus1]+1]-cur_data->rowptr[rpntr[iminus1]])/col_block_size;
           ll = indx[bpntr[iminus1]];
           /*loop over one blockrow*/
           for(kk = bpntr[iminus1]; kk < bpntr[i]; kk++)
           { 
             bindx[kk] = cur_data->columns[ll]/col_block_size;
             indx[kk+1] = indx[kk]+col_block_size*(rpntr[i] - rpntr[iminus1]);
             /*loop over a block and store data columns then rows*/
             for(jj = 0; jj < col_block_size; jj++)
             {
               for(ii = 0; ii < all_rpntr[k][j]; ii++) 
               {
                 vals[indx[kk]+ii+(cur_data->columns[kk*col_block_size+jj+(iminus1)*col_block_size]%col_block_size)*all_rpntr[k][j]] = cur_data->values[ll+jj+ii*(cur_data->rowptr[rpntr[iminus1]+1]-cur_data->rowptr[rpntr[iminus1]])];
               }
             }
             ll += col_block_size;
           }
         }
         k++;
       }
       ML_CSR_MSRdata_Destroy(cur_matrix->data);
       iplus1 = i+1;
       rpntr = (int*)realloc((void*)rpntr, iplus1*sizeof(int));
       bpntr = (int*)realloc((void*)bpntr, iplus1*sizeof(int));
       bindx = (int*)realloc((void*)bindx, bpntr[i]*sizeof(int));
       indx = (int*)realloc((void*)indx, (bpntr[i]+1)*sizeof(int));
        
       out_data = (struct ML_vbrdata *)ML_allocate(sizeof(struct ML_vbrdata)); 
       out_data->bindx = bindx;
       out_data->indx = indx;
       out_data->bpntr = bpntr;
       out_data->cpntr = cpntr;
       out_data->rpntr = rpntr;
       out_data->val = vals;
       cur_matrix->getrow->func_ptr = VBR_getrows;
       cur_matrix->getrow->data = (void *)out_data;
       cur_matrix->data = (void *)out_data;
       cur_matrix->blocks = cur_matrix->sub_matrix->blocks + bpntr[i];
       cur_matrix->data_destroy = ML_RECUR_VBRdata_Destroy; 
       for(j = 0; j < neighbors; j++)
         if (pre_comm->neighbors[j].N_rcv > 0)
           ML_free(all_rpntr[j]);
       ML_free(all_rpntr);
       ML_free(recv_lens);
       cur_data = NULL;
     }
     cur_matrix = in_matrix;
     if(cur_matrix->sub_matrix != NULL)
       cur_matrix->getrow->N_block_rows = i+cur_matrix->sub_matrix->getrow->N_block_rows;
   }
   else
   {
   /*settings to change since we now have a vbr matrix*/
   in_matrix->type = ML_TYPE_VBR_MATRIX;
   in_matrix->matvec->func_ptr = NULL;
   /*in_matrix->matvec->func_ptr = VBR_matvec;  this needs to be put back in at some point once the function exists*/

   /*find number of block rows and block columnsi and their location if using a fixed width*/
   if (row_block_size < 0)
   {
     pr_error("In function convert2vbr row_block_size is negative which is undefined.  Please see function header and pass in the appropriate value.\n");
   }
   if (row_block_size == 0)
   { 
     if(rpntr == NULL)
     {
       pr_error("In function convert2vbr rpntr is NULL when expecting row blocking data.\n");
     }
     i = rpntr[0];
     j = 0;
     while(i < in_matrix->outvec_leng)
     {
       blockrows++;
       i = rpntr[++j];
     }
   }
   else
   {
     blockrows = in_matrix->outvec_leng/row_block_size;
     if(in_matrix->outvec_leng != row_block_size*blockrows)
     {
       pr_error("In function convert2vbr the incoming row_block_size does not evenly divide the matrix rows.\n");
     }
     if(rpntr != NULL)
       ML_free(rpntr);
     rpntr = (int*)ML_allocate((blockrows+1) * sizeof(int));
     rpntr[0] = 0;
     for(i = 1; i <= blockrows; i++)
       rpntr[i] = rpntr[i-1] + row_block_size;
   }
   
   if (col_block_size < 0)
   {
     pr_error("In function convert2vbr col_block_size is negative which is undefined.  Please see function header and pass in the appropriate value.\n");
   }
   if (col_block_size == 0)
   { 
     Nghost_nodes = in_matrix->getrow->pre_comm->total_rcv_length;
     total_cols = in_matrix->invec_leng+Nghost_nodes;
     if(cpntr == NULL)
     {
       pr_error("In function convert2vbr cpntr is NULL when expecting column blocking data.\n");
     }
     i = cpntr[0];
     j = 0;
     while(i < in_matrix->invec_leng)
     {
       blockcolumns++;
       i = cpntr[++j];
     }
   }
   else
   { 
     ML_CommInfoOP_Compute_TotalRcvLength(in_matrix->getrow->pre_comm);
     Nghost_nodes = in_matrix->getrow->pre_comm->total_rcv_length;
     total_cols = in_matrix->invec_leng+Nghost_nodes;
     blockcolumns = total_cols/col_block_size;
     if(total_cols != col_block_size*blockcolumns)
     {
       pr_error("In function convert2vbr the incoming col_block_size does not evenly divide the matrix rows.\n");
     }
     if(cpntr != NULL)
       ML_free(cpntr);
     cpntr = (int*)ML_allocate((blockcolumns+1) * sizeof(int));
     cpntr[0] = 0;
     for(i = 1; i <= blockcolumns; i++)
       cpntr[i] = cpntr[i-1] + col_block_size;
   }

   /*at most we have the lesser of blockrows times blockcolumns blocks and the nnz in the matrix which might be better but we are not gareenteed to know the nnz this is a gross overestimate so finding a better prediction might be in order but nnz may not be right so we can't use that*/
   spaceneeded = blockcolumns*blockrows;
   /*if(in_matrix->N_nonzeros > 0)
   {
     if(in_matrix->N_nonzeros < spaceneeded)
       spaceneeded = in_matrix->N_nonzeros;
   }*/

   /*the number of blockrows is either exact or close enough however*/
   bpntr = (int*)ML_allocate((blockrows+1)*sizeof(int));


   /*Space for each point row being read in*/ 
   /*if(in_matrix->max_nz_per_row > 0)
     A_i_allocated = in_matrix->max_nz_per_row + 1; this is exact
   else*/
     A_i_allocated = total_cols; /*This is a very bad but safe estimate.*/
   A_i_cols  = (int    *) ML_allocate(A_i_allocated * sizeof(int) );
   accum_val = (double *) ML_allocate(A_i_allocated * sizeof(double));
   
   if(accum_val == NULL)
   { 
     pr_error("Not enough space to allocate %d doubles and ints in convert2vbr.\n", A_i_allocated);
   }


/*   if(in_matrix->N_nonzeros <= 0)*/ 
   /*We don't have a good estimate so lets get one this code is used since 
     N_nonzeros even when set is not always correct*/
   for(i=0; i<in_matrix->outvec_leng;i++)
   {  
     in_matrix->getrow->func_ptr(in_matrix, 1, &(i), A_i_allocated, A_i_cols, accum_val, &row_length);
     nnz += row_length;
   }
   spaceneeded = nnz;
   bindx = (int*)ML_allocate(spaceneeded*sizeof(int));
   indx = (int*)ML_allocate((spaceneeded+1)*sizeof(int));
   indx[0] = bpntr[0] = 0;
   /*else lets use our good estimate we already have
     nnz = in_matrix->N_nonzeros;*/

   /*10 is a complete guess.  One would hope the matrix resulting block matrix was more dense than this but there is no gareentee*/
   vals = (double *) ML_allocate(nnz*10*sizeof(double));
   i = 10;
   if(vals == NULL) /*if we don't have 10 times the space lets start trying smaller*/
   {
     for(i = 10; i > 1; i--)
     {
       vals = (double *) ML_allocate(nnz*i*sizeof(double));
       if(vals != NULL)
         break;
     }
   }
   /*entries avalible in the values array*/
   val_size = nnz*i;
   /*we don't have at least twice the nnz*/
   if(i == 1)
   {
     spacetight = 1;
     /*other arrays need to be made smaller since we're tight on space also its only going to work if the matrix blocks are fairly dense anyway*/
/*       if(in_matrix->N_nonzeros > 0)
       {
         if(nnz == spaceneeded)
           spaceneeded /= 
         else
           spaceneeded = sqrt(spaceneeded)*10;
           
       }*/
     /*keep making vals smaller till it fits*/
     x = 1.9;
     for(i = 9; i >= 0; i--)
     {
       temp_alc = (int) (nnz*x*sizeof(double));
       vals = (double *) ML_allocate(temp_alc);
       if(vals != NULL)
         break;
       x-=.1;
     }
     val_size = nnz*(10+i)/10;
     /*there just isn't enough space*/
     if(vals == NULL)
     {
       printf("Not enough space in ML_convert2vbr\n");
       printf("trying to allocate %d ints and %d doubles \n",A_i_allocated+spaceneeded*2+blockrows,in_matrix->N_nonzeros+A_i_allocated);
       printf("Matrix has %d rows \n", in_matrix->getrow->Nrows);
       pr_error("Matrix has %d nz\n", in_matrix->N_nonzeros);
     } 
   } 
   /*Converting functions start by looping over all block rows*/
   for(i = 0; i < blockrows; i++)
   {
     iplus1 = i+1;
     bpntr[iplus1] = bpntr[i];
     /*loop over point rows within blockrow*/
     for(j = rpntr[i]; j < rpntr[iplus1]; j++)
     {
       /*get row data*/
       in_matrix->getrow->func_ptr(in_matrix, 1, &(j), A_i_allocated, A_i_cols, accum_val, &row_length);
     
       /*loop over row data*/ 
       blockend = 1; 
       for(k = 0; k < row_length; k++)
       {
         /*loop over block data to find block we're in*/
         while(A_i_cols[k] >= cpntr[blockend]){
           blockend++;
         }
         blockstart = blockend-1;
         blockfound = 0;
         /*see if block already exists*/
         for (kk = bpntr[i]; kk < bpntr[iplus1]; kk++)
         {
           if(bindx[kk] == blockstart)
           {
             blockfound = 1;
             break;
           }
         }
         if(blockfound == 1) /*insert value*/
         {
           vals[indx[kk] + (A_i_cols[k]-cpntr[blockstart])*(rpntr[iplus1] - rpntr[i]) + j - rpntr[i]] = accum_val[k];
         }
         else /*create new block*/
         {
           bindx[bpntr[iplus1]] = blockend-1;
           bpntr[iplus1]++;
           indx[bpntr[iplus1]] = (cpntr[blockend]-cpntr[blockstart])*(rpntr[iplus1] - rpntr[i]) + indx[bpntr[iplus1]-1];
           if(indx[bpntr[iplus1]] >= val_size)/*if array is out of space try and make more if not fail*/
           {
             while(indx[bpntr[iplus1]] >= val_size)
             {
               val_size *= 3;
               val_size /= 2;
               temp_double = (double *)realloc(vals, val_size*sizeof(double));
               if(temp_double == NULL)
               {
                 printf("Not enough space in ML_convert2vbr to create a bigger values array\n");
                 printf("trying to allocate %d ints and %d doubles \n",A_i_allocated+spaceneeded*2+blockrows,in_matrix->N_nonzeros+A_i_allocated);
                 printf("Matrix has %d rows \n", in_matrix->getrow->Nrows);
                 pr_error("Matrix has %d nz\n", in_matrix->N_nonzeros);
               }
               else
                 vals = temp_double;
             }
           }
           /*initialize the new block to all zeros*/
           for(jj = indx[bpntr[iplus1]-1]; jj < indx[bpntr[iplus1]]; jj++)
           {
             vals[jj] = 0.0;
           }
           /*store value*/
           vals[indx[kk] + (A_i_cols[k]-cpntr[blockstart])*(rpntr[iplus1] - rpntr[i]) + j - rpntr[i]] = accum_val[k];
         }
       }
     }
   }

   /*we also need to set the block get row structure here*/
   in_matrix->getrow->func_ptr = VBR_getrows;
   

   /*We should reallocate to set the arrays to the actual sizes we know they now are*/
   bindx = (int*)realloc((void*)bindx, bpntr[iplus1]*sizeof(int)); 
   indx = (int*)realloc((void*)indx, (bpntr[iplus1]+1)*sizeof(int));
   vals = (double*)realloc((void*)vals, jj*sizeof(double));
 
   /*set real values to aliases*/
   out_data = (struct ML_vbrdata *)ML_allocate(sizeof(struct ML_vbrdata));
   out_data->bindx = bindx;
   out_data->indx = indx;
   out_data->bpntr = bpntr;
   out_data->cpntr = cpntr;
   out_data->rpntr = rpntr;
   out_data->val = vals;
   in_matrix->getrow->N_block_rows = blockrows;
   in_matrix->N_nonzeros = indx[bpntr[iplus1]];
   in_matrix->blocks = bpntr[iplus1];
   /*delete old csr matrix*/

   if ((in_matrix->data_destroy != NULL) && (in_matrix->data != NULL)) {
      in_matrix->data_destroy(in_matrix->data);
      in_matrix->data = NULL;
   }

   in_matrix->data = (void *)out_data;
   in_matrix->getrow->data = (void *)out_data;
   in_matrix->data_destroy = ML_RECUR_VBRdata_Destroy; 
   /*memory cleanup*/
   ML_free(A_i_cols);
   ML_free(accum_val);
 
   }
}

/************************************************************************/
/* For any accum_col[k] == accum_col[j] (k,j < Ncols), sum the          */
/* corresponding values and store them back in accum_val.               */
/* -------------------------------------------------------------------- */

void ML_sum_duplicates(int accum_col[], double accum_val[], int *Ncols)
{

   int i, new_length = 0;

   if ( *Ncols != 0) new_length++;
   
   for (i = 1; i < *Ncols ; i++ ) {
      if (accum_col[i] == accum_col[i-1]) 
         accum_val[new_length-1] += accum_val[i];
      else {
         accum_val[new_length] = accum_val[i];
         accum_col[new_length] = accum_col[i];
         new_length++;
      }
   }
   *Ncols = new_length;
}

/************************************************************************/
/* Routine used to expand the size of 'accum_val' and 'accum_col'.
 *
 * Parameters
 * ==========
 *   accum_size           On input, size of 'accum_col' and 'accum_val'
 *                        on exit to this routine.
 * 
 *   accum_col, accum_val On input, arrays contain 'Ncols' data values.
 *                        On output, arrays are expanded to contain 
 *                        'accum_size' elements preserving the 'Ncols' values.
 *
 *   Ncols                On input, number of data values stored in the
 *                        arrays 'accum_col' and 'accum_val'.
 * -------------------------------------------------------------------- */

void ML_expand_accum(int accum_size, int **accum_col, double **accum_val,
                     int Ncols)
{
   int    *itmp1, *itmp2, jj;
   double *dtmp1, *dtmp2;

   itmp1 = (int *) ML_allocate( accum_size * sizeof(int) );
   if (itmp1 == NULL) {
      printf("Out of space trying to expand accumulator\n");
      exit(1);
   }
   itmp2 = *accum_col;
   for (jj = 0; jj < Ncols; jj++) itmp1[jj] = itmp2[jj];
   ML_free(itmp2);
   *accum_col = itmp1;

   dtmp1 = (double *) ML_allocate( accum_size * sizeof(double) );
   if (dtmp1 == NULL) {
      printf("Out of space trying to expand accumulator\n");
      exit(1);
   }
   dtmp2 = *accum_val;
   for (jj = 0; jj < Ncols; jj++) dtmp1[jj] = dtmp2[jj];
   ML_free(dtmp2);
   *accum_val = dtmp1;
}

/* ******************************************************************** */
/* multiplying two matrices together                                    */
/* -------------------------------------------------------------------- */
void ML_2matmult(ML_Operator *Mat1, ML_Operator *Mat2,
                 ML_Operator *Result, int matrix_type)
{
   int N_input_vector, max_per_proc;
   ML_CommInfoOP *getrow_comm;
   ML_Operator   *Mat2comm, *Mat1Mat2, *tptr, *Mat1Mat2comm;
   ML_Comm       *comm;
   char          label1[80],label2[80];

   if (Mat1->invec_leng != Mat2->outvec_leng)
   {
     if (Mat1->label == NULL) sprintf(label1,"%s","mat1_not_labeled");
     else sprintf(label1,"%s",Mat1->label);
     if (Mat2->label == NULL) sprintf(label2,"%s","mat2_not_labeled");
     else sprintf(label2,"%s",Mat2->label);
     pr_error("In ML_2matmult: matrix dimensions do not agree:\n\tMat1->invec_leng = %d, Mat2->outvec_leng = %d, (%s & %s)\n", Mat1->invec_leng, Mat2->outvec_leng,label1,label2);
   }

   comm = Mat1->comm;
   N_input_vector = Mat2->invec_leng;
   getrow_comm    = Mat2->getrow->pre_comm;

   if (matrix_type != ML_EpetraCRS_MATRIX)
     ML_create_unique_col_id(N_input_vector, &(Mat2->getrow->loc_glob_map),
                             getrow_comm, &max_per_proc, comm);
   else
#ifdef ML_WITH_EPETRA
     ML_create_unique_col_id_exactoffset(N_input_vector, &(Mat2->getrow->loc_glob_map),
                                         getrow_comm, &max_per_proc, comm);
#else
     pr_error("ML_2matmult: ML_EpetraCRS_MATRIX requires epetra to be compiled in.\n");
#endif

   Mat2->getrow->use_loc_glob_map = ML_YES;

   if (matrix_type != ML_EpetraCRS_MATRIX)
     if (max_per_proc == 0 && comm->ML_mypid == 0) {
       pr_error("ERROR: In ML_2matmult, maximum number of local unknowns\n       on any processor (max_per_proc) is zero !\n");
     }

   if (Mat1->getrow->pre_comm != NULL)
      ML_exchange_rows( Mat2, &Mat2comm, Mat1->getrow->pre_comm);
   else Mat2comm = Mat2;
         
   ML_matmat_mult(Mat1, Mat2comm , &Mat1Mat2);

   ML_free(Mat2->getrow->loc_glob_map); Mat2->getrow->loc_glob_map = NULL;

   Mat2->getrow->use_loc_glob_map = ML_NO;
   if (Mat1->getrow->pre_comm != NULL) {
      tptr = Mat2comm;
      while ( (tptr!= NULL) && (tptr->sub_matrix != Mat2))
         tptr = tptr->sub_matrix;
      if (tptr != NULL) tptr->sub_matrix = NULL;
      ML_RECUR_CSR_MSRdata_Destroy(Mat2comm);
      ML_Operator_Destroy(&Mat2comm);
   }

   if (Mat1->getrow->post_comm != NULL) {
      ML_exchange_rows( Mat1Mat2, &Mat1Mat2comm, Mat1->getrow->post_comm);
   }
   else Mat1Mat2comm = Mat1Mat2;
   
   if (matrix_type == ML_CSR_MATRIX)
     ML_back_to_csrlocal(Mat1Mat2comm, Result, max_per_proc);
   else if (matrix_type == ML_MSR_MATRIX) {
     if (Mat1Mat2->invec_leng != Mat1Mat2->outvec_leng) 
       pr_error("ML_2matmult: MSR format only valid for square matrices.\n");
     ML_back_to_local(Mat1Mat2, Result, max_per_proc);
   }
   else if (matrix_type == ML_EpetraCRS_MATRIX)
#ifdef ML_WITH_EPETRA
     ML_back_to_epetraCrs(Mat1Mat2, Result, Mat1, Mat2);
#else
     pr_error("ML_2matmult: ML_EpetraCRS_MATRIX requires epetra to be compiled in.\n");
#endif
   else pr_error("ML_2matmult: Unknown matrix type\n");


   ML_RECUR_CSR_MSRdata_Destroy(Mat1Mat2comm);
   ML_Operator_Destroy(&Mat1Mat2comm);
}



/* ******************************************************************** */
/* multiplying two matrices together - w/ blocking                      */
/* -------------------------------------------------------------------- */
void ML_2matmult_block(ML_Operator *Mat1, ML_Operator *Mat2,
                 ML_Operator *Result, int matrix_type)
{
   int N_input_vector, max_per_proc;
   ML_CommInfoOP *getrow_comm;
   ML_Operator   *Mat2comm, *Mat1Mat2, *tptr, *Mat1Mat2comm;
   ML_Comm       *comm;
   char          label1[80],label2[80];

   if (Mat1->invec_leng != Mat2->outvec_leng)
   {
     if (Mat1->label == NULL) sprintf(label1,"%s","mat1_not_labeled");
     else sprintf(label1,"%s",Mat1->label);
     if (Mat2->label == NULL) sprintf(label2,"%s","mat2_not_labeled");
     else sprintf(label2,"%s",Mat2->label);
     pr_error("In ML_2matmult: matrix dimensions do not agree:\n\tMat1->invec_leng = %d, Mat2->outvec_leng = %d, (%s & %s)\n", Mat1->invec_leng, Mat2->outvec_leng,label1,label2);
   }

   comm = Mat1->comm;
   N_input_vector = Mat2->invec_leng;
   getrow_comm    = Mat2->getrow->pre_comm;

   if (matrix_type != ML_EpetraCRS_MATRIX)
     ML_create_unique_col_id(N_input_vector, &(Mat2->getrow->loc_glob_map),
                             getrow_comm, &max_per_proc, comm);
   else
#ifdef ML_WITH_EPETRA
     ML_create_unique_col_id_exactoffset(N_input_vector, &(Mat2->getrow->loc_glob_map),
                                         getrow_comm, &max_per_proc, comm);
#else
     pr_error("ML_2matmult: ML_EpetraCRS_MATRIX requires epetra to be compiled in.\n");
#endif

   Mat2->getrow->use_loc_glob_map = ML_YES;

   if (matrix_type != ML_EpetraCRS_MATRIX)
     if (max_per_proc == 0 && comm->ML_mypid == 0) {
       pr_error("ERROR: In ML_2matmult, maximum number of local unknowns\n       on any processor (max_per_proc) is zero !\n");
     }

   if (Mat1->getrow->pre_comm != NULL)
      ML_exchange_rows( Mat2, &Mat2comm, Mat1->getrow->pre_comm);
   else Mat2comm = Mat2;

   /*ML_convert2vbr(Mat2comm);*/
      
   ML_matmat_mult(Mat1, Mat2comm , &Mat1Mat2);

   ML_free(Mat2->getrow->loc_glob_map); Mat2->getrow->loc_glob_map = NULL;

   Mat2->getrow->use_loc_glob_map = ML_NO;
   if (Mat1->getrow->pre_comm != NULL) {
      tptr = Mat2comm;
      while ( (tptr!= NULL) && (tptr->sub_matrix != Mat2))
         tptr = tptr->sub_matrix;
      if (tptr != NULL) tptr->sub_matrix = NULL;
      ML_RECUR_CSR_MSRdata_Destroy(Mat2comm);
      ML_Operator_Destroy(&Mat2comm);
   }

   if (Mat1->getrow->post_comm != NULL) {
      ML_exchange_rows( Mat1Mat2, &Mat1Mat2comm, Mat1->getrow->post_comm);
   }
   else Mat1Mat2comm = Mat1Mat2;


   /* Copy over the num_PDEs and num_rigid info */
   Mat1Mat2comm->num_PDEs = Result->num_PDEs = Mat1->num_PDEs;
   Mat1Mat2comm->num_rigid = Result->num_rigid = Mat1->num_rigid;

   
   if (matrix_type == ML_CSR_MATRIX)
     ML_back_to_csrlocal(Mat1Mat2comm, Result, max_per_proc);
   else if (matrix_type == ML_MSR_MATRIX) {
     if (Mat1Mat2->invec_leng != Mat1Mat2->outvec_leng) 
       pr_error("ML_2matmult: MSR format only valid for square matrices.\n");
     ML_back_to_local(Mat1Mat2, Result, max_per_proc);
   }
   else if (matrix_type == ML_EpetraCRS_MATRIX)
#ifdef ML_WITH_EPETRA
     ML_back_to_epetraCrs(Mat1Mat2, Result, Mat1, Mat2);
#else
     pr_error("ML_2matmult: ML_EpetraCRS_MATRIX requires epetra to be compiled in.\n");
#endif
   else pr_error("ML_2matmult: Unknown matrix type\n");

   ML_RECUR_CSR_MSRdata_Destroy(Mat1Mat2comm);
   ML_Operator_Destroy(&Mat1Mat2comm);
}


/*
Memory improvements:
1) decide how many nonzeros to get in Bmatrix.
    if we can hold twice the size of B don't worry be happy.
    otherwise, take half of what we can hold.


2) estimate how many rows of B we can hold.
    just divide by the avg_per_row (maybe we can update the avg).

3) scan rows of A and make a list of corresponding rows of B
   until we hit the number of rows in B that we can hold. 
*/
int ML_determine_Brows(int start, int *end, ML_Operator *Amatrix,
		       int *irows[], int *irows_length, int *NBrows,
		       int *rows_that_fit, 
                void   (*Agetrow)(ML_Operator *,int,int *,int *,int **,
                       double **,int *,int))
{
  int i, j, rowi_N, hash_val = 0, N, *rows, rows_length,kk;
  int A_i_allocated = 0, *A_i_cols = NULL, hash_used;
  double *A_i_vals = NULL;

  rows = *irows;
  rows_length = *irows_length;
  N = Amatrix->getrow->Nrows;
  for (i = 0; i < rows_length; i++) rows[i] = -1;
  i = start;
  j = 0;
  rowi_N = 0;
  hash_used = 0;

  while ( *NBrows < *rows_that_fit ) {
    if (j < rowi_N) {
      ML_hash_it(A_i_cols[j], rows, rows_length, &hash_used, &hash_val);
      if (rows[hash_val] == -1) {
        (*NBrows)++;
        if ( *NBrows == *rows_that_fit) {
          if ( (j+1 < rowi_N) && (i-1 == start)) {
            (*rows_that_fit)++;
            if ( *rows_that_fit > rows_length) {
              (*irows_length) += 5;
              *irows = (int *) ML_allocate((*irows_length)*sizeof(int));
              if (*irows == NULL) pr_error("matmat: out of space\n");
              for (kk = 0; kk < rows_length; kk++) (*irows)[kk]=rows[kk];
              for (kk = rows_length; kk < *irows_length; kk++) (*irows)[kk]=-1;
              ML_free(rows);
              rows = *irows;
              rows_length = *irows_length;
            }
          }
        }
      }
      rows[hash_val] = A_i_cols[j++];
    } else {
      if (i == N) *rows_that_fit = -(*rows_that_fit);
      else {
        Agetrow(Amatrix,1, &i, &A_i_allocated, &A_i_cols, &A_i_vals, &rowi_N,0);
        i++;
        j = 0;
      }
    }
  }
  if (*rows_that_fit < 0) { *rows_that_fit = -(*rows_that_fit);}
  if (j != rowi_N) i--;
  *end = i;
  j = 0;
  for (i = 0; i < rows_length; i++) {
    if ( rows[i] != -1) rows[j++] = rows[i];
  }
  return 0;
} /* ML_determine_Brows() */

/* ******************************************************************** */
int ML_determine_Bblkrows(int start, int *end, ML_Operator *Amatrix,
               int *irows[], int *irows_length, int *NBrows,
               int *rows_that_fit, 
               void   (*Agetrow)(ML_Operator *,int,int *,int *,int **,
                       int **,int *,int))
/* -------------------------------------------------------------------- */
{
  int i, j, rowi_N, hash_val = 0, N, *rows, rows_length,kk;
  int A_i_allocated = 0, *A_i_cols = NULL, hash_used;
  int *A_i_vals = NULL;

  rows = *irows;
  rows_length = *irows_length;
  N = Amatrix->getrow->Nrows;
  for (i = 0; i < rows_length; i++) rows[i] = -1;
  i = start;
  j = 0;
  rowi_N = 0;
  hash_used = 0;

  while ( *NBrows < *rows_that_fit ) {
    if (j < rowi_N) {
      ML_hash_it(A_i_cols[j], rows, rows_length, &hash_used, &hash_val);
      if (rows[hash_val] == -1) {
        (*NBrows)++;
        if ( *NBrows == *rows_that_fit) {
          if ( (j+1 < rowi_N) && (i-1 == start)) {
            (*rows_that_fit)++;
            if ( *rows_that_fit > rows_length) {
              (*irows_length) += 5;
              *irows = (int *) ML_allocate((*irows_length)*sizeof(int));
              if (*irows == NULL) pr_error("matmat: out of space\n");
              for (kk = 0; kk < rows_length; kk++) (*irows)[kk]=rows[kk];
              for (kk = rows_length; kk < *irows_length; kk++) (*irows)[kk]=-1;
              ML_free(rows);
              rows = *irows;
              rows_length = *irows_length;
            }
          }
        }
      }
      rows[hash_val] = A_i_cols[j++];
    } else {
      if (i == N) *rows_that_fit = -(*rows_that_fit);
      else {
        Agetrow(Amatrix,1, &i, &A_i_allocated, &A_i_cols, &A_i_vals,&rowi_N,0);
        i++;
        j = 0;
      }
    }
  }
  if (*rows_that_fit < 0) { *rows_that_fit = -(*rows_that_fit);}
  if (j != rowi_N) i--;
  *end = i;
  j = 0;
  for (i = 0; i < rows_length; i++) {
    if ( rows[i] != -1) rows[j++] = rows[i];
  }
  return 0;
} /* ML_determine_Bblkrows() */
