/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ml_struct.h"
#include "ml_rap.h"
#include "ml_memory.h"

#define HASH_SIZE     1024   /* should be a power of two */
#define MASK          01777  /* should be bit pattern of hash size minus 1 */


/* ******************************************************************** */
/* matrix matrix multiplication                                         */ 
/* Cmatrix = Amatrix * Bmatrix                                          */
/*                                                                      */
/* Note: Matrices can be stored in chunks. See ml_rap.h for a           */
/*       description of the matrix structure ML_matrix.                 */
/* -------------------------------------------------------------------- */

void ML_matmat_mult(ML_Operator *Amatrix, ML_Operator *Bmatrix,
                    ML_Operator **Cmatrix)
{
   int    i,k, jj, next_nz, Ncols, N, Nnz_estimate, sub_i, accum_size, row;
   int    *C_ptr, *Ccol, *A_i_cols, rowi_N, *accum_col, row2_N;
   double *Cval, *A_i_vals, multiplier, *accum_val, dtemp;
   ML_Operator *current, *previous_matrix, *next;
   struct ML_CSR_MSRdata *temp;
   int    max_nz_row_new = 0, total_nz = 0, index_length = 0;
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
   /*
   double t1, t6;
   */
   int tcols, hash_used, j, *tptr;
   int *acc_col_ptr, *Bcol_ptr; double *acc_val_ptr, *Bval_ptr;
   /*
   t1 = GetClock();
   */
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
	  printf("%d: recomputing rcv length (%u %u)\n",Bmatrix->comm->ML_mypid,Bmatrix,current); fflush(stdout);
#endif
	  ML_CommInfoOP_Compute_TotalRcvLength(current->getrow->pre_comm);
	}

         Next_est += current->getrow->pre_comm->total_rcv_length;
#ifdef charles
	  printf("%d: Nghost = %d  %d\n",Bmatrix->comm->ML_mypid,Next_est,current->getrow->pre_comm->total_rcv_length); fflush(stdout);
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
   /*------------------------------------------------------------------------*/

   if ( Amatrix->getrow->Nrows > 0 )
   {
      row = 0;
      Agetrow(Amatrix,1, &row, &A_i_allocated , &A_i_cols, &A_i_vals, &i,0);
      row = (Amatrix->getrow->Nrows-1)/2;
      Agetrow(Amatrix,1, &row, &A_i_allocated , &A_i_cols, &A_i_vals, &k,0);
      row = Amatrix->getrow->Nrows-1;
      Agetrow(Amatrix,1, &row, &A_i_allocated , &A_i_cols, &A_i_vals, &jj,0);
      A_avg_nz_per_row = ((double) (i+k+jj))/3.0;
   } else A_avg_nz_per_row = 100;

   if ( Bmatrix->getrow->Nrows > 0 )
   {
      row = 0;
      Bgetrow(Bmatrix,1,&row, &accum_size, &accum_col, &accum_val, &i,0);
      row = (Bmatrix->getrow->Nrows-1)/2;
      Bgetrow(Bmatrix,1,&row, &accum_size, &accum_col, &accum_val, &k,0);
      row = Bmatrix->getrow->Nrows-1;
      Bgetrow(Bmatrix,1,&row, &accum_size, &accum_col, &accum_val,&jj,0);
      B_avg_nz_per_row = ((double) (i+k+jj))/3.0;
   } else B_avg_nz_per_row = i = k = jj = 100;

   if (i  > Bmatrix->max_nz_per_row) Bmatrix->max_nz_per_row =  i+10;
   if (k  > Bmatrix->max_nz_per_row) Bmatrix->max_nz_per_row =  k+10;
   if (jj > Bmatrix->max_nz_per_row) Bmatrix->max_nz_per_row = jj+10;

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

   if (Bmatrix->N_nonzeros <= 0) {
      B_total_Nnz = 0;
      for (i = 0; i < Bmatrix->getrow->Nrows; i++ ) {
         Bgetrow(Bmatrix,1,&i,&accum_size, &accum_col, &accum_val, &row2_N, 0);
         B_total_Nnz += row2_N;
      }
   }
   else B_total_Nnz = Bmatrix->N_nonzeros;

   index_length = 2*accum_size;
   dtemp = 2.*((double) accum_size);
   while (index_length < 0 ) { /* overflow */
     dtemp *= 1.7;
     index_length = (int) dtemp;
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
      /*
if ((lots_of_space < 7) ) Bvals = NULL; else 
      */
  /*
if ((lots_of_space < 4) && (B_allocated > 500)) Bvals = NULL; else 
  */
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
   for (i = 0; i < subB_Nnz; i++) {
     if (hash_used >= ((int) (.75 * index_length)) ) {
#ifdef charles
       printf("%d: running out of hashing space: row = %d, tcols=%d"
              " hash_length=%d, hash_used = %d\n",
	          Bmatrix->comm->ML_mypid,i,tcols,index_length,hash_used);
       fflush(stdout);
#endif
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
	 hash_val = ML_hash_it(Bcols[j], col_inds, index_length, &hash_used);
	 if (col_inds[hash_val] == -1) tcols++;
	 col_inds[hash_val] = Bcols[j];
	 Bcols[j] = hash_val;
       }
     }
       
     hash_val = ML_hash_it(Bcols[i], col_inds, index_length, &hash_used);
     if (col_inds[hash_val] == -1) tcols++;
     col_inds[hash_val] = Bcols[i];
     Bcols[i] = hash_val;
   }
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
#ifdef takeout
   for (k = 0; k < accum_size; k++) accum_val[k] = 0.;
#endif

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
#ifdef takeout
      for (k = 0; k < rowi_N; k++) {
	jj = Bptr[A_i_cols[k]];
	Bcol_ptr = &(Bcols[jj]);
	while (jj++ < Bptr[A_i_cols[k]+1]) {
	  *acc_col_ptr = *Bcol_ptr;  acc_col_ptr++;
	     accum_index[*Bcol_ptr++] = Ncols++;
	}
      }
      for (k = 0; k < rowi_N; k++) {
	multiplier = A_i_vals[k];
	if ( multiplier != 0.0 ) {
	  jj = Bptr[A_i_cols[k]];
	  Bval_ptr = &(Bvals[jj]);
	  Bcol_ptr = &(Bcols[jj]);
	  while (jj++ < Bptr[A_i_cols[k]+1]) {
	    accum_val[accum_index[*Bcol_ptr++]] += multiplier*(*Bval_ptr++);
	  }
	}
      }
      k = 0;
      for (jj = 0; jj < Ncols; jj++ ) {
	if (jj == accum_index[accum_col[jj]]) {
	  accum_val[k] = accum_val[jj];
	  accum_col[k++] = accum_col[jj];
	}
      }
      for (jj = k; jj < Ncols; jj++) accum_val[jj] = 0.;
      Ncols = k;
#endif
      /***********************************************************************/
      /* Convert back to the original column indices.                        */
      /*---------------------------------------------------------------------*/

      acc_col_ptr = accum_col;
      for (jj = 0; jj < Ncols; jj++ ) {
	accum_index[*acc_col_ptr] = -1;
	*acc_col_ptr = col_inds[*acc_col_ptr];
	acc_col_ptr++;
      }

      /* empty row. Let's just put a zero in the first column */

      if (Ncols == 0) {
         accum_col[Ncols] = 0;
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

      memcpy(&(Ccol[next_nz]),accum_col, sizeof(int)*Ncols);
      memcpy(&(Cval[next_nz]),accum_val, sizeof(double)*Ncols);
      next_nz += Ncols;

#ifdef out
      /* above code might be a bit faster??? */
      for (k = 0; k < Ncols; k++) {
	/* This 'if' might break some applications somewhere */
	/* but I can't remember who and where or why?        */
	/* For now, I want to reduce memory in alegra so I am*/
	/* putting it in. If we later want to take this out  */
	/* we should use the memcpy code above.              */
	if (accum_val[k] != 0.0) {
          Ccol[next_nz] = accum_col[k];
          Cval[next_nz++] = accum_val[k];
	}
      }
#endif
      /*      */
#ifdef takeout
for (jj = 0; jj < Ncols; jj++) accum_val[jj] = 0.;
#endif
      C_ptr[sub_i+1] = next_nz;
      sub_i++;
      if (Ncols > max_nz_row_new) max_nz_row_new = Ncols;
   }
   start = end;
   }
   if (Bvals != NULL) ML_free(Bvals);
   if (Bcols != NULL) ML_free(Bcols);
   if (Bptr != NULL) ML_free(Bptr);
   if (accum_index != NULL) ML_free(accum_index);
   if (col_inds != NULL) ML_free(col_inds);
   ML_free(accum_col);
   ML_free(accum_val);
   ML_free(A_i_vals);
   ML_free(A_i_cols);
   /*
   t6 = GetClock();
   printf("matmat  ==> %e   %d\n", t6-t1,Bmatrix->comm->ML_mypid); fflush(stdout);
   */

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

/********************************************************************/
/********************************************************************
 * Cmatrix = Amatrix * Bmatrix
 *
 * Note: Matrices can be stored in chunks. See ml_rap.h for a description
 *       of the matrix structure ML_matrix.
 * -------------------------------------------------------------------- */

void ML_oldmatmat_mult(ML_Operator *Amatrix, ML_Operator *Bmatrix,
                    ML_Operator **Cmatrix)
{

   int    i,k, jj, next_nz, Ncols, N, total, sub_i, accum_size, row;
   int    *C_ptr, *Ccol, *rowi_col, rowi_N, *accum_col, row2_N;
   double *Cval, *rowi_val, multiplier, *accum_val, *d2temp, dtemp;
   ML_Operator *current, *previous_matrix;
   struct ML_CSR_MSRdata *temp;
   int    max_nz_per_row, max_nz_row_new = 0, total_nz = 0;
   double A_avg_nz_per_row, B_avg_nz_per_row, estimated_nz_per_row;
   int allocated;
   *Cmatrix = NULL;

   /********************************************************/
   /* put here to alleviate a hidden error - C. Tong       */
   if ( Amatrix->max_nz_per_row > Amatrix->getrow->Nrows )
      Amatrix->max_nz_per_row = Amatrix->getrow->Nrows;
   if ( Bmatrix->max_nz_per_row > Bmatrix->getrow->Nrows )
      Bmatrix->max_nz_per_row = Bmatrix->getrow->Nrows;
   /********************************************************/

   previous_matrix = NULL;
   max_nz_per_row  = Bmatrix->max_nz_per_row + 1;
   N      = Amatrix->getrow->Nrows;
   sub_i  = 0;

   accum_size = 10*Bmatrix->max_nz_per_row+ 1;

   /* allocate space to hold rows of Amatrix and the accumulator */

   allocated = Amatrix->max_nz_per_row + 1;
   rowi_col  = (int    *) ML_allocate(allocated * sizeof(int) );
   rowi_val  = (double *) ML_allocate(allocated * sizeof(double));
   accum_col = (int    *) ML_allocate( accum_size * sizeof(int) );
   accum_val = (double *) ML_allocate( accum_size * sizeof(double) );
   if ( (rowi_val == NULL) || (accum_val == NULL) || (rowi_col==NULL) ||
        (accum_col == NULL) ) {
      printf("Not enough space in ML_matmatmult().\n");
      printf("trying to allocate %d %d elements \n",allocated,accum_size);
      printf("Left  matrix has %d rows \n", Amatrix->getrow->Nrows);
      printf("Right matrix has %d rows \n", Bmatrix->getrow->Nrows);
      printf("Left  matrix has %d nz per row\n", Amatrix->max_nz_per_row);
      printf("Right matrix has %d nz per row\n", Bmatrix->max_nz_per_row);
      exit(1);
   }

   /* Make conservative estimates as to the size of the resulting  */
   /* matrix and the size needed for the accumulator.              */
   /* NOTE: These arrays can be increased later in the computation */

   if ( Amatrix->getrow->Nrows > 0 )
   {
      row = 0;
      ML_get_matrix_row(Amatrix,1, &row, &allocated , &rowi_col, &rowi_val, &i,0);
      row = (Amatrix->getrow->Nrows-1)/2;
      ML_get_matrix_row(Amatrix,1, &row, &allocated , &rowi_col, &rowi_val, &k,0);
      row = Amatrix->getrow->Nrows-1;
      ML_get_matrix_row(Amatrix,1, &row, &allocated , &rowi_col, &rowi_val, &jj,0);
      A_avg_nz_per_row = ((double) (i+k+jj))/3.0;
   } else A_avg_nz_per_row = 100;

   if ( Bmatrix->getrow->Nrows > 0 )
   {
      row = 0;
      ML_get_matrix_row(Bmatrix,1,&row, &accum_size, &accum_col, &accum_val, &i,0);
      row = (Bmatrix->getrow->Nrows-1)/2;
      ML_get_matrix_row(Bmatrix,1,&row, &accum_size, &accum_col, &accum_val, &k,0);
      row = Bmatrix->getrow->Nrows-1;
      ML_get_matrix_row(Bmatrix,1,&row, &accum_size, &accum_col, &accum_val,&jj,0);
      B_avg_nz_per_row = ((double) (i+k+jj))/3.0;
   } else B_avg_nz_per_row = i = k = jj = 100;

   if (i  > Bmatrix->max_nz_per_row) Bmatrix->max_nz_per_row =  i+10;
   if (k  > Bmatrix->max_nz_per_row) Bmatrix->max_nz_per_row =  k+10;
   if (jj > Bmatrix->max_nz_per_row) Bmatrix->max_nz_per_row = jj+10;

   estimated_nz_per_row  = sqrt(A_avg_nz_per_row) + sqrt(B_avg_nz_per_row) - 1.;
   estimated_nz_per_row *= estimated_nz_per_row;
   total = (int) (((double) Amatrix->getrow->Nrows)*estimated_nz_per_row * .75) +
           100;

   if (total <= Bmatrix->max_nz_per_row) total = Bmatrix->max_nz_per_row + 1;

   /* Allocate space for the matrix. Since 'total' is just an */
   /* estimate of the space needed, we will reduce 'total' if */
   /* we are unsuccessful allocating space.                   */

   C_ptr     = (int    *) ML_allocate((N+1)* sizeof(int) );
   if (C_ptr == NULL) pr_error("ML_matmat_mult: No space for C_ptr\n");
   Cval = NULL; Ccol = NULL;
   while ( (Cval == NULL) && (total > Bmatrix->max_nz_per_row) ) {
      if (Ccol != NULL) ML_free(Ccol);
      Ccol  = (int    *) ML_allocate( total* sizeof(int) );
      Cval  = (double *) ML_allocate( total* sizeof(double));
      if (Cval == NULL) total = total/2;
   }
   if (Cval == NULL) {
      printf("Not enough space for new matrix in ML_matmatmult().\n");
      printf("trying to allocate %d elements \n",total);
      printf("Left  matrix has %d rows \n", Amatrix->getrow->Nrows);
      printf("Right matrix has %d rows \n", Bmatrix->getrow->Nrows);
      printf("Left  matrix has %d nz per row\n", Amatrix->max_nz_per_row);
      printf("Right matrix has %d nz per row\n", Bmatrix->max_nz_per_row);
      exit(1);
   }

   next_nz   = 0;
   C_ptr[0]  = next_nz;

   for (i = 0; i < N ; i++) {

      /* Compute a new row in the matrix */

      Ncols = 0;
      ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, 
                        &rowi_val, &rowi_N, 0);

      for (k = 0; k < rowi_N; k++) {
         multiplier = rowi_val[k];
         if ( multiplier == 0.0 ) continue; /* added by Tong */
         if (Ncols + max_nz_per_row > accum_size) {
            ML_az_sort(accum_col, Ncols, NULL, accum_val);
            ML_sum_duplicates(accum_col, accum_val, &Ncols);
            /*if (Ncols + max_nz_per_row > accum_size/3) { deleted by Tong*/
               accum_size += (Ncols + max_nz_per_row)*20;
               ML_expand_accum(accum_size, &accum_col, &accum_val, Ncols);
            /*} deleted by Tong */
         }

         ML_get_matrix_row(Bmatrix, 1, &(rowi_col[k]), &accum_size, &accum_col,
                           &accum_val, &row2_N, Ncols);

         d2temp = &(accum_val[Ncols]);
         for (jj = 0; jj < row2_N; jj++) d2temp[jj] *= multiplier;
         Ncols = Ncols + row2_N;
      }
      ML_az_sort(accum_col, Ncols, NULL, accum_val);
      ML_sum_duplicates(accum_col, accum_val, &Ncols);

      /* empty row. Let's just put a zero in the first column */
      if (Ncols == 0) {
         accum_col[Ncols] = 0;
         accum_val[Ncols++] = 0.0;
      }

      /* check if we have enough space to store the new matrix row */
      /* If not, place the current matrix into a submatrix and     */
      /* allocate new vectors (hopefully large enough to hold the  */
      /* rest of the matrix) and a 'parent' matrix to hold the     */
      /* remaining rows.                                           */

      if (next_nz+Ncols > total) {

         /* create sub_matrix object */

         total_nz += next_nz;

         temp = (struct ML_CSR_MSRdata *) 
                ML_allocate(sizeof(struct ML_CSR_MSRdata));
	 if (temp == NULL) pr_error("ML_matmat_mult: no space for temp\n");
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
         current->N_nonzeros     = total_nz; 
         current->sub_matrix   = previous_matrix;

         /* allocate space for new matrix */

         if (i != 0) {
            dtemp = ((double) N-i)/ ((double) i);
            dtemp *= (double) total_nz;
            total = (int)(1.1*dtemp);
            total += Ncols;
/*
            total = ((int) (1.1 * ((double) ((N-i)* total_nz))/
                                ((double) i))) + Ncols;
*/
         }
         else total = total*N + Ncols;

         C_ptr = (int    *) ML_allocate( (N-i+1)* sizeof(int) );
         Ccol  = (int    *) ML_allocate( total* sizeof(int) );
         Cval  = (double *) ML_allocate( total* sizeof(double));
         if ((Cval == NULL) || (Ccol == NULL))  {
            printf("Not enough space for matrix\n");
            exit(1);
         }
         next_nz   = 0;
         C_ptr[0]  = next_nz;

         previous_matrix = current;
         sub_i = 0;
      }

      /* store matrix row */

      for (k = 0; k < Ncols; k++) {
         Ccol[next_nz] = accum_col[k];
         Cval[next_nz++] = accum_val[k];
      }
      C_ptr[sub_i+1] = next_nz;
      sub_i++;
      if (Ncols > max_nz_row_new) max_nz_row_new = Ncols;
   }
   ML_free(accum_col);
   ML_free(accum_val);
   ML_free(rowi_val);
   ML_free(rowi_col);

   /* create 'parent' object corresponding to the resulting matrix */

   total_nz += next_nz;
   temp = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata));
   if (temp == NULL) pr_error("ML_matmat_mult: no space for temp2\n");
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
   (*Cmatrix)->N_nonzeros     = total_nz;
   (*Cmatrix)->sub_matrix     = previous_matrix;
   if (allocated-1 > Amatrix->max_nz_per_row) 
      Amatrix->max_nz_per_row = allocated;

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

   if (Mat1->invec_leng != Mat2->outvec_leng)
   {
     printf("In ML_2matmult: matrix dimensions do not agree:\n");
     printf("\tMat1->invec_leng = %d, Mat2->outvec_leng = %d\n",
                       Mat1->invec_leng, Mat2->outvec_leng);
      exit(1);
   }

   comm = Mat1->comm;
   N_input_vector = Mat2->invec_leng;
   getrow_comm    = Mat2->getrow->pre_comm;

   ML_create_unique_col_id(N_input_vector, &(Mat2->getrow->loc_glob_map),
                           getrow_comm, &max_per_proc, comm);
   Mat2->getrow->use_loc_glob_map = ML_YES;

   if (max_per_proc == 0 && comm->ML_mypid == 0) {
     printf("ERROR: In ML_2matmult, maximum number of local unknowns\n       on any processor (max_per_proc) is zero !\n");
     exit(1);
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
int ML_hash_init(int hash_list[], int hash_length, int *hash_used)
{
  int i;

  for (i = 0; i < hash_length; i++) hash_list[i] = -1;
  *hash_used = 0;
  return 0;
}

int ML_hash_it( int new_val, int hash_list[], int hash_length,int *hash_used) {

  int index;
#ifdef charles
  int origindex;
#endif

  index = new_val<<1;
  if (index < 0) index = new_val;
  index = index%hash_length;
#ifdef charles
  origindex = index;
#endif
  while ( hash_list[index] != new_val) {
    if (hash_list[index] == -1) { (*hash_used)++; break;}
    /* JJH */
    /*index = (++index)%hash_length;*/
    index++;
    index = (index)%hash_length;
    /* --JJH */
#ifdef charles
    if (origindex == index)
       fprintf(stderr,"ML_hash_it: looped around original index = %d, new_val = %d, hash_length = %d, hash_used = %d\n");
#endif
  }

  return(index);
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

#ifdef charles
   printf("%d: in ML_determine\n",Amatrix->comm->ML_mypid);
   fflush(stdout);
#endif
   while ( *NBrows < *rows_that_fit ) {
      if (j < rowi_N) {
         hash_val = ML_hash_it(A_i_cols[j], rows, rows_length, &hash_used);
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
      }
      else {
        if (i == N) *rows_that_fit = -(*rows_that_fit);
        else {
           Agetrow(Amatrix,1, &i, &A_i_allocated, &A_i_cols, &A_i_vals, 
                   &rowi_N, 0);
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
#ifdef charles
   printf("%d: leaving ML_determine\n",Amatrix->comm->ML_mypid);
   fflush(stdout);
#endif
   return 0;
} 
/*    
        

4) scan in B and change the col indices.
5) in the main for loop put an additional loop around it
   do those rows. 
6) repeat 3) 
7) compress the old B to the front of the array
8) starting from the back of the array ...either scan or copy
in the right rows of B.
*/

#ifdef previousmatmat

/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ml_struct.h"
#include "ml_rap.h"

#define HASH_SIZE     1024   /* should be a power of two */
#define MASK          01777  /* should be bit pattern of hash size minus 1 */


/* ******************************************************************** */
/* matrix matrix multiplication                                         */ 
/* Cmatrix = Amatrix * Bmatrix                                          */
/*                                                                      */
/* Note: Matrices can be stored in chunks. See ml_rap.h for a           */
/*       description of the matrix structure ML_matrix.                 */
/* -------------------------------------------------------------------- */

void ML_matmat_mult(ML_Operator *Amatrix, ML_Operator *Bmatrix,
                    ML_Operator **Cmatrix)
{
   int    i,k, jj, next_nz, Ncols, N, total, sub_i, accum_size, row;
   int    *C_ptr, *Ccol, *rowi_col, rowi_N, *accum_col, row2_N;
   int    *Tcol, array_flag;
   double *Cval, *rowi_val, multiplier, *accum_val, *d2temp, dtemp, *Tval;
   ML_Operator *current, *previous_matrix, *next;
   struct ML_CSR_MSRdata *temp;
   int    max_nz_row_new = 0, total_nz = 0;
   double A_avg_nz_per_row, B_avg_nz_per_row, estimated_nz_per_row;
   int    allocated;
   int    hash, index, flag;
   void   (*fun1)(ML_Operator *,int,int *,int *,int **,double **,int *,int);
   void   (*fun2)(ML_Operator *,int,int *,int *,int **,double **,int *,int);
   int    Ncols2, kk;
   int    last, *newwhere, tflag, start;
   double *val;

   fun1   = ML_get_matrix_row;
   fun2   = ML_get_matrix_row;

   flag = 1;
   next = Bmatrix;
   while ( (next != NULL) ) {
      if (next->getrow->func_ptr != CSR_getrow) flag = 0;
      next = next->sub_matrix;
   }
   if ( (flag == 1)) { 
      fun2 = ML_get_matrow_CSR;

      next = Bmatrix;
      while ( (next != NULL) ) {
         if (next->getrow->row_map!= NULL) flag = 0;
         next = next->sub_matrix;
      }
      if ( (flag == 1)) { fun2 = ML_get_row_CSR_norow_map; }
   }
   flag = 1;
   next = Amatrix;
   while ( (next != NULL) ) {
      if (next->getrow->func_ptr != CSR_getrow) flag = 0;
      next = next->sub_matrix;
   }
   if ( (flag == 1)) { fun1 = ML_get_matrow_CSR; }

   *Cmatrix = NULL;

    /********************************************************/
   /* put here to alleviate a hidden error - C. Tong       */
   if ( Amatrix->max_nz_per_row > Amatrix->getrow->Nrows )
      Amatrix->max_nz_per_row = Amatrix->getrow->Nrows;
   if ( Bmatrix->max_nz_per_row > Bmatrix->getrow->Nrows )
      Bmatrix->max_nz_per_row = Bmatrix->getrow->Nrows;
   /********************************************************/

   previous_matrix = NULL;
   N      = Amatrix->getrow->Nrows;
   sub_i  = 0;

   accum_size = 10*Bmatrix->max_nz_per_row+ 100;

   /* allocate space to hold rows of Amatrix and the accumulator */

   newwhere = (int *) ML_allocate(sizeof(int)*HASH_SIZE);
   if ( newwhere == NULL) {
      printf("Not enough space for hash table in ML_matmatmult().\n");
      exit(1);
   }
   for (i = 0; i < HASH_SIZE; i++) newwhere[i] = -1;

   allocated = Amatrix->max_nz_per_row + 1;
   rowi_col  = (int    *) ML_allocate(allocated * sizeof(int) );
   rowi_val  = (double *) ML_allocate(allocated * sizeof(double));
   accum_col = (int    *) ML_allocate( accum_size * sizeof(int) );
   accum_val = (double *) ML_allocate( accum_size * sizeof(double) );
   if ( (rowi_val == NULL) || (accum_val == NULL)) {
      printf("Not enough space in ML_matmatmult().\n");
      printf("trying to allocate %d %d elements \n",allocated,accum_size);
      printf("Left  matrix has %d rows \n", Amatrix->getrow->Nrows);
      printf("Right matrix has %d rows \n", Bmatrix->getrow->Nrows);
      printf("Left  matrix has %d nz per row\n", Amatrix->max_nz_per_row);
      printf("Right matrix has %d nz per row\n", Bmatrix->max_nz_per_row);
      exit(1);
   }

   /* Make conservative estimates as to the size of the resulting  */
   /* matrix and the size needed for the accumulator.              */
   /* NOTE: These arrays can be increased later in the computation */

   if ( Amatrix->getrow->Nrows > 0 )
   {
      row = 0;
      ML_get_matrix_row(Amatrix,1, &row, &allocated , &rowi_col, &rowi_val, &i,0);
      row = (Amatrix->getrow->Nrows-1)/2;
      ML_get_matrix_row(Amatrix,1, &row, &allocated , &rowi_col, &rowi_val, &k,0);
      row = Amatrix->getrow->Nrows-1;
      ML_get_matrix_row(Amatrix,1, &row, &allocated , &rowi_col, &rowi_val, &jj,0);
      A_avg_nz_per_row = ((double) (i+k+jj))/3.0;
   } else A_avg_nz_per_row = 100;

   if ( Bmatrix->getrow->Nrows > 0 )
   {
      row = 0;
      ML_get_matrix_row(Bmatrix,1,&row, &accum_size, &accum_col, &accum_val, &i,0);
      row = (Bmatrix->getrow->Nrows-1)/2;
      ML_get_matrix_row(Bmatrix,1,&row, &accum_size, &accum_col, &accum_val, &k,0);
      row = Bmatrix->getrow->Nrows-1;
      ML_get_matrix_row(Bmatrix,1,&row, &accum_size, &accum_col, &accum_val,&jj,0);
      B_avg_nz_per_row = ((double) (i+k+jj))/3.0;
   } else B_avg_nz_per_row = i = k = jj = 100;

   if (i  > Bmatrix->max_nz_per_row) Bmatrix->max_nz_per_row =  i+10;
   if (k  > Bmatrix->max_nz_per_row) Bmatrix->max_nz_per_row =  k+10;
   if (jj > Bmatrix->max_nz_per_row) Bmatrix->max_nz_per_row = jj+10;

   estimated_nz_per_row  = sqrt(A_avg_nz_per_row) + sqrt(B_avg_nz_per_row) - 1.;
   estimated_nz_per_row *= estimated_nz_per_row;
   total = (int) (((double) Amatrix->getrow->Nrows)*estimated_nz_per_row * .75) +
           100;

   if (total <= Bmatrix->max_nz_per_row) total = Bmatrix->max_nz_per_row + 1;

   /* Allocate space for the matrix. Since 'total' is just an */
   /* estimate of the space needed, we will reduce 'total' if */
   /* we are unsuccessful allocating space.                   */

   C_ptr     = (int    *) ML_allocate((N+1)* sizeof(int) );
   Cval = NULL; Ccol = NULL;
   while ( (Cval == NULL) && (total > Bmatrix->max_nz_per_row) ) {
      if (Ccol != NULL) ML_free(Ccol);
      Ccol  = (int    *) ML_allocate( total* sizeof(int) );
      Cval  = (double *) ML_allocate( total* sizeof(double));
      if (Cval == NULL) total = total/2;
   }
   if ( (Cval == NULL) && (N != 0)) {
      printf("Not enough space for new matrix in ML_matmatmult().\n");
      printf("trying to allocate %d elements \n",total);
      printf("Left  matrix has %d rows \n", Amatrix->getrow->Nrows);
      printf("Right matrix has %d rows \n", Bmatrix->getrow->Nrows);
      printf("Left  matrix has %d nz per row\n", Amatrix->max_nz_per_row);
      printf("Right matrix has %d nz per row\n", Bmatrix->max_nz_per_row);
      exit(1);
   }

   next_nz   = 0;
   C_ptr[0]  = next_nz;

   for (i = 0; i < N ; i++) {

      /* Compute a new row in the matrix */

      Ncols = 0;
      fun1(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val, &rowi_N, 0);

      for (k = 0; k < rowi_N; k++) {
         multiplier = rowi_val[k];
         if ( multiplier == 0.0 ) continue; 
         fun2(Bmatrix, 1, &(rowi_col[k]), &accum_size, &accum_col, 
              &accum_val, &row2_N, Ncols);

         d2temp = &(accum_val[Ncols]);
         for (jj = 0; jj < row2_N; jj++) *d2temp++ *= multiplier;
         Ncols = Ncols + row2_N;
      }

      /* lots of different ways to merge rows together */

      if (next_nz+Ncols < total) {
         Tcol = &(Ccol[next_nz]);
         Tval = &(Cval[next_nz]);
         array_flag = 1;
      }
      else { Tcol = accum_col; Tval = accum_val; array_flag = 0;}
      Ncols2 = Ncols;
      Ncols = 0;
      val  = accum_val;
      last = Ncols2;
      if (last > HASH_SIZE) last = HASH_SIZE;
      tflag = 1;
      start = 0;
      while (tflag) {

         for (jj = start; jj < last; jj++) {
            kk = accum_col[jj];
            hash= kk & MASK;
            index = newwhere[hash];
            if ( index == -1) {
               Tcol[Ncols] = kk;
               Tval[Ncols] = *val++;
               newwhere[hash]     = Ncols++;
            }
            else if ( Tcol[index] == kk ) {
               Tval[index] += *val++;
            }
            else {
               hash++;  hash = hash & MASK; index = newwhere[hash];
               while ( (index != -1) && (Tcol[index] != kk) ) {
                  hash++;  hash = hash & MASK; index = newwhere[hash];
               }
               if ( index == -1) { 
                  Tcol[Ncols] = kk;
                  Tval[Ncols] = *val++;
                  newwhere[hash]    = Ncols++;
               }
               else {
                  Tval[index] += *val++;
               }
            }
         }
         tflag = 0;
         if ((last != Ncols2) && (2*Ncols < HASH_SIZE)) {
            tflag = 1;
            start = last;
            last  = Ncols2;
            if (last > HASH_SIZE - Ncols+start) last = HASH_SIZE - Ncols+start;
         }
      }
      if (last != Ncols2) {
         if (array_flag) {
            for (jj = 0; jj < Ncols; jj++) {
               accum_col[jj] = Tcol[jj];
               accum_val[jj] = Tval[jj];
            }
            Tcol = accum_col; 
            Tval = accum_val;
            array_flag = 0;
         }
         for (jj = last; jj < Ncols2; jj++) {
            accum_col[Ncols+jj-last] = accum_col[jj];
            accum_val[Ncols+jj-last] = accum_val[jj];
         }
         Ncols += (Ncols2-last);

         ML_az_sort(accum_col, Ncols, NULL, accum_val);
         ML_sum_duplicates(accum_col, accum_val, &Ncols);
      }
      for (jj = 0; jj < Ncols; jj++) {
         hash= Tcol[jj] & MASK;
         newwhere[hash++] = -1;
         hash = hash & MASK;
         while ( newwhere[hash] != -1) {
            newwhere[hash++] = -1;
            hash = hash & MASK;
         }
      }

      /* empty row. Let's just put a zero in the first column */
      if (Ncols == 0) {
         Tcol[Ncols] = 0;
         Tval[Ncols++] = 0.0;
      }

      /* check if we have enough space to store the new matrix row */
      /* If not, place the current matrix into a submatrix and     */
      /* allocate new vectors (hopefully large enough to hold the  */
      /* rest of the matrix) and a 'parent' matrix to hold the     */
      /* remaining rows.                                           */

      if (next_nz+Ncols > total) {

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
         current->N_nonzeros     = total_nz; 
         current->sub_matrix   = previous_matrix;

         /* allocate space for new matrix */

         if (i != 0) {
            dtemp = ((double) N-i)/ ((double) i);
            dtemp *= (double) total_nz;
            total = (int)(1.1*dtemp);
            total += Ncols;
         }
         else total = total*N + Ncols;

         C_ptr = (int    *) ML_allocate( (N-i+1)* sizeof(int) );
         Ccol  = (int    *) ML_allocate( total* sizeof(int) );
         Cval  = (double *) ML_allocate( total* sizeof(double));
         if (Cval == NULL) {
            printf("Not enough space for matrix\n");
            exit(1);
         }
         next_nz   = 0;
         C_ptr[0]  = next_nz;

         previous_matrix = current;
         sub_i = 0;
      }

      /* store matrix row */

      if (array_flag) next_nz += Ncols;
      else {
         for (k = 0; k < Ncols; k++) {
            Ccol[next_nz] = Tcol[k];
            Cval[next_nz++] = Tval[k];
         }
      }
      C_ptr[sub_i+1] = next_nz;
      sub_i++;
      if (Ncols > max_nz_row_new) max_nz_row_new = Ncols;
   }
   ML_free(accum_col);
   ML_free(accum_val);
   ML_free(rowi_val);
   ML_free(rowi_col);
   ML_free(newwhere);

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
   (*Cmatrix)->N_nonzeros     = total_nz;
   (*Cmatrix)->sub_matrix     = previous_matrix;
   if (allocated-1 > Amatrix->max_nz_per_row) 
      Amatrix->max_nz_per_row = allocated;

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

/********************************************************************/
/********************************************************************
 * Cmatrix = Amatrix * Bmatrix
 *
 * Note: Matrices can be stored in chunks. See ml_rap.h for a description
 *       of the matrix structure ML_matrix.
 * -------------------------------------------------------------------- */

void ML_oldmatmat_mult(ML_Operator *Amatrix, ML_Operator *Bmatrix,
                    ML_Operator **Cmatrix)
{

   int    i,k, jj, next_nz, Ncols, N, total, sub_i, accum_size, row;
   int    *C_ptr, *Ccol, *rowi_col, rowi_N, *accum_col, row2_N;
   double *Cval, *rowi_val, multiplier, *accum_val, *d2temp, dtemp;
   ML_Operator *current, *previous_matrix;
   struct ML_CSR_MSRdata *temp;
   int    max_nz_per_row, max_nz_row_new = 0, total_nz = 0;
   double A_avg_nz_per_row, B_avg_nz_per_row, estimated_nz_per_row;
   int allocated;

   *Cmatrix = NULL;

   /********************************************************/
   /* put here to alleviate a hidden error - C. Tong       */
   if ( Amatrix->max_nz_per_row > Amatrix->getrow->Nrows )
      Amatrix->max_nz_per_row = Amatrix->getrow->Nrows;
   if ( Bmatrix->max_nz_per_row > Bmatrix->getrow->Nrows )
      Bmatrix->max_nz_per_row = Bmatrix->getrow->Nrows;
   /********************************************************/

   previous_matrix = NULL;
   max_nz_per_row  = Bmatrix->max_nz_per_row + 1;
   N      = Amatrix->getrow->Nrows;
   sub_i  = 0;

   accum_size = 10*Bmatrix->max_nz_per_row+ 1;

   /* allocate space to hold rows of Amatrix and the accumulator */

   allocated = Amatrix->max_nz_per_row + 1;
   rowi_col  = (int    *) ML_allocate(allocated * sizeof(int) );
   rowi_val  = (double *) ML_allocate(allocated * sizeof(double));
   accum_col = (int    *) ML_allocate( accum_size * sizeof(int) );
   accum_val = (double *) ML_allocate( accum_size * sizeof(double) );
   if ( (rowi_val == NULL) || (accum_val == NULL)) {
      printf("Not enough space in ML_matmatmult().\n");
      printf("trying to allocate %d %d elements \n",allocated,accum_size);
      printf("Left  matrix has %d rows \n", Amatrix->getrow->Nrows);
      printf("Right matrix has %d rows \n", Bmatrix->getrow->Nrows);
      printf("Left  matrix has %d nz per row\n", Amatrix->max_nz_per_row);
      printf("Right matrix has %d nz per row\n", Bmatrix->max_nz_per_row);
      exit(1);
   }

   /* Make conservative estimates as to the size of the resulting  */
   /* matrix and the size needed for the accumulator.              */
   /* NOTE: These arrays can be increased later in the computation */

   if ( Amatrix->getrow->Nrows > 0 )
   {
      row = 0;
      ML_get_matrix_row(Amatrix,1, &row, &allocated , &rowi_col, &rowi_val, &i,0);
      row = (Amatrix->getrow->Nrows-1)/2;
      ML_get_matrix_row(Amatrix,1, &row, &allocated , &rowi_col, &rowi_val, &k,0);
      row = Amatrix->getrow->Nrows-1;
      ML_get_matrix_row(Amatrix,1, &row, &allocated , &rowi_col, &rowi_val, &jj,0);
      A_avg_nz_per_row = ((double) (i+k+jj))/3.0;
   } else A_avg_nz_per_row = 100;

   if ( Bmatrix->getrow->Nrows > 0 )
   {
      row = 0;
      ML_get_matrix_row(Bmatrix,1,&row, &accum_size, &accum_col, &accum_val, &i,0);
      row = (Bmatrix->getrow->Nrows-1)/2;
      ML_get_matrix_row(Bmatrix,1,&row, &accum_size, &accum_col, &accum_val, &k,0);
      row = Bmatrix->getrow->Nrows-1;
      ML_get_matrix_row(Bmatrix,1,&row, &accum_size, &accum_col, &accum_val,&jj,0);
      B_avg_nz_per_row = ((double) (i+k+jj))/3.0;
   } else B_avg_nz_per_row = i = k = jj = 100;

   if (i  > Bmatrix->max_nz_per_row) Bmatrix->max_nz_per_row =  i+10;
   if (k  > Bmatrix->max_nz_per_row) Bmatrix->max_nz_per_row =  k+10;
   if (jj > Bmatrix->max_nz_per_row) Bmatrix->max_nz_per_row = jj+10;

   estimated_nz_per_row  = sqrt(A_avg_nz_per_row) + sqrt(B_avg_nz_per_row) - 1.;
   estimated_nz_per_row *= estimated_nz_per_row;
   total = (int) (((double) Amatrix->getrow->Nrows)*estimated_nz_per_row * .75) +
           100;

   if (total <= Bmatrix->max_nz_per_row) total = Bmatrix->max_nz_per_row + 1;

   /* Allocate space for the matrix. Since 'total' is just an */
   /* estimate of the space needed, we will reduce 'total' if */
   /* we are unsuccessful allocating space.                   */

   C_ptr     = (int    *) ML_allocate((N+1)* sizeof(int) );
   Cval = NULL; Ccol = NULL;
   while ( (Cval == NULL) && (total > Bmatrix->max_nz_per_row) ) {
      if (Ccol != NULL) ML_free(Ccol);
      Ccol  = (int    *) ML_allocate( total* sizeof(int) );
      Cval  = (double *) ML_allocate( total* sizeof(double));
      if (Cval == NULL) total = total/2;
   }
   if (Cval == NULL) {
      printf("Not enough space for new matrix in ML_matmatmult().\n");
      printf("trying to allocate %d elements \n",total);
      printf("Left  matrix has %d rows \n", Amatrix->getrow->Nrows);
      printf("Right matrix has %d rows \n", Bmatrix->getrow->Nrows);
      printf("Left  matrix has %d nz per row\n", Amatrix->max_nz_per_row);
      printf("Right matrix has %d nz per row\n", Bmatrix->max_nz_per_row);
      exit(1);
   }

   next_nz   = 0;
   C_ptr[0]  = next_nz;

   for (i = 0; i < N ; i++) {

      /* Compute a new row in the matrix */

      Ncols = 0;
      ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, 
                        &rowi_val, &rowi_N, 0);

      for (k = 0; k < rowi_N; k++) {
         multiplier = rowi_val[k];
         if ( multiplier == 0.0 ) continue; /* added by Tong */
         if (Ncols + max_nz_per_row > accum_size) {
            ML_az_sort(accum_col, Ncols, NULL, accum_val);
            ML_sum_duplicates(accum_col, accum_val, &Ncols);
            /*if (Ncols + max_nz_per_row > accum_size/3) { deleted by Tong*/
               accum_size += (Ncols + max_nz_per_row)*20;
               ML_expand_accum(accum_size, &accum_col, &accum_val, Ncols);
            /*} deleted by Tong */
         }

         ML_get_matrix_row(Bmatrix, 1, &(rowi_col[k]), &accum_size, &accum_col,
                           &accum_val, &row2_N, Ncols);

         d2temp = &(accum_val[Ncols]);
         for (jj = 0; jj < row2_N; jj++) d2temp[jj] *= multiplier;
         Ncols = Ncols + row2_N;
      }
      ML_az_sort(accum_col, Ncols, NULL, accum_val);
      ML_sum_duplicates(accum_col, accum_val, &Ncols);

      /* empty row. Let's just put a zero in the first column */
      if (Ncols == 0) {
         accum_col[Ncols] = 0;
         accum_val[Ncols++] = 0.0;
      }

      /* check if we have enough space to store the new matrix row */
      /* If not, place the current matrix into a submatrix and     */
      /* allocate new vectors (hopefully large enough to hold the  */
      /* rest of the matrix) and a 'parent' matrix to hold the     */
      /* remaining rows.                                           */

      if (next_nz+Ncols > total) {

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
         current->N_nonzeros     = total_nz; 
         current->sub_matrix   = previous_matrix;

         /* allocate space for new matrix */

         if (i != 0) {
            dtemp = ((double) N-i)/ ((double) i);
            dtemp *= (double) total_nz;
            total = (int)(1.1*dtemp);
            total += Ncols;
/*
            total = ((int) (1.1 * ((double) ((N-i)* total_nz))/
                                ((double) i))) + Ncols;
*/
         }
         else total = total*N + Ncols;

         C_ptr = (int    *) ML_allocate( (N-i+1)* sizeof(int) );
         Ccol  = (int    *) ML_allocate( total* sizeof(int) );
         Cval  = (double *) ML_allocate( total* sizeof(double));
         if (Cval == NULL) {
            printf("Not enough space for matrix\n");
            exit(1);
         }
         next_nz   = 0;
         C_ptr[0]  = next_nz;

         previous_matrix = current;
         sub_i = 0;
      }

      /* store matrix row */

      for (k = 0; k < Ncols; k++) {
         Ccol[next_nz] = accum_col[k];
         Cval[next_nz++] = accum_val[k];
      }
      C_ptr[sub_i+1] = next_nz;
      sub_i++;
      if (Ncols > max_nz_row_new) max_nz_row_new = Ncols;
   }
   ML_free(accum_col);
   ML_free(accum_val);
   ML_free(rowi_val);
   ML_free(rowi_col);

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
   (*Cmatrix)->N_nonzeros     = total_nz;
   (*Cmatrix)->sub_matrix     = previous_matrix;
   if (allocated-1 > Amatrix->max_nz_per_row) 
      Amatrix->max_nz_per_row = allocated;

}

/* ******************************************************************** */
/* multiplying two matrices together                                    */
/* -------------------------------------------------------------------- */

void ML_2matmult(ML_Operator *Mat1, ML_Operator *Mat2,
                 ML_Operator *Result, int matrix_type)
{
   int N_input_vector, max_per_proc;
   ML_CommInfoOP *getrow_comm;
   ML_Operator   *Mat2comm, *Mat1Mat2, *tptr;
   ML_Comm       *comm;

   comm = Mat1->comm;
   N_input_vector = Mat2->invec_leng;
   getrow_comm    = Mat2->getrow->pre_comm;

   ML_create_unique_col_id(N_input_vector, &(Mat2->getrow->loc_glob_map),
                           getrow_comm, &max_per_proc, comm);
   Mat2->getrow->use_loc_glob_map = ML_YES;

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

   if (matrix_type == ML_CSR_MATRIX)
     ML_back_to_csrlocal(Mat1Mat2, Result, max_per_proc);
   else if (matrix_type = ML_MSR_MATRIX) {
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

   ML_RECUR_CSR_MSRdata_Destroy(Mat1Mat2);
   ML_Operator_Destroy(&Mat1Mat2);
}

#endif
