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
/* functions defined outside this file                                  */
/* ******************************************************************** */

extern void ML_get_matrow_CSR(ML_Operator *input_matrix, int N_requested_rows,
        int requested_rows[], int *allocated_space, int **columns,
        double **values, int row_lengths[], int index);
extern void ML_get_row_CSR_norow_map(ML_Operator *input_matrix, 
        int N_requested_rows, int requested_rows[], int *allocated_space, 
        int **columns, double **values, int row_lengths[], int index);

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
      if (next->getrow->external != CSR_getrows) flag = 0;
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
      if (next->getrow->external != CSR_getrows) flag = 0;
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
				       Amatrix->outvec_leng, ML_EMPTY,temp,
				       i,NULL,0);

         ML_Operator_Set_Getrow(current, ML_EXTERNAL, i, CSR_getrows);

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
			        Amatrix->outvec_leng,ML_EMPTY,temp,N,NULL,0);
   ML_Operator_Set_Getrow(*Cmatrix, ML_EXTERNAL, N, CSR_getrows);
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
				       Amatrix->outvec_leng, ML_EMPTY,temp,
				       i,NULL,0);

         ML_Operator_Set_Getrow(current, ML_EXTERNAL, i, CSR_getrows);

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
			        Amatrix->outvec_leng,ML_EMPTY,temp,N,NULL,0);
   ML_Operator_Set_Getrow(*Cmatrix, ML_EXTERNAL, N, CSR_getrows);
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
                 ML_Operator *Result)
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
      ML_Operator_Destroy(Mat2comm);
   }

   ML_back_to_csrlocal(Mat1Mat2, Result, max_per_proc);

   ML_RECUR_CSR_MSRdata_Destroy(Mat1Mat2);
   ML_Operator_Destroy(Mat1Mat2);
}

