/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* ******************************************************************** */
/* Functions for manipulating matrices                                  */
/* ******************************************************************** */
/* ******************************************************************** */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ml_struct.h"

/*********************************************************************/
/* Perform a matrix-vector product:  ovec = matrix * vec             */
/*                                                                   */
/* Parameters                                                        */
/* ==========                                                        */
/*   matrix        On input, a sparse matrix.                        */
/*   vec           On input, a vector of length Nvec.                */
/*   Nvec          On input, length of 'vec'.                        */
/*   ovec          On output, ovec = matrix * vec                    */
/*   Novec         On output, length of 'ovec'.                      */
/*********************************************************************/

void ML_getrow_matvec(ML_Operator *matrix, double *vec, int Nvec, 
                      double *ovec, int *Novec)
{
   ML_Operator *temp, *temp2, *temp3, *temp4, *tptr;
   int *cols, i;
   int allocated, row_length;

   if (matrix->getrow->func_ptr == NULL) {
      printf("ML_getrow_matvec: empty object? \n");
      exit(1);
   }
   temp = ML_Operator_Create(matrix->comm);
   ML_Operator_Set_1Levels(temp, matrix->from, matrix->from);
   ML_Operator_Set_ApplyFuncData(temp,1,Nvec,vec,Nvec,NULL,0);

   ML_Operator_Set_Getrow(temp,Nvec, VECTOR_getrows);
   temp->max_nz_per_row = 1;
   temp->N_nonzeros     = Nvec;

   if (matrix->getrow->pre_comm != NULL) {
      ML_exchange_rows(temp, &temp2, matrix->getrow->pre_comm);
   }
   else temp2 = temp;

   ML_matmat_mult(matrix, temp2, &temp3);

   if (matrix->getrow->post_comm != NULL)
      ML_exchange_rows(temp3, &temp4, matrix->getrow->post_comm);
   else temp4 = temp3;

   allocated = temp4->getrow->Nrows + 1;
   cols = (int *) ML_allocate(allocated*sizeof(int));
   if (cols == NULL) {
      printf("no space in ML_getrow_matvec()\n");
      exit(1);
   }
   for (i = 0; i < temp4->getrow->Nrows; i++) {
      ML_get_matrix_row(temp4, 1, &i, &allocated , &cols, &ovec,
                   &row_length, i);
      if (allocated != temp4->getrow->Nrows + 1)
         printf("memory problems ... we can't reallocate here\n");
   }

   ML_free(cols);

   if ( *Novec != temp4->getrow->Nrows) {
     printf("Warning: The length of ML's output vector does not agree with\n");
     printf("         the user's length for the output vector (%d vs. %d).\n",
            *Novec, temp4->getrow->Nrows);
     printf("         indicate a problem.\n");
   }
   *Novec = temp4->getrow->Nrows;

   if (matrix->getrow->pre_comm != NULL) {
      tptr = temp2;
      while ( (tptr!= NULL) && (tptr->sub_matrix != temp))
         tptr = tptr->sub_matrix;
      if (tptr != NULL) tptr->sub_matrix = NULL;
      ML_RECUR_CSR_MSRdata_Destroy(temp2);
      ML_Operator_Destroy(&temp2);
   }
   if (matrix->getrow->post_comm != NULL) {
      tptr = temp4;
      while ( (tptr!= NULL) && (tptr->sub_matrix != temp3))
         tptr = tptr->sub_matrix;
      if (tptr != NULL) tptr->sub_matrix = NULL;
      ML_RECUR_CSR_MSRdata_Destroy(temp4);
      ML_Operator_Destroy(&temp4);
   }

   ML_Operator_Destroy(&temp);
   ML_RECUR_CSR_MSRdata_Destroy(temp3);
   ML_Operator_Destroy(&temp3);
}

/*********************************************************************/
/* check the correctness of RAP construction                         */
/*********************************************************************/

void ML_rap_check(ML *ml, ML_Operator *RAP, ML_Operator *R,
        ML_Operator *A, ML_Operator *P, int iNvec, int oNvec)
{
   int i,j;
   double *vec1, *vec2, *vec3, *vec4, *vec5;
   double norm1, norm2;

#ifdef DEBUG
   printf("ML_rap_check begins ...\n");
#endif

   if (RAP->getrow->ML_id != ML_ID_MATRIX) {
      if (ml->comm->ML_mypid == 0)
         printf("ML_rap_check: RAP is the wrong object (=%d). \n",
                         RAP->getrow->ML_id);
      exit(1);
   }
   if (R->getrow->ML_id != ML_ID_MATRIX) {
      if (ml->comm->ML_mypid == 0)
        printf("ML_rap_check: R is the wrong object (=%d). \n",
                         RAP->getrow->ML_id);
      exit(1);
   }
   if (P->getrow->ML_id != ML_ID_MATRIX) {
      if (ml->comm->ML_mypid == 0)
        printf("ML_rap_check: P is the wrong object (=%d). \n",
                         RAP->getrow->ML_id);
      exit(1);
   }
   if (A->getrow->ML_id != ML_ID_MATRIX) {
      if (ml->comm->ML_mypid == 0)
        printf("ML_rap_check: A is the wrong object (=%d). \n",
                         RAP->getrow->ML_id);
      exit(1);
   }

   /* j is the number of external variables */

   j = 0;
   for (i = 0; i < RAP->getrow->pre_comm->N_neighbors; i++)
      j += RAP->getrow->pre_comm->neighbors[i].N_rcv;

   vec1 = (double *) ML_allocate((iNvec+ 1 + j)*sizeof(double));
   vec2 = (double *) ML_allocate((P->getrow->Nrows+1)*sizeof(double));
   vec3 = (double *) ML_allocate((A->getrow->Nrows+1)*sizeof(double));
   vec4 = (double *) ML_allocate((oNvec + 1)*sizeof(double));
   vec5 = (double *) ML_allocate((oNvec + 1)*sizeof(double));

   for (i = 0; i < iNvec; i++)
      vec1[i] = (double) (ml->comm->ML_mypid*2301 + i*7 + 1);

   j = P->getrow->Nrows;
   ML_getrow_matvec(P, vec1, iNvec, vec2,&j);
   i = A->getrow->Nrows;
   ML_getrow_matvec(A, vec2, j, vec3,&i);
   ML_getrow_matvec(R, vec3, i, vec4,&oNvec);

   /* j is the number of variables sent in communication */

   j = 0;
   for (i = 0; i < RAP->getrow->pre_comm->N_neighbors; i++)
      j += RAP->getrow->pre_comm->neighbors[i].N_send;

   ML_restricted_MSR_mult( RAP, oNvec, vec1, vec5, j);

   norm1 = sqrt(ML_gdot(oNvec, vec5, vec5, ml->comm));
   for (i = 0; i < oNvec; i++) vec5[i] -= vec4[i];
   norm2 = sqrt(ML_gdot(oNvec, vec5, vec5, ml->comm));

   if (norm2 > norm1*1e-10) {
      norm2 = sqrt(ML_gdot(oNvec, vec4, vec4, ml->comm));
      if (ml->comm->ML_mypid == 0) {
          printf("***************************************\n");
          printf("RAP seems inaccurate:\n");
          printf("    ||    RAP v    ||_2 = %e\n\n", norm1);
          printf("    || R (A (P v)) ||_2 = %e\n",norm2);
          printf("***************************************\n");
      }
   }

   ML_free(vec5);
   ML_free(vec4);
   ML_free(vec3);
   ML_free(vec2);
   ML_free(vec1);

#ifdef DEBUG
   printf("ML_rap_check ends ...\n");
#endif
}

/*********************************************************************/
/* Get matrix rows (given in requested_rows[]) from the user's matrix*/
/* and return the information in (*columns)[index ...] and           */
/* (*values)[index ...].                                             */
/*                                                                   */
/* Parameters                                                        */
/* ==========                                                        */
/*   input_matrix      On input, data structure representing the     */
/*                     matrix whose rows will be returned            */
/*                     (see ml_rap.h).                               */
/*                                                                   */
/*   N_requested_rows  On input, number of rows that will are        */
/*                     requested.  Note: current implementation      */
/*                     allows for only 1 row.                        */
/*                                                                   */
/*   requested_rows    On input, requested_rows[i]                   */
/*                     ( 0 <= i < N_requested_rows) is the local row */
/*                     number of the ith requested row.              */
/*                                                                   */
/*   allocated_space   On input, length of *columns and *values. On  */
/*                     output, 'allocated_space' is the length of    */
/*                     *columns and *values on exit. NOTE: *columns  */
/*                     and *values may be reallocated within the     */
/*                     subroutine if there is not enough space.      */
/*                                                                   */
/*   columns      On input, an array of length 'allocated_space'. On */ 
/*                output, column numbers of the nonzero in retreived */
/*                rows are stored starting from (*columns)[index]    */
/*                (first row stored first). Note: if there is not    */
/*                enough available space, *columns is reallocated    */
/*                (and the contents of the first 'index' variables   */
/*                are copied to the new array.                       */
/*                                                                   */
/*   values       On input, an array of length 'allocated_space'. On */
/*                output, nonzero values in retreived rows           */
/*                are stored starting from (*values)[index] (first   */
/*                row stored first). Note: if there is not enough    */
/*                available space, *values is reallocated (and the   */
/*                contents of the first 'index' variables are copied */
/*                to the new array.                                  */
/*                                                                   */
/*   row_lengths[]  On input, an array of length 'N_requested_rows'. */
/*                  On output, row_lengths[i] is the of nonzeros in  */
/*                  'requested_row[i]'.                              */
/*                                                                   */
/*********************************************************************/

void ML_get_matrix_row(ML_Operator *input_matrix, int N_requested_rows,
        int requested_rows[], int *allocated_space, int **columns,
        double **values, int row_lengths[], int index)
{
   int    i, *mapper, *t1, row;
   ML_Operator *next;
   double *t2;
   void *data;
   int (*getfunction)(void *,int,int*,int,int*,double*,int*);

#ifdef DEBUG2
   if (N_requested_rows != 1) {
      printf("ML_get_matrix_row is currently implemented for only 1 row");
      printf(" at a time.\n");
      exit(1);
   }
#endif

   row = requested_rows[0]; 
#ifdef DEBUG2
   if ( (row >= input_matrix->getrow->Nrows) || (row < 0) ) {
      row_lengths[0] = 0;
      return;
   }
#endif

   if (input_matrix->getrow->row_map != NULL) {
      if (input_matrix->getrow->row_map[row] != -1) 
         row = input_matrix->getrow->row_map[row];
      else { 
	row_lengths[0] = 0; 
	ML_avoid_unused_param( (void *) &N_requested_rows);
	return;}
   }

   next = input_matrix->sub_matrix;

   while ( (next != NULL) && (row < next->getrow->Nrows) ) {
      input_matrix = next;
      next = next->sub_matrix;
   }
   if (next != NULL) row -= next->getrow->Nrows;

   data = (void *) input_matrix;
   getfunction = (int (*)(void *,int,int*,int,int*,double*,int*))
     input_matrix->getrow->func_ptr;

   while(getfunction(data,1,&row,*allocated_space-index,
               &((*columns)[index]), &((*values)[index]), row_lengths) == 0) {
      *allocated_space = 2*(*allocated_space) + 1;
      t1 = (int    *) ML_allocate(*allocated_space*sizeof(int   ));
      t2 = (double *) ML_allocate(*allocated_space*sizeof(double));
      if (t2 == NULL) {
         printf("Not enough space to get a matrix row. A row length of \n");
         printf("%d was not sufficient\n",(*allocated_space-1)/2);
	 fflush(stdout);
         exit(1);
      }
      for (i = 0; i < index; i++) t1[i] = (*columns)[i];
      for (i = 0; i < index; i++) t2[i] = (*values)[i];
      if (*columns != NULL) ML_free(*columns);  
      if (*values  != NULL) ML_free(*values);
      *columns = t1;
      *values  = t2;
   }

   if ( (input_matrix->getrow->use_loc_glob_map == ML_YES)) {
      mapper       = input_matrix->getrow->loc_glob_map;
      for (i = 0; i < row_lengths[0]; i++) 
         (*columns)[i+index] = mapper[(*columns)[index+i]];
   }
}

void ML_get_matrow_CSR(ML_Operator *input_matrix, int N_requested_rows,
        int requested_rows[], int *allocated_space, int **columns,
        double **values, int row_lengths[], int index)
{
   int    i, *mapper, *t1, row;
   ML_Operator *next;
   double *t2;
   struct ML_CSR_MSRdata *matrix;
   int    *rowptr, *bindx, *col_ptr, itemp, j;
   double *val, *val_ptr;


#ifdef DEBUG2
   if (N_requested_rows != 1) {
      printf("ML_get_matrix_row is currently implemented for only 1 row");
      printf(" at a time.\n");
      exit(1);
   }
#endif

   row = requested_rows[0]; 
#ifdef DEBUG2
   if ( (row >= input_matrix->getrow->Nrows) || (row < 0) ) {
      row_lengths[0] = 0;
      return;
   }
#endif

   if (input_matrix->getrow->row_map != NULL) {
      if (input_matrix->getrow->row_map[row] != -1) 
         row = input_matrix->getrow->row_map[row];
      else { row_lengths[0] = 0; 
	ML_avoid_unused_param( (void *) &N_requested_rows);
	return;}
   }

   next = input_matrix->sub_matrix;
   while ( (next != NULL) && (row < next->getrow->Nrows) ) {
      input_matrix = next;
      next = next->sub_matrix;
   }
   if (next != NULL) row -= next->getrow->Nrows;

   matrix = (struct ML_CSR_MSRdata *) input_matrix->data;
   rowptr = matrix->rowptr;
   itemp   = rowptr[row];
   bindx   = &(matrix->columns[itemp]);
   val     = &(matrix->values[itemp]);

   *row_lengths = rowptr[row+1] - itemp;

   if (*row_lengths+index > *allocated_space) {
      *allocated_space = 2*(*allocated_space) + 1;
      if (*row_lengths+index > *allocated_space) 
         *allocated_space = *row_lengths + 5 + index;
      t1 = (int    *) ML_allocate(*allocated_space*sizeof(int   ));
      t2 = (double *) ML_allocate(*allocated_space*sizeof(double));
      if (t2 == NULL) {
         printf("Not enough space to get a matrix row. A row length of \n");
         printf("%d was not sufficient\n",(*allocated_space-1)/2);
	 fflush(stdout);
         exit(1);
      }
      for (i = 0; i < index; i++) t1[i] = (*columns)[i];
      for (i = 0; i < index; i++) t2[i] = (*values)[i];
      ML_free(*columns);  ML_free(*values);
      *columns = t1;
      *values  = t2;
   }
   col_ptr = &((*columns)[index]);
   val_ptr = &((*values)[index]);

   for (j = 0 ; j < *row_lengths; j++) {
      *col_ptr++ = *bindx++;
   }
   for (j = 0 ; j < *row_lengths; j++) {
      *val_ptr++  = *val++;
   }

   if ( (input_matrix->getrow->use_loc_glob_map == ML_YES)) {
      mapper       = input_matrix->getrow->loc_glob_map;
      for (i = 0; i < row_lengths[0]; i++) 
         (*columns)[i+index] = mapper[(*columns)[index+i]];
   }
}
void ML_get_row_CSR_norow_map(ML_Operator *input_matrix, int N_requested_rows,
        int requested_rows[], int *allocated_space, int **columns,
        double **values, int row_lengths[], int index)
{
   int    i, *mapper, *t1, row;
   ML_Operator *next;
   double *t2;
   struct ML_CSR_MSRdata *matrix;
   int    *rowptr, *bindx, *col_ptr, itemp, j;
   double *val, *val_ptr;


#ifdef DEBUG2
   if (N_requested_rows != 1) {
      printf("ML_get_matrix_row is currently implemented for only 1 row");
      printf(" at a time.\n");
      exit(1);
   }
#endif

   row = requested_rows[0]; 
#ifdef DEBUG2
   if ( (row >= input_matrix->getrow->Nrows) || (row < 0) ) {
      row_lengths[0] = 0;
      return;
   }
#endif

   next = input_matrix->sub_matrix;
   while ( (next != NULL) && (row < next->getrow->Nrows) ) {
      input_matrix = next;
      next = next->sub_matrix;
   }
   if (next != NULL) row -= next->getrow->Nrows;

   matrix = (struct ML_CSR_MSRdata *) input_matrix->data;
   rowptr = matrix->rowptr;
   itemp   = rowptr[row];
   bindx   = &(matrix->columns[itemp]);
   val     = &(matrix->values[itemp]);

   *row_lengths = rowptr[row+1] - itemp;

   if (*row_lengths+index > *allocated_space) {
      *allocated_space = 2*(*allocated_space) + 1;
      if (*row_lengths+index > *allocated_space) 
         *allocated_space = *row_lengths + 5 + index;
      t1 = (int    *) ML_allocate(*allocated_space*sizeof(int   ));
      t2 = (double *) ML_allocate(*allocated_space*sizeof(double));
      if (t2 == NULL) {
         printf("Not enough space to get a matrix row. A row length of \n");
         printf("%d was not sufficient\n",(*allocated_space-1)/2);
	 fflush(stdout);
	 ML_avoid_unused_param( (void *) &N_requested_rows);
         exit(1);
      }
      for (i = 0; i < index; i++) t1[i] = (*columns)[i];
      for (i = 0; i < index; i++) t2[i] = (*values)[i];
      ML_free(*columns);  ML_free(*values);
      *columns = t1;
      *values  = t2;
   }
   col_ptr = &((*columns)[index]);
   val_ptr = &((*values)[index]);

   for (j = 0 ; j < *row_lengths; j++) {
      *col_ptr++ = *bindx++;
   }
   for (j = 0 ; j < *row_lengths; j++) {
      *val_ptr++  = *val++;
   }

   if ( (input_matrix->getrow->use_loc_glob_map == ML_YES)) {
      mapper       = input_matrix->getrow->loc_glob_map;
      for (i = 0; i < row_lengths[0]; i++) 
         (*columns)[i+index] = mapper[(*columns)[index+i]];
   }
}
