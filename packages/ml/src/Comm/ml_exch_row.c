/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "ml_struct.h"

/******************************************************************************/
/******************************************************************************
  Routine to exchange rows between neighboring processors. 
  Note: It is useful that matrix columns have a unique global id (that
        is column numbers are not local columns but instead unique global
        values that can be deciphered on other processors).
  The rows which are exchanged are given by comm_info. The exchange is done in 
  the following steps: 

  1) Count number of matrix rows that will be send and received by this node.
  2) Compute row lengths of rows to be sent. Communicate the length 
     information so processors know how big the rows they will receive are.
  3) Send column information for the new rows. Append received rows onto
     the end of the matrix.
  4) Send row information for the new rows. Append received rows onto the 
     end of the matrix.
  5) Create an ML_matrix structure and fill it so that it corresponds to
     the appended matrix (with the newly received rows).
  6) The appended matrix is now modified depending on whether a rcv_list
     is given or not. There are essentially three possible modifications:

        i) rcv_list given and received values over-write existing values

           We effectively move the appended rows into the matrix. This is done 
           by creating a 'map' array which indicates that the ith row is stored
           in map[i]. In this way, the appended matrix is not really altered 
           internally. However, ML_matrix_get_row() will know how to correctly 
           access the ith row. 

       ii) rows can be remapped. That is, row[k] can now be referred to by 
           row[j]. 

           We effectively move rows around. This is done in a fashion similar
           to 6i) using the array 'remap' passed in by the user.  
           Note: remap[i] = -1 indictes that the ith row is empty.
           
      iii) rcv_list is given and received values add with existing values.

           We effectively take the appended rows and add them in with the
           original rows according to the receive list. This is done by 
           creating a new little matrix and appending it to the original
           matrix given by the user. The new little matrix contains the
           sums corresponding to all the received rows. Once these sums are
           created the received rows (obtained in steps 3-5 are removed).
           Finally, a map array is created so that the correct row is referenced
           by ML_matrix_get_row(). 
     
  Author:          Ray S. Tuminaro, SNL, 9222
  =======

  Return code:     void
  ============

  Parameter list:
  ===============
  Pmatrix          On input, matrix whose rows will be communicated.

  Pappended        On output, matrix with rows received and possibly
                   remapped according to comm_info.

  comm_info        On input, structure (set using ML_CommInfoOP_Set_neighbors()
                   and ML_CommInfoOP_Set_exch_info()) containing communication 
                   info.

*******************************************************************************/

#define ML_MPI_MSG_NUM 2391

void ML_exchange_rows(ML_Operator *Pmatrix, ML_Operator **Pappended, 
                      ML_CommInfoOP *comm_info)
{

  /* local variables */

  int         *actual_send_length, *actual_recv_length, *start_send_proc;
  int         Nneighbors, *neighbor, *remap;
  int         Nrows_new, *cols_new, *rowptr_new, Nrows_send, Nrows, max_per_row;
  int         *ibuff, total_num_recv, total_send, total_recv;
  double      *vals_new, *dtemp, *dbuff, *dummy1;
  int         i, j, k, ii, jj, *newmap, *orig_map, nonNULL_rcv_list, *dummy2;
  static int  type = ML_MPI_MSG_NUM;
  struct      ML_CSR_MSRdata *temp;
  int         allocated_space, row_length;
  ML_Comm     *comm;
  int Nghost;
  int rcv_list_exists = 0, count = 0;

  /**************************** execution begins ****************************/

  comm        = Pmatrix->comm;
  Nrows       = Pmatrix->getrow->Nrows;
  Nneighbors  = comm_info->N_neighbors; 
  if (Pmatrix->getrow->pre_comm != NULL) {
    ML_CommInfoOP_Compute_TotalRcvLength(Pmatrix->getrow->pre_comm);

    if (Pmatrix->N_total_cols_est < Pmatrix->invec_leng + 
	Pmatrix->getrow->pre_comm->total_rcv_length) {
      Pmatrix->N_total_cols_est = Pmatrix->invec_leng + 
	Pmatrix->getrow->pre_comm->total_rcv_length;
    }
  }



  /* compute the total number of rows to send and receive */

  Nrows_new  = 0;
  Nrows_send = 0;
  Nghost = 0;
  nonNULL_rcv_list = 0;
  for (i = 0; i < Nneighbors; i++) 
  {
     Nrows_send += comm_info->neighbors[i].N_send; 
     Nrows_new  += comm_info->neighbors[i].N_rcv;
     if ((comm_info->neighbors[i].N_rcv != 0) && 
         (comm_info->neighbors[i].rcv_list != NULL)) nonNULL_rcv_list = 1;
     if (comm_info->neighbors[i].rcv_list != NULL) 
     {
        rcv_list_exists = 1;
        for (j = 0; j < comm_info->neighbors[i].N_rcv; j++) {
           if (comm_info->neighbors[i].rcv_list[j] > Nghost + Nrows - 1)
              Nghost = comm_info->neighbors[i].rcv_list[j] - Nrows + 1;
        }
     }
  }
  if (Nghost == 0) Nghost = Nrows_new;

  actual_recv_length = (int *) ML_allocate( (Nneighbors+1)*sizeof(int));
  actual_send_length = (int *) ML_allocate( (Nneighbors+1)*sizeof(int));
  start_send_proc    = (int *) ML_allocate( (Nneighbors+1)*sizeof(int));
  rowptr_new         = (int *) ML_allocate( (Nrows_new+1)*sizeof(int));
  if (rowptr_new == NULL) {
     fprintf(stderr,"Out of space in ML_exchange_rows\n");
     exit(1);
  }

  /* compute the row lengths of the external rows    */
  /* and the total number of nonzeros to be received */
  /* from each processor.                            */

  allocated_space = Pmatrix->max_nz_per_row+2;
  dummy1 = (double *) ML_allocate(allocated_space*sizeof(double));
  dummy2 = (int    *) ML_allocate(allocated_space*sizeof(   int));
  dtemp  = (double *) ML_allocate((Nrows+Nghost + 1)*sizeof(double));
  if (dtemp == NULL) 
  {
     printf("out of space in exch_row\n");  exit(1);
  }
   
  for (j = 0; j < Nrows+Nghost + 1; j++) dtemp[j] = -1.;
  total_send = 0;
  for (i = 0; i < Nneighbors; i++) 
  {
     for (ii = 0; ii < comm_info->neighbors[i].N_send; ii++ ) 
     {
        j = (comm_info->neighbors[i]).send_list[ii];
        if (j > Nrows) 
        {
           printf("Error: the %dth element sent to %d is greater than the\n",
                  ii+1,comm_info->neighbors[i].ML_id);
           printf("Error: total number of local elements: %d vs. %d\n",
                   j, Nrows);
           exit(1);
        }
        ML_get_matrix_row(Pmatrix, 1, &j, &allocated_space, &dummy2, 
                          &dummy1, &row_length, 0);

        dtemp[j] = (double) row_length;
        total_send += row_length;
     }
  }
  ML_cheap_exchange_bdry(dtemp, comm_info, Nrows, Nrows_send, comm);

/*
   ML_exchange_bdry(dtemp, comm_info, Nrows, comm, ML_OVERWRITE);
Unfortunately ML_exchange_bdry() has a problem when two processors
send to the same value `ML_ADD'. I'm not really sure if Kenneth Massey's
example (with nonutilized ghost variables still works
*/

  for (i = 0; i < Nneighbors; i++) 
  {
     for (ii = 0; ii < comm_info->neighbors[i].N_send; ii++ ) 
     {
       j = (comm_info->neighbors[i]).send_list[ii];
       dtemp[j] = -1.;
     }
  }

  j = Nrows;
/*
  j = 0;
*/
  for (i = 0 ; i < Nrows_new ; i++ ) 
  {
     while (dtemp[j] == -1.) j++;
     rowptr_new[i] = (int) dtemp[j];
     j++;
  }

  /* Set rowptr_new[i] so that it points to the beginning of row[i] */
  /* and compute maximum row length in the appended part of matrix  */

  j = rowptr_new[0];
  rowptr_new[0] = 0;
  if (allocated_space-2nz_per_row+ 1;
   row_map         = (int    *) ML_allocate(matrix->getrow->Nrows*sizeof(int) );
   accum_col       = (int    *) ML_allocate( accum_size * sizeof(int) );
   accum_val       = (double *) ML_allocate( accum_size * sizeof(double) );
   for (i = 0; i <  matrix->getrow->Nrows; i++) row_map[i] = -1;

   if (matrix->getrow->row_map != NULL) {
      for (i = 0; i <  matrix->getrow->Nrows; i++) {
         if (matrix->getrow->row_map[i] < orig_rows) 
            row_map[i] = matrix->getrow->row_map[i];
      }
   }
   else for (i = 0; i <  orig_rows; i++) row_map[i] = i;
         
   /* Allocate space for the matrix. Since 'total' is just an */
   /* estimate of the space needed, we will reduce 'total' if */
   /* we are unsuccessful allocating space.                   */

   C_ptr = (int    *) ML_allocate((N_append_rows+1)* sizeof(int) );
   Cval  = NULL; Ccol = NULL;
   total = appended_nzs;
   if (total <= matrix->max_nz_per_row) total = matrix->max_nz_per_row + 2;
   while ( (Cval == NULL) && (total > matrix->max_nz_per_row) ) {
      if (Ccol != NULL) ML_free(Ccol);
      Ccol  = (int    *) ML_allocate( total* sizeof(int) );
      Cval  = (double *) ML_allocate( total* sizeof(double));
      if (Cval == NULL) total = total/2;
   }
   if (Cval == NULL) {
      printf("Not enough space for new matrix in ML_add_appended_rows().\n");
      exit(1);
   }

   previous_matrix = matrix->sub_matrix;
   next_nz   = 0;
   C_ptr[0]  = next_nz;
   i = 0;
   row_count = orig_rows;
   for (ii = 0; ii < comm_info->N_neighbors; ii++) {
      for (jj = 0; jj < comm_info->neighbors[ii].N_rcv; jj++) {
         Ncols     = 0;

         /* determine row number of jjth received row from iith neighbor */

         row_location = comm_info->neighbors[ii].rcv_list[jj];

         if (row_location >= 0) { /* rows already processed are  */
                                  /* encoded as negative numbers */

            row_map[row_location] = row_count;

            /* get matrix row in the non-appended portion of matrix */

            ML_get_matrix_row(matrix, 1, &row_location, &accum_size, &accum_col,
                              &accum_val, &row_length, Ncols);
            Ncols += row_length;

            /* get matrix row in appended portion of matrix */


            itemp = matrix->getrow->row_map;  matrix->getrow->row_map = NULL; 
                                          /* kludge: we don't want to use    */
                                          /* row_map to get at appended rows */
            row_request = i+orig_rows;
            ML_get_matrix_row(matrix, 1, &row_request, &accum_size, &accum_col,
                              &accum_val, &row_length, Ncols);
            Ncols += row_length;
            comm_info->neighbors[ii].rcv_list[jj] = -1 - row_location;
            N_changed += t_changed;
            t_changed = 1;

            /* find any other rows in the appended portion of matrix that */
            /* are supposed to be added into to the row 'row_location'.   */

            start = jj + 1;
            i2    =  i + 1;
            for (iii = ii; iii < comm_info->N_neighbors; iii++) {
               for (jjj = start; jjj < comm_info->neighbors[iii].N_rcv; jjj++) {
                  new_row = comm_info->neighbors[iii].rcv_list[jjj];
                  if (new_row == row_location) {
                     comm_info->neighbors[iii].rcv_list[jjj]= -1 - new_row;
                     t_changed++;

                     if (Ncols + max_nz_per_row > accum_size) {
                        ML_az_sort(accum_col, Ncols, NULL, accum_val);
                        ML_sum_duplicates(accum_col, accum_val, &Ncols);
                        if (Ncols + max_nz_per_row > accum_size/2) {
                           accum_size += (Ncols + max_nz_per_row)*2;
                           ML_expand_accum(accum_size, &accum_col, &accum_val, 
                                           Ncols);
                        }
                     }
                     row_request = i2+orig_rows;
                     ML_get_matrix_row(matrix, 1, &row_request, &accum_size,
                                    &accum_col, &accum_val, &row_length, Ncols);
                     Ncols += row_length;
                  }
                  i2++;
               }
               start = 0;
            }

            /* store accumulator in new matrix */

            ML_az_sort(accum_col, Ncols, NULL, accum_val);
            ML_sum_duplicates(accum_col, accum_val, &Ncols);
            matrix->getrow->row_map = itemp;

            if (next_nz+Ncols > total) { /* not enough space in existing    */
                                         /* chunk of memory. Must close off */
                                         /* existing chunk and allocate a   */
                                         /* new chunk of memory.            */

               /* create sub_matrix object */

               total_nz += next_nz;
               temp = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct 
                                                            ML_CSR_MSRdata));
               temp->columns          = Ccol;
               temp->values           = Cval;
               temp->rowptr           = C_ptr;
               current = ML_Operator_Create(matrix->comm);
               ML_Operator_Set_1Levels(current, parent->from, parent->to);
               ML_Operator_Set_ApplyFuncData(current,row_count,row_count,
                                   ML_EMPTY,(void*)temp,row_count,NULL,0);
               ML_Operator_Set_Getrow(current, ML_EXTERNAL, row_count, 
				      CSR_getrows);
               current->max_nz_per_row = max_nz_row_new;
               if (matrix->sub_matrix->N_nonzeros >= 0)
                  current->N_nonzeros =total_nz+matrix->sub_matrix->N_nonzeros;
               current->sub_matrix = previous_matrix;

               /* allocate space for new matrix */

               if (row_count != orig_rows) {
		 dtemp =   1.1 * ((double) (total_rcvd-N_changed))/
		   ((double) N_changed);
		 total = ((int) ( ((double) total_nz) * dtemp)) + Ncols;
               }
               else total = total*total_rcvd + Ncols;

               C_ptr = (int    *) ML_allocate((N_append_rows-row_count+
                                                  orig_rows+1)*sizeof(int));
               Ccol  = (int    *) ML_allocate(total * sizeof(int) );
               Cval  = (double *) ML_allocate(total * sizeof(double));
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
            row_count++;
         }
         i++;
      }
   }
   ML_free(accum_col);
   ML_free(accum_val);

   /* restore the receive list (i.e. decode) */

   for (ii = 0; ii < comm_info->N_neighbors; ii++) {
      for (jj = 0; jj < comm_info->neighbors[ii].N_rcv; jj++) {
         comm_info->neighbors[ii].rcv_list[jj]= -1 -
                                        comm_info->neighbors[ii].rcv_list[jj];
      }
   }

   /* put together the final pieces */

   /* free the old appended matrix */

   temp = (struct ML_CSR_MSRdata *) matrix->data;
   if (temp->columns  != NULL) ML_free(temp->columns);
   if (temp->values   != NULL) ML_free(temp->values);
   if (temp->rowptr   != NULL) ML_free(temp->rowptr);

   total_nz += next_nz;
   matrix->getrow->Nrows  = row_count;
   matrix->getrow->external = CSR_getrows;
   matrix->getrow->loc_glob_map   = NULL;
   matrix->getrow->use_loc_glob_map = ML_NO;
   if (matrix->getrow->row_map != NULL) ML_free(matrix->getrow->row_map);
   matrix->getrow->row_map        = row_map;
   if (matrix->sub_matrix->N_nonzeros >= 0)
     matrix->N_nonzeros     = total_nz + matrix->sub_matrix->N_nonzeros;
   matrix->sub_matrix     = previous_matrix;
   matrix->max_nz_per_row = max_nz_row_new;
   temp = (struct ML_CSR_MSRdata *) matrix->data;
   temp->columns          = Ccol;
   temp->values           = Cval;
   temp->rowptr           = C_ptr;

}

void ML_globalcsr2localcsr(ML_Operator *imatrix, int max_per_proc)
{
  int    lower, upper, col, i, j, k, Nexternal;
   int    *bindx, *externals;
   double *val;
   struct ML_CSR_MSRdata *temp;
   int    allocated, row_length;
   ML_Comm *comm;

   comm  = imatrix->comm;
   lower = max_per_proc*comm->ML_mypid;
   upper = lower + max_per_proc;

   allocated = 100;
   bindx = (int    *)  ML_allocate( allocated*sizeof(int   ));
   val   = (double *)  ML_allocate( allocated*sizeof(double));
   if (val == NULL) 
     pr_error("ML_back_to_csrlocal: out of space\n");

   Nexternal = 0;
   for (i = 0 ; i < imatrix->getrow->Nrows; i++) {
      ML_get_matrix_row(imatrix, 1, &i, &allocated, &bindx, &val,
                        &row_length, Nexternal);
      k = 0;
      for (j = 0; j < row_length; j++) {
         if (  (bindx[Nexternal+j] < lower) || (bindx[Nexternal+j] >= upper)){
            bindx[Nexternal+k++] = bindx[Nexternal+j];
         }
      }
      Nexternal += k;
   }
   ML_az_sort(bindx, Nexternal, NULL, NULL);
   ML_rm_duplicates(bindx, &Nexternal);
   externals = (int *) ML_allocate( (Nexternal+1)*sizeof(int));
   for (i = 0; i < Nexternal; i++) externals[i] = bindx[i];
   ML_free(bindx); bindx = NULL;
   ML_free(val); val = NULL;

   temp = (struct ML_CSR_MSRdata *) imatrix->data;

   for (i = 0 ; i < temp->rowptr[imatrix->getrow->Nrows]; i++) {
     col   = temp->columns[i];
     if ( (col >= lower) && (col < upper) ) col -= lower;
     else {
       j = ML_find_index(col,externals,Nexternal);
       if (j == -1) {
	 printf("Column not found: %d\n",col);
	 exit(1);
       }
       else col = j + imatrix->invec_leng;
     }
     temp->columns[i] = col;
   }

   ML_set_message_info(Nexternal, externals, max_per_proc,imatrix);

   ML_free(externals);

}


/******************************************************************************/
/* Create a map between local variables on this processor and a unique
 * global number where local variables on different processors which
 * correspond to the same global variable have the same unique global number.
 *
 * Parameters
 * ==========
 *   N_local       On input, number of local variables assigned to this node.
 *
 *   map           On output, map[k] is the unique global id of the kth local
 *                 variable. Note: if the kth local variable on processor P0
 *                 corresponds to the jth local variable on processor P1, then
 *                 map[k] on P0 is equal to map[j] on P1.
 *
 *   comm_info     On input, communcation information (see ml_rap.h) which
 *                 indicates which local variables are sent to other processors
 *                 and where received information is stored locally.
 ******************************************************************************/

void ML_back_to_csrlocal(ML_Operator *imatrix, ML_Operator *omatrix,
                         int max_per_proc)
{
   int    lower, upper, next_nz, col, i, j, k, Nexternal, ii;
   int    *bindx, *externals, *rowptr;
   double *val, dtemp;
   struct ML_CSR_MSRdata *temp;
   int    allocated, row_length;
   ML_Comm *comm;

   comm  = imatrix->comm;
   lower = max_per_proc*comm->ML_mypid;
   upper = lower + max_per_proc;

   allocated = imatrix->N_nonzeros+2;
   rowptr = (int   *)  ML_allocate( (1+imatrix->getrow->Nrows)*sizeof(int));
   bindx = (int    *)  ML_allocate( allocated*sizeof(int   ));
   val   = (double *)  ML_allocate( allocated*sizeof(double));
   if (val == NULL) 
     pr_error("ML_back_to_csrlocal: out of space\n");

   Nexternal = 0;
   for (i = 0 ; i < imatrix->getrow->Nrows; i++) {
      ML_get_matrix_row(imatrix, 1, &i, &allocated, &bindx, &val,
                        &row_length, Nexternal);
      k = 0;
      for (j = 0; j < row_length; j++) {
         if (  (bindx[Nexternal+j] < lower) || (bindx[Nexternal+j] >= upper)){
            bindx[Nexternal+k++] = bindx[Nexternal+j];
         }
      }
      Nexternal += k;
   }
   ML_az_sort(bindx, Nexternal, NULL, NULL);
   ML_rm_duplicates(bindx, &Nexternal);
   externals = (int *) ML_allocate( (Nexternal+1)*sizeof(int));
   for (i = 0; i < Nexternal; i++) externals[i] = bindx[i];

   rowptr[0] = 0;
   next_nz   = 0;
   for (i = 0 ; i < imatrix->getrow->Nrows; i++) {
      ML_get_matrix_row(imatrix, 1, &i, &allocated, &bindx, &val,
                        &row_length, next_nz);
      ii = next_nz;

      for (k = 0; k < row_length; k++) {
         dtemp = val[ii];
         col   = bindx[ii++];
         if ( (col >= lower) && (col < upper) ) col -= lower;
         else {
            j = ML_find_index(col,externals,Nexternal);
            if (j == -1) {
               printf("Column not found: %d\n",col);
               exit(1);
            }
            /* else col = j + imatrix->getrow->Nrows; */
            else col = j + imatrix->invec_leng;
         }
         bindx[next_nz] = col; val[next_nz++] = dtemp;
      }
      rowptr[i+1] = next_nz;
   }

   temp = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata) );
   temp->columns       = bindx;
   temp->values        = val;
   temp->rowptr        = rowptr;

   /* Note: outvec_leng is garbage because I don't update it in RAP */
   /*       Use imatrix->getrow->Nrows                              */

   omatrix->data_destroy = ML_CSR_MSRdata_Destroy;
   ML_Operator_Set_1Levels(omatrix, imatrix->from, imatrix->to);
   ML_Operator_Set_ApplyFuncData(omatrix, imatrix->invec_leng,
                             imatrix->getrow->Nrows, ML_EMPTY, (void*)temp,
                             imatrix->getrow->Nrows, NULL, 0);
   ML_Operator_Set_Getrow(omatrix, ML_EXTERNAL, imatrix->getrow->Nrows,
                          CSR_getrows);
   omatrix->max_nz_per_row = imatrix->max_nz_per_row;
   omatrix->N_nonzeros     = next_nz;
   ML_Operator_Set_ApplyFunc (omatrix, ML_INTERNAL, CSR_matvec);

   ML_set_message_info(Nexternal, externals, max_per_proc,omatrix);

   ML_free(externals);

}

/******************************************************************************/
/* Create a map between local variables on this processor and a unique
 * global number where local variables on different processors which
 * correspond to the same global variable have the same unique global number.
 *
 * Parameters
 * ==========
 *   N_local       On input, number of local variables assigned to this node.
 * 
 *   map           On output, map[k] is the unique global id of the kth local
 *                 variable. Note: if the kth local variable on processor P0
 *                 corresponds to the jth local variable on processor P1, then
 *                 map[k] on P0 is equal to map[j] on P1.
 *   
 *   comm_info     On input, communcation information (see ml_rap.h) which
 *                 indicates which local variables are sent to other processors
 *                 and where received information is stored locally.
 ******************************************************************************/

void ML_back_to_local(ML_Operator *imatrix, ML_Operator *omatrix,
                      int max_per_proc)
{
   int    lower, upper, next_nz, col, i, j, k, Nexternal, ii;
   int    *bindx, *externals;
   double *val, dtemp;
   struct ML_CSR_MSRdata *temp;
   int    allocated, row_length;
   int    N_nonzeros = 0, max_per_row = 0, count, lowest, num_PDEs;
   ML_Comm *comm;

   comm  = imatrix->comm;
   num_PDEs = imatrix->num_rigid;
   omatrix->num_PDEs  = num_PDEs;
   omatrix->num_rigid = imatrix->num_rigid;
   lower = max_per_proc*comm->ML_mypid;
   upper = lower + max_per_proc;

   allocated = imatrix->N_nonzeros+2;
   bindx = (int    *)  ML_allocate( allocated*sizeof(int   ));
   val   = (double *)  ML_allocate( allocated*sizeof(double));

   Nexternal = 0;
   for (i = 0 ; i < imatrix->getrow->Nrows; i++) {
      ML_get_matrix_row(imatrix, 1, &i, &allocated, &bindx, &val, 
                        &row_length, Nexternal);
      k = 0;
      for (j = 0; j < row_length; j++) {
         if (  (bindx[Nexternal+j] < lower) || (bindx[Nexternal+j] >= upper)){
            bindx[Nexternal+k++] = bindx[Nexternal+j];
         }
      }
      if (row_length > max_per_row) max_per_row = row_length;
      Nexternal += k;
      N_nonzeros += row_length;
   }
   ML_az_sort(bindx, Nexternal, NULL, NULL);
   ML_rm_duplicates(bindx, &Nexternal);

   /* Make sure that for each block (for block matrices), we either have */
   /* no point entries or all the point entries in the external list.    */
   count = 0; i = 0;
   while (i < Nexternal) {
      lowest = (int) floor( .000001 + ((double)bindx[i])/ ((double) num_PDEs) );
      lowest *= num_PDEs;
      for (j = 0; j < num_PDEs; j++) {
         count++;      
         if (i < Nexternal) {
            if (bindx[i] == lowest+j) i++;
         }
         else i++;
      }
   }
   externals = (int *) ML_allocate( (count+1)*sizeof(int));
   count = 0; i = 0;
   while (i < Nexternal) {
      lowest = (int) floor( .000001 + ((double)bindx[i])/ ((double) num_PDEs) );
      lowest *= num_PDEs;
      for (j = 0; j < num_PDEs; j++) {
         externals[count++] = lowest+j;
         if (i < Nexternal) {
            if (bindx[i] == lowest+j) i++;
         }
         else i++;
      }
   }
   Nexternal = count;

/*
   for (i = 0; i < Nexternal; i++) externals[i] = bindx[i];
*/
   ML_free(val);
   ML_free(bindx);
   bindx = (int    *)  ML_allocate( (N_nonzeros+5)*sizeof(int   ));
   val   = (double *)  ML_allocate( (N_nonzeros+5)*sizeof(double));

      /* need extra room (2) for diagonal guy and wasted space */
  
   bindx[0] = imatrix->getrow->Nrows+1;
   next_nz = bindx[0];
   for (i = 0 ; i < imatrix->getrow->Nrows; i++) {
      ML_get_matrix_row(imatrix, 1, &i, &allocated, &bindx, &val, 
                        &row_length, next_nz);
      ii = next_nz;

      val[i] = 0.;
      for (k = 0; k < row_length; k++) {
         dtemp = val[ii];
         col   = bindx[ii++];
         if ( (col >= lower) && (col < upper) ) col -= lower;
         else {
            j = ML_find_index(col,externals,Nexternal);
            if (j == -1) {
               printf("Column not found: %d\n",col);
               exit(1);
            }
            else col = j + imatrix->getrow->Nrows;
         }
         if (col == i) { val[i] = dtemp; }
         else 
            if (dtemp != 0.) { bindx[next_nz] = col; val[next_nz++] = dtemp; }
      }
      bindx[i+1] = next_nz;
   }

   temp = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata) );
   temp->columns       = bindx;
   temp->values        = val;
   temp->rowptr        = NULL;

   /* Note: outvec_leng is garbage because I don't update it in RAP */
   /*       Use imatrix->getrow->Nrows                              */

   omatrix->data_destroy = ML_CSR_MSRdata_Destroy;
   ML_Operator_Set_1Levels(omatrix, imatrix->from, imatrix->to);
   ML_Operator_Set_ApplyFuncData(omatrix, imatrix->invec_leng, 
                             imatrix->getrow->Nrows, ML_EMPTY, (void*)temp, 
                             imatrix->getrow->Nrows, NULL, 0);
   ML_Operator_Set_Getrow(omatrix, ML_EXTERNAL, imatrix->getrow->Nrows, 
		          MSR_getrows);
   omatrix->max_nz_per_row = max_per_row;
   omatrix->N_nonzeros     = N_nonzeros;
   ML_Operator_Set_ApplyFunc (omatrix, ML_INTERNAL, MSR_matvec);
   ML_Operator_Set_Diag (omatrix, imatrix->getrow->Nrows, temp->values);

   ML_set_message_info(Nexternal, externals, max_per_proc,omatrix);

   ML_free(externals);

}

/*******************************************************************************
  Perform initializations so that local communication can occur to update the
  external elements. This initialization includes:

   1) determine the number of neighbors to which we send information as well as
      the number of neighbors from which we receive information.  These two
      should be the same. If they are not, we will set up things so that we 
      send 0 length messages to processors from which we receive (but have no
      information to send them).
   2) determine the total number of unknowns that we must send out.  Using this
      information we allocate space.
   3) Initialize the communication operator for 'matrix' so that it contains 
      the number of messages to be sent/received, the node numbers of the 
      processors to which we must send and receive, the length of the messages 
      that we must send and that we expect to receive from each of the neighbors, 
      and finally, a list of the indices of the elements that will be send to 
      other processors (in the order that they will be sent).

      NOTE: Implicitly the neighbors are numbered using the ordering of the
      external elements. In particular, the external elements are ordered such
      that those elements updated by the same processor are contiguous. In this
      way, the ordering of the external elements defines an ordering for the
      neighbors.


  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  N_external:      Number of external elements on this processor.

  N_rows:          Number of elements updated on this processor.

  external:        List (global indices) of external elements on this node.


*******************************************************************************/

void ML_set_message_info(int N_external, int external[], int max_per_proc, 
                         ML_Operator *matrix)

{
  /* local variables */

  int   num_recv_neighbors, num_send_neighbors, total_to_be_sent;
  int   i, j, ii, type, start, partner, proc, nprocs, neigh_proc;
  int  *neighbors, *tempneigh;
  USR_REQ     *request;
  int length;
  int *rcv_list, total;
  int *send_lengths, *rcv_lengths, *send2_ptr;
  int count, count2;
  ML_Comm *comm;

  char *yo = "ML_set_message_info: ";

  comm   = matrix->comm;
  proc   = comm->ML_mypid;
  nprocs = comm->ML_nprocs;

  neighbors = (int *) ML_allocate(nprocs*sizeof(int));
  tempneigh = (int *) ML_allocate(nprocs*sizeof(int));

  for (i = 0 ; i < nprocs ; i++ ) neighbors[i] = 0;

 /*
   * Count the number of neighbors from which we receive information to update
   * our external elements. Additionally, fill the array neighbors in the
   * following way:
   *      neighbors[i] = 0   ==>  No external elements are updated by
   *                              processor i.
   *      neighbors[i] = x   ==>  (x-1)/nprocs elements are updated from
   *                              processor i.
   */

  num_recv_neighbors = 0;

  for (i = 0; i < N_external; i++) {
    neigh_proc  = external[i]/max_per_proc;
    if ( neigh_proc >= nprocs ) {
       printf("%d : ML_set_message_info warning : index out of bound \
             = %d(%d) %d %d(%d)\n", proc, neigh_proc, nprocs, max_per_proc, 
             external[i], i);
       exit(1);
    }
    if (neighbors[neigh_proc] == 0) {
      num_recv_neighbors++;
      neighbors[neigh_proc] = 1;
    }
    neighbors[neigh_proc] += nprocs;
  }

  /* sum over all processors all the neighbors arrays */

  ML_gsum_vec_int(neighbors, tempneigh, comm->ML_nprocs, comm);

  /* decode the combined 'neighbors' array from all the processors */

  num_send_neighbors = neighbors[proc] % nprocs;

  /* decode 'neighbors[proc] to deduce total number of elements we must send */

  total_to_be_sent = (neighbors[proc] - num_send_neighbors) / nprocs;

  ML_free(neighbors);
  ML_free(tempneigh);


 /*
   * Make a list of the neighbors that will send information to update our
   * external elements (in the order that we will receive this information).
   */

  neighbors = (int *) ML_allocate((num_recv_neighbors+num_send_neighbors+1)*
                                  sizeof(int) );
  if (neighbors == NULL) {
    (void) fprintf(stderr, "%sERROR: Not enough dynamic space.\n", yo);
    exit(-1);
  }

  if (N_external > 0) {
     j = 0;
     neigh_proc = external[0]/max_per_proc;
     neighbors[j++] = neigh_proc;
     for (i = 1; i < N_external; i++) {
        ii = external[i]/max_per_proc;
        if ( neigh_proc != ii) neighbors[j++] = ii;
        neigh_proc = ii;
     }
  }

  type = 1992;

  request = (USR_REQ    *) ML_allocate((num_recv_neighbors+num_send_neighbors+1)
                                           *sizeof(USR_REQ    ));

  /* first post receives */

  length = 0;
  for (i = 0; i < num_send_neighbors; i++) {
    partner = -1;
    comm->USR_irecvbytes((void *) external, (unsigned int) length, &partner, 
                         &type, comm->USR_comm, request+i);
  }

  /* send messages */

  for (i = 0; i < num_recv_neighbors; i++) {
    comm->USR_sendbytes((void *) external, 0, neighbors[i],type,comm->USR_comm);
  }

    /*
   * Receive message from each send neighbor to construct send list.
   */

  length = 0;
  for (i = 0; i < num_send_neighbors; i++) {
    partner = -1;
    comm->USR_waitbytes((void *) external, length, &partner, &type, 
                        comm->USR_comm, request+i);
    if (partner != -1) neighbors[num_recv_neighbors++] = partner;
  }

  /* make one neighbor list based on recv_neighbors and send_neighbors */

  ML_az_sort(neighbors, num_recv_neighbors, NULL, NULL);
  ML_rm_duplicates(neighbors, &num_recv_neighbors);

  num_send_neighbors = num_recv_neighbors;

  send_lengths =(int *) ML_allocate( (num_recv_neighbors+1)*sizeof(int));
  rcv_lengths  =(int *) ML_allocate( (num_recv_neighbors+1)*sizeof(int));
  send2_ptr    =(int *) ML_allocate( (total_to_be_sent+1)*sizeof(int));
  if ( (rcv_lengths == NULL) || (send2_ptr   == NULL) ) {
    (void) fprintf(stderr, "%sERROR: Not enough dynamic space.\n", yo);
    exit(-1);
  }

  /*
  * Send each processor the global index list of the external elements in the
  * order that I will want to receive them when updating my external elements
  */

  type = 1993;

  for (i = 0; i < num_recv_neighbors; i++) {
    length     = sizeof(int);
    partner    = neighbors[i];
    comm->USR_irecvbytes((void *) &(send_lengths[i]), 
		   (unsigned int) length, 
                   &partner, &type, comm->USR_comm, request+i);
  }

  j = 0;
  for (i = 0; i < num_recv_neighbors; i++) {
    length = 0;

    /* go through list of external elements until updating processor changes */

    while ((j < N_external) && ( external[j]/max_per_proc == neighbors[i])) {
      length++; j++;
    }

    rcv_lengths[i] = length;

    comm->USR_sendbytes((void *) &length, sizeof(int), neighbors[i], type, 
                        comm->USR_comm);
  }

  /* receive from each neighbor the send length */

  for (i = 0; i < num_recv_neighbors; i++) {
    length     = sizeof(int);
    partner    = neighbors[i];
    comm->USR_waitbytes((void *) &(send_lengths[i]), length, 
                 &partner, &type, comm->USR_comm, request+i);
  }


  type = 1994;

  start = 0;
  for (i = 0; i < num_recv_neighbors; i++) {
    length     = send_lengths[i] * sizeof(int);
    partner    = neighbors[i];
    comm->USR_irecvbytes((void *) &(send2_ptr[start]),(unsigned int) length,
                                  &partner,&type, comm->USR_comm, request+i);
    start     += send_lengths[i];
  }

  start = 0;
  for (i = 0; i < num_recv_neighbors; i++) {
    length     = rcv_lengths[i] * sizeof(int);
    partner    = neighbors[i];
    comm->USR_sendbytes((void *) &(external[start]), (unsigned int) length, 
                        partner, type, comm->USR_comm);
    start     += rcv_lengths[i];
  }

  /* receive from each neighbor the global index list of external ele */

  start = 0;
  for (i = 0; i < num_recv_neighbors; i++) {
    length     = send_lengths[i] * sizeof(int);
    partner    = neighbors[i];
    comm->USR_waitbytes((void *) &(send2_ptr[start]), length, &partner, &type, 
                        comm->USR_comm, request+i);
    start += send_lengths[i];
  }

   ii = max_per_proc*proc;
   for (i = 0; i < total_to_be_sent; i++) send2_ptr[i] -= ii;

   ML_free(request);
   total = 0;
   for (i = 0; i < num_recv_neighbors; i++)
      total += rcv_lengths[i];
   rcv_list = (int *) ML_allocate( (total+1)*sizeof(int));
   for (i = 0; i < total; i++) rcv_list[i] = i + matrix->invec_leng;

   ML_CommInfoOP_Set_neighbors(&(matrix->getrow->pre_comm), 
                               num_recv_neighbors, neighbors, 
                               ML_OVERWRITE, NULL, 0);
   count  = 0;
   count2 = 0;
   for (i = 0; i < num_recv_neighbors; i++) {
      ML_CommInfoOP_Set_exch_info(matrix->getrow->pre_comm,
                                  neighbors[i], rcv_lengths[i], 
				  &(rcv_list[count2]),
                                  send_lengths[i], &(send2_ptr[count]));
      count  += send_lengths[i];
      count2 += rcv_lengths[i];
   }
   ML_free(neighbors);
   ML_free(send_lengths);
   ML_free(rcv_lengths);
   ML_free(send2_ptr);

   ML_free(rcv_list);

}
/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

void ML_convert_data_org(ML_Operator *matrix, int data_org[],
        int rcv_list[], int remap[], int leng, int add_or_not)
{
   int i, count, count2;



    ML_CommInfoOP_Set_neighbors( &(matrix->getrow->pre_comm), 
			   data_org[ML_N_neigh],&(data_org[ML_neighbors]),
                           add_or_not, remap, leng);

    count = ML_send_list;
    count2 = 0;

    if (rcv_list == NULL) {
       for (i = 0; i < data_org[ML_N_neigh]; i++) {
          ML_CommInfoOP_Set_exch_info(matrix->getrow->pre_comm, 
		    data_org[ML_neighbors+i], data_org[ML_rec_length+i], NULL, 
		    data_org[ML_send_length+i], &(data_org[count]));
          count += data_org[ML_send_length+i];
       }
    }
    else {
       for (i = 0; i < data_org[ML_N_neigh]; i++) {
          ML_CommInfoOP_Set_exch_info(matrix->getrow->pre_comm, 
		    data_org[ML_neighbors+i], data_org[ML_rec_length+i], 
		    &(rcv_list[count2]), data_org[ML_send_length+i],
		    &(data_org[count]));
          count2 += data_org[ML_rec_length+i];
          count += data_org[ML_send_length+i];
       }
    }
}

