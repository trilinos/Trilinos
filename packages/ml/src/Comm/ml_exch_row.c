/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "ml_struct.h"
#include "ml_memory.h"
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

#define ML_FUNCTION_NAME ML_exchange_rows
void ML_exchange_rows(ML_Operator *Pmatrix, ML_Operator **Pappended, 
                      ML_CommInfoOP *comm_info)
{

  /* local variables */

  int         *actual_send_length, *actual_recv_length, *start_send_proc;
  int         Nneighbors, *neighbor, *remap;
  int         Nrows_new, *cols_new, *rowptr_new, Nrows_send, Nrows, max_per_row;
  int         *ibuff, total_num_recv, total_send, total_recv, *iptr;
  double      *vals_new, *dtemp, *dbuff, *dummy1, *dptr;
  int         i, j, k, ii, jj, *newmap, *orig_map, nonNULL_rcv_list, *dummy2;
  static int  type = ML_MPI_MSG_NUM;
  struct      ML_CSR_MSRdata *temp;
  int         allocated_space, row_length;
  ML_Comm     *comm;
  int Nghost;
  int rcv_list_exists = 0, count = 0;
  int mypid = Pmatrix->comm->ML_mypid;

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
  if (Nghost < Nrows_new) Nghost = Nrows_new;

  actual_recv_length = (int *) ML_allocate( (Nneighbors+1)*sizeof(int));
  actual_send_length = (int *) ML_allocate( (Nneighbors+1)*sizeof(int));
  start_send_proc    = (int *) ML_allocate( (Nneighbors+1)*sizeof(int));
  rowptr_new         = (int *) ML_allocate( (Nrows_new+1)*sizeof(int));
  if (rowptr_new == NULL) {
     pr_error("(%d) %s, line %d: Out of space in %s\n   tried to allocate %d ints\n", mypid, __FILE__,__LINE__, ML_FUNCTION_NAME, Nrows_new+1);
  }

  /* compute the row lengths of the external rows    */
  /* and the total number of nonzeros to be received */
  /* from each processor.                            */

  allocated_space = Pmatrix->max_nz_per_row+2;
  dummy1 = (double *) ML_allocate(allocated_space*sizeof(double));
  dummy2 = (int    *) ML_allocate(allocated_space*sizeof(   int));
#if defined(HXVE_ML_PARMETIS_2x) || defined(HXVE_ML_PARMETIS_3x)
  dtemp  = (double *) ML_allocate(2*(Nrows+Nghost + 1)*sizeof(double));
#else
  dtemp  = (double *) ML_allocate((Nrows+Nghost + 1)*sizeof(double));
#endif
  
  if (dtemp == NULL) 
  {
     pr_error("(%d) %s, line %d: Out of space in %s\n   tried to allocate %d + %d doubles\n",mypid,__FILE__,__LINE__,ML_FUNCTION_NAME,Nrows,Nghost);
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
  ML_free(dummy2);
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
  if (allocated_space-2 > Pmatrix->max_nz_per_row) 
     Pmatrix->max_nz_per_row = allocated_space;

  max_per_row = Pmatrix->max_nz_per_row;

  for (i = 1 ; i < Nrows_new; i++) 
  {
     k = rowptr_new[i];
     rowptr_new[i] = rowptr_new[i-1] + j;
     j = k;
     if (j > max_per_row) max_per_row = j;
  }
  if (Nrows_new != 0) rowptr_new[Nrows_new] = rowptr_new[Nrows_new-1] + j;
  total_recv = rowptr_new[Nrows_new];

  /* compute number of nonzeros to receive from each neighbor */

  if (rcv_list_exists == 0) 
  {
     k = 0;
     for (j = 0; j < Nneighbors; j++) 
     {
        actual_recv_length[j] = rowptr_new[k + comm_info->neighbors[j].N_rcv] - 
                             rowptr_new[k];
        k += comm_info->neighbors[j].N_rcv;
     }
  }

  /* allocate integer space for the new matrix rows */
  /* allocate space to pack the messages to be sent */

  allocated_space = total_send+1;
  ibuff    = (int    *) ML_allocate( allocated_space*sizeof(int   ));
  if (ibuff == NULL) 
  {
     pr_error("(%d) %s, line %d: Out of space in %s\n   tried to allocate %d ints (Nneighbors = %d)\n", mypid, __FILE__,__LINE__, ML_FUNCTION_NAME, allocated_space,Nneighbors);
  }

  /* pack up the integer information to be sent and send it */

  i = 0; k = 0;
  for (ii = 0; ii < Nneighbors; ii++) 
  {
     start_send_proc[ii]  = i;
     for (jj = 0 ; jj < comm_info->neighbors[ii].N_send; jj++) 
     {
        j = comm_info->neighbors[ii].send_list[jj];
	iptr = &(ibuff[i]);
        ML_get_matrix_row(Pmatrix, 1, &j, &allocated_space, &iptr,
                          &dummy1, &row_length, 0);
        i += row_length;
     }
     actual_send_length[ii] = i - start_send_proc[ii];
  } 
  neighbor    = (int *) ML_allocate(Nneighbors*sizeof(int));
  if ((Nneighbors != 0) && (neighbor == NULL)) 
  {
     printf("Not enough space in ML_exchange_rows\n");
     printf("   tried to allocate %d ints for neighbor\n",Nneighbors);
     exit(1);
  }
  for (i = 0; i < Nneighbors; i++) 
     neighbor[i] = comm_info->neighbors[i].ML_id;

  if (rcv_list_exists) 
  {
     j = Nrows;
     i = 0; count = 0;
     for (ii = 0; ii < Nneighbors; ii++) 
     {
        actual_recv_length[ii] = 0;
        for (jj = 0 ; jj < comm_info->neighbors[ii].N_rcv; jj++) 
        {
           while (dtemp[j] == -1.) j++;
/*
           j = comm_info->neighbors[ii].rcv_list[jj];
*/
           actual_recv_length[ii] += (int) dtemp[j];
           rowptr_new[i++] = count;   count += (int) dtemp[j];
           j++;
        }
     }
     rowptr_new[i] = count;
  }
  ML_free(dtemp);
  cols_new = (int    *) ML_allocate( (total_recv+1)*sizeof(int));
  if (cols_new == NULL) {
     pr_error("(%d) %s, line %d: Out of space in %s\n   tried to allocate %d ints (Nneighbors = %d)\n", mypid, __FILE__,__LINE__, ML_FUNCTION_NAME, total_recv+1,Nneighbors);

  }

  type++;
  if (type > ML_MPI_MSG_NUM + 100) type = ML_MPI_MSG_NUM;
  ML_splitup_big_msg(Nneighbors,(char *) ibuff,(char *) cols_new, 
                     sizeof(int), start_send_proc, actual_send_length, 
                     actual_recv_length, neighbor, type, 
                     &total_num_recv, comm);
  ML_free(ibuff);
  dbuff    = (double *) ML_allocate( allocated_space*sizeof(double));
  if (dbuff == NULL) 
    {
     pr_error("(%d) %s, line %d: Out of space in %s\n   tried to allocate %d doubles (Nneighbors = %d)\n", mypid, __FILE__,__LINE__, ML_FUNCTION_NAME, allocated_space,Nneighbors);
    }

  /* pack up the float information to be sent and send it */

  i = 0; k = 0;
  for (ii = 0; ii < Nneighbors; ii++) 
  {
     for (jj = 0 ; jj < comm_info->neighbors[ii].N_send; jj++) 
     {
        j = comm_info->neighbors[ii].send_list[jj];
	dptr = &(dbuff[i]);
        ML_get_matrix_row(Pmatrix, 1, &j, &allocated_space, (int **) &dummy1,
                          &dptr, &row_length, 0);
        i += row_length;
     } 
  } 
  ML_free(dummy1);
  vals_new = (double *) ML_allocate( (total_recv+1)*sizeof(double));
  if (vals_new == NULL) {
    printf("out of space in ML_exchange_rows\n");
    printf("   tried to allocate %d doubles for vals_new\n",allocated_space);
    exit(1);
  }

  type++;
  if (type > ML_MPI_MSG_NUM + 100) type = ML_MPI_MSG_NUM;
  ML_splitup_big_msg(Nneighbors,(char *) dbuff,(char *) vals_new, 
                     sizeof(double), start_send_proc, actual_send_length, 
                     actual_recv_length, neighbor, type, 
                     &total_num_recv, comm);

  ML_free(neighbor);
  ML_free(dbuff);
  ML_free(start_send_proc);
  ML_free(actual_recv_length);
  ML_free(actual_send_length);

  temp = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata) );
  temp->columns       = cols_new;
  temp->values        = vals_new;
  temp->rowptr        = rowptr_new;

  *Pappended = ML_Operator_Create(comm);

  /* Estimate the total number of columns in Pappended  */
  /* This is done by taking the total number of columns */
  /* in Pmatrix and adding 1/4 the number of nonzeros   */
  /* received in the expanded matrix. This is a total   */
  /* guess.                                             */

  if (Pmatrix->N_total_cols_est != -1) 
    (*Pappended)->N_total_cols_est = Pmatrix->N_total_cols_est;
  else
    (*Pappended)->N_total_cols_est = Pmatrix->invec_leng;

  (*Pappended)->N_total_cols_est += count/4;

  ML_Operator_Set_1Levels(*Pappended, Pmatrix->from, Pmatrix->to);
  ML_Operator_Set_ApplyFuncData(*Pappended,Pmatrix->invec_leng, 
                             /* Nrows+Nrows_new, */ Nrows+Nrows_new,
                             (void*)temp,Nrows+Nrows_new,NULL,0);

  ML_Operator_Set_Getrow(*Pappended, Nrows + Nrows_new, CSR_getrow);
  (*Pappended)->max_nz_per_row = max_per_row;
  if (Pmatrix->N_nonzeros >= 0) 
     (*Pappended)->N_nonzeros     = Pmatrix->N_nonzeros + rowptr_new[Nrows_new];


  (*Pappended)->sub_matrix = Pmatrix; 

  if (Pmatrix->getrow->row_map == NULL)
     (*Pappended)->getrow->row_map = NULL;
  else 
  {
     (*Pappended)->getrow->row_map = 
           (int *) ML_allocate((*Pappended)->getrow->Nrows* sizeof(int));
     for (i =0; i < Nrows; i++) 
        (*Pappended)->getrow->row_map[i] = Pmatrix->getrow->row_map[i];
     for (i =Nrows; i < Nrows + Nrows_new; i++) 
        (*Pappended)->getrow->row_map[i] = i;
  }

  /* If the user wishes to remap the rows of the matrix corresponding to */
  /* Pmatrix, set or change the map to reflect this.                     */

  remap = comm_info->remap;
  if (remap != NULL ) 
  {
     newmap = (int *) ML_allocate(((*Pappended)->getrow->Nrows+1)*sizeof(int));
     if (newmap == NULL) 
     {
        printf("Not enough space in ML_exchange_rows for remap\n");
        printf("   tried to allocate %d ints for newmap\n",
               (*Pappended)->getrow->Nrows+1);
        exit(1);
     }
     for (i = 0; i < Nrows; i++) newmap[i] = -1;
     orig_map = (*Pappended)->getrow->row_map;
     if (orig_map == NULL) 
     {
        for (i = Nrows; i < (*Pappended)->getrow->Nrows; i++) newmap[i] = i;
        for (i = 0; i < Nrows; i++) {
           if (remap[i] != -1) newmap[remap[i]] = i;
        }
     }
     else 
     {
        for (i = Nrows; i < (*Pappended)->getrow->Nrows; i++) 
        newmap[i] = orig_map[i];
        for (i = 0; i < Nrows; i++)
           if (remap[i] != -1) newmap[remap[i]] = orig_map[i];
        if ( Pmatrix->getrow->row_map != orig_map) ML_free(orig_map);
     }
     (*Pappended)->getrow->row_map = newmap;
  }

  /* If there is a Rcv_list and the user has not set the add option, */
  /* then simply remap the received rows of the matrix to reflect    */
  /* the receive list.                                               */

  if ( (nonNULL_rcv_list == 1) && (comm_info->add_rcvd == 0)) 
  {
#if defined(HxVE_ML_PARMETIS_2x) || defined(HxVE_ML_PARMETIS_3x)
    j = (*Pappended)->getrow->Nrows;
     newmap = (int *) ML_allocate( j * sizeof(int));
     for (i = Nrows; i < j; i++) newmap[i] = -1;
#else
     newmap = (int *) ML_allocate( (Nrows + Nghost) * sizeof(int));
     for (i = Nrows; i < Nrows + Nghost; i++) newmap[i] = -1;
#endif
     
     orig_map = (*Pappended)->getrow->row_map;
     if (orig_map == NULL) 
     {
        for (i = 0; i < Nrows; i++) newmap[i] = i;
     }
     else 
     {
        for (i = 0; i < (*Pappended)->getrow->Nrows; i++) 
        {
           if (orig_map[i] < Nrows)  newmap[i] = orig_map[i];
        }
        if (Pmatrix->getrow->row_map != orig_map) ML_free(orig_map);
     }
     ii = Nrows;
     for (i = 0; i < comm_info->N_neighbors; i++) 
     {
        for (j = 0; j < comm_info->neighbors[i].N_rcv; j++) 
        {
           newmap[comm_info->neighbors[i].rcv_list[j]] = ii++;
        }
     }
     (*Pappended)->getrow->row_map = newmap;
  }

  if ( (nonNULL_rcv_list == 1) && (comm_info->add_rcvd == 1)) 
  {
     ML_add_appended_rows(comm_info,*Pappended, Nrows, Nrows_new,
                          rowptr_new[Nrows_new]);
  }

  /* Count the number of rows */
  if ((*Pappended)->getrow->row_map != NULL) 
  {
     for (i = (*Pappended)->getrow->Nrows - 1; i >= 0; i--) 
     {
        if ((*Pappended)->getrow->row_map[i] != -1) break; 
     }
     /* Should probably try and work in an ML call here */ 
     (*Pappended)->getrow->Nrows = i+1;
     (*Pappended)->outvec_leng   = i+1;
  }
} /* ML_exchange_rows */
#ifdef ML_FUNCTION_NAME
#undef ML_FUNCTION_NAME
#endif

/******************************************************************************/
/******************************************************************************/
/*
 * On input, 'matrix' has been appended with rows that have been received
 * from other processors. This routine will modify 'matrix' so that those
 * new rows are placed (actually added) into the matrix according to the 
 * row numbers given in the receive list. If a row already exists, then
 * the received row will be added to the already existing row.
 *
 * Parameters
 * ==========
 *
 *    comm_info       On input, structure indicating what information has
 *                    been received from what processor and where the appeneded
 *                    rows need to be stored. See ml_rap.h for more info.
 *
 *    matrix          On input, a matrix where rows received from other
 *                    processors are appended to the matrix. On output,
 *                    the appended rows are added into the matrix. That is,
 *                    each appended row is removed from the end of the matrix
 *                    and instead added to the row that corresponds to the
 *                    receive list.
 *
 *    orig_rows       On input, number of rows in the non-appended portion of 
 *                    the matrix.
 *
 *    total_rcvd      On input, the number of rows appended to the matirx.
 *
 *    appended_nzs    On input, the number of nonzeros in the appended portion
 *                    of the matrix.
 ******************************************************************************/

#define ML_FUNCTION_NAME ML_add_appended_rows
void ML_add_appended_rows(ML_CommInfoOP *comm_info, ML_Operator *matrix, 
                          int orig_rows, int total_rcvd, int appended_nzs)
{

   int i, ii, jj, iii, jjj, i2, start, k;
   int Ncols, row_location, new_row, row_count;
   int *accum_col, accum_size, max_nz_per_row, next_nz, total_nz = 0;
   int max_nz_row_new, N_changed = 0, t_changed = 0, sub_i = 0, total;
   int *Ccol, *C_ptr, *row_ids, N_append_rows, *row_map, *itemp;
   double *accum_val, *Cval, dtemp;
   ML_Operator *current, *previous_matrix, *parent;
   struct ML_CSR_MSRdata *temp = NULL;
   int row_request, row_length;

   parent = matrix;

   /* determine number of new rows (N_append_rows) in the resulting matrix */
 
   N_append_rows = total_rcvd;
   row_ids = (int *) ML_allocate(N_append_rows*sizeof(int));
   i = 0;
   for (ii = 0; ii < comm_info->N_neighbors; ii++)
      for (jj = 0; jj < comm_info->neighbors[ii].N_rcv; jj++)
         row_ids[i++] = comm_info->neighbors[ii].rcv_list[jj];
   ML_az_sort(row_ids, N_append_rows, NULL, NULL);
   ML_rm_duplicates(row_ids, &N_append_rows);
   ML_free(row_ids);
 
   /* allocate space */

   max_nz_per_row  = matrix->max_nz_per_row;
   max_nz_row_new  = max_nz_per_row;
   accum_size      = 5*matrix->max_nz_per_row+ 2;
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
   if (total <= matrix->max_nz_per_row) total = matrix->max_nz_per_row + 1;
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
                                   (void*)temp,row_count,NULL,0);
               ML_Operator_Set_Getrow(current, row_count, CSR_getrow);
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
   matrix->getrow->func_ptr = CSR_getrow;
   matrix->getrow->ML_id = ML_NONEMPTY;
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
#ifdef ML_FUNCTION_NAME
#undef ML_FUNCTION_NAME
#endif

#define ML_FUNCTION_NAME ML_globalcsr2localcsr
void ML_globalcsr2localcsr(ML_Operator *imatrix, int max_per_proc)
{
  int    lower, upper, col, i, j, k, Nexternal;
   int    *bindx, *externals;
   double *val;
   struct ML_CSR_MSRdata *temp;
   int    allocated, row_length;
   ML_Comm *comm;
   int mypid  = imatrix->comm->ML_mypid;

   comm  = imatrix->comm;
   lower = max_per_proc*comm->ML_mypid;
   upper = lower + max_per_proc;

   allocated = 100;
   bindx = (int    *)  ML_allocate( allocated*sizeof(int   ));
   val   = (double *)  ML_allocate( allocated*sizeof(double));
   if (val == NULL)
     pr_error("(%d) %s, line %d: Out of space in %s\n   tried to allocate %d doubles\n",mypid,__FILE__,__LINE__,ML_FUNCTION_NAME,allocated);

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

   ML_CommInfoOP_GenUsingGIDExternals(Nexternal, externals, max_per_proc,imatrix);

   ML_free(externals);

}
#ifdef ML_FUNCTION_NAME
#undef ML_FUNCTION_NAME
#endif


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

#define ML_FUNCTION_NAME ML_back_to_csrlocal
void ML_back_to_csrlocal(ML_Operator *imatrix, ML_Operator *omatrix,
                         int max_per_proc)
{
   int    lower, upper, next_nz, col, i, j, k, Nexternal, ii;
   int    *bindx, *externals, *rowptr;
   double *val, dtemp;
   struct ML_CSR_MSRdata *temp;
   int    allocated, row_length;
   ML_Comm *comm;
   int mypid = imatrix->comm->ML_mypid;

   comm  = imatrix->comm;
   lower = max_per_proc*comm->ML_mypid;
   upper = lower + max_per_proc;

   allocated = imatrix->N_nonzeros+2;
   rowptr = (int   *)  ML_allocate( (1+imatrix->getrow->Nrows)*sizeof(int));
   bindx = (int    *)  ML_allocate( allocated*sizeof(int   ));
   val   = (double *)  ML_allocate( allocated*sizeof(double));
   if (val == NULL)  {
     pr_error("(%d) %s, line %d: Out of space in %s\n   tried to allocate %d doubles\n",mypid,__FILE__,__LINE__,ML_FUNCTION_NAME,allocated);

   }

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
	 /* rst: changing this to weed out any zero values */
	 /* in the final matrix that is created. A long    */
	 /* time ago ... we used to do this. I'm not sure  */
	 /* if there is a reason that we stopped weeding   */
	 /* out zeros.                                     */
	 if (dtemp != 0.) {
	   bindx[next_nz] = col; val[next_nz++] = dtemp;
	 }
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
                             imatrix->getrow->Nrows, (void*)temp,
                             imatrix->getrow->Nrows, NULL, 0);
   ML_Operator_Set_Getrow(omatrix, imatrix->getrow->Nrows,
                          CSR_getrow);
   omatrix->max_nz_per_row = imatrix->max_nz_per_row;
   omatrix->N_nonzeros     = next_nz;
   ML_Operator_Set_ApplyFunc (omatrix, CSR_matvec);

   ML_CommInfoOP_GenUsingGIDExternals(Nexternal, externals, max_per_proc,omatrix);

   ML_free(externals);

}
#ifdef ML_FUNCTION_NAME
#undef ML_FUNCTION_NAME
#endif

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

#define ML_FUNCTION_NAME ML_back_to_local
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
                             imatrix->getrow->Nrows, (void*)temp, 
                             imatrix->getrow->Nrows, NULL, 0);
   ML_Operator_Set_Getrow(omatrix, imatrix->getrow->Nrows, 
		          MSR_getrows);
   omatrix->max_nz_per_row = max_per_row;
   omatrix->N_nonzeros     = N_nonzeros;
   ML_Operator_Set_ApplyFunc (omatrix, MSR_matvec);
   ML_Operator_Set_Diag (omatrix, imatrix->getrow->Nrows, temp->values);

   ML_CommInfoOP_GenUsingGIDExternals(Nexternal, externals, max_per_proc,omatrix);

   ML_free(externals);

}
#ifdef ML_FUNCTION_NAME
#undef ML_FUNCTION_NAME
#endif

/******************************************************************************
  Perform initializations so local communication can occur to update the
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
      that we must send and that we expect to receive from each neighbor, 
      and finally, a list of the indices of the elements that will be send to 
      other processors (in the order that they will be sent).


  Author:          Ray S. Tuminaro, SNL, 9214
  =======
  Return code:     void
  ============
  Parameter list:
  ===============

  N_external:      Number of external elements on this processor. 

  external:        List (global indices) of external elements on this node.
                   This list should be ordered so that all externals from
                   the same processor appear consecutively. It is assumed 
                   that the global elements are defined as 
                   gid = lid + max_per_proc*proc_id 
                   (normally computed by ML_create_unique_col_id).

  max_per_proc     Maximum number of local elements per processor
                   (normally computed by ML_create_unique_col_id).

  matrix           On input, matrix for which we wish to compute communication
                   information. On output, matrix->getrow->pre_comm is
                   set so that local communications can be performed on 
                   the matrix.
                  
*****************************************************************************/

#define ML_FUNCTION_NAME ML_CommInfoOP_GenUsingGIDExternals
void ML_CommInfoOP_GenUsingGIDExternals(int N_external, int external[], 
					int max_per_proc, ML_Operator *matrix)

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

  char *yo = "ML_CommInfoOP_GenUsingGIDExternals: ";

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
       printf("%d : ML_CommInfoOP_GenUsingGIDExternals warning : index out of bound \
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

  ML_gsum_vec_int(&neighbors, &tempneigh, comm->ML_nprocs, comm);

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
    comm->USR_cheapwaitbytes((void *) external, length, &partner, &type, 
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
    comm->USR_cheapwaitbytes((void *) &(send_lengths[i]), length, 
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
    comm->USR_cheapwaitbytes((void *) &(send2_ptr[start]), length, &partner, &type, 
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
#ifdef ML_FUNCTION_NAME
#undef ML_FUNCTION_NAME
#endif
/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

#define ML_FUNCTION_NAME ML_convert_data_org
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
#ifdef ML_FUNCTION_NAME
#undef ML_FUNCTION_NAME
#endif

