#include "ml_operator.h"
#include "ml_xyt.h"
#include "ml_utils.h"
#include "ml_memory.h"

void setup_henry(ML *my_ml, int grid0, int **imapper, int **separator,
        int **sep_size, int *Nseparators, int *Nlocal, int *Nghost,
        ML_Operator **matvec_data) {

   int *sep, *s_sizes, *mapper;
   int i, j, jj, kk, Nsend, Nrecv, N_procs, N_bits;
   int start, mask2;
   int Nrows;
   ML_Operator *Amat, *omatrix;
   ML_CommInfoOP *getrow_comm;
   int total_rows, sep_space, sub_mask, N_count, old_Ncount, mask, proc;
   double *dmapper;
   int N_nz, allocated, *bindx, *row_ptr, row_length, max_nz_per_row;
   int *neighbors;
   double *val;
   struct ML_CSR_MSRdata *temp;

   Amat = &(my_ml->Amat[grid0]);
   if (Amat->matvec->internal == NULL)
      perror("Can't get number of rows because matvec not set\n");

   if (Amat->getrow->internal == NULL)
      perror("Get row not set!!! Can't setup henry\n");


   Nrows = Amat->getrow->Nrows;
   getrow_comm = Amat->getrow->pre_comm;

   if (getrow_comm == NULL)
      perror("No communication information for getrow! Can't setup henry\n");

   N_procs = my_ml->comm->ML_nprocs;
   proc    = my_ml->comm->ML_mypid;
   *Nlocal = Nrows;
   /* generate a unique numbering */

   start = ML_gpartialsum_int(Nrows, my_ml->comm);

   Nrecv = 0; Nsend = 0;
   for (i = 0; i < getrow_comm->N_neighbors; i++) {
      Nrecv += getrow_comm->neighbors[i].N_rcv;
      Nsend += getrow_comm->neighbors[i].N_send;
   }
   *Nghost = Nrecv;

   dmapper = (double *) ML_allocate((Nrows+Nrecv)*sizeof(double));
   for (i = 0; i < Nrows; i++) dmapper[i] = (double) (i+start+1);
   ML_exchange_bdry(dmapper, getrow_comm, Nrows,my_ml->comm,ML_OVERWRITE,NULL);

   mapper  = (int    *) ML_allocate((Nrows+Nrecv)*sizeof(int));
   for (i = 0; i < Nrows+Nrecv; i++) mapper[i] = (int) dmapper[i];
   ML_free(dmapper);

   /* generate a bunch of separators */

   N_bits = 0;
   i      = N_procs;
   while( i > 1) {
     i = i/2;
     N_bits++;
   }

   /* This is an over estimate of what we really need */

   total_rows = Nrows;
   ML_gsum_scalar_int(&total_rows, &j,my_ml->comm);
   sep_space = (int) (10.*sqrt( (double) total_rows ));
   sep_space += Nrows;

   sep     = (int *) ML_allocate(sep_space*sizeof(int));
   s_sizes = (int *) ML_allocate( (N_bits+2)*sizeof(int) );
   sub_mask = 0;
   N_count = 0;
   *Nseparators = N_bits + 1;
   old_Ncount = N_count;
   mask2 = 0;
   for (j = N_bits-1 ; j >= 0; j--) {
     mask = 1<<j;
     if ( (mask & proc) == 0) {
        for (i = 0; i < getrow_comm->N_neighbors; i++) {
           if (((mask2 & proc) == (mask2 & getrow_comm->neighbors[i].ML_id)) &
               ((mask & getrow_comm->neighbors[i].ML_id) != 0) ) {
              for (jj = 0; jj < getrow_comm->neighbors[i].N_send; jj++) {
                kk = getrow_comm->neighbors[i].send_list[jj];
                if (mapper[kk] > 0) {
                   sep[N_count++] = mapper[kk];
                   mapper[kk] = -mapper[kk];
                }
              }
           }
        }
     }
     mask2 += mask;
     i = N_count - old_Ncount;
     ML_Comm_subGappendInt(my_ml->comm, &(sep[old_Ncount]),&i,
                           sep_space-old_Ncount, sub_mask);


     ML_sort(i, &(sep[old_Ncount]) );

     sub_mask += mask;
     s_sizes[N_bits-1-j] = i;
     old_Ncount += i;
     N_count = old_Ncount;

  }

  /* take the domain interior as the last separator */

  N_count = old_Ncount;
  for (i = 0; i < Nrows; i++) {
      if (mapper[i] > 0) sep[N_count++] = mapper[i];
  }
  s_sizes[N_bits] = N_count - old_Ncount;


  for (i = 0; i < Nrows+Nrecv; i++) mapper[i] = ML_dabs(mapper[i]);

   /* generate a matrix-vector product (with sub communication) */

   N_nz = 0;
   allocated = 10;
   bindx = (int    *) ML_allocate(allocated*sizeof(int));
   val   = (double *) ML_allocate(allocated*sizeof(double));

   for (i = 0; i < Nrows; i++) {
      ML_get_matrix_row(Amat, 1, &i, &allocated, &bindx, &val,
                        &row_length, 0);
      N_nz += row_length;
   }
   ML_free(bindx); ML_free(val);

   bindx = (int    *) ML_allocate((N_nz+1)*sizeof(int));
   val   = (double *) ML_allocate((N_nz+1)*sizeof(double));
   row_ptr=(int    *) ML_allocate((Nrows+1)*sizeof(int));

   N_nz = 0;
   row_ptr[0] = N_nz;
   max_nz_per_row = 0;
   for (i = 0; i < Nrows; i++) {
      ML_get_matrix_row(Amat, 1, &i, &allocated, &bindx, &val,
                        &row_length, N_nz);
/*
if (proc == 0) {
for(kkk = 0; kkk < row_length; kkk++) 
   printf("A(%d,%d) = %e\n",i,bindx[N_nz+kkk],val[N_nz+kkk]);
}
*/
      N_nz += row_length;
      if (row_length > max_nz_per_row) max_nz_per_row = row_length;
      row_ptr[i+1] = N_nz;
   }
   temp = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata) );
   temp->columns       = bindx;
   temp->values        = val;
   temp->rowptr        = row_ptr;

   omatrix = ML_Operator_Create(my_ml->comm);
   omatrix->data_destroy = ML_CSR_MSRdata_Destroy;
   ML_Operator_Set_1Levels(omatrix, Amat->from, Amat->to);
   ML_Operator_Set_ApplyFuncData(omatrix, Amat->invec_leng,
                             Amat->getrow->Nrows, (void*)temp,
                             Amat->getrow->Nrows, NULL, 0);
   ML_Operator_Set_Getrow(omatrix, omatrix->getrow->Nrows,CSR_getrow);

   omatrix->max_nz_per_row = max_nz_per_row;
   omatrix->N_nonzeros     = N_nz;
   ML_Operator_Set_ApplyFunc (omatrix, CSR_matvec);


   neighbors = (int *) ML_allocate(sizeof(int)*getrow_comm->N_neighbors);
   for (i = 0; i < getrow_comm->N_neighbors; i++)
      neighbors[i] = getrow_comm->neighbors[i].ML_id;

   ML_CommInfoOP_Set_neighbors( &(omatrix->getrow->pre_comm), 
			        getrow_comm->N_neighbors, neighbors,
	                        ML_OVERWRITE, NULL, 0);
   ML_free(neighbors);
   for (i = 0; i < getrow_comm->N_neighbors; i++)
      ML_CommInfoOP_Set_exch_info(omatrix->getrow->pre_comm, 
			getrow_comm->neighbors[i].ML_id,
                        getrow_comm->neighbors[i].N_rcv,
                        getrow_comm->neighbors[i].rcv_list,
                        getrow_comm->neighbors[i].N_send,
                        getrow_comm->neighbors[i].send_list);



  *matvec_data = omatrix;
  *separator = sep;
  *sep_size  = s_sizes;
  *imapper   = mapper;

}

int ML_Comm_subGappendInt(ML_Comm *com_ptr, int *vals, int *cur_length, 
                    int total_length,int sub_mask)
{
   int     mask, partner, hbit, msgtype, msgbase=145;
   int     nbytes, mypid, nprocs, sub_cube, nprocs_small;
   USR_REQ Request;

  /*********************** first executable statment *****************/

  /* ----------------------------------------------------------------- */
  /* get processor information                                         */
  /* ----------------------------------------------------------------- */

  mypid    = com_ptr->ML_mypid;
  nprocs   = com_ptr->ML_nprocs;
  sub_cube = sub_mask & mypid;


  msgtype = msgbase+1;

  /* Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;

  if (nprocs_small * 2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  partner = mypid ^ nprocs_small;

  if (mypid + nprocs_small < nprocs) {

    /* post receives on the hypercube portion of the machine partition */

if ((sub_mask & partner) == sub_cube) {
    if ( com_ptr->USR_irecvbytes((void *) &(vals[*cur_length]),
                       (total_length - *cur_length) * sizeof(int), &partner,
                       &msgtype, com_ptr->USR_comm,
#ifdef ML_CPP
                                    &Request ))
#else
                                    (void *)&Request ))
#endif
    {
      (void) fprintf(stderr, "ERROR on node %d\nread failed, message type = %d\n", mypid, msgtype);
      exit(-1);
    }
}
  }
  else if (mypid & nprocs_small) {

    /*
     * Send messages from the portion of the machine partition "above" the
     * largest hypercube to the hypercube portion.
     */

if ((sub_mask & partner) == sub_cube) {
    if (com_ptr->USR_sendbytes((void *) vals, (*cur_length)*sizeof(int), 
			       partner, msgtype, com_ptr->USR_comm)) {
      (void) fprintf(stderr, "ERROR on node %d\nwrite failed, message type = %d\n", mypid, msgtype);
      exit(-1);
    }
}
  }

  if (mypid + nprocs_small < nprocs) {

    /* wait to receive the messages */

if ((sub_mask & partner) == sub_cube) {
    nbytes = com_ptr->USR_waitbytes((void *) &(vals[*cur_length]),
                          (total_length - *cur_length)*sizeof(int), &partner,
                          &msgtype, com_ptr->USR_comm,
#ifdef ML_CPP
                                    &Request );
#else
                                    (void *)&Request );
#endif
    (*cur_length) += (nbytes / sizeof(int));
}
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(mypid & nprocs_small)) {
    for (mask = nprocs_small >> 1; mask; mask >>= 1) {
      partner = mypid ^ mask;
if ((sub_mask & partner) == sub_cube) {
      if (com_ptr->USR_irecvbytes((void *) &(vals[*cur_length]),
                        (total_length - *cur_length)*sizeof(int), &partner,
                        &msgtype, com_ptr->USR_comm,
#ifdef ML_CPP
                                    &Request ))
#else
                                    (void *)&Request ))
#endif
      {
        (void) fprintf(stderr, "ERROR on node %d\nread failed, message type = %d\n", mypid, msgtype);
        exit(-1);
      }

      if (com_ptr->USR_sendbytes((void *) vals, *cur_length*sizeof(int), 
			partner, msgtype, com_ptr->USR_comm)) {
        (void) fprintf(stderr, "ERROR on node %d\nwrite failed, message type = %d\n", mypid, msgtype);
        exit(-1);
      }

      nbytes = com_ptr->USR_waitbytes((void *) &(vals[*cur_length]),
                            (total_length - *cur_length)*sizeof(int), &partner,
                            &msgtype, com_ptr->USR_comm,
#ifdef ML_CPP
                                    &Request );
#else
                                    (void *)&Request );
#endif
      (*cur_length) += (nbytes / sizeof(int));
}
    }
  }

  /* Finally, send message from lower half to upper half. */

  partner = mypid ^ nprocs_small;
  if (mypid & nprocs_small) {
if ((sub_mask & partner) == sub_cube) {
    if (com_ptr->USR_irecvbytes((void *) vals, total_length*sizeof(int),
			&partner,&msgtype,com_ptr->USR_comm,
#ifdef ML_CPP
                                    &Request ))
#else
                                    (void *)&Request ))
#endif
    {
      (void) fprintf(stderr, "ERROR on node %d\nread failed, message type = %d\n", mypid, msgtype);
      exit(-1);
    }
}
  }

  else if (mypid+nprocs_small < nprocs ) {
if ((sub_mask & partner) == sub_cube) {
    if (com_ptr->USR_sendbytes((void *) vals, *cur_length*sizeof(int), partner,
				msgtype, com_ptr->USR_comm )) {
      (void) fprintf(stderr, "ERROR on node %d\nwrite failed, message type = %d\n", mypid, msgtype);
      exit(-1);
    }
}
  }

  if (mypid & nprocs_small) {
if ((sub_mask & partner) == sub_cube) {
    nbytes = com_ptr->USR_waitbytes((void *) vals, total_length*sizeof(int), 
			  &partner,&msgtype,com_ptr->USR_comm,
#ifdef ML_CPP
                                    &Request );
#else
                                    (void *)&Request );
#endif
    (*cur_length) = (nbytes / sizeof(int));
}
  }
  return 0;

} /* ML_gappend */


int CSR_submv(ML_Operator *Amat, double p[], double ap[])
{
  int i, k, Nrows, *bindx, total_send, total_rcv;
  double *p2, *val, sum;
  struct ML_CSR_MSRdata *temp;
  ML_CommInfoOP *getrow_comm;
  int *row_ptr;

   Nrows = Amat->matvec->Nrows;
   temp  = (struct ML_CSR_MSRdata *) Amat->data;
   val   = temp->values;
   bindx = temp->columns;
   row_ptr= temp->rowptr;

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) {
      total_send = 0;
      total_rcv = 0;
      for (i = 0; i < getrow_comm->N_neighbors; i++) {
         total_send += getrow_comm->neighbors[i].N_send;
         total_rcv += getrow_comm->neighbors[i].N_rcv;
      }
      p2 = (double *) ML_allocate((Nrows+total_rcv+1)*sizeof(double));
      for (i = 0; i < Nrows+total_rcv; i++) p2[i] = p[i];
   }
   else p2 = p;

  for (i = 0; i < Nrows; i++) {
     sum = 0;
     for (k = row_ptr[i]; k < row_ptr[i+1]; k++)
        sum  += val[k] * p2[bindx[k]];
     ap[i] = sum;
  }

  if (getrow_comm != NULL) {
     for (i = 0; i < Nrows; i++) p[i] = p2[i];
     ML_free(p2);
  }
  return(1);
}
int ML_submv(ML_Operator *Amat, double p[], double ap[])
{
   int i, j, Nrows, total_send, total_rcv, *cols, allocated_space;
   int length, col;
   double *p2, *vals, dtemp;
   ML_CommInfoOP *getrow_comm;

   Nrows = Amat->matvec->Nrows;
   allocated_space = Amat->max_nz_per_row+2;
   cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
   vals = (double *) ML_allocate(allocated_space*sizeof(double));
   if (vals == NULL) pr_error("Error in ML_submatvec(): Not enough space\n");

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) {
      total_send = 0;
      total_rcv = 0;
      for (i = 0; i < getrow_comm->N_neighbors; i++) {
         total_send += getrow_comm->neighbors[i].N_send;
         total_rcv += getrow_comm->neighbors[i].N_rcv;
      }
      p2 = (double *) ML_allocate((Nrows+total_rcv+1)*sizeof(double));
      for (i = 0; i < Nrows+total_rcv; i++) p2[i] = p[i];
   }
   else p2 = p;

  for (i = 0; i < Nrows; i++) {
     ML_get_matrix_row(Amat, 1, &i , &allocated_space , &cols, &vals,
                           &length, 0);
     dtemp = 0.;
     for (j = 0; j < length; j++)
     {
        col = cols[j];
        dtemp += vals[j]*p2[col];
     }
     ap[i] = dtemp;


  }

  if (getrow_comm != NULL) {
     for (i = 0; i < Nrows; i++) p[i] = p2[i];
     ML_free(p2);
  }
  return(1);
}
int ML_submatvec(ML_Operator *Amat, double p[], double ap[], int mask)
{
   int i, j, Nrows, total_send, total_rcv, *cols, allocated_space;
   int col, length;
   double *p2, *vals, dtemp;
   ML_CommInfoOP *getrow_comm;

   Nrows = Amat->matvec->Nrows;
   allocated_space = Amat->max_nz_per_row+2;
   cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
   vals = (double *) ML_allocate(allocated_space*sizeof(double));
   if (vals == NULL) pr_error("Error in ML_submatvec(): Not enough space\n");

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) {
      total_send = 0;
      total_rcv = 0;
      for (i = 0; i < getrow_comm->N_neighbors; i++) {
         total_send += getrow_comm->neighbors[i].N_send;
         total_rcv += getrow_comm->neighbors[i].N_rcv;
      }
      p2 = (double *) ML_allocate((Nrows+total_rcv+1)*sizeof(double));
      for (i = 0; i < Nrows; i++) p2[i] = p[i];
      for (i = Nrows; i < Nrows+total_rcv; i++) p2[i] = 0.0;
      ML_subexchange_bdry(p2,getrow_comm, Nrows,total_send,Amat->to->comm,mask);
   }
   else p2 = p;

  for (i = 0; i < Nrows; i++) {
     ML_get_matrix_row(Amat, 1, &i , &allocated_space , &cols, &vals,
                           &length, 0);
     dtemp = 0.;
     for (j = 0; j < length; j++)
     {
        col = cols[j];
        dtemp += vals[j]*p2[col];
     }
     ap[i] = dtemp;


  }

  if (getrow_comm != NULL) {
     for (i = 0; i < Nrows; i++) p[i] = p2[i];
     ML_free(p2);
  }
  return(1);
}
int CSR_submatvec(ML_Operator *Amat, double p[], double ap[], int mask)
{
  int i, k, Nrows, *bindx, total_send, total_rcv;
   double *p2, *val, sum;
   struct ML_CSR_MSRdata *temp;
   ML_CommInfoOP *getrow_comm;
   int *row_ptr;

   Nrows = Amat->matvec->Nrows;
   temp  = (struct ML_CSR_MSRdata *) Amat->data;
   val   = temp->values;
   bindx = temp->columns;
   row_ptr= temp->rowptr;

   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) {
      total_send = 0;
      total_rcv = 0;
      for (i = 0; i < getrow_comm->N_neighbors; i++) {
         total_send += getrow_comm->neighbors[i].N_send;
         total_rcv += getrow_comm->neighbors[i].N_rcv;
      }
      p2 = (double *) ML_allocate((Nrows+total_rcv+1)*sizeof(double));
      for (i = 0; i < Nrows; i++) p2[i] = p[i];
      for (i = Nrows; i < Nrows+total_rcv; i++) p2[i] = 0.0;
      ML_subexchange_bdry(p2,getrow_comm, Nrows,total_send,Amat->to->comm,mask);
   }
   else p2 = p;

  for (i = 0; i < Nrows; i++) {
     sum = 0;
     for (k = row_ptr[i]; k < row_ptr[i+1]; k++)
        sum  += val[k] * p2[bindx[k]];
     ap[i] = sum;
  }

  if (getrow_comm != NULL) {
     for (i = 0; i < Nrows; i++) p[i] = p2[i];
     ML_free(p2);
  }
  return(1);
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void ML_subexchange_bdry(double x[], ML_CommInfoOP *comm_info,
                      int start_location, int total_send, ML_Comm *comm,
                      int mask)

/*******************************************************************************

  Routine to locally exchange components of the vector "x". This routine
  differs from ML_exchange_bdry() in that processors only send to neighbors
  whose node id has some subset of bits (indicated by mask) which match
  the processors node id.

  Author:          Ray Tuminaro, SNL, 9222
  =======

  Parameter list:
  ===============

  x                On input, vector of unknowns defined on current processor.
                   On output, x is appended with information received from
                   other processors (dictated by 'comm_info').

  comm_info        On input, comm_info contains information on what needs to be
                   exchanged with other processors. See ml_rap.h.

  start_location   On input, starting location in 'x' where received information
                   will be placed.

  total_send       On input, total number of components in 'x' to be sent.

*******************************************************************************/

{
  register double *ptrd;
  double          *ptr_send_list, *ptr_recv_list, *orig_ptr;
  int              type, N_neighbors, *temp, i, j, rtype;
  int             proc, sub_proc;
  USR_REQ         *request;
  ML_NeighborList *neighbor;

  /**************************** execution begins ******************************/
  N_neighbors              = comm_info->N_neighbors;
  proc  = comm->ML_mypid;
  sub_proc = proc & mask;
  if (N_neighbors == 0) return;

  /* Set up send messages: Gather send unknowns from "x" vector */

  request = (USR_REQ     *)  ML_allocate(N_neighbors*sizeof(USR_REQ    ));
  ptrd = (double *) ML_allocate( total_send*sizeof(double));
  if (ptrd == NULL) {
     printf("Out of space in ML_exchange_bdry\n");
     exit(1);
  }
  ptr_send_list = ptrd;
  orig_ptr      = ptrd;


  for (i = 0; i < N_neighbors; i++) {
     temp = comm_info->neighbors[i].send_list;
     for (j = 0; j < comm_info->neighbors[i].N_send; j++) {
        *ptrd++ = x[ temp[j] ];
     }
  }

  type = 1991;

  /* post receives for all messages */

  ptr_recv_list = &x[start_location];
  for (i = 0; i < N_neighbors; i++) {
    neighbor = &(comm_info->neighbors[i]);
    rtype = type;   j = sizeof(double)* neighbor->N_rcv;
    if ((neighbor->ML_id & mask) == sub_proc ) {
       comm->USR_irecvbytes((void *) ptr_recv_list, j, &(neighbor->ML_id), &rtype,
                         comm->USR_comm, request+i);
    }
    ptr_recv_list         += neighbor->N_rcv;
  }

  /* write out all messages */

  for (i = 0; i < N_neighbors; i++) {
    neighbor = &(comm_info->neighbors[i]);
    j = sizeof(double)* neighbor->N_send;
    if ((neighbor->ML_id & mask) == sub_proc ) {
       comm->USR_sendbytes((void *) ptr_send_list, j, neighbor->ML_id,
                          rtype, comm->USR_comm);
    }
    ptr_send_list         += neighbor->N_send;
  }

  /* wait for all messages */

  ptr_recv_list = &x[start_location];
  for (i = 0; i < N_neighbors; i++) {
    neighbor = &(comm_info->neighbors[i]);
    rtype = type;   j = sizeof(double)* neighbor->N_rcv;
    if ((neighbor->ML_id & mask) == sub_proc ) {
       comm->USR_cheapwaitbytes((void *) ptr_recv_list, j, &(neighbor->ML_id),
                        &rtype, comm->USR_comm, request+i);
    }
    ptr_recv_list         += neighbor->N_rcv;
  }

  ML_free(orig_ptr);
  ML_free(request);


} /* ML_exchange_bdry */

