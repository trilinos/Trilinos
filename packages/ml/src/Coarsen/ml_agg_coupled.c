/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person    */
/* and disclaimer.                                                           */
/* ************************************************************************* */

/* ************************************************************************* */
/* ************************************************************************* */
/* Functions to create tentative prolongators (coupled aggregation)          */
/* ************************************************************************* */
/* Author        : Charles Tong (LLNL)                                       */
/* Date          : August, 2000                                              */
/* ************************************************************************* */
/* Local Functions :                                                         */
/*    ML_Aggregate_CoarsenCoupled                                            */
/*    ML_Aggregate_CoarsenCoupledCore                                        */
/*    ML_Aggregate_Compress_Matrix                                           */
/*    ML_Aggregate_ExchangeStatus                                            */
/*    ML_Aggregate_ComposeExpandedCommInfo                                   */
/*    ML_Aggregate_ComposeRecvFromSend                                       */
/*    ML_Aggregate_Form_Aggregates                                           */
/*    ML_Aggregate_PutInto_Aggregates                                        */
/* ************************************************************************* */

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ml_aggregate.h"
#include "ml_lapack.h"
#include "ml_utils.h"

/* ************************************************************************* */
/* local defines                                                             */
/* ------------------------------------------------------------------------- */

#define ML_AGGR_READY      -11
#define ML_AGGR_NOTSEL     -12
#define ML_AGGR_SELECTED   -13
#define ML_AGGR_SELECTED2  -14
#define ML_AGGR_BDRY       -15
#define ML_AGGR_MINRANK      1
#define ML_AGGR_MAXLINK      2
#define ML_AGGR_INTERIOR     0
#define ML_AGGR_BORDER       1

/* ************************************************************************* */
/* ML_Aggregate_CoarsenCoupled subroutine.                                   */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_CoarsenCoupled( ML_Aggregate *ml_ag,
       ML_Operator *Amatrix, ML_Operator **Pmatrix, ML_Comm *comm)
{
   int     i, j, k, m, jj, index, index3, index4, offset, count; 
   int     max_count, nbytes, length, level, diff_level;
   int     Nrows, exp_Nrows, *mat_indx=NULL, *amal_mat_indx;
   int     N_neighbors, *neighbors, *recv_leng, *send_leng, *send_list;
   int     total_recv_leng, total_send_leng, msgtype, mypid, new_N_send;
   int     *new_send_neighbors, *new_send_list, *new_send_leng;
   int     new_N_recv, *new_recv_leng, *new_recv_neighbors, *int_buf;
   int     *int_buf2, *recv_list, nprocs;
   double  printflag;
   int     aggr_count, *aggr_index, *aggr_index2;
   int     *aggr_cnt_array, max_agg_size, **rows_in_aggs;
   int     Ncoarse, exp_Ncoarse, *new_ia, *new_ja, new_Nrows;
   int     num_PDE_eqns, nullspace_dim, lwork, info, *bc_array;
   double  dcompare1, *new_val=NULL;
   double  epsilon, *dble_buf=NULL, *nullspace_vect=NULL, *qr_tmp=NULL;
   double  *tmp_vect=NULL, *work=NULL, *new_null=NULL, *comm_val=NULL;
   double  *dble_buf2=NULL, largest, thesign, dtemp;
   struct ML_CSR_MSRdata *csr_data;
   ML_Aggregate_Comm     *aggr_comm;
   ML_GetrowFunc         *getrow_obj;

   /* ============================================================= */
   /* get machine and matrix information                            */
   /* ============================================================= */

   mypid          = comm->ML_mypid;
   nprocs         = comm->ML_nprocs;
   nullspace_dim  = ml_ag->nullspace_dim;
   nullspace_vect = ml_ag->nullspace_vect;
   printflag      = ml_ag->print_flag;
   num_PDE_eqns   = ml_ag->num_PDE_eqns;
   Nrows          = Amatrix->outvec_leng;
   getrow_obj     = Amatrix->getrow;

   /* ============================================================= */
   /* initialize and update the threshold                           */
   /* ============================================================= */

   diff_level = ml_ag->begin_level - ml_ag->cur_level;
   if ( diff_level == 0 ) ml_ag->curr_threshold = ml_ag->threshold;
   epsilon = ml_ag->curr_threshold;
   ml_ag->curr_threshold *= 0.5;

   if ( mypid == 0 && printflag < ML_Get_PrintLevel())
   {
      printf("ML_Aggregate_CoarsenCoupled : current level = %d\n",
                           ml_ag->cur_level);
      printf("ML_Aggregate_CoarsenCoupled : current eps = %e\n",epsilon);
   }
   epsilon = epsilon * epsilon;

   /* ============================================================= */
   /* form a matrix graph after taking out weak edges -> mat_indx   */
   /* ============================================================= */

   ML_Graph_CreateFromMatrix(ml_ag,Amatrix,&mat_indx,comm,epsilon,&exp_Nrows,
                             &bc_array);

   /* ============================================================= */
   /* compress the matrix using num_PDE_eqns information            */
   /* ============================================================= */

   ML_Aggregate_Compress_Matrix(getrow_obj,mat_indx,num_PDE_eqns, 
                comm,&amal_mat_indx,&N_neighbors,&neighbors,&recv_leng,
                &send_leng, &send_list, &recv_list, bc_array);
   if ( mat_indx != amal_mat_indx ) ML_memory_free((void**) &mat_indx);
   mat_indx = NULL;

   /* ============================================================= */
   /* perform coarsening on the compressed matrix                   */
   /* ============================================================= */

   ML_Aggregate_CoarsenCoupledCore(ml_ag, comm, amal_mat_indx,
       &aggr_count, &aggr_index2, N_neighbors, neighbors, recv_leng,
       send_leng, send_list, recv_list, &aggr_cnt_array, bc_array);
   if ( amal_mat_indx != NULL ) ML_memory_free( (void**) &amal_mat_indx );
   if ( neighbors     != NULL ) ML_memory_free( (void**) &neighbors );
   if ( recv_leng     != NULL ) ML_memory_free( (void**) &recv_leng );
   if ( send_leng     != NULL ) ML_memory_free( (void**) &send_leng );
   if ( send_list     != NULL ) ML_memory_free( (void**) &send_list );
   if ( recv_list     != NULL ) ML_memory_free( (void**) &recv_list );
   if ( bc_array      != NULL ) ML_memory_free( (void**) &bc_array );

   /* ============================================================= */
   /* compose a new communication info for A in view of num_PDEs    */
   /* This is needed to send null spaces around between processors) */
   /* ============================================================= */

   ML_Aggregate_ComposeExpandedCommInfo(getrow_obj, num_PDE_eqns, 
           comm, &N_neighbors, &neighbors, &recv_leng, &send_leng, 
            &send_list, &recv_list);

   /* ============================================================= */
   /* decode aggr_index2 to recover the block information           */
   /* after this, nodes that are aggregated locally are given a     */
   /* non-negative integer, and nodes that are not aggregated       */
   /* locally are given a negative number with PID on it.           */
   /* ============================================================= */

   exp_Nrows = Nrows;
   for ( i = 0; i < N_neighbors; i++ ) exp_Nrows += recv_leng[i];
   
   for ( i = 0; i < aggr_count; i++ ) aggr_cnt_array[i] = 0;
   nbytes = exp_Nrows * sizeof(int);
   if ( nbytes > 0 ) ML_memory_alloc((void**) &aggr_index, nbytes, "ACJ");
   else              aggr_index = NULL;
   for ( i = 0; i < exp_Nrows; i++ )
   {
      aggr_index[i] = aggr_index2[i/num_PDE_eqns];
      if ( aggr_index[i] >= 0 && aggr_index[i] < aggr_count )
          aggr_cnt_array[aggr_index[i]]++;
   }   
   if ( aggr_index2 != NULL ) ML_memory_free( (void**) &aggr_index2 );
   aggr_index2 = NULL;

   /* ============================================================= */
   /* Form tentative prolongator                                    */
   /* ============================================================= */

   Ncoarse = aggr_count;

   /* ------------------------------------------------------------- */
   /* send the remote node index back to remote processors, with    */
   /* added information on which remote nodes have been aggregated  */
   /* by the local aggregates (and also the aggregate numbers).     */
   /* (in aggr_index[0:Nrows-1], a postive integer means this node  */
   /* has been aggregated locally; while a negative integer (-1XX)  */
   /* means it has been aggregated by remote processors.  As to     */
   /* aggr_index[Nrows:exp_Nrows-1], a positive integer means the   */
   /* remote node has been aggregated locally, and negative         */
   /* otherwise) ===> int_buf                                       */
   /* ------------------------------------------------------------- */

   total_send_leng = 0;
   for ( i = 0; i < N_neighbors; i++ ) total_send_leng += send_leng[i];
   total_recv_leng = 0;
   for ( i = 0; i < N_neighbors; i++ ) total_recv_leng += recv_leng[i];

   nbytes = total_send_leng * sizeof(int);
   if ( nbytes > 0 ) ML_memory_alloc((void**)&int_buf, nbytes, "Ag1");
   else              int_buf = NULL;
   nbytes = total_recv_leng * sizeof(int);
   if ( nbytes > 0 ) ML_memory_alloc((void**)&int_buf2, nbytes, "Ag2");
   else              int_buf2 = NULL;

   offset = 0;
   for ( i = 0; i < N_neighbors; i++ )
   {
      for ( j = 0; j < recv_leng[i]; j++ )
      {
         if (aggr_index[recv_list[offset+j]]<0) int_buf2[offset+j] = -1;
         else int_buf2[offset+j] = aggr_index[recv_list[offset+j]];
      }
      offset += recv_leng[i];
   }
   msgtype = 15963;
   ML_Aggregate_ExchangeStatus((char*) int_buf, (char*) int_buf2,
      N_neighbors,neighbors,send_leng,recv_leng,NULL,Nrows,msgtype, 
      ML_INT, comm);

   if ( int_buf2 != NULL ) ML_memory_free( (void**) &int_buf2 );

   /* ------------------------------------------------------------- */
   /* if int_buf[i] > 0, this means that aggr_index[send_list[i]]   */
   /* has been aggregated by aggregate int_buf[i] of a remote       */
   /* processor.  The next step uses this information to update     */
   /* aggr_index (from a negative integer for remote aggregate to   */
   /* a positive local index for remote aggregate) and exp_Ncoarse  */
   /* ------------------------------------------------------------- */

   offset = 0;
   m      = 0; /* store the local index offset for remote processors */
   for ( i = 0; i < N_neighbors; i++ )
   {
      /* ---------------------------------------------------------- */
      /* find out how large an array to allocate for int_buf2,      */
      /* which is used to count the number of distinct aggregates   */
      /* remote processor neighbor[i] used on my local nodes        */
      /* ---------------------------------------------------------- */

      max_count = -1;
      for ( j = 0; j < send_leng[i]; j++ )
      {
         index = int_buf[offset+j];
         max_count = (index > max_count ) ? index : max_count;
      }
      nbytes = ( max_count + 2 ) * sizeof(int);
      if ( nbytes > 0 ) ML_memory_alloc((void**)&int_buf2, nbytes, "Ag3");
      else              int_buf2 = NULL;

      /* ---------------------------------------------------------- */
      /* see how many distinct remote aggregates are referenced by  */
      /* local fine nodes in aggregation in neighbor[i]             */
      /* ---------------------------------------------------------- */

      for ( j = 0; j <= max_count; j++ ) int_buf2[j] = 0;
      for ( j = 0; j < send_leng[i]; j++ )
      {
         index = int_buf[offset+j];
         if ( index >= 0 ) int_buf2[index]++;
         if (index >= 0 && index > max_count)
            {printf("int_buf2 error : maxcount\n");exit(1);}
      }
      count = 0;
      for ( j = 0; j <= max_count; j++ )
      {
         if (int_buf2[j] > 0)
         {
            count++; int_buf2[j] = 1;
         }
      }
      for ( j = max_count; j > 0; j-- ) int_buf2[j] = int_buf2[j-1];
      int_buf2[0] = 0;
      for ( j = 0; j < max_count; j++ ) int_buf2[j+1] += int_buf2[j];

      /* ---------------------------------------------------------- */
      /* now assign local aggregate indices to local nodes that are */
      /* aggregated by remote processors                            */
      /* ---------------------------------------------------------- */

      for ( j = 0; j < send_leng[i]; j++ )
      {
         index = send_list[offset+j];

         /* ------------------------------------------------------- */
         /* The first condition indicates that the local node has   */
         /* been registered to have been aggregated by remote       */
         /* aggregates.  The second condition is needed in case     */
         /* the local node is linked to more than 1 remote          */
         /* processor (but only to one aggregate though)            */
         /* int_buf2 contains local indices of remote aggregates    */
         /* ------------------------------------------------------- */

         if ( aggr_index[index] <= -100 && int_buf[offset+j] >= 0 )
         {
            k = int_buf[offset+j];
            aggr_index[index] = int_buf2[k] + Ncoarse + m;
         }
      }
      if ( int_buf2 != NULL ) ML_memory_free((void**)&int_buf2);
      m += count;
      offset += send_leng[i];
   }
   exp_Ncoarse = Ncoarse + m;

   if ( int_buf != NULL ) ML_memory_free( (void**) &int_buf );

   /* ============================================================= */
   /* check and copy aggr_index to ml_ag (for block smoothing)      */
   /* ------------------------------------------------------------- */

   level = ml_ag->cur_level;
   nbytes = Nrows * sizeof( int );
   if ( nbytes > 0 )
      ML_memory_alloc((void**) &(ml_ag->aggr_info[level]), nbytes, "ACK");
   else
      ml_ag->aggr_info[level] = NULL;
   count = aggr_count;
   for ( i = 0; i < Nrows; i++ )
   {
      ml_ag->aggr_info[level][i] = aggr_index[i];
      if (aggr_index[i] >= count) count = aggr_index[i] + 1;
   }
   ml_ag->aggr_count[level] = count; /* for relaxing boundary points */

   /* ============================================================= */
   /* find out how many local coarse aggregates are needed by       */
   /* remote processors for interpolation (to construct the         */
   /* communicator - send info - for P)                             */
   /* ------------------------------------------------------------- */

   new_N_send = new_N_recv = 0;
   nbytes = N_neighbors * sizeof(int);
   if ( nbytes > 0 ) int_buf = (int *) ML_allocate( nbytes );
   else              int_buf = NULL;
   nbytes = Ncoarse * sizeof(int);
   if ( nbytes > 0 ) int_buf2 = (int *) ML_allocate( nbytes );
   else              int_buf2 = NULL;
   for ( i = 0; i < N_neighbors; i++ ) int_buf[i] = 0;

   /* ------------------------------------------------------------- */
   /* count which remote fine nodes belong to local aggregates      */
   /* in order to generate the communication pattern for            */
   /* the interpolation operator.                                   */
   /* ------------------------------------------------------------- */

   offset = 0;
   for ( i = 0; i < N_neighbors; i++ )
   {
      for ( j = 0; j < Ncoarse; j++ ) int_buf2[j] = 0;
      for ( j = 0; j < recv_leng[i]; j++ )
      {
         index = aggr_index[recv_list[offset++]];
         if ( index >= 0 ) int_buf2[index]++;
      }
      count = 0;
      for (j = 0; j < Ncoarse; j++) if ( int_buf2[j] > 0 ) count++;
      int_buf[i] = count * nullspace_dim;
      if ( int_buf[i] > 0 ) new_N_send++;
   }

   /* ------------------------------------------------------------- */
   /* now the number of neighbors for P has been found, the next    */
   /* step is to find the send_list and send_leng for the matvec    */
   /* function for interpolation                                    */
   /* ------------------------------------------------------------- */

   nbytes = new_N_send * sizeof(int);
   if ( nbytes > 0 )
   {
      ML_memory_alloc((void**) &new_send_leng, nbytes, "ACL");
      ML_memory_alloc((void**) &new_send_neighbors, nbytes, "ACM");
      new_N_send = 0;
      for ( i = 0; i < N_neighbors; i++ )
      {
         if ( int_buf[i] > 0 )
         {
            new_send_leng[new_N_send] = int_buf[i];
            new_send_neighbors[new_N_send] = neighbors[i];
            new_N_send++;
         }
      }
      count = 0;
      for ( i = 0; i < new_N_send; i++ ) count += new_send_leng[i];
      nbytes = count * sizeof(int);
      ML_memory_alloc((void**) &new_send_list, nbytes, "ACN");
      offset = 0;
      m = count;
      count = 0;
      for ( i = 0; i < N_neighbors; i++ )
      {
         for ( j = 0; j < Ncoarse; j++ ) int_buf2[j] = 0;
         for ( j = 0; j < recv_leng[i]; j++ )
         {
            index = aggr_index[recv_list[offset++]];
            if ( index >= 0 ) int_buf2[index]++;
         }
         for ( j = 0; j < Ncoarse; j++ )
         {
            if ( int_buf2[j] > 0 )
            {
               for ( jj = 0; jj < nullspace_dim; jj++ )
                  new_send_list[count++] = j * nullspace_dim + jj;
            }
         }
      }
      if ( m != count )
      {
         printf("ML_Aggregate_Coupled : internal error (1).\n");
         exit(-1);
      }
   }
   else
   {
      new_send_leng = NULL;
      new_send_neighbors = NULL;
      new_send_list = NULL;
   } 
   if ( int_buf  != NULL ) ML_free (int_buf);
   if ( int_buf2 != NULL ) ML_free (int_buf2);
   ML_Aggregate_ComposeRecvFromSend(nprocs, mypid, new_N_send, new_send_leng, 
          new_send_neighbors, &new_N_recv, &new_recv_leng, &new_recv_neighbors, 
          comm);

   /* ============================================================= */
   /* set up the new operator                                       */
   /* ------------------------------------------------------------- */

   new_Nrows = Nrows;
   for ( i = 0; i < new_Nrows; i++ )
   {
      if ( aggr_index[i] >= exp_Ncoarse )
         printf("WARNING : index out of bound %d = %d(%d)\n",i,aggr_index[i],
                exp_Ncoarse);
   }
   nbytes = ( new_Nrows + 1 ) * sizeof(int);
   ML_memory_alloc((void**)&(new_ia), nbytes, "ACO");
   nbytes = new_Nrows * nullspace_dim * sizeof(int);
   if ( nbytes > 0 ) ML_memory_alloc((void**)&(new_ja), nbytes, "ACP");
   else              new_ja = NULL;
   nbytes = new_Nrows * nullspace_dim * sizeof(double);
   if ( nbytes > 0 ) ML_memory_alloc((void**)&(new_val), nbytes, "ACQ");
   else              new_ja = NULL;
   for ( i = 0; i < new_Nrows*nullspace_dim; i++ ) new_val[i] = 0.0;
   for ( i = 0; i < new_Nrows*nullspace_dim; i++ ) new_ja[i] = 0;

   /* ------------------------------------------------------------- */
   /* set up the space for storing the new null space               */
   /* ------------------------------------------------------------- */

   nbytes = Ncoarse * nullspace_dim * nullspace_dim * sizeof(double);
   if ( nbytes > 0 )
      ML_memory_alloc((void**)&(new_null),nbytes,"ACR");
   else
      new_null = NULL;
   for (i = 0; i < Ncoarse*nullspace_dim*nullspace_dim; i++)
      new_null[i] = 0.0;

   /* ------------------------------------------------------------- */
   /* initialize the row pointer for the CSR prolongation operator  */
   /* (each row will have at most nullspace_dim nonzero entries)    */
   /* ------------------------------------------------------------- */

   for (i = 0; i <= Nrows; i++) new_ia[i] = i * nullspace_dim;

   /* ------------------------------------------------------------- */
   /* generate an array to store which aggregate has which rows.Then*/
   /* loop through the rows of A checking which aggregate each row  */
   /* is in, and adding it to the appropriate spot in rows_in_aggs  */
   /* ------------------------------------------------------------- */

   if ( aggr_count > 0 )
      ML_memory_alloc((void**)&rows_in_aggs,aggr_count*sizeof(int*),"ACS");
   else
      rows_in_aggs = NULL;

   for (i = 0; i < aggr_count; i++)
   {
      rows_in_aggs[i] = (int *) ML_allocate(aggr_cnt_array[i]*sizeof(int));
      aggr_cnt_array[i] = 0;
      if (rows_in_aggs[i] == NULL)
      {
         printf("ERROR: couldn't allocate memory in CoarsenCoupled\n");
         printf("       requested = %d\n",aggr_cnt_array[i]);
         exit(1);
      }
   }
   for (i = 0; i < exp_Nrows; i++)
   {
      if ( aggr_index[i] >= 0 && aggr_index[i] < aggr_count)
      {
         index = aggr_cnt_array[aggr_index[i]]++;
         rows_in_aggs[aggr_index[i]][index] = i;
      }
   }

   /* ------------------------------------------------------------- */
   /* allocate work arrays for QR factorization                     */
   /* work and lwork are needed for lapack's QR routine.  These     */
   /* settings seemed easiest since I don't quite understand        */
   /* what they do, but may want to do something better here later  */
   /* ------------------------------------------------------------- */

   nbytes = total_recv_leng * nullspace_dim * sizeof(double);
   if ( nbytes > 0 ) ML_memory_alloc((void**)&comm_val, nbytes, "ACT");
   else              comm_val = NULL;
   for (i = 0; i < total_recv_leng*nullspace_dim; i++) comm_val[i] = 0.0;
   max_agg_size = 0;
   for (i = 0; i < aggr_count; i++)
   {
      if (aggr_cnt_array[i] > max_agg_size) max_agg_size = aggr_cnt_array[i];
   }
   nbytes = max_agg_size * nullspace_dim * sizeof(double);
   if ( nbytes > 0 ) ML_memory_alloc((void**)&qr_tmp, nbytes, "ACU");
   else              qr_tmp = NULL;
   nbytes = nullspace_dim * sizeof(double);
   if ( nbytes > 0 ) ML_memory_alloc((void**)&tmp_vect, nbytes, "ACV");
   else              tmp_vect = NULL;

   lwork  = nullspace_dim;
   nbytes = nullspace_dim * sizeof(double);
   if ( nbytes > 0 ) ML_memory_alloc((void**)&work, nbytes, "ACW");
   else              work = NULL;

   /* ------------------------------------------------------------- */
   /* ship the null space information to other processors           */
   /* ------------------------------------------------------------- */

   if (nullspace_vect != NULL)
   {
      nbytes = total_send_leng * nullspace_dim * sizeof(double);
      if ( nbytes > 0 ) ML_memory_alloc((void**) &dble_buf, nbytes,"ACX");
      else              dble_buf = NULL;
      nbytes = total_recv_leng * nullspace_dim * sizeof(double);
      if ( nbytes > 0 ) ML_memory_alloc((void**) &dble_buf2, nbytes,"ACY");
      else              dble_buf2 = NULL;
      length = total_send_leng * nullspace_dim;
      for ( i = 0; i < total_send_leng; i++ )
      {
         index = send_list[i];
         for ( j = 0; j < nullspace_dim; j++ )
            dble_buf[i*nullspace_dim+j] = nullspace_vect[j*Nrows+index];
      }
      msgtype = 12093;
      length = sizeof(double) * nullspace_dim;
      ML_Aggregate_ExchangeStatus((char*)dble_buf2,(char*) dble_buf,
            N_neighbors, neighbors, recv_leng, send_leng, recv_list, 
            Nrows, msgtype,length,comm);
      if ( dble_buf != NULL ) ML_memory_free((void**) &dble_buf);
   } 
   else dble_buf2 = NULL; 

   /* ------------------------------------------------------------- */
   /* perform block QR decomposition                                */
   /* ------------------------------------------------------------- */

   work[0] = lwork;
   for (i = 0; i < aggr_count; i++)
   {
      /* ---------------------------------------------------------- */
      /* set up the matrix we want to decompose into Q and R:       */
      /* ---------------------------------------------------------- */

      length = aggr_cnt_array[i];
      if (nullspace_vect == NULL)
      {
         for (j = 0; j < length; j++)
         {
            index = rows_in_aggs[i][j];
            count = num_PDE_eqns;
            index = 0;
            for (k = 0; k < nullspace_dim; k++)
            {
               if (index % count == k) qr_tmp[k*length+j] = 1.0;
               else                    qr_tmp[k*length+j] = 0.0;
               index++;
            }
         }
      }
      else
      {
         for (k = 0; k < nullspace_dim; k++)
         {
            for (j = 0; j < length; j++)
            {
               index = rows_in_aggs[i][j];
               if (index < Nrows)
               {
                  qr_tmp[k*length+j] = nullspace_vect[k*Nrows+index];
               }
               else
               {
                  qr_tmp[k*length+j] =
                        dble_buf2[(index-Nrows)*nullspace_dim+k];
               }
            }
         }
      }

      /* ---------------------------------------------------------- */
      /* now calculate QR using an LAPACK routine                   */
      /* ---------------------------------------------------------- */

      if ( nullspace_dim == 1 )
      {
         dtemp = 0.0;
         for (j = 0; j < aggr_cnt_array[i]; j++)
            dtemp += ( qr_tmp[j] * qr_tmp[j] ); 
         dtemp = sqrt( dtemp );
         tmp_vect[0] = qr_tmp[0];
         qr_tmp[0] = dtemp;
      }
      else
      {
         DGEQRF_F77(&(aggr_cnt_array[i]), &nullspace_dim, qr_tmp,
                  &(aggr_cnt_array[i]), tmp_vect, work, &lwork, &info);
         if (info != 0)
            pr_error("CoarsenCoupled ERROR : dgeqrf returned a non-zero\n");
      }

      if (work[0] > lwork)
      {
         lwork=(int) work[0];
         ML_memory_free((void**) &work);
         ML_memory_alloc((void**) &work, sizeof(double)*lwork, "ACa");
      }
      else lwork=(int) work[0];

      /* ---------------------------------------------------------- */
      /* the upper triangle of qr_tmp is now R, so copy that into   */
      /* the new nullspace                                          */
      /* ---------------------------------------------------------- */

      for (j = 0; j < nullspace_dim; j++)
         for (k = j; k < nullspace_dim; k++)
            new_null[i*nullspace_dim+j+k*Ncoarse*nullspace_dim] =
               qr_tmp[j+aggr_cnt_array[i]*k];

      /* ---------------------------------------------------------- */
      /* to get this block of P, need to run qr_tmp through another */
      /* LAPACK function:                                           */
      /* ---------------------------------------------------------- */

      if ( aggr_cnt_array[i] < nullspace_dim )
         printf("ERROR : performing QR on a MxN matrix where M<N.\n");

      if ( nullspace_dim == 1 )
      {
         dtemp = qr_tmp[0];
         qr_tmp[0] = tmp_vect[0];
         dtemp = 1.0 / dtemp;
         for (j = 0; j < aggr_cnt_array[i]; j++)
            qr_tmp[j] *= dtemp;
      }
      else
      {
         DORGQR_F77(&(aggr_cnt_array[i]),&nullspace_dim,&nullspace_dim,
                 qr_tmp, &(aggr_cnt_array[i]), tmp_vect, work, &lwork, &info);
         if (info != 0)
            pr_error("CoarsenCoupled ERROR : dorgqr returned a non-zero\n");
      }

      if (work[0] > lwork)
      {
         lwork=(int) work[0];
         ML_memory_free((void**) &work);
         ML_memory_alloc((void**) &work, sizeof(double)*lwork, "ACb");
      }
      else lwork=(int) work[0];

      /* ---------------------------------------------------------- */
      /* now copy Q over into the appropriate part of P:            */
      /* The rows of P get calculated out of order, so I assume the */
      /* Q is totally dense and use what I know of how big each Q   */
      /* will be to determine where in ia, ja, etc each nonzero in  */
      /* Q belongs.  If I did not assume this, I would have to keep */
      /* all of P in memory in order to determine where each entry  */
      /* should go                                                  */
      /* ---------------------------------------------------------- */

      largest = 0.; thesign = 1.;
      for (j = 0; j < aggr_cnt_array[i]; j++)
      {
         index = rows_in_aggs[i][j];
         for (k = 0; k < nullspace_dim; k++)
         {
            dtemp = qr_tmp[ k*aggr_cnt_array[i]+j];
            if ( ML_dabs( dtemp ) > largest )
            {
               largest = ML_dabs( dtemp );
               if ( dtemp < 0.0) thesign = -1.;
               else thesign = 1.;
            }
         }
      }
      for (j = 0; j < aggr_cnt_array[i]; j++)
      {
         index = rows_in_aggs[i][j];
         if ( index < Nrows )
         {
            index3 = new_ia[index];
            for (k = 0; k < nullspace_dim; k++)
            {
               new_ja [index3+k] = i * nullspace_dim + k;
               new_val[index3+k] = thesign * qr_tmp[ k*aggr_cnt_array[i]+j];
            }
         }
         else
         {
            index3 = (index - Nrows) * nullspace_dim;
            for (k = 0; k < nullspace_dim; k++)
               comm_val[index3+k] = thesign * qr_tmp[ k*aggr_cnt_array[i]+j];
         }
      }
   }

   ML_Aggregate_Set_NullSpace(ml_ag, num_PDE_eqns, nullspace_dim,
                              new_null, Ncoarse*nullspace_dim);
   
   ML_memory_free( (void **) &new_null);
   if (dble_buf2 != NULL) ML_memory_free( (void **) &dble_buf2);

   /* ------------------------------------------------------------- */
   /* send the P rows back to its parent processor                  */
   /* ------------------------------------------------------------- */

   nbytes = total_send_leng * nullspace_dim * sizeof(double);
   if ( nbytes > 0 ) ML_memory_alloc((void**) &dble_buf, nbytes,"ACc");
   else              dble_buf = NULL;
   msgtype = 24945;
   length = sizeof(double) * nullspace_dim;
   ML_Aggregate_ExchangeStatus((char*)dble_buf,(char*) comm_val,
         N_neighbors, neighbors, send_leng, recv_leng, NULL, 
         Nrows, msgtype,length,comm);
   for ( i = 0; i < total_send_leng; i++ )
   {
      index = send_list[i];
      if ( aggr_index[index] >= aggr_count )
      {
         dcompare1 = 0.0;
         for ( j = 0; j < nullspace_dim; j++ )
         {
            index4 = i * nullspace_dim + j;
            dcompare1 += dble_buf[index4];
         }
         if ( dcompare1 != 0.0 )
         {
            index4 = i * nullspace_dim;
            k      = index * nullspace_dim;
            for ( j = 0; j < nullspace_dim; j++ )
            {
               new_val[k+j] = dble_buf[index4+j];
               new_ja[k+j]  = aggr_index[index]*nullspace_dim+j;
            }
         }
      }
   }
   if ( comm_val != NULL ) ML_memory_free( (void **) &comm_val);
   if ( dble_buf != NULL ) ML_memory_free( (void **) &dble_buf);

   /* ------------------------------------------------------------- */
   /* check P (row sum = 1)                                         */
   /* ------------------------------------------------------------- */

   for ( i = 0; i < Nrows; i++ )
   {
      dcompare1 = 0.0;
      for (j = new_ia[i]; j < new_ia[i+1]; j++)
      {
         dcompare1 += new_val[j];
      }
/*
      printf("%d : CoarsenCoupled : rowsum(P(%d)) = %e (%d)\n",
             mypid, i, dcompare1, aggr_index[i]);
*/
      if ( dcompare1 == 0.0 && aggr_index[i] >= 0 )
         printf("%d : CoarsenCoupled WARNING : rowsum(P(%d)) = 0 (%d)\n",
                 mypid, i, aggr_index[i]);
   }

   /* ------------------------------------------------------------- */
   /* set up the csr_data data structure                            */
   /* ------------------------------------------------------------- */

   ML_memory_alloc((void**)&csr_data,sizeof(struct ML_CSR_MSRdata),"CSR");
   csr_data->rowptr  = new_ia;
   csr_data->columns = new_ja;
   csr_data->values  = new_val;
   ML_Operator_Set_ApplyFuncData(*Pmatrix, nullspace_dim*Ncoarse, Nrows,
                                 csr_data, Nrows, NULL, 0);
   (*Pmatrix)->data_destroy = ML_CSR_MSR_ML_memorydata_Destroy;
   ML_memory_alloc((void**) &aggr_comm, sizeof(ML_Aggregate_Comm),"ACO");
   aggr_comm->comm = comm;
   aggr_comm->N_send_neighbors = new_N_send;
   aggr_comm->N_recv_neighbors = new_N_recv;
   aggr_comm->send_neighbors = new_send_neighbors;
   aggr_comm->recv_neighbors = new_recv_neighbors;
   aggr_comm->send_leng = new_send_leng;
   aggr_comm->recv_leng = new_recv_leng;
   aggr_comm->send_list = new_send_list;
   aggr_comm->local_nrows = Ncoarse * nullspace_dim;

   m = exp_Ncoarse - Ncoarse;
   ML_CommInfoOP_Generate( &((*Pmatrix)->getrow->pre_comm),
                           ML_Aggregate_ExchangeBdry, aggr_comm,
                           comm, Ncoarse*nullspace_dim, m*nullspace_dim);
   ML_Operator_Set_Getrow((*Pmatrix), Nrows, CSR_getrow);
   ML_Operator_Set_ApplyFunc((*Pmatrix), CSR_matvec);
   (*Pmatrix)->max_nz_per_row = nullspace_dim;

   /* ------------------------------------------------------------- */
   /* clean up                                                      */
   /* ------------------------------------------------------------- */

   ML_memory_free((void**) &neighbors);
   ML_memory_free((void**) &recv_leng);
   ML_memory_free((void**) &send_leng);
   ML_memory_free((void**) &send_list);
   ML_memory_free((void**) &recv_list);
   ML_memory_free((void**) &aggr_index);
   ML_memory_free((void**) &aggr_cnt_array);
   for (i = 0; i < aggr_count; i++) ML_free(rows_in_aggs[i]);
   ML_memory_free((void**)&rows_in_aggs);
   ML_memory_free((void**)&qr_tmp);
   ML_memory_free((void**)&tmp_vect);
   ML_memory_free((void**)&work);
   if ( new_N_send > 0 )
   {
      ML_memory_free((void**) &new_send_leng);
      ML_memory_free((void**) &new_send_list);
      ML_memory_free((void**) &new_send_neighbors);
   }
   if ( new_N_recv > 0 )
   {
      ML_memory_free((void**) &new_recv_leng);
      ML_memory_free((void**) &new_recv_neighbors);
   }
   ML_memory_free((void**) &aggr_comm);
   return Ncoarse*nullspace_dim;
}

/* ************************************************************************* */
/* This algorithm goes as follow :                                           */
/*                                                                           */
/* 1) At the beginning of the aggregation phase, each processor picks a      */
/*    point on the interprocessor boundary. Each processor then computes     */
/*    the graph distance between every other interprocessor boundary         */
/*    point and this first point. I was thinking that when computing         */
/*    these distances we only consider the graph corresponding to the        */
/*    interprocessor points.  We would need to handle the case where         */
/*    the interprocessor boundary is not completely connected. I've          */
/*    attached some matlab code that will hopefully make this clearer.       */
/* 2) When we create aggregates on the interprocessor boundary, we take      */
/*    the point with the smallest distance among all valid points and        */
/*    we use this point to build the next aggregate.                         */
/* 3) When finished with the interprocessor boundary, we do the same         */
/*    thing in the interior. That is, we pick a point in the interior        */
/*    and then compute a graph distance between it and every other           */
/*    interior point.  We might want to choose this first point to be        */
/*    close to the already computed aggregates on the interprocessor         */
/*    boundary ... if this is not hard.                                      */
/* 4) When creating interior aggregates, we take the point with the          */
/*    smallest distance among valid points.                                  */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_CoarsenCoupledCore(ML_Aggregate *ml_ag, ML_Comm *comm,
       int *mat_indx, int *aggr_count_out, int **aggr_index_out, 
       int N_neighbors, int *neighbors, int *recv_leng, int *send_leng, 
       int *send_list, int *recv_list, int **cnt_array, int *bc_array)
{
   int     i, j, k, inode, jnode, nbytes, length, Nrows;
   int     aggr_count, mypid, *order_array;
   int     *aggr_index, *order_array2;
   int     *aggr_stat, aggr_cnt_leng, *aggr_cnt_array;
   int     seed_node, *node_type, *node_dist, node_left, *dist_array;
   int     exp_Nrows, max_dist, *sendlist_proc;
   int     phaseAFlag;
   double  printflag;
   int     max_length=0;
   int     min_agg_size, max_neigh_selected;
   ML_Node       *node_head, *node_tail, *new_node;

   /* ============================================================= */
   /* get the machine information and matrix references             */
   /* ============================================================= */

   mypid = comm->ML_mypid;
   Nrows = mat_indx[0] - 1;
   exp_Nrows = Nrows;
   for ( i = 0; i < N_neighbors; i++ ) exp_Nrows += recv_leng[i];
   printflag = ml_ag->print_flag;

   /* ============================================================= */
   /* construct an array indicating boundary or interior nodes      */
   /* ============================================================= */

   nbytes = Nrows * sizeof( int );
   if ( nbytes > 0 ) ML_memory_alloc((void**) &node_type, nbytes, "Ag4");
   else              node_type = NULL;
   for ( i = 0; i < Nrows; i++ ) node_type[i] = 0; /* all interior */
   for ( i = 0; i < Nrows; i++ ) 
   {
      for ( j = mat_indx[i]; j < mat_indx[i+1]; j++ ) 
         if ( mat_indx[j] >= Nrows ) {node_type[i] = 1; break;}
   }

   /* ============================================================= */
   /* construct an array indicating isolated nodes                  */
   /* ============================================================= */

   for ( i = 0; i < Nrows; i++ ) 
   {
      length = mat_indx[i+1] - mat_indx[i];
      if ( length+1 > max_length ) max_length = length + 1;
   }

   /* ============================================================= */
   /* Pick a seed node in each processor and compute node distances */
   /* from seed point to every other node. The outer loop is here   */
   /* to deal with unconnected regions.                             */
   /* ============================================================= */

   ML_memory_alloc((void**) &node_dist, Nrows*sizeof(int), "Ag9");
   for ( i = 0; i < Nrows; i++ ) node_dist[i] = -1;
   node_left = Nrows - 1;

   while ( node_left > 0 )
   {
      /* ---------------------------------------------------------- */
      /* pick a seed node (if boundary and has not been visited)    */
      /* ---------------------------------------------------------- */

      seed_node = -1;
      for ( inode = 0; inode < Nrows; inode++ ) 
      { 
         if ( node_type[inode] == 1 && node_dist[inode] == -1 ) 
         {
            seed_node = inode; 
            break;
         } 
      } 
      if ( seed_node == -1 )
      {
         for ( inode = 0; inode < Nrows; inode++ ) 
            if ( node_dist[inode] == -1 ) { seed_node = inode; break; } 
      } 

      /* ---------------------------------------------------------- */
      /* initialize search queue                                    */
      /* ---------------------------------------------------------- */

      node_head = NULL;
      new_node = (ML_Node *) ML_allocate(sizeof(ML_Node));      
      new_node->node_id = seed_node;
      node_head = new_node;
      node_tail = new_node;
      new_node->next = NULL;
      node_dist[seed_node] = 0; 
   
      /* ---------------------------------------------------------- */
      /* process the subgraph                                       */
      /* ---------------------------------------------------------- */

      while ( node_head != NULL )
      {
         new_node = node_head;
         inode = new_node->node_id;
         node_head = new_node->next;
         ML_free(new_node);
         for ( j = mat_indx[inode]; j < mat_indx[inode+1]; j++ ) 
         {
            jnode = mat_indx[j];
            if ( jnode < Nrows && node_dist[jnode] == -1 ) 
            {
               node_dist[jnode] = node_dist[inode] + 1;
               new_node = (ML_Node *) ML_allocate(sizeof(ML_Node));
               new_node->node_id = jnode;
               new_node->next = NULL;
               if ( node_head == NULL ) {node_head = node_tail = new_node;}
               else { node_tail->next = new_node; node_tail = new_node; }
            }      
         }
      }

      node_left = 0;
      for ( j = 0; j < Nrows; j++ )
         if ( node_dist[j] == -1 ) node_left++;
   }

   /* ============================================================= */
   /* based on the node_dist information, build an reorder array    */
   /* ============================================================= */

   if ( ml_ag->ordering == 0 )
   {
      ML_memory_alloc((void**) &order_array, Nrows*sizeof(int), "Ag6");
      for (i = 0; i < Nrows; i++) order_array[i] = i;
   } 
   else if ( ml_ag->ordering == 1 )
   {
      ML_memory_alloc((void**) &order_array, Nrows*sizeof(int), "Ag6");
      for (i = 0; i < Nrows; i++) order_array[i] = i;
      ML_randomize(Nrows, order_array);
   } 
   else
   {
      max_dist = 0;
      for (i = 0; i < Nrows; i++) 
         if ( node_dist[i] > max_dist ) max_dist = node_dist[i];
      max_dist++;
      ML_memory_alloc((void**) &dist_array, max_dist*sizeof(int), "Ag0");
      for (i = 0; i < max_dist; i++) dist_array[i] = 0;
      for (i = 0; i < Nrows; i++) dist_array[node_dist[i]]++;
      k = dist_array[0];
      dist_array[0] = 0;
      for (i = 1; i < max_dist; i++) 
      {
         j = dist_array[i]; 
         dist_array[i] = k;
         k += j;
      }
      ML_memory_alloc((void**) &order_array, Nrows*sizeof(int), "Ag6");
      for (i = 0; i < Nrows; i++) 
         order_array[i] = dist_array[node_dist[i]]++;
      ML_memory_free( (void**) &dist_array );
      ML_memory_alloc((void**) &order_array2, Nrows*sizeof(int), "Ag6");
      for (i = 0; i < Nrows; i++) order_array2[i] = i; 
      ML_az_sort(order_array, Nrows, order_array2, NULL);
      ML_memory_free( (void**) &order_array );
      order_array = order_array2;
   }

   /* ============================================================= */
   /* prepare for aggregation                                       */
   /* ============================================================= */

   /* ------------------------------------------------------------- */
   /* sendlist_proc is used to find out, in the aggregation process,*/
   /* which processor holds neighbors of my local nodes             */
   /* ------------------------------------------------------------- */

   nbytes = (N_neighbors + 1 ) * sizeof(int);
   ML_memory_alloc((void**) &sendlist_proc, nbytes, "Aga");
   sendlist_proc[0] = 0;
   for ( i = 1; i <= N_neighbors; i++ )
      sendlist_proc[i] = sendlist_proc[i-1] + recv_leng[i-1];

   /* ------------------------------------------------------------- */
   /* set up bookkeeping arrays for aggregation                     */
   /* ------------------------------------------------------------- */

   aggr_count = 0;
   aggr_cnt_leng = Nrows / 5 + 2;
   nbytes = aggr_cnt_leng * sizeof(int);
   if (nbytes > 0) ML_memory_alloc((void**) &aggr_cnt_array,nbytes,"Agb");
   else            aggr_cnt_array = NULL;

   for ( i = 0; i < aggr_cnt_leng; i++ ) aggr_cnt_array[i] = 0;

   /* ------------------------------------------------------------- */
   /* Construct an initial status array to store whether the nodes  */
   /* have been aggregated. aggr_index stores the aggregate number  */
   /* where this node has been aggregated into.                     */
   /* ------------------------------------------------------------- */

   nbytes = exp_Nrows * sizeof( int );
   if ( nbytes > 0 )
   {
      ML_memory_alloc((void**) &aggr_index, nbytes, "Ag7");
      ML_memory_alloc((void**) &aggr_stat, nbytes, "Ag8");
   } else aggr_index = aggr_stat = NULL;

   for ( i = 0; i < Nrows; i++ )
   {
      if (bc_array[i] == 1) aggr_stat[i] = ML_AGGR_BDRY;
      else                  aggr_stat[i] = ML_AGGR_READY;
   }
   for ( i = 0; i < exp_Nrows; i++ ) aggr_index[i] = -1;
   for ( i = Nrows; i < exp_Nrows; i++ ) aggr_stat[i] = 0;

   /* ============================================================= */
   /* Phase 1 :                                                     */
   /*    This consists of two parts - aggregate border nodes first  */
   /*    followed by aggregating interior nodes.  This goes on      */
   /*    until all nodes are either selected or not selected.       */
   /* ============================================================= */

   min_agg_size = 2;
   max_neigh_selected = 0;
   phaseAFlag = 1;
   ML_Aggregate_Form_Aggregates('1', phaseAFlag, Nrows, mat_indx, 
        aggr_index, aggr_stat, node_type, NULL, order_array, &aggr_count, 
        &aggr_cnt_leng, &aggr_cnt_array, max_length, min_agg_size, 
        max_neigh_selected, N_neighbors, neighbors, send_leng, send_list, 
        recv_leng, recv_list, sendlist_proc, comm, printflag);

   /* ============================================================= */
   /* Phase 2 :                                                     */
   /*    This phase consists of forming aggregates where possible   */
   /*    if the aggregate size is large enough.  This phase is      */
   /*    different from phase 1 in that the seed node does not have */
   /*    to be next to no selected nodes.  This phase consists of   */ 
   /*    two parts - aggregate border nodes first followed by       */
   /*    the interior nodes.  This goes on until all nodes are      */
   /*    either selected or not selected.                           */
   /* ============================================================= */

   max_neigh_selected = 10000;
   min_agg_size = 10000;
   for ( i = 0; i < aggr_count; i++ )
      if ( aggr_cnt_array[i] < min_agg_size ) 
         min_agg_size = aggr_cnt_array[i];

   min_agg_size = - ML_gmax_int(-min_agg_size, comm) - 1;
   if ( min_agg_size <= 1 ) min_agg_size = 2;
   min_agg_size = 3;

   ML_Aggregate_Form_Aggregates('2',phaseAFlag,Nrows, mat_indx, 
        aggr_index, aggr_stat, node_type, NULL, order_array, &aggr_count, 
        &aggr_cnt_leng, &aggr_cnt_array, max_length, min_agg_size, 
        max_neigh_selected, N_neighbors, neighbors, send_leng, send_list, 
        recv_leng, recv_list, sendlist_proc, comm, printflag);

   /* ============================================================= */
   /* Phase 3 :                                                     */
   /*    for all nodes, see if it can be aggregated into one of the */
   /*    existing LOCAL aggregates.                                 */
   /* ============================================================= */

/*
   ML_Aggregate_PutInto_Aggregates('3', attach_scheme, mat_indx, 
        aggr_index, aggr_stat, &aggr_count, &aggr_cnt_array, 
        max_length, N_neighbors, neighbors, send_leng, send_list, 
        recv_leng, recv_list, comm, printflag);
*/

   /* ============================================================= */
   /* Phase 4 :                                                     */
   /*    for all remaining nodes, form new aggregates               */
   /* ------------------------------------------------------------- */

   min_agg_size = 1;
   max_neigh_selected = 10000;
   phaseAFlag = 0;
   ML_Aggregate_Form_Aggregates('4',phaseAFlag,Nrows,mat_indx, 
        aggr_index, aggr_stat, node_type, NULL, order_array, &aggr_count, 
        &aggr_cnt_leng, &aggr_cnt_array, max_length, min_agg_size, 
        max_neigh_selected, N_neighbors, neighbors, send_leng, send_list, 
        recv_leng, recv_list, sendlist_proc, comm, printflag);

   /* ============================================================= */
   /* final checking                                                */
   /* ============================================================= */

   for (i = 0; i < Nrows; i++)
   {
      if (aggr_index[i] >= aggr_count)
      {
         printf("WARNING (CC) : index out of range (%d,%d,%d,%d)\n",
                    mypid,i,aggr_index[i], aggr_count);
         break;
      }
      if (aggr_stat[i] != ML_AGGR_SELECTED && aggr_stat[i] != ML_AGGR_BDRY)
      {
         printf("%d : ERROR (CC) : node %d not aggregated (%d)\n", mypid,
                i, aggr_stat[i]);
         for ( j = mat_indx[i]; j < mat_indx[i+1]; j++ )
            printf("%d : neighbors = %d %d %d\n", mypid, mat_indx[j], 
                   aggr_stat[mat_indx[j]], aggr_index[mat_indx[j]]);
/*
         exit(1);
*/
         aggr_stat[i] = ML_AGGR_SELECTED;
         aggr_stat[i] = aggr_count++;
         aggr_cnt_array[aggr_count-1] = 1;
      }
   }

   /* ============================================================= */
   /* final clean up                                                */
   /* ============================================================= */

   ML_memory_free((void**) &sendlist_proc);
   ML_memory_free((void**) &node_dist);
   ML_memory_free((void**) &aggr_stat);
   ML_memory_free((void**) &node_type);
   ML_memory_free((void**) &order_array);
   (*cnt_array) = aggr_cnt_array;
   (*aggr_count_out) = aggr_count;
   (*aggr_index_out) = aggr_index;
   return 0;
}

/* ************************************************************************* */
/* ************************************************************************* */
/*          Support subroutines                                              */
/* ************************************************************************* */
/* ************************************************************************* */

/* ************************************************************************* */
/* ************************************************************************* */
/* Compress a matrix into a block matrix                                     */
/*    getrow_obj : for extracting communication information                  */
/*    mat_indx   : pruned matrix                                             */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Compress_Matrix(ML_GetrowFunc *getrow_obj, int *mat_indx, 
                int num_PDE_eqns, ML_Comm *comm, int **new_mat_indx, 
                int *N_neighbors, int **neighbors, int **recv_leng,
                int **send_leng, int **send_list, int **recv_list,
                int *bc_array)
{
   int    i, j, k, Nrows, nz_cnt, nbytes, *amal_mat_indx, LN_neighbors;
   int    *Lneighbors, *Lsend_leng, *Lrecv_leng, *Lsend_list, *Lrecv_list;
   int    *Aneighbors, *Asend_leng, *Arecv_leng, *Asend_list, *Arecv_list;
   int    AN_neighbors, total_send_leng, count, label, lcount;
   int    nblocks, row, amal_count, index;
   int    total_recv_leng, chk_bdry;
   char   *col_entered;
   double *dbuf;
   ML_CommInfoOP *getrow_comm;

   /* ------------------------------------------------------------- */
   /* retrieve matrix parameters                                    */
   /* ------------------------------------------------------------- */
 
   Nrows  = mat_indx[0] - 1;
   nz_cnt = mat_indx[Nrows];

   /* ------------------------------------------------------------- */
   /* retrieve communication information for the original A         */
   /* ------------------------------------------------------------- */

   AN_neighbors = getrow_obj->pre_comm->N_neighbors;
   nbytes = AN_neighbors * sizeof( int );
   if ( nbytes > 0 )
   {
      ML_memory_alloc((void**) &Aneighbors, nbytes, "AgA");
      ML_memory_alloc((void**) &Asend_leng, nbytes, "AgB");
      ML_memory_alloc((void**) &Arecv_leng, nbytes, "AgC");
   }
   else
   {
      Aneighbors = Arecv_leng = Asend_leng = NULL;
   }
   for ( i = 0; i < AN_neighbors; i++ )
   {
      Aneighbors[i] = getrow_obj->pre_comm->neighbors[i].ML_id;
      Arecv_leng[i] = getrow_obj->pre_comm->neighbors[i].N_rcv;
      Asend_leng[i] = getrow_obj->pre_comm->neighbors[i].N_send;
   }
   total_send_leng = 0;
   for ( i = 0; i < AN_neighbors; i++ ) total_send_leng += Asend_leng[i];
   nbytes = total_send_leng * sizeof( int );
   if ( nbytes > 0 ) ML_memory_alloc((void**) &Asend_list, nbytes, "AgD");
   else              Asend_list = NULL;
   count = 0;
   for ( i = 0; i < AN_neighbors; i++ )
   {
      for (j = 0; j < Asend_leng[i]; j++)
         Asend_list[count++] =
            getrow_obj->pre_comm->neighbors[i].send_list[j];
   }
   total_recv_leng = 0;
   for ( i = 0; i < AN_neighbors; i++ ) total_recv_leng += Arecv_leng[i];
   nbytes = total_recv_leng * sizeof( int );
   if ( nbytes > 0 ) ML_memory_alloc((void**) &Arecv_list, nbytes, "AgE");
   else              Arecv_list = NULL;
   count = 0;
   for ( i = 0; i < AN_neighbors; i++ )
   {
      for (j = 0; j < Arecv_leng[i]; j++)
         Arecv_list[count++] =
            getrow_obj->pre_comm->neighbors[i].rcv_list[j];
   }

   /* ------------------------------------------------------------- */
   /* if the system is scalar, the graph is the same                */
   /* ------------------------------------------------------------- */

   if ( num_PDE_eqns == 1 )
   {
      (*new_mat_indx) = mat_indx;
      (*N_neighbors) = AN_neighbors;
      (*neighbors) = Aneighbors;
      (*send_list) = Asend_list;
      (*recv_list) = Arecv_list;
      (*recv_leng) = Arecv_leng;
      (*send_leng) = Asend_leng;
      return 0;
   }

   /* ------------------------------------------------------------- */
   /* copy neighbor processor information - LN_neighbors, Lneighbors*/
   /* ------------------------------------------------------------- */

   LN_neighbors = AN_neighbors;
   nbytes = LN_neighbors * sizeof( int );
   if ( nbytes > 0 )
   {
      ML_memory_alloc((void**) &Lneighbors, nbytes, "AgA");
      ML_memory_alloc((void**) &Lsend_leng, nbytes, "AgB");
      ML_memory_alloc((void**) &Lrecv_leng, nbytes, "AgC");
   }
   else
   {
      Lneighbors = Lrecv_leng = Lsend_leng = NULL;
   }
   for ( i = 0; i < LN_neighbors; i++ ) Lneighbors[i] = Aneighbors[i];

   /* ------------------------------------------------------------- */
   /* compress send stuff (Lsend_leng, Lsend_list)                  */
   /* ------------------------------------------------------------- */

   count = 0;
   for ( i = 0; i < LN_neighbors; i++ )
   {
      for (j = 0; j < Asend_leng[i]; j++) Asend_list[j] /= num_PDE_eqns;
      lcount = 0;
      if ( Asend_leng[i] > 0 ) 
      {
         label = Asend_list[count];
         lcount++;
      }
      for ( j = 1; j < Asend_leng[i]; j++ )
      {
         index = Asend_list[count+j]; 
         if ( index != label )
         {
            lcount++; 
            label = index;
         }
      }
      Lsend_leng[i] = lcount;
      count += Asend_leng[i];
      for (j = 0; j < Asend_leng[i]; j++) Asend_list[j] *= num_PDE_eqns;
   }
   total_send_leng = 0;
   for ( i = 0; i < LN_neighbors; i++ ) total_send_leng += Lsend_leng[i];
   nbytes = total_send_leng * sizeof( int );
   if ( nbytes > 0 ) ML_memory_alloc((void**) &Lsend_list, nbytes, "AgD");
   else              Lsend_list = NULL;

   count = lcount = 0;
   for ( i = 0; i < LN_neighbors; i++ )
   {
      for (j = 0; j < Asend_leng[i]; j++) Asend_list[j] /= num_PDE_eqns;
      if ( Asend_leng[i] > 0 ) 
      {
         label = Asend_list[count];
         Lsend_list[lcount++] = label;
      }
      for ( j = 1; j < Asend_leng[i]; j++ )
      {
         index = Asend_list[count+j]; 
         if ( index != label ) 
         { 
            label = index;
            Lsend_list[lcount++] = label;
         }
      }
      count += Asend_leng[i];
      for (j = 0; j < Asend_leng[i]; j++) Asend_list[j] *= num_PDE_eqns;
   }
   
   /* ------------------------------------------------------------- */
   /* compress recv stuff (Lrecv_leng)                              */
   /* ------------------------------------------------------------- */

   total_recv_leng = 0;
   for ( i = 0; i < AN_neighbors; i++ ) total_recv_leng += Arecv_leng[i];
   nbytes = (Nrows + total_recv_leng) * sizeof( int );
   if ( nbytes > 0 ) dbuf = (double *) ML_allocate( nbytes );
   else              dbuf = NULL;
   for ( i = 0; i < Nrows; i++ ) dbuf[i] = i / num_PDE_eqns;
   getrow_comm = getrow_obj->pre_comm;
   if ( getrow_comm != NULL )
      ML_exchange_bdry(dbuf, getrow_comm, Nrows, comm, ML_OVERWRITE,NULL);

   count = Nrows;
   total_recv_leng = 0; 
   for ( i = 0; i < AN_neighbors; i++ )
   {
      lcount = 0;
      if ( Arecv_leng[i] > 0 ) 
      {
         label = (int) dbuf[count];
         lcount++;
      }
      for ( j = 1; j < Arecv_leng[i]; j++ )
      {
         index = (int) dbuf[count+j];
         if ( index != label )
         {
            lcount++; 
            label = index;
         }
      }
      Lrecv_leng[i] = lcount;
      count += Arecv_leng[i];
      total_recv_leng += lcount;
   }

   if ( Arecv_list != NULL )
   {
      count = 0;
      for ( i = 0; i < AN_neighbors; i++ )
      {
         for ( j = 0; j < Arecv_leng[i]; j++ )
         {
            if ( Arecv_list[count] != count )
            {
               printf("Arecv_list cannot be processed.\n");
               exit(-1);
            }
            count++;
         }
      }

      nbytes = total_recv_leng * sizeof( int );
      if ( nbytes > 0 ) ML_memory_alloc((void**) &Lrecv_list, nbytes, "AgE");
      else              Lrecv_list = NULL;

      count = lcount = 0;
      for ( i = 0; i < LN_neighbors; i++ )
      {
         if ( Arecv_leng[i] > 0 ) 
         {
            label = (int) dbuf[Arecv_list[count]];
            Lrecv_list[lcount] = lcount;
            dbuf[count] = (double) lcount;
            lcount++;
         }
         for ( j = 1; j < Arecv_leng[i]; j++ )
         {
            index = (int) dbuf[Arecv_list[count+j]]; 
            if ( index != label ) 
            { 
               label = index;
               /* JJH */
               /*Lrecv_list[lcount++] = lcount;*/
               Lrecv_list[lcount] = lcount;
               lcount++;
               /* --JJH */
            }
            dbuf[count+j] = (double) lcount;
         }
         count += Arecv_leng[i];
      }
   } else Lrecv_list = NULL;

   /* ------------------------------------------------------------- */
   /* allocate storage for block matrix                             */
   /* ------------------------------------------------------------- */

   nbytes = (nz_cnt + Nrows + 1) * sizeof( int ); /* probably excessive */
   ML_memory_alloc((void**) &amal_mat_indx, nbytes, "AgF");

   /* ------------------------------------------------------------- */
   /* start compressing                                             */
   /* ------------------------------------------------------------- */

   nblocks = Nrows / num_PDE_eqns;
   amal_count = nblocks + 1;
   amal_mat_indx[0] = amal_count;
   row = 0;
   col_entered = (char *) ML_allocate(sizeof(char)*(total_recv_leng+nblocks) );
   if (col_entered == NULL) 
   {
      printf("Not enough space in ML_aggregate\n");
      exit(1);
   }
   for ( i = 0; i < nblocks+total_recv_leng; i++) col_entered[i] = 'F';

   for ( i = 0; i < nblocks; i++) 
   {
      col_entered[i] = 'T';
      chk_bdry = 0;
      for ( j = 0; j < num_PDE_eqns; j++) 
      {
         if ( mat_indx[row+1] - mat_indx[row] == 0 ) chk_bdry++;
         if ( mat_indx[row+1] - mat_indx[row] == 1 && 
              mat_indx[mat_indx[row]] == row ) chk_bdry++;

         for ( k = mat_indx[row]; k < mat_indx[row+1]; k++) 
         {
            if ( mat_indx[k] < Nrows ) index = mat_indx[k] / num_PDE_eqns;
            else                       index = (int) dbuf[mat_indx[k]-Nrows];
            if (col_entered[index] == 'F') 
            {
               amal_mat_indx[ amal_count++] = index;
               col_entered[index] = 'T';
            }
         }
         row++;
      }
      if ( chk_bdry == num_PDE_eqns ) bc_array[i] = 1;
      else                            bc_array[i] = 0;
      if ( amal_count > nz_cnt )
         printf("WARNING in Compress_Matrix : not enough memory %d %d %d %d %d\n",
                amal_count, nz_cnt, comm->ML_mypid, i, nblocks);

      amal_mat_indx[i+1] = amal_count;
      col_entered[i] = 'F';
      for ( j = amal_mat_indx[i]; j < amal_mat_indx[i+1]; j++)
         col_entered[ amal_mat_indx[j]] = 'F';
   }
   if ( dbuf != NULL ) ML_free( dbuf );

   /* ------------------------------------------------------------- */
   /* shuffle pointers                                              */
   /* ------------------------------------------------------------- */

   (*new_mat_indx) = amal_mat_indx;
   (*N_neighbors) = LN_neighbors;
   (*neighbors) = Lneighbors;
   (*send_list) = Lsend_list;
   (*recv_list) = Lrecv_list;
   (*send_leng) = Lsend_leng;
   (*recv_leng) = Lrecv_leng;

   /* ------------------------------------------------------------- */
   /* clean up                                                      */
   /* ------------------------------------------------------------- */

   ML_free(col_entered);
   ML_free(Aneighbors);
   ML_free(Asend_leng);
   ML_free(Arecv_leng);
   ML_free(Asend_list);
   ML_free(Arecv_list);

   return 0;
}

/* ************************************************************************* */
/* Exchange data between processors given communication information          */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_ExchangeStatus(char *recvbuf,char *sendbuf,int N_neighbors,
              int *neighbors,int *recv_leng,int *send_leng,int *recv_list,
              int Nrows, int msgid, int datatype, ML_Comm *comm)
{
   int     i, nbytes, fromproc, length, typeleng, msgtype, offset;
   int     total_recv_leng, *int_array, *iarray;
   char    *char_array, *carray;
   double  *dble_array, *darray;
   USR_REQ *Request;

   switch ( datatype ) 
   {
      case ML_CHAR    : typeleng = sizeof(char);   break;
      case ML_INT     : typeleng = sizeof(int);    break;
      case ML_DOUBLE  : typeleng = sizeof(double); break;
      default :         typeleng = datatype;       break;
   }

   nbytes = N_neighbors * sizeof(USR_REQ);
   if ( nbytes > 0 ) Request = (USR_REQ *) ML_allocate(nbytes);
   else              Request = NULL;
   offset = 0;
   msgtype = msgid;
   for ( i = 0; i < N_neighbors; i++ ) 
   {
      fromproc = neighbors[i];
      length = recv_leng[i] * typeleng;
      if ( length > 0 )
         comm->USR_irecvbytes(&recvbuf[offset*typeleng],length,&fromproc,
#ifdef ML_CPP
                           &msgtype, comm->USR_comm, &Request[i] );
#else
                           &msgtype, comm->USR_comm, (void *) &Request[i] );
#endif
      offset += recv_leng[i];
   }
   offset = 0;
   msgtype = msgid;
   for ( i = 0; i < N_neighbors; i++ ) 
   {
      length = send_leng[i] * typeleng;
      if ( length > 0 )
         comm->USR_sendbytes((void*) &sendbuf[offset*typeleng], length,
                             neighbors[i], msgtype, comm->USR_comm );
      offset += send_leng[i];
   }
   offset = 0;
   for ( i = 0; i < N_neighbors; i++ ) 
   {
      fromproc = neighbors[i];
      length = recv_leng[i] * typeleng;
      msgtype = msgid;
      if ( length > 0 )
         comm->USR_cheapwaitbytes(&recvbuf[offset*typeleng], length, &fromproc,
#ifdef ML_CPP
                             &msgtype, comm->USR_comm, &Request[i] );
#else
                             &msgtype, comm->USR_comm, (void *) &Request[i] );
#endif
      offset += recv_leng[i];
   }
   if ( Request != NULL ) ML_free( Request );

   /* ----------------------------------------------------------------- */
   /* if a receive list is given, then permute the incoming data        */
   /* ----------------------------------------------------------------- */

   if ( recv_list != NULL )
   {
      total_recv_leng = 0;
      for ( i = 0; i < N_neighbors; i++ ) total_recv_leng += recv_leng[i];
      switch ( datatype ) 
      {
         case ML_CHAR    : nbytes = total_recv_leng * sizeof(char);
                           char_array = (char *) ML_allocate(nbytes);
                           carray = (char *) recvbuf;
                           for ( i = 0; i < total_recv_leng; i++ )
                              char_array[recv_list[i]-Nrows] = carray[i];
                           for ( i = 0; i < total_recv_leng; i++ )
                              carray[i] = char_array[i];
                           ML_free(char_array);
                           break;
         case ML_INT     : nbytes = total_recv_leng * sizeof(int);
                           int_array = (int *) ML_allocate(nbytes);
                           iarray = (int *) recvbuf;
                           for ( i = 0; i < total_recv_leng; i++ )
                              int_array[recv_list[i]-Nrows] = iarray[i];
                           for ( i = 0; i < total_recv_leng; i++ )
                              iarray[i] = int_array[i];
                           ML_free(int_array);
                           break;
         case ML_DOUBLE  : nbytes = total_recv_leng * sizeof(double);
                           dble_array = (double *) ML_allocate(nbytes);
                           darray = (double *) recvbuf;
                           for ( i = 0; i < total_recv_leng; i++ )
                              dble_array[recv_list[i]-Nrows] = darray[i];
                           for ( i = 0; i < total_recv_leng; i++ )
                              darray[i] = dble_array[i];
                           ML_free(dble_array);
                           break;
      }
   }
   return 0;
}

/* ************************************************************************* */
/* compose send and receive information in view of block                     */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_ComposeExpandedCommInfo(ML_GetrowFunc *getrow_obj, 
           int num_PDEs, ML_Comm *comm, 
           int *out_N_neighbors, int **out_neighbors, 
           int **out_recv_leng, int **out_send_leng, int **out_send_list,
           int **out_recv_list)
{
   int      i, j, k, m, N_neighbors, *send_leng, *recv_leng, *send_list;
   int      *neighbors, nbytes, total_send_leng, total_recv_leng;
   int      index, count, count2, label, nprocs, msgtype;
   int      new_N_neighbors, *new_send_leng, *new_send_list;
   int      *new_recv_list, *new_recv_leng, *new_neighbors;
   int      toproc, fromproc, *recv_list; 
   USR_REQ  *Request;

   /* ----------------------------------------------------------------- */
   /* get machine information                                           */
   /* ----------------------------------------------------------------- */

   nprocs = comm->ML_nprocs;
   if ( nprocs == 1 )
   {
      (*out_N_neighbors) = 0;
      (*out_neighbors) = NULL;
      (*out_send_leng) = NULL;
      (*out_recv_leng) = NULL;
      (*out_send_list) = NULL;
      return 0;
   }

   /* ----------------------------------------------------------------- */
   /* allocate storage for the communication information                */
   /* ----------------------------------------------------------------- */

   N_neighbors = getrow_obj->pre_comm->N_neighbors;
   nbytes = N_neighbors * sizeof( int );
   if ( nbytes > 0 ) 
   {
      ML_memory_alloc((void**) &neighbors,  nbytes, "CI1");
      ML_memory_alloc((void**) &send_leng,  nbytes, "CI2");
      ML_memory_alloc((void**) &recv_leng,  nbytes, "CI3");
   } 
   else neighbors = send_leng = recv_leng = NULL;

   for ( i = 0; i < N_neighbors; i++ ) 
   {
      neighbors[i] = getrow_obj->pre_comm->neighbors[i].ML_id;
      send_leng[i] = getrow_obj->pre_comm->neighbors[i].N_send;
      recv_leng[i] = getrow_obj->pre_comm->neighbors[i].N_rcv;
   }
   total_send_leng = 0;
   for ( i = 0; i < N_neighbors; i++ ) total_send_leng += send_leng[i];

   nbytes = total_send_leng * sizeof( int );
   if ( nbytes > 0 ) ML_memory_alloc((void**) &send_list,nbytes,"CI4");
   else              send_list = NULL;

   count = 0;
   for ( i = 0; i < N_neighbors; i++ )
   {
      for (j = 0; j < send_leng[i]; j++)
         send_list[count++] =
            getrow_obj->pre_comm->neighbors[i].send_list[j];
   }

   total_recv_leng = 0;
   for ( i = 0; i < N_neighbors; i++ ) total_recv_leng += recv_leng[i];
   nbytes = total_recv_leng * sizeof( int );
   if ( nbytes > 0 ) ML_memory_alloc((void**) &recv_list,nbytes,"CI5");
   else              recv_list = NULL;
   count = 0;
   for ( i = 0; i < N_neighbors; i++ )
   {
      for (j = 0; j < recv_leng[i]; j++)
         recv_list[count++] =
            getrow_obj->pre_comm->neighbors[i].rcv_list[j];
   }

   /* ----------------------------------------------------------------- */
   /* check to see if the block sizes are all 1's.                      */
   /* ----------------------------------------------------------------- */

   if ( num_PDEs == 1 )
   {
      (*out_N_neighbors) = N_neighbors;
      (*out_neighbors) = neighbors;
      (*out_send_leng) = send_leng;
      (*out_recv_leng) = recv_leng;
      (*out_send_list) = send_list;
      (*out_recv_list) = recv_list;
      return 0;
   }

   /* ----------------------------------------------------------------- */
   /* compose send information based on block information               */
   /* ----------------------------------------------------------------- */

   new_N_neighbors = N_neighbors;
   nbytes = new_N_neighbors * sizeof( int );
   if ( nbytes > 0 ) 
   {
      ML_memory_alloc((void**) &new_neighbors,  nbytes, "CI1");
      ML_memory_alloc((void**) &new_send_leng,  nbytes, "CI2");
   } 
   else new_neighbors = new_send_leng = NULL;

   count = 0;
   for ( i = 0; i < N_neighbors; i++ )
   {
      count2 = 0;
      if ( send_leng[i] > 0 )
      {
         index = send_list[count];
         label = index / num_PDEs;
         count2 += num_PDEs; 
      }
      for ( j = 1; j < send_leng[i]; j++ )
      {
         index = send_list[count+j];
         m     = index / num_PDEs;
         if ( m != label )
         {
            label = m;
            count2 += num_PDEs;
         }
      }
      count += send_leng[i];
      new_send_leng[i] = count2;
   }

   total_send_leng = 0;
   for (i = 0; i < new_N_neighbors; i++) total_send_leng += new_send_leng[i];
   nbytes = total_send_leng * sizeof( int );
   if ( nbytes > 0 ) ML_memory_alloc((void**) &new_send_list,nbytes,"CI6");

   count = count2 = 0;
   for ( i = 0; i < N_neighbors; i++ )
   {
      new_neighbors[i] = neighbors[i];
      if ( send_leng[i] > 0 )
      {
         index = send_list[count];
         label = index / num_PDEs;
         for ( k = 0; k < num_PDEs; k++ )
            new_send_list[count2++] = label * num_PDEs + k;
      }
      for ( j = 1; j < send_leng[i]; j++ )
      {
         index = send_list[count+j];
         m     = index / num_PDEs;
         if ( m != label )
         {
            label = m;
            for ( k = 0; k < num_PDEs; k++ )
               new_send_list[count2++] = label * num_PDEs + k;
         }
      }
      count += send_leng[i];
   }
      
   /* ----------------------------------------------------------------- */
   /* compute receive information                                       */
   /* ----------------------------------------------------------------- */

   nbytes = N_neighbors * sizeof(int);
   if ( nbytes > 0 ) ML_memory_alloc((void**)&new_recv_leng, nbytes, "RI1");
   else              new_recv_leng = NULL;
   nbytes = N_neighbors * sizeof(USR_REQ);
   if ( nbytes > 0 ) Request = (USR_REQ *) ML_allocate(nbytes);
   else              Request = NULL;
   msgtype = 97531;
   for ( i = 0; i < new_N_neighbors; i++ ) 
   {
      fromproc = new_neighbors[i];
      comm->USR_irecvbytes(&new_recv_leng[i],sizeof(int),&fromproc,
#ifdef ML_CPP
                &msgtype, comm->USR_comm, &Request[i] );
#else
                &msgtype, comm->USR_comm, (void *) &Request[i] );
#endif
   }
   for ( i = 0; i < new_N_neighbors; i++ ) 
   {
      toproc = new_neighbors[i];
      comm->USR_sendbytes((void*) &new_send_leng[i], sizeof(int), 
                   toproc, msgtype, comm->USR_comm);
   }
   for ( i = 0; i < new_N_neighbors; i++)
   {
      fromproc = new_neighbors[i];
      comm->USR_cheapwaitbytes((void*) &new_recv_leng[i], sizeof(int), 
#ifdef ML_CPP
               &fromproc,&msgtype,comm->USR_comm,&Request[i]);
#else
               &fromproc,&msgtype,comm->USR_comm,(void *)&Request[i]);
#endif
   }

   /* ----------------------------------------------------------------- */
   /* generate the new receive list                                     */
   /* ----------------------------------------------------------------- */

   total_recv_leng = 0;
   for (i = 0; i < N_neighbors; i++) total_recv_leng += recv_leng[i];
   nbytes = total_recv_leng * sizeof( int );
   if ( nbytes > 0 ) ML_memory_alloc((void**)&new_recv_list, nbytes, "RI3");
   else              new_recv_leng = NULL;
   for ( i = 0; i < total_recv_leng; i++ ) new_recv_list[i] = i;

   /* ----------------------------------------------------------------- */
   /* final clean up                                                    */
   /* ----------------------------------------------------------------- */

   ML_memory_free((void**) &send_leng);
   ML_memory_free((void**) &send_list);
   ML_memory_free((void**) &recv_leng);
   ML_memory_free((void**) &neighbors);

   (*out_N_neighbors) = new_N_neighbors;
   (*out_neighbors) = new_neighbors;
   (*out_send_leng) = new_send_leng;
   (*out_recv_leng) = new_recv_leng;
   (*out_send_list) = new_send_list;
   (*out_recv_list) = new_recv_list;
   return 0;
}

/* ************************************************************************* */
/* compose receive information from send information                         */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_ComposeRecvFromSend(int nprocs, int mypid, int new_N_send,
       int *new_send_leng, int *new_send_neighbors, int *N_rcv, 
       int **recv_leng, int **recv_neighbors, ML_Comm *comm)
{
   int     i, nbytes, *int_buf, *int_buf2, new_N_rcv, *new_recv_neighbors;
   int     msgtype, fromproc, *new_recv_leng;
   USR_REQ *Request;

   if ( nprocs > 1 )
   {
      /* -------------------------------------------------------------- */
      /* compute new_N_rcv                                              */
      /* -------------------------------------------------------------- */

      int_buf  = (int *) ML_allocate(nprocs * sizeof(int));
      int_buf2 = (int *) ML_allocate(nprocs * sizeof(int));
      for (i = 0; i < nprocs; i++) int_buf[i] = 0;
      for (i = 0; i < new_N_send; i++) int_buf[new_send_neighbors[i]] = 1;
      /*ML_Comm_GsumVecInt(comm, int_buf, int_buf2, nprocs);*/
      ML_gsum_vec_int(&int_buf, &int_buf2, nprocs, comm);
      ML_free( int_buf2 );
      new_N_rcv = int_buf[mypid];
      ML_free( int_buf );

      /* -------------------------------------------------------------- */
      /* get receive processor and the receive length                   */
      /* -------------------------------------------------------------- */

      nbytes = new_N_rcv * sizeof(int);
      if ( nbytes > 0 )
      {
         ML_memory_alloc((void**)&new_recv_leng, nbytes, "AGM");
         ML_memory_alloc((void**)&new_recv_neighbors, nbytes, "AGM");
      } else new_recv_leng = new_recv_neighbors = NULL;
      nbytes = new_N_rcv * sizeof(USR_REQ);
      if ( nbytes > 0 ) Request = (USR_REQ *) ML_allocate(nbytes);
      else              Request = NULL;
      msgtype = 97531;
      for ( i = 0; i < new_N_rcv; i++ ) 
      {
         fromproc = -1;
         comm->USR_irecvbytes(&new_recv_leng[i],sizeof(int),&fromproc,
#ifdef ML_CPP
                   &msgtype, comm->USR_comm, &Request[i] );
#else
                   &msgtype, comm->USR_comm, (void *)&Request[i] );
#endif
      }
      for ( i = 0; i < new_N_send; i++ ) 
      {
         comm->USR_sendbytes((void*) &new_send_leng[i], sizeof(int), 
                      new_send_neighbors[i], msgtype, comm->USR_comm);
      }
      for ( i = 0; i < new_N_rcv; i++)
      {
         fromproc = -1;
         comm->USR_cheapwaitbytes((void*) &new_recv_leng[i], sizeof(int), 
#ifdef ML_CPP
                  &fromproc,&msgtype,comm->USR_comm,&Request[i]);
#else
                  &fromproc,&msgtype,comm->USR_comm,(void *)&Request[i]);
#endif
         new_recv_neighbors[i] = fromproc;
      }
      ML_az_sort( new_recv_neighbors, new_N_rcv, new_recv_leng, NULL);
      if ( new_N_rcv > 0 ) ML_free( Request );

      (*recv_leng) = new_recv_leng;
      (*recv_neighbors) = new_recv_neighbors;
      (*N_rcv) = new_N_rcv;
   }
   else
   {
      (*recv_leng) = NULL;
      (*recv_neighbors) = NULL;
      (*N_rcv) = 0;
   } 
   return 0;
} 

/* ************************************************************************* */
/* form new aggregates                                                       */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Form_Aggregates(char phaseID, int phaseAFlag, int Nrows, 
        int *mat_indx, int *aggr_index, int *aggr_stat, 
        int *node_type, int *node_type2, int *order_array, int *aggr_count_in, 
        int *aggr_cnt_leng_in,
        int **aggr_cnt_array_in, int max_row_nnz, int min_agg_size, 
        int max_neigh_selected, int N_neighbors, int *neighbors, 
        int *send_leng, int *send_list, int *recv_leng, int *recv_list, 
        int *sendlist_proc, ML_Comm *comm, double printflag)
{
   int           i, j, k, index, inode, inode2, jnode, mdiff, count;
   int           mypid, msgtype, procnum, *com_buf, *com_buf2, nready;
   int           nbytes, total_send_leng, total_recv_leng, select_flag;
   int           aggr_count, *aggr_cnt_array, loop_flag, *itmp_array;
   int           nneigh_selected, node_type_scalar;
   int           nselected, total_nselected, total_Nrows, nwaiting;
   int           total_nwaiting, total_aggr_count, aggr_cnt_leng; 
   int           sqrtnprocs;
   ML_SuperNode  *supernode;

   /* ---------------------------------------------------------------- */
   /* get inputs and compute local parameters                          */
   /* ---------------------------------------------------------------- */

   mypid = comm->ML_mypid;
   sqrtnprocs = comm->ML_nprocs;
   sqrtnprocs = (int) sqrt((double) sqrtnprocs);
   aggr_count = (*aggr_count_in);
   aggr_cnt_leng = (*aggr_cnt_leng_in);
   aggr_cnt_array = (*aggr_cnt_array_in);
   nbytes = aggr_cnt_leng * sizeof(int);
   if ( nbytes > 0 ) ML_memory_alloc((void**)&aggr_cnt_array,nbytes,"AGF"); 
   for ( i = 0; i < aggr_count; i++ ) 
      aggr_cnt_array[i] = (*aggr_cnt_array_in)[i];
   if ( nbytes > 0 )
      ML_memory_free((void**)aggr_cnt_array_in);

   total_send_leng = 0;
   for ( i = 0; i < N_neighbors; i++ ) total_send_leng += send_leng[i];
   total_recv_leng = 0;
   for ( i = 0; i < N_neighbors; i++ ) total_recv_leng += send_leng[i];

   for ( inode = 0; inode < Nrows; inode++ )
      if (aggr_stat[inode] != ML_AGGR_SELECTED &&
          aggr_stat[inode] != ML_AGGR_BDRY) aggr_stat[inode] = ML_AGGR_READY;
   nbytes = total_send_leng * sizeof(int);
   if ( nbytes > 0 ) com_buf = (int *) ML_allocate(nbytes);
   else              com_buf = com_buf2 = NULL;
   for (i = 0; i < total_send_leng; i++) com_buf[i] = aggr_stat[send_list[i]];
   msgtype = 18445;
   ML_Aggregate_ExchangeStatus((char*)&aggr_stat[Nrows],(char*) com_buf,
                N_neighbors, neighbors, recv_leng, send_leng, recv_list,
                Nrows, msgtype, ML_INT, comm);

   supernode = (ML_SuperNode *) ML_allocate(sizeof(ML_SuperNode));
   supernode->list = (int*) ML_allocate(max_row_nnz*sizeof(int));
   supernode->maxlength = max_row_nnz;

   /* ---------------------------------------------------------------- */
   /* Phase A handles the border nodes                                 */
   /* ---------------------------------------------------------------- */

   if ( phaseAFlag )
   {
      /* ------------------------------------------------------------- */
      /* label all READY nodes that have to wait for its neighbors     */
      /* (If any of my neighbors reside in processors with processor   */
      /*  ID lower than my processor ID, then my node has to wait)     */
      /* ------------------------------------------------------------- */

      for ( inode = 0; inode < Nrows; inode++ )
      {
         if ( aggr_stat[inode] == ML_AGGR_READY && 
              node_type[inode] == ML_AGGR_BORDER) 
         {
            for (jnode=mat_indx[inode]; jnode<mat_indx[inode+1]; jnode++)
            {
               index = mat_indx[jnode];
               if (aggr_stat[index] >= 0 || aggr_stat[index] == ML_AGGR_READY)
               {
                  mdiff = index - Nrows;

                  /* ------------------------------------------------- */
                  /* search for the processor the node is coming from  */
                  /* ------------------------------------------------- */

                  for ( k = 0; k <= N_neighbors; k++ )
                     if ( mdiff < sendlist_proc[k] ) break;

                  /* ------------------------------------------------- */
                  /* if the processor number < mypid, tag it with the  */
                  /* neighbor processor with the smallest rank         */
                  /* ------------------------------------------------- */

                  if ( k != 0 && neighbors[k-1] < mypid &&
                       aggr_stat[inode] < 0 )
                     aggr_stat[inode] = neighbors[k-1];
                  else if ( k != 0 && neighbors[k-1] < mypid &&
                            neighbors[k-1] < aggr_stat[inode] )
                     aggr_stat[inode] = neighbors[k-1];
               }
            }
         }
      }

      /* ------------------------------------------------------------- */
      /* send my status information to remote processors               */
      /* ------------------------------------------------------------- */

      if ( total_send_leng > total_recv_leng )
         nbytes = total_send_leng * sizeof(int);
      else
         nbytes = total_recv_leng * sizeof(int);
      if ( nbytes > 0 ) 
      {
         com_buf  = (int *) ML_allocate(nbytes);
         com_buf2 = (int *) ML_allocate(nbytes);
      } else com_buf = com_buf2 = NULL;

      for ( i = 0; i < total_send_leng; i++ )
      {
         com_buf[i] = aggr_stat[send_list[i]];
      }
      msgtype = 13445;
      ML_Aggregate_ExchangeStatus((char*)&aggr_stat[Nrows],(char*) com_buf,
                   N_neighbors, neighbors, recv_leng, send_leng, recv_list,
                   Nrows, msgtype, ML_INT, comm);

      /* ------------------------------------------------------------- */
      /* begin aggregating the border nodes                            */
      /* ------------------------------------------------------------- */

      loop_flag = 1;

      while ( loop_flag != 0 )
      {
	/*
         old_aggr_count = aggr_count;
	*/

         for ( inode2 = 0; inode2 < Nrows; inode2++ )
         {
            inode = order_array[inode2];

            /* ------------------------------------------------------- */
            /* if it is a READY and boundary node, do the following    */
            /* ------------------------------------------------------- */

            if ( node_type[inode] == ML_AGGR_BORDER && 
                 aggr_stat[inode] == ML_AGGR_READY)
            {  
               /* ---------------------------------------------------- */
               /* first put the nodes in the supernode list            */
               /* ---------------------------------------------------- */

               supernode->length = 1;
               supernode->list[0] = inode;
               nneigh_selected = 0;
               select_flag = 1;
               nready = 0;

               /* ---------------------------------------------------- */
               /* examine all of its neighbors                         */
               /* ---------------------------------------------------- */

               for (jnode=mat_indx[inode]; jnode<mat_indx[inode+1]; jnode++)
               {
                  index = mat_indx[jnode];
                  if (aggr_stat[index] == ML_AGGR_SELECTED) nneigh_selected++; 
                  else if (aggr_stat[index] >= mypid && index >= Nrows)
                     supernode->list[supernode->length++] = index;
                  else if (aggr_stat[index] == ML_AGGR_READY && index < Nrows)
                  {
                     supernode->list[supernode->length++] = index;
                     nready++;
                  }
                  else if (aggr_stat[index] == ML_AGGR_NOTSEL && index < Nrows) 
                     supernode->list[supernode->length++] = index;
                  else if (aggr_stat[index] < mypid && index < Nrows) 
                     select_flag = 0;
                  else
                  {
                     select_flag = 0;
                     /*printf("%d : anomalous behavior %d %d.\n", mypid,
                             index, aggr_stat[index]);
                     */
                  }
               }

               /* ---------------------------------------------------- */
               /* if the number of neighbors having been selected      */ 
               /* exceeds the threshold, or the size of the aggregate  */
               /* is too small, then my node is not to be used as seed */
               /* node for aggregation.  If I do not have neighbors    */
               /* that can aggregate me, then do the following :       */
               /* if all my off-processor neighbors are all in         */
               /* processor with lower rank, set my status to be       */
               /* NOT_SEL, otherwise, set my status to be the          */
               /* processor closest to my processor rank.              */
               /* ---------------------------------------------------- */

               if ( nneigh_selected > max_neigh_selected )
               {
                  select_flag = 0;
                  if ( nready == 0 )
                  {
                     procnum = mypid;
                     for (j = mat_indx[inode]; j< mat_indx[inode+1];j++)
                     {
                        if (aggr_stat[index] >= 0 || 
                            aggr_stat[index] == ML_AGGR_READY)
                        {
                           index = mat_indx[j] - Nrows;
                           if ( index >= 0)
                           {
                              count = 0;
                              for (k = 0; k < N_neighbors; k++ )
                              {
                                 if ( index < (count+recv_leng[k]) ) break;
                                 count += recv_leng[k];
                              }
                              if (procnum == mypid && neighbors[k] > procnum) 
                                 procnum = neighbors[k];
                              else if (procnum != mypid && 
                                       neighbors[k] > mypid && 
                                       neighbors[k] < procnum) 
                                 procnum = neighbors[k];
                           }
                        }
                     }
                     if ( procnum == mypid ) 
                        aggr_stat[inode] = ML_AGGR_NOTSEL;
                     else 
                        aggr_stat[inode] = procnum;
                  }
               }
               else if ( supernode->length < min_agg_size || select_flag == 0)
               {
                  select_flag = 0;
                  aggr_stat[inode] = ML_AGGR_NOTSEL;
               }

               /* ---------------------------------------------------- */
               /* aggregation successful                               */
               /* ---------------------------------------------------- */

               if ( select_flag == 1 )
               {
                  /* ------------------------------------------------- */
                  /* label the aggregate members SELECTED, and assign  */
                  /* an aggregate number for the node (aggr_index)     */
                  /* ------------------------------------------------- */

#ifdef ML_AGGR_OUTPUT
printf("%d : Phase %cA - aggr %4d has ", mypid, phaseID, aggr_count);
#endif
                  for ( j = 0; j < supernode->length; j++ )
                  {
                     jnode = supernode->list[j];
                     aggr_stat[jnode] = ML_AGGR_SELECTED;
                     aggr_index[jnode] = aggr_count;
#ifdef ML_AGGR_OUTPUT
printf("%d ", jnode);
#endif
                     if ( jnode < Nrows )
                     {
                        for (k=mat_indx[jnode]; k<mat_indx[jnode+1]; k++)
                           if (aggr_stat[mat_indx[k]] == ML_AGGR_READY)
                              aggr_stat[mat_indx[k]] = ML_AGGR_NOTSEL;
                     }
                  }
#ifdef ML_AGGR_OUTPUT
printf("\n");
#endif
                  /* ------------------------------------------------- */
                  /* stretch aggr_cnt_array, if needed                 */
                  /* ------------------------------------------------- */

                  aggr_cnt_array[aggr_count++] = supernode->length;
                  if ( aggr_count >= aggr_cnt_leng )
                  {
                     itmp_array = aggr_cnt_array;
                     aggr_cnt_leng = aggr_cnt_leng * 6 / 5 + 1;
                     nbytes = aggr_cnt_leng * sizeof( int );
                     ML_memory_alloc((void**)&aggr_cnt_array,nbytes,"Agc");
                     for ( k = 0; k < aggr_count; k++ )
                        aggr_cnt_array[k] = itmp_array[k];
                     ML_memory_free((void**)&itmp_array);
                  }
               }
            }
         }

         /* ---------------------------------------------------------- */
         /* communicate remote node info back to remote processors     */
         /* (tell remote processor that some of their nodes have been  */
         /*  aggregated by this processor)                             */
         /* ---------------------------------------------------------- */

         msgtype = 33945 + loop_flag;
         ML_Aggregate_ExchangeStatus((char*)com_buf2,(char*)&aggr_stat[Nrows],
            N_neighbors,neighbors,send_leng,recv_leng,NULL,Nrows,msgtype,
            ML_INT,comm);

         /* ---------------------------------------------------------- */
         /* after my processor obtains information from other          */
         /* processors about my nodes being selected, update my local  */
         /* node status array aggr_stat (mark ML_AGGR_SELECTED for the */
         /* local nodes that have been aggregated by remote processors)*/
         /* ---------------------------------------------------------- */

         count = 0;
         for ( i = 0; i < N_neighbors; i++ )
         {
            for ( j = 0; j < send_leng[i]; j++ )
            {
               inode = send_list[count];
               if ( com_buf2[count] == ML_AGGR_SELECTED &&
                    aggr_stat[inode] != ML_AGGR_SELECTED )
               {
                  aggr_stat[inode]  = ML_AGGR_SELECTED;
                  aggr_index[inode] = - 100 - neighbors[i];
               }
               count++;
            }
         }

         /* ---------------------------------------------------------- */
         /* now my aggr_stat contains latest information about the     */
         /* status of the nodes I own.  Next, send this updated info   */
         /* to other processors                                        */
         /* ---------------------------------------------------------- */

         for ( i = 0; i < total_send_leng; i++ )
         {
            com_buf[i] = aggr_stat[send_list[i]];
         }
         msgtype = 13945 + loop_flag;
         ML_Aggregate_ExchangeStatus((char*)&aggr_stat[Nrows],(char*)com_buf,
            N_neighbors,neighbors,recv_leng,send_leng,recv_list,Nrows,
            msgtype,ML_INT,comm);

         /* ---------------------------------------------------------- */
         /* update my waiting nodes' status                            */
         /* ---------------------------------------------------------- */

         for ( inode = 0; inode < Nrows; inode++ )
         {
            if ( aggr_stat[inode] >= 0 )
            {
               procnum = mypid;
               for (jnode=mat_indx[inode];jnode<mat_indx[inode+1];jnode++)
               {
                  index = mat_indx[jnode];
                  mdiff = index - Nrows;
                  if ( mdiff >= 0 )
                  {
                     for ( k = 0; k <= N_neighbors; k++ )
                        if ( mdiff < sendlist_proc[k] ) break;
   
                     if ( aggr_stat[index] == ML_AGGR_READY &&
                          neighbors[k-1] < procnum )
                        procnum = neighbors[k-1];
                     else if (aggr_stat[index] >= 0 &&
                              aggr_stat[index] < mypid)
                     {
                        if ( neighbors[k-1] < procnum )
                           procnum = neighbors[k-1];
                     }
                  }
               }
               if ( procnum == mypid ) aggr_stat[inode] = ML_AGGR_READY;
               else                    aggr_stat[inode] = procnum;
            }
         }

         /* ---------------------------------------------------------- */
         /* now my aggr_stat contains latest information about the     */
         /* status of the nodes I own.  Next, send this updated info   */
         /* to other processors                                        */
         /* ---------------------------------------------------------- */

         for ( i = 0; i < total_send_leng; i++ )
         {
            com_buf[i] = aggr_stat[send_list[i]];
         }
         msgtype = 13965 + loop_flag;
         ML_Aggregate_ExchangeStatus((char*)&aggr_stat[Nrows],(char*)com_buf,
            N_neighbors,neighbors,recv_leng,send_leng,recv_list,Nrows,
            msgtype, ML_INT,comm);

         /* ---------------------------------------------------------- */
         /* output information about aggregation progress              */
         /* ---------------------------------------------------------- */

         nselected = 0;
         for (i = 0; i < Nrows; i++)
            if (aggr_stat[i] == ML_AGGR_SELECTED) nselected++;
         total_nselected = ML_Comm_GsumInt( comm, nselected);
         total_Nrows = ML_Comm_GsumInt( comm, Nrows);
         total_aggr_count = ML_Comm_GsumInt( comm, aggr_count );
         nwaiting = 0;
         for (i = 0; i < Nrows; i++) if (aggr_stat[i] >= 0) nwaiting++;
         total_nwaiting = ML_Comm_GsumInt( comm, nwaiting );
         if ( mypid == 0 && printflag < ML_Get_PrintLevel() )
         {
            printf("Aggregation(CC) : Phase %cA - Iteration        = %d\n",
                   phaseID, loop_flag);
            printf("Aggregation(CC) : Phase %cA - nodes aggregated = %d(%d)\n",
                   phaseID, total_nselected, total_Nrows);
            printf("Aggregation(CC) : Phase %cA - nodes waiting    = %d\n",
                   phaseID, total_nwaiting);
            printf("Aggregation(CC) : Phase %cA - total aggregates = %d\n",
                   phaseID, total_aggr_count);
         }

         /* ---------------------------------------------------------- */
         /* check to see if further loop is needed                     */
         /* ---------------------------------------------------------- */

/*
total_nwaiting = 0;
*/
         if ( total_nwaiting == 0 ) loop_flag = 0;
         else                       loop_flag++;
         if ( loop_flag > sqrtnprocs*sqrtnprocs ) loop_flag = 0; 
/*
         else
         {
            k = aggr_count - old_aggr_count;
            m = ML_Comm_GsumInt( comm, k);
            if ( m == 0 ) loop_flag = 0;
            else          loop_flag++;
         }
*/
      }
      if ( com_buf  != NULL ) ML_free( com_buf ); 
      if ( com_buf2 != NULL ) ML_free( com_buf2 ); 
   }

   /* ---------------------------------------------------------------- */
   /* turn the waiting nodes into READY nodes                          */
   /* ---------------------------------------------------------------- */

   for (i = 0; i < Nrows; i++)
      if (aggr_stat[i] >= 0) aggr_stat[i] = ML_AGGR_READY;

   /* ================================================================ */
   /* Phase B :                                                        */
   /*    aggregate interior nodes                                      */
   /* ---------------------------------------------------------------- */

   for ( inode2 = 0; inode2 < Nrows; inode2++ )
   {
      inode = order_array[inode2];

      if ( node_type2 != NULL ) node_type_scalar = node_type2[inode];
      else                      node_type_scalar = ML_AGGR_INTERIOR;

      /* ------------------------------------------------------------- */
      /* choose the ready nodes only                                   */
      /* ------------------------------------------------------------- */

      if ( aggr_stat[inode] == ML_AGGR_READY && 
           node_type_scalar == ML_AGGR_INTERIOR)
      {
         supernode->length = 1;
         supernode->list[0] = inode;
         nneigh_selected = 0;

         /* ---------------------------------------------------------- */
         /* examine all of its neighbors                               */
         /* ---------------------------------------------------------- */

         for (jnode=mat_indx[inode]; jnode<mat_indx[inode+1]; jnode++)
         {
            index = mat_indx[jnode];
            if ( index < Nrows )
            {
               if ( aggr_stat[index] == ML_AGGR_READY ||
                    aggr_stat[index] == ML_AGGR_NOTSEL )
                  supernode->list[supernode->length++] = index;
               else if ( aggr_stat[index] == ML_AGGR_SELECTED ) 
                  nneigh_selected++;
            }
         }

         /* ---------------------------------------------------------- */
         /* if critical mass is there, aggregation is successful       */
         /* ---------------------------------------------------------- */

         if ( nneigh_selected <= max_neigh_selected &&
              supernode->length >= min_agg_size )
         {
#ifdef ML_AGGR_OUTPUT
printf("%d : Phase %cB - aggr %4d has ", mypid, phaseID, aggr_count);
#endif
            for ( j = 0; j < supernode->length; j++ )
            {
               jnode = supernode->list[j];
#ifdef ML_AGGR_OUTPUT
printf("%d ", jnode);
#endif
               aggr_stat[jnode] = ML_AGGR_SELECTED;
               aggr_index[jnode] = aggr_count;
            }
#ifdef ML_AGGR_OUTPUT
printf("\n");
#endif
            aggr_cnt_array[aggr_count++] = supernode->length;
            if ( aggr_count >= aggr_cnt_leng )
            {
               itmp_array = aggr_cnt_array;
               aggr_cnt_leng = aggr_cnt_leng * 6 / 5 + 1;
               nbytes = aggr_cnt_leng * sizeof( int );
               ML_memory_alloc((void**)&aggr_cnt_array,nbytes,"Agf");
               for ( k = 0; k < aggr_count; k++ )
                  aggr_cnt_array[k] = itmp_array[k];
               ML_memory_free((void**)&itmp_array);
            }
         }
      }
   }
   ML_free( supernode->list );
   ML_free( supernode );

   /* ---------------------------------------------------------------- */
   /* send my aggr_stat information to other processors                */
   /* ---------------------------------------------------------------- */

   if ( total_send_leng > total_recv_leng )
      nbytes = total_send_leng * sizeof(int);
   else
      nbytes = total_recv_leng * sizeof(int);
   if ( nbytes > 0 ) com_buf = (int *) ML_allocate(nbytes);
   else              com_buf = NULL;

   for ( i = 0; i < total_send_leng; i++ )
   {
      com_buf[i] = aggr_stat[send_list[i]];
   }
   msgtype = 23965;
   ML_Aggregate_ExchangeStatus((char*)&aggr_stat[Nrows],(char*) com_buf,
      N_neighbors, neighbors,recv_leng,send_leng,recv_list,Nrows,msgtype,
      ML_INT, comm);
   if ( com_buf != NULL ) ML_free( com_buf );

   /* ---------------------------------------------------------------- */
   /* output information about aggregation progress                    */
   /* ---------------------------------------------------------------- */

   nselected = 0;
   for (i = 0; i < Nrows; i++)
      if (aggr_stat[i] == ML_AGGR_SELECTED) nselected++;
   total_nselected = ML_Comm_GsumInt( comm, nselected);
   total_Nrows = ML_Comm_GsumInt( comm, Nrows);
   total_aggr_count = ML_Comm_GsumInt( comm, aggr_count );
   if ( mypid == 0 && printflag < ML_Get_PrintLevel())
   {
      printf("Aggregation(CC) : Phase %cB - nodes aggregated = %d(%d)\n",
             phaseID, total_nselected, total_Nrows);
      printf("Aggregation(CC) : Phase %cB - total aggregates = %d\n",
             phaseID, total_aggr_count);
   }

   /* ---------------------------------------------------------------- */
   /* return the output parameters                                     */
   /* ---------------------------------------------------------------- */

   (*aggr_count_in)     = aggr_count;
   (*aggr_cnt_leng_in)  = aggr_cnt_leng;
   (*aggr_cnt_array_in) = aggr_cnt_array;

#ifdef ML_MPI
   MPI_Barrier(MPI_COMM_WORLD);
#endif

   return 0;
}

/* ************************************************************************* */
/* put nodes into existing local aggregates                                  */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_PutInto_Aggregates(char phaseID, int attach_scheme, 
        int *mat_indx, int *aggr_index, int *aggr_stat, int *aggr_count_in, 
        int **aggr_cnt_array_in, int N_neighbors, 
        int *neighbors, int *send_leng, int *send_list, int *recv_leng, 
        int *recv_list, ML_Comm *comm, double printflag )
{
   int          i, k, m, inode, jnode, index, mincount, select_flag;
   int          length, *int_array, *int_array2, *com_buf;
   int          total_send_leng, msgtype, mypid, aggr_count, *aggr_cnt_array;
   int          nselected, total_nselected, total_Nrows, total_aggr_count; 
   int          Nrows, nbytes;

   /* ----------------------------------------------------------------- */
   /* get inputs and compute local parameters                           */
   /* ----------------------------------------------------------------- */

   if ( mat_indx != NULL ) Nrows = mat_indx[0] - 1;
   else                    return 0;

   for (i = 0; i < Nrows; i++)
      if (aggr_stat[i] == ML_AGGR_NOTSEL ||
          aggr_stat[i] >= 0 ) aggr_stat[i] = ML_AGGR_READY;

   mypid = comm->ML_mypid;
   aggr_count = (*aggr_count_in);
   aggr_cnt_array = (*aggr_cnt_array_in);

   /* ----------------------------------------------------------------- */
   /* loop over all local nodes                                         */
   /* ----------------------------------------------------------------- */

   for ( inode = 0; inode < Nrows; inode++ )
   {
      /* -------------------------------------------------------------- */
      /* if the node in question is either READY or NOTSEL              */
      /* -------------------------------------------------------------- */

      if ( aggr_stat[inode] == ML_AGGR_READY || 
           aggr_stat[inode] == ML_AGGR_NOTSEL ) 
      { 
         select_flag = 0;

         /* ----------------------------------------------------------- */
         /* search for a neighboring aggregate that has the fewest nodes*/
         /* ----------------------------------------------------------- */

         if ( attach_scheme == ML_AGGR_MINRANK )
         {
            mincount = 100000;
            for (jnode=mat_indx[inode]; jnode<mat_indx[inode]; jnode++)
            {
               k = mat_indx[jnode];
               if ( k < Nrows )
               {
                  if (aggr_stat[k] == ML_AGGR_SELECTED &&
                      aggr_index[k] >= 0) /* locally aggregated */
                  {
                     select_flag = 1;
                     m = aggr_index[k];
                     if ( aggr_cnt_array[m] < mincount )
                     {
                        mincount = aggr_cnt_array[m];
                        index = k;
                     }
                  }
               }
            }
         }

         /* ----------------------------------------------------------- */
         /* search for a neighboring aggregate that has the most        */
         /* connection to my node                                       */
         /* ----------------------------------------------------------- */

         else if ( attach_scheme == ML_AGGR_MAXLINK )
         {
	   /*
            maxcount = 0;
	   */
            length = mat_indx[inode+1] - mat_indx[inode];
            if (length>0) int_array  = (int *) ML_allocate(length * sizeof(int));
            if (length>0) int_array2 = (int *) ML_allocate(length * sizeof(int));
            for ( i = 0; i < length; i++ ) int_array2[i] = i;
            length = 0;

            for (jnode=mat_indx[inode]; jnode<mat_indx[inode+1]; jnode++)
            {
               k = mat_indx[jnode];
               if (aggr_index[k] >= 0) 
               {
                  int_array2[length]  = k;
                  int_array[length++] = aggr_index[k];
               }
            }
            if ( length > 0 ) 
            {
               if (length > 1) ML_az_sort(int_array,length,int_array2,NULL);
               index = aggr_index[int_array2[length-1]];
               select_flag = 1;
            }
            length = mat_indx[inode+1] - mat_indx[inode];
            if ( length > 0 ) ML_free( int_array );
            if ( length > 0 ) ML_free( int_array2 );
         }

         /* ----------------------------------------------------------- */
         /* if search is successful, put thid node in the aggregate     */
         /* ----------------------------------------------------------- */

         if ( select_flag == 1 )
         {
            aggr_cnt_array[index]++;
            aggr_index[inode] = index;
            aggr_stat[inode] = ML_AGGR_SELECTED2;
         }
      }
   }

   /* ----------------------------------------------------------------- */
   /* restore the selected status (modified above to prevent chain)     */
   /* ----------------------------------------------------------------- */

   for ( inode = 0; inode < Nrows; inode++ )
      if ( aggr_stat[inode] == ML_AGGR_SELECTED2 )
         aggr_stat[inode] = ML_AGGR_SELECTED;

   /* ----------------------------------------------------------------- */
   /* communicate the information                                       */
   /* ----------------------------------------------------------------- */

   total_send_leng = 0;
   for ( i = 0; i < N_neighbors; i++ ) total_send_leng += send_leng[i];
   nbytes = total_send_leng * sizeof(int);
   if ( nbytes > 0 ) com_buf  = (int *) ML_allocate( nbytes );
   else              com_buf = NULL;
   for ( i = 0; i < total_send_leng; i++ ) 
   {
      com_buf[i] = aggr_stat[send_list[i]];
   }
   msgtype = 48934;
   ML_Aggregate_ExchangeStatus((char*)&aggr_stat[Nrows],(char*) com_buf,
        N_neighbors,neighbors,recv_leng,send_leng,recv_list,Nrows,msgtype,
        ML_INT,comm);
   if ( com_buf != NULL ) ML_free( com_buf );

   /* ----------------------------------------------------------------- */
   /* output information about aggregation progress                     */
   /* ----------------------------------------------------------------- */

   nselected = 0;
   for (i = 0; i < Nrows; i++)
      if (aggr_stat[i] == ML_AGGR_SELECTED) nselected++;
   total_nselected = ML_Comm_GsumInt( comm, nselected);
   total_Nrows = ML_Comm_GsumInt( comm, Nrows);
   total_aggr_count = ML_Comm_GsumInt( comm, aggr_count );
   if ( mypid == 0 && printflag < ML_Get_PrintLevel())
   {
      printf("Aggregation(CC) : Phase %c  - nodes aggregated = %d(%d)\n",
             phaseID, total_nselected, total_Nrows);
      printf("Aggregation(CC) : Phase %c  - total aggregates = %d\n",
             phaseID, total_aggr_count);
   }

   /* ---------------------------------------------------------------- */
   /* return the output parameters                                     */
   /* ---------------------------------------------------------------- */

   (*aggr_count_in)     = aggr_count;
   (*aggr_cnt_array_in) = aggr_cnt_array;

#ifdef ML_MPI
   MPI_Barrier(MPI_COMM_WORLD);
#endif

   return 0;
}

/* ************************************************************************* */
/* put nodes into existing local aggregates                                  */
/* ------------------------------------------------------------------------- */

int ML_Graph_CreateFromMatrix(ML_Aggregate *ml_ag, ML_Operator *Amatrix,
                              int **mat_indx_out, ML_Comm *comm, 
                              double epsilon, int *exp_Nrows_out, 
                              int **bc_array_out)
{
   int           i, j, k, m, jnode, length, count, nz_cnt, *mat_indx;
   int           Nrows, exp_Nrows, nbytes, maxnnz_per_row=500, *col_ind;
   int           (*getrowfunc)(ML_Operator *,int,int*,int,int*,double*,int*);
   void          *getrowdata;
   int           *bc_array;
   double        *diagonal, *dble_buf, *col_val, dcompare1, dcompare2;
   ML_GetrowFunc *getrow_obj;
   ML_CommInfoOP *getrow_comm;

   /* ============================================================= */
   /* fetch the getrow function for the incoming matrix             */
   /* ============================================================= */

   getrow_obj = Amatrix->getrow;
   getrowfunc = getrow_obj->internal;
   getrowdata = (void *) Amatrix;
   if ( getrowfunc == NULL )
   {
      printf("ML_Graph_CreateFromMatrix ERROR : null getrowfunc.\n");
      exit(-1);
   }

   /* ============================================================= */
   /* allocate initial temporary storage space for getrow           */
   /* also allocate space for storing the diagonal                  */
   /* ============================================================= */

   nbytes  = maxnnz_per_row * sizeof( int );
   col_ind = (int *) ML_allocate( nbytes );
   nbytes  = maxnnz_per_row * sizeof( double );
   col_val = (double *) ML_allocate( nbytes );
   Nrows   = Amatrix->outvec_leng;
   if ( Nrows > 0 ) diagonal = (double *) ML_allocate(Nrows * sizeof(double));
   else             diagonal = NULL;

   /* ============================================================= */
   /* fill in the diagonal array, also find out about the size of   */
   /* the incoming matrix (for allocation purpose)                  */
   /* ============================================================= */

   exp_Nrows = Nrows - 1;
   count = 0;
   for ( i = 0; i < Nrows; i++ )
   {
      diagonal[i]     = 0.0;
      while (getrowfunc(getrowdata,1,&i,maxnnz_per_row,col_ind,
                        col_val, &m) == 0 )
      {
         ML_free( col_ind );
         ML_free( col_val );
         maxnnz_per_row = maxnnz_per_row * 2 + 1;
         col_ind = (int *)    ML_allocate(maxnnz_per_row * sizeof(int));
         col_val = (double *) ML_allocate(maxnnz_per_row * sizeof(double));
      }
      for ( j = 0; j < m; j++ )
      {
         if ( col_ind[j] > exp_Nrows ) exp_Nrows = col_ind[j];
         if ( col_ind[j] == i )        diagonal[i] = col_val[j];
      }
      count += m;
      if ( diagonal[i] == 0.0 )
      {
         printf("%d : ML_Graph_CreateFromMatrix WARNING - diag %d is 0.\n",
                comm->ML_mypid,i);
         count++;
      }
   }
   exp_Nrows++;

   /* ============================================================= */
   /* build the diagonals of the expanded rows if epsilon != 0      */
   /* (diagonal elements are needed for pruning weak edges)         */
   /* ============================================================= */

   if ( epsilon == 0.0 && diagonal != NULL )
   {
      ML_free( diagonal );
      diagonal = NULL;
   }
   if ( epsilon != 0.0 && exp_Nrows > 0 )
   {
      dble_buf = diagonal;
      nbytes = exp_Nrows * sizeof(double);
      if ( nbytes > 0 ) diagonal = (double *) ML_allocate(nbytes);
      else              diagonal = NULL;
      for ( i = 0; i < Nrows; i++ ) diagonal[i] = dble_buf[i];
      for ( i = Nrows; i < exp_Nrows; i++ ) diagonal[i] = 0.0;
      if ( dble_buf != NULL ) ML_free(dble_buf);
      getrow_comm = getrow_obj->pre_comm;
      if ( getrow_comm != NULL )
         ML_exchange_bdry(diagonal,getrow_comm,Nrows,comm,ML_OVERWRITE,NULL);
   }

   /* ============================================================= */
   /* allocate temporary storage space for getrow                   */
   /* ============================================================= */

   nbytes = Nrows * sizeof( int );
   ML_memory_alloc((void**) &bc_array, nbytes, "ABC");
   nbytes = (count + Nrows + 1) * sizeof( int );
   ML_memory_alloc((void**) &mat_indx, nbytes, "ACG");
   k = ML_Comm_GsumInt( comm, Nrows);
   m = ML_Comm_GsumInt( comm, count);

   k = ML_Comm_GsumInt( comm, Nrows);
   m = ML_Comm_GsumInt( comm, count);

   if ( comm->ML_mypid == 0 && ml_ag->print_flag < ML_Get_PrintLevel())
      printf("Aggregation(CVB) : Total nnz = %d (Nrows=%d)\n",m,k);

   if ( ml_ag->operator_complexity == 0.0 )
   {
      ml_ag->fine_complexity = 1.0 * m;
      ml_ag->operator_complexity = 1.0 * m;
   }
   else
   {
      ml_ag->operator_complexity += 1.0 * m;
   }

   /* ============================================================= */
   /* extract and prune the matrix using the getrow function        */
   /* ============================================================= */

   nz_cnt = Nrows + 1;
   mat_indx[0] = nz_cnt;
   for ( i = 0; i < Nrows; i++ )
   {
      getrowfunc(getrowdata,1,&i,maxnnz_per_row,col_ind,col_val,&m);
      length = 0;
      for (j = 0; j < m; j++)
      {
         jnode = col_ind[j];
         if ( jnode != i ) length++;
         if ( jnode != i && epsilon > 0.0 )
         {
            dcompare1 = col_val[j] * col_val[j];
            if ( dcompare1 > 0.0 )
            {
               dcompare2 = ML_dabs((diagonal[i] * diagonal[jnode]));
               if ( dcompare1 >= epsilon * dcompare2 )
                  mat_indx[nz_cnt++] = col_ind[j];
            }
         }
         else if ( jnode != i && col_val[j] != 0.0 )
         {
            mat_indx[nz_cnt++] = col_ind[j];
         }
      }
      if ( m == 0 || (m == 1 && col_ind[0] == i)) bc_array[i] = 1;
      else                                        bc_array[i] = 0; 
      mat_indx[i+1] = nz_cnt;
      ML_sort(mat_indx[i+1]-mat_indx[i], mat_indx+mat_indx[i]);
   }
   if ( col_ind  != NULL ) ML_free(col_ind);
   if ( col_val  != NULL ) ML_free(col_val);
   if ( diagonal != NULL ) ML_free(diagonal);

   (*mat_indx_out)  = mat_indx;
   (*exp_Nrows_out) = exp_Nrows;
   (*bc_array_out)  = bc_array;
   return 0;
}

