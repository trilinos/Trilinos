/* ******************************************************************** */
/* ******************************************************************** */
/* Functions to create tentative prolongators (coupled aggregation)     */
/* ******************************************************************** */
/* Author        : Charles Tong                                         */
/* Organization  : Sandia National Laboratories                         */
/* Date          : August, 2000                                         */
/* ******************************************************************** */
/* Local Functions :                                                    */
/*    ML_Aggregate_CoarsenCoupledVBlock                                 */
/*    ML_Aggregate_CoarsenCoupledCore                                   */
/*    ML_Aggregate_Compress_Matrix                                      */
/*    ML_Aggregate_ExchangeData2                                        */
/*    ML_Aggregate_ComposeVBlockInfo                                    */
/*    ML_Aggregate_ComposeBlockCommInfo                                 */
/*    ML_Aggregate_ComposeRecvInfo                                      */
/*    ML_Aggregate_Form_Aggregates                                      */
/*    ML_Aggregate_PutInto_Aggregates                                   */
/* ******************************************************************** */

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "ml_aggregate.h"
#include "ml_lapack.h"

#define abs(x) (((x) > 0) ? x : (-(x)))

/* ******************************************************************** */
/* internal function defined later on in this file                      */
/* -------------------------------------------------------------------- */

int ML_Aggregate_Compress_Matrix(ML_GetrowFunc *getrow_obj, int *mat_indx, 
           int N_info, int *vinfo, ML_Comm *comm, int **new_mat_indx, 
           int *N_neighbors, int **neighbors, int **recv_leng,
           int **send_leng, int **send_list, int **recv_list);
int ML_Aggregate_CoarsenCoupledCore(ML_Aggregate *, ML_Comm *comm,
           int *amal_mat_indx, int *aggr_count, int **aggr_index2, 
           int N_neighbors, int *neighbors, int *recv_leng, int *send_leng,
           int *send_list,int *,int **); 
int  ML_Aggregate_ExchangeData2(char *recvbuf, char *sendbuf, 
           int N_neighbors, int *neighbors, int *recv_leng, 
           int *send_leng, int *recv_list, int Nrows, int msgid, 
           int datatype, ML_Comm *comm);
int ML_Aggregate_ComposeVBlockInfo(int Nrows, int exp_Nrows, 
           int *new_nvblocks, int **new_vblock_info, ML_Aggregate *ml_ag, 
           int mypid, ML_GetrowFunc *getrow_obj, ML_Comm *comm);
int ML_Aggregate_ComposeVBlockCommInfo(ML_GetrowFunc *getrow_obj, 
           int *mat_indx, int nvblocks, int *vblock_info, ML_Comm *comm, 
           int *new_N_neighbors, int **new_neighbors, int **new_recv_leng, 
           int **new_send_leng, int **new_send_list, int **new_recv_list);
int ML_Aggregate_ComposeRecvInfo(int nprocs, int mypid, int new_N_send,
           int *new_send_leng, int *new_send_neighbors, int *N_rcv, 
           int **recv_leng, int **recv_neighbors, ML_Comm *comm);
int ML_Aggregate_Form_Aggregates(char phaseID, int phaseAFlag, int Nrows, 
           int exp_Nrows, int *mat_indx, int *aggr_index, int *aggr_stat, 
           int *node_type, int *node_type2, int *order_array, 
           int *aggr_count_in, int *aggr_cnt_leng_in,
           int **aggr_cnt_array_in, int max_row_nnz, int min_agg_size, 
           int max_neigh_selected, int N_neighbors, int *neighbors, 
           int *send_leng, int *send_list, int *recv_leng, int *recv_list, 
           int *sendlist_proc, ML_Comm *comm);
int ML_Aggregate_PutInto_Aggregates(char phaseID, int attach_scheme, 
           int *mat_indx, int *aggr_index, int *aggr_stat, 
           int *aggr_count_in, int **aggr_cnt_array_in, int max_row_nnz, 
           int N_neighbors, int *neighbors, int *send_leng, int *send_list, 
           int *recv_leng, int *recv_list, ML_Comm *comm);

/* ******************************************************************** */
/* external functions called from this file                             */
/* -------------------------------------------------------------------- */

extern void ML_CSR_MSR_ML_memorydata_Destroy(void *data);
extern int  ML_randomize(int nlist, int *list);

/* ******************************************************************** */
/* local defines                                                        */
/* -------------------------------------------------------------------- */

#define ML_AGGR_READY      -11
#define ML_AGGR_NOTSEL     -12
#define ML_AGGR_SELECTED   -13
#define ML_AGGR_SELECTED2  -14
#define ML_AGGR_BDRY       -15
#define ML_AGGR_MINRANK      1
#define ML_AGGR_MAXLINK      2
#define ML_AGGR_INTERIOR     0
#define ML_AGGR_BORDER       1

/* ******************************************************************** */
/* ******************************************************************** */
/* ML_Aggregate_CoarsenCoupledVBlock subroutine.                        */
/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */

int ML_Aggregate_CoarsenCoupledVBlock( ML_Aggregate *ml_ag,
       ML_Operator *Amatrix, ML_Operator **Pmatrix, ML_Comm *comm)
{
   int     i, j, k, m, jj, jnode, index, index3, index4, offset, count; 
   int     max_count, nz_cnt, nbytes, length, level, diff_level;
   int     Nrows, exp_Nrows, *mat_indx=NULL, *amal_mat_indx, nvblocks;
   int     maxnnz_per_row=500, *vblock_info, *col_ind;
   int     N_neighbors, *neighbors, *recv_leng, *send_leng, *send_list;
   int     total_recv_leng, total_send_leng, msgtype, mypid, new_N_send;
   int     *new_send_neighbors, *new_send_list, *new_send_leng;
   int     new_N_recv, *new_recv_leng, *new_recv_neighbors, *int_buf;
   int     *int_buf2, *recv_list, nprocs, label;
   int     aggr_count, *aggr_index, *aggr_index2;
   int     *aggr_cnt_array, max_agg_size, **rows_in_aggs;
   int     Ncoarse, exp_Ncoarse, *new_ia, *new_ja, new_Nrows;
   int     num_PDE_eqns, nullspace_dim, lwork, info;
   double  *col_val, *diagonal=NULL, dcompare1, dcompare2, *new_val=NULL;
   double  epsilon, *dble_buf=NULL, *nullspace_vect=NULL, *qr_tmp=NULL;
   double  *tmp_vect=NULL, *work=NULL, *new_null=NULL, *comm_val=NULL;
   double  *dble_buf2=NULL;
   int     (*getrowfunc)(void *,int,int*,int,int*,double*,int*);
   void    *getrowdata;
   struct ML_CSR_MSRdata *csr_data;
   ML_Aggregate_Comm     *aggr_comm;
   ML_GetrowFunc         *getrow_obj;
   ML_CommInfoOP         *getrow_comm;

   /* ============================================================= */
   /* get machine and matrix information                            */
   /* ============================================================= */

   mypid          = comm->ML_mypid;
   nprocs         = comm->ML_nprocs;
   nullspace_dim  = ml_ag->nullspace_dim;
   nullspace_vect = ml_ag->nullspace_vect;
   Nrows          = Amatrix->outvec_leng;

   /* ============================================================= */
   /* initialize and update the threshold                           */
   /* ============================================================= */

   diff_level = ml_ag->begin_level - ml_ag->cur_level;
   if ( diff_level == 0 ) ml_ag->curr_threshold = ml_ag->threshold;
   epsilon = ml_ag->curr_threshold;
   ml_ag->curr_threshold *= 0.5;

#ifdef AGG_OUTPUT
   if ( mypid == 0 )
   {
      printf("ML_Aggregate_CoarsenCoupledVBlock : current level = %d\n",
                           ml_ag->cur_level);
      printf("ML_Aggregate_CoarsenCoupledVBlock : current eps = %e\n",
                           epsilon);
   }
#endif
   epsilon = epsilon * epsilon;

   /* ============================================================= */
   /* fetch the getrow function for the incoming matrix             */
   /* ============================================================= */

   getrow_obj = Amatrix->getrow;
   getrowfunc = getrow_obj->func_ptr;
   getrowdata = (void *) Amatrix;
   if ( getrowfunc == NULL )
   {
      printf("ML_Aggregate_CoarsenCoupledVBlock ERROR : null getrowfunc.\n");
      exit(-1);
   }

   /* ============================================================= */
   /* allocate initial temporary storage space for getrow           */
   /* also allocate space for storing the diagonal                  */
   /* ============================================================= */

   nbytes = maxnnz_per_row * sizeof( int );
   ML_memory_alloc((void**) &col_ind, nbytes, "ACA");
   nbytes = maxnnz_per_row * sizeof( double );
   ML_memory_alloc((void**) &col_val, nbytes, "ACB");
   if ( Nrows > 0 )
   {
      nbytes = Nrows * sizeof( double );
      ML_memory_alloc((void**) &diagonal, nbytes, "ACC");
   }
   else diagonal = NULL;

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
         ML_memory_free((void**) &col_ind);
         ML_memory_free((void**) &col_val);
         maxnnz_per_row = maxnnz_per_row * 2 + 1;
         nbytes = maxnnz_per_row * sizeof( int );
         ML_memory_alloc((void**) &col_ind, nbytes, "ACD");
         nbytes = maxnnz_per_row * sizeof( double );
         ML_memory_alloc((void**) &col_val,  nbytes, "ACE");
      }
      for ( j = 0; j < m; j++ )
      {
         if ( col_ind[j] > exp_Nrows ) exp_Nrows = col_ind[j];
         if ( col_ind[j] == i )        diagonal[i] = col_val[j];
      }
      count += m;
      if ( diagonal[i] == 0.0 )
      {
         printf("%d : CoarsenCoupledVBlock WARNING - diag %d is 0.\n",
                mypid,i);
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
      ML_memory_free((void**) &diagonal);
      diagonal = NULL;
   }
   if ( epsilon != 0.0 && exp_Nrows > 0 )
   {
      dble_buf = diagonal;
      nbytes = exp_Nrows * sizeof(double);
      if ( nbytes > 0 ) ML_memory_alloc((void**) &diagonal, nbytes, "ACF");
      else              diagonal = NULL;
      for ( i = 0; i < Nrows; i++ ) diagonal[i] = dble_buf[i];
      for ( i = Nrows; i < exp_Nrows; i++ ) diagonal[i] = 0.0;
      if ( dble_buf != NULL ) ML_memory_free((void**) &dble_buf);
      getrow_comm = getrow_obj->pre_comm;
      if ( getrow_comm != NULL )
         ML_exchange_bdry(diagonal,getrow_comm,Nrows,comm,ML_OVERWRITE);
   }

   /* ============================================================= */
   /* allocate temporary storage space for getrow                   */
   /* ============================================================= */

   nbytes = (count + 1) * sizeof( int );
   ML_memory_alloc((void**) &mat_indx, nbytes, "ACG");
   k = ML_Comm_GsumInt( comm, Nrows);
   m = ML_Comm_GsumInt( comm, count);

   k = ML_Comm_GsumInt( comm, Nrows);
   m = ML_Comm_GsumInt( comm, count);
#ifdef AGG_OUTPUT
   if ( mypid == 0 )
      printf("Aggregation(CVB) : Total nnz = %d (Nrows=%d)\n",m,k);
#endif
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
      for (j = 0; j < m; j++)
      {
         jnode = col_ind[j];
         if ( jnode != i && epsilon > 0.0 )
         {
            dcompare1 = col_val[j] * col_val[j];
            if ( dcompare1 > 0.0 )
            {
               dcompare2 = abs((diagonal[i] * diagonal[jnode]));
               if ( dcompare1 >= epsilon * dcompare2 )
                  mat_indx[nz_cnt++] = col_ind[j];
            }
         }
         else if ( jnode != i && col_val[j] != 0.0 )
         {
            mat_indx[nz_cnt++] = col_ind[j];
         }
      }
      mat_indx[i+1] = nz_cnt;
      ML_sort(mat_indx[i+1]-mat_indx[i], mat_indx+mat_indx[i]);
   }
   if ( col_ind  != NULL ) ML_memory_free((void**) &col_ind);
   if ( col_val  != NULL ) ML_memory_free((void**) &col_val);
   if ( diagonal != NULL ) ML_memory_free((void**) &diagonal);

   /* ============================================================= */
   /* get the off-processor variable block information              */
   /* ==> nvblocks, vblock_info (of length exp_Nrows)               */
   /*  (vblock_info[i] contains block number of node i)             */
   /* ============================================================= */

   ML_Aggregate_ComposeVBlockInfo(Nrows, exp_Nrows, &nvblocks,
                      &vblock_info, ml_ag, mypid, getrow_obj, comm);

   /* ============================================================= */
   /* compress the matrix using vblock information                  */
   /* ============================================================= */

   ML_Aggregate_Compress_Matrix(getrow_obj,mat_indx,nvblocks,vblock_info, 
                comm,&amal_mat_indx,&N_neighbors,&neighbors,&recv_leng,
                &send_leng, &send_list, &recv_list);

   /* ============================================================= */
   /* perform coarsening on the compressed matrix                   */
   /* ============================================================= */

   ML_Aggregate_CoarsenCoupledCore(ml_ag, comm, amal_mat_indx,
       &aggr_count, &aggr_index2, N_neighbors, neighbors, recv_leng,
       send_leng, send_list, recv_list, &aggr_cnt_array);
   if ( amal_mat_indx != NULL ) ML_memory_free( (void**) &amal_mat_indx );
   if ( neighbors     != NULL ) ML_memory_free( (void**) &neighbors );
   if ( recv_leng     != NULL ) ML_memory_free( (void**) &recv_leng );
   if ( send_leng     != NULL ) ML_memory_free( (void**) &send_leng );
   if ( send_list     != NULL ) ML_memory_free( (void**) &send_list );
   if ( recv_list     != NULL ) ML_memory_free( (void**) &recv_list );

   /* ------------------------------------------------------------- */
   /* compose a new communication info for A in view of blocks      */
   /* This is needed to send null spaces around between processors) */
   /* ------------------------------------------------------------- */

   ML_Aggregate_ComposeVBlockCommInfo(getrow_obj, mat_indx, nvblocks, 
           vblock_info, comm, &N_neighbors, &neighbors, 
           &recv_leng, &send_leng, &send_list, &recv_list);

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
      aggr_index[i] = aggr_index2[vblock_info[i]];
      if ( aggr_index[i] >= 0 && aggr_index[i] < aggr_count )
          aggr_cnt_array[aggr_index[i]]++;
   }   
   if ( aggr_index2 != NULL ) ML_memory_free( (void**) &aggr_index2 );
   aggr_index2 = NULL;
/*
for (i=0; i<exp_Nrows; i++)
   printf("%d : aggr_index[%4d] = %d\n", mypid, i, aggr_index[i]);
*/

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
   ML_Aggregate_ExchangeData2((char*) int_buf, (char*) int_buf2,
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
   if ( N_neighbors > 0 )
   {
      nbytes = N_neighbors * sizeof(int);
      if ( nbytes > 0 ) int_buf = (int *) ML_allocate( nbytes );
      else              int_buf = NULL;
      nbytes = Ncoarse * sizeof(int);
      if ( nbytes > 0 ) int_buf2 = (int *) ML_allocate( nbytes );
      else              int_buf2 = NULL;
      for ( i = 0; i < N_neighbors; i++ ) int_buf[i] = 0;

      /* ---------------------------------------------------------- */
      /* count which remote fine nodes belong to local aggregates   */
      /* in order to generate the communication pattern for         */
      /* the interpolation operator.                                */
      /* ---------------------------------------------------------- */

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

      /* ---------------------------------------------------------- */
      /* now the number of neighbors for P has been found, the next */
      /* step is to find the send_list and send_leng for the matvec */
      /* function for interpolation                                 */
      /* ---------------------------------------------------------- */

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
            printf("ML_Aggregate_CoupledVBlock : internal error (1).\n");
            exit(-1);
         }
      }
      else
      {
         new_send_leng = NULL;
         new_send_neighbors = NULL;
         new_send_list = NULL;
      } 
      if ( int_buf  != NULL ) free (int_buf);
      if ( int_buf2 != NULL ) free (int_buf2);
      ML_Aggregate_ComposeRecvInfo(nprocs, mypid, new_N_send, new_send_leng, 
          new_send_neighbors, &new_N_recv, &new_recv_leng, &new_recv_neighbors, 
          comm);
   }
   else
   {
      new_send_leng = NULL;
      new_send_neighbors = NULL;
      new_send_list = NULL;
   }

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
         printf("       requested = %d\n",aggr_cnt_array[i]*sizeof(int));
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
      ML_Aggregate_ExchangeData2((char*)dble_buf2,(char*) dble_buf,
            N_neighbors, neighbors, recv_leng, send_leng, recv_list, 
            Nrows, msgtype,length,comm);
      if ( dble_buf != NULL ) ML_memory_free((void**) &dble_buf);
   } 
   else dble_buf2 = NULL; 

   /* ------------------------------------------------------------- */
   /* compose the vblock information for block QR decomposition     */
   /* coming into this part, vblock_info[0:exp_Nrows-1] contains    */
   /* block numbers of individual equations.  Would like to change  */
   /* vblock_info[i] to be the blocksize node i is in.              */ 
   /* ------------------------------------------------------------- */

   if ( vblock_info != NULL ) ML_memory_free( (void**) &vblock_info );
   nvblocks = exp_Nrows;
   nbytes = nvblocks * sizeof(int);
   if ( nbytes > 0 ) ML_memory_alloc((void**) &vblock_info, nbytes, "ACZ");
   else              vblock_info = NULL;

   if ( ml_ag->nvblocks == 0 )
   {
      for ( i = 0; i < exp_Nrows; i++ ) vblock_info[i] = 1;
   } 
   else
   {
      nbytes = total_send_leng * sizeof(int);
      if ( nbytes > 0 ) int_buf = (int *) ML_allocate(nbytes);
      else              int_buf = NULL;
      count = 0;
      for ( i = 0; i < ml_ag->nvblocks; i++ )
         for ( j = 0; j < ml_ag->vblock_info[i]; j++ )
            vblock_info[count++] = i;

      count = 0;
      for ( i = 0; i < N_neighbors; i++ )
      {
         for ( j = 0; j < send_leng[i]; j++ )
         {
            int_buf[count] = mypid * 10000 + send_list[count];
            count++;
         }
      }
      msgtype = 12094;
      ML_Aggregate_ExchangeData2((char*)&vblock_info[Nrows],
            (char*) int_buf, N_neighbors, neighbors, recv_leng, 
            send_leng, recv_list, Nrows, msgtype,sizeof(int),comm);

      nvblocks = exp_Nrows;
      label = vblock_info[0];
      count = 1;
      index = 0;
      for ( i = 1; i < exp_Nrows; i++ )
      {
         if ( vblock_info[i] != label )
         {
            for ( j = index; j < i; j++ ) vblock_info[j] = count;
            index = i;
            label = vblock_info[i];
            count = 1;
         } else count++;
      }
      for ( j = index; j < exp_Nrows; j++ ) vblock_info[j] = count;
      if ( int_buf != NULL ) free( int_buf );
   }

   /* ------------------------------------------------------------- */
   /* perform block QR decomposition                                */
   /* ------------------------------------------------------------- */

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
            count = vblock_info[index];
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

      MLFORTRAN(dgeqrf)(&(aggr_cnt_array[i]), &nullspace_dim, qr_tmp,
               &(aggr_cnt_array[i]), tmp_vect, work, &lwork, &info);
      if (info != 0)
         pr_error("CoarsenCoupled ERROR : dgeqrf returned a non-zero\n");

      if (work[0] > lwork)
      {
         lwork=(int) work[0];
         ML_memory_free((void**) &work);
         ML_memory_alloc((void**) &work, sizeof(double)*lwork, "ACa");
      }
      else lwork=work[0];

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
      MLFORTRAN(dorgqr)(&(aggr_cnt_array[i]),&nullspace_dim,&nullspace_dim,
              qr_tmp, &(aggr_cnt_array[i]), tmp_vect, work, &lwork, &info);
      if (info != 0)
         pr_error("CoarsenCoupled ERRO : dorgqr returned a non-zero\n");

      if (work[0] > lwork)
      {
         lwork=(int) work[0];
         ML_memory_free((void**) &work);
         ML_memory_alloc((void**) &work, sizeof(double)*lwork, "ACb");
      }
      else lwork=work[0];

      /* ---------------------------------------------------------- */
      /* now copy Q over into the appropriate part of P:            */
      /* The rows of P get calculated out of order, so I assume the */
      /* Q is totally dense and use what I know of how big each Q   */
      /* will be to determine where in ia, ja, etc each nonzero in  */
      /* Q belongs.  If I did not assume this, I would have to keep */
      /* all of P in memory in order to determine where each entry  */
      /* should go                                                  */
      /* ---------------------------------------------------------- */

      for (j = 0; j < aggr_cnt_array[i]; j++)
      {
         index = rows_in_aggs[i][j];
         if ( index < Nrows )
         {
            index3 = new_ia[index];
            for (k = 0; k < nullspace_dim; k++)
            {
               new_ja [index3+k] = i * nullspace_dim + k;
               new_val[index3+k] = qr_tmp[ k*aggr_cnt_array[i]+j];
            }
         }
         else
         {
            index3 = (index - Nrows) * nullspace_dim;
            for (k = 0; k < nullspace_dim; k++)
               comm_val[index3+k] = qr_tmp[ k*aggr_cnt_array[i]+j];
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
   ML_Aggregate_ExchangeData2((char*)dble_buf,(char*) comm_val,
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
      if ( dcompare1 == 0.0 )
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
   (*Pmatrix) = ML_Operator_Create();
   ML_Operator_Set_ApplyFuncData(*Pmatrix, nullspace_dim*Ncoarse, Nrows,
                                 ML_EMPTY, csr_data, Nrows, NULL, 0);
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

   ML_memory_free((void**) &vblock_info );
   ML_memory_free((void**) &comm_val);
   ML_memory_free((void**) &mat_indx);
   ML_memory_free((void**) &neighbors);
   ML_memory_free((void**) &recv_leng);
   ML_memory_free((void**) &send_leng);
   ML_memory_free((void**) &send_list);
   ML_memory_free((void**) &recv_list);
   ML_memory_free((void**) &aggr_index);
   ML_memory_free((void**) &aggr_cnt_array);
   for (i = 0; i < aggr_count; i++) free(rows_in_aggs[i]);
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
   if ( N_neighbors > 0 )
   {
      ML_memory_free((void**) &new_recv_leng);
      ML_memory_free((void**) &new_recv_neighbors);
   }
   ML_memory_free((void**) &aggr_comm);
   return Ncoarse*nullspace_dim;
}

/* ******************************************************************** */
/* Coarsening with coupling according to Tuminaro                       */
/* -------------------------------------------------------------------- */
/* This algorithm goes as follow :                                      */
/*                                                                      */
/* 1) At the beginning of the aggregation phase, each processor picks a */
/*    point on the interprocessor boundary. Each processor then computes*/
/*    the graph distance between every other interprocessor boundary    */
/*    point and this first point. I was thinking that when computing    */
/*    these distances we only consider the graph corresponding to the   */
/*    interprocessor points.  We would need to handle the case where    */
/*    the interprocessor boundary is not completely connected. I've     */
/*    attached some matlab code that will hopefully make this clearer.  */
/* 2) When we create aggregates on the interprocessor boundary, we take */
/*    the point with the smallest distance among all valid points and   */
/*    we use this point to build the next aggregate.                    */
/* 3) When finished with the interprocessor boundary, we do the same    */
/*    thing in the interior. That is, we pick a point in the interior   */
/*    and then compute a graph distance between it and every other      */
/*    interior point.  We might want to choose this first point to be   */
/*    close to the already computed aggregates on the interprocessor    */
/*    boundary ... if this is not hard.                                 */
/* 4) When creating interior aggregates, we take the point with the     */
/*    smallest distance among valid points.                             */
/* -------------------------------------------------------------------- */

int ML_Aggregate_CoarsenCoupledCore(ML_Aggregate *ml_ag, ML_Comm *comm,
       int *mat_indx, int *aggr_count_out, int **aggr_index_out, 
       int N_neighbors, int *neighbors, int *recv_leng, int *send_leng, 
       int *send_list, int *recv_list, int **cnt_array)
{
   int     i, j, k, inode, jnode, nbytes, length, Nrows;
   int     aggr_count, mypid, *order_array;
   int     *aggr_index, *order_array2;
   int     *aggr_stat, aggr_cnt_leng, *aggr_cnt_array;
   int     seed_node, *node_type, *node_dist, node_left, *dist_array;
   int     exp_Nrows, *trackbc, max_dist, *sendlist_proc;
   int     phaseAFlag;
   int     max_length=0;
   int     attach_scheme, min_agg_size, max_neigh_selected;
   ML_Node       *node_head, *node_tail, *new_node;

   /* ============================================================= */
   /* get the machine information and matrix references             */
   /* ============================================================= */

   mypid = comm->ML_mypid;
   Nrows = mat_indx[0] - 1;
   exp_Nrows = Nrows;
   for ( i = 0; i < N_neighbors; i++ ) exp_Nrows += recv_leng[i];
   attach_scheme = ml_ag->attach_scheme;

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

   nbytes = Nrows * sizeof(int);
   if ( nbytes > 0 ) ML_memory_alloc((void**) &trackbc, nbytes, "Ag5");
   else              trackbc = NULL;
   for ( i = 0; i < Nrows; i++ ) 
   {
      length = mat_indx[i+1] - mat_indx[i];
      if ( length == 0 ) trackbc[i] = 1; else trackbc[i] = 0;
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
         free(new_node);
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
      if (trackbc[i] == 1) aggr_stat[i] = ML_AGGR_BDRY;
      else                 aggr_stat[i] = ML_AGGR_READY;
   }
   ML_memory_free( (void**) &trackbc );
   for ( i = 0; i < exp_Nrows; i++ ) aggr_index[i] = -1;
   for ( i = Nrows; i < exp_Nrows; i++ ) aggr_stat[i] = 0;

   /* ============================================================= */
   /* Phase 1 :                                                     */
   /*    This consists of two parts - aggregate border nodes first  */
   /*    followed by aggregating interior nodes.  This goes on      */
   /*    until all nodes are either selected or not selected.       */
   /* ============================================================= */

   min_agg_size = 0;
   max_neigh_selected = 0;
   phaseAFlag = 1;
   ML_Aggregate_Form_Aggregates('1', phaseAFlag, Nrows, exp_Nrows, mat_indx, 
        aggr_index, aggr_stat, node_type, node_type, order_array, &aggr_count, 
        &aggr_cnt_leng, &aggr_cnt_array, max_length, min_agg_size, 
        max_neigh_selected, N_neighbors, neighbors, send_leng, send_list, 
        recv_leng, recv_list, sendlist_proc, comm);

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

   ML_Aggregate_Form_Aggregates('2',phaseAFlag,Nrows,exp_Nrows,mat_indx, 
        aggr_index, aggr_stat, node_type, node_type, order_array, &aggr_count, 
        &aggr_cnt_leng, &aggr_cnt_array, max_length, min_agg_size, 
        max_neigh_selected, N_neighbors, neighbors, send_leng, send_list, 
        recv_leng, recv_list, sendlist_proc, comm);

   /* ============================================================= */
   /* Phase 3 :                                                     */
   /*    for all nodes, see if it can be aggregated into one of the */
   /*    existing LOCAL aggregates.                                 */
   /* ============================================================= */

   ML_Aggregate_PutInto_Aggregates('3', attach_scheme, mat_indx, 
        aggr_index, aggr_stat, &aggr_count, &aggr_cnt_array, 
        max_length, N_neighbors, neighbors, send_leng, send_list, 
        recv_leng, recv_list, comm);

   /* ============================================================= */
   /* Phase 4 :                                                     */
   /*    for all remaining nodes, form new aggregates               */
   /* ------------------------------------------------------------- */

   min_agg_size = 1;
   max_neigh_selected = 10000;
   phaseAFlag = 0;
   ML_Aggregate_Form_Aggregates('4',phaseAFlag,Nrows,exp_Nrows,mat_indx, 
        aggr_index, aggr_stat, node_type, NULL, order_array, &aggr_count, 
        &aggr_cnt_leng, &aggr_cnt_array, max_length, min_agg_size, 
        max_neigh_selected, N_neighbors, neighbors, send_leng, send_list, 
        recv_leng, recv_list, sendlist_proc, comm);

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
      if ( aggr_stat[i] != ML_AGGR_SELECTED )
      {
         printf("%d : ERROR (CC) : node %d not aggregated (%d)\n", mypid,
                i, aggr_stat[i]);
         for ( j = mat_indx[i]; j < mat_indx[i+1]; j++ )
            printf("%d : neighbors = %d %d %d\n", mypid, mat_indx[j], 
                   aggr_stat[mat_indx[j]], aggr_index[mat_indx[j]]);
         exit(1);
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

/* ******************************************************************** */
/* ******************************************************************** */
/*          Support subroutines                                         */
/* ******************************************************************** */
/* ******************************************************************** */

/* ******************************************************************** */
/* ******************************************************************** */
/* Compress a matrix into a block matrix                                */
/*    getrow_obj : for extracting communication information             */
/*    mat_indx   : pruned matrix                                        */
/*    N_info     : = Nrows                                              */
/*    vinfo      : vinfo[i] - block number of node i                    */
/* -------------------------------------------------------------------- */

int ML_Aggregate_Compress_Matrix(ML_GetrowFunc *getrow_obj, int *mat_indx, 
                int N_info, int *vinfo, ML_Comm *comm, int **new_mat_indx, 
                int *N_neighbors, int **neighbors, int **recv_leng,
                int **send_leng, int **send_list, int **recv_list)
{
   int   i, j, k, Nrows, nz_cnt, nbytes, *amal_mat_indx, LN_neighbors;
   int   *Lneighbors, *Lsend_leng, *Lrecv_leng, *Lsend_list, *Lrecv_list;
   int   *Aneighbors, *Asend_leng, *Arecv_leng, *Asend_list, *Arecv_list;
   int   AN_neighbors, total_send_leng, count, label, lcount, nvblk_local;
   int   nvblocks, *vblock_info, *vblock_info2, row, amal_count, index;
   int   total_recv_leng, *Lsend_count, mcount;
   char  *col_entered;

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
      Aneighbors = (int *) ML_allocate( nbytes );
      Arecv_leng = (int *) ML_allocate( nbytes );
      Asend_leng = (int *) ML_allocate( nbytes );
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
   if ( nbytes > 0 ) Asend_list = (int *) ML_allocate( nbytes );
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
   if ( nbytes > 0 ) Arecv_list = (int *) ML_allocate( nbytes );
   else              Arecv_list = NULL;
   count = 0;
   for ( i = 0; i < AN_neighbors; i++ )
   {
      for (j = 0; j < Arecv_leng[i]; j++)
         Arecv_list[count++] =
            getrow_obj->pre_comm->neighbors[i].rcv_list[j];
   }

   /* ------------------------------------------------------------- */
   /* compress neighbor processor information - Lneighbors          */
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
   /* compress send stuff                                           */
   /* ------------------------------------------------------------- */

   count = 0;
   for ( i = 0; i < LN_neighbors; i++ )
   {
      lcount = 0;
      if ( Asend_leng[i] > 0 ) 
      {
         label = vinfo[Asend_list[count++]];
         lcount++;
      }
      for ( j = 1; j < Asend_leng[i]; j++ )
      {
         index = Asend_list[count++]; 
         if ( vinfo[index] != label )
         {
            lcount++; 
            label = vinfo[index];
         }
      }
      Lsend_leng[i] = lcount;
   }
   total_send_leng = 0;
   for ( i = 0; i < LN_neighbors; i++ ) total_send_leng += Lsend_leng[i];
   nbytes = total_send_leng * sizeof( int );
   if ( nbytes > 0 ) ML_memory_alloc((void**) &Lsend_list, nbytes, "AgD");
   else              Lsend_list = NULL;
   if ( nbytes > 0 ) Lsend_count = (int *) ML_allocate( nbytes );
   else              Lsend_count = NULL;

   count = lcount = 0;
   for ( i = 0; i < LN_neighbors; i++ )
   {
      if ( Asend_leng[i] > 0 ) 
      {
         label = vinfo[Asend_list[count++]];
         Lsend_list[lcount++] = label;
         mcount = 1;
      }
      for ( j = 1; j < Asend_leng[i]; j++ )
      {
         index = Asend_list[count++]; 
         if ( vinfo[index] != label ) 
         { 
            label = vinfo[index];
            Lsend_count[lcount-1] = mcount;
            mcount = 1;
            Lsend_list[lcount++] = label;
         } else mcount++;
      }
      if ( Asend_leng[i] > 0 ) Lsend_count[lcount-1] = mcount;
   }
   free( Lsend_count );
   
   /* ------------------------------------------------------------- */
   /* compress recv stuff                                           */
   /* ------------------------------------------------------------- */

   count = Nrows;
   for ( i = 0; i < AN_neighbors; i++ )
   {
      lcount = 0;
      if ( Arecv_leng[i] > 0 ) 
      {
         label = vinfo[count++];
         lcount++;
      }
      for ( j = 1; j < Arecv_leng[i]; j++ )
      {
         if ( vinfo[count] != label )
         {
            lcount++; label = vinfo[count];
         }
         count++;
      }
      Lrecv_leng[i] = lcount;
   }
   total_recv_leng = 0;
   for ( i = 0; i < LN_neighbors; i++ ) total_recv_leng += Lrecv_leng[i];
   nbytes = total_recv_leng * sizeof( int );
   if ( nbytes > 0 ) ML_memory_alloc((void**) &Lrecv_list, nbytes, "AgE");
   else              Lrecv_list = NULL;

   count = lcount = 0;
   for ( i = 0; i < LN_neighbors; i++ )
   {
      if ( Arecv_leng[i] > 0 ) 
      {
         label = vinfo[Arecv_list[count++]];
         Lrecv_list[lcount++] = label;
      }
      for ( j = 1; j < Arecv_leng[i]; j++ )
      {
         index = Arecv_list[count++]; 
         if ( vinfo[index] != label ) 
         { 
            label = vinfo[index];
            Lrecv_list[lcount++] = label;
         }
      }
   }

   /* ------------------------------------------------------------- */
   /* allocate storage for block matrix                             */
   /* ------------------------------------------------------------- */

   nbytes = (nz_cnt + 1) * sizeof( int ); /* probably excessive */
   ML_memory_alloc((void**) &amal_mat_indx, nbytes, "AgF");

   /* ------------------------------------------------------------- */
   /* reformat block information                                    */
   /* ------------------------------------------------------------- */

   nvblocks = vinfo[N_info-1] + 1;
   vblock_info = (int *) ML_allocate(nvblocks * sizeof(int));
   for ( i = 0; i < nvblocks; i++ ) vblock_info[i] = 0;
   for ( i = 0; i < N_info; i++ ) vblock_info[vinfo[i]]++; 
  
   /* ------------------------------------------------------------- */
   /* allocate temporary storage for block information              */
   /* ------------------------------------------------------------- */

   vblock_info2 = (int *) ML_allocate(nvblocks * sizeof(int));
   vblock_info2[0] = vblock_info[0];
   for ( i = 1; i < nvblocks; i++ )
      vblock_info2[i] = vblock_info2[i-1] + vblock_info[i];
   for ( i = nvblocks-1; i >= 0; i-- )
      if ( vblock_info2[i] == Nrows ) {nvblk_local = i + 1; break;}

   /* ------------------------------------------------------------- */
   /* start compressing                                             */
   /* ------------------------------------------------------------- */

   amal_count = nvblk_local + 1;
   amal_mat_indx[0] = amal_count;
   row = 0;
   col_entered = (char *) ML_allocate(sizeof(char)*(1+ nvblocks) );
   if (col_entered == NULL) 
   {
      printf("Not enough space in ML_aggregate\n");
      exit(1);
   }
   for ( i = 0; i < nvblocks; i++) col_entered[i] = 'F';

   for ( i = 0; i < nvblk_local; i++) 
   {
      col_entered[i] = 'T';
      for ( j = 0; j < vblock_info[i]; j++) 
      {
         for ( k = mat_indx[row]; k < mat_indx[row+1]; k++) 
         {
            if ( mat_indx[k] < vblock_info2[0] ) index = 0;
            else
            {
               index=ML_sorted_search(mat_indx[k],nvblocks,vblock_info2);
               if ( index < 0 ) index = - index;
               else             index++;
            }
            if ( index < 0 || index >= nvblocks )
               printf("ERROR : in almalgamation %d => %d(%d).\n",mat_indx[k],
                       index,nvblocks);
            if (col_entered[index] == 'F') {
               amal_mat_indx[ amal_count++] = index;
               col_entered[index] = 'T';
            }
         }
         row++;
      }
      amal_mat_indx[i+1] = amal_count;
      col_entered[i] = 'F';
      for ( j = amal_mat_indx[i]; j < amal_mat_indx[i+1]; j++)
         col_entered[ amal_mat_indx[j]] = 'F';
   }

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

   free(col_entered);
   free(vblock_info);
   free(vblock_info2);
   free(Aneighbors);
   free(Asend_leng);
   free(Arecv_leng);
   free(Asend_list);
   free(Arecv_list);

   return 0;
}

/* ******************************************************************** */
/* Exchange data between processors given communication information     */
/* -------------------------------------------------------------------- */

int ML_Aggregate_ExchangeData2(char *recvbuf,char *sendbuf,int N_neighbors,
              int *neighbors,int *recv_leng,int *send_leng,int *recv_list,
              int Nrows, int msgid, int datatype, ML_Comm *comm)
{
   int     i, nbytes, fromproc, length, typeleng, msgtype, offset;
   int     total_recv_leng, *int_array, *iarray;
   char    *char_array, *carray;
   double  *dble_array, *darray;
   USR_REQ *Request;

/*
static int flag=0;
int j, count;
if (flag == 0)
{
   flag = 1;
   count = 0;
   for (i=0; i <N_neighbors; i++)
      for (j=0; j <recv_leng[i]; j++)
         printf("recv_list[%4d,%4d] = %d\n", i, j, recv_list[count++]);
}
*/
   
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
                           &msgtype, comm->USR_comm, (void *) &Request[i] );
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
                             &msgtype, comm->USR_comm, (void *) &Request[i] );
      offset += recv_leng[i];
   }
   if ( Request != NULL ) free( Request );

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
                           free(char_array);
                           break;
         case ML_INT     : nbytes = total_recv_leng * sizeof(int);
                           int_array = (int *) ML_allocate(nbytes);
                           iarray = (int *) recvbuf;
                           for ( i = 0; i < total_recv_leng; i++ )
                              int_array[recv_list[i]-Nrows] = iarray[i];
                           for ( i = 0; i < total_recv_leng; i++ )
                              iarray[i] = int_array[i];
                           free(int_array);
                           break;
         case ML_DOUBLE  : nbytes = total_recv_leng * sizeof(double);
                           dble_array = (double *) ML_allocate(nbytes);
                           darray = (double *) recvbuf;
                           for ( i = 0; i < total_recv_leng; i++ )
                              dble_array[recv_list[i]-Nrows] = darray[i];
                           for ( i = 0; i < total_recv_leng; i++ )
                              darray[i] = dble_array[i];
                           free(dble_array);
                           break;
      }
   }
   return 0;
}

/* ******************************************************************** */
/* get the off-processor variable block information                     */
/* ==> nvblocks, vblock_info (of length exp_Nrows)                      */
/*     vblock_info[i] - block number of equation i                      */
/* -------------------------------------------------------------------- */

int ML_Aggregate_ComposeVBlockInfo(int Nrows, int exp_Nrows, 
          int *new_nvblocks, int **new_vblock_info, ML_Aggregate *ml_ag, 
          int mypid, ML_GetrowFunc *getrow_obj, ML_Comm *comm)
{
   int           i, j, nvblocks, nbytes, count, *vblock_info;
   double        *dble_buf;
   ML_CommInfoOP *getrow_comm;

   nvblocks = exp_Nrows;
   nbytes   = nvblocks * sizeof(int);
   ML_memory_alloc((void**) &vblock_info, nbytes, "VB1");
   if ( ml_ag->nvblocks == 0 )
   {
      for ( i = 0; i < exp_Nrows; i++ ) vblock_info[i] = i;
   } 
   else
   {
      nbytes = exp_Nrows * sizeof(double);
      dble_buf = (double *) ML_allocate(nbytes);
      count = 0;
      for ( i = 0; i < ml_ag->nvblocks; i++ )
      {
         for ( j = 0; j < ml_ag->vblock_info[i]; j++ )
         {
            vblock_info[count] = i;
            dble_buf[count++] = mypid * 10000 + i;
         }
      }
      getrow_comm = getrow_obj->pre_comm;
      if ( getrow_comm != NULL )
         ML_exchange_bdry(dble_buf,getrow_comm,Nrows,comm,ML_OVERWRITE);
      count = ml_ag->nvblocks - 1;
      for ( i = Nrows; i < exp_Nrows; i++ )
      {
         if ( Nrows == 0 || dble_buf[i] != dble_buf[i-1] ) count++;
         vblock_info[i] = count;
      }
      count++;
      free( dble_buf );
   }
   (*new_nvblocks) = nvblocks;
   (*new_vblock_info) = vblock_info;
   return 0;
}

/* ******************************************************************** */
/* compose send and receive information in view of block                */
/* -------------------------------------------------------------------- */

int ML_Aggregate_ComposeVBlockCommInfo(ML_GetrowFunc *getrow_obj, 
           int *mat_indx, int nvblocks, int *vblock_info, ML_Comm *comm, 
           int *out_N_neighbors, int **out_neighbors, 
           int **out_recv_leng, int **out_send_leng, int **out_send_list,
           int **out_recv_list)
{
   int      i, j, k, N_neighbors, *send_leng, *recv_leng, *send_list;
   int      *neighbors, nbytes, total_send_leng, total_recv_leng;
   int      index, count, count2, label, mypid, nprocs, msgtype;
   int      procnum, length, offset, *ind_array2, *recv_list; 
   int      new_N_neighbors, *new_send_leng, *new_send_list;
   int      *new_recv_list, *ext_mat_indx;
   int      *new_recv_leng, *new_neighbors, Nrows, *ind_array;
   int      proc_begin, toproc, fromproc, unitflag;
   USR_REQ  *Request;

   /* ----------------------------------------------------------------- */
   /* get machine information                                           */
   /* ----------------------------------------------------------------- */

   Nrows = mat_indx[0] - 1;
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
   mypid  = comm->ML_mypid;

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
   else neighbors = send_leng = NULL;

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

   unitflag = 1;
   for ( i = 0; i < nvblocks; i++ )
      if ( vblock_info[i] != i ) {unitflag = 0; break;}
   if ( unitflag == 1 )
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
      for ( j = 0; j < send_leng[i]; j++ )
      {
         index = send_list[count];
         label = vblock_info[index];

         /* -------------------------------------- */
         /* count number of equations in the group */
         /* -------------------------------------- */

         if ( label >= 0 )
         {
            while (label == vblock_info[index] && index < Nrows) 
            {
               vblock_info[index] = - vblock_info[index] - 1;
               index++;
               count2++;
            }
            index = send_list[count++] - 1;
            while (label == vblock_info[index] && index >= 0) 
            {
               vblock_info[index] = - vblock_info[index] - 1;
               index--; count2++;
            }
         }
      }
      new_send_leng[i] = count2;
      for ( k = 0; k < nvblocks; k++ ) 
         if ( vblock_info[k] < 0 )
            vblock_info[k] = - vblock_info[k] - 1; 
   }

   nbytes = total_send_leng * sizeof( int );
   if ( nbytes > 0 ) ML_memory_alloc((void**) &ind_array,nbytes,"CI7");

   total_send_leng = 0;
   for ( i = 0; i < new_N_neighbors; i++ )
      total_send_leng += new_send_leng[i];
   nbytes = total_send_leng * sizeof( int );
   if ( nbytes > 0 ) ML_memory_alloc((void**) &new_send_list,nbytes,"CI6");

   count = count2 = proc_begin = 0;
   for ( i = 0; i < N_neighbors; i++ )
   {
      new_neighbors[i] = neighbors[i];
      for ( j = 0; j < send_leng[i]; j++ )
      {
         index = send_list[count];
         label = vblock_info[index];
         if ( label >= 0 ) 
         {
            while (label == vblock_info[index] && index >= 0) { index--; }
            index++;
            while (label == vblock_info[index] && index < Nrows) 
            {
               vblock_info[index] = - vblock_info[index] - 1;
               if ( index == send_list[count] ) 
                  ind_array[count] = count2 - proc_begin;
               new_send_list[count2++] = index++;
            }
         } 
         else
         {
            for ( k = proc_begin; k < count2; k++ )
            {
               if (new_send_list[k] == index) 
               {
                  ind_array[count] = k - proc_begin;
                  break;
               }
            }
         }
         count++;
      }
      proc_begin += send_leng[i];
      for ( k = 0; k < nvblocks; k++ ) 
         if ( vblock_info[k] < 0 )
            vblock_info[k] = - vblock_info[k] - 1; 
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
                &msgtype, comm->USR_comm, (void *) &Request[i] );
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
               &fromproc,&msgtype,comm->USR_comm,(void *) &Request[i]);
   }

   /* ----------------------------------------------------------------- */
   /* set up true external indices to be shipped to receive processors  */
   /* (true in view of that num_PDE_eqns can be > 1)                    */
   /* ----------------------------------------------------------------- */

   total_recv_leng = 0;
   for (i = 0; i < N_neighbors; i++) total_recv_leng += recv_leng[i];
   nbytes = total_recv_leng * sizeof(int);
   if ( nbytes > 0 ) ML_memory_alloc((void**)&ind_array2, nbytes, "RI1");
   else              ind_array2 = NULL;

   offset = 0;
   for (i = 0; i < N_neighbors; i++) 
   {
      msgtype = 2000;
      length = recv_leng[i] * sizeof( int );
      procnum = neighbors[i];
      comm->USR_irecvbytes((void *) &(ind_array2[offset]),length,&procnum,
                           &msgtype, comm->USR_comm, Request+i);
      offset += recv_leng[i];
   }
   offset = 0;
   for (i = 0; i < N_neighbors; i++) 
   {
      msgtype = 2000;
      length = send_leng[i] * sizeof( int );
      procnum = neighbors[i];
      comm->USR_sendbytes((void *) &(ind_array[offset]),length,procnum, 
                              msgtype, comm->USR_comm);
      offset += send_leng[i];
   }
   offset = 0;
   for (i = 0; i < N_neighbors; i++) 
   {
      msgtype = 2000;
      length = recv_leng[i] * sizeof( int );
      procnum = neighbors[i];
      comm->USR_cheapwaitbytes((void *) &(ind_array2[offset]),length,&procnum,
                           &msgtype, comm->USR_comm, Request+i);
      for (j = 0; j < recv_leng[i]; j++) ind_array2[offset+j] += offset; 
      offset += recv_leng[i];
   }
   if ( N_neighbors > 0 ) free( Request );
   if ( ind_array != NULL ) ML_memory_free( (void**) &ind_array);

   /* ----------------------------------------------------------------- */
   /* adjust the matrix to reflect the index remapping                  */
   /* ----------------------------------------------------------------- */

   count = Nrows + 1;
   for ( i = 0; i < Nrows; i++ )
   {
      for ( j = mat_indx[i]; j < mat_indx[i+1]; j++ )
         if ( mat_indx[j] >= Nrows ) count++;
   }
   nbytes = count * sizeof(int);
   if ( nbytes > 0 ) ML_memory_alloc((void**)&ext_mat_indx, nbytes, "RI4");
   else              ext_mat_indx = NULL;
   if ( nbytes > 0 ) ext_mat_indx[0] = Nrows + 1;

   count = 0;
   for ( i = 0; i < Nrows; i++ )
   {
      for ( j = mat_indx[i]; j < mat_indx[i+1]; j++ )
         if ( mat_indx[j] >= Nrows ) 
            ext_mat_indx[j] = ind_array2[mat_indx[j]-Nrows] + Nrows;
   }

   nbytes = total_recv_leng * sizeof( int );
   if ( nbytes > 0 ) ML_memory_alloc((void**)&new_recv_list, nbytes, "RI3");
   else              new_recv_leng = NULL;
   for ( i = 0; i < total_recv_leng; i++ ) new_recv_list[i] = i;

   ML_memory_free((void**) &ext_mat_indx);
   ML_memory_free((void**) &ind_array2);
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

/* ******************************************************************** */
/* compose receive information from send information                    */
/* -------------------------------------------------------------------- */

int ML_Aggregate_ComposeRecvInfo(int nprocs, int mypid, int new_N_send,
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
      free( int_buf2 );
      new_N_rcv = int_buf[mypid];
      free( int_buf );

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
                   &msgtype, comm->USR_comm, (void *) &Request[i] );
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
                  &fromproc,&msgtype,comm->USR_comm,(void *) &Request[i]);
         new_recv_neighbors[i] = fromproc;
      }
      ML_az_sort( new_recv_neighbors, new_N_rcv, new_recv_leng, NULL);
      if ( new_N_rcv > 0 ) free( Request );

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

/* ******************************************************************** */
/* form new aggregates                                                  */
/* -------------------------------------------------------------------- */

int ML_Aggregate_Form_Aggregates(char phaseID, int phaseAFlag, int Nrows, 
        int exp_Nrows, int *mat_indx, int *aggr_index, int *aggr_stat, 
        int *node_type, int *node_type2, int *order_array, int *aggr_count_in, 
        int *aggr_cnt_leng_in,
        int **aggr_cnt_array_in, int max_row_nnz, int min_agg_size, 
        int max_neigh_selected, int N_neighbors, int *neighbors, 
        int *send_leng, int *send_list, int *recv_leng, int *recv_list, 
        int *sendlist_proc, ML_Comm *comm)
{
   int           i, j, k, m, index, inode, inode2, jnode, mdiff, count;
   int           mypid, msgtype, procnum, *com_buf, *com_buf2, nready;
   int           nbytes, total_send_leng, total_recv_leng, select_flag;
   int           aggr_count, *aggr_cnt_array, loop_flag, *itmp_array;
   int           old_aggr_count, nneigh_selected, node_type_scalar;
   int           nselected, total_nselected, total_Nrows, nwaiting;
   int           total_nwaiting, total_aggr_count, aggr_cnt_leng; 
   ML_SuperNode  *supernode;

   /* ---------------------------------------------------------------- */
   /* get inputs and compute local parameters                          */
   /* ---------------------------------------------------------------- */

   mypid = comm->ML_mypid;
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
      if ( aggr_stat[inode] != ML_AGGR_SELECTED ) 
         aggr_stat[inode] = ML_AGGR_READY;

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
               mdiff = index - Nrows;

               /* ---------------------------------------------------- */
               /* search for the processor the node is coming from     */
               /* ---------------------------------------------------- */

               for ( k = 0; k <= N_neighbors; k++ )
                  if ( mdiff < sendlist_proc[k] ) break;

               /* ---------------------------------------------------- */
               /* if the processor number < mypid, tag it with the     */
               /* neighbor processor with the smallest rank            */
               /* ---------------------------------------------------- */

               if ( k != 0 && neighbors[k-1] < mypid )
                     aggr_stat[inode] = neighbors[k-1];
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
      ML_Aggregate_ExchangeData2((char*)&aggr_stat[Nrows],(char*) com_buf,
                   N_neighbors, neighbors, recv_leng, send_leng, recv_list,
                   Nrows, msgtype, ML_INT, comm);

      /* ------------------------------------------------------------- */
      /* begin aggregating the border nodes                            */
      /* ------------------------------------------------------------- */

      loop_flag = 1;

      while ( loop_flag != 0 )
      {
         old_aggr_count = aggr_count;

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
                        index = mat_indx[j] - Nrows;
                        if ( index >= 0)
                        {
                           count = 0;
                           for (k = 0; k < N_neighbors; k++ )
                           {
                              if ( index < (count+recv_leng[k]) ) break;
                              count += recv_leng[k];
                           }
                           if (neighbors[k] > procnum) procnum = neighbors[k];
                        }
                     }
                     if (procnum==mypid) aggr_stat[inode] = ML_AGGR_NOTSEL;
                     else                aggr_stat[inode] = procnum;
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

#ifdef AGG_OUTPUT2
printf("%d : Phase %cA - aggr %4d has ", mypid, phaseID, aggr_count);
#endif
                  for ( j = 0; j < supernode->length; j++ )
                  {
                     jnode = supernode->list[j];
                     aggr_stat[jnode] = ML_AGGR_SELECTED;
                     aggr_index[jnode] = aggr_count;
#ifdef AGG_OUTPUT2
printf("%d ", jnode);
#endif
                     if ( jnode < Nrows )
                     {
                        for (k=mat_indx[jnode]; k<mat_indx[jnode+1]; k++)
                           if (aggr_stat[mat_indx[k]] == ML_AGGR_READY)
                              aggr_stat[mat_indx[k]] = ML_AGGR_NOTSEL;
                     }
                  }
#ifdef AGG_OUTPUT2
printf("\n");
#endif
                  /* ------------------------------------------------- */
                  /* stretch aggr_cnt_array, if needed                 */
                  /* ------------------------------------------------- */

                  aggr_cnt_array[aggr_count++] = supernode->length;
                  if ( aggr_count >= aggr_cnt_leng )
                  {
                     itmp_array = aggr_cnt_array;
printf("%d : AGGR_LENG : %d", mypid, aggr_cnt_leng);
                     aggr_cnt_leng = aggr_cnt_leng * 6 / 5 + 1;
printf(" %d \n", aggr_cnt_leng);
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
         ML_Aggregate_ExchangeData2((char*)com_buf2,(char*)&aggr_stat[Nrows],
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
         ML_Aggregate_ExchangeData2((char*)&aggr_stat[Nrows],(char*)com_buf,
            N_neighbors,neighbors,recv_leng,send_leng,recv_list,Nrows,
            msgtype,ML_INT,comm);

         /* ---------------------------------------------------------- */
         /* update my waiting nodes' status                            */
         /* ---------------------------------------------------------- */

         for ( inode = 0; inode < Nrows; inode++ )
         {
            if ( aggr_stat[inode] >= 0 )
            {
               procnum = 100000;
               for (jnode=mat_indx[inode];jnode<mat_indx[inode+1];jnode++)
               {
                  index = mat_indx[jnode];
                  mdiff = index - Nrows;
                  if ( mdiff >= 0 )
                  {
                     for ( k = 0; k <= N_neighbors; k++ )
                        if ( mdiff < sendlist_proc[k] ) break;
   
                     if ( aggr_stat[index] == ML_AGGR_READY )
                        procnum = neighbors[k-1];
                     else if (aggr_stat[index] >= 0 &&
                              aggr_stat[index] < mypid)
                     {
                        if ( neighbors[k-1] < procnum )
                           procnum = neighbors[k-1];
                     }
                  }
               }
               if ( procnum == 100000 ) aggr_stat[inode] = ML_AGGR_READY;
               else                     aggr_stat[inode] = procnum;
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
         ML_Aggregate_ExchangeData2((char*)&aggr_stat[Nrows],(char*)com_buf,
            N_neighbors,neighbors,recv_leng,send_leng,recv_list,Nrows,
            msgtype, ML_INT,comm);

#ifdef AGG_OUTPUT
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
         if ( mypid == 0 )
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
#endif

         /* ---------------------------------------------------------- */
         /* check to see if further loop is needed                     */
         /* ---------------------------------------------------------- */

         if ( total_nwaiting == 0 ) loop_flag = 0;
         else
         {
            k = aggr_count - old_aggr_count;
            m = ML_Comm_GsumInt( comm, k);
            if ( m == 0 ) loop_flag = 0;
            else          loop_flag++;
         }
      }
      if ( com_buf  != NULL ) free( com_buf ); 
      if ( com_buf2 != NULL ) free( com_buf2 ); 
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
#ifdef AGG_OUTPUT2
printf("%d : Phase %cB - aggr %4d has ", mypid, phaseID, aggr_count);
#endif
            for ( j = 0; j < supernode->length; j++ )
            {
               jnode = supernode->list[j];
#ifdef AGG_OUTPUT2
printf("%d ", jnode);
#endif
               aggr_stat[jnode] = ML_AGGR_SELECTED;
               aggr_index[jnode] = aggr_count;
            }
#ifdef AGG_OUTPUT2
printf("\n");
#endif
            aggr_cnt_array[aggr_count++] = supernode->length;
            if ( aggr_count >= aggr_cnt_leng )
            {
               itmp_array = aggr_cnt_array;
printf("%d : AGGR_LENG : %d", mypid, aggr_cnt_leng);
               aggr_cnt_leng = aggr_cnt_leng * 6 / 5 + 1;
printf(" %d \n", aggr_cnt_leng);
               nbytes = aggr_cnt_leng * sizeof( int );
               ML_memory_alloc((void**)&aggr_cnt_array,nbytes,"Agf");
               for ( k = 0; k < aggr_count; k++ )
                  aggr_cnt_array[k] = itmp_array[k];
               ML_memory_free((void**)&itmp_array);
            }
         }
      }
   }
   free( supernode->list );
   free( supernode );

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
   ML_Aggregate_ExchangeData2((char*)&aggr_stat[Nrows],(char*) com_buf,
      N_neighbors, neighbors,recv_leng,send_leng,recv_list,Nrows,msgtype,
      ML_INT, comm);
   if ( com_buf != NULL ) free( com_buf );

#ifdef AGG_OUTPUT
   /* ---------------------------------------------------------------- */
   /* output information about aggregation progress                    */
   /* ---------------------------------------------------------------- */

   nselected = 0;
   for (i = 0; i < Nrows; i++)
      if (aggr_stat[i] == ML_AGGR_SELECTED) nselected++;
   total_nselected = ML_Comm_GsumInt( comm, nselected);
   total_Nrows = ML_Comm_GsumInt( comm, Nrows);
   total_aggr_count = ML_Comm_GsumInt( comm, aggr_count );
   if ( mypid == 0 )
   {
      printf("Aggregation(CC) : Phase %cB - nodes aggregated = %d(%d)\n",
             phaseID, total_nselected, total_Nrows);
      printf("Aggregation(CC) : Phase %cB - total aggregates = %d\n",
             phaseID, total_aggr_count);
   }
#endif

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

/* ******************************************************************** */
/* put nodes into existing local aggregates                             */
/* -------------------------------------------------------------------- */

int ML_Aggregate_PutInto_Aggregates(char phaseID, int attach_scheme, 
        int *mat_indx, int *aggr_index, int *aggr_stat, int *aggr_count_in, 
        int **aggr_cnt_array_in, int max_row_nnz, int N_neighbors, 
        int *neighbors, int *send_leng, int *send_list, int *recv_leng, 
        int *recv_list, ML_Comm *comm)
{
   int          i, k, m, inode, jnode, index, mincount, select_flag;
   int          maxcount, length, *int_array, *int_array2, *com_buf;
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
            maxcount = 0;
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
            if ( length > 0 ) free( int_array );
            if ( length > 0 ) free( int_array2 );
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
   ML_Aggregate_ExchangeData2((char*)&aggr_stat[Nrows],(char*) com_buf,
        N_neighbors,neighbors,recv_leng,send_leng,recv_list,Nrows,msgtype,
        ML_INT,comm);
   if ( com_buf != NULL ) free( com_buf );

#ifdef AGG_OUTPUT

   /* ----------------------------------------------------------------- */
   /* output information about aggregation progress                     */
   /* ----------------------------------------------------------------- */

   nselected = 0;
   for (i = 0; i < Nrows; i++)
      if (aggr_stat[i] == ML_AGGR_SELECTED) nselected++;
   total_nselected = ML_Comm_GsumInt( comm, nselected);
   total_Nrows = ML_Comm_GsumInt( comm, Nrows);
   total_aggr_count = ML_Comm_GsumInt( comm, aggr_count );
   if ( mypid == 0 )
   {
      printf("Aggregation(CC) : Phase %c  - nodes aggregated = %d(%d)\n",
             phaseID, total_nselected, total_Nrows);
      printf("Aggregation(CC) : Phase %c  - total aggregates = %d\n",
             phaseID, total_aggr_count);
   }
#endif

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

#ifdef PHASE5
   /* ============================================================= */
   /* Phase 4 :                                                     */
   /*    for all remaining nodes, see if they can be aggregated     */
   /*    into one of the existing remote aggregates.                */
   /* ============================================================= */

   /* ------------------------------------------------------------- */
   /* communicate the index information for this phase              */
   /* ------------------------------------------------------------- */

   nbytes = total_send_leng * sizeof(int);
   if ( nbytes > 0 ) com_buf  = (int *) ML_allocate( nbytes );
   nbytes = total_recv_leng * sizeof(int);
   if ( nbytes > 0 ) com_buf2 = (int *) ML_allocate( nbytes );
   for ( i = 0; i < total_send_leng; i++ ) 
   {
      com_buf[i] = aggr_index[send_list[i]];
   }
   msgtype = 49934;
   ML_Aggregate_ExchangeData2((char*)com_buf2,(char*) com_buf,
        N_neighbors,neighbors,recv_leng,send_leng,recv_list,Nrows,
        msgtype,ML_INT,comm);

   /* ------------------------------------------------------------- */
   /* reset unselected nodes to be READY                            */
   /* ------------------------------------------------------------- */

   for ( inode = 0; inode < Nrows; inode++ )
   {
      if ( aggr_stat[inode] != ML_AGGR_SELECTED &&
           aggr_stat[inode] != ML_AGGR_BDRY ) 
         aggr_stat[inode] = ML_AGGR_READY;
   }

   /* ------------------------------------------------------------- */
   /* begin phase 4 aggregation                                     */
   /* ------------------------------------------------------------- */

   for ( inode = 0; inode < Nrows; inode++ )
   {
      /* ---------------------------------------------------------- */
      /* if the node in question is either READY or NOTSEL          */
      /* ---------------------------------------------------------- */

      if ((aggr_stat[inode] == ML_AGGR_READY || 
          aggr_stat[inode] == ML_AGGR_NOTSEL) && node_type[inode] == 1) 
      { 
         for (jnode=mat_indx[inode]; jnode<mat_indx[inode]; jnode++)
         {
            index = mat_indx[inode];
            if (index >= Nrows && aggr_stat[index]==ML_AGGR_SELECTED &&
                com_buf2[index-Nrows] >= 0)
            {
               mdiff = index - Nrows;
               for ( k = 0; k <= N_neighbors; k++ )
               if ( mdiff < sendlist_proc[k] ) break;
               aggr_stat[inode]  = mdiff - sendlist_proc[k-1];
               aggr_index[inode] = - 100 - neighbors[k-1];
               break;
            }
         }
      }
   }
 
   /* ------------------------------------------------------------- */
   /* communicate the index information for this phase              */
   /* ------------------------------------------------------------- */

   for ( i = 0; i < total_send_leng; i++ ) 
   {
      com_buf[i] = aggr_stat[send_list[i]];
   }
   msgtype = 49935;
   ML_Aggregate_ExchangeData2((char*)&aggr_stat[Nrows],(char*) com_buf,
        N_neighbors,neighbors,recv_leng,send_leng,recv_list,Nrows,
        msgtype,ML_INT,comm);

   for ( i = 0; i < total_send_leng; i++ ) 
   {
      com_buf[i] = aggr_index[send_list[i]];
   }
   msgtype = 49936;
   ML_Aggregate_ExchangeData2((char*)com_buf2,(char*) com_buf,
        N_neighbors,neighbors,recv_leng,send_leng,recv_list,Nrows,
        msgtype,ML_INT,comm);

   free( com_buf );

   /* ------------------------------------------------------------- */
   /* process the incoming data                                     */
   /* ------------------------------------------------------------- */

   length = 0;
   for ( i = 0; i < N_neighbors; i++ )
   {
      for ( j = 0; j < recv_leng[i]; j++ )
      {
         if ( aggr_stat[Nrows+length+j] >= 0 &&
              (-com_buf2[length+j]+100) == mypid )
         {
            index = aggr_stat[Nrows+length+j]; /* intra-proc offset */
            for ( k = 0; k < i; k++ ) index += send_leng[k];
            index = aggr_index[send_list[index]];
            aggr_cnt_array[index]++;
            aggr_index[Nrows+length+j] = index;
         }
      }
      length += recv_leng[i];
   } 
   free( com_buf2 );

   /* ------------------------------------------------------------- */
   /* set the status of the selected node to the SELECTED state     */
   /* ------------------------------------------------------------- */

   for ( inode = 0; inode < Nrows; inode++ )
   {
      if ( aggr_stat[inode] >= 0 ) aggr_stat[inode] = ML_AGGR_SELECTED;
   }

#ifdef AGG_OUTPUT

   /* ------------------------------------------------------------- */
   /* output information about aggregation progress                 */
   /* ------------------------------------------------------------- */

   m = 0;
   for (i = 0; i < Nrows; i++) if (aggr_stat[i] == ML_AGGR_SELECTED) m++;
   k = ML_Comm_GsumInt( comm, m);
   m = ML_Comm_GsumInt( comm, Nrows);
   j = ML_Comm_GsumInt( comm, aggr_count );
   if ( mypid == 0 )
   {
      printf("Aggregation(CC) : Phase 4  - nodes aggregated = %d (%d)\n",k,m);
      printf("Aggregation(CC) : Phase 4  - total aggregates = %d \n",j);
   }
#endif
#endif


