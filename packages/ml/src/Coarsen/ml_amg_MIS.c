/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* ******************************************************************** */
/* Functions to create AMG prolongators                                 */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL)                                  */
/* Date          : September, 2000                                      */
/* ******************************************************************** */
/* ******************************************************************** */

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "ml_operator.h"
#include "ml_amg.h"
#include "ml_aggregate.h"
#include "ml_lapack.h"

int ML_AMG_LabelVertices(int vlist_cnt, int *vlist, char Vtype,
                          char *vertex_state, char *vertex_type,
                          int nvertices, int *rptr, int *cptr, 
                          int myrank, int **proclist, int send_cnt, 
                          int **send_buf, int *send_proc, int *send_leng,
                          int recv_cnt, int **recv_buf, int *recv_proc, 
                          int *recv_leng, int **recv_list, int msgtype, 
                          ML_Comm *comm, int aggr_index[]);

#define dabs(x) (((x) >= 0) ? x : -x)

/* ******************************************************************** */
/* construct the tentative prolongator allowing aggregate to cross      */
/* processor boundaries                                                 */
/* -------------------------------------------------------------------- */

int ML_AMG_CoarsenMIS( ML_AMG *ml_amg, ML_Operator *Amatrix, 
                       ML_Operator **Pmatrix, ML_Comm *comm)
{
   int     i, j, k, m, nbytes, offset, msgtype, Nrows, index, mypid, nprocs;
   int     num_PDE_eqns, count, Ncoarse, *int_buf, *int_buf2, max_element;
   int     level, printflag, diff_level, exp_Nrows, exp_Ncoarse;
   int     N_neighbors, *neighbors=NULL, *send_list=NULL, *sendlist_proc=NULL;
   int     *recv_list=NULL, total_recv_leng=0, total_send_leng=0;
   int     total_nnz, *new_ia = NULL, *new_ja = NULL, *CF_array;
   int     nvertices, *vlist, *rowptr, *column, ncoarse;
   int     nrows1d, nproc1d, mypidx, mypidy;
   int     A_ntotal, A_Nneigh, *A_rcvleng= NULL, *A_sndleng= NULL;
   int     *A_neigh=NULL, Nghost, **A_sndbuf=NULL, **A_rcvbuf = NULL;
   int     **A_sendlist=NULL, **A_recvlist=NULL, **proclist, *templist;
   int     *send_leng = NULL, *recv_leng = NULL;
   int     allocated = 0, *rowi_col = NULL, rowi_N;
   double  *new_val=NULL, epsilon, *dble_buf=NULL, *rowi_val=NULL, *dtemp;
   double  rowmax;
   char    *vtype, *state, *bdry;
   struct ML_CSR_MSRdata *temp;
   struct ML_CSR_MSRdata *csr_data;
   ML_GetrowFunc         *getrow_obj;
   ML_CommInfoOP         *getrow_comm;
   ML_CommInfoOP         *mat_comm;

   /* ============================================================= */
   /* get the machine information and matrix references             */
   /* ============================================================= */

   mypid        = comm->ML_mypid;
   nprocs       = comm->ML_nprocs;
   epsilon      = ml_amg->threshold;
   num_PDE_eqns = ml_amg->num_PDE_eqns;
   printflag    = ml_amg->print_flag;
   Nrows        = Amatrix->outvec_leng;

   /* ============================================================= */
   /* check the system size                                         */
   /* ============================================================= */

   if ( Nrows % num_PDE_eqns != 0 )
   {
      printf("ML_AMG_CoarsenMIS ERROR : Nrows must be multiples");
      printf(" of num_PDE_eqns.\n");
      exit(1);
   }
   if ( mypid == 0 && printflag )
   {
      printf("ML_AMG_CoarsenMIS : current level = %d\n", ml_amg->cur_level);
      printf("ML_AMG_CoarsenMIS : current eps   = %e\n", epsilon);
   }
   if ( num_PDE_eqns > 1 ) 
      ML_Operator_AmalgamateAndDropWeak(Amatrix, num_PDE_eqns, 0.0);
   Nrows /= num_PDE_eqns;
   nvertices = Amatrix->outvec_leng;

   /* ============================================================= */
   /* fetch getrow/comm information                                 */
   /* ============================================================= */

   mat_comm    = Amatrix->getrow->pre_comm;
   A_Nneigh    = ML_CommInfoOP_Get_Nneighbors(mat_comm);
   A_neigh     = ML_CommInfoOP_Get_neighbors(mat_comm);
   A_sendlist  = (int **) malloc(sizeof(int *)*A_Nneigh);
   A_recvlist  = (int **) malloc(sizeof(int *)*A_Nneigh);
   A_rcvbuf    = (int **) malloc(sizeof(int *)*A_Nneigh);
   A_sndbuf    = (int **) malloc(sizeof(int *)*A_Nneigh);
   A_rcvleng   = (int  *) malloc(sizeof(int  )*A_Nneigh);
   A_sndleng   = (int  *) malloc(sizeof(int  )*A_Nneigh);

   max_element = nvertices - 1;
   for (i = 0; i < A_Nneigh; i++) 
   {
      A_recvlist[i] = ML_CommInfoOP_Get_rcvlist(mat_comm, A_neigh[i]);
      A_rcvleng[i]  = ML_CommInfoOP_Get_Nrcvlist (mat_comm, A_neigh[i]);
      A_sendlist[i] = ML_CommInfoOP_Get_sendlist (mat_comm, A_neigh[i]);
      A_sndleng[i]  = ML_CommInfoOP_Get_Nsendlist(mat_comm, A_neigh[i]);
      A_rcvbuf[i]   = (int *) malloc(sizeof(int)*(A_rcvleng[i]+1));
      A_sndbuf[i]   = (int *) malloc(sizeof(int)*(A_sndleng[i]+1));
                           /* +1 needed inside ML_Aggregate_LabelVertices */
      for (j = 0; j < A_rcvleng[i]; j++) 
         if (A_recvlist[i][j] > max_element ) max_element = A_recvlist[i][j];
   }
   Nghost = max_element - nvertices + 1;
   A_ntotal = nvertices + Nghost;

   templist = (int *) malloc(sizeof(int)*nvertices);
   for ( i = 0; i < nvertices; i++ ) templist[i] = 0;
   for ( i = 0; i < A_Nneigh; i++ ) 
   {
      for ( j = 0; j < A_sndleng[i]; j++ ) 
      {
         index = A_sendlist[i][j];
         if ( index >= nvertices || index < 0 ) 
         {
            printf("%d : ERROR : in ML_AMG_CoarsenMIS.\n", mypid);
            exit(0);
         }
         templist[index]++;
      }
   }

   /* ============================================================= */
   /* Allocate proclist to record the processors and indices each   */
   /* of my local vertices are to send.  The first element of the   */
   /* array is a counter of how many processors, followed by a      */
   /* number of processor and index pairs.                          */
   /* ============================================================= */

   proclist = (int **) malloc(A_ntotal * sizeof( int *));
   for ( i = 0; i < nvertices; i++ ) 
   {
      if ( templist[i] > 0 )
      {
         proclist[i] = (int *) malloc( (2*templist[i]+1) * sizeof( int ) );
         proclist[i][0] = 0;
         templist[i]    = 0;
      }
      else proclist[i] = NULL;
   }
   for ( i = 0; i < A_Nneigh; i++ ) 
   {
      for ( j = 0; j < A_sndleng[i]; j++ ) 
      {
         index = A_sendlist[i][j];
         proclist[index][templist[index]+1] = i;
         proclist[index][templist[index]+2] = j;
         templist[index] += 2;
         proclist[index][0]++;
      }
   }
   for ( i = nvertices; i < A_ntotal; i++ ) 
   {
      proclist[i] = (int *) malloc( sizeof( int ) );
   }
   for ( i = 0; i < A_Nneigh; i++ ) {
      for ( j = 0; j < A_rcvleng[i]; j++ ) {
         index = A_recvlist[i][j];
         proclist[index][0] = A_neigh[i];
      }
   }
   free(templist);

   /* ============================================================= */
   /* record the Dirichlet boundary                                 */
   /* ============================================================= */

   bdry = (char *) malloc(sizeof(char)*(A_ntotal + 1));
   total_nnz = 0;
   for (i = 0; i < Nrows; i++) 
   {
      bdry[i] = 'T';
      ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
                        &rowi_N, 0);
      if (rowi_N > 1) bdry[i] = 'F';
      total_nnz += rowi_N;
   }
   rowptr = (int *) malloc( (Nrows + 1)* sizeof(int ) );
   column = (int *) malloc( total_nnz * sizeof(int ) );
   total_nnz = 0;
   for (i = 0; i < Nrows; i++) 
   {
      ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
                        &rowi_N, 0);
      rowmax = 0.0;
      for (j = 0; j < rowi_N; j++) 
         if ( dabs(rowi_val[j]) < rowmax ) rowmax = dabs( rowi_val[j] );
      rowmax *= epsilon;
      for (j = 0; j < rowi_N; j++) 
         if ( dabs(rowi_val[j]) > rowmax ) column[total_nnz++] = rowi_col[j];
      rowptr[i+1] = total_nnz;
   }

   /* ============================================================= */
   /* communicate the boundary information                          */
   /* ============================================================= */

   dtemp = (double *) ML_allocate(sizeof(double)*(A_ntotal+1));
   for (i = 0; i < nvertices; i++) 
   {
      if (bdry[i] == 'T') dtemp[i] = 1.;
      else  dtemp[i] = 0.;
   }
   ML_exchange_bdry(dtemp,Amatrix->getrow->pre_comm,nvertices,comm,ML_OVERWRITE);
   for (i = nvertices; i < A_ntotal; i++) 
   {
      if (dtemp[i] == 1.) bdry[i] = 'T';
      else bdry[i] = 'F';
   }
   ML_free(dtemp);

   CF_array = (int *) malloc(sizeof(int)* A_ntotal);
   for (i = 0; i < A_ntotal; i++) CF_array[i] = -1;

   temp   = (struct ML_CSR_MSRdata *) Amatrix->data;
   vlist  = (int *) malloc(sizeof(int)* nvertices);
   state  = (char *) malloc(sizeof(char)* A_ntotal);
   vtype  = (char *) malloc(sizeof(char)* A_ntotal);
   for (i = 0; i < nvertices   ; i++) vlist[i] = i;
   for (i = 0; i < A_ntotal; i++)     state[i] = 'F';
   for (i = 0; i < A_ntotal; i++)     vtype[i] = 'x';

   /* ============================================================= */
   /* delete nodes that are just isolated Dirichlet points          */
   /* then label the vertices as selected(S) or deleted (D)         */
   /* ============================================================= */

   for (i = 0; i < nvertices ; i++) if (bdry[i] == 'T') state[i] = 'D'; 

   k = ML_AMG_LabelVertices(nvertices, vlist, 'x', state, vtype, 
                      nvertices, rowptr, column, mypid, proclist, 
                      A_Nneigh,A_sndbuf,A_neigh, A_sndleng, A_Nneigh,
                      A_rcvbuf, A_neigh, A_rcvleng, A_recvlist, 1532, 
                      comm, CF_array);
   ncoarse = 0;
   for (i = 0; i < nvertices ; i++) if ( state[i] == 'S' ) ncoarse++;

   /* ============================================================= */
   /* check to see which fine nodes do not have a C neighbor        */
   /* ============================================================= */

   total_send_leng = 0;
   for (i = 0; i < A_Nneigh ; i++) total_send_leng += A_sndleng[i]; 
   total_recv_leng = 0;
   for (i = 0; i < A_Nneigh ; i++) total_recv_leng += A_rcvleng[i]; 
   int_buf  = (int *) malloc(total_recv_leng * sizeof(int) );
   int_buf2 = (int *) malloc(total_send_leng * sizeof(int) );

   offset = 0;
   for (i = 0; i < A_Nneigh ; i++) 
   {
      for (j = 0; j < A_sndleng[i] ; j++) 
         if ( state[A_sendlist[i][j]] == 'S' ) int_buf2[offset++] = 1; 
         else                                  int_buf2[offset++] = 0; 
   }

   ML_Aggregate_ExchangeData((char*) int_buf, (char*) int_buf2,
      A_Nneigh, A_neigh, A_rcvleng, A_sndleng, 94566, ML_INT, comm);
   for (i = 0; i < nvertices ; i++) 
   {
      if ( state[i] == 'D' )
      {
         count = 0;
         for (j = rowptr[i]; j < rowptr[i+1] ; j++) 
            if ( column[j] != i ) count += int_buf[column[j]];
         if ( count == 0 )
            printf("%4d : node %4d does not have a C neighbor\n", mypid, i);
         if ( (count+1) < (rowptr[i+1]-rowptr[i]) )
            printf("%4d : node %4d has F neighbors\n", mypid, i);
      }
   }
   free(int_buf);
   free(int_buf2);
   free(rowptr);
   free(column);

   /* ============================================================= */
   /* debug (true only for the square PDE problem)                  */
   /* ============================================================= */

   nproc1d = (int) sqrt( (double) nprocs ); 
   nrows1d = (int) sqrt( (double) nvertices );
   mypidy  = mypid / nproc1d;
   mypidx  = mypid % nproc1d;
   count   = nproc1d;
   k       = ML_gmax_int(k, comm);
   while ( count > 0 )
   {
      count--;
      m = nrows1d - 1;
      for ( i = 0; i < nrows1d; i++ )
      {
         for ( j = 0; j < nproc1d; j++ )
         {
            if ( j == mypidx && mypidy == count )
            {
               for ( k = 0; k < nrows1d; k++ )
                  printf("(%d)%c ", mypid, state[m*nrows1d+k]);
               fflush(stdout);
            }
            k = ML_gmax_int(k, comm);
         }
         if ( mypidx == nproc1d-1 && mypidy == count ) 
            {printf("\n"); fflush(stdout);}
         k = ML_gmax_int(k, comm);
         m--;
      }
   }
   k = ML_gmax_int(k, comm);

   /* ============================================================= */
   /* free memory used for doing the MIS stuff                      */
   /* ============================================================= */

   for ( i = 0; i < A_ntotal; i++ ) 
      if ( proclist[i] != NULL ) free(proclist[i]);
   free(proclist);
   free(vlist); free(state); free(vtype);
   for (i = 0; i < A_Nneigh; i++) 
   {
      free(A_recvlist[i]);
      free(A_sendlist[i]);
      free(A_rcvbuf[i]);
      free(A_sndbuf[i]);
   }
   free(A_sndleng); free(A_rcvleng);  free(A_sndbuf);
   free(A_rcvbuf);  free(A_recvlist); free(A_sendlist);
   free(A_neigh);

   /* ============================================================= */
   /* fetch communication information for A                         */
   /* ============================================================= */

   getrow_obj  = Amatrix->getrow;
   N_neighbors = getrow_obj->pre_comm->N_neighbors;
   nbytes = N_neighbors * sizeof( int );
   if ( nbytes > 0 ) 
   {
      neighbors = (int *) malloc( nbytes );
      recv_leng = (int *) malloc( nbytes );
      send_leng = (int *) malloc( nbytes );
   } 
   else neighbors = recv_leng = send_leng = NULL;

   for ( i = 0; i < N_neighbors; i++ ) 
   {
      neighbors[i] = getrow_obj->pre_comm->neighbors[i].ML_id;
      recv_leng[i] = getrow_obj->pre_comm->neighbors[i].N_rcv;
      send_leng[i] = getrow_obj->pre_comm->neighbors[i].N_send;
   }
   total_recv_leng = total_send_leng = 0;
   for ( i = 0; i < N_neighbors; i++ ) 
   {
      total_recv_leng += recv_leng[i];
      total_send_leng += send_leng[i];
   }
   nbytes = total_send_leng * sizeof( int );
   if ( nbytes > 0 ) send_list = (int *) malloc(nbytes);
   else              send_list = NULL;
   if ( total_recv_leng+Nrows != A_ntotal ) 
   {
      printf("%d : ML_AMG_CoarsenMIS - internal error.\n",mypid);
      printf("     lengths = %d %d \n",total_recv_leng+Nrows,A_ntotal);
      exit(-1);
   }
   count = 0;
   for ( i = 0; i < N_neighbors; i++ ) 
   {
      for (j = 0; j < send_leng[i]; j++)
         send_list[count++] = getrow_obj->pre_comm->neighbors[i].send_list[j];
   }
   if ( count > total_send_leng ) 
   {
      printf("%d : ML_AMG_CoarsenMIS ERROR : count < total_send_leng\n",mypid);
      exit(1);
   }

   /* ============================================================= */
   /* make sure that boundary nodes are not any C points            */
   /* ============================================================= */

   for (i = 0; i < A_ntotal; i++) if (bdry[i] == 'T') CF_array[i] = -1;

   /* ============================================================= */
   /* recover unamalgamated data                                    */
   /* ============================================================= */

   ML_Operator_UnAmalgamateAndDropWeak(Amatrix, num_PDE_eqns, 0.0);
   Nrows *= num_PDE_eqns;
   nvertices *= num_PDE_eqns;
   exp_Nrows = A_ntotal * num_PDE_eqns;
   getrow_obj   = Amatrix->getrow;
   getrow_comm  = getrow_obj->pre_comm;
   N_neighbors  = getrow_obj->pre_comm->N_neighbors;

   /* ============================================================= */
   /* Form prolongator                                              */
   /* ============================================================= */

   Ncoarse = 0;

   /* ============================================================= */
   /* update CF_array to find out which local fine point is mapped  */
   /* to which coarse point in remote processors                    */
   /* ------------------------------------------------------------- */

   /* ------------------------------------------------------------- */
   /* allocate communication buffer                                 */
   /* ------------------------------------------------------------- */

   nbytes = total_send_leng * sizeof(int);
   if ( nbytes > 0 ) int_buf = (int *) malloc( nbytes );
   else              int_buf = NULL;
   nbytes = total_recv_leng * sizeof(int);
   if ( nbytes > 0 ) int_buf2 = (int *) malloc( nbytes );
   else              int_buf2 = NULL;

   /* ------------------------------------------------------------- */
   /* send the remote node index back to remote processors, with    */
   /* added information on which remote nodes have been aggregated  */
   /* by the local aggregates (and also the aggregate numbers).     */
   /* ------------------------------------------------------------- */

   offset = 0;
   for ( i = 0; i < N_neighbors; i++ ) 
   {
      for ( j = 0; j < recv_leng[i]; j++ ) 
      {
         if ( CF_array[Nrows+offset+j] < 0 ) int_buf2[offset+j] = -1;
         else int_buf2[offset+j] = CF_array[Nrows+offset+j];
      }
      offset += recv_leng[i];
   }
   msgtype = 15963;
   ML_Aggregate_ExchangeData((char*) int_buf, (char*) int_buf2,
      N_neighbors, neighbors, send_leng, recv_leng, msgtype, ML_INT, comm);

   if ( int_buf2 != NULL ) ML_memory_free((void**) &int_buf2);


   return Ncoarse;
}

/* ******************************************************************** */
/* A subroutine to label vertices of a particular type                  */
/* -------------------------------------------------------------------- */

int ML_AMG_LabelVertices(int vlist_cnt, int *vlist, char Vtype,
                          char *vertex_state, char *vertex_type,
                          int nvertices, int *rptr, int *cptr, 
                          int myrank, int **proclist, int send_cnt, 
                          int **send_buf, int *send_proc, int *send_leng,
                          int recv_cnt, int **recv_buf, int *recv_proc, 
                          int *recv_leng, int **recv_list, int msgtype, 
                          ML_Comm *comm, int aggr_index[])
{
   int     i, j, k, m, N_remaining_vertices, index, select_flag, fproc, col;
   int     NremainingRcvProcs, change_flag, *proc_flag, send_flag,nselected;
   int     *pref_list, col2, loop_cnt, nbytes, *tlist, pref_cnt;
   int     pref_flag, pref_index;
   char    *in_preflist;
   USR_REQ *Request;
   int msg_type = 1041;

   N_remaining_vertices = vlist_cnt;
   NremainingRcvProcs   = recv_cnt;
   send_flag            = 0;
   if ( recv_cnt > 0 )
   {
      nbytes = recv_cnt * sizeof( USR_REQ );
      ML_memory_alloc((void**) &Request, nbytes, "ggp" );
      nbytes = recv_cnt * sizeof( int );
      ML_memory_alloc((void**) &proc_flag, nbytes, "ggq" );
      for ( i = 0; i < recv_cnt; i++ ) proc_flag[i] = 0;
   }
   for ( j = 0; j < send_cnt; j++ )
      for ( k = 0; k <= send_leng[j]; k++ ) send_buf[j][k] = 0;

   /* First clear out any already deleted vertices (for example corresponding */
   /* to Dirichlet BC points. */

   for ( i = 0; i < vlist_cnt; i++ ) 
   {
      index = vlist[i];
      if (vertex_state[index] == 'D') 
      {
         N_remaining_vertices--;
         if ( proclist[index] != NULL )
         {
            for ( k = 0; k < proclist[index][0]; k++ ) 
            {
               fproc = proclist[index][2*k+1];
               m     = proclist[index][2*k+2];
               send_buf[fproc][m] = 2;
            }
         }
      }
   }

   ML_Aggregate_UpdateVertexStates(N_remaining_vertices, vertex_state,
	recv_cnt, recv_proc, recv_leng, recv_buf,
	recv_list, proc_flag, &NremainingRcvProcs,
	send_cnt, send_proc, send_leng, send_buf,
	&send_flag, Request, comm, msg_type);

   /* ---------------------------------------------------------- */
   /* give the vertices adjacent to deleted vertices preferences */
   /* ---------------------------------------------------------- */

   in_preflist = (char *) malloc(vlist_cnt*sizeof(char) );
   for (i = 0; i < vlist_cnt; i++) in_preflist[i] = 'f';
   if ( vlist_cnt > 0 )
   {
      nbytes = vlist_cnt * sizeof( int );
      ML_memory_alloc((void**) &tlist, nbytes, "ggn" );
      for ( i = 0; i < vlist_cnt; i++ ) tlist[i] = vlist[i];
      for ( i = 0; i < vlist_cnt; i++ )
      {
         index = tlist[i];
         for ( j = rptr[index]; j < rptr[index+1]; j++ )
         {
            col = cptr[j];
            if ( vertex_state[col] == 'D' )
            {
               tlist[i] = - index;
               break;
            }
         }
      }
      m = 0;
      for ( i = 0; i < vlist_cnt; i++ )
      {
         if ( tlist[i] < 0 ) vlist[m++] = - tlist[i];
      }
      for ( i = 0; i < vlist_cnt; i++ )
      {
         if ( tlist[i] >= 0 ) vlist[m++] = tlist[i];
      }
      ML_memory_free( (void**) &tlist );
   }
   if ( nvertices > 0 )
   {
      nbytes = nvertices * sizeof( int );
      ML_memory_alloc((void**) &pref_list, nbytes, "ggo" );
   }   
   pref_cnt = 0;
   
   /* -------------------------------------------------------- */
   /* get ready for the coarsening                             */
   /* -------------------------------------------------------- */

   nselected = 0;

   /* -------------------------------------------------------- */
   /* let's actually do coarsening                             */
   /* -------------------------------------------------------- */

   change_flag = 1;
   loop_cnt = 0;
   pref_index = 0;     /* pointer to a stack of vertex numbers */

   do {
      pref_index = 0;     /* pointer to a stack of vertex numbers */
      /* loop_cnt is to monitor the performance of coarsening */

      loop_cnt++;

      /* reset all buffers to zero only if it has been changed */

      if ( change_flag == 1 )
      {
         for ( j = 0; j < send_cnt; j++ )
            for ( k = 0; k <= send_leng[j]; k++ ) send_buf[j][k] = 0;
         change_flag = 0;
      }

      /* examine the vertices in vlist */

      for ( i = 0; i < vlist_cnt; i++ )
      {

         /* handle the preference list first, if there is any */
         /* Note : we want to fetch the pref_list from the    */
         /*        front                                      */

         index = vlist[i];
         pref_flag = 0;
         if ( pref_cnt > pref_index ) {
            do {
               index = pref_list[pref_index];    
               in_preflist[index] = 'f';
               for (j = pref_index+1; j < pref_cnt; j++) 
                  pref_list[j-1] = pref_list[j];
               pref_cnt--;
            } while ( (vertex_state[index] != 'F') && (pref_cnt > pref_index));
            pref_flag = 1;
            i--;
         }

         /* if the vertex in question has not been considered F(ree) */

         if ( vertex_state[index] == 'F' )
         {
            select_flag = 1;
            for ( j = rptr[index]; j < rptr[index+1]; j++ )
            {
               /* if its neighbor is selected, delete this vertex */

               col = cptr[j];
               if ( vertex_state[col] == 'S' )
               {
                  if (vertex_state[index] == 'F') {
                     vertex_state[index] = 'D';
                     N_remaining_vertices--;
                     if ( proclist[index] != NULL )
                     {
                        for ( k = 0; k < proclist[index][0]; k++ ) 
                        {
                           fproc = proclist[index][2*k+1];
                           m     = proclist[index][2*k+2];
                           send_buf[fproc][m] = 2;
                           change_flag = 1;
                        }
                     }

                     /* also, put the next set of vertices into the  */
                     /* preference list (try to mimic the sequential */
                     /* maximally independent set algorithm          */

                        for ( k = rptr[index]; k < rptr[index+1]; k++ ) {
                           col2 = cptr[k];
                           if (col2 < nvertices && 
                               vertex_state[col2] == 'F' &&
                               vertex_type[col2] == Vtype ) {
                              if (in_preflist[col2] != 't') {
                                 pref_list[pref_cnt++] = col2;
                                 in_preflist[col2] = 't';
                              }
                           }
                        }
                  }
                  select_flag = 0;
                  break;
               }
               
               /* If its neighbor is of the same type and not been   */
               /* considered. Furthermore, if it is a remote vertex  */
               /* and its owner processor has rank smaller than mine,*/
               /* my processor should wait(thus turn off select_flag)*/

               else if ( vertex_type[col] == Vtype && 
                         vertex_state[col] == 'F')
               {
                  if ( col >= nvertices )
                  {
                     if ( proclist[col][0] < myrank )
                     {
                        select_flag = 0;
                        break;
                     }
                  }
               }
            }

            /* if the vertex in question is not any of those considered */
            /* above, select this vertex.                               */

            if ( select_flag == 1 )
            {
               if ((vertex_state[index] == 'F') &&(index < nvertices)) 
                  N_remaining_vertices--;
               vertex_state[index] = 'S';
               aggr_index[index] = nselected;
               nselected++;

               /* set the flag that this vertex has been selected in */
               /* the buffer which is to be sent to other processors */

               if ( proclist[index] != NULL )
               {
                  for ( k = 0; k < proclist[index][0]; k++ )
                  {
                     fproc = proclist[index][2*k+1];
                     m     = proclist[index][2*k+2];
                     send_buf[fproc][m] = 1;
                     change_flag = 1;
                  }
               }

               /* delete all vertices adjacent to this vertex and */
               /* indicate that also in the communication buffer  */

               for ( j = rptr[index]; j < rptr[index+1]; j++ )
               {
                  col = cptr[j];
                  if (vertex_state[col] == 'F') {
                     vertex_state[col] = 'D';
                     if ( col < nvertices ) {
                        N_remaining_vertices--;
                        if ( proclist[col] != NULL )
                        {
                           for ( k = 0; k < proclist[col][0]; k++ ) 
                           {
                              fproc = proclist[col][2*k+1];
                              m     = proclist[col][2*k+2];
                              send_buf[fproc][m] = 2;
                              change_flag = 1;
                           }
                        }

                        /* also, put the next set of vertices into the  */
                        /* preference list (try to mimic the sequential */
                        /* maximally independent set algorithm          */

                        for ( k = rptr[col]; k < rptr[col+1]; k++ ) {
                           col2 = cptr[k];
                           if (col2 < nvertices && 
                               vertex_state[col2] == 'F' &&
                               vertex_type[col2] == Vtype ) {
                              if (in_preflist[col2] != 't') {
                                 pref_list[pref_cnt++] = col2;
                                 in_preflist[col2] = 't';
                              }
                           }
                        }
                     }
                  }
               } 
            }
         }

         /* if after the steps above, the vertex is still not */
         /* selected.  Well, do something about it.           */

         if ( vertex_state[index] == 'F' )
         {
            /* if a vertex in the pref_list has been considered */
            /* but not selected, need to put the vertex back to */
            /* the list, and move on to consider the next one   */
            /* (i.e. advance pref_index)                        */

            if ( pref_flag == 1 )
            {
               for (j = pref_cnt-1; j >= pref_index; j--) 
                  pref_list[j+1] = pref_list[j];
               pref_list[pref_index] = index;
               in_preflist[index] = 't';
               pref_index++; pref_cnt++;
            }
         }
      }
      ML_Aggregate_UpdateVertexStates(N_remaining_vertices, vertex_state,
	recv_cnt, recv_proc, recv_leng, recv_buf,
	recv_list, proc_flag, &NremainingRcvProcs,
	send_cnt, send_proc, send_leng, send_buf,
	&send_flag, Request, comm, msg_type);

#ifdef ML_AGGR_out
      /* update the states to/from other processors */

      msgtype += 131;
      for ( j = 0; j < recv_cnt; j++ )
      {
         if ( proc_flag[j] == 0 )
         {
            fproc = recv_proc[j];
            nbytes = (recv_leng[j] + 1) * sizeof( int );
            comm->USR_irecvbytes((char*) recv_buf[j], nbytes, &fproc,
                    &msgtype, comm->USR_comm, (void *) &Request[j] );
         }
      }
      if ( send_flag == 0 ) {
      for ( j = 0; j < send_cnt; j++ )
      {
         nbytes = (send_leng[j] + 1) * sizeof( int );
         if ( N_remaining_vertices <= 0 ) { 
            send_buf[j][send_leng[j]] = 1; send_flag = 1; 
         }
         comm->USR_sendbytes((void*) send_buf[j], nbytes,
                             send_proc[j], msgtype, comm->USR_comm );
      }
      }
      for ( j = 0; j < recv_cnt; j++ )
      {
         if ( proc_flag[j] == 0 )
         {
            fproc = recv_proc[j];
            nbytes = (recv_leng[j] + 1) * sizeof( int );
            comm->USR_waitbytes((char*) recv_buf[j], nbytes, &fproc,
                     &msgtype, comm->USR_comm, (void *) &Request[j] );
            for ( k = 0; k < recv_leng[j]; k++ )
            {
               kkk = recv_list[j][k];
               if      (recv_buf[j][k] == 1)
               {
                  vertex_state[kkk] = 'S';
               }
               else if (recv_buf[j][k] == 2)
               {
                  vertex_state[kkk] = 'D';
               }
            }
            if ( recv_buf[j][recv_leng[j]] == 1 )
            {
               proc_flag[j] = 1;
               NremainingRcvProcs--;
            }
         }
      }
#endif
   } while ( NremainingRcvProcs > 0 || N_remaining_vertices > 0 );

/* #ifdef ML_AGGR_DEBUG */
   printf("%d : loop_count = %d \n", myrank, loop_cnt ); fflush(stdout);
/* #endif */
   if ( recv_cnt > 0 )
   {
      ML_memory_free( (void **) &proc_flag );
      ML_memory_free( (void **) &Request );
   }
   if ( nvertices > 0 ) ML_memory_free( (void **) &pref_list );
   free(in_preflist);
   return nselected;
}
