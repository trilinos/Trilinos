/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* ******************************************************************** */
/* Functions to create tentative prolongators                           */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Ray Tuminaro (SNL)           */
/* Date          : May, 2000                                            */
/* ******************************************************************** */
/* ******************************************************************** */

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ml_aggregate.h"
#include "ml_lapack.h"
#include "ml_agg_genP.h"

/* ******************************************************************** */
/* variables used for parallel debugging  (Ray)                         */
/* -------------------------------------------------------------------- */

#ifdef ML_AGGR_PARTEST
extern int **global_mapping = NULL, global_nrows, global_ncoarse;
#endif

/* ******************************************************************** */
/* local defines                                                        */
/* -------------------------------------------------------------------- */

#define AGGR_READY      -11
#define AGGR_NOTSEL     -12
#define AGGR_SELECTED   -13
#define AGGR_SELECTED2  -14
#define AGGR_READY2     -15
#define AGGR_BDRY       -16
#define ML_AGGR_MINRANK   1
#define ML_AGGR_MAXLINK   2
#define ML_AGGR_UNCOUPLED 1
#define ML_AGGR_COUPLED   2

#if defined(OUTPUT_AGGREGATES) || defined(INPUT_AGGREGATES) || (ML_AGGR_INAGGR) || (ML_AGGR_OUTAGGR) || (ML_AGGR_MARKINAGGR)
#ifndef MAXWELL
#ifndef ALEGRA
extern int *update_index, *update, *extern_index, *external;
#endif
#else
extern int *reordered_glob_nodes, *global_node_inds,
           *reordered_node_externs, *global_node_externs;
#endif /*ifdef MAXWELL */
extern int ML_gpartialsum_int(int val, ML_Comm *comm);
#endif


/* ******************************************************************** */
/* ******************************************************************** */
/* The following subroutine is useful to ensure coupling between        */
/* processors.  However, it has not been updated as good as the         */
/* ML_Aggregate_CoarsenUncoupled subroutine.                            */
/*               --- Charles Tong, December 29, 1999                    */
/* ******************************************************************** */
/* ******************************************************************** */
/* ******************************************************************** */
/* construct the tentative prolongator allowing aggregate to cross      */
/* processor boundaries                                                 */
/* -------------------------------------------------------------------- */

int ML_Aggregate_CoarsenMIS( ML_Aggregate *ml_ag, ML_Operator *Amatrix, 
                                 ML_Operator **Pmatrix, ML_Comm *comm)
{
  unsigned int nbytes, length;
   int     i, j, jj, k, m, Nrows, exp_Nrows;
   int     N_neighbors, *neighbors = NULL, diff_level;
   double  printflag;
   int     *Asqrd_rcvleng= NULL, *Asqrd_sndleng= NULL, *send_list = NULL;
   int     total_recv_leng = 0, total_send_leng = 0, offset, msgtype;
   int     aggr_count, index, mypid, num_PDE_eqns;
   int     *aggr_index = NULL, nullspace_dim;
   int     *dof_accounted_for = NULL;
   int     *sendlist_proc = NULL, Ncoarse, count, *int_buf = NULL;
   int     *int_buf2, *aggr_stat = NULL;
   int     index4, *new_send_leng = NULL, new_N_send;
   int     *new_send_neighbors = NULL, *new_send_list = NULL;
   int     max_count, *new_ia = NULL, *new_ja = NULL, new_Nrows;
   int     *new_recv_leng=NULL, exp_Ncoarse, new_N_recv;
   int     *new_recv_neighbors=NULL, *aggr_cnt_array = NULL;
   int     level, index3, count3, *recv_list = NULL, max_agg_size;
   int     **rows_in_aggs = NULL, lwork, info;
   double  dcompare1, *new_val = NULL, epsilon;
   double  *dble_buf = NULL, *nullspace_vect = NULL, *qr_tmp = NULL;
   double  *tmp_vect = NULL, *work = NULL, *new_null = NULL, *comm_val = NULL;
   double  *dble_buf2;
   ML_SuperNode          *aggr_head = NULL, *aggr_curr, *supernode;
   struct ML_CSR_MSRdata *csr_data;
   ML_Aggregate_Comm     *aggr_comm;
   ML_GetrowFunc         *getrow_obj;
   ML_Operator           *Asqrd = NULL, *tmatrix;
   struct ML_CSR_MSRdata *temp;
   char                  *vtype, *state, *bdry, *unamalg_bdry;
   int                   nvertices, *vlist, Asqrd_ntotal, Asqrd_Nneigh;
   int                   *Asqrd_neigh = NULL, max_element, Nghost;
   int                   **Asqrd_sndbuf = NULL, **Asqrd_rcvbuf = NULL;
   int                   **Asqrd_sendlist = NULL, **Asqrd_recvlist = NULL;
   ML_CommInfoOP         *mat_comm;
   int                   allocated = 0, *rowi_col = NULL, rowi_N, current;
   double                *rowi_val = NULL, *dtemp;
   int                   *templist, **proclist, *temp_index;
   int                   *temp_leng, *tem2_index, *tempaggr_index = NULL;
   int                   *send_leng = NULL, *recv_leng = NULL;
   int                   total_nz = 0;
   int                   count2;
   ML_agg_indx_comm      agg_indx_comm;

   /*   int kk, old_upper, nnzs, count2, newptr; */

#ifdef CLEAN_DEBUG
   char tlabel[80];
#endif

#ifdef DDEBUG
   int curagg,myagg,*good,*bad, kk;
#endif

#if defined(OUTPUT_AGGREGATES) || defined(INPUT_AGGREGATES) || (ML_AGGR_INAGGR) || (ML_AGGR_OUTAGGR) || (ML_AGGR_MARKINAGGR)
FILE *fp;
char fname[80];
static int level_count = 0;
double *d2temp;
int agg_offset, vertex_offset;
#endif

   /* ============================================================= */
   /* get the machine information and matrix references             */
   /* ============================================================= */

   mypid                   = comm->ML_mypid;
   epsilon                 = ml_ag->threshold;
   num_PDE_eqns            = ml_ag->num_PDE_eqns;
   nullspace_dim           = ml_ag->nullspace_dim;
   nullspace_vect          = ml_ag->nullspace_vect;
   Nrows                   = Amatrix->outvec_leng;
   printflag               = ml_ag->print_flag;

   /* ============================================================= */
   /* check the system size versus null dimension size              */
   /* ============================================================= */

   if ( Nrows % num_PDE_eqns != 0 )
   {
      printf("ML_Aggregate_CoarsenMIS ERROR : Nrows must be multiples");
      printf(" of num_PDE_eqns.\n");
      exit(1);
   }

    /*= the following two lines are for the solution of elasticity
     problems on fine grid */

   diff_level = ml_ag->max_levels - ml_ag->cur_level - 1;
   if ( diff_level > 0 ) num_PDE_eqns = nullspace_dim; /* ## 12/20/99 */


   /* ============================================================= */
   /* Figure out where the Dirichlet points are on the fine grid of */ 
   /* the unamalgmated system.                                      */ 
   /* ============================================================= */

   Nghost = Amatrix->getrow->pre_comm->total_rcv_length;
   unamalg_bdry = (char *) ML_allocate(sizeof(char)*(Nrows+Nghost+1));
   dtemp = (double *) ML_allocate(sizeof(double)*(Nrows+Nghost+1));
   if (dtemp == NULL) pr_error("ml_agg_MIS: out of space.\n");

   for (i = 0; i < Nrows+Nghost; i++) dtemp[i] = 0.;

   for (i = 0; i < Nrows; i++) {
      ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
                        &rowi_N, 0);
      count2 = 0;
      for (j = 0; j < rowi_N; j++) if (rowi_val[j] != 0.) count2++;
      if (count2 <= 1) dtemp[i] = 1.;
   }
   ML_free(rowi_col); ML_free(rowi_val);
   rowi_col = NULL; rowi_val = NULL;
   allocated = 0; 

   ML_exchange_bdry(dtemp,Amatrix->getrow->pre_comm,Amatrix->outvec_leng,
                    comm, ML_OVERWRITE,NULL);
   for (i = 0; i < Nrows+Nghost; i++) {
      if (dtemp[i] == 1.) unamalg_bdry[i] = 'T';
      else unamalg_bdry[i] = 'F';
   }
   ML_free(dtemp);

   /* ============================================================= */
   /* set up the threshold for weight-based coarsening              */
   /* ============================================================= */

   /*= the following four lines define some stuff for dropping */

   diff_level = ml_ag->begin_level - ml_ag->cur_level;
   if ( diff_level == 0 ) ml_ag->curr_threshold = ml_ag->threshold;
   epsilon = ml_ag->curr_threshold;
   ml_ag->curr_threshold *= 0.5;

   if ( mypid == 0 && printflag < ML_Get_PrintLevel())
   {
      printf("ML_Aggregate_CoarsenMIS : current level = %d\n",
                                            ml_ag->cur_level);
      printf("ML_Aggregate_CoarsenMIS : current eps = %e\n",epsilon);
   }
/*
   epsilon = epsilon * epsilon;
*/

#ifdef CLEAN_DEBUG
   sprintf(tlabel,"before amalg %d",comm->ML_mypid);
   ML_CommInfoOP_Print(Amatrix->getrow->pre_comm, tlabel);
#endif
   /*= `amalgamate' refers to the conversion from block matrices
     to point matrices */
   ML_Operator_AmalgamateAndDropWeak(Amatrix, num_PDE_eqns, epsilon);
   Nrows /= num_PDE_eqns;
#ifdef CLEAN_DEBUG
   sprintf(tlabel,"after amalg %d",comm->ML_mypid);
   ML_CommInfoOP_Print(Amatrix->getrow->pre_comm, tlabel);
#endif

   nvertices = Amatrix->outvec_leng;

   /* ============================================================= */
   /* square the matrix and fetch getrow/comm information           */
   /* ============================================================= */

   Asqrd    = ML_Operator_Create(comm);
   tmatrix  = ML_Operator_halfClone(Amatrix);

   /*= construction of the square of matrix A */

   ML_2matmult(tmatrix, Amatrix, Asqrd, ML_CSR_MATRIX );
   ML_Operator_halfClone_Destroy(tmatrix);

   mat_comm        = Asqrd->getrow->pre_comm;
   Asqrd_Nneigh    = ML_CommInfoOP_Get_Nneighbors(mat_comm);
   Asqrd_neigh     = ML_CommInfoOP_Get_neighbors(mat_comm);
   Asqrd_sendlist  = (int **) ML_allocate(sizeof(int *)*Asqrd_Nneigh);
   Asqrd_recvlist  = (int **) ML_allocate(sizeof(int *)*Asqrd_Nneigh);
   Asqrd_rcvbuf    = (int **) ML_allocate(sizeof(int *)*Asqrd_Nneigh);
   Asqrd_sndbuf    = (int **) ML_allocate(sizeof(int *)*Asqrd_Nneigh);
   Asqrd_rcvleng   = (int  *) ML_allocate(sizeof(int  )*Asqrd_Nneigh);
   Asqrd_sndleng   = (int  *) ML_allocate(sizeof(int  )*Asqrd_Nneigh);

   max_element = nvertices - 1;
   for (i = 0; i < Asqrd_Nneigh; i++) {
      Asqrd_recvlist[i]  = ML_CommInfoOP_Get_rcvlist(mat_comm, Asqrd_neigh[i]);
      Asqrd_rcvleng[i]   = ML_CommInfoOP_Get_Nrcvlist (mat_comm, Asqrd_neigh[i]);
      Asqrd_sendlist[i]  = ML_CommInfoOP_Get_sendlist (mat_comm, Asqrd_neigh[i]);
      Asqrd_sndleng[i]  = ML_CommInfoOP_Get_Nsendlist(mat_comm, Asqrd_neigh[i]);
      Asqrd_rcvbuf[i]    = (int *) ML_allocate(sizeof(int)*(Asqrd_rcvleng[i]+1));
      Asqrd_sndbuf[i]    = (int *) ML_allocate(sizeof(int)*(Asqrd_sndleng[i]+1));
                           /* +1 needed inside ML_Aggregate_LabelVertices */
      for (j = 0; j < Asqrd_rcvleng[i]; j++) 
         if (Asqrd_recvlist[i][j] > max_element ) 
            max_element = Asqrd_recvlist[i][j];
   }
   Nghost = max_element - nvertices + 1;
   Asqrd_ntotal = nvertices + Nghost;

   templist = (int *) ML_allocate(sizeof(int)*nvertices);
   for ( i = 0; i < nvertices; i++ ) templist[i] = 0;
   for ( i = 0; i < Asqrd_Nneigh; i++ ) {
      for ( j = 0; j < Asqrd_sndleng[i]; j++ ) {
         index = Asqrd_sendlist[i][j];
         if ( index >= nvertices || index < 0 ) {
            printf("%d : Error : in Coarsening.\n", mypid);
            exit(0);
         }
         templist[index]++;
	 /*= templist[j] is the number of processors who need `j' */	 
      }
   }

   /* Allocate proclist to record the processors and indices each of */
   /* my local vertices are to send.  The first element of the array */
   /* is a counter of how many processors, followed by a number of   */
   /* processor and index pairs.                                     */

   /*= other stuff needed by MIS */
   proclist = (int **) ML_allocate(Asqrd_ntotal * sizeof( int *));
   for ( i = 0; i < nvertices; i++ ) {
      proclist[i]    = (int *) ML_allocate( (2*templist[i]+1) * sizeof( int ) );
      proclist[i][0] = 0;
      templist[i]    = 0;
   }
   for ( i = 0; i < Asqrd_Nneigh; i++ ) {
      for ( j = 0; j < Asqrd_sndleng[i]; j++ ) {
         index = Asqrd_sendlist[i][j];
         proclist[index][templist[index]+1] = i;
         proclist[index][templist[index]+2] = j;
         templist[index] += 2;
         proclist[index][0]++;
      }
   }
   for ( i = nvertices; i < Asqrd_ntotal; i++ ) {
      proclist[i] = (int *) ML_allocate( sizeof( int ) );
   }
   for ( i = 0; i < Asqrd_Nneigh; i++ ) {
      for ( j = 0; j < Asqrd_rcvleng[i]; j++ ) {
         index = Asqrd_recvlist[i][j];
         proclist[index][0] = Asqrd_neigh[i];
      }
   }
   ML_free(templist);

   /* Need to set up communication for the matrix A (as opposed to Asqrd) */

   getrow_obj   = Amatrix->getrow;
   N_neighbors  = getrow_obj->pre_comm->N_neighbors;
   total_recv_leng = 0;
   for (i = 0; i < N_neighbors; i++) {
      total_recv_leng += getrow_obj->pre_comm->neighbors[i].N_rcv;
   }
   recv_list   = (int *) ML_allocate(sizeof(int *)*total_recv_leng);
   max_element = nvertices - 1;
   count = 0;
   for (i = 0; i < N_neighbors; i++) {
      for (j = 0; j < getrow_obj->pre_comm->neighbors[i].N_rcv; j++) {
         recv_list[count] = getrow_obj->pre_comm->neighbors[i].rcv_list[j];
         if (recv_list[count] > max_element ) max_element = recv_list[count];
         count++;
      }
   }
   Nghost = max_element - nvertices + 1;
   exp_Nrows = nvertices + Nghost;

   /* record the Dirichlet boundary */
   /*= amalgamate boundary */

   bdry = (char *) ML_allocate(sizeof(char)*(exp_Nrows + 1));
   for (i = Nrows ; i < exp_Nrows; i++) bdry[i] = 'F';
   for (i = 0; i < Nrows; i++) {
      bdry[i] = 'T';
      ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
                        &rowi_N, 0);
      if (rowi_N > 1) bdry[i] = 'F';
   }

   /* communicate the boundary information */

   dtemp = (double *) ML_allocate(sizeof(double)*(exp_Nrows+1));
   for (i = nvertices; i < exp_Nrows; i++) dtemp[i] = 0;
   for (i = 0; i < nvertices; i++) {
      if (bdry[i] == 'T') dtemp[i] = 1.;
      else  dtemp[i] = 0.;
   }
   ML_exchange_bdry(dtemp,Amatrix->getrow->pre_comm,nvertices,comm,
                    ML_OVERWRITE,NULL);
   for (i = nvertices; i < exp_Nrows; i++) {
      if (dtemp[i] == 1.) bdry[i] = 'T';
      else bdry[i] = 'F';
   }
   ML_free(dtemp);


   aggr_index = (int *) ML_allocate(sizeof(int)* Asqrd_ntotal*num_PDE_eqns);
   for (i = 0; i < Asqrd_ntotal; i++) aggr_index[i] = -1;

   temp   = (struct ML_CSR_MSRdata *) Asqrd->data;
   vlist  = (int *) ML_allocate(sizeof(int)* nvertices);
   state  = (char *) ML_allocate(sizeof(char)* Asqrd_ntotal);
   vtype  = (char *) ML_allocate(sizeof(char)* Asqrd_ntotal);
   for (i = 0; i < nvertices   ; i++)  vlist[i] =   i;
   for (i = 0; i < Asqrd_ntotal; i++)  state[i] = 'F';
   for (i = 0; i < Asqrd_ntotal; i++)  vtype[i] = 'x';

   /* delete nodes that are just isolated Dirichlet points */

   for (i = 0; i < nvertices ; i++) {
     if (bdry[i] == 'T') state[i] = 'D'; 
                               /* remember communication for A and Asqrd */
                               /* not the same.  We will communicate     */
                               /* this stuff inside ML_Aggregate_LabelVertices*/
   }

   aggr_count = ML_Aggregate_LabelVertices(nvertices, vlist,'x',state, vtype, 
                      nvertices, temp->rowptr, temp->columns, mypid, proclist, 
                      Asqrd_Nneigh,Asqrd_sndbuf,Asqrd_neigh, Asqrd_sndleng,
                      Asqrd_Nneigh,Asqrd_rcvbuf, Asqrd_neigh, Asqrd_rcvleng,
                      Asqrd_recvlist, 1532, comm, aggr_index);

   /* free memory used for doing the MIS stuff */

   for ( i = 0; i < Asqrd_ntotal; i++ ) ML_free(proclist[i]);
   ML_free(proclist);
   ML_free(vlist); ML_free(state); ML_free(vtype);
   for (i = 0; i < Asqrd_Nneigh; i++) {
      ML_free(Asqrd_recvlist[i]);
      ML_free(Asqrd_sendlist[i]);
      ML_free(Asqrd_rcvbuf[i]);
      ML_free(Asqrd_sndbuf[i]);
   }
   ML_free(Asqrd_sndleng); ML_free(Asqrd_rcvleng);  ML_free(Asqrd_sndbuf);
   ML_free(Asqrd_rcvbuf);  ML_free(Asqrd_recvlist); ML_free(Asqrd_sendlist);
   ML_free(Asqrd_neigh);
   ML_Operator_Destroy(&Asqrd);


   /* Make a new aggr_index[] corresponding to Amatrix as opposed to the one  */
   /* that already exists for Asqrd. These are not entirely the same because  */
   /* the two have different receive lists and in fact Amatrix could have a   */
   /* longer receive list than Asqrd (due to dropping). The good news is that */
   /* at this point we only need information in regular variables (not ghost) */

   tempaggr_index = (int *) ML_allocate(sizeof(int)* (exp_Nrows+1)*num_PDE_eqns);
   for (i = 0;         i < nvertices ; i++) tempaggr_index[i] = aggr_index[i];
   for (i = nvertices; i < exp_Nrows ; i++) tempaggr_index[i] = -1;
   ML_free(aggr_index);
   aggr_index = tempaggr_index;

   N_neighbors = getrow_obj->pre_comm->N_neighbors;
   nbytes = N_neighbors * sizeof( int );
   if ( nbytes > 0 ) {
         ML_memory_alloc((void**) &neighbors,  nbytes, "AGL");
         ML_memory_alloc((void**) &recv_leng,  nbytes, "AGM");
         ML_memory_alloc((void**) &send_leng,  nbytes, "AGN");
   } 
   else {
         neighbors = recv_leng = send_leng = NULL;
   }
   for ( i = 0; i < N_neighbors; i++ ) {
      neighbors[i] = getrow_obj->pre_comm->neighbors[i].ML_id;
      recv_leng[i] = getrow_obj->pre_comm->neighbors[i].N_rcv;
      send_leng[i] = getrow_obj->pre_comm->neighbors[i].N_send;
   }
   total_recv_leng = total_send_leng = 0;
   for ( i = 0; i < N_neighbors; i++ ) {
         total_recv_leng += recv_leng[i];
         total_send_leng += send_leng[i];
   }
   nbytes = total_send_leng * sizeof( int );
   if ( nbytes > 0 ) ML_memory_alloc((void**) &send_list,nbytes,"AGO");
   else              send_list = NULL;
   if ( total_recv_leng+Nrows != exp_Nrows ) {
         printf("%d : ML_Aggregate_CoarsenMIS - internal error.\n",mypid);
         printf("     lengths = %d %d \n",total_recv_leng+Nrows,exp_Nrows);
         exit(-1);
   }
   count = 0;
   for ( i = 0; i < N_neighbors; i++ ) {
     /*= build one big send list for A */
         for (j = 0; j < send_leng[i]; j++)
            send_list[count++] = 
               getrow_obj->pre_comm->neighbors[i].send_list[j];
   }
   if ( count > total_send_leng ) {
         printf("%d : CoarsenMIS ERROR : count < total_send_leng\n",mypid);
         exit(1);
   }

   /* take MIS points and make phase1 style aggregates (i.e. just take */
   /* the neighbors of each MIS point as an aggregate).                */

   total_nz = 0;
   count    = 0;
   for (i = 0; i < nvertices; i++) {
      if (aggr_index[i] >= 0) {
	/*= encode aggregate numbers */
         current = -aggr_index[i]-2;
         aggr_index[i] = current;
         ML_get_matrix_row(Amatrix,1,&i,&allocated,&rowi_col,&rowi_val,&rowi_N,0);
         for (j = 0; j < rowi_N; j++) aggr_index[rowi_col[j]] = current;
      }
      else {
         /* still get the rows to count nonzeros */
         ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
			   &rowi_N, 0);
      }
      total_nz += rowi_N;
   }
   total_nz = ML_Comm_GsumInt( comm, total_nz);
   i = ML_Comm_GsumInt( comm, nvertices);

   if ( mypid == 0 && printflag < ML_Get_PrintLevel())
      printf("Aggregation(MIS) : Total nonzeros = %d (Nrows=%d)\n",total_nz,i);

   if ( ml_ag->operator_complexity == 0.0 ) {
      ml_ag->fine_complexity = total_nz;
      ml_ag->operator_complexity = total_nz;
   }
   else ml_ag->operator_complexity += total_nz;

   /*= decode aggregate numbers */
   for (i = 0; i < exp_Nrows; i++) 
      aggr_index[i] = -aggr_index[i]-2;

   /* make sure that boundary nodes are not in any aggregate */

   for (i = 0; i < exp_Nrows; i++)
      if (bdry[i] == 'T') aggr_index[i] = -1;

   /* communicate the aggregate information.  Use temp_index and tem2_index as */
   /* communication buffer arrays.                                             */
   temp_index = (int *) ML_allocate(sizeof(int)*(1+total_send_leng+total_recv_leng));
   tem2_index = (int *) ML_allocate(sizeof(int)*(1+total_send_leng+total_recv_leng));
   temp_leng  = (int *) ML_allocate(sizeof(int)*(N_neighbors+1));

    /*= update aggregate ghost nodes */
   ML_aggr_index_communicate(N_neighbors, temp_leng, send_leng, recv_leng, send_list, 
			     recv_list, nvertices, comm, aggr_index, 1563, tem2_index, 
			     neighbors, temp_index,1);

   agg_indx_comm.N_neighbors = N_neighbors;
   agg_indx_comm.temp_leng = temp_leng;
   agg_indx_comm.send_leng = send_leng;
   agg_indx_comm.recv_leng = recv_leng;
   agg_indx_comm.send_list = send_list;
   agg_indx_comm.recv_list = recv_list;
   agg_indx_comm.tem2_index = tem2_index;
   agg_indx_comm.neighbors = neighbors;
   agg_indx_comm.temp_index = temp_index;

   ML_Aggregate_Phase2_3_Cleanup(ml_ag, Amatrix, &aggr_count, nvertices, aggr_index,
				 exp_Nrows, comm, bdry, "MIS", &agg_indx_comm);

   /* make sure that boundary nodes are not in any aggregate */
   /* Also, make sure that everyone else is aggregated.      */

   for (i = 0; i < exp_Nrows; i++) {
     if (bdry[i] == 'T') aggr_index[i] = -1;
   }


   /* communicate the phase two/three information  */

   if ( mypid == 0 && 8 < ML_Get_PrintLevel())
   {
      printf("Communicating phase 2/3 info\n");
      fflush(stdout);
   }
   ML_aggr_index_communicate(N_neighbors, temp_leng, send_leng, recv_leng, send_list, 
			     recv_list, nvertices, comm, aggr_index, 1963, tem2_index, 
			     neighbors, temp_index,1);

   ML_free(temp_leng);
   ML_free(tem2_index);
   ML_free(temp_index);
   ML_free(rowi_col); rowi_col = NULL;
   ML_free(rowi_val); rowi_val = NULL;
   allocated = 0;
   ML_free(recv_list);

#ifdef ML_AGGR_MARKINAGGR

   for (i = 0; i < exp_Nrows; i++) aggr_index[i] = -1;
   sprintf(fname,"agg_%d",level_count); level_count++;
   fp = fopen(fname,"r");
   aggr_count = 0;
   for (i = 0; i <nvertices; i++) {
      fscanf(fp,"%d%d",&k,&j);
      aggr_index[j] = k;
      if (k >= aggr_count) aggr_count = k+1;
   }
   fclose(fp);
#endif /*ifdef ML_AGGR_MARKINAGGR*/

   /* make sure that boundary nodes are not in any aggregate */

   for (i = 0; i < exp_Nrows; i++) {
     if (bdry[i] == 'T') aggr_index[i] = -1;
     else if (aggr_index[i] == -1)
	 {
	   printf("for node %d, bndry[.] = %c\n",i,bdry[i]);
       printf("ML_agg_MIS: I'm not sure who takes care of this guy\n");
	   printf("probably need to use the other aggregate file...\n");
     }
   }

   ML_free(bdry);

   /*= unamalgamate point matrices to block matrices */
   for (i = exp_Nrows - 1; i >= 0; i-- ) {
      for (j = num_PDE_eqns-1; j >= 0; j--) {
         aggr_index[i*num_PDE_eqns+j] = aggr_index[i];
      }
   }

   if ( mypid == 0 && 8 < ML_Get_PrintLevel())
   {
      printf("Calling ML_Operator_UnAmalgamateAndDropWeak\n");
      fflush(stdout);
   }

   ML_Operator_UnAmalgamateAndDropWeak(Amatrix, num_PDE_eqns, epsilon);

   Nrows      *= num_PDE_eqns;
   nvertices  *= num_PDE_eqns;
   exp_Nrows  *= num_PDE_eqns;
   getrow_obj  = Amatrix->getrow;
   N_neighbors = getrow_obj->pre_comm->N_neighbors;

#ifdef ML_AGGR_INAGGR

   for (i = 0; i < exp_Nrows; i++) aggr_index[i] = -1;
   sprintf(fname,"agg_%d",level_count); level_count++;
   fp = fopen(fname,"r");
   if (fp == NULL)
   {
      printf("Cannot open aggregate file %s for reading.\n",fname);
      printf("exp_Nrows = %d\n",exp_Nrows);
      exit(1);
   }
   else
   {
      printf("Reading aggregate for level %d from file %s\n",
             level_count-1,fname);
      fflush(stdout);
   }
   aggr_count = 0;
   for (i = 0; i <nvertices; i++) {
     if ( fscanf(fp,"%d%d",&j,&k) != 2) break;
      aggr_index[j] = k;
      if (k >= aggr_count) aggr_count = k+1;
   }
   fclose(fp);
   printf("Read in %d aggregates\n\n",aggr_count);
#endif /*ifdef ML_AGGR_INAGGR*/

   /* We have to recompute the send and receive lists. In particular, */
   /* we must assure that all degrees of freedom at a node point are  */
   /* in the send and receive lists (not just those in Amat).         */

   if ( num_PDE_eqns != 1 )
   {
      ML_memory_free((void**) &neighbors);
      ML_memory_free((void**) &recv_leng);
      ML_memory_free((void**) &send_leng);
      ML_memory_free((void**) &send_list);
      /* ---------------------------------------------------------- */
      /* allocate storage for the communication information         */
      /* ---------------------------------------------------------- */

      N_neighbors = getrow_obj->pre_comm->N_neighbors;
      nbytes = N_neighbors * sizeof( int );
      if ( nbytes > 0 )
      {
         ML_memory_alloc((void**) &neighbors,  nbytes, "AGL");
         ML_memory_alloc((void**) &recv_leng,  nbytes, "AGM");
         ML_memory_alloc((void**) &send_leng,  nbytes, "AGN");
      }
      else
      {
         neighbors = recv_leng = send_leng = NULL;
      }
      for ( i = 0; i < N_neighbors; i++ )
      {
         neighbors[i] = getrow_obj->pre_comm->neighbors[i].ML_id;
         recv_leng[i] = getrow_obj->pre_comm->neighbors[i].N_rcv;
         send_leng[i] = getrow_obj->pre_comm->neighbors[i].N_send;
      }
      total_recv_leng = 0;
      for ( i = 0; i < N_neighbors; i++ ) total_recv_leng += recv_leng[i];

      /* fix the send lists so that all degree of freedom associated with a */
      /* vertex in the original send_list is in the new send_list */

      /* first we need to count things up */
      nbytes = (Nrows + total_recv_leng*num_PDE_eqns) * sizeof( int );
               /* This needs to be big enough to hold all the receives */
      if ( nbytes > 0 ) ML_memory_alloc((void**) &dof_accounted_for,nbytes,"AGQ");

      count3 = 0;
      for ( i = 0; i < N_neighbors; i++ ) {
         for ( j = 0; j < Nrows; j++ ) dof_accounted_for[j] = -1;
         for (j = 0; j < send_leng[i]; j++) {
            index3 = getrow_obj->pre_comm->neighbors[i].send_list[j];
            k      = index3 / num_PDE_eqns;
            index3 = k * num_PDE_eqns;
            for (k = 0; k < num_PDE_eqns; k++) {
               if ( dof_accounted_for[index3+k] < 0 ) 
		 dof_accounted_for[index3+k] = count3++;
	    }
	 }
      }
      total_send_leng = count3;

      count3 = 0;
      for ( i = 0; i < N_neighbors; i++ ) {
         for ( j = nvertices-1; j < nvertices + total_recv_leng*num_PDE_eqns; j++ ) 
	   dof_accounted_for[j] = -1;
         for (j = 0; j < recv_leng[i]; j++) {
            index3 = getrow_obj->pre_comm->neighbors[i].rcv_list[j];
            k      = index3 / num_PDE_eqns;
            index3 = k * num_PDE_eqns;
            for (k = 0; k < num_PDE_eqns; k++) {
               if ( dof_accounted_for[index3+k] < 0 ) 
		 dof_accounted_for[index3+k] = count3++;
	    }
	 }
      }
      total_recv_leng = count3;

      nbytes = total_send_leng * sizeof( int );
      if ( nbytes > 0 ) ML_memory_alloc((void**) &send_list,nbytes,"AGO");
      else              send_list = NULL;
      nbytes = total_recv_leng * sizeof( int );
      if ( nbytes > 0 ) ML_memory_alloc((void**) &recv_list,nbytes,"AGP");
      else              recv_list = NULL;


      count3 = 0;
      for ( i = 0; i < N_neighbors; i++ ) {
	count = count3;
         for ( j = 0; j < Nrows; j++ ) dof_accounted_for[j] = -1;
         for (j = 0; j < send_leng[i]; j++) {
            index3 = getrow_obj->pre_comm->neighbors[i].send_list[j];
            k      = index3 / num_PDE_eqns;
            index3 = k * num_PDE_eqns;
            for (k = 0; k < num_PDE_eqns; k++) {
	      if ( dof_accounted_for[index3+k] < 0 ) {
		send_list[count3] = index3 + k;
		dof_accounted_for[index3+k] = count3++;
	      }
	    }
	 }
	 send_leng[i] = count3 - count;
      }
      count3 = 0;
      for ( i = 0; i < N_neighbors; i++ ) {
	count = count3;
	for ( j = nvertices-1; j < nvertices + total_recv_leng*num_PDE_eqns; j++ ) 
	  dof_accounted_for[j] = -1;
         for (j = 0; j < recv_leng[i]; j++) {
            index3 = getrow_obj->pre_comm->neighbors[i].rcv_list[j];
            k      = index3 / num_PDE_eqns;
            index3 = k * num_PDE_eqns;
            for (k = 0; k < num_PDE_eqns; k++) {
	      if ( dof_accounted_for[index3+k] < 0 ) {
		recv_list[count3] = index3 + k;
		dof_accounted_for[index3+k] = count3++;
	      }
	    }
	 }
	 recv_leng[i] = count3 - count;
      }
      ML_memory_free((void**) &dof_accounted_for);
      exp_Nrows = Nrows + total_recv_leng;
   }


   /* count the size of each aggregate */

   aggr_cnt_array = (int *) ML_allocate(sizeof(int)*aggr_count);
   for (i = 0; i < aggr_count ; i++) aggr_cnt_array[i] = 0;
   for (i = 0; i < exp_Nrows; i++) 
      if (aggr_index[i] >= 0) 
         aggr_cnt_array[aggr_index[i]]++;

/*
for (i = 0; i < exp_Nrows; i++) printf("%d: AGG_INDX %d %d\n",comm->ML_mypid, i,aggr_index[i]);
for (i = 0; i < aggr_count ; i++) printf("counts %d %d\n",i,aggr_cnt_array[i]);
*/

#ifdef ML_AGGR_OUTAGGR
   sprintf(fname,"agg%d_%d",comm->ML_mypid,level_count);
   fp = fopen(fname,"w");
   agg_offset = ML_gpartialsum_int(aggr_count, comm);
   vertex_offset = ML_gpartialsum_int(nvertices, comm);
   for (i = 0; i < nvertices ; i++)
   {
#ifndef MAXWELL
#ifdef ALEGRA
      if (level_count == 0) { j = i; k = i+vertex_offset;}
#else
      if (level_count == 0) { j = update_index[i]; k = update[i];}
#endif /*ifdef ALEGRA*/
#else
      if (level_count == 0) { j = reordered_glob_nodes[i]; k = global_node_inds[i];}
#endif /* ifndef MAXWELL */
      else                  { j = i              ; k = i+vertex_offset;}
      if (aggr_index[j] >= 0)
         fprintf(fp,"%5d %5d\n",k, aggr_index[j]+agg_offset);
   }

   dtemp = (double *) ML_allocate(sizeof(double)*(exp_Nrows+1));
   for (i = 0; i < nvertices; i++) dtemp[i] = (double) (i + vertex_offset);
   ML_exchange_bdry(dtemp,Amatrix->getrow->pre_comm, nvertices, comm,
                    ML_OVERWRITE,NULL);
   for (i = 0; i < exp_Nrows-nvertices; i++)
   {
#ifndef MAXWELL
#ifdef ALEGRA
      if (level_count == 0) { j = i+nvertices; k =  (int) dtemp[i+nvertices];}
#else
      if (level_count == 0) { j = extern_index[i]; k = external[i];} 
#endif /*ifdef ALEGRA*/
#else
      if (level_count == 0) { j = reordered_node_externs[i]; k =
	  global_node_externs[i];}
#endif /* ifndef MAXWELL */

      else                 { j = i+nvertices    ; k = (int) dtemp[i+nvertices];}
      if (aggr_index[j] >= 0)
         fprintf(fp,"%5d %5d\n", k, aggr_index[j]+agg_offset);
   }
   fclose(fp);
   level_count++;
   ML_free(dtemp);
#else
#ifdef INPUT_AGGREGATES
   agg_offset = ML_gpartialsum_int(aggr_count, comm);
   vertex_offset = ML_gpartialsum_int(nvertices, comm);
#endif
#endif /*ifdef ML_AGGR_OUTAGGR*/

   /* ============================================================= */
   /* Form tentative prolongator                                    */
   /* ============================================================= */

   Ncoarse = aggr_count;

   /* ============================================================= */
   /* update aggr_index to find out which local fine node is mapped */
   /* to which coarse aggregate in remote processors                */
   /* ------------------------------------------------------------- */

/****************************************************************************
*****************************************************************************

We want to assign local ghost aggregate numbers to local nodes that are 
within another processors aggregate. We also need to update communication 
information so that it now corresponds to aggregates.

Here is how we do all this:
  1) int_buf2[.] corresponds to rcv_list[] and contains the local aggr_index[]
     or -1 (if not assigned to any local aggregate).
  2) int_buf corresponds to send_list[] and contains remote aggr_index[] or 
     -1 (if not assigned to the corresponding remote processor). 'int_buf' is 
     filled via  ML_Aggregate_ExchangeData().
  3) for each neighbor ...
       a) find largest remote aggregate index in int_buf
       b) reallocate int_buf2 to hold largest remote index
       c) compute int_buf2 such that int_buf2[i] = 
          # of local nodes assigned to remote aggregate i.
       d) count number of remote aggregates for this neighbor
          and set int_buf2[i] = 1 or 0 depending on whether there
          is a remote aggregate with this neighbor.
       e) redefine int_buf2[i] in the following way
              int_buf2[i+1] - int_buf2[i] = 1 indicates that there is
              a remote aggregate i and int_buf2[i] is the local ghost 
              number associated with that remote aggregate.
       f) record the appropriate recv_length, recv_neighbors, etc.
       g) Now jam in the right ghost number into aggr_index[]. I'm
          not really sure why we need the first condition in the
          double condition if statement?

*****************************************************************************
****************************************************************************/

   nbytes = total_send_leng * sizeof(int);
   if ( nbytes > 0 ) ML_memory_alloc((void**) &int_buf, nbytes, "AGg");
   else              int_buf = NULL;
   nbytes = total_recv_leng * sizeof(int);
   if ( nbytes > 0 ) ML_memory_alloc((void**) &int_buf2, nbytes, "AGh");
   else              int_buf2 = NULL;

   /* ------------------------------------------------------------- */
   /* send the remote node index back to remote processors, with   */
   /* added information on which remote nodes have been aggregated */
   /* by the local aggregates (and also the aggregate numbers).    */
   /* ------------------------------------------------------------- */

   offset = 0;
   for ( i = 0; i < N_neighbors; i++ ) 
   {
      for ( j = 0; j < recv_leng[i]; j++ ) 
      {
         if ( aggr_index[Nrows+offset+j] < 0 ) int_buf2[offset+j] = -1;
         else int_buf2[offset+j] = aggr_index[Nrows+offset+j];
      }
      offset += recv_leng[i];
   }
   msgtype = 15963;
   ML_Aggregate_ExchangeData((char*) int_buf, (char*) int_buf2,
      N_neighbors, neighbors, send_leng, recv_leng, msgtype, ML_INT, comm);

   if ( int_buf2 != NULL ) ML_memory_free((void**) &int_buf2);

   /* ------------------------------------------------------------- */
   /* if int_buf[i] > 0, this means that aggr_index[send_list[i]]   */ 
   /* has been aggregated by a remote processor                     */
   /* int_buf2 is used to tabulate how many distinct aggregates     */
   /* in remote processors are used.                                */
   /* ------------------------------------------------------------- */

   offset = 0;
   m      = 0; /* store the local index offset for remote processors */ 
   new_N_recv = 0;
   nbytes = N_neighbors * sizeof(int);
   if ( nbytes > 0 ) 
   {
      ML_memory_alloc((void**) &new_recv_leng, nbytes, "AGi");
      ML_memory_alloc((void**) &new_recv_neighbors, nbytes, "AGj");
   } 
   else 
   {
      new_recv_leng = new_recv_neighbors = NULL;
   }
   for ( i = 0; i < N_neighbors; i++ ) 
   {
      /* ---------------------------------------------------------- */
      /* find out how large an array to allocate for int_buf2       */
      /* ---------------------------------------------------------- */

      max_count = -1;
      for ( j = 0; j < send_leng[i]; j++ ) 
      {
         index = int_buf[offset+j];
         max_count = (index > max_count ) ? index : max_count;
      }
      nbytes = ( max_count + 2 ) * sizeof(int);
      if (nbytes > 0) ML_memory_alloc((void **) &int_buf2, nbytes, "AGk");

      /* ---------------------------------------------------------- */
      /* see how many distinct remote aggregates are referenced by  */
      /* local fine nodes in aggregation in proc i ==> count        */
      /* ---------------------------------------------------------- */

      for ( j = 0; j <= max_count; j++ ) int_buf2[j] = 0;
      for ( j = 0; j < send_leng[i]; j++ ) 
      {
         index = int_buf[offset+j];
         if ( index >= 0 ) int_buf2[index]++;
	 /*= int_blk2[k] = # of pts we own that are in aggregate k
	   on processor i */
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

      if ( count > 0 ) 
      {
         new_recv_leng[new_N_recv] = count * nullspace_dim;
         new_recv_neighbors[new_N_recv] = neighbors[i];
         new_N_recv++;
      } 

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

	 /* I have replaced the double condition with just a single
	    one. I hope I haven't made a mistake, rst 
	    if ( aggr_index[index] <= -100 && int_buf[offset+j] >= 0 )*/ 
         if ( int_buf[offset+j] >= 0 ) 
         {
            k = int_buf[offset+j];
            aggr_index[index] = int_buf2[k] + Ncoarse + m;
 	    /*= this is the ghost number for remote aggregates */
         } 
      }
      if (nbytes > 0) ML_memory_free((void **) &int_buf2);
      m += count;
      offset += send_leng[i];
   }
   exp_Ncoarse = Ncoarse + m;
 
   if ( int_buf  != NULL ) ML_memory_free((void**) &int_buf);

   /* ============================================================= */
   /* check and copy aggr_index                                     */
   /* ------------------------------------------------------------- */

   /*= record aggregates in ml_ag */
   level = ml_ag->cur_level;
   nbytes = Nrows * sizeof( int );
   ML_memory_alloc((void**) &(ml_ag->aggr_info[level]), nbytes, "AGl");
   count = aggr_count;
   for ( i = 0; i < Nrows; i+=num_PDE_eqns ) 
   {
      if ( aggr_index[i] >= 0 )
      {
         for ( j = 0; j < num_PDE_eqns; j++ ) 
            ml_ag->aggr_info[level][i+j] = aggr_index[i];
         if (aggr_index[i] >= count) count = aggr_index[i] + 1;
      }
      /*else
       *{
       *   printf("%d : CoarsenMIS error : aggr_index[%d] < 0\n",
       *          mypid,i);
       *   exit(1);
       *}*/
   }
   ml_ag->aggr_count[level] = count; /* for relaxing boundary points */
  
   /* ============================================================= */
   /* find out how many local coarse aggregates are needed by       */
   /* remote processors for interpolation (to construct the         */
   /* communicator - send info - for P)                             */
   /* ------------------------------------------------------------- */

   new_N_send = 0;
   if ( N_neighbors > 0 ) 
   {
      nbytes = N_neighbors * sizeof(int);
      ML_memory_alloc((void**) &int_buf, nbytes, "AGm");
      nbytes = Ncoarse * sizeof(int);
      ML_memory_alloc((void**) &int_buf2, nbytes, "AGn");
      for ( i = 0; i < N_neighbors; i++ ) int_buf[i] = 0;

      /* ---------------------------------------------------------- */
      /* count which remote fine nodes belong to local aggregates   */
      /* in order to generate the communication pattern for         */
      /* the interpolation operator.                                */
      /* ---------------------------------------------------------- */

      offset = Nrows; 
      for ( i = 0; i < N_neighbors; i++ ) 
      {
         for ( j = 0; j < Ncoarse; j++ ) int_buf2[j] = 0;
         for ( j = 0; j < recv_leng[i]; j++ ) 
         {
            index = aggr_index[offset++];
            if ( index >= 0 ) int_buf2[index]++;
         }
         count = 0;
         for ( j = 0; j < Ncoarse; j++ ) if ( int_buf2[j] > 0 ) count++;
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
         ML_memory_alloc((void**) &new_send_leng, nbytes, "AGo");
         ML_memory_alloc((void**) &new_send_neighbors, nbytes, "AGp");
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
         ML_memory_alloc((void**) &new_send_list, nbytes, "AGq");
         offset = Nrows;
         m = count;
         count = 0;
         for ( i = 0; i < N_neighbors; i++ ) 
         {
            for ( j = 0; j < Ncoarse; j++ ) int_buf2[j] = 0;
            for ( j = 0; j < recv_leng[i]; j++ ) 
            {
               index = aggr_index[offset++];
               if ( index >= 0 ) int_buf2[index]++;
            }
            for ( j = 0; j < Ncoarse; j++ ) 
            {
               if ( int_buf2[j] > 0 ) 
               {
                  for ( jj = 0; jj < nullspace_dim; jj++ ) 
                     new_send_list[count++] = j * nullspace_dim + jj;
		  /*= new_send_list is the processor's send list */
               } 
            } 
         } 
         if ( m != count ) 
         {
            printf("ML_Aggregate_CoarsenMIS : internal error (1).\n");
            exit(-1);
         }
      } 
      else 
      {
         new_send_leng = NULL;
         new_send_neighbors = NULL;
         new_send_list = NULL;
      }  
      ML_memory_free((void**) &int_buf);
      ML_memory_free((void**) &int_buf2);
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
      {
         printf("WARNING : index out of bound %d = %d(%d)\n", i, aggr_index[i], 
                exp_Ncoarse);
/*
         for ( j = 0; j < new_Nrows; j++ ) 
            if ( aggr_index[j] >= exp_Ncoarse )
               printf("%d : aggr_index[%5d] = %5d *\n", mypid, j, aggr_index[j]); 
            else
               printf("%d : aggr_index[%5d] = %5d \n", mypid, j, aggr_index[j]); 
         exit(1);
*/
      }
   }
   nbytes = ( new_Nrows + 1 ) * sizeof(int); 
   ML_memory_alloc((void**)&(new_ia), nbytes, "AIA");
   nbytes = new_Nrows * nullspace_dim * sizeof(int); 
   ML_memory_alloc((void**)&(new_ja), nbytes, "AJA");
   nbytes = new_Nrows * nullspace_dim * sizeof(double); 
   ML_memory_alloc((void**)&(new_val), nbytes, "AVA");
   for ( i = 0; i < new_Nrows*nullspace_dim; i++ ) new_val[i] = 0.0;

   /* ------------------------------------------------------------- */
   /* set up the space for storing the new null space               */
   /* ------------------------------------------------------------- */

   nbytes = Ncoarse * nullspace_dim * nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&(new_null),nbytes,"AGr");
   for (i = 0; i < Ncoarse*nullspace_dim*nullspace_dim; i++) 
      new_null[i] = 0.0;

   /* ------------------------------------------------------------- */
   /* initialize the row pointer for the CSR prolongation operator  */
   /* (each row will have at most nullspace_dim nonzero entries)    */
   /* ------------------------------------------------------------- */

   for (i = 0; i <= Nrows; i++) new_ia[i] = i * nullspace_dim;

   /* trying this when a Dirichlet row is taken out */
   j = 0;
   new_ia[0] = 0;
   for (i = 0; i < Nrows; i++) {
      if (aggr_index[i] != -1) j += nullspace_dim;
      new_ia[i+1] = j;
   }

   /* ------------------------------------------------------------- */
   /* generate an array to store which aggregate has which rows.Then*/
   /* loop through the rows of A checking which aggregate each row  */
   /* is in, and adding it to the appropriate spot in rows_in_aggs  */
   /* ------------------------------------------------------------- */

   ML_memory_alloc((void**)&rows_in_aggs,aggr_count*sizeof(int*),"MLs");
   for (i = 0; i < aggr_count; i++) 
   {
      rows_in_aggs[i] = (int *) ML_allocate(aggr_cnt_array[i]*sizeof(int));
      aggr_cnt_array[i] = 0;
      if (rows_in_aggs[i] == NULL) 
      {
         printf("Error: couldn't allocate memory in CoarsenMIS\n");
         exit(1);
      }
   }
   for (i = 0; i < exp_Nrows; i+=num_PDE_eqns) 
   {
      if ( aggr_index[i] >= 0 && aggr_index[i] < aggr_count)
      {
         for (j = 0; j < num_PDE_eqns; j++)
         {
            index = aggr_cnt_array[aggr_index[i]]++; 
            rows_in_aggs[aggr_index[i]][index] = i + j;
         }
      }
   }

   /* ------------------------------------------------------------- */
   /* allocate work arrays for QR factorization                     */
   /* work and lwork are needed for lapack's QR routine.  These     */
   /* settings seemed easiest since I don't quite understand        */
   /* what they do, but may want to do something better here later  */
   /* ------------------------------------------------------------- */

   nbytes = total_recv_leng * nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&comm_val, nbytes, "AGt");
   for (i = 0; i < total_recv_leng*nullspace_dim; i++) comm_val[i] = 0.0; 
   max_agg_size = 0;
   for (i = 0; i < aggr_count; i++) 
   {
      if (aggr_cnt_array[i] > max_agg_size) max_agg_size = aggr_cnt_array[i];
   }
   nbytes = max_agg_size * nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&qr_tmp, nbytes, "AGu");
   nbytes = nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&tmp_vect, nbytes, "AGv");

   lwork  = nullspace_dim;
   nbytes = nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&work, nbytes, "AGw");

   /* ------------------------------------------------------------- */
   /* ship the null space information to other processors           */
   /* ------------------------------------------------------------- */
 
   if (nullspace_vect != NULL) 
   {
      nbytes = total_send_leng * nullspace_dim * sizeof(double);
      ML_memory_alloc((void**) &dble_buf, nbytes,"AG1");
      nbytes = total_recv_leng * nullspace_dim * sizeof(double);
      ML_memory_alloc((void**) &dble_buf2, nbytes,"AG2");
      length = total_send_leng * nullspace_dim;
      for ( i = 0; i < total_send_leng; i++ ) 
      {
         index = send_list[i];
         for ( j = 0; j < nullspace_dim; j++ ) 
            dble_buf[i*nullspace_dim+j] = nullspace_vect[j*Nrows+index];
      }
      msgtype = 12093;
      length = sizeof(double) * nullspace_dim;
      ML_Aggregate_ExchangeData((char*)dble_buf2,(char*) dble_buf,
            N_neighbors, neighbors, recv_leng, send_leng,msgtype,
				(int) length,comm);
      ML_memory_free((void**) &dble_buf);
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
            for (k = 0; k < nullspace_dim; k++)
            {
              if ( unamalg_bdry[index] == 'T') qr_tmp[k*length+j] = 0.;
               else
               {
                  if (index % num_PDE_eqns == k) qr_tmp[k*length+j] = 1.0;
                  else                           qr_tmp[k*length+j] = 0.0;
               }
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
               if ( unamalg_bdry[index] == 'T') qr_tmp[k*length+j] = 0.;
               else {
                  if (index < Nrows) {
                     qr_tmp[k*length+j] = nullspace_vect[k*Nrows+index];
                  }
                  else {
                     qr_tmp[k*length+j] = 
                        dble_buf2[(index-Nrows)*nullspace_dim+k];
                  }
               }
            }
         }
      }

      /* ---------------------------------------------------------- */
      /* now calculate QR using an LAPACK routine                   */
      /* ---------------------------------------------------------- */
      if (aggr_cnt_array[i] >= nullspace_dim) {

	MLFORTRAN(dgeqrf)(&(aggr_cnt_array[i]), &nullspace_dim, qr_tmp, 
			  &(aggr_cnt_array[i]), tmp_vect, work, &lwork, &info);
	if (info != 0)
	  pr_error("ErrOr in CoarsenMIS : dgeqrf returned a non-zero %d %d\n",
		   aggr_cnt_array[i],i);

	if (work[0] > lwork) 
	  {
	    lwork=(int) work[0]; 
	    ML_memory_free((void**) &work);
	    ML_memory_alloc((void**) &work, sizeof(double)*lwork, "AGx");
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

	if ( aggr_cnt_array[i] < nullspace_dim ){
	  printf("Error in dorgqr on %d row (dims are %d, %d)\n",i,aggr_cnt_array[i],
                 nullspace_dim);
	  printf("ERROR : performing QR on a MxN matrix where M<N.\n");
	}
	MLFORTRAN(dorgqr)(&(aggr_cnt_array[i]), &nullspace_dim, &nullspace_dim, 
			  qr_tmp, &(aggr_cnt_array[i]), tmp_vect, work, &lwork, &info);
	if (info != 0) {
	  printf("Error in dorgqr on %d row (dims are %d, %d)\n",i,aggr_cnt_array[i],
                 nullspace_dim);
	  pr_error("Error in CoarsenMIS: dorgqr returned a non-zero\n");
	}

	if (work[0] > lwork) 
	  {
	    lwork=(int) work[0]; 
	    ML_memory_free((void**) &work);
	    ML_memory_alloc((void**) &work, sizeof(double)*lwork, "AGy");
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
      else {
	/* We have a small aggregate such that the QR factorization can not */
	/* be performed. Instead let us copy the null space from the fine   */
        /* into the coarse grid nullspace and put the identity for the      */
	/* prolongator????                                                  */
	for (j = 0; j < nullspace_dim; j++)
	  for (k = 0; k < nullspace_dim; k++)
            new_null[i*nullspace_dim+j+k*Ncoarse*nullspace_dim] = 
	      qr_tmp[j+aggr_cnt_array[i]*k];
	for (j = 0; j < aggr_cnt_array[i]; j++) {
	  index = rows_in_aggs[i][j];
	  index3 = new_ia[index];
	  for (k = 0; k < nullspace_dim; k++) {
	    new_ja [index3+k] = i * nullspace_dim + k;
	    if (k == j) new_val[index3+k] = 1.;
	    else new_val[index3+k] = 0.;
	  }
	}
      }


   }
	 
   ML_Aggregate_Set_NullSpace(ml_ag, num_PDE_eqns, nullspace_dim, 
                              new_null, Ncoarse*nullspace_dim);
   ML_memory_free( (void **) &new_null);
   if (nullspace_vect != NULL) ML_memory_free( (void **) &dble_buf2);

   /* ------------------------------------------------------------- */
   /* send the P rows back to its parent processor                  */
   /* ------------------------------------------------------------- */
 
   nbytes = total_send_leng * nullspace_dim * sizeof(double);
   ML_memory_alloc((void**) &dble_buf, nbytes,"AGz");
   msgtype = 24945;
   length = sizeof(double) * nullspace_dim;
   ML_Aggregate_ExchangeData((char*)dble_buf,(char*) comm_val,
         N_neighbors, neighbors, send_leng, recv_leng,msgtype,(int)length,comm);
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
               /*new_val[k+j] = dble_buf[index4+j];*/
               new_val[ new_ia[index]+j] = dble_buf[index4+j];
               /* new_ja = column index */
               /*new_ja[k+j]  = aggr_index[index]*nullspace_dim+j;*/
               new_ja[ new_ia[index]+j ]  = aggr_index[index]*nullspace_dim+j;
            }
         }
      }
   }
   ML_memory_free( (void **) &comm_val);
   ML_memory_free( (void **) &dble_buf);
 
   /* ------------------------------------------------------------- */
   /* check P (row sum = 1)                                         */
   /* ------------------------------------------------------------- */

#ifdef DDEBUG
   /* count up the number of connections/edges that leave aggregate */
   /* versus those that are within the aggregate.                   */
   good = (int *) ML_allocate(Nrows*sizeof(int));
   bad  = (int *) ML_allocate(Nrows*sizeof(int));
   for (i = 0; i < Nrows; i++) { good[i] = 0; bad[i] = 0; }
   count = 0; 
   for (i = 0; i < Nrows; i++) {
      /* figure out my aggregate */
      myagg = -1;
      for (kk = new_ia[i]; kk < new_ia[i+1]; kk++) {
         if (myagg != new_ja[kk]/6) {
            if (myagg == -1) myagg = new_ja[kk]/6;
            else {
               printf("a:something is wrong %d %d in row %d\n",
                       myagg, new_ja[kk]/6, i);
               /*exit(1);*/
            } 
          }
      }

      ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
                        &rowi_N, 0);
      count2 = 0;
      for (j = 0; j < rowi_N; j++) {
         if (rowi_col[j]/num_PDE_eqns != -500 - i/num_PDE_eqns) {
            /* for each column figure out what aggregate it corresponds */
	    /* to. If it is the same as myagg, add 1 to good, otherwise */
	    /* add 1 to bad */
            curagg = -1;
	    for (kk = new_ia[rowi_col[j]]; kk < new_ia[rowi_col[j]+1]; kk++) {
	      if (curagg != new_ja[kk]/6) {
		if (curagg == -1) curagg = new_ja[kk]/6;
		else {
		  printf("b:Something is wrong %d %d in row %d\n",
			 curagg, new_ja[kk]/6, rowi_col[j]);
		  /*exit(1);*/
		} 
	      }
	    }
            if ((curagg != -1) && (myagg != -1)) {
               if (curagg == myagg) good[myagg]++;
               else bad[myagg]++;
            }
            
         }
      }
   }
   myagg = 0;
   sprintf(fname,"goodbad%d",level_count);
   fp = fopen(fname,"w");
   for (i = 0; i < Nrows; i++) { 
      if ((good[i] != 0) || (bad[i] != 0)) {
         myagg += good[i]; 
         myagg += bad[i]; 
         fprintf(fp,"%d (%d,%d)\n",i,good[i],bad[i]);
      }
   }
   fclose(fp);
   printf("total number of connections counted is %d\n",myagg);
   ML_free(bad);
   ML_free(good);
   ML_free(rowi_col); rowi_col = NULL;
   ML_free(rowi_val); rowi_val = NULL;
   allocated = 0;
#endif
   /* ------------------------------------------------------------- */
   /* set up the csr_data data structure                            */
   /* ------------------------------------------------------------- */

   ML_memory_alloc((void**) &csr_data, sizeof(struct ML_CSR_MSRdata),"CSR");
   csr_data->rowptr  = new_ia;
   csr_data->columns = new_ja;
   csr_data->values  = new_val;

   ML_Operator_Set_ApplyFuncData( *Pmatrix, nullspace_dim*Ncoarse, Nrows, 
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
#ifdef LEASTSQ_SERIAL
   ML_Operator_Set_Getrow((*Pmatrix), ML_EXTERNAL, Nrows, CSR_get_ones_rows);
#else
   ML_Operator_Set_Getrow((*Pmatrix), ML_EXTERNAL, Nrows, CSR_getrows);
#endif
   ML_Operator_Set_ApplyFunc((*Pmatrix), ML_INTERNAL, CSR_matvec);
   (*Pmatrix)->max_nz_per_row = 1;

   /* ------------------------------------------------------------- */
   /* clean up                                                      */
   /* ------------------------------------------------------------- */

   ML_free(unamalg_bdry);
   ML_memory_free((void**) &comm_val);
   ML_memory_free((void**) &neighbors);
   ML_memory_free((void**) &recv_leng);
   ML_memory_free((void**) &send_leng);
   ML_memory_free((void**) &send_list);
   ML_free(aggr_index);
   ML_memory_free((void**) &aggr_stat);
   ML_memory_free((void**) &sendlist_proc);
   ML_free(aggr_cnt_array);
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
   if ( N_neighbors > 0 ) 
   {
      ML_memory_free((void**) &new_recv_leng);
      ML_memory_free((void**) &new_recv_neighbors);
   }
   ML_memory_free((void**) &aggr_comm);
   aggr_curr = aggr_head;
   while ( aggr_curr != NULL ) 
   {
      supernode = aggr_curr;
      aggr_curr = aggr_curr->next;
      if ( supernode->length > 0 ) ML_free( supernode->list );
      ML_free( supernode );
   }
#if defined(OUTPUT_AGGREGATES) || defined(INPUT_AGGREGATES)
   /* Print Pmatrix*v (where v is constructed using global indices) */

   dtemp = (double *) ML_allocate(sizeof(double)*(*Pmatrix)->invec_leng);
   d2temp = (double *) ML_allocate(sizeof(double)*(*Pmatrix)->outvec_leng);
   for (i = 0; i < (*Pmatrix)->outvec_leng; i++) d2temp[i] = 0.;
   for (i = 0; i < (*Pmatrix)->invec_leng; i++)
      dtemp[i] = (double) (i+agg_offset);

   sprintf(fname,"PP%d_%d",comm->ML_mypid,level_count);
   fp = fopen(fname,"w");
   ML_Operator_Apply(*Pmatrix, (*Pmatrix)->invec_leng, dtemp, 
                     (*Pmatrix)->outvec_leng, d2temp);
   for (i = 0; i < nvertices; i++) {
#ifndef MAXWELL
#ifdef ALEGRA
      if (level_count == 1) { j = i; k = i;} 
#else
      if (level_count == 1) { j = update_index[i]; k = update[i];} 
#endif
#else
      if (level_count == 1) { j = reordered_glob_nodes[i]; k = global_node_inds[i];}
#endif /* ifndef MAXWELL */
      else                  { j = i              ; k = i+vertex_offset;}
      fprintf(fp,"PP%d(%d) (loc=%d) = %e\n",level_count,k,j, d2temp[j]);
   }
   fclose(fp);
   ML_free(dtemp);
   ML_free(d2temp);

   /*
   csr_data = (struct ML_CSR_MSRdata *) (*Pmatrix)->data;
   if (comm->ML_mypid == 1)
   {
      printf("%d : row_ptr = %d\nrow_ptr = %d\ncol = %d\nval = %e\n\n\n",
             comm->ML_mypid,
             csr_data->rowptr[14], csr_data->rowptr[15],
             csr_data->columns[14], csr_data->values[14]);
   }
   sprintf(fname,"comm%d",comm->ML_mypid);
   ML_CommInfoOP_Print((*Pmatrix)->getrow->pre_comm, fname);
   fflush(stdout);
   */

#endif
   return Ncoarse*nullspace_dim;
}
#include "ml_ggraph.h"

/* ******************************************************************** */
/* A subroutine to label vertices of a particular type                  */
/* -------------------------------------------------------------------- */

int ML_Aggregate_LabelVertices(int vlist_cnt, int *vlist, char Vtype,
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
   int     *pref_list, col2, loop_cnt, *tlist, pref_cnt;
   int     pref_flag, pref_index;
   char    *in_preflist;
   USR_REQ *Request;
   int msg_type = 1041;
   unsigned int nbytes;
#ifdef RANDMISROOTS
   double *tvec2;
   int *tvec, *neworder;
#endif

   ML_avoid_unused_param((void *) &msgtype);
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

   for ( i = 0; i < vlist_cnt; i++ ) {
      index = vlist[i];
      if (vertex_state[index] == 'D') {
         N_remaining_vertices--;
         for ( k = 0; k < proclist[index][0]; k++ ) {
            fproc = proclist[index][2*k+1];
            m     = proclist[index][2*k+2];
            send_buf[fproc][m] = 2;
#ifdef ML_AGGR_DEBUG
            printf("%d: Recording %d (%d)\n",comm->ML_mypid,fproc,m); 
            fflush(stdout);
#endif
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

   in_preflist = (char *) ML_allocate(vlist_cnt*sizeof(char) );
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
   /* Why is pref_cnt 0????? shouldn't it be m????? */
   /* perhaps it is because nothing is aggregated yet? */
   
   /* -------------------------------------------------------- */
   /* get ready for the coarsening                             */
   /* -------------------------------------------------------- */

#ifdef ML_AGGR_DEBUG
printf("%d: N_remaining_vertices is %d |%d\n",comm->ML_mypid, 
        N_remaining_vertices,recv_cnt); fflush(stdout);
#endif
   nselected = 0;

   /* -------------------------------------------------------- */
   /* let's actually do coarsening                             */
   /* -------------------------------------------------------- */

   change_flag = 1;
   loop_cnt = 0;
   pref_index = 0;     /* pointer to a stack of vertex numbers */

#ifdef RANDMISROOTS
   tvec2 = (double *) ML_allocate( (vlist_cnt+1)*sizeof(double));
   tvec  = (int *) ML_allocate( (vlist_cnt+1)*sizeof(int));
   ML_random_vec(tvec2, vlist_cnt, comm);
   for (i = 0; i < vlist_cnt; i++) tvec[i] = (int) (100000.*tvec2[i]);
   ML_free(tvec2);
   neworder  = (int *) ML_allocate( (vlist_cnt+1)*sizeof(int));
   for (i = 0; i < vlist_cnt; i++) neworder[i] = i;
   ML_az_sort(tvec, vlist_cnt, neworder, NULL);
   ML_free(tvec);
#endif


   do {
      pref_index = 0;     /* pointer to a stack of vertex numbers */
#ifdef RANDMISROOTS
   pref_index = vlist_cnt + 1000;
#endif

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
#ifdef RANDMISROOTS
         index = vlist[neworder[i]];
#endif
#ifdef ML_AGGR_DEBUG
printf("%d: near the top%d, (%d %d) | %d %d %d\n",comm->ML_mypid,
        N_remaining_vertices,i,index,pref_index,pref_cnt,pref_list[pref_index]); 
fflush(stdout);
#endif
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
#ifdef ML_AGGR_DEBUG
printf("%d: my(%d) neighs(%d) are %c\n",comm->ML_mypid,index,col,
       vertex_state[col]); fflush(stdout);
#endif
               if ( vertex_state[col] == 'S' )
               {
                  if (vertex_state[index] == 'F') {
                     vertex_state[index] = 'D';
                     N_remaining_vertices--;
                     for ( k = 0; k < proclist[index][0]; k++ ) {
                        fproc = proclist[index][2*k+1];
                        m     = proclist[index][2*k+2];
                        send_buf[fproc][m] = 2;
#ifdef ML_AGGR_DEBUG
printf("%d: recording %d (%d)\n",comm->ML_mypid,fproc,m); fflush(stdout);
#endif
                        change_flag = 1;
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
#ifdef ML_AGGR_DEBUG
                        printf("%d : %d not selected due to N %d\n",myrank,
                               index,nvertices);
#endif
                        select_flag = 0;
                        break;
                     }
                  }
               }
            }
#ifdef ML_AGGR_DEBUG
printf("%d: out of loop %d | %d\n",comm->ML_mypid,N_remaining_vertices,
       select_flag); fflush(stdout);
#endif

            /* if the vertex in question is not any of those considered */
            /* above, select this vertex.                               */

            if ( select_flag == 1 )
            {
#ifdef ML_AGGR_DEBUG
               printf("%d : %d selected %d\n", myrank, index,nvertices); 
               fflush(stdout);
#endif
               if ((vertex_state[index] == 'F') &&(index < nvertices)) 
                  N_remaining_vertices--;
               vertex_state[index] = 'S';
               aggr_index[index] = nselected;
               nselected++;

               /* set the flag that this vertex has been selected in */
               /* the buffer which is to be sent to other processors */

               for ( k = 0; k < proclist[index][0]; k++ )
               {
                  fproc = proclist[index][2*k+1];
                  m     = proclist[index][2*k+2];
                  send_buf[fproc][m] = 1;
                  change_flag = 1;
               }

               /* delete all vertices adjacent to this vertex and */
	       /* INDICATE that also in the communication buffer  */

               for ( j = rptr[index]; j < rptr[index+1]; j++ )
               {
                  col = cptr[j];
#ifdef ML_AGGR_DEBUG
printf("%d: deleting %d\n",comm->ML_mypid,col); fflush(stdout);
#endif
                  if (vertex_state[col] == 'F') {
                     vertex_state[col] = 'D';
                     if ( col < nvertices ) {
                        N_remaining_vertices--;
                        for ( k = 0; k < proclist[col][0]; k++ ) {
                           fproc = proclist[col][2*k+1];
                           m     = proclist[col][2*k+2];
                           send_buf[fproc][m] = 2;
#ifdef ML_AGGR_DEBUG
printf("%d: recording %d (%d)\n",comm->ML_mypid,fproc,m); fflush(stdout);
#endif
                           change_flag = 1;
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
#ifdef ML_AGGR_DEBUG
printf("%d: near the bottom %d\n",comm->ML_mypid,N_remaining_vertices); 
fflush(stdout);
#endif

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
#ifdef ML_AGGR_DEBUG
printf("%d: N_remaining_vertices is %d  | %d\n",comm->ML_mypid, 
        N_remaining_vertices,send_leng[j]); fflush(stdout);
printf("%d: a couple %d %d \n", comm->ML_mypid, send_buf[j][0], 
        send_buf[j][1]);
#endif
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
            comm->USR_cheapwaitbytes((char*) recv_buf[j], nbytes, &fproc,
                     &msgtype, comm->USR_comm, (void *) &Request[j] );
#ifdef ML_AGGR_DEBUG
printf("%d: rcvd is %d\n",comm->ML_mypid, recv_buf[0][0]); fflush(stdout);
#endif
            for ( k = 0; k < recv_leng[j]; k++ )
            {
               kkk = recv_list[j][k];
               if      (recv_buf[j][k] == 1)
               {
                  vertex_state[kkk] = 'S';
#ifdef ML_AGGR_DEBUG
                  printf("%d : incoming state %d = S \n",myrank,kkk);
#endif
               }
               else if (recv_buf[j][k] == 2)
               {
                  vertex_state[kkk] = 'D';
#ifdef ML_AGGR_DEBUG
                  printf("%d : incoming state %d = D \n",myrank,kkk);
#endif
               }
            }
#ifdef ML_AGGR_DEBUG
printf("%d: what are these %d %d\n",
comm->ML_mypid, recv_leng[j],recv_buf[j][recv_leng[j]]); fflush(stdout);
#endif
            if ( recv_buf[j][recv_leng[j]] == 1 )
            {
               proc_flag[j] = 1;
               NremainingRcvProcs--;
#ifdef ML_AGGR_DEBUG
printf("the new NremainingRcvProcs is %d\n",NremainingRcvProcs); fflush(stdout);
#endif
            }
         }
      }
#endif
   } while ( NremainingRcvProcs > 0 || N_remaining_vertices > 0 );

#ifdef ML_AGGR_DEBUG
   printf("%d : loop_count = %d \n", myrank, loop_cnt ); fflush(stdout);
#endif
   if ( recv_cnt > 0 )
   {
      ML_memory_free( (void **) &proc_flag );
      ML_memory_free( (void **) &Request );
   }
   if ( nvertices > 0 ) ML_memory_free( (void **) &pref_list );
   ML_free(in_preflist);
#ifdef RANDMISROOTS
   ML_free(neworder);
#endif
   return nselected;
}

/* ******************************************************************** */
/* A subroutine to update vertex states                                 */
/* -------------------------------------------------------------------- */

int ML_Aggregate_UpdateVertexStates(int N_remaining_vertices, 
        char vertex_state[], int recv_cnt, int recv_proc[], int recv_leng[], 
        int **recv_buf, int **recv_list, int proc_flag[], 
        int *NremainingRcvProcs, int send_cnt, int send_proc[], int send_leng[], 
        int **send_buf, int *send_flag, USR_REQ *Request, ML_Comm *comm, 
        int msgtype) 
{
   int    j, k, kkk, fproc;
   unsigned int nbytes;

      /* update the states to/from other processors */

      msgtype += 131;
      for ( j = 0; j < recv_cnt; j++ )
      {
         if ( proc_flag[j] == 0 )
         {
            fproc = recv_proc[j];
            nbytes = (recv_leng[j] + 1) * sizeof( int );
            comm->USR_irecvbytes((char*) recv_buf[j], nbytes, &fproc,
#ifdef ML_CPP
                    &msgtype, comm->USR_comm, &Request[j] );
#else
                    &msgtype, comm->USR_comm, (void *) &Request[j] );
#endif
         }
      }
      if ( *send_flag == 0 ) {
         for ( j = 0; j < send_cnt; j++ ) {
            nbytes = (send_leng[j] + 1) * sizeof( int );
#ifdef ML_AGGR_DEBUG
printf("%d: N_remaining_vertices is %d  | %d\n",comm->ML_mypid, 
        N_remaining_vertices,send_leng[j]); fflush(stdout);
printf("%d: a couple %d %d \n",
comm->ML_mypid, send_buf[j][0], send_buf[j][1]);
#endif
            if ( N_remaining_vertices <= 0 ) { 
               send_buf[j][send_leng[j]] = 1; 
               *send_flag = 1;
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
            comm->USR_cheapwaitbytes((char*) recv_buf[j], nbytes, &fproc,
#ifdef ML_CPP
                     &msgtype, comm->USR_comm, &Request[j] );
#else
                     &msgtype, comm->USR_comm, (void *) &Request[j] );
#endif
#ifdef ML_AGGR_DEBUG
printf("%d: rcvd is %d %d (%d)\n",comm->ML_mypid, recv_buf[0][0],fproc,
        recv_leng[j]+1); fflush(stdout);
#endif
            for ( k = 0; k < recv_leng[j]; k++ )
            {
               kkk = recv_list[j][k];
               if      (recv_buf[j][k] == 1)
               {
                  vertex_state[kkk] = 'S';
#ifdef ML_AGGR_DEBUG
                  printf("%d : incoming state %d = S \n",comm->ML_mypid,kkk);
#endif
               }
               else if (recv_buf[j][k] == 2)
               {
                  vertex_state[kkk] = 'D';
#ifdef ML_AGGR_DEBUG
                  printf("%d : incoming state %d = D \n",comm->ML_mypid,kkk);
#endif
               }
            }
#ifdef ML_AGGR_DEBUG
printf("%d: what are these %d %d\n",
comm->ML_mypid, recv_leng[j],recv_buf[j][recv_leng[j]]); fflush(stdout);
#endif
            if ( recv_buf[j][recv_leng[j]] == 1 )
            {
               proc_flag[j] = 1;
               (*NremainingRcvProcs)--;
#ifdef ML_AGGR_DEBUG
printf("the new NremainingRcvProcs is %d\n",*NremainingRcvProcs); fflush(stdout);
#endif
            }
         }
      }
      return 0;
}


int ML_Aggregate_Phase2_3_Cleanup(ML_Aggregate *ml_ag, ML_Operator *Amatrix,
				  int *aggr_count, int nvertices, 
				  int *aggr_index, int exp_Nrows, 
				  ML_Comm *comm, char *input_bdry, char *label,
                                  ML_agg_indx_comm *agg_indx_comm)
{
  int Nphase1_agg, phase_one_aggregated, i, j, aggd_neighbors, nonaggd_neighbors;
  int current_agg, total_aggs, *agg_incremented = NULL, *connect_type = NULL, 
    *number_connections = NULL;
  int best_score, score, best_agg, best_connect, mypid;
  double printflag;
  int                   allocated = 0, *rowi_col = NULL, rowi_N,total_vertices;
  double                *rowi_val = NULL, factor = 1.;
  char *bdry = NULL;
  int send_count, k, kk, flag, *temp_aggr_index;
  int Nleftover = 0, Nsingle = 0;

  if (input_bdry == NULL) {
   bdry = (char *) ML_allocate(sizeof(char)*(exp_Nrows + 1));
   for (i = nvertices ; i < exp_Nrows; i++) input_bdry[i] = 'F';
   for (i = 0; i < nvertices; i++) {
      bdry[i] = 'T';
      ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
                        &rowi_N, 0);
      if (rowi_N > 1) bdry[i] = 'F';
   }
  }
  else bdry = input_bdry;

    
   mypid       = comm->ML_mypid;
   printflag   = ml_ag->print_flag;
   Nphase1_agg = *aggr_count;

   /* Among unaggregated points, see if we can make a reasonable size    */
   /* aggregate out of it. We do this by looking at neighbors and seeing */
   /* how many are unaggregated and on my processor.                     */
   /*      NOTE: these aggregates must lie entirely on a processor.      */
   /*            The code/coordination between processors is just too    */
   /*            hard to work out the details.                           */

   phase_one_aggregated = 0;
   for (i = 0; i < nvertices; i++) {
     if ((aggr_index[i] != -1) && (bdry[i] == 'F')) {
        phase_one_aggregated++;
     }
   }
   phase_one_aggregated = ML_Comm_GsumInt( comm, phase_one_aggregated);
   total_vertices       = ML_Comm_GsumInt( comm, nvertices);

   /* base the number of new aggregates created on the percentage of */
   /* unaggregated nodes.                                            */

   factor = ((double) phase_one_aggregated)/((double)(total_vertices + 1));
   factor = pow(factor, ml_ag->phase3_agg_creation);

   for (i = 0; i < nvertices; i++) {
     if ((aggr_index[i] == -1) && (bdry[i] != 'T')) {
       ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val, &rowi_N, 0);
       aggd_neighbors = 0;
       nonaggd_neighbors = 0;
       for (j = 0; j < rowi_N; j++) {
	 current_agg = aggr_index[rowi_col[j]];
	 if (current_agg != -1) aggd_neighbors++;
	 else if (rowi_col[j] < nvertices) nonaggd_neighbors++;
       }
       if ((rowi_N > 3) && 
	   (((double) nonaggd_neighbors)/((double) rowi_N) > factor)) {
	 aggr_index[i] = (*aggr_count)++;
	 for (j = 0; j < rowi_N; j++) {
	   if ( (aggr_index[ rowi_col[j]] == -1) && (rowi_col[j] < nvertices))
	     aggr_index[rowi_col[j]] = aggr_index[i];
	 }
       }
     }
   }


   if ( printflag < ML_Get_PrintLevel()) {

     total_aggs    = ML_Comm_GsumInt( comm, Nphase1_agg);
     if (mypid == 0) {
       printf("Aggregation(%s) : Phase 1 - nodes aggregated = %d \n",label,
             phase_one_aggregated);
       printf("Aggregation(%s) : Phase 1 - total aggregates = %d\n",label, total_aggs);
     }
     i = *aggr_count - Nphase1_agg;
     i = ML_Comm_GsumInt( comm, i);
     if ( mypid == 0 ) {
       printf("Aggregation(%s) : Phase 2a- additional aggregates = %d\n",label, i);
     }
   }

   agg_incremented = (int *) ML_allocate(sizeof(int)* (exp_Nrows+1));
                                /* This guarantees that array will be long */
                                /* enough even for newly created aggregates. */
   for (i = 0; i < exp_Nrows; i++) agg_incremented[i] = 0;
   connect_type = (int *) ML_allocate(sizeof(int)* (exp_Nrows+1));
   for (i = 0; i < exp_Nrows; i++) connect_type[i] = 100;
   number_connections = (int *) ML_allocate(sizeof(int)*(exp_Nrows+1));
   for (i = 0; i < exp_Nrows; i++) number_connections[i] = 0;


   /* Try to stick unaggregated nodes into a neighboring aggegrate (the */
   /* smallest on processor) if they are not already too big. Otherwise */
   /* make a new aggregate.                                             */

   for (kk = 0; kk < 2; kk++) {
   for (i = 0; i < nvertices; i++) {
     if ((aggr_index[i] == -1) && (bdry[i] != 'T')) {
       ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
			 &rowi_N, 0);
       for (j = 0; j < rowi_N; j++) {
	 current_agg = aggr_index[rowi_col[j]];

	 /*********** do we neeeeed the extra condition????????? */
	 if ((current_agg >= 0) && /*(current_agg < Nphase1_agg) &&*/ (rowi_col[j] < nvertices)) {
	   number_connections[current_agg] += connect_type[rowi_col[j]];
	 }
       }
       best_score = -1000000;
       best_agg = -1; best_connect = -1;
       for (j = 0; j < rowi_N; j++) {
	 current_agg = aggr_index[rowi_col[j]];
 if (( current_agg >= 0) && /* (current_agg < Nphase1_agg) && do we need this???? */
	     (number_connections[current_agg] != 0) ) {
	   score = number_connections[current_agg]-agg_incremented[current_agg];
	   if (score > best_score) {
	     best_agg = current_agg; 
	     best_score = score;
	     best_connect = connect_type[rowi_col[j]];
	   }
	   else if ((best_agg == current_agg) && (connect_type[rowi_col[j]] > best_connect))
	     best_connect = connect_type[rowi_col[j]];

	   number_connections[current_agg] = 0;
	 }
       }
       if (best_score >= 0) { 
	 aggr_index[i] = best_agg;
	 agg_incremented[best_agg]++;
	 connect_type[i] = best_connect - 10;
       }
     }
   }
   }

   /* only MIS can have aggregates that span processors */

   if (agg_indx_comm != NULL) {
     ML_aggr_index_communicate(agg_indx_comm->N_neighbors, agg_indx_comm->temp_leng, 
			     agg_indx_comm->send_leng, agg_indx_comm->recv_leng, 
			     agg_indx_comm->send_list, agg_indx_comm->recv_list, 
			     nvertices, comm, aggr_index, 1573, agg_indx_comm->tem2_index, 
			     agg_indx_comm->neighbors, agg_indx_comm->temp_index,1);


     /* Find unaggregated points that I own. Assign them to an aggregate defined */
     /* on another processor. The main difficulty is that we do not have the     */
     /* aggregate index on the other processor .... so .... set temp_aggr_index  */
     /* to aggr_index[] and communicate without stuffing in the '-100-mypid'     */
     /* stuff. Now, look at neighbors of each unaggregated point, i,  until we   */
     /* find one assigned to an aggregate on another processor (we must also     */
     /* check that the ith point is a ghost node of the other processor). Set    */
     /* aggr_index[i] = aggr_index[neighbor] and communicate. This gives all     */
     /* processors who do not own the aggregate the correct '-100-mypid value.Set*/
     /* temp_aggr_index[i] = temp_aggr_index[neighbor] and communicate. This puts*/
     /* the correct aggregate number in the temp_aggr_index. Find all k such that*/
     /* aggr_index[k] = -100 - mypid and do aggr_index[k] = temp_aggr_index[k]   */

     temp_aggr_index = (int *) ML_allocate(sizeof(int)* (exp_Nrows+1));
     for (i = 0; i < exp_Nrows; i++) {
       if (aggr_index[i] >= 0) temp_aggr_index[i] = aggr_index[i];
       else temp_aggr_index[i] = -1;
     }
     /* communicate 'temp_aggr_index' twice to make sure that it gets there */
     /* even if stuff is communicated via an intermediate processor.        */

     ML_aggr_index_communicate(agg_indx_comm->N_neighbors, agg_indx_comm->temp_leng, 
			     agg_indx_comm->send_leng, agg_indx_comm->recv_leng, 
			     agg_indx_comm->send_list, agg_indx_comm->recv_list, 
			     nvertices, comm, 
			     temp_aggr_index, 1593, 
                             agg_indx_comm->tem2_index, 
			     agg_indx_comm->neighbors, agg_indx_comm->temp_index,0);

     ML_aggr_index_communicate(agg_indx_comm->N_neighbors, agg_indx_comm->temp_leng, 
			     agg_indx_comm->send_leng, agg_indx_comm->recv_leng, 
			     agg_indx_comm->send_list, agg_indx_comm->recv_list, 
			     nvertices, comm, 
			     temp_aggr_index, 1594, 
                             agg_indx_comm->tem2_index, 
			     agg_indx_comm->neighbors, agg_indx_comm->temp_index,0);


     for (i = 0; i < nvertices; i++) {
       if ((aggr_index[i] == -1) && (bdry[i] != 'T')) {
#ifdef CLEAN_DEBUG
	 printf("%d: Trying to assign %d to someone\n",comm->ML_mypid,i);
#endif
	 ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
			   &rowi_N, 0);
	 flag = 1;
	 for (j = 0; j < rowi_N; j++) {
	   current_agg = aggr_index[rowi_col[j]];
#ifdef CLEAN_DEBUG
	   printf("%d: a neighbor %d %d | %d\n",
		  comm->ML_mypid,rowi_col[j],current_agg,temp_aggr_index[rowi_col[j]]);
#endif

	   /* This node needs to be a ghost node of the processor that */
	   /* owns the aggregate!!!!!                                  */
	   if  ((flag==1)&&(current_agg < -1) && (temp_aggr_index[rowi_col[j]] >= 0)) {
	     current_agg = -100 - current_agg;
#ifdef CLEAN_DEBUG
	     printf("%d: lets go for it %d | proc = %d \n",comm->ML_mypid,i,current_agg); fflush(stdout);
#endif

	     /* Do I send this node to the processor 'current_agg'  */
	     send_count = 0;
	     for ( k = 0; k < agg_indx_comm->N_neighbors; k++ ) {
#ifdef CLEAN_DEBUG
	       printf("%d: my neighbor %d, length %d\n",comm->ML_mypid,
		      agg_indx_comm->neighbors[k],
		      agg_indx_comm->send_leng[k]);
#endif

	       if (agg_indx_comm->neighbors[k] == current_agg) {
		 for (kk = 0; kk < agg_indx_comm->send_leng[k]; kk++) {
#ifdef CLEAN_DEBUG
	       printf("%d: send %d %d\n",comm->ML_mypid,
		      agg_indx_comm->send_list[send_count+kk],send_count);
#endif
		   if (agg_indx_comm->send_list[send_count+kk] == i) {
#ifdef CLEAN_DEBUG
		     printf("%d: yes this is in %d on %d\n",
			    comm->ML_mypid,temp_aggr_index[rowi_col[j]],current_agg);
#endif
		     flag = 0;
		     aggr_index[i] = aggr_index[rowi_col[j]];
		     temp_aggr_index[i] = -10 - temp_aggr_index[rowi_col[j]];
		   }
		 }
	       }
	       send_count += agg_indx_comm->send_leng[k];
	     }
	   }
	 }
       }
     }

     /* communicate twice 'aggr_index' and 'temp_aggr_index' to make sure that */
     /* everyone gets everything. I think this is needed because some data     */
     /* might be communicated through another processor.                       */

     ML_aggr_index_communicate(agg_indx_comm->N_neighbors, agg_indx_comm->temp_leng, 
			     agg_indx_comm->send_leng, agg_indx_comm->recv_leng, 
			     agg_indx_comm->send_list, agg_indx_comm->recv_list, 
			     nvertices, comm, 
			     aggr_index, 
			     1578, agg_indx_comm->tem2_index, 
			     agg_indx_comm->neighbors, agg_indx_comm->temp_index,1);
     ML_aggr_index_communicate(agg_indx_comm->N_neighbors, agg_indx_comm->temp_leng, 
			     agg_indx_comm->send_leng, agg_indx_comm->recv_leng, 
			     agg_indx_comm->send_list, agg_indx_comm->recv_list, 
			     nvertices, comm, 
			     temp_aggr_index, 
			     1579, agg_indx_comm->tem2_index, 
			     agg_indx_comm->neighbors, agg_indx_comm->temp_index,0);
     ML_aggr_index_communicate(agg_indx_comm->N_neighbors, agg_indx_comm->temp_leng, 
			     agg_indx_comm->send_leng, agg_indx_comm->recv_leng, 
			     agg_indx_comm->send_list, agg_indx_comm->recv_list, 
			     nvertices, comm, 
			     aggr_index, 
			     1580, agg_indx_comm->tem2_index, 
			     agg_indx_comm->neighbors, agg_indx_comm->temp_index,1);
     ML_aggr_index_communicate(agg_indx_comm->N_neighbors, agg_indx_comm->temp_leng, 
			     agg_indx_comm->send_leng, agg_indx_comm->recv_leng, 
			     agg_indx_comm->send_list, agg_indx_comm->recv_list, 
			     nvertices, comm, 
			     temp_aggr_index, 
			     1581, agg_indx_comm->tem2_index, 
			     agg_indx_comm->neighbors, agg_indx_comm->temp_index,0);

     for (i= nvertices; i < exp_Nrows; i++) {
       if (aggr_index[i] == -100 - comm->ML_mypid) {
	 aggr_index[i] = -10 - temp_aggr_index[i];
#ifdef CLEAN_DEBUG
	 printf("%d: assigning %d to %d\n", comm->ML_mypid,i,aggr_index[i]);
#endif
       }
     }
     ML_free(temp_aggr_index);
   }

   for (i = 0; i < nvertices; i++) { 
     if ((aggr_index[i] == -1) && (bdry[i] != 'T')) {
       Nleftover++;
#ifdef CLEAN_DEBUG
       printf("%d: this guy is left %d    %d %d\n",comm->ML_mypid,i,nvertices,Nphase1_agg);
#endif
       ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
			 &rowi_N, 0);
	 /* We don't want a singleton. So lets see if there is an unaggregated */
	 /* neighbor that we can also but with this point.                     */

       flag = 1;
	 for (j = 0; j < rowi_N; j++) {
	   if ((rowi_col[j] != i) && (rowi_col[j] < nvertices)) {
	     current_agg = aggr_index[rowi_col[j]];
	     if ( current_agg == -1) {
	       flag = 0;
#ifdef CLEAN_DEBUG
	       printf("%d:      neighboring aggs %d %d\n",comm->ML_mypid,current_agg,rowi_col[j]);
#endif
	       aggr_index[rowi_col[j]] = *aggr_count;
	       connect_type[rowi_col[j]] = 100;
	     }
	   }
	 }
	 if (flag == 0) {
	   connect_type[i] = 100;
	   aggr_index[i] = (*aggr_count)++;
	 }
	 else {

	   /* We don't want a singleton. So lets see if we can connect it to */
	   /* any higher phase aggregate.                                    */
	   for (j = 0; j < rowi_N; j++) {
	     current_agg = aggr_index[rowi_col[j]];
	     if (( current_agg >= 0) && (rowi_col[j] < nvertices) ) break;
	   }
	   if (j < rowi_N) {
	     aggr_index[i] = current_agg;
	     (agg_incremented[current_agg])++;
	     connect_type[i] = connect_type[rowi_col[j]] - 11;
	   }
	   else {
	     Nsingle++;
#ifdef CLEAN_DEBUG
	     printf("%d: oh well  %d\n",comm->ML_mypid,i);
#endif
	     /* Singletons are problematics when the nullspace dimension */
	     /* is larger than the number of pde equations!!! This is    */
	     /* due to the QR factorization that occurs later.           */
	     if ((ml_ag->nullspace_dim > ml_ag->num_PDE_eqns) && ((*aggr_count)>0)) {
	       aggr_index[i] = 0; connect_type[i] = 1;
	     }
	     else {
	       aggr_index[i] = (*aggr_count)++;
	       connect_type[i] = 100;
	     }
	   }
	 } 
     }
   }



   if (printflag < ML_Get_PrintLevel()) {
     Nsingle = ML_Comm_GsumInt( comm, Nsingle);
     Nleftover = ML_Comm_GsumInt( comm, Nleftover);

     total_aggs = ML_Comm_GsumInt( comm, *aggr_count);
     j = 0;
     for (i = 0; i < nvertices; i++) if (bdry[i] == 'T') j++;
     j = ML_Comm_GsumInt( comm, j);

     if ( mypid == 0 ) {
       printf("Aggregation(%s) : Phase 2 - total aggregates = %d\n",label, total_aggs);
       printf("Aggregation(%s) : Phase 2 - boundary nodes   = %d\n",label, j);
       printf("Aggregation(%s) : Phase 3 - leftovers = %d and singletons = %d\n",label ,Nleftover, Nsingle);
     }
   }
  if (input_bdry == NULL) ML_free(bdry);
  if (rowi_col != NULL) ML_free(rowi_col);
  if (rowi_val != NULL) ML_free(rowi_val);
  if (agg_incremented != NULL) ML_free(agg_incremented);
  if (connect_type != NULL) ML_free(connect_type);
  if (number_connections != NULL) ML_free(number_connections);
  return 0;
}


int ML_aggr_index_communicate(int N_neighbors, int temp_leng[], int send_leng[],
			   int recv_leng[], int send_list[], int recv_list[],
			   int nvertices,
			   ML_Comm *comm, int aggr_index[], int msgtype,
			   int tem2_index[], int neighbors[], int temp_index[],
			      int flag)
{
   int count = 0, send_count = 0, recv_count = 0, i, j;
   int index;

   for ( i = 0; i < N_neighbors; i++ ) {
      temp_leng[i] = send_leng[i] + recv_leng[i];
      for ( j = 0; j < send_leng[i]; j++ ) {
         index = send_list[send_count++];
         if (index > nvertices) {
	   printf("ML_agg_MIS%d: sending something in rcv list %d %d\n",
		  comm->ML_mypid, index,nvertices);
	   exit(1);
         }
         if (aggr_index[index] == -1) temp_index[count] = -1;
         else if (aggr_index[index] <= -100) temp_index[count]=aggr_index[index];
         else if (flag == 1) temp_index[count] = -100 - comm->ML_mypid;
	 else  temp_index[count]=aggr_index[index];
         count++;
      }

      for ( j = 0; j < recv_leng[i]; j++ ) {
         index = recv_list[recv_count++];
         if (aggr_index[index] == -1) temp_index[count] = -1;
         else if (aggr_index[index] <= -100) temp_index[count]=aggr_index[index];
         else if (flag == 1) temp_index[count] = -100 - comm->ML_mypid;
	 else  temp_index[count]=aggr_index[index];
         count++;
      }
   }
   ML_Aggregate_ExchangeData((char*) tem2_index, (char*) temp_index,
      N_neighbors, neighbors, temp_leng, temp_leng, msgtype, ML_INT, comm);
   count = 0;
   send_count = 0;
   recv_count = 0;
   for ( i = 0; i < N_neighbors; i++ ) {
      for ( j = 0; j < recv_leng[i]; j++ ) {
         index = recv_list[recv_count++];
         if ((aggr_index[index] == -1) &&(tem2_index[count] != -1))
	   aggr_index[index] = tem2_index[count];
         count++;
      }
      for ( j = 0; j < send_leng[i]; j++ ) {
         index = send_list[send_count++];
         if ((aggr_index[index] == -1) && (tem2_index[count] != -1))
	   aggr_index[index] = tem2_index[count];
         count++;
      }
   }
   return 1;
}
