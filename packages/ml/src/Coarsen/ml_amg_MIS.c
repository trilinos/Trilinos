/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person    */
/* and disclaimer.                                                           */
/* ************************************************************************* */

/* ************************************************************************* */
/* ************************************************************************* */
/* Functions to create AMG prolongators                                      */
/* ************************************************************************* */
/* Author        : Charles Tong (LLNL)                                       */
/* Date          : October, 2000                                             */
/* ************************************************************************* */
/* ************************************************************************* */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ml_operator.h"
#include "ml_mat_formats.h"
#include "ml_amg.h"
#include "ml_aggregate.h"
#include "ml_smoother.h"
#include "ml_lapack.h"
#include "ml_utils.h"

/* ************************************************************************* */
/* construct the prolongator using Mark Adam's MIS algorithm                 */
/* ------------------------------------------------------------------------- */

int ML_AMG_CoarsenMIS( ML_AMG *ml_amg, ML_Operator *Amatrix, 
                       ML_Operator **Pmatrix, ML_Comm *comm)
{
  unsigned int nbytes;
   int     i, j, k, m, offset, count, index, ind2, col, *vlist;
   int     num_PDE_eqns, Nrows, exp_Nrows, Ncoarse, exp_Ncoarse, total_nnz;
   int     *rowptr, *column, *new_ia=NULL, *new_ja=NULL, *short_list;
   int     *CF_array, *sort_array, sortleng, bitindex,intindex, short_leng;
   double  printflag;
   int     sizeint, logsizeint, offnrows, *offibuffer, *offmap;
   int     *offmap2, *offlengths, numCi, *Ci_array,/* *int_array,*/ *int_buf;
   int     allocated=0, *rowi_col=NULL, rowi_N, nnz_per_row, *sys_array;
   int     msgtype, mypid, *send_leng=NULL, *recv_leng=NULL;
   int     N_neighbors, *neighbors=NULL, *send_list=NULL, sys_unk_filter;
   int     *recv_list=NULL, total_recv_leng=0, total_send_leng=0, idiag;
   int     A_ntotal, A_Nneigh, *A_rcvleng=NULL, *A_sndleng=NULL, *sys_info;
   int     *A_neigh=NULL, Nghost, **A_sndbuf=NULL, **A_rcvbuf=NULL;
   int     **A_sendlist=NULL, **A_recvlist=NULL, **proclist, *templist, nbdry;
   int     new_Nsend, new_Nrecv, *new_send_leng, *new_recv_leng, A_nnz;
   int     *new_send_list, *new_send_neigh, *new_recv_neigh, tmpcnt;
   int     *short_size;
   double  *new_val=NULL, epsilon, *rowi_val=NULL, *darray, rowmax, diag;
   double  *values, *offdbuffer, *dsumCij, dsumCi, rowsum, dtemp;
   char    *vtype, *state, *bdry, *border_flag;
   ML_GetrowFunc *getrow_obj;
   ML_CommInfoOP *mat_comm;
   struct ML_CSR_MSRdata *csr_data;
   ML_Aggregate_Comm *aggr_comm;

   /* ============================================================= */
   /* get the machine information and matrix references             */
   /* ============================================================= */

   mypid          = comm->ML_mypid;
   epsilon        = ml_amg->threshold;
   num_PDE_eqns   = ml_amg->num_PDE_eqns;
   printflag      = ml_amg->print_flag;
   Nrows          = Amatrix->outvec_leng;
   sys_unk_filter = 0;
   mat_comm       = Amatrix->getrow->pre_comm;

   /* ============================================================= */
   /* if system AMG (unknown approach) is requested, communicate    */
   /* the degree of freedom information                             */
   /* ============================================================= */

   if (ml_amg->amg_scheme == ML_AMG_SYSTEM_UNKNOWN && num_PDE_eqns > 1) 
   {
      sys_unk_filter = 1;
      count    = Nrows + 1;
      if ( mat_comm != NULL) count += mat_comm->total_rcv_length;
      darray = (double *) ML_allocate(sizeof(double) * count);
      for (i = 0; i < Nrows; i++) darray[i] = (double) ml_amg->blk_info[i];
      if ( mat_comm != NULL )
         ML_exchange_bdry(darray,mat_comm,Nrows,comm,ML_OVERWRITE,NULL);
      sys_info = (int *) ML_allocate(sizeof(double) * count);
      for (i = 0; i < count; i++) sys_info[i] = (int) darray[i];
      ML_free(darray);
   } else sys_info = NULL;

   /* ============================================================= */
   /* check the system size                                         */
   /* ============================================================= */

   if ( Nrows % num_PDE_eqns != 0 )
   {
      printf("ML_AMG_CoarsenMIS ERROR : Nrows must be multiples");
      printf(" of num_PDE_eqns.\n");
      exit(1);
   }
   if ( mypid == 0 && printflag < ML_Get_PrintLevel())
   {
      printf("ML_AMG_CoarsenMIS : current level = %d\n", ml_amg->cur_level);
      printf("ML_AMG_CoarsenMIS : current eps   = %e\n", epsilon);
   }

   /* ============================================================= */
   /* fetch getrow/comm information of the Amat                     */
   /* ============================================================= */

   ML_AMG_GetCommInfo(mat_comm, Nrows, &A_Nneigh, &A_neigh, &A_sendlist, 
                      &A_recvlist, &A_sndbuf, &A_rcvbuf, &A_sndleng, 
                      &A_rcvleng, &Nghost);

   A_ntotal = Nrows + Nghost;

   /* ============================================================= */
   /* set up the proclist array for labeling                        */
   /* ============================================================= */
   /* templist[i] tabulates how many processors vertex i is to be   */
   /* sent to.  Then allocate proclist to record the processors and */
   /* indices each of my local vertices are to send.  The first     */
   /* element of the array proclist[i][0] is a counter of how many  */
   /* processors, followed by a number of processor (not processor  */
   /* number, but processor osfset) and index pairs.                */
   /* For the receive vertices, proclist[i][0] indicates which      */
   /* processor (offset) the data is coming from                    */
   /* ============================================================= */

   if ( Nrows > 0 ) templist = (int *) ML_allocate(sizeof(int)*Nrows);
   else             templist = NULL;
   for ( i = 0; i < Nrows; i++ ) templist[i] = 0;
   for ( i = 0; i < A_Nneigh; i++ ) 
   {
      for ( j = 0; j < A_sndleng[i]; j++ ) 
      {
         index = A_sendlist[i][j];
         if ( index >= Nrows || index < 0 ) 
         {
            printf("%d : SYSTEM ERROR (1) in ML_AMG_CoarsenMIS.\n", mypid);
            exit(0);
         }
         templist[index]++;
      }
   }
   if (A_ntotal > 0) proclist = (int **) ML_allocate(A_ntotal * sizeof(int *));
   else              proclist = NULL;
   for ( i = 0; i < Nrows; i++ ) 
   {
      if ( templist[i] > 0 )
      {
         proclist[i] = (int *) ML_allocate( (2*templist[i]+1) * sizeof( int ) );
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
   for ( i = Nrows; i < A_ntotal; i++ ) 
   {
      proclist[i] = (int *) ML_allocate( sizeof( int ) );
   }
   for ( i = 0; i < A_Nneigh; i++ ) {
      for ( j = 0; j < A_rcvleng[i]; j++ ) {
         index = A_recvlist[i][j];
         proclist[index][0] = A_neigh[i];
      }
   }
   if ( templist != NULL ) ML_free(templist);

   /* ============================================================= */
   /* record the Dirichlet boundary and processor boundaries        */
   /* ============================================================= */

   if (A_ntotal > 0) bdry = (char *) ML_allocate(sizeof(char)*(A_ntotal + 1));
   else              bdry = NULL;
   if ( Nrows > 0 ) border_flag = (char *) ML_allocate(sizeof(char)*Nrows);
   else             border_flag = NULL;
   total_nnz = 0;
   nbdry     = 0;
   for (i = 0; i < Nrows; i++) 
   {
      bdry[i] = 'T';
      border_flag[i] = 'F';
      ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
                        &rowi_N, 0);
      rowsum = diag = 0.0;
      count = 0;
      for (j = 0; j < rowi_N; j++) 
      {
         if ( rowi_col[j] == i ) diag = ML_dabs( rowi_val[j] );
         else                    rowsum += ML_dabs( rowi_val[j] );
         if ( rowi_col[j] >= Nrows ) border_flag[i] = 'T';
         if ( rowi_val[j] != 0.0 ) count++;
      }
      if (count > 1) bdry[i] = 'F';
      if ( diag != 0.0 )
      {
         if ( (rowsum / diag) < 0.0e-8 ) bdry[i] = 'T';
      }
      total_nnz += rowi_N;
      if ( bdry[i] == 'T' ) nbdry++;
   }
   m = ML_Comm_GsumInt( comm, nbdry );
   if ( mypid == 0 && printflag < ML_Get_PrintLevel())
      printf("AMG Phase 1  - total number of bndry points  = %6d \n",m); 

   /* ============================================================= */
   /* create the strength matrix in (rowptr, column)                */
   /* ============================================================= */

   A_nnz = total_nnz;
   rowptr = (int *) ML_allocate( (Nrows + 1)* sizeof(int) );
   if ( total_nnz > 0 ) 
   {
      column = (int *)    ML_allocate( total_nnz * sizeof(int) );
      values = (double *) ML_allocate( total_nnz * sizeof(double) );
   }
   else
   {
      column = NULL;
      values = NULL;
   }
   total_nnz = 0;
   rowptr[0] = 0;
   for (i = 0; i < Nrows; i++)
   {
      ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
                        &rowi_N, 0);
      if ( sys_unk_filter )
      {
         for (j = 0; j < rowi_N; j++) 
            if (sys_info[rowi_col[j]] != sys_info[i]) rowi_val[j] = 0.0; 
      } 
      diag = 0.0;
      for (j = 0; j < rowi_N; j++) 
         if ( rowi_col[j] == i ) diag = rowi_val[j];
      rowmax = 0.0;
      if ( diag >= 0. )
      {
         for (j = 0; j < rowi_N; j++) 
            if (rowi_col[j] != i) rowmax = ML_min(rowmax, rowi_val[j]); 
      }
      else
      {
         for (j = 0; j < rowi_N; j++) 
            if (rowi_col[j] != i) rowmax = ML_max(rowmax, rowi_val[j]); 
      }
      rowmax *= epsilon;
      if ( diag >= 0. )
      {
         for (j = 0; j < rowi_N; j++) 
         {
            if ( rowi_col[j] != i && rowi_val[j] != 0 && rowi_val[j] < rowmax ) 
            {
               values[total_nnz]   = rowi_val[j];
               column[total_nnz++] = rowi_col[j];
            }
         }
      }
      else
      {
         for (j = 0; j < rowi_N; j++) 
         {
            if ( rowi_col[j] != i && rowi_val[j] != 0 && rowi_val[j] > rowmax ) 
            {
               values[total_nnz]   = rowi_val[j];
               column[total_nnz++] = rowi_col[j];
            }
         }
      }
      rowptr[i+1] = total_nnz;
   }
   ML_free( rowi_col );
   ML_free( rowi_val );
   dtemp = A_nnz;
   dtemp = ML_gsum_double(dtemp, comm);
   if ( ml_amg->operator_complexity == 0.0 )
   {
      ml_amg->fine_complexity = dtemp;
      ml_amg->operator_complexity = dtemp;
   }
   else
   {
      ml_amg->operator_complexity += dtemp;
   }

   /* ============================================================= */
   /* communicate the boundary information                          */
   /* ============================================================= */

   darray = (double *) ML_allocate(sizeof(double)*(A_ntotal+1));
   for (i = 0; i < Nrows; i++) 
   {
      if (bdry[i] == 'T') darray[i] = 1.;
      else  darray[i] = 0.;
   }
   ML_exchange_bdry(darray,Amatrix->getrow->pre_comm,Nrows,comm,
                    ML_OVERWRITE,NULL);
   for (i = Nrows; i < A_ntotal; i++) 
   {
      if (darray[i] == 1.) bdry[i] = 'T';
      else bdry[i] = 'F';
   }
   ML_free(darray);

   /* ============================================================= */
   /* get ready for coarsening (allocate temporary arrays)          */
   /* ============================================================= */

   if ( A_ntotal > 0 ) CF_array = (int *) ML_allocate(sizeof(int)* A_ntotal);
   else                CF_array = NULL;
   for (i = 0; i < A_ntotal; i++) CF_array[i] = -1;

   if ( Nrows > 0 ) vlist = (int *) ML_allocate(sizeof(int)* Nrows);
   else             vlist = NULL;
   if ( A_ntotal > 0 ) state = (char *) ML_allocate(sizeof(char)* A_ntotal);
   else                state = NULL;
   if ( A_ntotal > 0 ) vtype = (char *) ML_allocate(sizeof(char)* A_ntotal);
   else                vtype = NULL;
   for (i = 0; i < Nrows   ; i++) vlist[i] = i;
   for (i = 0; i < A_ntotal; i++) state[i] = 'F';
   for (i = 0; i < A_ntotal; i++) vtype[i] = 'x';

   /* ============================================================= */
   /* Phase 1 : compute an initial MIS                              */
   /* delete nodes that are just isolated Dirichlet points (B)      */
   /* then label the vertices as selected(S) or deleted (D)         */
   /* ============================================================= */

   m = ML_Comm_GsumInt( comm, Nrows );
   k = ML_Comm_GsumInt( comm, A_nnz );
   if ( mypid == 0 && printflag < ML_Get_PrintLevel())
   {
      printf("AMG Phase 1 begins, total_nrows, total_nnz   = %d %d\n", m, k); 
      fflush(stdout);
   }
   for (i = 0; i < A_ntotal ; i++) if (bdry[i] == 'T') state[i] = 'B'; 

   k = ML_AMG_LabelVertices(Nrows, vlist, 'x', state, vtype, 
                      Nrows, rowptr, column, mypid, proclist, 
                      A_Nneigh,A_sndbuf,A_neigh, A_sndleng, A_Nneigh,
                      A_rcvbuf, A_neigh, A_rcvleng, A_recvlist, 
                      comm, CF_array);

   Ncoarse = 0;
   for (i = 0; i < Nrows ; i++) if ( state[i] == 'S' ) Ncoarse++;
   m = ML_Comm_GsumInt( comm, Ncoarse );
   if ( mypid == 0 && printflag < ML_Get_PrintLevel())
      printf("AMG Phase 1  - total number of coarse points = %6d (%3d)\n",m,k); 
   if ( printflag < ML_Get_PrintLevel() )
   {
      printf("%4d : Phase 1 - number of coarse points = %6d\n",mypid,Ncoarse); 
      fflush(stdout);
   }

   /* ============================================================= */
   /* Phase 2 : if any F point has a neighbor F point strongly      */
   /*           coupled to each other, make sure they share at      */
   /*           least one common C point                            */
   /*   This phase consists of 2 passes :                           */
   /*   (1) take care of the border nodes first                     */
   /*   (2) take care of the interior nodes                         */
   /* ============================================================= */

   /* ============================================================= */
   /* that rows for neighbors that are in other processors          */
   /* ------------------------------------------------------------- */

   ML_Smoother_ComposeOverlappedMatrix(Amatrix, comm, &offnrows,
              &offlengths, &offibuffer, &offdbuffer, &offmap, &offmap2, &k);
   offset = 0;
   exp_Nrows = Nrows + offnrows;
   for ( i = 0; i < offnrows; i++ )
   {
      diag = 0.0;
      for ( j = offset; j < offset+offlengths[i]; j++ )
      {
         index = offibuffer[j];
         if ( index >= k && index < k+Nrows ) offibuffer[j] = index - k;
         else
         {
            count = ML_sorted_search(index, exp_Nrows-Nrows, offmap);
            if ( count >= 0 ) offibuffer[j] = offmap2[count] + Nrows;
            else              offibuffer[j] = -1;
         }
         if ( index == (i+Nrows) ) diag = offdbuffer[j];
      }
      if ( diag >= 0.0 ) rowmax =  1.0e20;
      else               rowmax = -1.0e20;
      for ( j = offset; j < offset+offlengths[i]; j++ )
      {
         index = offibuffer[j];
         if ( index != (i+Nrows) )
         {
            if ( diag >= 0. ) rowmax = ML_min(rowmax, offdbuffer[j]); 
            else              rowmax = ML_max(rowmax, offdbuffer[j]); 
         }
      }
      rowmax *= epsilon;
      for ( j = offset; j < offset+offlengths[i]; j++ )
      {
         index = offibuffer[j];
         if ( index != (i+Nrows) )
         {
            if ( diag >= 0. )
            {
               if ( offdbuffer[j] != 0. && offdbuffer[j] >= rowmax ) 
                  offibuffer[j] = - index - 2;
            }
            else
            {
               if ( offdbuffer[j] != 0. && offdbuffer[j] <= rowmax ) 
                  offibuffer[j] = - index - 2;
            }
         }
         else offibuffer[j] = - index - 2;
      }
      offset += offlengths[i];
   }
   for ( i = 1; i < offnrows; i++ ) offlengths[i] += offlengths[i-1];
   if ( offmap  != NULL ) ML_free( offmap );
   if ( offmap2 != NULL ) ML_free( offmap2 );

   /* ============================================================= */
   /* communicate of offprocessor C/F information                   */
   /* ------------------------------------------------------------- */

   total_send_leng = 0;
   for (i = 0; i < A_Nneigh ; i++) total_send_leng += A_sndleng[i]; 
   nbytes = total_send_leng * sizeof(int);
   if ( nbytes > 0 ) int_buf = (int *) ML_allocate( nbytes );
   else              int_buf = NULL;
   offset = 0;
   for ( i = 0; i < A_Nneigh; i++ ) 
   {
      for ( j = 0; j < A_sndleng[i]; j++ ) 
         int_buf[offset+j] = CF_array[A_sendlist[i][j]];
      offset += A_sndleng[i];
   }
   msgtype = 35733;
   ML_Aggregate_ExchangeData((char*) &(CF_array[Nrows]), (char*) int_buf,
      A_Nneigh, A_neigh, A_rcvleng, A_sndleng, msgtype, ML_INT, comm);

   if ( int_buf != NULL ) ML_free(int_buf);
/*
for ( i = 0; i < Nrows; i++ )
   if ( CF_array[i] >= 0 ) printf("%d : C point = %d\n",mypid,i);
*/

   /* ============================================================= */
   /* get ready to examine all local F points                       */
   /* ------------------------------------------------------------- */

   sortleng  = A_ntotal / (8 * sizeof(int)) + 1;
   sort_array = (int *) ML_allocate( sortleng * sizeof(int) );
   for ( i = 0; i < sortleng; i++ ) sort_array[i] = 0;
   sizeint = sizeof(int) * 8;
   if ( sizeint == 16 )      logsizeint = 4;
   else if ( sizeint == 32 ) logsizeint = 5;
   else if ( sizeint == 64 ) logsizeint = 6;
   else                      logsizeint = 5;

   /* ============================================================= */
   /* Pass one : take care of the border nodes first                */
   /* ============================================================= */
   /* search for common C point between two F points                */
   /* five cases to deal with :                                     */
   /* 1. F1 on P0, F2 on P0, C on P0                                */
   /* 2. F1 on P0, F2 on P0, C on P1                                */
   /* 3. F1 on P0, F2 on P1, C on P0                                */
   /* 4. F1 on P0, F2 on P1, C on P1                                */
   /* 5. F1 on P0, F2 on P1, C on P2                                */
   /* ------------------------------------------------------------- */

   for ( i = 0; i < Nrows; i++ )
   {
      if (state[i] != 'B' && CF_array[i] < 0 && border_flag[i] == 'T')
      { /* -- border F -- */
         /* ----- register my C neighbors ----- */

         for ( j = rowptr[i]; j < rowptr[i+1] ; j++) 
         {
            col = column[j];
            if ( col != i && CF_array[col] >= 0 ) 
            {
               intindex = col >> logsizeint;
               bitindex = col % sizeint;
               sort_array[intindex] |= ( 1 << bitindex );
            }
         }

         /* ----- examine my strongly connected neighbors ----- */

         for (j = rowptr[i]; j < rowptr[i+1] ; j++) 
         {
            col = column[j];
            if ( col != i && CF_array[col] < 0 ) /* --- F-F --- */
            {
               if ( col < Nrows )  /* --- case 1 and 2 --- */
               {
                  for (k = rowptr[col]; k < rowptr[col+1] ; k++) 
                  {
                     ind2 = column[k];
                     if ( ind2 != col && CF_array[ind2] >= 0 ) 
                     {
                        intindex = ind2 >> logsizeint;
                        bitindex = ind2 % sizeint;
                        if (sort_array[intindex] & (1 << bitindex)) break;
                     }
                  }
                  if ( k == rowptr[col+1] ) /* --- not found --- */
                  {
                     CF_array[i] = Ncoarse++;
                     break;
                  }
               } 
               else /* --- case 3, 4, and 5 --- */
               {
                  if ( (col-Nrows) == 0 ) index = 0;
                  else                    index = offlengths[col-Nrows-1];

                  for (k = index; k < offlengths[col-Nrows] ; k++) 
                  {
                     ind2 = offibuffer[k];
                     if (ind2 >= 0 && CF_array[ind2] >= 0) 
                     {
                        intindex = ind2 >> logsizeint;
                        bitindex = ind2 % sizeint;
                        if (sort_array[intindex] & (1 << bitindex)) break;
                     }
                  }
                  if ( k == offlengths[col-Nrows] ) /* --- not found --- */
                  {
                     /* first search to see if equation col has my    */
                     /* equation (i) as strong connection.  If not, I */
                     /* declare myself as a C point, else whoever is  */
                     /* in processor with rank is the C point         */

                     count = offlengths[col-Nrows];
                     ind2  = -1;
                     for ( m = index; m < count; m++ )
                        if ( offibuffer[m] == i ) { ind2 = m - index; break;}
                     if ( ind2 < 0 ) 
                     {
                        CF_array[i] = Ncoarse++; 
                        break;
                     }
                     else
                     {
                        count = Nrows;
                        for ( m = 0; m < A_Nneigh; m++ )
                           if ( col < count ) count += A_rcvleng[m];

/* ############################################################### */
/* don't decide on coarse point on other processors        
                        if ( mypid < A_neigh[m-1] ) 
                        {
                           CF_array[col] = Ncoarse++;
                        }
                        else 
*/
/* ############################################################### */
/* taking this if statement out improves convergence dramatically,
   but also increase the operator complexity
*/
                        if ( mypid > A_neigh[m-1] ) 
/* ############################################################### */
                        {
                           CF_array[i] = Ncoarse++; 
                           break;
                        }
                     }
                  }
               } 
            }  /* if a neighbor of F point i is also a F point */
         }  /* for all neighbors of an F point i */ 

         /* ----- reset the C neighbor registers ----- */

         for ( j = rowptr[i]; j < rowptr[i+1] ; j++ ) 
         {
            col = column[j];
            if ( col != i && CF_array[col] >= 0 ) 
            {
               intindex = col >> logsizeint;
               sort_array[intindex] = 0;
            }
         }
      } /* if CF_array[i] < 0 - i a F point */ 
   } 
   m = ML_Comm_GsumInt( comm, Ncoarse );
   if ( mypid == 0 && printflag < ML_Get_PrintLevel())
   {
      printf("AMG Phase 2a - total number of coarse points = %6d\n", m); 
   }

   /* ============================================================= */
   /* communicate of offprocessor C/F information                   */
   /* ------------------------------------------------------------- */

   nbytes = total_send_leng * sizeof(int);
   if ( nbytes > 0 ) int_buf = (int *) ML_allocate( nbytes );
   else              int_buf = NULL;
   offset = 0;
   for ( i = 0; i < A_Nneigh; i++ ) 
   {
      for ( j = 0; j < A_sndleng[i]; j++ ) 
         int_buf[offset+j] = CF_array[A_sendlist[i][j]];
      offset += A_sndleng[i];
   }
   msgtype = 35734;
   ML_Aggregate_ExchangeData((char*) &(CF_array[Nrows]), (char*) int_buf,
      A_Nneigh, A_neigh, A_rcvleng, A_sndleng, msgtype, ML_INT, comm);

   if ( int_buf != NULL ) ML_free(int_buf);

   /* ============================================================= */
   /* Pass 2 : handle all interior nodes                            */
   /* ============================================================= */
   /* search for common C point between two F points                */
   /* only 1 case to deal with :                                    */
   /* 1. F1 on P0, F2 on P0, C on P0                                */
   /* ------------------------------------------------------------- */

   short_leng = 0;
   short_list = (int *) ML_allocate( Nrows * sizeof(int) );
   short_size = (int *) ML_allocate( Nrows * sizeof(int) );
   for ( i = 0; i < Nrows; i++ )
   {
      if (state[i] != 'B' && CF_array[i] < 0 && border_flag[i] == 'F') 
      {
         short_list[short_leng] = i;
         short_size[short_leng] = 0;
         for ( j = rowptr[i]; j < rowptr[i+1]; j++ )
            if ( CF_array[column[j]] < 0 ) short_size[short_leng]++;
         short_leng++;
      }
   }
   ML_az_sort( short_size, short_leng, short_list, NULL );
   for ( i = 0; i < short_leng/2; i++ )
   {
      j = short_list[i];
      short_list[i] = short_list[short_leng-1-i];
      short_list[short_leng-1-i] = j;
      j = short_size[i];
      short_size[i] = short_size[short_leng-1-i];
      short_size[short_leng-1-i] = j;
   }
   if ( short_leng > 0 ) index = short_list[0];
   if ( short_leng > 0 ) j     = short_size[0];
   count = 0;
   for ( i = 1; i < short_leng; i++ )
   {
      if ( short_size[i] != j )
      {
         ML_az_sort( &short_list[count], i-count, NULL, NULL );
         index = short_list[i];
         j     = short_size[i];
         count = i; 
      }
   }
   for ( i = 0; i < short_leng; i++ )
   {
      index = short_list[i];
      for ( j = rowptr[index]; j < rowptr[index+1] ; j++) 
      {
         col = column[j];
         if ( col != index && CF_array[col] >= 0 ) 
         {
            intindex = col >> logsizeint;
            bitindex = col % sizeint;
            sort_array[intindex] |= ( 1 << bitindex );
         }
      }
      for (j = rowptr[index]; j < rowptr[index+1] ; j++) 
      {
         col = column[j];
         if ( col != index && CF_array[col] < 0 ) /* --- F-F --- */
         {
            for (k = rowptr[col]; k < rowptr[col+1] ; k++) 
            {
               ind2 = column[k];
               if ( ind2 != col && CF_array[ind2] >= 0 ) 
               {
                  intindex = ind2 >> logsizeint;
                  bitindex = ind2 % sizeint;
                  if (sort_array[intindex] & (1 << bitindex)) break;
               }
            }
            if ( k == rowptr[col+1] ) /* --- shared C not found --- */
            {
               CF_array[index] = Ncoarse++;
               break;
            }
         }
      }
      for ( j = rowptr[index]; j < rowptr[index+1] ; j++ ) 
      {
         col = column[j];
         if ( col != index && CF_array[col] >= 0 ) 
         {
            intindex = col >> logsizeint;
            sort_array[intindex] = 0;
         }
      }
   } 

/* #ifdef ML_DEBUG_AMG */
   for ( i = 0; i < short_leng; i++ )
   {
      index = short_list[i];
      if ( CF_array[index] < 0 )
      {
         for ( j = rowptr[index]; j < rowptr[index+1] ; j++) 
         {
            col = column[j];
            if ( col != index && CF_array[col] >= 0 ) 
            {
               intindex = col >> logsizeint;
               bitindex = col % sizeint;
               sort_array[intindex] |= ( 1 << bitindex );
            }
         }
         for (j = rowptr[index]; j < rowptr[index+1] ; j++) 
         {
            col = column[j];
            if ( col != index && CF_array[col] < 0 ) /* --- F-F --- */
            {
               for (k = rowptr[col]; k < rowptr[col+1] ; k++) 
               {
                  ind2 = column[k];
                  if ( ind2 != col && CF_array[ind2] >= 0 ) 
                  {
                     intindex = ind2 >> logsizeint;
                     bitindex = ind2 % sizeint;
                     if (sort_array[intindex] & (1 << bitindex)) break;
                  }
               }
               if ( k == rowptr[col+1] ) /* --- shared C not found --- */
               {
                  printf("%d : AMG Rule C1 violated.\n", mypid);
                  break;
               }
            }
         }
         for ( j = rowptr[index]; j < rowptr[index+1] ; j++ ) 
         {
            col = column[j];
            if ( col != index && CF_array[col] >= 0 ) 
            {
               intindex = col >> logsizeint;
               sort_array[intindex] = 0;
            }
         }
      } 
   } 
/* #endif */

   ML_free( short_list );
   ML_free( short_size );

#ifdef obsolete
   /* This is obsolete because I want to give preference to vertices */
   /* with the highest number of F neighbors                         */
   for ( i = 0; i < Nrows; i++ )
   {
      if (state[i] != 'B' && CF_array[i] < 0 && border_flag[i] == 'F') 
      {  /* -- interior F -- */
         /* ----- register my C neighbors ----- */

         for ( j = rowptr[i]; j < rowptr[i+1] ; j++) 
         {
            col = column[j];
            if ( col != i && CF_array[col] >= 0 ) 
            {
               intindex = col >> logsizeint;
               bitindex = col % sizeint;
               sort_array[intindex] |= ( 1 << bitindex );
            }
         }

         /* ----- examine my strongly connected neighbors ----- */

         for (j = rowptr[i]; j < rowptr[i+1] ; j++) 
         {
            col = column[j];
            if ( col != i && CF_array[col] < 0 ) /* --- F-F --- */
            {
               for (k = rowptr[col]; k < rowptr[col+1] ; k++) 
               {
                  ind2 = column[k];
                  if ( ind2 != col && CF_array[ind2] >= 0 ) 
                  {
                     intindex = ind2 >> logsizeint;
                     bitindex = ind2 % sizeint;
                     if (sort_array[intindex] & (1 << bitindex)) break;
                  }
               }
               if ( k == rowptr[col+1] ) /* --- shared C not found --- */
               {
                  CF_array[i] = Ncoarse++;
                  break;
               }
            }
         }

         /* ----- reset the C neighbor registers ----- */

         for ( j = rowptr[i]; j < rowptr[i+1] ; j++ ) 
         {
            col = column[j];
            if ( col != i && CF_array[col] >= 0 ) 
            {
               intindex = col >> logsizeint;
               sort_array[intindex] = 0;
            }
         }
      }
   } 
#endif
   m = ML_Comm_GsumInt( comm, Ncoarse );
   if ( mypid == 0 && printflag < ML_Get_PrintLevel())
   {
      printf("AMG Phase 2b - total number of coarse points = %6d\n", m); 
   }
   if ( printflag < ML_Get_PrintLevel())
   {
      printf("%4d : Phase 2 - number of coarse points = %6d\n",mypid,Ncoarse); 
      fflush(stdout);
   }
/*
for ( i = 0; i < Nrows; i++ )
   if ( CF_array[i] >= 0 ) printf("%d : C point = %d\n",mypid,i);
*/

   /* ============================================================= */
   /* checking                                                      */
   /* ============================================================= */

   count = 0;
   for ( i = 0; i < Nrows; i++ )
   {
      m = 1;
      if ( bdry[i] == 'T' ) m = 0;
      if ( (rowptr[i+1]-rowptr[i]) == 0 ) m = 0;
      if ( m == 1 && CF_array[i] < 0 )
      {
         for ( j = rowptr[i]; j < rowptr[i+1]; j++ )
            if ( CF_array[column[j]] >= 0 ) {m = 0; break;}
      } else m = 0;
      count += m;
      if ( mypid == 0  && m == 1 && i == 0 ) 
      {
         printf("BAD ROW = %d (%d)\n", i, Nrows);
         allocated = 0;
         rowi_col = NULL;
         rowi_val = NULL;
         ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
                           &rowi_N, 0);
         for ( j = 0; j < rowi_N; j++ )
         {
            col = rowi_col[j];
            if ( col < Nrows )
            {
               for (k = rowptr[col]; k < rowptr[col+1] ; k++) 
                  printf("   COL,VAL = %7d %e\n", column[k], values[k]);
            }
            else
            {
               if ( (col-Nrows) == 0 ) ind2 = 0;
               else                    ind2 = offlengths[col-Nrows-1];
            }
         }
      }
   }
   k = ML_Comm_GsumInt( comm, count );
   if ( mypid == 0 && printflag < ML_Get_PrintLevel())
   {
      printf("AMG Coarsen  - total number of bad points    = %6d\n", k); 
   }
 
   /* ============================================================= */
   /* compatible relaxation                                         */
   /* ------------------------------------------------------------- */

/* Haven't proved it to be useful yet, but may continue to try
   ML_AMG_CompatibleRelaxation(CF_array,Amatrix,&Ncoarse,Ncoarse/10+1);
*/

   m = ML_Comm_GsumInt( comm, Ncoarse );
   if ( mypid == 0 && printflag < ML_Get_PrintLevel())
   {
      printf("AMG Phase 2c - total number of coarse points = %6d\n", m); 
   }

   /* ============================================================= */
   /* communicate of offprocessor C/F information                   */
   /* ------------------------------------------------------------- */

   if ( count > 0 )
   {
      for ( i = 0; i < Nrows; i++ )
      {
         m = 1;
         if ( bdry[i] == 'T' ) m = 0;
         if ( (rowptr[i+1]-rowptr[i]) == 0 ) m = 0;
         if ( m == 1 && CF_array[i] < 0 )
         {
            for ( j = rowptr[i]; j < rowptr[i+1]; j++ )
               if ( CF_array[column[j]] >= 0 ) {m = 0; break;}
         } else m = 0;
         if ( m == 1 ) CF_array[i] = Ncoarse++; 
      }
   }
   if ( k > 0 )
   {
      nbytes = total_send_leng * sizeof(int);
      if ( nbytes > 0 ) int_buf = (int *) ML_allocate( nbytes );
      else              int_buf = NULL;
      offset = 0;
      for ( i = 0; i < A_Nneigh; i++ ) 
      {
         for ( j = 0; j < A_sndleng[i]; j++ ) 
            int_buf[offset+j] = CF_array[A_sendlist[i][j]];
         offset += A_sndleng[i];
      }
      msgtype = 35736;
      ML_Aggregate_ExchangeData((char*) &(CF_array[Nrows]), (char*) int_buf,
         A_Nneigh, A_neigh, A_rcvleng, A_sndleng, msgtype, ML_INT, comm);
      if ( int_buf != NULL ) ML_free(int_buf);
      m = ML_Comm_GsumInt( comm, Ncoarse );
      if ( mypid == 0 && printflag < ML_Get_PrintLevel())
      {
         printf("AMG Phase 2d - total number of coarse points = %6d\n", m); 
      }
   }

   /* ============================================================= */
   /* register system information for the next level, as needed by  */
   /* the unknown approach.                                         */
   /* ------------------------------------------------------------- */

   if ( ml_amg->amg_scheme == ML_AMG_SYSTEM_UNKNOWN &&
        ml_amg->blk_info != NULL && num_PDE_eqns > 1 )
   {
      sys_array = ml_amg->blk_info;
      if ( Ncoarse > 0 ) 
      {
         nbytes = Ncoarse * sizeof(int);
         ML_memory_alloc((void**)&(ml_amg->blk_info), nbytes, "AM2");
      }
      else ml_amg->blk_info = NULL;
      m = 0;
      for ( i = 0; i < Nrows; i++ )
         if ( CF_array[i] >= 0 ) ml_amg->blk_info[m++] = sys_array[i];
      ML_memory_free((void**) &sys_array);
   }
       
   /* ============================================================= */
   /* clean up                                                      */
   /* ------------------------------------------------------------- */

   if ( sort_array != NULL ) ML_free( sort_array );

   /* ============================================================= */
   /* free memory used for doing the MIS stuff                      */
   /* ============================================================= */

   for ( i = 0; i < A_ntotal; i++ ) 
      if ( proclist[i] != NULL ) ML_free(proclist[i]);
   if ( proclist    != NULL ) ML_free(proclist);
   if ( vlist       != NULL ) ML_free(vlist); 
   if ( state       != NULL ) ML_free(state); 
   if ( vtype       != NULL ) ML_free(vtype);
   if ( bdry        != NULL ) ML_free( bdry );
   if ( border_flag != NULL ) ML_free(border_flag);
   for (i = 0; i < A_Nneigh; i++) 
   {
      ML_free(A_recvlist[i]);
      ML_free(A_sendlist[i]);
      ML_free(A_rcvbuf[i]);
      ML_free(A_sndbuf[i]);
   }
   if ( A_sndleng  != NULL ) ML_free(A_sndleng); 
   if ( A_rcvleng  != NULL ) ML_free(A_rcvleng);  
   if ( A_sndbuf   != NULL ) ML_free(A_sndbuf);
   if ( A_rcvbuf   != NULL ) ML_free(A_rcvbuf);  
   if ( A_recvlist != NULL ) ML_free(A_recvlist); 
   if ( A_sendlist != NULL ) ML_free(A_sendlist);
   if ( A_neigh    != NULL ) ML_free(A_neigh);

   /* ============================================================= */
   /* recover unamalgamated matrix and fetch communication info     */
   /* ============================================================= */

   /*
   Nrows     *= num_PDE_eqns;
   exp_Nrows = A_ntotal * num_PDE_eqns;
   if ( num_PDE_eqns > 1 ) 
   {
      ML_Operator_UnAmalgamateAndDropWeak(Amatrix, num_PDE_eqns, 0.0);
      if (exp_Nrows > 0) int_array = (int *) ML_allocate(exp_Nrows*sizeof(int));
      else               int_array = NULL;
      for ( i = 0; i < exp_Nrows; i+= num_PDE_eqns )
         for ( j = 0; j < num_PDE_eqns; j++ )
            int_array[i+j] = CF_array[i/num_PDE_eqns]; 
      if ( CF_array != NULL ) ML_free( CF_array );
      CF_array = int_array;
   }
   */
   exp_Nrows = A_ntotal;
   getrow_obj  = Amatrix->getrow;
   N_neighbors = getrow_obj->pre_comm->N_neighbors;
   nbytes      = N_neighbors * sizeof( int );
   if ( nbytes > 0 ) 
   {
      neighbors = (int *) ML_allocate( nbytes );
      recv_leng = (int *) ML_allocate( nbytes );
      send_leng = (int *) ML_allocate( nbytes );
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
   if ( nbytes > 0 ) send_list = (int *) ML_allocate(nbytes);
   else              send_list = NULL;
   if ( total_recv_leng+Nrows != exp_Nrows ) 
   {
      printf("%d : ML_AMG_CoarsenMIS - internal error.\n",mypid);
      printf("     lengths = %d %d \n",total_recv_leng+Nrows,exp_Nrows);
      exit(-1);
   }
   count = 0;
   for ( i = 0; i < N_neighbors; i++ ) 
   {
      for (j = 0; j < send_leng[i]; j++)
         send_list[count++] = getrow_obj->pre_comm->neighbors[i].send_list[j];
   }
   nbytes = total_recv_leng * sizeof( int );
   if ( nbytes > 0 ) recv_list = (int *) ML_allocate(nbytes);
   else              recv_list = NULL;
   count = 0;
   for ( i = 0; i < N_neighbors; i++ ) 
   {
      for (j = 0; j < recv_leng[i]; j++)
      {
         if ( getrow_obj->pre_comm->neighbors[i].rcv_list != NULL )
            recv_list[count] = getrow_obj->pre_comm->neighbors[i].rcv_list[j];
         else
            recv_list[count] = Nrows + count;
         count++;
      }
   }

   /* ============================================================= */
   /* Form prolongator                                              */
   /* ============================================================= */

   /* ============================================================= */
   /* update CF_array to find out how my neighbor's labelings       */
   /* ------------------------------------------------------------- */

   nbytes = total_send_leng * sizeof(int);
   if ( nbytes > 0 ) int_buf = (int *) ML_allocate( nbytes );
   else              int_buf = NULL;
   offset = 0;
   for ( i = 0; i < N_neighbors; i++ ) 
   {
      for ( j = 0; j < send_leng[i]; j++ ) 
         int_buf[offset+j] = CF_array[send_list[offset+j]];
      offset += send_leng[i];
   }
   msgtype = 35763;
   ML_Aggregate_ExchangeData((char*) &(CF_array[Nrows]), (char*) int_buf,
      N_neighbors, neighbors, recv_leng, send_leng, msgtype, ML_INT, comm);

   if ( int_buf != NULL ) ML_free(int_buf);
   Ncoarse = exp_Ncoarse = 0;;
   for (i = 0; i < Nrows; i++) if (CF_array[i] >= 0) CF_array[i] = Ncoarse++;
   exp_Ncoarse = Ncoarse;
   for ( i = Nrows; i < exp_Nrows; i++ )
      if ( CF_array[i] >= 0 ) CF_array[i] = exp_Ncoarse++;
   
   /* ============================================================= */
   /* allocate memory to hold the interpolation operator            */
   /* ------------------------------------------------------------- */

   nbytes = ( Nrows + 1 ) * sizeof(int);
   if ( nbytes > 0 ) ML_memory_alloc((void**)&new_ia, nbytes, "AM1");
   allocated = total_nnz = 0;
   rowi_col  = NULL;
   rowi_val  = NULL;
   for (i = 0; i < Nrows; i++) 
   {
      ML_get_matrix_row(Amatrix,1,&i,&allocated,&rowi_col,&rowi_val,&rowi_N,0);
      total_nnz += rowi_N;
   }
   nbytes = total_nnz * sizeof(int);
   if ( nbytes > 0 ) ML_memory_alloc((void**)&new_ja, nbytes, "AM2");
   else              new_ja = NULL;
   nbytes = total_nnz * sizeof(double);
   if ( nbytes > 0 ) ML_memory_alloc((void**)&new_val, nbytes, "AM3");
   else              new_val = NULL;
   
   /* ============================================================= */
   /* for each of the coarse grid point, create interpolant         */
   /* ------------------------------------------------------------- */

   sortleng  = A_ntotal / (8 * sizeof(int)) + 1;
   sort_array = (int *) ML_allocate( sortleng * sizeof(int) );
   for ( i = 0; i < sortleng; i++ ) sort_array[i] = 0;
   sizeint = sizeof(int) * 8;
   if ( sizeint == 16 )      logsizeint = 4;
   else if ( sizeint == 32 ) logsizeint = 5;
   else if ( sizeint == 64 ) logsizeint = 6;
   else                      logsizeint = 5;

   new_ia[0]   = 0;
   count       = 0;
   nnz_per_row = 0;
   for (i = 0; i < Nrows; i++) 
   {
      if ( CF_array[i] >= 0 ) /* ----- C point ----- */
      {
         new_val[count]  = 1.0;
         new_ja[count++] = CF_array[i];
         new_ia[i+1]     = count;
      }
      else /* ----- F points ----- */ 
      {
         /* ----- fetch the row i ----- */

         diag = 0.0;
         ML_get_matrix_row(Amatrix,1,&i,&allocated,&rowi_col,&rowi_val,
                           &rowi_N,0);

         /* ----- look for diagonal element ----- */

         if ( sys_unk_filter )
         {
            for (j = 0; j < rowi_N; j++) 
            {
               if (sys_info[rowi_col[j]] != sys_info[i]) 
               {
                  diag += rowi_val[j];
                  rowi_val[j] = 0.0; 
               } 
            } 
         } 
         else
         {
            for (j = 0; j < rowi_N; j++) 
               if ( rowi_col[j] == i ) diag = rowi_col[j]; 
         }
         if ( diag >= 0.0 ) idiag = 0; else idiag = -1;

         /* ----- compute row max and find number of C's connected ----- */

         rowmax = 0.0;
         numCi  = 0;
         for (j = 0; j < rowi_N; j++) 
         {
            col = rowi_col[j];
            if ( idiag >= 0 ) rowmax = ML_min(rowmax, rowi_val[j]);
            else              rowmax = ML_max(rowmax, rowi_val[j]);
            if ( CF_array[col] >= 0 ) numCi++;
         }
         rowmax *= epsilon;

         /* ----- get a list for C_i ----- */

         Ci_array = (int *)    ML_allocate( numCi * sizeof(int) );
         dsumCij  = (double *) ML_allocate( numCi * sizeof(double) );
         numCi    = 0;
         for (j = 0; j < rowi_N; j++) 
         {
            col = rowi_col[j];
            if ( CF_array[col] >= 0 )
            {
               if ( (idiag >= 0 && (rowi_val[j] < rowmax)) ||
                    (idiag <  0 && (rowi_val[j] > rowmax)) )  
               {
                  Ci_array[numCi++] = col;
                  intindex = col >> logsizeint;
                  bitindex = col % sizeint;
                  sort_array[intindex] |= ( 1 << bitindex );
               }
            }
         }
         ML_sort( numCi, Ci_array );
         for (j = 0; j < numCi; j++) 
         {
            col = Ci_array[j];
            new_ja[count+j] = CF_array[col];
            new_val[count+j] = 0.0;
         }

         /* ----- examine each of i's neighbors ----- */

         diag = 0.0;
         for (j = 0; j < rowi_N; j++) 
         {
            col = rowi_col[j];
            if ( col == i ) diag += rowi_val[j];
            else
            {
               /* ----- weak couplings, add to diagonal ----- */

               if ( (idiag >= 0 && ( rowi_val[j] >= rowmax ) ) || 
                    (idiag <  0 && ( rowi_val[j] <= rowmax ) ) ) 
                  diag += rowi_val[j];

               /* ----- strong coupling to 'C', put into array ----- */ 

               else if ( CF_array[col] >= 0 ) 
               {
                  k = ML_sorted_search(col, numCi, Ci_array);
                  new_val[count+k] += rowi_val[j];
               }

               /* ----- strong coupling to 'F' ----- */ 

               else if ( CF_array[col] < 0 ) 
               {
                  dsumCi = 0.0;
                  for (k = 0; k < numCi; k++) dsumCij[k] = 0.0; 

                  /* --- the F point is local to my processor ----- */

                  if ( col < Nrows )
                  {
                     for (k = rowptr[col]; k < rowptr[col+1] ; k++) 
                     {
                        ind2 = column[k];
                        if ( ind2 != col && CF_array[ind2] >= 0 ) 
                        {
                           intindex = ind2 >> logsizeint;
                           bitindex = ind2 % sizeint;
                           if ( sort_array[intindex] & ( 1 << bitindex ) ) 
                           {
                              m = ML_sorted_search(ind2, numCi, Ci_array);
                              dsumCij[m] += ( rowi_val[j] * values[k] );
                              dsumCi  += values[k];
                           }
                        }
                     }
                  }
                  else

                  /* --- the F point is in remote processor ----- */
                  {
                     if ( (col-Nrows) == 0 ) index = 0;
                     else                    index = offlengths[col-Nrows-1];
                     for (k = index; k < offlengths[col-Nrows] ; k++) 
                     {
                        ind2 = offibuffer[k];
                        if (ind2 >= 0 && CF_array[ind2] >= 0) 
                        {
                           intindex = ind2 >> logsizeint;
                           bitindex = ind2 % sizeint;
                           if ( sort_array[intindex] & ( 1 << bitindex ) ) 
                           {
                              m = ML_sorted_search(ind2, numCi, Ci_array);
                              dsumCij[m] += ( rowi_val[j] * offdbuffer[k] );
                              dsumCi  += offdbuffer[k];
                           }
                        }
                     }
                  }
                  for ( k = 0; k < numCi; k++ ) 
                  {
                     if ( dsumCij[k] != 0.0 )
                        new_val[count+k] = new_val[count+k] + dsumCij[k]/dsumCi; 
                  }
               }
            }
         }

         /* ----- put the element into the matrix ----- */

         for ( j = 0; j < numCi; j++ ) 
            new_val[count+j] = - ( new_val[count+j] / diag );
         count += numCi;
         ML_free( dsumCij );
         ML_free( Ci_array );

         /* ----- reset the sort_array ----- */

         for (j = 0; j < rowi_N; j++) 
         {
            col = rowi_col[j];
            if ( CF_array[col] >= 0 )  
            {
               intindex = col >> logsizeint;
               sort_array[intindex] = 0;
            }
         }

         /* ----- update the row pointer of P ----- */

         new_ia[i+1] = count;
         if ( (new_ia[i+1] - new_ia[i]) > nnz_per_row ) 
            nnz_per_row = new_ia[i+1] - new_ia[i];
      }
   }
   ML_free( rowi_col );
   ML_free( rowi_val );
/*
if ( mypid == 1 )
   for (i = 0; i < Nrows; i++) 
   {
      for (j = new_ia[i]; j < new_ia[i+1]; j++) 
         printf("P%d(%6d,%6d) = %e;\n", ml_amg->cur_level, i+1, 
                 new_ja[j]+1, new_val[j]);
   } 
MPI_Barrier(MPI_COMM_WORLD);
*/

   /* ------------------------------------------------------------- */
   /* compress the send and receive information                     */
   /* ------------------------------------------------------------- */

   offset = count = new_Nsend = 0;
   for ( i = 0; i < N_neighbors; i++ )
   {
      for ( j = 0; j < send_leng[i]; j++ )
      {
         index = send_list[offset+j];
         if ( CF_array[index] >= 0 ) send_list[count++] = CF_array[index];
      }
      offset += send_leng[i];
      if ( i == 0 ) send_leng[i] = count;
      else          send_leng[i] = count - tmpcnt;
      if ( send_leng[i] > 0 ) new_Nsend++;
      tmpcnt = count;
   }
   nbytes = new_Nsend * sizeof(int);
   if ( nbytes > 0 ) 
   {
      ML_memory_alloc((void**) &new_send_leng, nbytes, "SL1");
      ML_memory_alloc((void**) &new_send_neigh, nbytes, "SN1");
      ML_memory_alloc((void**) &new_send_list, count*sizeof(int), "SN1");
   } 
   else new_send_leng = new_send_list = new_send_neigh = NULL;
   new_Nsend = 0;
   for ( i = 0; i < N_neighbors; i++ )
   {
      if ( send_leng[i] > 0 )
      {
         new_send_leng[new_Nsend]    = send_leng[i];
         new_send_neigh[new_Nsend++] = neighbors[i];
      }
   }
   for ( i = 0; i < count; i++ ) new_send_list[i] = send_list[i];
   if ( send_list != NULL ) ML_free( send_list );
   if ( send_leng != NULL ) ML_free( send_leng );
   
   offset = count = new_Nrecv = 0;
   for ( i = 0; i < N_neighbors; i++ )
   {
      for ( j = 0; j < recv_leng[i]; j++ )
      {
         index = offset + j;
         if ( CF_array[Nrows+index] >= 0 ) count++;
      }
      offset += recv_leng[i];
      if ( i == 0 ) recv_leng[i] = count;
      else          recv_leng[i] = count - tmpcnt;
      if ( recv_leng[i] > 0 ) new_Nrecv++;
      tmpcnt = count;
   }
   nbytes = new_Nrecv * sizeof(int);
   if ( nbytes > 0 ) 
   {
      ML_memory_alloc((void**) &new_recv_leng, nbytes, "RL1");
      ML_memory_alloc((void**) &new_recv_neigh, nbytes, "RN1");
   } 
   else new_recv_leng = new_recv_neigh = NULL;
   new_Nrecv = 0;
   for ( i = 0; i < N_neighbors; i++ )
   {
      if ( recv_leng[i] > 0 )
      {
         new_recv_leng[new_Nrecv]    = recv_leng[i];
         new_recv_neigh[new_Nrecv++] = neighbors[i];
      }
   }
   if ( neighbors != NULL ) ML_free( neighbors );
   if ( recv_leng != NULL ) ML_free( recv_leng );
   if ( recv_list != NULL ) ML_free( recv_list );

   /* ------------------------------------------------------------- */
   /* set up the csr_data data structure                            */
   /* ------------------------------------------------------------- */

   ML_memory_alloc((void**)&csr_data,sizeof(struct ML_CSR_MSRdata),"CSR");
   csr_data->rowptr  = new_ia;
   csr_data->columns = new_ja;
   csr_data->values  = new_val;
/*
for ( i = 0; i < Nrows; i++ )
   for ( j = new_ia[i]; j < new_ia[i+1]; j++ )
      printf("P(%4d,%4d) = %e\n", i, new_ja[j], new_val[j]);
*/
   ML_Operator_Set_ApplyFuncData(*Pmatrix, Ncoarse, Nrows, ML_EMPTY, csr_data, 
                                 Nrows, NULL, 0);
   (*Pmatrix)->data_destroy = ML_CSR_MSR_ML_memorydata_Destroy;
   ML_memory_alloc((void**) &aggr_comm, sizeof(ML_Aggregate_Comm),"ACO");
   aggr_comm->comm = comm;
   aggr_comm->N_send_neighbors = new_Nsend;
   aggr_comm->N_recv_neighbors = new_Nrecv;
   aggr_comm->send_neighbors = new_send_neigh;
   aggr_comm->recv_neighbors = new_recv_neigh;
   aggr_comm->send_leng = new_send_leng;
   aggr_comm->recv_leng = new_recv_leng;
   aggr_comm->send_list = new_send_list;
   aggr_comm->local_nrows = Ncoarse;

   k = exp_Ncoarse - Ncoarse;
   ML_CommInfoOP_Generate( &((*Pmatrix)->getrow->pre_comm),
                           ML_Aggregate_ExchangeBdry, aggr_comm,
                           comm, Ncoarse, k);
   ML_Operator_Set_Getrow((*Pmatrix), ML_INTERNAL, Nrows, CSR_getrow);
   ML_Operator_Set_ApplyFunc((*Pmatrix), ML_INTERNAL, CSR_matvec);
   (*Pmatrix)->max_nz_per_row = nnz_per_row;

   /* ============================================================= */
   /* clean up                                                      */
   /* ------------------------------------------------------------- */

   if ( offibuffer != NULL ) ML_free( offibuffer );
   if ( offlengths != NULL ) ML_free( offlengths );
   if ( offdbuffer != NULL ) ML_free( offdbuffer );
   if ( rowptr     != NULL ) ML_free( rowptr );
   if ( column     != NULL ) ML_free( column );
   if ( values     != NULL ) ML_free( values );
   if ( CF_array   != NULL ) ML_free( CF_array );
   if ( sys_info   != NULL ) ML_free( sys_info );
   if ( new_Nsend > 0 )
   {
      ML_memory_free((void**) &new_send_leng);
      ML_memory_free((void**) &new_send_list);
      ML_memory_free((void**) &new_send_neigh);
   }
   if ( new_Nrecv > 0 )
   {
      ML_memory_free((void**) &new_recv_leng);
      ML_memory_free((void**) &new_recv_neigh);
   }
   ML_memory_free((void**) &aggr_comm);

   return Ncoarse;
}

/* ************************************************************************* */
/* A subroutine to label vertices of a particular type                       */
/* ------------------------------------------------------------------------- */

int ML_AMG_LabelVertices(int vlist_cnt2, int *vlist2, char Vtype,
                          char *vertex_state, char *vertex_type,
                          int nvertices, int *rptr, int *cptr, 
                          int myrank, int **proclist, int send_cnt, 
                          int **send_buf, int *send_proc, int *send_leng,
                          int recv_cnt, int **recv_buf, int *recv_proc, 
                          int *recv_leng, int **recv_list, 
                          ML_Comm *comm, int amg_index[])
{
   int     i, j, k, m, N_remaining_vertices, index, select_flag, fproc, col;
   int     NremainingRcvProcs, change_flag, *proc_flag, send_flag,nselected;
   int     *pref_list=NULL, col2, loop_cnt, pref_cnt;
   unsigned int nbytes;
   int     pref_index, *vlist, vlist_cnt;
   int     *pref_rank, vlist_ind;
   char    *in_preflist=NULL;
   USR_REQ *Request;
   int msg_type = 1041;

   /* ------------------------------------------------------------- */
   /* take out boundary vertices from the vertex list               */
   /* ------------------------------------------------------------- */

   nbytes = vlist_cnt2 * sizeof(int);
   if ( nbytes > 0 ) 
   {
      vlist = (int *) ML_allocate( nbytes );
      if ( vlist == NULL ) printf("MALLOC ERROR (LabelVertices) : vlist\n");
   }
   else              vlist = NULL;
   vlist_cnt = 0;
   for ( i = 0; i < vlist_cnt2; i++ )
   {
      index = vlist2[i];
      if ( vertex_state[index] == 'F' ) 
         vlist[vlist_cnt++] = index;
   }
   N_remaining_vertices = vlist_cnt;

   /* ------------------------------------------------------------- */
   /* set up communication buffers                                  */
   /* ------------------------------------------------------------- */

   NremainingRcvProcs = recv_cnt;
   send_flag          = 0;
   if ( recv_cnt > 0 )
   {
      nbytes    = recv_cnt * sizeof( USR_REQ );
      Request   = (USR_REQ *) ML_allocate( nbytes );
      nbytes    = recv_cnt * sizeof( int );
      proc_flag = (int *) ML_allocate( nbytes );
      if (proc_flag == NULL) printf("MALLOC ERROR (LabelVertices) : pflag.\n");
      for ( i = 0; i < recv_cnt; i++ ) proc_flag[i] = 0;
   }
   for ( j = 0; j < send_cnt; j++ )
      for ( k = 0; k <= send_leng[j]; k++ ) send_buf[j][k] = 0;

   /* ------------------------------------------------------------- */
   /* each processor updates the state of nodes belonging to other  */
   /* processors especially to take node of the boundary nodes      */
   /* ------------------------------------------------------------- */

   for ( i = 0; i < vlist_cnt2; i++ ) 
   {
      index = vlist2[i];
      if (vertex_state[index] == 'B' )
      {
         if ( proclist[index] != NULL )
         {
            for ( k = 0; k < proclist[index][0]; k++ ) 
            {
               fproc = proclist[index][2*k+1];
               m     = proclist[index][2*k+2];
               send_buf[fproc][m] = 3;
            }
         }
      }
   }
   ML_AMG_UpdateVertexStates(N_remaining_vertices, vertex_state,
	                     recv_cnt, recv_proc, recv_leng, recv_buf,
                     	     recv_list, proc_flag, &NremainingRcvProcs,
                     	     send_cnt, send_proc, send_leng, send_buf,
                     	     &send_flag, Request, comm, msg_type);

   /* ------------------------------------------------------------- */
   /* give the vertices with higher degree preferences (10/19/00)   */
   /* ------------------------------------------------------------- */

   if ( vlist_cnt > 0 )
   {
      nbytes      = vlist_cnt2 * sizeof(char);
      in_preflist = (char *) ML_allocate( nbytes );
      if (in_preflist == NULL) printf("MALLOC ERROR (LabelVertices) : inplist\n");
      for (i = 0; i < vlist_cnt2; i++) in_preflist[i] = 'f';

      nbytes    = vlist_cnt * sizeof(int);
      pref_rank = (int *) ML_allocate( nbytes );
      if ( pref_rank == NULL) printf("MALLOC ERROR (LabelVertices) : prank\n");
      for ( i = 0; i < vlist_cnt; i++ )
      {
         index = vlist[i];
         pref_rank[i] = rptr[index+1] - rptr[index];
      }
      ML_az_sort(pref_rank, vlist_cnt, vlist, NULL);
      ML_free( pref_rank );

      for ( i = 0; i < vlist_cnt/2; i++ )
      {
         j = vlist[i];
         vlist[i] = vlist[vlist_cnt-1-i];
         vlist[vlist_cnt-1-i] = j;
      }
      nbytes    = vlist_cnt * sizeof( int );
      pref_list = (int *) ML_allocate( nbytes );
      if (pref_list == NULL) printf("MALLOC ERROR (LabelVertices) : plist\n");
      pref_list[0] = vlist[0];
      in_preflist[vlist[0]] = 't';
      pref_cnt = 1;
   }   
   else 
   {
      pref_list = NULL;
      pref_cnt = 0;
   }

   /* ------------------------------------------------------------- */
   /* let's actually do coarsening                                  */
   /* ------------------------------------------------------------- */

   nselected   = 0;
   change_flag = 1;
   loop_cnt    = 0;     /* used to monitor the performance of coarsening */

   do {
      pref_index = 0;
      loop_cnt++;
      vlist_ind = 0;

      /* ---------------------------------------------------------- */
      /* reset all buffers to zero only if it has been changed      */
      /* ---------------------------------------------------------- */

      if ( change_flag == 1 )
      {
         for ( j = 0; j < send_cnt; j++ )
            for ( k = 0; k <= send_leng[j]; k++ ) send_buf[j][k] = 0;
         change_flag = 0;
      }

      /* ---------------------------------------------------------- */
      /* examine the vertices in vlist                              */
      /* ---------------------------------------------------------- */

      for ( i = 0; i < vlist_cnt; i++ )
      {
         /* ------------------------------------------------------- */
         /* handle the preference list first, if there is any       */
         /* Note : we want to fetch the pref_list from the front    */
         /* ------------------------------------------------------- */

         index = -1;
         if ( pref_cnt > pref_index ) 
         {
            /* ---------------------------------------------------- */
            /* search for the first F(ree) vertex                   */
            /* ---------------------------------------------------- */

            for (j = pref_index; j < pref_cnt; j++) 
            {
               index = pref_list[j];    
               if ( index < 0 || index >= nvertices )
               {
                  printf("LabelVertices ERROR : in pref_list %d %d\n",j,index);
                  exit(1);
               }
               if ( vertex_state[index] == 'F' ) break;
            }

            /* ---------------------------------------------------- */
            /* if found, take it off from the pref_list             */
            /* ---------------------------------------------------- */

            if ( j != pref_cnt ) pref_index = j + 1;
            else 
            {
               index = -1;
               pref_index = pref_cnt;
            }
         }

         /* ------------------------------------------------------- */
         /* if not found from the pref_list, look for one in the    */
         /* incoming vertex list (sorted)                           */
         /* ------------------------------------------------------- */

         if ( index == -1 )
         {
            for ( j = vlist_ind; j < vlist_cnt; j++ )
            {
               index = vlist[j];
               if ( vertex_state[index] == 'F' ) break;
            }
            if ( j == vlist_cnt ) break;
            else                  vlist_ind = j + 1;
         }

         /* ------------------------------------------------------- */
         /* at this point it is assumed that vertex_state[index]    */
         /* is F(ree)                                               */
         /* ------------------------------------------------------- */

         select_flag = 1;
         for ( j = rptr[index]; j < rptr[index+1]; j++ )
         {
            /* if any of its neighbor is selected, delete this node */

            col = cptr[j];
            if ( vertex_state[col] == 'S' )
            {
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

               /* ------------------------------------------------- */
               /* put the next set of vertices into the preference  */
               /* list (try to mimic the sequential MIS algorithm)  */
               /* ------------------------------------------------- */

               for ( k = rptr[index]; k < rptr[index+1]; k++ ) 
               {
                  col2 = cptr[k];
                  if (col2 < nvertices && vertex_state[col2] == 'F' &&
                      vertex_type[col2] == Vtype ) 
                  {
                     if (in_preflist[col2] != 't') 
                     {
                        if ( pref_cnt >= vlist_cnt ) 
                        {
                           printf("LabelVertices ERROR : pref_cnt too long\n");
                           exit(1);
                        }
                        pref_list[pref_cnt++] = col2;
                        in_preflist[col2] = 't';
                     }
                  }
               }
               select_flag = 0;
               break;
            }
               
            /* ---------------------------------------------------- */
            /* If its neighbor is of the same type and not been     */
            /* considered. Furthermore, if it is a remote vertex    */
            /* and its owner processor has rank smaller than mine,  */
            /* my processor should wait(thus turn off select_flag)  */
            /* ---------------------------------------------------- */

            else if (vertex_type[col] == Vtype && vertex_state[col] == 'F')
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

         /* ------------------------------------------------------- */
         /* if the vertex in question is not any of those           */
         /* considered above, select this vertex.                   */
         /* ------------------------------------------------------- */

         if ( select_flag == 1 )
         {
            if ((vertex_state[index] == 'F') && (index < nvertices)) 
               N_remaining_vertices--;
            vertex_state[index] = 'S';
            amg_index[index] = nselected;
            nselected++;

            /* ---------------------------------------------------- */
            /* set the flag that this vertex has been selected in   */
            /* the buffer which is to be sent to other processors   */
            /* ---------------------------------------------------- */

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

            /* ---------------------------------------------------- */
            /* put the next set of vertices into the preference     */
            /* list (try to mimic the sequential MIS algorithm)     */
            /* ---------------------------------------------------- */

            for ( j = rptr[index]; j < rptr[index+1]; j++ ) 
            {
               col = cptr[j];
               if (col < nvertices && vertex_state[col] == 'F' &&
                   vertex_type[col] == Vtype ) 
               {
                  if (in_preflist[col] != 't') 
                  {
                     if ( pref_cnt >= vlist_cnt ) 
                     {
                        printf("LabelVertices ERROR : pref_cnt too long\n");
                        exit(1);
                     }
                     pref_list[pref_cnt++] = col;
                     in_preflist[col] = 't';
                  }
               }
/*
               if (col2 < nvertices && vertex_state[col2] == 'F' &&
                   vertex_type[col2] == Vtype ) 
               {
                  for ( k = rptr[col]; k < rptr[col+1]; k++ ) 
                  {
                     col2 = cptr[k];
                     if (col2 < nvertices && vertex_state[col2] == 'F' &&
                         vertex_type[col2] == Vtype ) 
                     {
                        if (in_preflist[col2] != 't') 
                        {
                           if ( pref_cnt >= vlist_cnt ) 
                           {
                              printf("LabelVertices ERROR : pref_cnt too long\n");
                              exit(1);
                           }
                           pref_list[pref_cnt++] = col2;
                           in_preflist[col2] = 't';
                        }
                     }
                  }
               }
*/
            }
         }
      }

      /* ---------------------------------------------------------- */
      /* communicate the state information to other processors      */ 
      /* ---------------------------------------------------------- */

      ML_AMG_UpdateVertexStates(N_remaining_vertices, vertex_state,
	                  recv_cnt, recv_proc, recv_leng, recv_buf,
	                  recv_list, proc_flag, &NremainingRcvProcs,
	                  send_cnt, send_proc, send_leng, send_buf,
	                  &send_flag, Request, comm, msg_type);

   } while ( NremainingRcvProcs > 0 || N_remaining_vertices > 0 );

   if ( recv_cnt > 0 )
   {
      ML_free( proc_flag );
      ML_free( Request );
   }
   if ( vlist_cnt  > 0 ) ML_free( pref_list );
   if ( vlist_cnt2 > 0 ) ML_free( vlist );
   if ( vlist_cnt  > 0 ) ML_free( in_preflist );
   /* return nselected; */
   return loop_cnt;
}

/* ************************************************************************* */
/* A subroutine to update vertex states                                      */
/* ------------------------------------------------------------------------- */

int ML_AMG_UpdateVertexStates(int N_remaining_vertices, char vertex_state[], 
                 int recv_cnt, int recv_proc[], int recv_leng[], 
                 int **recv_buf, int **recv_list, int proc_flag[], 
                 int *NremainingRcvProcs, int send_cnt, int send_proc[], 
                 int send_leng[], int **send_buf, int *send_flag, 
                 USR_REQ *Request, ML_Comm *comm, int msgtype) 
{
   int j, k, kkk, fproc;
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
   if ( *send_flag == 0 ) 
   {
      for ( j = 0; j < send_cnt; j++ ) 
      {
         nbytes = (send_leng[j] + 1) * sizeof( int );
         if ( N_remaining_vertices <= 0 ) 
         { 
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
         for ( k = 0; k < recv_leng[j]; k++ )
         {
            kkk = recv_list[j][k];
            if      (recv_buf[j][k] == 1) vertex_state[kkk] = 'S';
            else if (recv_buf[j][k] == 2) vertex_state[kkk] = 'D';
            else if (recv_buf[j][k] == 3) vertex_state[kkk] = 'B';
         }
         if ( recv_buf[j][recv_leng[j]] == 1 )
         {
            proc_flag[j] = 1;
            (*NremainingRcvProcs)--;
         }
      }
   }
   return 0;
}

/* ************************************************************************* */
/* A subroutine to get communication information                             */
/* ------------------------------------------------------------------------- */

int ML_AMG_GetCommInfo(ML_CommInfoOP *mat_comm, int Nrows, int *A_Nneigh, 
                       int **A_neigh, int ***A_sendlist, int ***A_recvlist, 
                       int ***A_sndbuf, int ***A_rcvbuf, int **A_sndleng, 
                       int **A_rcvleng, int *Nghost)
{
   int i, j, max_element;

   if ( mat_comm == NULL )
   {
      (*A_Nneigh)   = 0;
      (*A_neigh)    = NULL;
      (*A_sendlist) = NULL;
      (*A_recvlist) = NULL;
      (*A_sndbuf)   = NULL;
      (*A_rcvbuf)   = NULL;
      (*A_sndleng)  = NULL;
      (*A_rcvleng)  = NULL;
      (*Nghost)     = 0;
      return 0;
   }
      
   (*A_Nneigh) = ML_CommInfoOP_Get_Nneighbors(mat_comm);
   if ( (*A_Nneigh) > 0 )
   {
      (*A_neigh)    = ML_CommInfoOP_Get_neighbors(mat_comm);
      (*A_sendlist) = (int **) ML_allocate(sizeof(int *)*(*A_Nneigh));
      (*A_recvlist) = (int **) ML_allocate(sizeof(int *)*(*A_Nneigh));
      (*A_sndbuf)   = (int **) ML_allocate(sizeof(int *)*(*A_Nneigh));
      (*A_rcvbuf)   = (int **) ML_allocate(sizeof(int *)*(*A_Nneigh));
      (*A_sndleng)  = (int  *) ML_allocate(sizeof(int  )*(*A_Nneigh));
      (*A_rcvleng)  = (int  *) ML_allocate(sizeof(int  )*(*A_Nneigh));
   }
   else
   {
      (*A_sendlist) = (*A_recvlist) = (*A_rcvbuf) = (*A_sndbuf) = NULL;
      (*A_neigh) = (*A_rcvleng) = (*A_sndleng) = NULL;
   }

   max_element = Nrows - 1;
   for (i = 0; i < (*A_Nneigh); i++) 
   {
      (*A_recvlist)[i] = ML_CommInfoOP_Get_rcvlist(mat_comm,(*A_neigh)[i]);
      (*A_rcvleng)[i]  = ML_CommInfoOP_Get_Nrcvlist (mat_comm,(*A_neigh)[i]);
      (*A_sendlist)[i] = ML_CommInfoOP_Get_sendlist (mat_comm,(*A_neigh)[i]);
      (*A_sndleng)[i]  = ML_CommInfoOP_Get_Nsendlist(mat_comm,(*A_neigh)[i]);
      (*A_rcvbuf)[i]   = (int *) ML_allocate(sizeof(int)*((*A_rcvleng)[i]+1));
      (*A_sndbuf)[i]   = (int *) ML_allocate(sizeof(int)*((*A_sndleng)[i]+1));
      for (j = 0; j < (*A_rcvleng)[i]; j++) 
         if ((*A_recvlist)[i][j] > max_element ) 
            max_element = (*A_recvlist)[i][j];
   }
   (*Nghost) = max_element - Nrows + 1;
   return 0;
}

/* ************************************************************************* */
/* newly added to generate coarse grid based on compatible relaxation        */
/* using symmetric Gauss Seidel                                              */
/* ------------------------------------------------------------------------- */

int ML_AMG_CompatibleRelaxation(int *CF_array,
                                ML_Operator *Amat, int *Ncoarse, int limit)
{
   int           i, j, iter, Nrows, length, allocated, *cols, col;
   int           *indices, count;
   double        *initsol, *rhs, *sol;
   double        dtemp, diag_term, *vals, omega=1.0;
   ML_Comm       *comm;
   ML_CommInfoOP *getrow_comm;

   /* ------------------------------------------------------------- */
   /* fetch Amatrix information                                     */
   /* ------------------------------------------------------------- */

   comm        = Amat->comm;
   Nrows       = Amat->getrow->Nrows;
   getrow_comm = Amat->getrow->pre_comm;

   /* ------------------------------------------------------------- */
   /* set up arrays (rhs and initial guesses)                       */
   /* ------------------------------------------------------------- */

   indices = (int *)    ML_allocate( Nrows * sizeof(int) );
   initsol = (double *) ML_allocate( Nrows * sizeof(double) );
   sol     = (double *) ML_allocate((Nrows+getrow_comm->total_rcv_length+1)*
                               sizeof(double));
   rhs     = (double *) ML_allocate( Nrows * sizeof(double) );
   for ( i = 0; i < Nrows; i++ ) rhs[i] = 0.0;
   ML_random_vec(initsol, Nrows, comm);
   for ( i = 0; i < Nrows; i++ ) initsol[i] += 1.0;
   for ( i = 0; i < Nrows; i++ ) sol[i] = initsol[i];
   for ( i = 0; i < Nrows; i++ ) 
      if ( CF_array[i] >= 0 ) sol[i] = 0.0;

   /* ------------------------------------------------------------- */
   /* allocate matrix buffers                                       */
   /* ------------------------------------------------------------- */

   allocated = Amat->max_nz_per_row+1;
   cols = (int    *) ML_allocate(allocated * sizeof(int   ));
   vals = (double *) ML_allocate(allocated * sizeof(double));

   /* ------------------------------------------------------------- */
   /* iterate                                                       */
   /* ------------------------------------------------------------- */

   for (iter = 0; iter < 2; iter++)
   {
      if (getrow_comm != NULL)
         ML_exchange_bdry(sol, getrow_comm, Nrows, comm, ML_OVERWRITE,NULL);

      for (i = 0; i < Nrows; i++)
      {
         if ( CF_array[i] < 0 )
         {
            dtemp = 0.0;
            diag_term = 0.0;
            ML_get_matrix_row(Amat,1,&i,&allocated,&cols,&vals,&length,0);
            for (j = 0; j < length; j++)
            {
               col = cols[j];
               dtemp += vals[j]*sol[col];
               if (col == i) diag_term = vals[j];
            }
            if (diag_term != 0.0)
               sol[i] += omega*(rhs[i] - dtemp)/diag_term;
         }
      }

      if (getrow_comm != NULL)
         ML_exchange_bdry(sol, getrow_comm, Nrows, comm, ML_OVERWRITE,NULL);

      for (i = Nrows-1; i >= 0; i--)
      {
         if ( CF_array[i] < 0 )
         {
            dtemp = 0.0;
            diag_term = 0.0;
            ML_get_matrix_row(Amat,1,&i,&allocated,&cols,&vals,&length,0);
            for (j = 0; j < length; j++)
            {
               col = cols[j];
               dtemp += vals[j]*sol[col];
               if (col == i) diag_term = vals[j];
            }
            if (diag_term != 0.0)
               sol[i] += omega*(rhs[i] - dtemp)/diag_term;
         }
      }
   }
   ML_free(vals); ML_free(cols);

   /* ------------------------------------------------------------- */
   /* sort the solution                                             */
   /* ------------------------------------------------------------- */

   for ( i = 0; i < Nrows; i++ ) indices[i] = i;
   for ( i = 0; i < Nrows; i++ ) sol[i] = ML_dabs(sol[i]) / initsol[i];
   ML_split_dsort( sol, Nrows, indices, limit );
   count = 0;
   for ( i = 0; i < Nrows; i++ )
   {
      if ( CF_array[indices[i]] < 0 ) 
      {
         CF_array[indices[i]] = (*Ncoarse)++;
         count++;
         if ( count >= limit ) break;
      }
   }
   ML_free(indices);
   ML_free(initsol);
   ML_free(sol);
   ML_free(rhs);
   return 0;
}

