/* ************************************************************************* */
/* ************************************************************************* */
/* Functions to create tentative prolongators (domain decomposition)         */
/* ************************************************************************* */
/* Author        : Charles Tong (LLNL)                                       */
/* Date          : August, 2000                                              */
/* ************************************************************************* */
/* Local Function :                                                          */
/*    ML_Aggregate_CoarsenDomainDecomp                                       */
/* ************************************************************************* */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ml_aggregate.h"
#include "ml_lapack.h"
#include "ml_utils.h"


/* ************************************************************************* */
/* ************************************************************************* */
/* ML_Aggregate_CoarsenDomainDecomp subroutine.                              */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_CoarsenDomainDecomp( ML_Aggregate *ml_ag,
       ML_Operator *Amatrix, ML_Operator **Pmatrix, ML_Comm *comm)
{
   int     i, j, k, m, jnode, index, index3, count; 
   int     nz_cnt, nbytes, length, level, diff_level;
   int     Nrows, *mat_indx=NULL;
   int     maxnnz_per_row=500, *col_ind;
   int     mypid;
   int     aggr_count, *aggr_index;
   int     *aggr_cnt_array, max_agg_size, **rows_in_aggs;
   int     Ncoarse, *new_ia, *new_ja, new_Nrows;
   int     num_PDE_eqns, nullspace_dim, lwork, info;
   double  *col_val, *diagonal=NULL, dcompare1, dcompare2, *new_val=NULL;
   double  epsilon, *nullspace_vect=NULL, *qr_tmp=NULL;
   double  *tmp_vect=NULL, *work=NULL, *new_null=NULL;
   int     (*getrowfunc)(ML_Operator *,int,int*,int,int*,double*,int*);
   void    *getrowdata;
   struct ML_CSR_MSRdata *csr_data;
   ML_Aggregate_Comm     *aggr_comm;
   ML_GetrowFunc         *getrow_obj;

   /* ============================================================= */
   /* get machine and matrix information                            */
   /* ============================================================= */

   mypid          = comm->ML_mypid;
   nullspace_dim  = ml_ag->nullspace_dim;
   nullspace_vect = ml_ag->nullspace_vect;
   Nrows          = Amatrix->outvec_leng;
   num_PDE_eqns   = ml_ag->num_PDE_eqns;

   /* ============================================================= */
   /* initialize and update the threshold                           */
   /* ============================================================= */

   diff_level = ml_ag->begin_level - ml_ag->cur_level;
   if ( diff_level == 0 ) ml_ag->curr_threshold = ml_ag->threshold;
   epsilon = ml_ag->curr_threshold;
   ml_ag->curr_threshold *= 0.5;

#ifdef ML_AGGR_DEBUG
   if ( mypid == 0 )
   {
      printf("ML_Aggregate_DomainDecomp : current level = %d\n",
                           ml_ag->cur_level);
      printf("ML_Aggregate_DomainDecomp : current eps = %e\n", epsilon);
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
      printf("ML_Aggregate_DomainDecomp ERROR : null getrowfunc.\n");
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

   count = 0;
   for ( i = 0; i < Nrows; i++ )
   {
      diagonal[i]     = 0.0;
      while (getrowfunc((ML_Operator *) getrowdata,1,&i,maxnnz_per_row,col_ind,
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
         if ( col_ind[j] == i ) diagonal[i] = col_val[j];
      }
      count += m;
      if ( diagonal[i] == 0.0 )
      {
         printf("%d : DomainDecomp WARNING - diag %d is 0.\n",
                mypid,i);
         count++;
      }
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
#ifdef ML_AGGR_DEBUG
   if ( mypid == 0 )
      printf("Aggregation(DD) : Total nnz = %d (Nrows=%d)\n",m,k);
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
      getrowfunc((ML_Operator *) getrowdata,1,&i,maxnnz_per_row,col_ind,col_val,&m);
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
      if ( nz_cnt == mat_indx[i] && length > 0 ) mat_indx[nz_cnt++] = i;
      mat_indx[i+1] = nz_cnt;
      ML_sort(mat_indx[i+1]-mat_indx[i], mat_indx+mat_indx[i]);
   }
   if ( col_ind  != NULL ) ML_memory_free((void**) &col_ind);
   if ( col_val  != NULL ) ML_memory_free((void**) &col_val);
   if ( diagonal != NULL ) ML_memory_free((void**) &diagonal);

   /* ============================================================= */
   /* create a domain decomposed aggregate                          */
   /* ============================================================= */

   aggr_count = 1;
   ML_memory_alloc((void**) &aggr_cnt_array, sizeof(int), "ACJ");
   aggr_cnt_array[0] = Nrows;
   nbytes = Nrows * sizeof(int);
   if ( nbytes > 0 ) ML_memory_alloc((void**) &aggr_index, nbytes, "ACJ");
   else              aggr_index = NULL;
   for ( i = 0; i < Nrows; i++ ) aggr_index[i] = 0;

   /* ============================================================= */
   /* Form tentative prolongator                                    */
   /* ============================================================= */

   Ncoarse = aggr_count;

   /* ============================================================= */
   /* check and copy aggr_index to ml_ag (for block smoothing)      */
   /* ------------------------------------------------------------- */

   level = ml_ag->cur_level;
   nbytes = Nrows * sizeof( int );
   if ( nbytes > 0 )
      ML_memory_alloc((void**) &(ml_ag->aggr_info[level]), nbytes, "ACK");
   else
      ml_ag->aggr_info[level] = NULL;

   for ( i = 0; i < Nrows; i++ )
   {
      ml_ag->aggr_info[level][i] = aggr_index[i];
   }
   ml_ag->aggr_count[level] = aggr_count; 

   /* ============================================================= */
   /* set up the new operator                                       */
   /* ------------------------------------------------------------- */

   new_Nrows = Nrows;
   for ( i = 0; i < new_Nrows; i++ )
   {
      if ( aggr_index[i] >= Ncoarse )
         printf("WARNING : index out of bound %d = %d(%d)\n",i,aggr_index[i],
                Ncoarse);
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
         printf("       requested = %d\n",aggr_cnt_array[i]);
         exit(1);
      }
   }
   for (i = 0; i < Nrows; i++)
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
            for (k = 0; k < nullspace_dim; k++) qr_tmp[k*length+j] = 1.0;
      }
      else
      {
         for (k = 0; k < nullspace_dim; k++)
         {
            for (j = 0; j < length; j++)
            {
               index = rows_in_aggs[i][j];
               qr_tmp[k*length+j] = nullspace_vect[k*Nrows+index];
            }
         }
      }

      /* ---------------------------------------------------------- */
      /* now calculate QR using an LAPACK routine                   */
      /* ---------------------------------------------------------- */

      DGEQRF_F77(&(aggr_cnt_array[i]), &nullspace_dim, qr_tmp,
               &(aggr_cnt_array[i]), tmp_vect, work, &lwork, &info);
      if (info != 0)
         pr_error("CoarsenCoupled ERROR : dgeqrf returned a non-zero\n");

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
      DORGQR_F77(&(aggr_cnt_array[i]),&nullspace_dim,&nullspace_dim,
              qr_tmp, &(aggr_cnt_array[i]), tmp_vect, work, &lwork, &info);
      if (info != 0)
         pr_error("CoarsenCoupled ERRO : dorgqr returned a non-zero\n");

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

      for (j = 0; j < aggr_cnt_array[i]; j++)
      {
         index = rows_in_aggs[i][j];
         index3 = new_ia[index];
         for (k = 0; k < nullspace_dim; k++)
         {
            new_ja [index3+k] = i * nullspace_dim + k;
            new_val[index3+k] = qr_tmp[ k*aggr_cnt_array[i]+j];
         }
      }
   }

   ML_Aggregate_Set_NullSpace(ml_ag, num_PDE_eqns, nullspace_dim,
                              new_null, Ncoarse*nullspace_dim);
   
   ML_memory_free( (void **) &new_null);

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
         printf("%d : CoarsenDD WARNING : rowsum(P(%d)) = 0 (%d)\n",
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
   aggr_comm->N_send_neighbors = 0;
   aggr_comm->N_recv_neighbors = 0;
   aggr_comm->send_neighbors = NULL;
   aggr_comm->recv_neighbors = NULL;
   aggr_comm->send_leng = NULL;
   aggr_comm->recv_leng = NULL;
   aggr_comm->send_list = NULL;
   aggr_comm->local_nrows = Ncoarse * nullspace_dim;

   ML_CommInfoOP_Generate( &((*Pmatrix)->getrow->pre_comm),
                           ML_Aggregate_ExchangeBdry, aggr_comm,
                           comm, Ncoarse*nullspace_dim, 0);
   ML_Operator_Set_Getrow((*Pmatrix), Nrows, CSR_getrow);
   ML_Operator_Set_ApplyFunc((*Pmatrix), CSR_matvec);
   (*Pmatrix)->max_nz_per_row = nullspace_dim;

   /* ------------------------------------------------------------- */
   /* clean up                                                      */
   /* ------------------------------------------------------------- */

   ML_memory_free((void**) &mat_indx);
   ML_memory_free((void**) &aggr_index);
   ML_memory_free((void**) &aggr_cnt_array);
   for (i = 0; i < aggr_count; i++) ML_free(rows_in_aggs[i]);
   ML_memory_free((void**)&rows_in_aggs);
   ML_memory_free((void**)&qr_tmp);
   ML_memory_free((void**)&tmp_vect);
   ML_memory_free((void**)&work);
   ML_memory_free((void**) &aggr_comm);
   return Ncoarse*nullspace_dim;
}

