/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
/* ************************************************************************* */
/* ************************************************************************* */
/* Functions to create tentative prolongators using METIS to create the      */
/* aggregate.                                                                */
/* ************************************************************************* */
/* Author        : Marzio Sala (SNL)                                         */
/* Date          : December 2004                                             */
/* ************************************************************************* */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ml_aggregate.h"
#include "ml_lapack.h"
#include "ml_utils.h"
#include "ml_agg_METIS.h" 
#include "ml_agg_user.h"
#include "ml_viz_stats.h"
#include "ml_agg_info.h"

/* ---------------------------------- */
/* Set label for user's defined label */
/* ---------------------------------- */

static char*(*UserLabel_)() = NULL;

char* ML_GetUserLabel()
{
  assert (UserLabel_ != NULL);
  return((*UserLabel_)());
}

int ML_SetUserLabel(char*(user)())
{
  UserLabel_ = user;
  return(0);
}

/* --------------------------------- */
/* function to define the aggregates */
/* --------------------------------- */

static int (*UserPartitions_)(ML_Operator* Amat, char* bdry_nodes,
                              double epsilon,
                              double* x,double* y,double* z,
                              int* partitions, int *LocalNonzeros) = NULL;

int ML_GetUserPartitions(ML_Operator* Amat, char* bdry_nodes,
                         double epsilon,
                         double* x,double* y,double* z,
                         int* partitions, int* LocalNonzeros)
{
  assert (UserPartitions_ != NULL);
  return((*UserPartitions_)(Amat,bdry_nodes,epsilon,
                            x,y,z,partitions,LocalNonzeros));
}

int ML_SetUserPartitions(int (user)(ML_Operator* Amat, char* bdry_nodes,
                                    double epsilon,
                                    double* x,double* y,double* z,
                                    int* partitions, int* LocalNonzeros))
{
  UserPartitions_ = user;
  return(0);
}

/* ------------------------------------------------------------------------ */

int ML_Aggregate_CoarsenUser(ML_Aggregate *ml_ag, ML_Operator *Amatrix, 
			      ML_Operator **Pmatrix, ML_Comm *comm)
{
  unsigned int nbytes, length;
  int     i, j,  k, Nrows, exp_Nrows;
  int     diff_level;
  int     aggr_count, index, mypid, num_PDE_eqns;
  int     *aggr_index = NULL, nullspace_dim;
  int     Ncoarse, count;
  int     *new_ia = NULL, *new_ja = NULL, new_Nrows;
  int     exp_Ncoarse;
  int     *aggr_cnt_array = NULL;
  int     level, index3, max_agg_size;
  int     **rows_in_aggs = NULL, lwork, info;
  double  *new_val = NULL, epsilon;
  double  *nullspace_vect = NULL, *qr_tmp = NULL;
  double  *tmp_vect = NULL, *work = NULL, *new_null = NULL;
  ML_SuperNode          *aggr_head = NULL, *aggr_curr, *supernode;
  struct ML_CSR_MSRdata *csr_data;
  int                   total_nz = 0;
  char str[80];

  int * graph_decomposition = NULL;
  ML_Aggregate_Viz_Stats * aggr_viz_and_stats;
  ML_Aggregate_Viz_Stats * grid_info;
  int Nprocs;
  char * unamalg_bdry = NULL;
  char* label;
  int N_dimensions;
  double* x_coord = NULL;
  double* y_coord = NULL;
  double* z_coord = NULL;

  /* ------------------- execution begins --------------------------------- */

  label =  ML_GetUserLabel();
  sprintf(str, "%s (level %d) :", label, ml_ag->cur_level);

  /* ============================================================= */
  /* get the machine information and matrix references             */
  /* ============================================================= */

  mypid                   = comm->ML_mypid;
  Nprocs                  = comm->ML_nprocs;
  epsilon                 = ml_ag->threshold;
  num_PDE_eqns            = ml_ag->num_PDE_eqns;
  nullspace_dim           = ml_ag->nullspace_dim;
  nullspace_vect          = ml_ag->nullspace_vect;
  Nrows                   = Amatrix->outvec_leng;

  if (mypid == 0 && 5 < ML_Get_PrintLevel()) {
    printf("%s num PDE eqns = %d\n",
           str,
           num_PDE_eqns);
  }

  /* ============================================================= */
  /* check the system size versus null dimension size              */
  /* ============================================================= */

  if ( Nrows % num_PDE_eqns != 0 )
  {
    printf("ML_Aggregate_CoarsenUser ERROR : Nrows must be multiples");
    printf(" of num_PDE_eqns.\n");
    exit(EXIT_FAILURE);
  }
  diff_level = ml_ag->max_levels - ml_ag->cur_level - 1;
  if ( diff_level > 0 ) num_PDE_eqns = nullspace_dim; /* ## 12/20/99 */

  /* ============================================================= */
  /* set up the threshold for weight-based coarsening              */
  /* ============================================================= */

  diff_level = ml_ag->begin_level - ml_ag->cur_level;
  if (diff_level == 0) 
    ml_ag->curr_threshold = ml_ag->threshold;
  epsilon = ml_ag->curr_threshold;
  ml_ag->curr_threshold *= 0.5;

  if (mypid == 0 && 7 < ML_Get_PrintLevel())
    printf("%s current eps = %e\n", str, epsilon);

  epsilon = epsilon * epsilon;

  ML_Operator_AmalgamateAndDropWeak(Amatrix, num_PDE_eqns, epsilon);
  Nrows /= num_PDE_eqns;

  exp_Nrows = Nrows;

  /* ********************************************************************** */
  /* allocate memory for aggr_index, which will contain the decomposition   */
  /* ********************************************************************** */

  nbytes = (Nrows*num_PDE_eqns) * sizeof(int);

  if ( nbytes > 0 ) {
    ML_memory_alloc((void**) &aggr_index, nbytes, "ACJ");
    if( aggr_index == NULL ) {
      fprintf( stderr,
              "*ML*ERR* not enough memory for %d bytes\n"
              "*ML*ERR* (file %s, line %d)\n",
              nbytes,
              __FILE__,
              __LINE__ );
      exit( EXIT_FAILURE );
    }
  }
  else              aggr_index = NULL;

  for( i=0 ; i<Nrows*num_PDE_eqns ; i++ ) aggr_index[i] = -1;

  unamalg_bdry = (char *) ML_allocate( sizeof(char) * (Nrows+1) );

  if( unamalg_bdry == NULL ) {
    fprintf( stderr,
            "*ML*ERR* on proc %d, not enough space for %d bytes\n"
            "*ML*ERR* (file %s, line %d)\n",
            mypid,
            (int)sizeof(char) * Nrows,
            __FILE__,
            __LINE__ );
    exit( EXIT_FAILURE );
  }

  N_dimensions = ml_ag->N_dimensions;
  grid_info = (ML_Aggregate_Viz_Stats*) Amatrix->to->Grid->Grid;
  x_coord = grid_info->x;

  if (N_dimensions > 1 && x_coord)
    y_coord = grid_info->y;
  else
    y_coord = 0;
  if (N_dimensions > 2 && x_coord)
    z_coord = grid_info->z;
  else
    z_coord = 0;

  aggr_count = ML_GetUserPartitions(Amatrix,unamalg_bdry,
                                    epsilon,
                                    x_coord,y_coord,z_coord,
                                    aggr_index,&total_nz);

#ifdef ML_MPI
  MPI_Allreduce( &Nrows, &i, 1, MPI_INT, MPI_SUM, Amatrix->comm->USR_comm );
  MPI_Allreduce( &aggr_count, &j, 1, MPI_INT, MPI_SUM, Amatrix->comm->USR_comm );
#else
  i = Nrows;
  j = aggr_count;
#endif

  if( mypid == 0 && 7 < ML_Get_PrintLevel() ) {
    printf("%s Using %d (block) aggregates (globally)\n",
           str,
           j );
    printf("%s # (block) aggre/ # (block) rows = %8.5f %% ( = %d / %d)\n",
           str,
           100.0*j/i,
           j, i);
  }

  j = ML_gsum_int( aggr_count, comm );
  if (mypid == 0 && 7 < ML_Get_PrintLevel())  {
    printf("%s %d (block) aggregates (globally)\n",
           str, j );
  }   

  /* ********************************************************************** */
  /* I allocate room to copy aggr_index and pass this value to the user,    */
  /* who will be able to analyze and visualize this after the construction  */
  /* of the levels. This way, the only price we have to pay for stats and   */
  /* viz is essentially a little bit of memory.                             */
  /* this memory will be cleaned with the object ML_Aggregate ml_ag.        */
  /* I set the pointers using the ML_Aggregate_Info structure. This is      */
  /* allocated using ML_Aggregate_Info_Setup(ml,MaxNumLevels)               */
  /* ********************************************************************** */

  if (Amatrix->to->Grid->Grid != NULL) {

    graph_decomposition = (int *)ML_allocate(sizeof(int)*(Nrows+1));
    if( graph_decomposition == NULL ) {
      fprintf( stderr,
              "*ML*ERR* Not enough memory for %d bytes\n"
              "*ML*ERR* (file %s, line %d)\n",
              (int)sizeof(int)*Nrows,
              __FILE__,
              __LINE__ );
      exit( EXIT_FAILURE );
    }

    for( i=0 ; i<Nrows ; i++ ) graph_decomposition[i] = aggr_index[i];

    aggr_viz_and_stats = (ML_Aggregate_Viz_Stats *) (Amatrix->to->Grid->Grid);
    aggr_viz_and_stats->graph_decomposition = graph_decomposition;
    aggr_viz_and_stats->Nlocal = Nrows;
    aggr_viz_and_stats->Naggregates = aggr_count;
    aggr_viz_and_stats->local_or_global = ML_LOCAL_INDICES;
    aggr_viz_and_stats->is_filled = ML_YES;
    aggr_viz_and_stats->Amatrix = Amatrix;
  }

  /* ********************************************************************** */
  /* take the decomposition as created by METIS and form the aggregates     */
  /* ********************************************************************** */

  total_nz = ML_Comm_GsumInt( comm, total_nz);
  i = ML_Comm_GsumInt( comm, Nrows);

  if ( mypid == 0 && 7 < ML_Get_PrintLevel())
    printf("%s Total (block) nnz = %d ( = %5.2f/(block)row)\n",
           str,
           total_nz,1.0*total_nz/i);

  if ( ml_ag->operator_complexity == 0.0 ) {
    ml_ag->fine_complexity = total_nz;
    ml_ag->operator_complexity = total_nz;
  }
  else ml_ag->operator_complexity += total_nz;

  /* fix aggr_index for num_PDE_eqns > 1 */

  for (i = Nrows - 1; i >= 0; i-- ) {
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
  exp_Nrows  *= num_PDE_eqns;

  /* count the size of each aggregate */

  aggr_cnt_array = (int *) ML_allocate(sizeof(int)*(aggr_count+1));
  for (i = 0; i < aggr_count ; i++) aggr_cnt_array[i] = 0;
  for (i = 0; i < exp_Nrows; i++) {
    if (aggr_index[i] >= 0) {
      if( aggr_index[i] >= aggr_count ) {
        fprintf( stderr,
                "*ML*WRN* on process %d, something weird happened...\n"
                "*ML*WRN* node %d belong to aggregate %d (#aggr = %d)\n"
                "*ML*WRN* (file %s, line %d)\n",
                comm->ML_mypid,
                i,
                aggr_index[i],
                aggr_count,
                __FILE__,
                __LINE__ );
      } else {
        aggr_cnt_array[aggr_index[i]]++;
      }
    }
  }

  /* ============================================================= */
  /* Form tentative prolongator                                    */
  /* ============================================================= */

  Ncoarse = aggr_count;

  /* ============================================================= */
  /* check and copy aggr_index                                     */
  /* ------------------------------------------------------------- */

  level = ml_ag->cur_level;
  nbytes = (Nrows+1) * sizeof( int );
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
  /* set up the new operator                                       */
  /* ------------------------------------------------------------- */

  new_Nrows = Nrows;
  exp_Ncoarse = Nrows;

  for ( i = 0; i < new_Nrows; i++ ) 
  {
    if ( aggr_index[i] >= exp_Ncoarse ) 
    {
      printf("*ML*WRN* index out of bound %d = %d(%d)\n",
             i, aggr_index[i], 
             exp_Ncoarse);
    }
  }
  nbytes = ( new_Nrows+1 ) * sizeof(int); 
  ML_memory_alloc((void**)&(new_ia), nbytes, "AIA");
  nbytes = ( new_Nrows+1)  * nullspace_dim * sizeof(int); 
  ML_memory_alloc((void**)&(new_ja), nbytes, "AJA");
  nbytes = ( new_Nrows+1)  * nullspace_dim * sizeof(double); 
  ML_memory_alloc((void**)&(new_val), nbytes, "AVA");
  for ( i = 0; i < new_Nrows*nullspace_dim; i++ ) new_val[i] = 0.0;

  /* ------------------------------------------------------------- */
  /* set up the space for storing the new null space               */
  /* ------------------------------------------------------------- */

  nbytes = (Ncoarse+1) * nullspace_dim * nullspace_dim * sizeof(double);
  ML_memory_alloc((void**)&(new_null),nbytes,"AGr");
  if( new_null == NULL ) {
    fprintf( stderr,
            "*ML*ERR* on process %d, not enough memory for %d bytes\n"
            "*ML*ERR* (file %s, line %d)\n",
            mypid,
            nbytes,
            __FILE__,
            __LINE__ );
    exit( EXIT_FAILURE );
  }

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
  for (i = 0; i < aggr_count; i++) {
    nbytes = aggr_cnt_array[i]+1;
    rows_in_aggs[i] = (int *) ML_allocate(nbytes*sizeof(int));
    aggr_cnt_array[i] = 0;
    if (rows_in_aggs[i] == NULL)  {
      printf("*ML*ERR* couldn't allocate memory in CoarsenMETIS\n");
      exit(1);
    }
  }
  for (i = 0; i < exp_Nrows; i+=num_PDE_eqns) {
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
      for (j = 0; j < (int) length; j++)
      {
        index = rows_in_aggs[i][j];

        for (k = 0; k < nullspace_dim; k++)
        {
          if ( unamalg_bdry[index/num_PDE_eqns] == 'T')
            qr_tmp[k*length+j] = 0.;
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
        for (j = 0; j < (int) length; j++)
        {
          index = rows_in_aggs[i][j];
          if ( unamalg_bdry[index/num_PDE_eqns] == 'T')
            qr_tmp[k*length+j] = 0.;
          else {
            if (index < Nrows) {
              qr_tmp[k*length+j] = nullspace_vect[k*Nrows+index];
            }
            else {
              fprintf( stderr,
                      "*ML*ERR* in QR\n"
                      "*ML*ERR* (file %s, line %d)\n",
                      __FILE__,
                      __LINE__ );
              exit( EXIT_FAILURE );
            }
          }
        }
      }
    }

    /* ---------------------------------------------------------- */
    /* now calculate QR using an LAPACK routine                   */
    /* ---------------------------------------------------------- */

    if (aggr_cnt_array[i] >= nullspace_dim) {

      DGEQRF_F77(&(aggr_cnt_array[i]), &nullspace_dim, qr_tmp, 
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

      if ( aggr_cnt_array[i] < nullspace_dim ){
        printf("Error in dorgqr on %d row (dims are %d, %d)\n",i,aggr_cnt_array[i],
               nullspace_dim);
        printf("ERROR : performing QR on a MxN matrix where M<N.\n");
      }
      DORGQR_F77(&(aggr_cnt_array[i]), &nullspace_dim, &nullspace_dim, 
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
          fprintf( stderr,
                  "*ML*ERR* in QR: index out of bounds (%d - %d)\n",
                  index,
                  Nrows );
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

  /* ------------------------------------------------------------- */
  /* set up the csr_data data structure                            */
  /* ------------------------------------------------------------- */

  ML_memory_alloc((void**) &csr_data, sizeof(struct ML_CSR_MSRdata),"CSR");
  csr_data->rowptr  = new_ia;
  csr_data->columns = new_ja;
  csr_data->values  = new_val;

  ML_Operator_Set_ApplyFuncData( *Pmatrix, nullspace_dim*Ncoarse, Nrows, 
                                csr_data, Nrows, NULL, 0);
  (*Pmatrix)->data_destroy = ML_CSR_MSR_ML_memorydata_Destroy;
  (*Pmatrix)->getrow->pre_comm = ML_CommInfoOP_Create();
  (*Pmatrix)->max_nz_per_row = 1;

  ML_Operator_Set_Getrow((*Pmatrix), Nrows, CSR_getrow);
  ML_Operator_Set_ApplyFunc((*Pmatrix), CSR_matvec);
  (*Pmatrix)->max_nz_per_row = 1;
  /* this must be set so that the hierarchy generation does not abort early
     in adaptive SA */
  (*Pmatrix)->num_PDEs = nullspace_dim;

  /* ------------------------------------------------------------- */
  /* clean up                                                      */
  /* ------------------------------------------------------------- */

  ML_free(unamalg_bdry);
  ML_memory_free((void**)&aggr_index);
  ML_free(aggr_cnt_array);
  for (i = 0; i < aggr_count; i++) ML_free(rows_in_aggs[i]);
  ML_memory_free((void**)&rows_in_aggs);
  ML_memory_free((void**)&qr_tmp);
  ML_memory_free((void**)&tmp_vect);
  ML_memory_free((void**)&work);

  aggr_curr = aggr_head;
  while ( aggr_curr != NULL ) 
  {
    supernode = aggr_curr;
    aggr_curr = aggr_curr->next;
    if ( supernode->length > 0 ) ML_free( supernode->list );
    ML_free( supernode );
  }

  return Ncoarse*nullspace_dim;

} /* ML_Aggregate_CoarsenUser */
