/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ************************************************************************* */
/* Author        : Marzio Sala (SNL)                                         */
/* Date          : September 2004                                            */
/* ************************************************************************* */
/* ************************************************************************* */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ml_aggregate.h"
#include "ml_lapack.h"
#include "ml_utils.h"
#include "ml_agg_Zoltan.h"
#include "ml_agg_METIS.h"
#include "ml_agg_ParMETIS.h"
#include "ml_viz_stats.h"
#include "ml_op_utils.h"

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" 
{
#endif
#endif

int ML_BuildReorderedDecomposition( int starting_decomposition[],
					   int reordered_decomposition[],
					   int Nrows, int Naggregates,
					   int nodes_per_aggre[],
					   USR_COMM comm );
extern ML_Operator * ML_BuildQ( int StartingNumElements,
				int ReorderedNumElements,
				int num_PDE_eqns, int,
				int * reordered_decomposition,
				double * StartingNullSpace,
				double * ReorderedNullSpace,
				int ComputeNewNullSpace,
				double * StartingBdry, double * ReorderedBdry,
				USR_COMM mpi_communicator,
				ML_Comm *ml_communicator );

extern int ML_ApplyQ(int StartingNumElements,
		     int ReorderedNumElements,
		     int NumVectors,
		     double* StartingVectors,
		     double* ReorderedVectors);

extern void ML_DestroyQ(void);

extern int OPTIMAL_VALUE;
extern int PARMETIS_DEBUG_LEVEL;
extern int OPTIMAL_LOCAL_COARSE_SIZE;

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#if defined(HAVE_ML_ZOLTAN) && defined(HAVE_MPI)
#include "zoltan.h"

/*
 *  PROTOTYPES for load-balancer interface functions.
 */

ZOLTAN_NUM_OBJ_FN ML_get_num_entries;
ZOLTAN_OBJ_LIST_FN ML_get_entries;

ZOLTAN_NUM_GEOM_FN ML_get_num_geom;
ZOLTAN_GEOM_MULTI_FN ML_get_geom_multi;

#ifdef USE_GRAPH
/* These two functions are needed to use ParMETIS thru Zoltan. */
/* They are not needed for geometric methods. */
ZOLTAN_NUM_EDGES_MULTI_FN ML_get_num_edges_multi;
ZOLTAN_EDGE_LIST_MULTI_FN ML_get_edge_list_multi;
#endif

/* Hypergraph functions */
#ifdef HAVE_ML_ZOLTAN_THREE
ZOLTAN_HG_SIZE_CS_FN ML_zoltan_hg_size_cs_fn;
ZOLTAN_HG_CS_FN ML_zoltan_hg_cs_fn;
ZOLTAN_OBJ_SIZE_MULTI_FN ML_zoltan_obj_size_multi_fn;
#endif

/* "...crappy code..."
 *
 * I use a bunch of static variables, with suffix `MLZ':
 * - MLZ_dim specifies the number of dimensions (1,2 or 3);
 * - MLZ_x contains the x-coord;
 * - MLZ_y contains the y-coord;
 * - MLZ_y contains the z-coord;
 * - ML_offset contains the offest of local row i. Global
 *   numbering of local row i is simply i + offset
 */
static int MLZ_dim;
static double* MLZ_x;
static double* MLZ_y;
static double* MLZ_z;
static int MLZ_offset;

static int* MLZ_indices;


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int setup_zoltan(struct Zoltan_Struct *zz, ML_Operator* A, int zoltan_type, int zoltan_estimated_its,
                        int zoltan_timers, int smoothing_steps, int rows_per_amalgamated_row)
{
  /* Fix the output level */
  char str[80];
  if (ML_Get_PrintLevel() > 9) strcpy(str,"1");
  else strcpy(str,"0");   

  if (Zoltan_Set_Param(zz, "DEBUG_LEVEL", str) == ZOLTAN_FATAL) {
    printf("fatal(0)  error returned from Zoltan_Set_Param(LB_METHOD)\n");
    return 0;
  }
  
  /* Set the load-balance method */
  if(zoltan_type == ML_ZOLTAN_TYPE_RCB) strcpy(str,"RCB");
#ifdef HAVE_ML_ZOLTAN_THREE  
  else if(zoltan_type == ML_ZOLTAN_TYPE_HYPERGRAPH) strcpy(str,"hypergraph");
  else if(zoltan_type == ML_ZOLTAN_TYPE_FAST_HYPERGRAPH) strcpy(str,"hypergraph");  
  else {
    printf("fatal(1) no valid zoltan type set\n");
    return 0;
  }
#else
  if(!A->comm->ML_mypid && zoltan_type != ML_ZOLTAN_TYPE_RCB) printf("ML-Zoltan: Zoltan 3.0 support not enabled, resetting method to RCB\n");
  strcpy(str,"RCB");
#endif
  
  if (Zoltan_Set_Param(zz, "LB_METHOD", str) == ZOLTAN_FATAL) {
    printf("fatal(1)  error returned from Zoltan_Set_Param(LB_METHOD)\n");
    return 0;
  }  

  /* Zoltan Timers & Output */
  if(zoltan_timers == 1){
    if (Zoltan_Set_Param(zz, "USE_TIMERS", "0") == ZOLTAN_FATAL) {
      printf("fatal(1)  error returned from Zoltan_Set_Param(USE_TIMERS)\n");
      return 0;
    }  
    if (Zoltan_Set_Param(zz, "FINAL_OUTPUT", "1") == ZOLTAN_FATAL) {
      printf("fatal(1)  error returned from Zoltan_Set_Param(FINAL_OUTPUT)\n");
      return 0;
    }  
  }


  /* Hypergraph parameters */
#ifdef HAVE_ML_ZOLTAN_THREE
  if(zoltan_type == ML_ZOLTAN_TYPE_HYPERGRAPH){    
    /* Set the repartitioning flag */
    if (Zoltan_Set_Param(zz, "LB_APPROACH", "repartition") == ZOLTAN_FATAL) {
      printf("fatal(1)  error returned from Zoltan_Set_Param(LB_APPROACH)\n");
      return 0;
    }  
  }
  else if (zoltan_type == ML_ZOLTAN_TYPE_FAST_HYPERGRAPH){
    /* Set the fast repartitioning flag */
    if (Zoltan_Set_Param(zz, "LB_APPROACH", "fast_repartition") == ZOLTAN_FATAL) {
      printf("fatal(1)  error returned from Zoltan_Set_Param(LB_APPROACH)\n");
      return 0;
    }  
  }
  
  if(zoltan_type == ML_ZOLTAN_TYPE_HYPERGRAPH || zoltan_type == ML_ZOLTAN_TYPE_FAST_HYPERGRAPH){
    /* Set the repartitioning multiplier (aka alpha) */    
    sprintf(str,"%d",(int)(zoltan_estimated_its*rows_per_amalgamated_row*(smoothing_steps+1)*sizeof(double)));
    if (Zoltan_Set_Param(zz, "PHG_REPART_MULTIPLIER",str) == ZOLTAN_FATAL) {
      printf("fatal(1)  error returned from Zoltan_Set_Param(PHG_REPART_MULTIPLIER)\n");
      return 0;
    }
  }

#endif
  
  /*
   * Set the callback functions
   */

  if (Zoltan_Set_Fn(zz, ZOLTAN_NUM_OBJ_FN_TYPE, (void (*)()) ML_get_num_entries,
                    (void *) A) == ZOLTAN_FATAL) {
    printf("fatal(2)  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }

  if (Zoltan_Set_Fn(zz, ZOLTAN_OBJ_LIST_FN_TYPE,
                    (void (*)()) ML_get_entries,
                    (void *) A) == ZOLTAN_FATAL) {
    printf("fatal(3)  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }

  /* Functions for geometry based algorithms */
  if (Zoltan_Set_Fn(zz, ZOLTAN_NUM_GEOM_FN_TYPE, (void (*)()) ML_get_num_geom,
                    (void *) A) == ZOLTAN_FATAL) {
    printf("fatal(4)  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }

  if (Zoltan_Set_Fn(zz, ZOLTAN_GEOM_MULTI_FN_TYPE,
                    (void (*)()) ML_get_geom_multi,
                    (void *) A) == ZOLTAN_FATAL) {
    printf("fatal(5)  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }

#ifdef USE_GRAPH
  /* Functions for graph based algorithms */
  /* These two functions are needed to use ParMETIS thru Zoltan. */
  /* They are not needed for geometric methods. */
  if (Zoltan_Set_Fn(zz, ZOLTAN_NUM_EDGES_MULTI_FN_TYPE,
                    (void (*)()) ML_get_num_edges_multi,
                    (void *) A) == ZOLTAN_FATAL) {
    printf("fatal(6)  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }
  if (Zoltan_Set_Fn(zz, ZOLTAN_EDGE_LIST_MULTI_FN_TYPE,
                    (void (*)()) ML_get_edge_list_multi,
                    (void *) A) == ZOLTAN_FATAL) {
    printf("fatal(7)  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }
#endif /* USE_GRAPH */

  /* Hypergraph Functions */
#ifdef HAVE_ML_ZOLTAN_THREE
  if(zoltan_type == ML_ZOLTAN_TYPE_HYPERGRAPH || zoltan_type == ML_ZOLTAN_TYPE_FAST_HYPERGRAPH){
    if(Zoltan_Set_Fn(zz, ZOLTAN_HG_SIZE_CS_FN_TYPE,
                      (void (*)()) ML_zoltan_hg_size_cs_fn,
                      (void *) A) == ZOLTAN_FATAL) {
      printf("fatal(8)  error returned from Zoltan_Set_Fn()\n");
      return 0;
    }
    if(Zoltan_Set_Fn(zz, ZOLTAN_HG_CS_FN_TYPE,
                      (void (*)()) ML_zoltan_hg_cs_fn,
                      (void *) A) == ZOLTAN_FATAL) {
      printf("fatal(9)  error returned from Zoltan_Set_Fn()\n");
      return 0;
    }
    if(Zoltan_Set_Fn(zz, ZOLTAN_OBJ_SIZE_MULTI_FN_TYPE,
                     (void (*)()) ML_zoltan_obj_size_multi_fn,
                     (void *) A) == ZOLTAN_FATAL) {
      printf("fatal(10)  error returned from Zoltan_Set_Fn()\n");
      return 0;
    }       
  }
#endif
  return(1);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int run_zoltan(int N_parts, struct Zoltan_Struct *zz, ML_Operator* A,
                      int Nrows, int graph_decomposition[],
                      int MyPID)
{
  /* Variables returned by Zoltan */
  ZOLTAN_ID_PTR import_gids = NULL;  /* Global nums of objs to be imported   */
  ZOLTAN_ID_PTR import_lids = NULL;  /* Local indices to objs to be imported */
  int   *import_procs = NULL;        /* Proc IDs of procs owning objs to be
                                        imported.                            */
  int   *import_to_part = NULL;      /* Partition #s to which imported objs 
                                        should be assigned.                  */
  ZOLTAN_ID_PTR export_gids = NULL;  /* Global nums of objs to be exported   */
  ZOLTAN_ID_PTR export_lids = NULL;  /* local indices to objs to be exported */
  int   *export_procs = NULL;        /* Proc IDs of destination procs for objs
                                        to be exported.                      */
  int   *export_to_part = NULL;      /* Partition #s for objs to be exported.*/
  int num_imported;              /* Number of objs to be imported.          */
  int num_exported;              /* Number of objs to be exported.          */
  int new_decomp;                /* Flag indicating whether the decomposition
                                    has changed                              */
  int num_gid_entries;           /* Number of array entries in a global ID.  */
  int num_lid_entries;           /* Number of array entries in a local ID.   */

  int i;
  char value[80];

  /* sets the number of partitions */
  sprintf(value,"%d",N_parts);
  Zoltan_Set_Param(zz,"num_global_partitions", value);

  /*
   * Call Zoltan to compute a new decomposition. 
   * The result is given as export and/or import lists.
   * If you need only the export lists, perform the 
   * following call before this function:
   *   Zoltan_Set_Param(zz, "RETURN_LISTS", "EXPORT");
   *
   * In this example, we ignore the import_to_part/export_to_part
   * arguments, as they are only useful if we distinguish between
   * partitions and processors.
   */
  if (Zoltan_LB_Partition(zz, &new_decomp, &num_gid_entries, &num_lid_entries,
               &num_imported, &import_gids,
               &import_lids, &import_procs, &import_to_part,
               &num_exported, &export_gids,
               &export_lids, &export_procs, &export_to_part) == ZOLTAN_FATAL){
    printf("fatal(8)  error returned from Zoltan_LB_Partition()\n");
    return 0;
  }

  for (i = 0 ; i < Nrows ; ++i)
    graph_decomposition[i] = MyPID;
  if (1 || new_decomp){
    for (i = 0 ; i < num_exported ; ++i) {
#ifdef DEBUG
      assert (export_gids[i] - MLZ_offset >= 0);
      assert (export_gids[i] - MLZ_offset < Nrows);
      assert (export_to_part[i] >= 0);
      assert (export_to_part[i] < value);
#endif
      graph_decomposition[export_gids[i] - MLZ_offset] = export_to_part[i];
    }
#if 0
  for (i = 0 ; i < Nrows ; ++i)
    printf("graph_decomposition[%d] = %d\n", i, graph_decomposition[i]);
#endif

#if 0
    printf("[Proc %1d] My data to export are:\n", myrank);
    for (k=0; k<num_exported; k++){
      ptr = (int *) &export_gids[num_gid_entries*k];
      printf("[Proc %1d] Export (%d,%d) to proc %d\n", myrank, 
         *ptr, *(ptr+1), export_procs[k]);
    }
#endif
  }

  /* Clean up */
  Zoltan_LB_Free_Part(&import_gids, &import_lids,
                      &import_procs, &import_to_part);
  Zoltan_LB_Free_Part(&export_gids, &export_lids,
                      &export_procs, &export_to_part);

  return(1);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int ML_get_num_entries (void *data, int *ierr)
{

  ML_Operator *A;

  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return 0;
  }
  /* cast the data pointer to correct data type */
  A = (ML_Operator *) data;

  *ierr = ZOLTAN_OK; /* set error code */

  return(A->getrow->Nrows);  /* # local rows */
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void ML_get_entries(void *data, int num_gid_entries, int num_lid_entries,
                 ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
                 int wdim, float *wgt, int *ierr)
{

  ML_Operator* A;
  int k;
  int allocated, *bindx, row_length;
  double *val;

  *ierr = ZOLTAN_OK; 

  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return;
  }
  A = (ML_Operator *) data;

  /* We should be using (at least) one int for each GID. */
  if (num_gid_entries < 1) {
    *ierr = ZOLTAN_FATAL;
    return;
  }
  allocated = 100;
  bindx = (int    *)  ML_allocate( allocated*sizeof(int   ));
  val   = (double *)  ML_allocate( allocated*sizeof(double));

  A = (ML_Operator*) data;
  for (k = 0; k < A->getrow->Nrows; k++) {
    global_id[k] = (ZOLTAN_ID_TYPE) (k + MLZ_offset);

    /* Add (optional) local ids and/or weights here if desired.  */

    ML_get_matrix_row(A, 1, &k, &allocated, &bindx, &val,
                        &row_length, 0);
    wgt[k] = (float) row_length;
    
  }
  ML_free(bindx);
  ML_free(val);


}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int ML_get_num_geom(void *data, int *ierr)
{
  *ierr = ZOLTAN_OK; /* set error flag */
  return(MLZ_dim);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void ML_get_geom_multi(void *data, int num_gid_entries, int num_lid_entries,
              int num_obj, ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
              int num_dim, double *coor, int *ierr)
{
  ML_Operator* A; 
  int k;

  *ierr = ZOLTAN_OK;

  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return;
  }

  A = (ML_Operator*) data;
  for (k = 0; k < A->getrow->Nrows; ++k) {
    coor[MLZ_dim * k]   = (double) MLZ_x[k];
    if (MLZ_dim > 1) {
      coor[MLZ_dim * k +1] = (double) MLZ_y[k];
      if (MLZ_dim > 2)
        coor[MLZ_dim * k + 2] = (double) MLZ_z[k];
    }
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#ifdef HAVE_ML_ZOLTAN_THREE  
void ML_zoltan_hg_size_cs_fn(void *data, int *num_lists, int *num_pins, int *format, int *ierr){
  ML_Operator* A; 
  *ierr = ZOLTAN_OK;
  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return;
  }

  A = (ML_Operator*) data;
  *num_lists = A->getrow->Nrows;             /* Local # of rows */
  *num_pins  = ML_Operator_ComputeNumNzs(A); /* Local nnz       */
  *format    = ZOLTAN_COMPRESSED_VERTEX;
}


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void ML_zoltan_hg_cs_fn(void *data, int num_gid_entries, int num_vtx_edge, int num_pins, int format, ZOLTAN_ID_PTR vtxedge_GID, int *vtxedge_ptr, ZOLTAN_ID_PTR pin_GID, int *ierr){
  ML_Operator* A;
  struct ML_CSR_MSRdata * CSR_Data;
  double *values;
  int i,N,maxnz,rv,rowlength,rowtotal=0;
  int *pin_ids,j;

  *ierr = ZOLTAN_OK;
  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return;
  }
  
  A = (ML_Operator*) data;
  N= A->getrow->Nrows; 

  /* Use getrow to pull all the data over */
  values =(double*) ML_allocate(A->max_nz_per_row*sizeof(double));  
  pin_ids =(int*) ML_allocate(A->max_nz_per_row*sizeof(int));  

  for(i=0;i<N;i++) {
    vtxedge_GID[i] = (ZOLTAN_ID_TYPE) (i + MLZ_offset);   
    vtxedge_ptr[i]=rowtotal;
    rv=(*A->getrow->func_ptr)(A,1,&i,A->max_nz_per_row,pin_ids,values,&rowlength);
    if(rv==0) {printf("ML: Out of space in getrow i=%d/%d",i,N);fflush(stdout);*ierr=ZOLTAN_FATAL;}
    for(j=0;j<rowlength;j++) pin_GID[rowtotal+j] = (ZOLTAN_ID_TYPE) pin_ids[j];
    /* ADD: Safety catch if A->max_nz_per_row is not correct.  Use ML_get_matrix_row and do a copy-back. */
    rowtotal+=rowlength;
  }
  /* Note: vtxedge_ptr does not have the extra "N+1" element that regular CSR
     matrices have. */
  
  /* Reindex to global IDs */
  for(i=0;i<num_pins;i++) pin_GID[i]=MLZ_indices[pin_GID[i]];    
  
  ML_free(values);
  ML_free(pin_ids);
}


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void ML_zoltan_obj_size_multi_fn(void * data,int num_gid_entries,int num_lid_entries,int num_ids,
                                 ZOLTAN_ID_PTR global_ids,ZOLTAN_ID_PTR local_ids,int *sizes,int *ierr){
  ML_Operator* A;
  struct ML_CSR_MSRdata * CSR_Data;
  int *indices;
  double *values;
  int i,N,maxnz,rowlength,rv;
  struct ML_CSR_MSRdata *input_matrix;
  
  *ierr = ZOLTAN_OK;
  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return;
  }
  
  A = (ML_Operator*) data;
  N= A->getrow->Nrows; 
  /* Get the row lengths --- use special purpose code for MSR matrices,
     just call the getrow routine for everything else */
  /*  if(A->getrow->func_ptr==MSR_getrows){
      input_matrix = (struct ML_CSR_MSRdata *) ML_Get_MyGetrowData(A);
      for(i=0;i<N;i++) 
      sizes[i] = (input_matrix->columns[i+1] - input_matrix->columns[i] + 1)* (sizeof(double)+sizeof(int));
      }
      else{*/
    indices =(int*) ML_allocate(A->max_nz_per_row*sizeof(int));  
    values  =(double*) ML_allocate(A->max_nz_per_row*sizeof(double));  
    for(i=0;i<N;i++) {
      rv=(*A->getrow->func_ptr)(A,1,&i,A->max_nz_per_row,indices,values,&rowlength);
      if(rv==0) {printf("ML: Out of space in getrow i=%d/%d",i,N);fflush(stdout);*ierr=ZOLTAN_FATAL;}
      sizes[i]=rowlength*(sizeof(double)+sizeof(int)); /* add its cost here */
    }
    ML_free(indices);
    ML_free(values);
    /*  }*/
}
#endif

#endif

/* ======================================================================== */
/*!
 \brief Reorder the local graph for the coarser level matrix using a Zoltan.

 \param Amatrix (In) -
 pointer to the amalgamated ML_Operator to decompose

 \param N_parts (In) -
 total number of parts in which the graph has to be decomposed

 \param graph_decomposition (Out) -
 integer vector, of size Amatrix->getrow->Nrows, that in output
 will contain the global partition number of each local row.

 \param bdry_nodes (In) -

 \param N_nonzeros () -

 \param current_level (In) -
 
*/
/* ------------------------------------------------------------------------ */

int ML_DecomposeGraph_with_Zoltan(ML_Operator *Amatrix,
				  int N_parts,
				  int graph_decomposition[],
				  double bdry_nodes[], double old_x[], 
				  double old_y[], double old_z[],
				  int current_level,
                                  int zoltan_type, int zoltan_estimated_its,
				  int zoltan_timers,
                                  int smoothing_steps, int rows_per_amalgamated_row){

  int i, Nrows;
#if defined(HAVE_ML_ZOLTAN) && defined(HAVE_MPI)
  ML_Comm * comm = Amatrix->comm;
  int mypid = Amatrix->comm->ML_mypid;
  double t0;
  struct Zoltan_Struct *zz;
  float version;
  int error;
  MPI_Comm Zoltan_SubComm = comm->USR_comm;
  int color = MPI_UNDEFINED; 

  /* ------------------- execution begins --------------------------------- */

  t0 = GetClock();
  
  Nrows = Amatrix->getrow->Nrows;
  
  /* ********************************************************************** */
  /* some general variables                                                 */
  /* ********************************************************************** */
  
  if (old_z)
    MLZ_dim = 3;
  else if (old_y)
    MLZ_dim = 2;
  else if (old_x) 
    MLZ_dim = 1;
  else 
    MLZ_dim = 0;

  MLZ_dim = ML_Comm_GmaxInt(Amatrix->comm, MLZ_dim);

  if (MLZ_dim == 0) {
    if (Amatrix->comm->ML_mypid == 0 && ML_Get_PrintLevel() > 0) 
      printf("ML*WRN* ML_DecomposeGraph_with_Zoltan: No coordinates given.\nML*WRN* Giving up on repartitioning.\n");
    return(-1);
  }

  /* FIXME */
  for (i = 0 ; i < Nrows * 0; ++i) {
    printf("[%d] = %e %e\n", i, old_x[i], old_y[i]);
  }
  /* junk */
  MLZ_x = old_x;
  MLZ_y = old_y;
  MLZ_z = old_z;

#ifdef ML_MPI
  MPI_Scan(&Nrows, &MLZ_offset, 1, MPI_INT, MPI_SUM, comm->USR_comm);
  MLZ_offset -= Nrows;
#else
  MLZ_offset = 0;
#endif

  /* Build a Global Column Numbering --- needed for the HG partitioner */
  ML_build_global_numbering(Amatrix, &MLZ_indices,"cols");

  if (N_parts <= 0) N_parts = 1;

  /* ********************************************************************** */
  /* no need to call Zoltan if only one global aggregate is required.       */
  /* ********************************************************************** */

  if (N_parts == 1) {
    for (i = 0 ; i < Nrows ; i++) {
      graph_decomposition[i] = 0;
    }
    return 1;
  }


  /* Define Subcommunicator */
  if (Nrows > 0) color = 0;
  MPI_Comm_split(comm->USR_comm,color,comm->ML_mypid,&Zoltan_SubComm);  

  if (Nrows > 0){
    /* BEGIN OF ZOLTAN CALL */
    
    /*  Initialize Zoltan. It will start MPI if we haven't already. */
    /*  Do this only once. */

    if ((error = Zoltan_Initialize(0, NULL, &version)) != ZOLTAN_OK) {
      printf("fatal(10) Zoltan_Initialize returned error code, %d", error);
      goto End;
    }



    /*
     *  Create a Zoltan structure.
     */
    if ((zz = Zoltan_Create(Zoltan_SubComm)) == NULL) {    
      printf("fatal(11)  NULL returned from Zoltan_Create()\n");
      goto End;
    }

    /*
     *  Tell Zoltan what kind of local/global IDs we will use.
     *  In our case, each GID is two ints and there are no local ids.
     *  One can skip this step if the IDs are just single ints.
     */
    Zoltan_Set_Param(zz, "num_gid_entries", "1");
    Zoltan_Set_Param(zz, "num_lid_entries", "0");
    Zoltan_Set_Param(zz, "obj_weight_dim", "1");
    
    /*
     *  Set up Zoltan query functions for our Matrix data structure.
     */
    
    if (!setup_zoltan(zz, Amatrix,zoltan_type, zoltan_estimated_its,zoltan_timers,
		      smoothing_steps, rows_per_amalgamated_row)){
      printf("fatal(12) Error returned from setup_zoltan\n");
      goto End;
    }
  
    /*
     * Run Zoltan to compute a new load balance.
     * Data migration may also happen here.
     */
    /* Uncomment the next line to produce debugging information.*/
    
    /*Zoltan_Generate_Files(zz, "ZoltanDebugging", 1, 1, 0, 0);*/
    if (!run_zoltan(N_parts, zz, Amatrix, Nrows, graph_decomposition,
		    comm->ML_mypid)) {
      printf("fatal(13) Error returned from run_zoltan\n");
      goto End;
    }
    
    
  End:
    /* Destroy Zoltan structure */
    Zoltan_Destroy(&zz);
    
    /* END OF ZOLTAN CALL */
    /* ------------------- that's all folks --------------------------------- */
    
  }
  
  t0 = GetClock() - t0;

  if ( mypid == 0 &&  ML_Get_PrintLevel() > 7 ) {
   
    printf("Zoltan (level %d) : time required = %e\n",
	   current_level,
	   t0 );
  }

  
  /* Cleanup */
  free(MLZ_indices);

  if(Nrows > 0)
    MPI_Comm_free(&Zoltan_SubComm);
  
  /* returns the *global* number of partitions */
  return(N_parts);
#else

  puts("*ML*ERR* You must configure ml with Zoltan support, using");
  puts("*ML*ERR* parameter --with-ml_zoltan in your configuration line.");
  puts("*ML*ERR* You also need --enable-mpi to use Zoltan");
  puts("*ML*ERR* Now inserting all local nodes in the same aggregate...");

  Nrows = Amatrix->getrow->Nrows;
  for (i = 0 ; i < Nrows ; ++i)
    graph_decomposition[i] = 0;

  return(1);
#endif
  
} /* ML_DecomposeGraph_with_Zoltan */

/* ======================================================================== */
/*!
 \brief create non-smoothed aggregates using Zoltan. 

*/
/* ------------------------------------------------------------------------ */

int ML_Aggregate_CoarsenZoltan(ML_Aggregate *ml_ag, ML_Operator *Amatrix, 
			       ML_Operator **Pmatrix, ML_Comm *comm)
{

   unsigned int nbytes, length;
   int     i, j, jj, k, Nrows, exp_Nrows,  N_bdry_nodes;
   int     diff_level, Nrows_global;
   int     aggr_count, index, mypid, num_PDE_eqns;
   int     *aggr_index = NULL, nullspace_dim;
   int     Ncoarse, count;
   int     *new_ia = NULL, *new_ja = NULL, new_Nrows;
   int     *aggr_cnt_array = NULL;
   int     level, index3, max_agg_size;
   int     **rows_in_aggs = NULL, lwork, info;
   double  *new_val = NULL, epsilon;
   double  *nullspace_vect = NULL, *qr_tmp = NULL;
   double  *tmp_vect = NULL, *work = NULL, *new_null = NULL;
   ML_SuperNode          *aggr_head = NULL, *aggr_curr, *supernode;
   struct ML_CSR_MSRdata *csr_data;
   double                *starting_amalg_bdry, *reordered_amalg_bdry;
   int                   Nghost;
   int                   allocated = 0, *rowi_col = NULL, rowi_N;
   double                *rowi_val = NULL;
   double Nnonzeros2 = 0;
   int optimal_value;
   ML_Operator * Pmatrix2 = NULL;
   ML_Aggregate_Viz_Stats *grid_info;
   
#if defined(OUTPUT_AGGREGATES) || defined(DDEBUG) || defined(INPUT_AGGREGATES) || (ML_AGGR_INAGGR) || (ML_AGGR_OUTAGGR) || (ML_AGGR_MARKINAGGR)
   FILE *fp;
   char fname[80];
   static int level_count = 0;
   double *d2temp;
   int agg_offset, vertex_offset;
#endif
   ML_Aggregate_Viz_Stats * aggr_viz_and_stats;
   ML_Aggregate_Options * aggr_options;
   int Nprocs;
   int * starting_offset = NULL, * reordered_offset = NULL;
   int desired_aggre_per_proc;
   int * nodes_per_aggre = NULL;
   int * starting_decomposition = NULL;
   int * reordered_decomposition = NULL;
   ML_Operator * QQ = NULL;
   ML_Operator *Pstart = NULL;
   int starting_aggr_count;
   char str[80];
   double * new_nullspace_vect = NULL;
   int * graph_decomposition = NULL;
   int N_dimensions = 0;
   double* old_x = NULL;
   double* old_y = NULL;
   double* old_z = NULL;
   
   /* ------------------- execution begins --------------------------------- */

   sprintf( str, "Zoltan (level %d) :", ml_ag->cur_level );
   
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

   if( mypid == 0 && 5 < ML_Get_PrintLevel() ) {
     printf("%s num_PDE_eqns = %d\n",
	    str,
	    num_PDE_eqns);
   }

#ifdef ML_MPI
   MPI_Allreduce( &Nrows, &Nrows_global, 1, MPI_INT, MPI_SUM, Amatrix->comm->USR_comm );
#else
   Nrows_global = Nrows;
#endif
   
   /* ============================================================= */
   /* check the system size versus null dimension size              */
   /* ============================================================= */
     
   if ( Nrows % num_PDE_eqns != 0 )
   {
      printf("*ML*ERR* : Nrows must be multiples of num_PDE_eqns.\n");
      exit(EXIT_FAILURE);
   }
   diff_level = ml_ag->max_levels - ml_ag->cur_level - 1;
   if ( diff_level > 0 ) num_PDE_eqns = nullspace_dim; /* ## 12/20/99 */

   Nghost = Amatrix->getrow->pre_comm->total_rcv_length;
   
   /* ============================================================= */
   /* set up the threshold for weight-based coarsening              */
   /* ============================================================= */

   diff_level = ml_ag->begin_level - ml_ag->cur_level;
   if ( diff_level == 0 ) ml_ag->curr_threshold = ml_ag->threshold;
   epsilon = ml_ag->curr_threshold;
   ml_ag->curr_threshold *= 0.5;

   if ( mypid == 0 && 8 < ML_Get_PrintLevel() )
   {
      printf("%s current eps = %e\n",
	     str,
	     epsilon);
      if( epsilon != 0.0 ) {
	fprintf( stderr,
		 "WARNING: Zoltan may not work with dropping!\n"
		 "WARNING: Now proceeding -- with fingers crossed\n" );
      }
   }
   
   ML_Operator_AmalgamateAndDropWeak(Amatrix, num_PDE_eqns, epsilon);
   Nrows /= num_PDE_eqns;
   Nrows_global /= num_PDE_eqns;
   Nghost /= num_PDE_eqns;
   exp_Nrows = Nrows+Nghost;

   /* record the Dirichlet boundary. unamalg_bdry[i] == 1 ==> the nodes is */
   /* a boundary node, unamalg_bdry[i] == 0 ==> the node is not            */

   nbytes = sizeof(double)*(exp_Nrows + 1);
   starting_amalg_bdry = (double *) ML_allocate(nbytes);
   if( starting_amalg_bdry == NULL ) {
     fprintf( stderr,
	      "*ML*ERR* Not enough space to allocate %d bytes\n"
	      "*ML*ERR* (file %s, line %d)\n",
	      nbytes,
	      __FILE__,
	      __LINE__ );
     exit( EXIT_FAILURE );
   }
   for (i = Nrows ; i < exp_Nrows; i++) starting_amalg_bdry[i] = 0.0;

   Nnonzeros2 = 0, N_bdry_nodes = 0;
   for (i = 0; i < Nrows; i++) {
      ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
                        &rowi_N, 0);

      if (rowi_N > 1) {
        starting_amalg_bdry[i] = 0.0;
        Nnonzeros2 += rowi_N;
      } else {
	starting_amalg_bdry[i] = 1.0;
	N_bdry_nodes++;
      }
   }

   if( rowi_col != NULL ) ML_free(rowi_col );
   if( rowi_val != NULL ) ML_free(rowi_val );
   
   i = ML_Comm_GsumInt(comm, N_bdry_nodes);
   
   if( mypid == 0 && 5 < ML_Get_PrintLevel() ) {
     printf("%s # bdry (block) nodes = %d, # (block) nodes = %d\n",
	    str,
	    i, Nrows_global);
   }
   
   /* communicate the boundary information */
   
   ML_exchange_bdry(starting_amalg_bdry,Amatrix->getrow->pre_comm,
		    Nrows,comm, ML_OVERWRITE,NULL);

   /* ********************************************************************** */
   /* allocate memory for starting_decomposition and call Zoltan to          */
   /* decompose the local                                                    */
   /* graph into the number of parts specified by the user with a call       */
   /* ML_Aggregate_Set_LocalNumber( ml, ag, level, Nparts)                   */
   /* ********************************************************************** */

   /* FIXME: only Nrows??? */
   nbytes = (Nrows+Nghost) * sizeof(int);

   if ( nbytes > 0 ) starting_decomposition = (int *)ML_allocate(nbytes);
   else              starting_decomposition = NULL;
   
   if( starting_decomposition == NULL && nbytes > 0 ) {
     
     fprintf( stderr,
	      "*ML*ERR* not enough memory to allocated %d bytes\n"
	      "*ML*ERR* (file %s, line %d)\n",
	      nbytes,
	      __FILE__,
	      __LINE__ );
     exit( EXIT_FAILURE );
   }

   /* ********************************************************************** */
   /* Retrive the user's defined data. If not set, pick up default settings  */
   /* ********************************************************************** */
     
   aggr_options = (ML_Aggregate_Options *) ml_ag->aggr_options;
   
   if( aggr_options == NULL ) {

     if( mypid == 0 && 8 < ML_Get_PrintLevel() ) {
       printf("%s Using default values\n",
	      str );
     }
    
     optimal_value = ML_Aggregate_Get_OptimalNumberOfNodesPerAggregate();
     starting_aggr_count = (int)1.0*Nrows_global/optimal_value;
     if( starting_aggr_count < 1 ) starting_aggr_count = 1;
     /* 'sto xxxcio e` un po` piu` difficile, non ne ho idea, piglio
	a caso... */
     desired_aggre_per_proc = ML_max( 128, Nrows );
       
   } else {

     /* ******************************************************************** */
     /* Retrive the user's defined choice to define the number of aggregates */
     /* For local number, it is ok.                                          */
     /* If global number of aggregates or nodes per aggregate have been      */
     /* specified, compute the local one (evenly dividing this global number)*/
     /* For those two latter cases, I suppose that the input value is the    */
     /* same on all processes (but I don't check ... )                       */
     /* ******************************************************************** */

     switch( aggr_options[ml_ag->cur_level].choice ) {

     case ML_NUM_LOCAL_AGGREGATES:
       i = aggr_options[ml_ag->cur_level].Naggregates_local;
       if( mypid == 0 && 7 < ML_Get_PrintLevel() ) {
	 printf( "%s Requested %d local aggregates (on proc 0)\n",
		 str,
		 i );
       }
#ifdef ML_MPI
       MPI_Allreduce(&i,&starting_aggr_count,1,MPI_INT,MPI_SUM,comm->USR_comm);
#else
       starting_aggr_count = i;
#endif
       break;
       
     case ML_NUM_GLOBAL_AGGREGATES:
       
       starting_aggr_count = aggr_options[ml_ag->cur_level].Naggregates_global;
       if( mypid == 0 && 7 < ML_Get_PrintLevel() ) {
	 printf( "%s Requested %d global aggregates\n",
		 str,
		 starting_aggr_count );
       }

       break;
       
     case ML_NUM_NODES_PER_AGGREGATE:

       starting_aggr_count = aggr_options[ml_ag->cur_level].Nnodes_per_aggregate;

       if( mypid == 0 && 7 < ML_Get_PrintLevel() ) {
	 printf( "%s Requested %d nodes per aggregate\n",
		 str,
		 starting_aggr_count );
       }
       
       if( starting_aggr_count >= Nrows_global) {

	 i = starting_aggr_count;

	 starting_aggr_count = Nrows/OPTIMAL_VALUE;
	 if( starting_aggr_count == 0) starting_aggr_count = 1;
	 
	 if( mypid == 0 && 7 < ML_Get_PrintLevel() ) {
	   fprintf( stderr,
		    "%s WARNING : # nodes per aggregate (%d) > # nodes (%d)\n"
		    "%s WARNING : now proceeding with # aggregates = %d\n",
		    str,
		    i,
		    Nrows,
		    str,
		    starting_aggr_count);
	 }

       } else {

	 starting_aggr_count = Nrows_global/starting_aggr_count;
	 if( starting_aggr_count == 0 ) starting_aggr_count = 1;

       }
       
       break;
       
     } /* switch */

     /* reorder_flag = aggr_options[ml_ag->cur_level].reordering_flag; */

     desired_aggre_per_proc = aggr_options[ml_ag->cur_level].desired_aggre_per_proc;

     if( desired_aggre_per_proc <= 0 )
       desired_aggre_per_proc = OPTIMAL_LOCAL_COARSE_SIZE;
     
   } /* if( aggr_options == NULL )*/

   if( mypid == 0 && 5 < ML_Get_PrintLevel() ) {
     printf("%s Objective : %d aggregates over %d nodes\n",
	    str,
	    starting_aggr_count,
	    Nrows_global);
     printf("%s Objective : %d aggregates on each process\n",
	    str,
	    desired_aggre_per_proc );
   } 
   
   grid_info = (ML_Aggregate_Viz_Stats *) Amatrix->to->Grid->Grid;
   N_dimensions = grid_info->Ndim;
   old_x = grid_info->x;
   old_y = grid_info->y;
   old_z = grid_info->z;

   if (N_dimensions >= 1 && old_x == NULL)
   {
     fprintf(stderr,
             "*ML*ERR* null x pointer!\n"
             "*ML*ERR* (file %s, line %d)\n",
             __FILE__, __LINE__);
     exit(EXIT_FAILURE);
   }

   if (N_dimensions >= 2 && old_y == NULL)
   {
     fprintf(stderr,
             "*ML*ERR* null y pointer!\n"
             "*ML*ERR* (file %s, line %d)\n",
             __FILE__, __LINE__);
     exit(EXIT_FAILURE);
   }
   if (N_dimensions == 3 && old_z == NULL)
   {
     fprintf(stderr,
             "*ML*ERR* null z pointer!\n"
             "*ML*ERR* (file %s, line %d)\n",
             __FILE__, __LINE__);
     exit(EXIT_FAILURE);
   }

   /* Amatrix is the *amalgamated* matrix.
    *
    * starting_aggr_count is the number of *global* partition
    * I would like to create. The function returns the number
    * that it could actually create.
    *
    * starting_decomposition[i] will contain the *global*
    * partition ID for *local* row i. `local' means that it still
    * refers to the non-reordered (starting) row layout.
    *
    * starting_amalg_bdry is a vector containing the bc.
    *
    * old_x, old_y, and old_z are passed to the function,
    * so that they can be used with (or in substitution to) the graph.
    *
    * See also notes in the function itself. */
   starting_aggr_count =
     ML_DecomposeGraph_with_Zoltan(Amatrix, starting_aggr_count,
				   starting_decomposition,
				   starting_amalg_bdry,
				   old_x, old_y, old_z,
				   ml_ag->cur_level,
                                   grid_info->zoltan_type,
                                   grid_info->zoltan_estimated_its,
				   grid_info->zoltan_timers,
                                   grid_info->smoothing_steps,
                                   nullspace_dim);
   
   /* From now on the code should not change because of Zoltan
    * until the next Zoltan comment (marked with `RAY')... */

   if( starting_aggr_count <= 0 ) {
     fprintf( stderr,
	      "*ML*ERR* Something went *very* wrong in Zoltan...\n"
	      "*ML*ERR* (file %s, line %d)\n",
	      __FILE__,
	      __LINE__ );
     exit( EXIT_FAILURE );
   }
   
   if( mypid == 0 && 7 < ML_Get_PrintLevel() ) 
     printf("%s Using %d aggregates (globally)\n",
	    str,
	    starting_aggr_count );
   
   if( mypid == 0 && 7 < ML_Get_PrintLevel() ) {
     printf("%s # aggre/ # (block) rows = %7.3f %%  (= %d/%d)\n",
	    str,
	    100.0*starting_aggr_count/Nrows_global,
	    starting_aggr_count,
	    Nrows_global);
   }
   
   /* ********************************************************************** */
   /* compute operator complexity                                            */
   /* ********************************************************************** */
   
   Nnonzeros2 = ML_Comm_GsumDouble( comm, Nnonzeros2);

   if ( mypid == 0 && 7 < ML_Get_PrintLevel())
     printf("%s Total (block) nnz = %g ( = %5.2f/(block)row)\n",
	    str,
	    Nnonzeros2,1.0*Nnonzeros2/Nrows_global);
   
   if ( ml_ag->operator_complexity == 0.0 ) {
      ml_ag->fine_complexity = Nnonzeros2;
      ml_ag->operator_complexity = Nnonzeros2;
   }
   else ml_ag->operator_complexity += Nnonzeros2;

   /* FIXME: erase meeeeeeeee
      fix aggr_index for num_PDE_eqns > 1 
   
   for (i = Nrows - 1; i >= 0; i-- ) {
      for (j = num_PDE_eqns-1; j >= 0; j--) {
         aggr_index[i*num_PDE_eqns+j] = aggr_index[i];
      }
   }
   */
   
   /* ********************************************************************** */
   /* I allocate room to copy aggr_index and pass this value to the user,    */
   /* who will be able to analyze and visualize this after the construction  */
   /* of the levels. This way, the only price we have to pay for stats and   */
   /* viz is essentially a little bit of memory.                             */
   /* this memory will be cleaned with the object ML_Aggregate ml_ag.        */
   /* I set the pointers using the ML_Aggregate_Info structure. This is      */
   /* allocated using ML_Aggregate_Info_Setup(ml,MaxNumLevels)               */
   /* ********************************************************************** */
   
   if( Amatrix->to->Grid->Grid != NULL ) {

     graph_decomposition = (int *) ML_allocate(sizeof(int)*Nrows );

     if( graph_decomposition == NULL ) {
       fprintf( stderr,
		"*ML*ERR* Not enough memory for %d bytes\n"
		"*ML*ERR* (file %s, line %d)\n",
		(int)sizeof(int)*Nrows,
		__FILE__,
	      __LINE__ );
       exit( EXIT_FAILURE );
     }

     for( i=0 ; i<Nrows ; i++ ) {
       graph_decomposition[i] = starting_decomposition[i];
     }

     aggr_viz_and_stats = (ML_Aggregate_Viz_Stats *) (Amatrix->to->Grid->Grid);
     aggr_viz_and_stats->graph_decomposition = graph_decomposition;
     aggr_viz_and_stats->Nlocal = Nrows;
     aggr_viz_and_stats->Naggregates = starting_aggr_count;
     aggr_viz_and_stats->local_or_global = ML_GLOBAL_INDICES;
     aggr_viz_and_stats->is_filled = ML_YES;
     
   }

   /* ********************************************************************** */
   /* Compute the new distribution, so that `desired_aggre_per_proc' aggre   */
   /* are stored on each processor (up to the maximum number of aggregates). */
   /* - starting_offset : decomposition of the unknowns for the finer grid   */
   /*                     before redistribution                              */
   /* - reordered_offset : decomposition of the unknowns for the finer grid  */
   /*                      as will be after redistribution, as done by       */
   /*                      operator QQ                                       */
   /* - Nrows, new_Nrows : number of local rows for the finer grid before    */
   /*                      and after redistribution                          */
   /* ********************************************************************** */

   starting_offset  = (int *)ML_allocate( sizeof(int) * (Nprocs+1));
   reordered_offset = (int *)ML_allocate( sizeof(int) * (Nprocs+1));
   nodes_per_aggre = (int *) ML_allocate( sizeof(int) * starting_aggr_count );

   if( starting_offset == NULL || reordered_offset == NULL
       || nodes_per_aggre == NULL ) {
     fprintf( stderr,
	      "*ML*ERR* Not enough memory\n"
	      "*ML*ERR* (file %s, line %d)\n",
	      __FILE__,
	      __LINE__ );
     exit( EXIT_FAILURE );
   }
   
   ML_DecomposeGraph_BuildOffsets( Nrows, starting_offset, Nprocs,
				   Amatrix->comm->USR_comm);
   
   /* ********************************************************************** */
   /* Compute how many nodes are contained in each aggregate. This will be   */
   /* done for all the aggregates (so some communications will occur).       */
   /* ********************************************************************** */
   
   ML_CountNodesPerAggre( Nrows, starting_decomposition,
			  starting_aggr_count, nodes_per_aggre,
			  Amatrix->comm->USR_comm );

   /* ********************************************************************** */
   /* Compute how many aggregates will be stored on this process. This is    */
   /* based on the `desired_aggre_per_proc', so that the first processes will*/
   /* have about this number (and then maybe some processes will have none). */
   /* This is used to determine a reorderd offset, so that each processor    */
   /* will hold the rows of the matrix required to form the given aggregates */
   /* This new row decomposition is hold in `reordered_decomposition'        */
   /* ********************************************************************** */

   aggr_count = ML_BuildReorderedOffset(starting_offset,
					desired_aggre_per_proc,
					Nprocs, nodes_per_aggre,
					starting_aggr_count,
					reordered_offset, mypid );

   new_Nrows = reordered_offset[mypid+1] - reordered_offset[mypid];
   
   i = 0;
   if( new_Nrows > 0 ) i++;

#ifdef ML_MPI
   MPI_Reduce( &i, &j, 1, MPI_INT, MPI_SUM, 0, comm->USR_comm);
#else
   j = i;
#endif

   if( mypid == 0 && 5 < ML_Get_PrintLevel() ) {
     printf( "%s Processes with at least 1 row at next level = %d\n",
	     str,
	     j );
   } 
   
   reordered_decomposition = (int *) ML_allocate( sizeof(int) * (Nrows+1) );
   if( reordered_decomposition == NULL ) {
     fprintf( stderr,
	      "*ML*ERR* Not enough memory to allocate %d bytes\n"
	      "*ML*ERR* (file %s, line %d)\n",
	      (int)sizeof(int) * (Nrows+1),
	      __FILE__,
	      __LINE__ );
     exit( EXIT_FAILURE );
   }

   ML_BuildReorderedDecomposition(starting_decomposition,
				  reordered_decomposition, Nrows,
				  starting_aggr_count, nodes_per_aggre,
				  Amatrix->comm->USR_comm );
   
   /* ********************************************************************** */
   /* Finally, built the operator QQ, moving from the reorderd decomposition */
   /* to the starting one. So, QQ will be applied to vectors (or matrices)   */
   /* whose row decomposition is the reordered one. I need QQ because we have*/
   /* P_tent = QQ * \hat{P}_tent, where P_tent is the tentative prolongator  */
   /* as used by ML (from the next level to this level, in the starting      */
   /* decomposition), and \hat{P}_tent the one based on the reordered dec.   */
   /* ********************************************************************** */

   if( nullspace_vect == NULL /*&& diff_level == 0*/ ) {
     new_nullspace_vect = nullspace_vect;
     i = 0;
   } else {
     nbytes = sizeof(double) * (new_Nrows * num_PDE_eqns * nullspace_dim );

     if( nbytes == 0 ) new_nullspace_vect = NULL;
     else {
       new_nullspace_vect = (double *) ML_allocate( nbytes );
       if( new_nullspace_vect == NULL ) {
	 fprintf( stderr,
		  "*ML*ERR* Not enough memory to allocate %d bytes\n"
		  "*ML*ERR* (file %s, line %d)\n",
		  nbytes,
		  __FILE__,
		  __LINE__ );
	 exit( EXIT_FAILURE );
       }
     }
       
     i = 1;
   }

   reordered_amalg_bdry = (double *) ML_allocate(sizeof(double)*(new_Nrows+1));
   
#ifdef ML_WITH_EPETRA
   if( nullspace_dim != num_PDE_eqns ) {
     printf("---------> Never tested with nullspace_dim != num_PDE_eqns\n"
	    "---------> Memory allocation within  ML_BuildQ to be checked...\n" );
   }
   QQ = ML_BuildQ(Nrows, new_Nrows, num_PDE_eqns, nullspace_dim,
		  reordered_decomposition,
		  nullspace_vect, new_nullspace_vect, i,
		  starting_amalg_bdry, reordered_amalg_bdry,
		  Amatrix->comm->USR_comm,
		  comm );
#else
   if( mypid == 0 ) 
     fprintf( stderr,
	      "*ML*ERR* Sorry, you cannot redistribute matrices within the Zoltan\n"
	      "*ML*ERR* aggregation without epetra. Please recompile using epetra...\n" );
   exit( EXIT_FAILURE );
#endif

   if (starting_decomposition != NULL) {
     ML_free( starting_decomposition );
     starting_decomposition = NULL;
   }
   if (reordered_decomposition != NULL) {
     ML_free( reordered_decomposition );
     reordered_decomposition = NULL;
   }
   if (starting_amalg_bdry != NULL) {
     ML_free( starting_amalg_bdry );
     starting_amalg_bdry = NULL;
   }
   if (starting_offset != NULL) {
     ML_free( starting_offset );
     starting_offset = NULL;
   }
   if (reordered_offset != NULL) {
     ML_free( reordered_offset );
     reordered_offset = NULL;
   }

   /* ********************************************************************** */
   /* Now reallocating aggr_index so that we can build the prolongator       */
   /* as if all the aggregates were local. Need some reallocations here.     */
   /* ********************************************************************** */

   nbytes = sizeof(int) * new_Nrows * num_PDE_eqns;
   
   if ( nbytes > 0 ) ML_memory_alloc((void**) &aggr_index, nbytes, "ACJ");
   else              aggr_index = NULL;

   /* k is the ID of the first aggregate on this proc */
#ifdef ML_MPI
   MPI_Scan( &aggr_count, &k, 1, MPI_INT, MPI_SUM, Amatrix->comm->USR_comm );
   k -= aggr_count;
#else
   k = 0;
#endif

   if( new_Nrows != 0 ) {

     jj = 0;
     for( i=0 ; i<aggr_count ; i++ ) {
       for( j=0 ; j<nodes_per_aggre[i+k] ; j++ ) {
	 aggr_index[jj] = i;
	 jj++;
       }
     }

     if( new_Nrows != jj ) {
       fprintf( stderr,
		"*ML*ERR* something went wrong in coarsening with Zoltan:\n"
		"*ML*ERR* new_Nrows = %d, jj = %d\n"
		"*ML*ERR* (file %s, line %d)\n",
		new_Nrows, jj,
		__FILE__,
		__LINE__ );
       exit( EXIT_FAILURE );
     }
   }
   
   if (nodes_per_aggre != NULL) {
     ML_free( nodes_per_aggre );
     nodes_per_aggre = NULL;
   }

   /* *************************** */
   /* TO ADD: OPERATOR COMPLEXITY */
   /* *************************** */

   for (i = new_Nrows - 1; i >= 0; i-- ) {
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
   
   new_Nrows  *= num_PDE_eqns;
   Nrows_global  *= num_PDE_eqns;

   /* count the size of each aggregate. Now all aggregates are local */
   
   aggr_cnt_array = (int *) ML_allocate(sizeof(int)*aggr_count);
   for (i = 0; i < aggr_count ; i++) aggr_cnt_array[i] = 0;
   for (i = 0; i < new_Nrows; i++) 
      if (aggr_index[i] >= 0) 
         aggr_cnt_array[aggr_index[i]]++;

   /* ============================================================= */
   /* Form tentative prolongator                                    */
   /* ============================================================= */
   
   Ncoarse = aggr_count;
   
   /* ============================================================= */
   /* check and copy aggr_index                                     */
   /* ------------------------------------------------------------- */

   level = ml_ag->cur_level;
   nbytes = new_Nrows * sizeof( int );
   ML_memory_alloc((void**) &(ml_ag->aggr_info[level]), nbytes, "AGl");
   count = aggr_count;
   for ( i = 0; i < new_Nrows; i+=num_PDE_eqns ) 
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

   for ( i = 0; i < new_Nrows; i++ ) 
   {
      if ( aggr_index[i] >= Ncoarse ) 
      {
         printf("*ML*WRN* index out of bound %d = %d (%d)\n"
		"*ML*WRN* (file %s, line %d)\n",
		i, aggr_index[i], 
                Ncoarse,
		__FILE__,
		__LINE__ );
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
   if( nbytes != 0 ) {
     ML_memory_alloc((void**)&(new_null),nbytes,"AGr");
     for (i = 0; i < Ncoarse*nullspace_dim*nullspace_dim; i++) 
       new_null[i] = 0.0;
   } else {
     new_null = NULL;
   }
   
   /* ------------------------------------------------------------- */
   /* initialize the row pointer for the CSR prolongation operator  */
   /* (each row will have at most nullspace_dim nonzero entries)    */
   /* ------------------------------------------------------------- */

   for (i = 0; i <= new_Nrows; i++) new_ia[i] = i * nullspace_dim;

   /* trying this when a Dirichlet row is taken out */
   j = 0;
   new_ia[0] = 0;
   for (i = 0; i < new_Nrows; i++) {
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
         printf("*ML*ERR* couldn't allocate memory in CoarsenZoltan\n");
         exit(1);
      }
   }
   for (i = 0; i < new_Nrows; i+=num_PDE_eqns) 
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

   max_agg_size = 0;
   for (i = 0; i < aggr_count; i++) 
   {
      if (aggr_cnt_array[i] > max_agg_size) max_agg_size = aggr_cnt_array[i];
   }
   nbytes = max_agg_size * nullspace_dim * sizeof(double);
   if( nbytes > 0 ) ML_memory_alloc((void**)&qr_tmp, nbytes, "AGu");
   else             qr_tmp = NULL;
   nbytes = nullspace_dim * sizeof(double);
   if( nbytes > 0 ) ML_memory_alloc((void**)&tmp_vect, nbytes, "AGv");
   else             tmp_vect = NULL;
   
   lwork  = nullspace_dim;
   nbytes = nullspace_dim * sizeof(double);
   if( nbytes > 0 ) ML_memory_alloc((void**)&work, nbytes, "AGw");
   else             work = NULL;
  
   /* ------------------------------------------------------------- */
   /* perform block QR decomposition                                */
   /* ------------------------------------------------------------- */
     
   for (i = 0; i < aggr_count; i++) 
   {
      /* ---------------------------------------------------------- */
      /* set up the matrix we want to decompose into Q and R:       */
      /* ---------------------------------------------------------- */

      length = aggr_cnt_array[i];

      if (new_nullspace_vect == NULL) 
      {
         for (j = 0; j < (int) length; j++)
         {
            index = rows_in_aggs[i][j];
            for (k = 0; k < nullspace_dim; k++)
            {
              if ( reordered_amalg_bdry[index/num_PDE_eqns] == 1.0) qr_tmp[k*length+j] = 0.;
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
	       
               if ( reordered_amalg_bdry[index/num_PDE_eqns] == 1.0) qr_tmp[k*length+j] = 0.;
               else {
                  if (index < new_Nrows) {
		    qr_tmp[k*length+j] = new_nullspace_vect[k*new_Nrows+index];
                  }
                  else {
		    fprintf( stderr,
			     "*ML*ERR* error in QR factorization within Zoltan aggregation\n"
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
	  pr_error("ErrOr in CoarsenZoltan : "
		   "dgeqrf returned a non-zero %d %d\n",
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
	  printf("*ML*ERR* in dorgqr on %d row (dims are %d, %d)\n",
		 i,aggr_cnt_array[i],
                 nullspace_dim);
	  printf("*ML*ERR* performing QR on a MxN matrix where M<N.\n");
	}
	DORGQR_F77(&(aggr_cnt_array[i]), &nullspace_dim,
			  &nullspace_dim, qr_tmp, &(aggr_cnt_array[i]),
			  tmp_vect, work, &lwork, &info);
	if (info != 0) {
	  printf("*ML*ERR* in dorgqr on %d row (dims are %d, %d)\n",
		 i,aggr_cnt_array[i],
                 nullspace_dim);
	  pr_error("Error in CoarsenZoltan: dorgqr returned a non-zero\n");
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
	    if ( index < new_Nrows )
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
			 "*ML*ERR* error in QR factorization within Zoltan\n" );
		exit( EXIT_FAILURE );
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

   } /* for( i over aggregates ) */
   
   ML_Aggregate_Set_NullSpace(ml_ag, num_PDE_eqns, nullspace_dim, 
                              new_null, Ncoarse*nullspace_dim);
   if( new_null != NULL ) {
     ML_memory_free( (void **) &new_null);
     new_null = NULL;
   }
   
   /* ------------------------------------------------------------- */
   /* set up the csr_data data structure                            */
   /* ------------------------------------------------------------- */

   ML_memory_alloc((void**) &csr_data, sizeof(struct ML_CSR_MSRdata),"CSR");
   csr_data->rowptr  = new_ia;
   csr_data->columns = new_ja;
   csr_data->values  = new_val;

   Pstart = ML_Operator_Create( Amatrix->comm );
			       
   ML_Operator_Set_ApplyFuncData( Pstart, nullspace_dim*Ncoarse, new_Nrows, 
                                  csr_data, new_Nrows, NULL, 0);
   Pstart->data_destroy = ML_CSR_MSR_ML_memorydata_Destroy;

   Pstart->getrow->pre_comm = ML_CommInfoOP_Create();
   
   ML_Operator_Set_Getrow((Pstart), new_Nrows, CSR_getrow);
   ML_Operator_Set_ApplyFunc(Pstart, CSR_matvec);
   Pstart->max_nz_per_row = 1;

   Pmatrix2 = ML_Operator_Create( Amatrix->comm );
   
   ML_2matmult(QQ, Pstart, Pmatrix2, ML_CSR_MATRIX );
   

   ML_Operator_Set_1Levels(Pmatrix2, (*Pmatrix)->from, (*Pmatrix)->to);
   ML_Operator_Set_BdryPts(Pmatrix2, (*Pmatrix)->BCs);
   /* JJH I've observed that (*Pmatrix)->label is sometimes null.
      JJH Not sure if this is a problem. */
   if ((*Pmatrix)->label) ML_Operator_Set_Label(Pmatrix2,(*Pmatrix)->label);
   else                   ML_Operator_Set_Label(Pmatrix2,"unknown");
/* this must be set so that the hierarchy generation does not abort early
   in adaptive SA */
   Pmatrix2->num_PDEs = nullspace_dim;   


   ML_Operator_Clean( *Pmatrix );

   memcpy((void *) *Pmatrix, (void *)Pmatrix2, sizeof(ML_Operator));
   /* FIXME : am I ok  ????? */
   ML_free(Pmatrix2);
      
   /* ********************************************************************** */
   /* I have to destroy the tentative local matrix, and the redistribution   */
   /* matrix QQ. This is actually an ML_Operator on the top of an Epetra     */
   /* object. So, I call ML_DestroyQ, which is a CPP function, to delete the */
   /* memory interally used by Epetra.                                       */
   /* ********************************************************************** */

   ML_Operator_Destroy( &Pstart ); 
   ML_Operator_Destroy( &QQ );

#ifdef ML_WITH_EPETRA
   ML_DestroyQ( );
#endif
   
   /* ------------------------------------------------------------- */
   /* clean up                                                      */
   /* ------------------------------------------------------------- */

   ML_memory_free((void**) &aggr_index);
   if( aggr_cnt_array != NULL ) ML_free(aggr_cnt_array);
   for (i = 0; i < aggr_count; i++) ML_free(rows_in_aggs[i]);
   ML_memory_free((void**)&rows_in_aggs);
   if( qr_tmp != NULL ) {
     ML_memory_free((void**)&qr_tmp);
     qr_tmp = NULL;
   }
   if( tmp_vect != NULL ) {
     ML_memory_free((void**)&tmp_vect);
     tmp_vect = NULL;
   }
   if( work != NULL ) {
     ML_memory_free((void**)&work);
     work = NULL;
   }

   aggr_curr = aggr_head;
   while ( aggr_curr != NULL ) 
   {
      supernode = aggr_curr;
      aggr_curr = aggr_curr->next;
      if ( supernode->length > 0 ) ML_free( supernode->list );
      ML_free( supernode );
   }

   if( reordered_amalg_bdry != NULL ) {
     ML_free( reordered_amalg_bdry );
     reordered_amalg_bdry = NULL;
   }
   if( new_nullspace_vect != NULL ) {
     ML_free( new_nullspace_vect );
     new_nullspace_vect = NULL;
   }

   /* ------------------- that's all folks --------------------------------- */

   return Ncoarse*nullspace_dim;

} /* ML_Aggregate_CoarsenZoltan */
