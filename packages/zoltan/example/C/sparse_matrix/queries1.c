/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     * 
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    Revision: 1.79.2.8 $
 ****************************************************************************/

#include <stdio.h>
#include "zoltan.h"
#include "matrix.h"


/* Simple example of how to use Zoltan.
 * This file contains the query functions
 * for a sparse matrix.
 */

/*
 *  PROTOTYPES for load-balancer interface functions.
 */

ZOLTAN_NUM_OBJ_FN get_num_entries;
ZOLTAN_OBJ_LIST_FN get_entries;

ZOLTAN_NUM_GEOM_FN get_num_geom;
ZOLTAN_GEOM_MULTI_FN get_geom_multi;

#ifdef USE_GRAPH
/* These two functions are needed to use ParMETIS thru Zoltan. */
/* They are not needed for geometric methods. */
ZOLTAN_NUM_EDGES_MULTI_FN get_num_edges_multi;
ZOLTAN_EDGE_LIST_MULTI_FN get_edge_list_multi;
#endif


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int setup_zoltan(struct Zoltan_Struct *zz, Matrix A)
{
  /* Set the load-balance method */
  /* You can change "RCB" to any Zoltan method below. */

  if (Zoltan_Set_Param(zz, "LB_METHOD", "RCB") == ZOLTAN_FATAL) {
    printf("fatal:  error returned from Zoltan_Set_Param(LB_METHOD)\n");
    return 0;
  }

  /*
   * Set the callback functions
   */

  if (Zoltan_Set_Fn(zz, ZOLTAN_NUM_OBJ_FN_TYPE, (void (*)()) get_num_entries,
                    (void *) &A) == ZOLTAN_FATAL) {
    printf("fatal:  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }

  if (Zoltan_Set_Fn(zz, ZOLTAN_OBJ_LIST_FN_TYPE,
                    (void (*)()) get_entries,
                    (void *) &A) == ZOLTAN_FATAL) {
    printf("fatal:  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }

  /* Functions for geometry based algorithms */
  if (Zoltan_Set_Fn(zz, ZOLTAN_NUM_GEOM_FN_TYPE, (void (*)()) get_num_geom,
                    (void *) &A) == ZOLTAN_FATAL) {
    printf("fatal:  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }

  if (Zoltan_Set_Fn(zz, ZOLTAN_GEOM_MULTI_FN_TYPE,
                    (void (*)()) get_geom_multi,
                    (void *) &A) == ZOLTAN_FATAL) {
    printf("fatal:  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }

#ifdef USE_GRAPH
  /* Functions for graph based algorithms */
  /* These two functions are needed to use ParMETIS thru Zoltan. */
  /* They are not needed for geometric methods. */
  if (Zoltan_Set_Fn(zz, ZOLTAN_NUM_EDGES_MULTI_FN_TYPE,
                    (void (*)()) get_num_edges_multi,
                    (void *) &A) == ZOLTAN_FATAL) {
    printf("fatal:  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }
  if (Zoltan_Set_Fn(zz, ZOLTAN_EDGE_LIST_MULTI_FN_TYPE,
                    (void (*)()) get_edge_list_multi,
                    (void *) &A) == ZOLTAN_FATAL) {
    printf("fatal:  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }
#endif /* USE_GRAPH */

  return 1;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int run_zoltan(struct Zoltan_Struct *zz, Matrix A)
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

  int *ptr;
  int myrank;
  int k;

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
    printf("fatal:  error returned from Zoltan_LB_Partition()\n");
    return 0;
  }

  /* Print results. */
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (new_decomp){
    printf("[Proc %1d] My data to export are:\n", myrank);
    for (k=0; k<num_exported; k++){
      ptr = (int *) &export_gids[num_gid_entries*k];
      printf("[Proc %1d] Export (%d,%d) to proc %d\n", myrank, 
         *ptr, *(ptr+1), export_procs[k]);
    }
  }
  else
    printf("[Proc %1d] No changes over here.\n", myrank);

  /* Data migration goes here. */

  /* Clean up */
  Zoltan_LB_Free_Part(&import_gids, &import_lids,
                      &import_procs, &import_to_part);
  Zoltan_LB_Free_Part(&export_gids, &export_lids,
                      &export_procs, &export_to_part);

  return 1;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Each Zoltan "object" is a nonzero in the matrix */
int get_num_entries (void *data, int *ierr)
{

  Matrix *A;

  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return 0;
  }
  /* cast the data pointer to correct data type */
  A = (Matrix *) data;

  *ierr = ZOLTAN_OK; /* set error code */

  return(A->num_mynz);  /* # local nonzeroes */
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void get_entries(void *data, int num_gid_entries, int num_lid_entries,
                  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
                  int wdim, float *wgt, int *ierr)
{
/* Return global ids for each matrix entry. Must be unique.
 * We use (i,j) as the gid for a nonzero in position (i,j).
 * Most applications use integer gids but Zoltan allows several ints.
 */
  Matrix *A;
  int k;

  *ierr = ZOLTAN_OK; 

  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return;
  }

  /* We should be using (at least) two ints for each GID. */
  if (num_gid_entries < 2) {
    *ierr = ZOLTAN_FATAL;
    return;
  }
  
  A = (Matrix *) data;
  for (k = 0; k < A->num_mynz; k++) {
    global_id[k*num_gid_entries] =   (ZOLTAN_ID_TYPE) A->entries[k].i;
    global_id[k*num_gid_entries+1] = (ZOLTAN_ID_TYPE) A->entries[k].j;

    /* Add (optional) local ids and/or weights here if desired.  */
  }

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int get_num_geom(void *data, int *ierr)
{
  /* Our coordinates are the (i,j) indices, so return 2. */

  *ierr = ZOLTAN_OK; /* set error flag */
  return(2);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void get_geom_multi(void *data, int num_gid_entries, int num_lid_entries,
              int num_obj, ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
              int num_dim, double *coor, int *ierr)
{
/*  Return in coor the coordinates for each nonzero
 *  in the global_id/local_id arrays.  This example uses
 *  only global ids. (Local ids are optional.)
 */
Matrix *A; 
int k;

  *ierr = ZOLTAN_OK;

  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return;
  }
  A = (Matrix *) data;

  for (k = 0; k < A->num_mynz; k++) {
    /*
     * The coordinates of matrix entry (i,j)
     * are simply (i,j).
     */
      coor[2*k]   = (double) A->entries[k].i;
      coor[2*k+1] = (double) A->entries[k].j;
  }
}

#ifdef GRAPH_QUERY
/* Graph query functions not yet implemented. */
#endif

