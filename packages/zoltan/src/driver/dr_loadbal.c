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
 *    $Revision$
 ****************************************************************************/


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdlib.h>
#include <stdio.h>
#include <strings.h>

#include <mpi.h>

#include "dr_const.h"
#include "dr_err_const.h"
#include "dr_loadbal_const.h"
#include "dr_eval_const.h"

static int Num_GID = 1, Num_LID = 1;

/*--------------------------------------------------------------------------*/
/* Purpose: Call Zoltan to determine a new load balance.                    */
/*          Contains all of the callback functions that Zoltan needs        */
/*          for the load balancing.                                         */
/*--------------------------------------------------------------------------*/

/*
 *  PROTOTYPES for load-balancer interface functions.
 */

ZOLTAN_NUM_OBJ_FN get_num_elements;
/* not used right now --->
ZOLTAN_OBJ_LIST_FN get_elements;
*/
ZOLTAN_FIRST_OBJ_FN get_first_element;
ZOLTAN_NEXT_OBJ_FN get_next_element;

ZOLTAN_NUM_GEOM_FN get_num_geom;
ZOLTAN_GEOM_FN get_geom;

ZOLTAN_NUM_EDGES_FN get_num_edges;
ZOLTAN_EDGE_LIST_FN get_edge_list;

ZOLTAN_NUM_CHILD_FN get_num_child;
ZOLTAN_CHILD_LIST_FN get_child_elements;
ZOLTAN_FIRST_COARSE_OBJ_FN get_first_coarse_element;
ZOLTAN_NEXT_COARSE_OBJ_FN get_next_coarse_element;

ZOLTAN_PARTITION_FN get_partition;

ZOLTAN_HG_EDGE_LIST_FN get_hg_edge_list;
ZOLTAN_NUM_HG_EDGES_FN get_num_hg_edges;
ZOLTAN_NUM_HG_PINS_FN get_num_hg_pins;

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int setup_zoltan(struct Zoltan_Struct *zz, int Proc, PROB_INFO_PTR prob,
                 MESH_INFO_PTR mesh)
{
/* Local declarations. */
  char *yo = "run_zoltan";

#define MAXPROC 20
  int i;                         /* Loop index */
  float psize[MAXPROC];          /* Partition size */
  int part[MAXPROC];             /* Partition numbers */
  int idx[MAXPROC];              /* Index numbers */
  int ierr;                      /* Error code */
  char errmsg[128];              /* Error message */

  DEBUG_TRACE_START(Proc, yo);

  /* Set the user-specified parameters */
  for (i = 0; i < prob->num_params; i++) {
    ierr = Zoltan_Set_Param(zz, prob->params[i][0], prob->params[i][1]);
    if (ierr == ZOLTAN_FATAL) {
      sprintf(errmsg,"fatal: error in Zoltan_Set_Param when setting parameter %s\n",
              prob->params[i][0]);
      Gen_Error(0, errmsg);
      return 0;
    }
    if (strcasecmp(prob->params[i][0], "NUM_GID_ENTRIES") == 0) 
      Num_GID = atoi(prob->params[i][1]);
    else if (strcasecmp(prob->params[i][0], "NUM_LID_ENTRIES") == 0) 
      Num_LID = atoi(prob->params[i][1]);
  }

  /* Set the method */
  if (Zoltan_Set_Param(zz, "LB_METHOD", prob->method) == ZOLTAN_FATAL) {
    Gen_Error(0, "fatal:  error returned from Zoltan_Set_Param(LB_METHOD)\n");
    return 0;
  }

  if (Test_Local_Partitions == 1) {
    /* Compute Proc partitions for each processor */
    char s[8];
    sprintf(s, "%d", Proc);
    if (Zoltan_Set_Param(zz, "NUM_LOCAL_PARTITIONS", s) == ZOLTAN_FATAL) {
      Gen_Error(0, "fatal:  error returned from Zoltan_Set_Param()\n");
      return 0;
    }
  }
  else if (Test_Local_Partitions == 2) {
    /* Compute Proc partitions for odd-ranked processors; let remaining
     * partitions be in even-ranked processors. */
    if (Proc%2) {
      char s[8];
      sprintf(s, "%d", Proc);
      if (Zoltan_Set_Param(zz, "NUM_LOCAL_PARTITIONS", s) == ZOLTAN_FATAL) {
        Gen_Error(0, "fatal:  error returned from Zoltan_Set_Param()\n");
        return 0;
      }
    }
  }
  else if (Test_Local_Partitions == 3) {
    /* Variable partition sizes, but one partition per proc */
    i = 0;
    psize[0] = (float) Proc;              /* Partition size = proc number */
    Zoltan_LB_Set_Part_Sizes(zz, 0, 1, &Proc, &i, psize);
  }
  else if (Test_Local_Partitions == 4) {
    /* Variable number of partitions per proc and partition sizes. */
    /* Compute Proc partitions for each processor */
    char s[8];
    sprintf(s, "%d", Proc);
    if (Zoltan_Set_Param(zz, "NUM_LOCAL_PARTITIONS", s) == ZOLTAN_FATAL) {
      Gen_Error(0, "fatal:  error returned from Zoltan_Set_Param()\n");
      return 0;
    }
    /* Each partition size is inverse to the no. of partitions on a proc. */ 
    for (i=0; i<Proc; i++){
      part[i] = Proc*(Proc-1)/2 + i;  /* Global part number */
      idx[i] = 0;
      psize[i] = 1.0/Proc; 
    }
    if (Proc) Zoltan_LB_Set_Part_Sizes(zz, 0, Proc, part, idx, psize);
  }


  /*
   * Set the callback functions
   */

  if (Zoltan_Set_Fn(zz, ZOLTAN_NUM_OBJ_FN_TYPE, (void (*)()) get_num_elements,
                    (void *) mesh) == ZOLTAN_FATAL) {
    Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }

  if (Zoltan_Set_Fn(zz, ZOLTAN_FIRST_OBJ_FN_TYPE,
                    (void (*)()) get_first_element,
                    (void *) mesh) == ZOLTAN_FATAL) {
    Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }

  if (Zoltan_Set_Fn(zz, ZOLTAN_NEXT_OBJ_FN_TYPE, (void (*)()) get_next_element,
                    (void *) mesh) == ZOLTAN_FATAL) {
    Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }

  /* Functions for geometry based algorithms */
  if (Zoltan_Set_Fn(zz, ZOLTAN_NUM_GEOM_FN_TYPE, (void (*)()) get_num_geom,
                    (void *) mesh) == ZOLTAN_FATAL) {
    Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }

  if (Zoltan_Set_Fn(zz, ZOLTAN_GEOM_FN_TYPE, (void (*)()) get_geom,
                    (void *) mesh) == ZOLTAN_FATAL) {
    Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }

  /* Functions for graph based algorithms */
  if (Zoltan_Set_Fn(zz, ZOLTAN_NUM_EDGES_FN_TYPE, (void (*)()) get_num_edges,
                    (void *) mesh) == ZOLTAN_FATAL) {
    Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }

  if (Zoltan_Set_Fn(zz, ZOLTAN_EDGE_LIST_FN_TYPE, (void (*)()) get_edge_list,
                    (void *) mesh)== ZOLTAN_FATAL) {
    Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }

  /* Functions for tree-based algorithms */
  if (Zoltan_Set_Fn(zz, ZOLTAN_NUM_COARSE_OBJ_FN_TYPE,
                    (void (*)()) get_num_elements,
                    (void *) mesh) == ZOLTAN_FATAL) {
    Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }

  if (Zoltan_Set_Fn(zz, ZOLTAN_FIRST_COARSE_OBJ_FN_TYPE,
                    (void (*)()) get_first_coarse_element,
                    (void *) mesh) == ZOLTAN_FATAL) {
    Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }

  if (Zoltan_Set_Fn(zz, ZOLTAN_NEXT_COARSE_OBJ_FN_TYPE,
                    (void (*)()) get_next_coarse_element,
                    (void *) mesh) == ZOLTAN_FATAL) {
    Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }

  if (Zoltan_Set_Fn(zz, ZOLTAN_NUM_CHILD_FN_TYPE,
                    (void (*)()) get_num_child,
                    (void *) mesh) == ZOLTAN_FATAL) {
    Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }

  if (Zoltan_Set_Fn(zz, ZOLTAN_CHILD_LIST_FN_TYPE,
                    (void (*)()) get_child_elements,
                    (void *) mesh) == ZOLTAN_FATAL) {
    Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }

  if (mesh->data_type == HYPERGRAPH) {
    if (Zoltan_Set_Fn(zz, ZOLTAN_NUM_HG_EDGES_FN_TYPE, 
                      (void (*)()) get_num_hg_edges,
                      (void *) mesh) == ZOLTAN_FATAL) {
      Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
      return 0;
    }

    if (Zoltan_Set_Fn(zz, ZOLTAN_NUM_HG_PINS_FN_TYPE, 
                      (void (*)()) get_num_hg_pins,
                      (void *) mesh) == ZOLTAN_FATAL) {
      Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
      return 0;
    }

    if (Zoltan_Set_Fn(zz, ZOLTAN_HG_EDGE_LIST_FN_TYPE, 
                      (void (*)()) get_hg_edge_list,
                      (void *) mesh) == ZOLTAN_FATAL) {
      Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
      return 0;
    }
  }

  /* Functions for partitions */
  if (Zoltan_Set_Fn(zz, ZOLTAN_PARTITION_FN_TYPE, (void (*)()) get_partition,
                (void *) mesh) == ZOLTAN_FATAL) {
    Gen_Error(0, "fatal:  error returned from Zoltan_Set_Fn()\n");
    return 0;
  }

  DEBUG_TRACE_END(Proc, yo);
  return 1;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int run_zoltan(struct Zoltan_Struct *zz, int Proc, PROB_INFO_PTR prob,
               MESH_INFO_PTR mesh, PARIO_INFO_PTR pio_info)
{
/* Local declarations. */
  char *yo = "run_zoltan";

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

  int i;                         /* Loop index                               */
  int *order;			 /* Ordering vector(s) */
  ZOLTAN_ID_PTR order_gids = NULL;  /* List of all gids for ordering */
  ZOLTAN_ID_PTR order_lids = NULL;  /* List of all lids for ordering */
  double stime = 0.0, mytime = 0.0, maxtime = 0.0;

/***************************** BEGIN EXECUTION ******************************/

  DEBUG_TRACE_START(Proc, yo);

  if (Driver_Action & 1){

    /* Load balancing part */
  
    /* Evaluate the old balance */
    if (Debug_Driver > 0) {
      if (Proc == 0) printf("\nBEFORE load balancing\n");
      driver_eval(mesh);
      i = Zoltan_LB_Eval(zz, 1, NULL, NULL, NULL, NULL, NULL, NULL);
      if (i) printf("Warning: Zoltan_LB_Eval returned code %d\n", i);
    }
  
    /*
     * Call Zoltan
     */
    MPI_Barrier(MPI_COMM_WORLD);   /* For timings only */
    stime = MPI_Wtime();
    if (Zoltan_LB_Partition(zz, &new_decomp, &num_gid_entries, &num_lid_entries,
                 &num_imported, &import_gids,
                 &import_lids, &import_procs, &import_to_part,
                 &num_exported, &export_gids,
                 &export_lids, &export_procs, &export_to_part) == ZOLTAN_FATAL){
      Gen_Error(0, "fatal:  error returned from Zoltan_LB_Partition()\n");
      return 0;
    }
    mytime = MPI_Wtime() - stime;
    MPI_Allreduce(&mytime, &maxtime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (Proc == 0)
      printf("DRIVER:  Zoltan_LB_Partition time = %g\n", maxtime);

    {int mine[2], gmax[2], gmin[2];
    mine[0] = num_imported;
    mine[1] = num_exported;
    MPI_Allreduce(mine, gmax, 2, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(mine, gmin, 2, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (Proc == 0) {
      printf("DRIVER:  Min/Max Import: %d %d\n", gmin[0], gmax[0]);
      printf("DRIVER:  Min/Max Export: %d %d\n", gmin[1], gmax[1]);
    }
    }

#ifdef ZOLTAN_NEMESIS
    if (pio_info->file_type == NEMESIS_FILE && Nemesis_Output) {
      i = write_elem_vars(Proc, mesh, pio_info, num_exported, export_gids,
                          export_procs, export_to_part);
    }
#endif
  
    if (Test_Drops) {
      double x[] = {0.0L, 0.0L, 0.0L} ;
      int proc ;
      int status ;
      double xlo, ylo, zlo ;
      double xhi, yhi, zhi ;
      int procs[1000] ;
      int count ;
  
      status = Zoltan_LB_Point_Assign (zz, x, &proc) ;
      if (status != ZOLTAN_OK) printf ("Point_Assign returned an error\n") ;

      xlo = ylo = zlo = 0.0L ;
      xhi = yhi = zhi = 1.0L ;
      status = Zoltan_LB_Box_Assign (zz, xlo, ylo, zlo, xhi, yhi, zhi, procs, 
                                     &count) ;
      if (status != ZOLTAN_OK) printf ("Box_Assign returned an error\n") ;
    }
  
    /*
     * Call another routine to perform the migration
     */
    MPI_Barrier(MPI_COMM_WORLD);   /* For timings only */
    stime = MPI_Wtime();
    if (new_decomp) {
      if (!migrate_elements(Proc, mesh, zz, num_gid_entries, num_lid_entries,
          num_imported, import_gids, import_lids, import_procs, import_to_part,
          num_exported, export_gids, export_lids, export_procs, export_to_part))
      {
        Gen_Error(0, "fatal:  error returned from migrate_elements()\n");
        return 0;
      }
    }
    mytime = MPI_Wtime() - stime;
    MPI_Allreduce(&mytime, &maxtime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (Proc == 0)
      printf("DRIVER:  Total migration time = %g\n", maxtime);
  
    /* Evaluate the new balance */
    if (Debug_Driver > 0) {
      if (Proc == 0) printf("\nAFTER load balancing\n");
      driver_eval(mesh);
      i = Zoltan_LB_Eval(zz, 1, NULL, NULL, NULL, NULL, NULL, NULL);
      if (i) printf("Warning: Zoltan_LB_Eval returned code %d\n", i);
    }
  
    /* Clean up */
    Zoltan_LB_Free_Part(&import_gids, &import_lids,
                        &import_procs, &import_to_part);
    Zoltan_LB_Free_Part(&export_gids, &export_lids,
                        &export_procs, &export_to_part);
  }

  if (Driver_Action & 2){
    /* Only do ordering if this was specified in the driver input file */
    order = (int *) malloc (2*(mesh->num_elems) * sizeof(int));
    order_gids = (ZOLTAN_ID_PTR) malloc(mesh->num_elems * sizeof(int));
    order_lids = (ZOLTAN_ID_PTR) malloc(mesh->num_elems * sizeof(int));

    /* Evaluate the old ordering */
    if (Debug_Driver > 0) {
      if (Proc == 0) printf("\nBEFORE ordering\n");
      /* Not yet impl. */
    }

    if (Zoltan_Order(zz, &num_gid_entries, &num_lid_entries,
        mesh->num_elems, order_gids, order_lids,
        order, &order[mesh->num_elems], NULL) == ZOLTAN_FATAL) {
      Gen_Error(0, "fatal:  error returned from Zoltan_Order()\n");
      return 0;
    }

    /* Evaluate the new ordering */
    if (Debug_Driver > 0) {
      if (Proc == 0) printf("\nAFTER ordering\n");
      /* Not yet impl. */
    }

    /* Copy ordering permutation into mesh structure */
    for (i = 0; i < mesh->num_elems; i++){
      mesh->elements[i].perm_value = order[i];
      mesh->elements[i].invperm_value = order[(mesh->num_elems)+i];
    }

    /* Free order data */
    free(order);
    free(order_gids);
    free(order_lids);
  }

  DEBUG_TRACE_END(Proc, yo);
  return 1;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int get_num_elements(void *data, int *ierr)
{
MESH_INFO_PTR mesh;

  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return 0;
  }
  mesh = (MESH_INFO_PTR) data;

  *ierr = ZOLTAN_OK; /* set error code */

  return(mesh->num_elems);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int get_first_element(void *data, int num_gid_entries, int num_lid_entries,
                      ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
                      int wdim, float *wgt, int *ierr)
{
  MESH_INFO_PTR mesh;
  ELEM_INFO *elem;
  ELEM_INFO *current_elem;
  int i;
  int gid = num_gid_entries-1;
  int lid = num_lid_entries-1;

 *ierr = ZOLTAN_OK; 

  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return 0;
  }
  
  mesh = (MESH_INFO_PTR) data;
  if (mesh->num_elems == 0) {
    /* No elements on this processor */
    return 0;
  }

  elem = mesh->elements;
  current_elem = &elem[0];
  if (num_lid_entries) local_id[lid] = 0;
  global_id[gid] = (ZOLTAN_ID_TYPE) current_elem->globalID;

  if (wdim>0){
    for (i=0; i<wdim; i++){
      *wgt++ = current_elem->cpu_wgt[i];
      /* printf("Debug: In query function, object = %d, weight no. %1d = %f\n",
             global_id[gid], i, current_elem->cpu_wgt[i]); */
    }
  }

  return 1;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int get_next_element(void *data, int num_gid_entries, int num_lid_entries,
                     ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
                     ZOLTAN_ID_PTR next_global_id, ZOLTAN_ID_PTR next_local_id, 
                     int wdim, float *next_wgt, int *ierr)
{
  int found = 0;
  ELEM_INFO *elem;
  ELEM_INFO *next_elem;
  MESH_INFO_PTR mesh;
  int i, idx;
  int gid = num_gid_entries-1;
  int lid = num_lid_entries-1;

  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return 0;
  }
  
  mesh = (MESH_INFO_PTR) data;
  elem = mesh->elements;

  if (num_lid_entries) {
    idx = local_id[lid];
  }
  else {
    /* testing zero-length local IDs; search by global ID for current elem */
    (void) search_by_global_id(mesh, global_id[gid], &idx);
  }

  if (idx+1 < mesh->num_elems) { 
    found = 1;
    if (num_lid_entries) next_local_id[lid] = idx + 1;
    next_elem = &elem[idx+1];
    next_global_id[gid] = next_elem->globalID;

    if (wdim>0){
      for (i=0; i<wdim; i++){
        *next_wgt++ = next_elem->cpu_wgt[i];
      }
    }

    *ierr = ZOLTAN_OK; 
  }

  return(found);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int get_num_geom(void *data, int *ierr)
{
  MESH_INFO_PTR mesh;

  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return 0;
  }
  mesh = (MESH_INFO_PTR) data;

  *ierr = ZOLTAN_OK; /* set error flag */

  return(mesh->num_dims);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void get_geom(void *data, int num_gid_entries, int num_lid_entries,
              ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
              double *coor, int *ierr)
{
  ELEM_INFO *elem;
  ELEM_INFO *current_elem;
  int i, j, idx;
  double tmp;
  MESH_INFO_PTR mesh;
  int gid = num_gid_entries-1;
  int lid = num_lid_entries-1;

  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return;
  }
  mesh = (MESH_INFO_PTR) data;
  elem = mesh->elements;
  current_elem = (num_lid_entries 
                    ? &elem[local_id[lid]] 
                    : search_by_global_id(mesh, global_id[gid], &idx));

  if (mesh->eb_nnodes[current_elem->elem_blk] == 0) {
    /* No geometry info was read. */
    *ierr = ZOLTAN_FATAL;
    return;
  }
  
  /*
   * calculate the geometry of the element by averaging
   * the coordinates of the nodes in its connect table
   */
  for (i = 0; i < mesh->num_dims; i++) {
    tmp = 0.0;
    for (j = 0; j < mesh->eb_nnodes[current_elem->elem_blk]; j++)
      tmp += current_elem->coord[j][i];

    coor[i] = tmp / mesh->eb_nnodes[current_elem->elem_blk];
  }

  *ierr = ZOLTAN_OK;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int get_num_edges(void *data, int num_gid_entries, int num_lid_entries,
                  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *ierr)
{
  MESH_INFO_PTR mesh;
  ELEM_INFO *elem, *current_elem;
  int gid = num_gid_entries-1;
  int lid = num_lid_entries-1;
  int idx;

  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return 0;
  }
  mesh = (MESH_INFO_PTR) data;
  elem = mesh->elements;

  *ierr = ZOLTAN_OK;

  current_elem = (num_lid_entries 
                    ? &elem[local_id[lid]] 
                    : search_by_global_id(mesh, global_id[gid], &idx));

  return(current_elem->nadj);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void get_edge_list (void *data, int num_gid_entries, int num_lid_entries, 
                   ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
                   ZOLTAN_ID_PTR nbor_global_id, int *nbor_procs,
                   int get_ewgts, float *nbor_ewgts, int *ierr)
{
  MESH_INFO_PTR mesh;
  ELEM_INFO *elem;
  ELEM_INFO *current_elem;
  int i, j, proc, local_elem, idx;
  int gid = num_gid_entries-1;
  int lid = num_lid_entries-1;

  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return;
  }

  mesh = (MESH_INFO_PTR) data;
  elem = mesh->elements;
  current_elem = (num_lid_entries
                   ? &elem[local_id[lid]] 
                   : search_by_global_id(mesh, global_id[gid], &idx));

  /* get the processor number */
  MPI_Comm_rank(MPI_COMM_WORLD, &proc);

  j = 0;
  for (i = 0; i < current_elem->adj_len; i++) {

    /* Skip NULL adjacencies (sides that are not adjacent to another elem). */
    if (current_elem->adj[i] == -1) continue;

    if (current_elem->adj_proc[i] == proc) {
      local_elem = current_elem->adj[i];
      nbor_global_id[gid+j*num_gid_entries] = elem[local_elem].globalID;
    }
    else { /* adjacent element on another processor */
      nbor_global_id[gid+j*num_gid_entries] = current_elem->adj[i];
    }
    nbor_procs[j] = current_elem->adj_proc[i];

    if (get_ewgts) {
      if (current_elem->edge_wgt == NULL)
        nbor_ewgts[j] = 1.0; /* uniform weights is default */
      else
        nbor_ewgts[j] = current_elem->edge_wgt[i];
    }
    j++;
  }

  *ierr = ZOLTAN_OK;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int get_first_coarse_element(void *data, int num_gid_entries, 
                      int num_lid_entries,
                      ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
                      int *assigned, int *num_vert, ZOLTAN_ID_PTR vertices,
                      int *in_order, ZOLTAN_ID_PTR in_vertex, 
                      ZOLTAN_ID_PTR out_vertex, int *ierr)
{

MESH_INFO_PTR mesh;
ELEM_INFO *elem;
ELEM_INFO *current_elem;
int gid = num_gid_entries-1;
int lid = num_lid_entries-1;
int idx, i;
int ok;

  *ierr = ZOLTAN_OK;

  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return 0;
  }


  /* 
   * Assumption:  data is same for get_first_coarse_element and 
   * get_first_element 
   */
  ok = get_first_element(data, num_gid_entries, num_lid_entries, 
                         global_id, local_id, 0, NULL, ierr);

  if (ok) {
    mesh = (MESH_INFO_PTR) data;
    elem = mesh->elements;
    current_elem = (num_lid_entries
                     ? &elem[local_id[lid]]
                     : search_by_global_id(mesh, global_id[gid], &idx));

    *assigned = 1;
    *in_order = 0;
    *num_vert = mesh->eb_nnodes[current_elem->elem_blk];
    for (i = 0; i < *num_vert; i++)
      vertices[i*num_gid_entries + gid] = current_elem->connect[i];
  }
  return ok;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int get_next_coarse_element(void *data, int num_gid_entries, 
                      int num_lid_entries,
                      ZOLTAN_ID_PTR prev_global_id, ZOLTAN_ID_PTR prev_local_id,
                      ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
                      int *assigned, int *num_vert, ZOLTAN_ID_PTR vertices,
                      ZOLTAN_ID_PTR in_vertex, ZOLTAN_ID_PTR out_vertex,
                      int *ierr)
{

MESH_INFO_PTR mesh;
ELEM_INFO *elem;
ELEM_INFO *current_elem;
int gid = num_gid_entries-1;
int lid = num_lid_entries-1;
int idx, i;
int ok;

  *ierr = ZOLTAN_OK;

  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return 0;
  }


  /* 
   * Assumption:  data is same for get_first_coarse_element and 
   * get_first_element 
   */
  ok = get_next_element(data, num_gid_entries, num_lid_entries, 
                        prev_global_id, prev_local_id,
                        global_id, local_id,
                        0, NULL, ierr);

  if (ok) {
    mesh = (MESH_INFO_PTR) data;
    elem = mesh->elements;
    current_elem = (num_lid_entries
                     ? &elem[local_id[lid]]
                     : search_by_global_id(mesh, global_id[gid], &idx));

    *assigned = 1;
    *num_vert = mesh->eb_nnodes[current_elem->elem_blk];
    for (i = 0; i < *num_vert; i++)
      vertices[i*num_gid_entries + gid] = current_elem->connect[i];
  }
  return ok;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int get_num_child(void *data, int num_gid_entries, int num_lid_entries, 
                  ZOLTAN_ID_PTR global_id,
                  ZOLTAN_ID_PTR local_id, int *ierr)
{
  *ierr = ZOLTAN_OK;
  return 0;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void get_child_elements(void *data, int num_gid_entries, int num_lid_entries,
                   ZOLTAN_ID_PTR parent_gid, ZOLTAN_ID_PTR parent_lid, 
                   ZOLTAN_ID_PTR child_gids, ZOLTAN_ID_PTR child_lids, 
                   int *assigned, int *num_vert, ZOLTAN_ID_PTR vertices, 
                   ZOLTAN_REF_TYPE *ref_type,
                   ZOLTAN_ID_PTR in_vertex, ZOLTAN_ID_PTR out_vertex, int *ierr)
{
  *ierr = ZOLTAN_OK;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int get_partition(void *data, int num_gid_entries, int num_lid_entries,
                  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *ierr)
{
  ELEM_INFO *elem;
  ELEM_INFO *current_elem;
  MESH_INFO_PTR mesh;
  int idx;
  int gid = num_gid_entries-1;
  int lid = num_lid_entries-1;

  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return -1;
  }

  mesh = (MESH_INFO_PTR) data;
  elem = mesh->elements;
  current_elem = (num_lid_entries 
                    ? &elem[local_id[lid]] 
                    : search_by_global_id(mesh, global_id[gid], &idx));


  *ierr = ZOLTAN_OK;
  return current_elem->my_part;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int get_num_hg_edges(
  void *data, 
  int *ierr)
{
  MESH_INFO_PTR mesh;
  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return -1;
  }

  mesh = (MESH_INFO_PTR) data;
  *ierr = ZOLTAN_OK;
  return mesh->nhedges;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int get_num_hg_pins(
  void *data, 
  int *ierr)
{
  MESH_INFO_PTR mesh;

  if (data == NULL) {
    *ierr = ZOLTAN_FATAL;
    return -1;
  }

  mesh = (MESH_INFO_PTR) data;

  *ierr = ZOLTAN_OK;
  return mesh->hindex[mesh->nhedges];
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int get_hg_edge_list(
  void *data, 
  int num_gid_entries,
  int ewgt_dim,
  int nedges,
  int max_size,
  int *edge_sizes,
  ZOLTAN_ID_PTR edge_verts,
  int *edge_procs,
  float *edge_weights
)
{
  MESH_INFO_PTR mesh;
  int i, j;
  int tmp;
  int gid = num_gid_entries-1;
  int *hindex;
  int Proc;
  int ierr = ZOLTAN_OK;

  MPI_Comm_rank(MPI_COMM_WORLD, &Proc);

  if (data == NULL) {
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  mesh = (MESH_INFO_PTR) data;
  hindex = mesh->hindex;
  if (nedges != mesh->nhedges) {
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  for (i = 0; i < nedges; i++) {
    tmp = hindex[i+1] - hindex[i];
    edge_sizes[i] = tmp;
    for (j = hindex[i]; j < hindex[i+1]; j++) {
      edge_verts[gid + j*num_gid_entries] = mesh->hvertex[j];
      edge_procs[j] = Proc;  /* KDD -- temporary; will have to do something
                                KDD -- different in parallel. */
    }
    for (j = 0; j < ewgt_dim; j++) {
      if (mesh->hewgt_dim >= j)
        edge_weights[j + i*ewgt_dim] = mesh->hewgts[j + i*mesh->hewgt_dim];
      else
        edge_weights[j + i*ewgt_dim] = 1.0;
    }
  }

End:
  return ierr;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
ELEM_INFO *search_by_global_id(MESH_INFO *mesh, int global_id, int *idx)
{
/*
 * Function that searchs for an element based upon its global ID.
 * This function does not provide the most efficient implementation of
 * the query functions; more efficient implementation uses local IDs
 * to directly access element info.  However, this function is useful
 * for testing Zoltan when the number of entries in a local ID 
 * (NUM_LID_ENTRIES) is zero.
 */

int i;
ELEM_INFO *elem, *found_elem = NULL;


  elem = mesh->elements;

  for (i = 0; i < mesh->elem_array_len; i++)
    if (elem[i].globalID == global_id) {
      found_elem = &elem[i];
      *idx = i;
      break;
    }
  
  return(found_elem);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
