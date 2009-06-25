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


#include <limits.h>
#include <ctype.h>
#include "zoltan_mem.h"    
#include "zz_const.h"
#include "coloring.h"
#include "g2l_hash.h"
#include "params_const.h"
#include "zz_util_const.h"
#include "third_library_const.h"
#include "all_allo_const.h"


/*****************************************************************************/
/*  Parameters structure for Color method.  Used in  */
/*  Zoltan_Color_Set_Param and Zoltan_Color.         */
static PARAM_VARS Color_params[] = {
                  { "DISTANCE", NULL, "INT", 0 },
                  { "SUPERSTEP_SIZE", NULL, "INT", 0},
                  { "COMM_PATTERN", NULL, "CHAR", 0 },
                  { "COLOR_ORDER", NULL, "CHAR", 0 },
                  { "COLORING_METHOD", NULL, "CHAR", 0},
                  { NULL, NULL, NULL, 0 } };
       
/*****************************************************************************/
/* Interface routine for Graph Coloring Testing */

int Zoltan_Color_Test(
    ZZ *zz,                   /* Zoltan structure */
    int *num_gid_entries,     /* # of entries for a global id */
    int *num_lid_entries,     /* # of entries for a local id */
    int num_obj,              /* Input: number of objects */
    ZOLTAN_ID_PTR global_ids, /* Input: global ids of the vertices */
                              /* The application must allocate enough space */    
    ZOLTAN_ID_PTR local_ids,  /* Input: local ids of the vertices */
                              /* The application must allocate enough space */
    int *color_exp            /* Output: Colors assigned to local vertices */
                              /* The application must allocate enough space */
) 
{
  static char *yo = "color_test_fn";
  indextype *vtxdist, *xadj, *adjncy; /* arrays to store the graph structure */
  int *adjproc;                     
  int *input_parts;                 /* Initial partitions for objects. */
  int nvtx = num_obj;               /* number of vertices */
  float *ewgts, *float_vwgt;        /* weights - not used */
  int obj_wgt_dim, edge_wgt_dim;    /* weight dimensions - not used */
  int graph_type, check_graph;
  int i, j;
  int ierr = ZOLTAN_OK;
  int ferr = ZOLTAN_OK;             /* final error signal */
  int distance=1;                   /* distance for coloring. This test is only implemeted for distance=1 */
  int ss=100;
  char comm_pattern='S', color_order='I', color_method='F';  
  int comm[2],gcomm[2];
  int *color=NULL, *reccnt=NULL;
  
  /* PARAMETER SETTINGS */
  Zoltan_Bind_Param(Color_params, "DISTANCE", (void *) &distance);
  Zoltan_Bind_Param(Color_params, "SUPERSTEP_SIZE", (void *) &ss);
  Zoltan_Bind_Param(Color_params, "COMM_PATTERN", (void *) &comm_pattern);
  Zoltan_Bind_Param(Color_params, "COLOR_ORDER", (void *) &color_order);
  Zoltan_Bind_Param(Color_params, "COLORING_METHOD", (void *) &color_method);
  
  Zoltan_Assign_Param_Vals(zz->Params, Color_params, zz->Debug_Level, zz->Proc,
                           zz->Debug_Proc);

  /* Check validity of parameters - they should be consistent with Zoltan_Color */
  if (distance != 1 && distance != 2) {
      distance = 1;
  }
  if (ss == 0) {
      ss = 100;
  }
  if (comm_pattern != 'S' && comm_pattern != 'A') {
      comm_pattern = 'S';
  }
  if (comm_pattern == 'A' && distance == 2) {
      comm_pattern = 'S';
  }
  if (color_order != 'I' && color_order != 'B' && color_order != 'U') {
      color_order = 'I';
  }  
  if (color_order == 'U' && distance == 2) {
      color_order = 'I';
  }
  if (color_method !='F') {
      color_method = 'F';
  }
  
  /* Compute Max number of array entries per ID over all processors.
     This is a sanity-maintaining step; we don't want different
     processors to have different values for these numbers. */
  comm[0] = zz->Num_GID;
  comm[1] = zz->Num_LID;
  MPI_Allreduce(comm, gcomm, 2, MPI_INT, MPI_MAX, zz->Communicator);
  zz->Num_GID = *num_gid_entries = gcomm[0];
  zz->Num_LID = *num_lid_entries = gcomm[1];
  
  /* Return if this processor is not in the Zoltan structure's
     communicator. */
  if (ZOLTAN_PROC_NOT_IN_COMMUNICATOR(zz))
      return ZOLTAN_OK;


  /* BUILD THE GRAPH */
  /* Check that the user has allocated space for the return args. */
  if (!color_exp)
      ZOLTAN_COLOR_ERROR(ZOLTAN_FATAL, "Output argument is NULL. Please allocate all required arrays before calling this routine.");
  
  /* Initialize all local pointers to NULL. This is necessary
     because we free all non-NULL pointers upon errors. */
  vtxdist = xadj = adjncy = adjproc = NULL;
  ewgts = float_vwgt = NULL;
  input_parts = NULL;
 
  /* Default graph type is GLOBAL. */
  SET_GLOBAL_GRAPH(&graph_type);
  check_graph = 1;
  obj_wgt_dim = 0; /* We do not use weights */
  edge_wgt_dim = 0;

  /* Get object ids and part information */
  ierr = Zoltan_Get_Obj_List(zz, &nvtx, &global_ids, &local_ids,
                             obj_wgt_dim, &float_vwgt, &input_parts);
  if (ierr) { /* Return error */      
      ZOLTAN_COLOR_ERROR(ierr, "Get_Obj_List returned error.");
  }

  /* Build ParMetis data structures, or just get vtxdist. */
  ierr = Zoltan_Build_Graph(zz, &graph_type, check_graph, nvtx,
         global_ids, local_ids, obj_wgt_dim, edge_wgt_dim,
         &vtxdist, &xadj, &adjncy, &ewgts, &adjproc);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_COLOR_ERROR(ierr, "Zoltan_Build_Graph returned error.");
  }

  if (nvtx!=num_obj) {
      ierr = ZOLTAN_FATAL;
      ZOLTAN_COLOR_ERROR(ierr, "Zoltan_Build_Graph returned different number of vertices.");
  }


  /* Exchange global color information */
  color = (int *) ZOLTAN_MALLOC(vtxdist[zz->Num_Proc] * sizeof(int));
  reccnt = (int *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int));
  if (!color || !reccnt)
      MEMORY_ERROR;

  if (distance != 1) 
      ZOLTAN_COLOR_ERROR(ZOLTAN_WARN, "Zoltan_Color_Test is only implemented for distance-1 coloring. Skipping verification.");


  for (i=0; i<zz->Num_Proc; i++)
      reccnt[i] = vtxdist[i+1]-vtxdist[i];
  MPI_Allgatherv(color_exp, nvtx, MPI_INT, color, reccnt, vtxdist, MPI_INT, zz->Communicator);

  /* Check if there is an error in coloring */
  for (i=0; i<nvtx; i++) {
      int gno = i + vtxdist[zz->Proc];
      if (color[gno] == 0) { /* object i is not colored */
          ierr = ZOLTAN_FATAL;
          break;
          /* printf("Error in coloring! u:%d, cu:%d\n", gno, color[gno]); */
      }
      for (j = xadj[i]; j < xadj[i+1]; j++) {
          int v = adjncy[j];
          if (color[gno] == color[v]) { /* neighbors have the same color */
              ierr = ZOLTAN_FATAL;
              break;
              /* printf("Error in coloring! u:%d, v:%d, cu:%d, cv:%d\n", gno, v, color[gno], color[v]); */
          }
      }
      if (ierr == ZOLTAN_FATAL)
          break;
  }

 End:
  if (ierr==ZOLTAN_FATAL)
      ierr = 2;
  MPI_Allreduce(&ierr, &ferr, 1, MPI_INT, MPI_MAX, zz->Communicator);  
  if (ferr == 2)
      ierr = ZOLTAN_FATAL;
  else
      ierr = ZOLTAN_OK;
  
  Zoltan_Multifree(__FILE__, __LINE__, 7, &vtxdist, &xadj, &adjncy, &input_parts, &adjproc, &color, &reccnt);
  
  return ierr;
}

