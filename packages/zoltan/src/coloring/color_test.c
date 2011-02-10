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
#include "graph.h"
#include "all_allo_const.h"
#include "third_library_const.h"


/*****************************************************************************/
/*  Parameters structure for Color method.  Used in  */
/*  Zoltan_Color_Set_Param and Zoltan_Color.         */
static PARAM_VARS Color_params[] = {
		  { "COLORING_PROBLEM", NULL, "STRING", 0 },
		  { "SUPERSTEP_SIZE",   NULL, "INT", 0},
		  { "COMM_PATTERN",     NULL, "CHAR", 0 },
		  { "COLORING_ORDER",   NULL, "CHAR", 0 },
		  { "COLORING_METHOD",  NULL, "CHAR", 0},
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
    int *color_exp            /* Input: Colors assigned to local vertices */
)
{
  static char *yo = "color_test_fn";
  int nvtx = num_obj;               /* number of vertices */
  int i, j;
  int ierr = ZOLTAN_OK;
  int ferr = ZOLTAN_OK;             /* final error signal */
  char coloring_problem;   /* Input: which coloring to perform;
			   currently only supports D1, D2 coloring and variants */
  char coloring_problemStr[MAX_PARAM_STRING_LEN]; /* string version coloring problem name */
  int ss=100;
  char comm_pattern='S', coloring_order='I', coloring_method='F';
  int comm[2],gcomm[2];
  int *color=NULL, *reccnt=NULL, *tmp=NULL;

  int *adjproc=NULL, *xadj=NULL;
  ZOLTAN_GNO_TYPE gvtx;                         /* number of global vertices */
  ZOLTAN_GNO_TYPE *vtxdist=NULL, *adjncy=NULL;
  ZG graph;
  int request_GNOs = 0;                    /* Flag indicating calling code 
                                              needs translation of extra GIDs
                                              to GNOs; partial 2D coloring
                                              needs this feature. */
  int num_requested = 0;                   /* Local # of GIDs needing 
                                              translation to GNOs. */
  ZOLTAN_ID_PTR requested_GIDs = NULL;     /* Calling code requests the 
                                              GNOs for these GIDs */
  ZOLTAN_GNO_TYPE *requested_GNOs = NULL;  /* Return GNOs of 
                                              the requested GIDs.  */

  memset (&graph, 0, sizeof(ZG));

  /* PARAMETER SETTINGS */
  Zoltan_Bind_Param(Color_params, "COLORING_PROBLEM", (void *) &coloring_problemStr);
  Zoltan_Bind_Param(Color_params, "SUPERSTEP_SIZE", (void *) &ss);
  Zoltan_Bind_Param(Color_params, "COMM_PATTERN", (void *) &comm_pattern);
  Zoltan_Bind_Param(Color_params, "COLORING_ORDER", (void *) &coloring_order);
  Zoltan_Bind_Param(Color_params, "COLORING_METHOD", (void *) &coloring_method);

  strncpy(coloring_problemStr, "distance-1", MAX_PARAM_STRING_LEN);

  Zoltan_Assign_Param_Vals(zz->Params, Color_params, zz->Debug_Level, zz->Proc,
			   zz->Debug_Proc);

  /* Check validity of parameters - they should be consistent with Zoltan_Color */
  if (!strcasecmp(coloring_problemStr, "distance-1"))
      coloring_problem = '1';
  else if (!strcasecmp(coloring_problemStr, "distance-2"))
      coloring_problem = '2';
  else if (!strcasecmp(coloring_problemStr, "partial-distance-2")
      || !strcasecmp(coloring_problemStr, "bipartite"))
      coloring_problem = 'P';
  else {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Unknown coloring requested. Using Distance-1 coloring.");
      coloring_problem = '1';
  }
  if (ss == 0) {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Invalid superstep size. Using default value 100.");
      ss = 100;
  }
  if (comm_pattern != 'S' && comm_pattern != 'A') {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Invalid communication pattern. Using synchronous communication (S).");
      comm_pattern = 'S';
  }
  if (comm_pattern == 'A' && (coloring_problem == '2' || coloring_problem == 'P')) {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Asynchronous communication pattern is not implemented for distance-2 coloring and its variants. Using synchronous communication (S).");
      comm_pattern = 'S';
  }
  if (coloring_order != 'I' && coloring_order != 'B' && coloring_order != 'U' && coloring_order != 'L' && coloring_order != 'N' && coloring_order != 'S') {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Invalid coloring order. Using internal first coloring order (I).");
      coloring_order = 'I';
  }
  if (coloring_order == 'U' && (coloring_problem == '2' || coloring_problem == 'P')) {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Interleaved coloring order is not implemented for distance-2 coloring and its variants. Using internal first coloring order (I).");
      coloring_order = 'I';
  }
  if (coloring_method !='F') {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Invalid coloring method. Using first fit method (F).");
      coloring_method = 'F';
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
      ZOLTAN_COLOR_ERROR(ZOLTAN_FATAL, "Color argument is NULL. Please give colors of local vertices.");


  if (coloring_problem == 'P') {
    /* Partial distance-2 coloring requires additional translate of input
     * GIDs to GNOs, using same mapping as in graph.  Pass additional arguments
     * to Zoltan_ZG_Build. 
     */
    request_GNOs = 1;
    num_requested = num_obj;
    requested_GIDs = global_ids;
    requested_GNOs = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC(num_requested *
                                                       sizeof(ZOLTAN_GNO_TYPE));
  }

  Zoltan_ZG_Build (zz, &graph, 0, request_GNOs, num_requested, 
                   requested_GIDs, requested_GNOs);

  Zoltan_ZG_Export (zz, &graph,
		    &gvtx, &nvtx, NULL, NULL, 
                    &vtxdist, &xadj, &adjncy, &adjproc,
		    NULL, NULL);

  if (gvtx > (ZOLTAN_GNO_TYPE)INT_MAX){
    if (zz->Proc == 0){
      fprintf(stderr,
      "Zoltan_Color_Test assumes number of vertices (%ld) is less than INT_MAX\n",gvtx);
    }
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  /* Exchange global color information */
  color = (int *) ZOLTAN_MALLOC(vtxdist[zz->Num_Proc] * sizeof(int));
  reccnt = (int *) ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int));
  tmp = (int *)ZOLTAN_MALLOC(sizeof(int) * zz->Num_Proc);
  if (!color || !reccnt || !tmp)
      MEMORY_ERROR;

  for (i=0; i<zz->Num_Proc; i++){
      reccnt[i] = (int)(vtxdist[i+1]-vtxdist[i]);
      tmp[i] = (int)vtxdist[i];
  }

  /*
   * MPI_Allgatherv requires that offsets be ints. So this won't work if global number
   *   of vertices exceeds 2^32.
   */

  MPI_Allgatherv(color_exp, nvtx, MPI_INT, color, reccnt, tmp, MPI_INT, zz->Communicator);

  ZOLTAN_FREE(&tmp);

  /* Check if there is an error in coloring */
  if (coloring_problem == '1') {
      for (i=0; i<nvtx; i++) {
          int gno = i + (int)vtxdist[zz->Proc];
          if (color[gno] <= 0) { /* object i is not colored properly */
              ierr = ZOLTAN_FATAL;
              break;
              /* printf("Error in coloring! u:%d, cu:%d\n", gno, color[gno]); */
          }
          for (j = xadj[i]; j < xadj[i+1]; ++j) {
              int v = (int)adjncy[j];
              if (color[gno] == color[v]) { /* neighbors have the same color */
                  ierr = ZOLTAN_FATAL;
                  break;
                  /* printf("Error in coloring! u:%d, v:%d, cu:%d, cv:%d\n", gno, v, color[gno], color[v]); */
              }
          }
          if (ierr == ZOLTAN_FATAL)
              break;
      }
  } else if (coloring_problem == '2') {
      for (i=0; i<nvtx; i++) {
          int gno = i + (int)vtxdist[zz->Proc];
          if (color[gno] <= 0) { /* object i is not colored properly */
              ierr = ZOLTAN_FATAL;
              break;
              /* printf("Error in coloring! u:%d, cu:%d\n", gno, color[gno]); */
          }
          for (j = xadj[i]; j < xadj[i+1]; ++j) {
              int v = (int)adjncy[j], k;
              if (color[gno] == color[v]) { /* d-1 neighbors have the same color */
                  ierr = ZOLTAN_FATAL;
                  break;
                  /* printf("Error in coloring! d1-neigh: u:%d, v:%d, cu:%d, cv:%d\n", gno, v, color[gno], color[v]); */
              }
              for (k = j+1; k < xadj[i+1]; ++k) {
                  int w = (int)adjncy[k];
                  if (color[v] == color[w]) { /* d-2 neighbors have the same color */
                      ierr = ZOLTAN_FATAL;
                      break;
                      /* printf("Error in coloring! d2-neigh: v:%d, w:%d, cv:%d, cw:%d\n", v, w, color[v], color[w]); */
                  }

              }
          }
          if (ierr == ZOLTAN_FATAL)
              break;
      }      
  } else 
      ZOLTAN_COLOR_ERROR(ZOLTAN_WARN, "Zoltan_Color_Test is only implemented for distance-1 and distance-2 coloring. Unknown coloring, skipping verification.");      
  
 End:
  if (ierr==ZOLTAN_FATAL)
      ierr = 2;
  MPI_Allreduce(&ierr, &ferr, 1, MPI_INT, MPI_MAX, zz->Communicator);
  if (ferr == 2)
      ierr = ZOLTAN_FATAL;
  else
      ierr = ZOLTAN_OK;

  Zoltan_ZG_Free (&graph);
  ZOLTAN_FREE(&adjproc);
  ZOLTAN_FREE(&color);
  ZOLTAN_FREE(&reccnt);
  ZOLTAN_FREE(&requested_GNOs);

  return ierr;
}
