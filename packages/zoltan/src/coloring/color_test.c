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


/*****************************************************************************/
/*  Parameters structure for Color method.  Used in  */
/*  Zoltan_Color_Set_Param and Zoltan_Color.         */
static PARAM_VARS Color_params[] = {
		  { "COLORING_PROBLEM",   NULL, "STRING", 0 },
		  { "SUPERSTEP_SIZE",     NULL, "INT", 0},
		  { "COMM_PATTERN",       NULL, "CHAR", 0 },
		  { "VERTEX_VISIT_ORDER", NULL, "CHAR", 0 },
		  { "COLORING_METHOD",    NULL, "CHAR", 0},
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
  int *color=NULL;

  int *adjproc=NULL, *xadj=NULL;
  ZOLTAN_GNO_TYPE gvtx;                         /* number of global vertices */
  ZOLTAN_GNO_TYPE *vtxdist=NULL, *adjncy=NULL;
  ZG graph;
  ZOLTAN_GNO_TYPE *requested_GNOs = NULL;  /* Return GNOs of 
                                              the requested GIDs.  */
  int *loc_partialD2 = NULL;    /* local binary array showing which vertices to be colored */
  int *partialD2 = NULL;        /* global binary array showing which vertices to be colored */
  struct Zoltan_DD_Struct *dd_color;  /* DDirectory for colors */
  ZOLTAN_GNO_TYPE *local_GNOs = NULL;
  ZOLTAN_GNO_TYPE *global_GNOs = NULL;
  
  memset (&graph, 0, sizeof(ZG));

  /* PARAMETER SETTINGS */
  Zoltan_Bind_Param(Color_params, "COLORING_PROBLEM",   (void *) &coloring_problemStr);
  Zoltan_Bind_Param(Color_params, "SUPERSTEP_SIZE",     (void *) &ss);
  Zoltan_Bind_Param(Color_params, "COMM_PATTERN",       (void *) &comm_pattern);
  Zoltan_Bind_Param(Color_params, "VERTEX_VISIT_ORDER", (void *) &coloring_order);
  Zoltan_Bind_Param(Color_params, "COLORING_METHOD",    (void *) &coloring_method);

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


  requested_GNOs = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC(num_obj *
                                                     sizeof(ZOLTAN_GNO_TYPE));
  Zoltan_ZG_Build (zz, &graph, 0, 1, num_obj, 
                   global_ids, requested_GNOs);

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
  if (vtxdist[zz->Num_Proc] && !(color = (int *) ZOLTAN_CALLOC(vtxdist[zz->Num_Proc], sizeof(int))))
      MEMORY_ERROR;

  if (nvtx && !(local_GNOs = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC(nvtx * sizeof(ZOLTAN_GNO_TYPE))))
      MEMORY_ERROR;
  for (i=0; i<nvtx; ++i)
      local_GNOs[i] = i+vtxdist[zz->Proc];
  if (vtxdist[zz->Num_Proc] && !(global_GNOs = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC(vtxdist[zz->Num_Proc] * sizeof(ZOLTAN_GNO_TYPE))))
      MEMORY_ERROR;
  for (i=0; i<vtxdist[zz->Num_Proc]; ++i)
      global_GNOs[i] = i;

  ierr = Zoltan_DD_Create (&dd_color, zz->Communicator, 
                           sizeof(ZOLTAN_GNO_TYPE)/sizeof(ZOLTAN_ID_TYPE), 0, 0, 0, 0);
  if (ierr != ZOLTAN_OK)
      ZOLTAN_COLOR_ERROR(ierr, "Cannot construct DDirectory.");
  /* Put req obs with 1 but first inialize the rest with 0 */
  ierr = Zoltan_DD_Update (dd_color, (ZOLTAN_ID_PTR)local_GNOs, NULL,
                           NULL, color, nvtx);
  if (ierr != ZOLTAN_OK)
      ZOLTAN_COLOR_ERROR(ierr, "Cannot update DDirectory.");
  ierr = Zoltan_DD_Update (dd_color, (ZOLTAN_ID_PTR)requested_GNOs, NULL,
                           NULL, color_exp, num_obj);
  if (ierr != ZOLTAN_OK)
      ZOLTAN_COLOR_ERROR(ierr, "Cannot update DDirectory.");
  /* Get requested colors from the DDirectory. */
  ierr = Zoltan_DD_Find (dd_color, (ZOLTAN_ID_PTR)global_GNOs, NULL, NULL,
                         color, vtxdist[zz->Num_Proc], NULL);

  if (ierr != ZOLTAN_OK)
      ZOLTAN_COLOR_ERROR(ierr, "Cannot find object in DDirectory.");
  /* Free DDirectory */
  Zoltan_DD_Destroy(&dd_color);
  ZOLTAN_FREE(&local_GNOs);
  ZOLTAN_FREE(&global_GNOs);


  if (coloring_problem == 'P' || coloring_problem == '2') {
      if (vtxdist[zz->Num_Proc] && !(partialD2 = (int *) ZOLTAN_CALLOC(vtxdist[zz->Num_Proc], sizeof(int))))
	  MEMORY_ERROR;
      if (vtxdist[zz->Num_Proc] && !(loc_partialD2 = (int *) ZOLTAN_CALLOC(vtxdist[zz->Num_Proc], sizeof(int))))
	  MEMORY_ERROR;

      if (coloring_problem == 'P') {
          for (i=0; i<num_obj; ++i) {
              int gno=requested_GNOs[i];
              loc_partialD2[gno] = 1;           
          }

          MPI_Allreduce(loc_partialD2, partialD2, vtxdist[zz->Num_Proc], MPI_INT, MPI_LOR, zz->Communicator);
      } else {
          for (i=0; i<vtxdist[zz->Num_Proc]; ++i)
              partialD2[i] = 1;
      }      
  }

  
  /* Check if there is an error in coloring */
  if (coloring_problem == '1') {
      for (i=0; i<nvtx; i++) {
          int gno = i + (int)vtxdist[zz->Proc];
          if (color[gno] <= 0) { /* object i is not colored properly */
              ierr = ZOLTAN_FATAL;
              printf("Error in coloring! u:%d, cu:%d\n", gno, color[gno]);               
              break;
          }
          for (j = xadj[i]; j < xadj[i+1]; ++j) {
              int v = (int)adjncy[j];
              if (color[gno] == color[v]) { /* neighbors have the same color */
                  ierr = ZOLTAN_FATAL;
                  printf("Error in coloring! u:%d, v:%d, cu:%d, cv:%d\n", gno, v, color[gno], color[v]); 
                  break; 
              }
          }
          if (ierr == ZOLTAN_FATAL)
              break;
      }
  } else if (coloring_problem == '2' || coloring_problem == 'P') {
      for (i=0; i<nvtx; i++) {
          int gno = i + (int)vtxdist[zz->Proc];
          if (partialD2[gno] && color[gno] <= 0) { /* object i is not colored properly */
              ierr = ZOLTAN_FATAL;
              printf("Error in coloring! u:%d, cu:%d\n", gno, color[gno]); 
              break;
          }
          for (j = xadj[i]; j < xadj[i+1]; ++j) {
              int v = (int)adjncy[j], k;
              if (partialD2[v] && color[v] <= 0) {
                  ierr = ZOLTAN_FATAL;
                  printf("Error in coloring! d1-neigh: u:%d, v:%d, cu:%d, cv:%d  pu:%d pv:%d\n", gno, v, color[gno], color[v], partialD2[gno], partialD2[v]);
                  break;
              }
              if (partialD2[gno] && partialD2[v] && color[gno] == color[v]) { /* d-1 neighbors have the same color */
                  ierr = ZOLTAN_FATAL;
                  printf("Error in coloring! d1-neigh: u:%d, v:%d, cu:%d, cv:%d  pu:%d pv:%d\n", gno, v, color[gno], color[v], partialD2[gno], partialD2[v]); 
                  break;
              }
              for (k = j+1; k < xadj[i+1]; ++k) {
                  int w = (int)adjncy[k];
                  if (partialD2[v] && partialD2[w] && color[v] == color[w]) { /* d-2 neighbors have the same color */
                      ierr = ZOLTAN_FATAL;
                      printf("Error in coloring! d2-neigh: v:%d, w:%d, cv:%d, cw:%d   pv:%d pw:%d\n", v, w, color[v], color[w], partialD2[v], partialD2[w]); 
                      break;

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
  ZOLTAN_FREE(&requested_GNOs);
  ZOLTAN_FREE(&partialD2);
  ZOLTAN_FREE(&loc_partialD2);

  return ierr;
}
