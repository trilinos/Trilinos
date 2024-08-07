// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include <ctype.h>
#include "zz_const.h"
#include "zz_util_const.h"
#include "all_allo_const.h"
#include "params_const.h"
#include "graph.h"

int Zoltan_Build_Graph(ZZ *zz, int *graph_type, int check_graph,
		       int num_obj, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
		       int obj_wgt_dim, int * edge_wgt_dim,
		       ZOLTAN_GNO_TYPE **vtxdist, int **xadj, ZOLTAN_GNO_TYPE **adjncy, float **ewgts,
		       int **adjproc)
{
  int ierr = ZOLTAN_OK;
  int local;

  int my_num_obj;
  ZOLTAN_GNO_TYPE glb_obj;
  int my_obj_wgt_dim;
  ZG graph;

  local = IS_LOCAL_GRAPH(*graph_type);
  ierr = Zoltan_ZG_Build (zz, &graph, local, 0,0,NULL,NULL); /* Normal graph */
  ierr = Zoltan_ZG_Export (zz, &graph,
			   &glb_obj, &my_num_obj, &my_obj_wgt_dim, edge_wgt_dim,
			   vtxdist, xadj, adjncy, adjproc,
			   ewgts, NULL);

  graph.mtx.mtx.pinwgt = NULL;
  graph.mtx.mtx.ystart = NULL;
  graph.mtx.mtx.yend = NULL;
  graph.mtx.mtx.pinGNO = NULL;
  graph.mtx.dist_y = NULL;
  Zoltan_ZG_Free(&graph);

  return (ierr);
}


/**************************************************************************/

int Zoltan_Get_Num_Edges_Per_Obj(
  ZZ *zz,
  int num_obj,
  ZOLTAN_ID_PTR global_ids,
  ZOLTAN_ID_PTR local_ids,
  int **edges_per_obj,
  int *max_edges,
  int *num_edges
)
{
/* Calls ZOLTAN_NUM_EDGE_FN or ZOLTAN_NUM_EDGE_MULTI_FN to obtain number
 * of edges per object.
 * Returns number of edges per object in array edges_per_obj.
 * Computes max edges per obj and total edges per obj.
 */
  static char *yo = "Zoltan_Get_Num_Edges_Per_Obj";
int ierr = ZOLTAN_OK;
int i;
int nedges;
int num_gid_entries = zz->Num_GID;
int num_lid_entries = zz->Num_LID;
ZOLTAN_ID_PTR lid;

  *max_edges = *num_edges = 0;
  if (num_obj) {

    *edges_per_obj = (int *) ZOLTAN_MALLOC(num_obj * sizeof(int));
    if (*edges_per_obj == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }

    if (zz->Get_Num_Edges_Multi) {
      zz->Get_Num_Edges_Multi(zz->Get_Num_Edges_Multi_Data,
                              num_gid_entries, num_lid_entries, num_obj,
                              global_ids, local_ids, *edges_per_obj, &ierr);
      if (ierr) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error in Get_Num_Edges_Multi.");
        goto End;
      }

      for (i = 0; i < num_obj; i++) {
        nedges = (*edges_per_obj)[i];
        *num_edges += nedges;
        if (nedges > *max_edges) *max_edges = nedges;
      }
    }
    else {
      for (i=0; i< num_obj; i++) {
        lid = (num_lid_entries ? &(local_ids[i*num_lid_entries]) : NULL);
        nedges = zz->Get_Num_Edges(zz->Get_Num_Edges_Data,
                                   num_gid_entries, num_lid_entries,
                                   &(global_ids[i*num_gid_entries]),
                                   lid, &ierr);
        if (ierr) {
          /* Return error */
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error in Get_Num_Edges.");
          goto End;
        }
        *num_edges += nedges;
        if (nedges > *max_edges) *max_edges = nedges;
        (*edges_per_obj)[i] = nedges;
      }
    }
  }

End:
  return ierr;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
