/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/
#ifndef lint
static char *cvs_ch_driver_util_id = "$Id$";
#endif

#include <stdio.h>
#include <stdlib.h>
#include "lbi_const.h"
#include "all_allo_const.h"
#include "ch_input_const.h"
#include "ch_driver_const.h"

/*
 * Zoltan query functions for the Chaco/ParMetis graph structure.
 * Local id is the local number of a vertex on a processor,
 * while global id is the global vertex number (dense on all procs).
 */

int num_obj_fn (void *data, int *ierr)
{
  Graph *graph;
  int n;

  if (data == NULL){
    *ierr = LB_FATAL;
    n = 0;
  }
  else{
    graph = (Graph *)data;
    n = graph->lnvtxs;
    *ierr = LB_OK;
  }
  return n;
}


void list_obj_fn(void *data, LB_GID *global_ids, LB_LID *local_ids,
                            int *ierr)
{
  Graph *graph;
  int i;

  if (data == NULL)
    *ierr = LB_FATAL;
  else{
    graph = (Graph *)data;
    for (i=0; i<graph->lnvtxs; i++){
      local_ids[i]  = i;
      global_ids[i] = i + graph->vtxdist[graph->myproc];
    }
    *ierr = LB_OK;
  }
}

float obj_weight_fn (void *data, LB_GID global_id, LB_LID local_id,
                                int *ierr)
{
  Graph *graph;
  float w;

  if (data == NULL)
    *ierr = LB_FATAL;
  else{
    graph = (Graph *)data;
    w = graph->vwgts[local_id];
    *ierr = LB_OK;
  }
  return w;
}

int num_edges_fn(void *data, LB_GID global_id, LB_LID local_id,
                            int *ierr)
{
  Graph *graph;
  int nedges;

  if (data == NULL)
    *ierr = LB_FATAL;
  else{
    graph = (Graph *)data;
    nedges = graph->start[local_id+1] - graph->start[local_id];
    *ierr = LB_OK;
  }
  return nedges;
}

void edge_list_fn (void *data, LB_GID global_id, LB_LID local_id,
                             LB_GID *nbor_global_id, int *nbor_procs,
                             int get_ewgts, int *nbor_ewgts, int *ierr)
{
  Graph *graph;
  int i,j,p;
  int offset = 1; /* Lowest numbered vertex in input file. Ugly hack... */

  if (data == NULL)
    *ierr = LB_FATAL;
  else{
    graph = (Graph *)data;
    *ierr = LB_OK;
    for (i=0,j=graph->start[local_id]; j<graph->start[local_id+1]; i++,j++){
      nbor_global_id[i] = graph->adjncy[j] - offset;
      if (get_ewgts)
        nbor_ewgts[i] = graph->ewgts[j];
      /* Determine host proc of current neighbor vertex */
      for (p=0; nbor_global_id[i] >= graph->vtxdist[p+1]; p++);
      nbor_procs[i] = p; 
    }
    *ierr = LB_OK;
  }
}


