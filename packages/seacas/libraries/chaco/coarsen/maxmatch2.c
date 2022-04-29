/*
 * Copyright(C) 1999-2021 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "defs.h"    // for FALSE, TRUE
#include "smalloc.h" // for sfree, smalloc
#include "structs.h" // for vtx_data

/* Find a maximal matching in a graph using simple greedy algorithm. */
/* Randomly permute vertices, and then have each select first unmatched */
/* neighbor. */

int maxmatch2(struct vtx_data **graph,      /* array of vtx data for graph */
              int               nvtxs,      /* number of vertices in graph */
              int              *mflag,      /* flag indicating vtx selected or not */
              int               using_ewgts /* are edge weights being used? */
)
{
  extern int HEAVY_MATCH; /* encourage heavy matching edges? */
  void       randomize(int *array, int n);

  /* First, randomly permute the vertices. */
  int *order = smalloc((nvtxs + 1) * sizeof(int));
  int *iptr  = order;
  int *jptr  = mflag;
  for (int i = 1; i <= nvtxs; i++) {
    *(++iptr) = i;
    *(++jptr) = 0;
  }
  randomize(order, nvtxs);

  int nmerged = 0;

  if (!using_ewgts || !HEAVY_MATCH) { /* Simple greedy approach. */
    for (int i = 1; i <= nvtxs; i++) {
      int vtx = order[i];
      if (mflag[vtx] == 0) { /* Not already matched. */
        /* Find first unmatched neighbor of vtx. */
        int matched = FALSE;
        for (int j = 1; j < graph[vtx]->nedges && !matched; j++) {
          int neighbor = graph[vtx]->edges[j];
          if (mflag[neighbor] == 0) { /* Found one! */
            mflag[vtx]      = neighbor;
            mflag[neighbor] = vtx;
            matched         = TRUE;
            nmerged++;
          }
        }
      }
    }
  }

  else { /* Find heavy edge to match */
    for (int i = 1; i <= nvtxs; i++) {
      int vtx = order[i];
      if (mflag[vtx] == 0) { /* Not already matched. */
        /* Find heaviest unmatched neighbor of vtx. */
        int   jsave    = 0;
        float ewgt_max = 0;
        for (int j = 1; j < graph[vtx]->nedges; j++) {
          int neighbor = graph[vtx]->edges[j];
          if (mflag[neighbor] == 0 && graph[vtx]->ewgts[j] > ewgt_max) {
            ewgt_max = graph[vtx]->ewgts[j];
            jsave    = j;
          }
        }
        if (jsave > 0) {
          int neighbor    = graph[vtx]->edges[jsave];
          mflag[vtx]      = neighbor;
          mflag[neighbor] = vtx;
          nmerged++;
        }
      }
    }
  }

  sfree(order);
  return (nmerged);
}
