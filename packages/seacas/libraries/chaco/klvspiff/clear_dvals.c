/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "structs.h" // for vtx_data

void clear_dvals(struct vtx_data **graph,      /* data structure for graph */
                 int               nvtxs,      /* number of vtxs in graph */
                 int *             ldvals,     /* d-values for each transition */
                 int *             rdvals,     /* d-values for each transition */
                 int *             bspace,     /* list of activated vertices */
                 int               list_length /* number of activated vertices */
)
{
  int *edges;    /* loops through edge list */
  int  vtx;      /* vertex in bspace */
  int  neighbor; /* neighbor of vtx */
  int  i, j;     /* loop counters */

  if (list_length > .05 * nvtxs) { /* Do it directly. */
    for (i = 1; i <= nvtxs; i++) {
      ldvals[i] = rdvals[i] = 0;
    }
  }
  else { /* Do it more carefully */
    for (i = 0; i < list_length; i++) {
      vtx = bspace[i];
      if (vtx < 0) {
        vtx = -vtx;
      }
      ldvals[vtx] = rdvals[vtx] = 0;
      edges                     = graph[vtx]->edges;
      for (j = graph[vtx]->nedges - 1; j; j--) {
        neighbor         = *(++edges);
        ldvals[neighbor] = rdvals[neighbor] = 0;
      }
    }
  }
}
