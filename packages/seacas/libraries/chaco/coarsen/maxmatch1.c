/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "defs.h"    // for FALSE, TRUE
#include "structs.h" // for vtx_data

/* Find a maximal matching in a graph using simple greedy algorithm. */
/* Randomly permute vertices, and then have each select an unmatched */
/* neighbor. */

int maxmatch1(struct vtx_data **graph,      /* array of vtx data for graph */
              int               nvtxs,      /* number of vertices in graph */
              int              *mflag,      /* flag indicating vtx selected or not */
              int               using_ewgts /* are edge weights being used? */
)
{
  extern int HEAVY_MATCH; /* choose heavy edges in matching? */
  float      ewgt_max;    /* largest edge weight seen so far */
  int       *jptr;        /* loops through integer arrays */
  int        vtx;         /* vertex to process next */
  int        neighbor;    /* neighbor of a vertex */
  int        nmerged;     /* number of edges in matching */
  int        matched;     /* is a vertex matched yet? */
  int        jsave;       /* best matching edge found so far */
  int        i, j;        /* loop counters */
  double     drandom(void);

  /* Initialize mflag array. */
  jptr = mflag;
  for (i = nvtxs; i; i--) {
    *(++jptr) = 0;
  }

  nmerged = 0;

  /* Select random starting point in list of vertices. */
  vtx = 1 + drandom() * nvtxs;

  if (!using_ewgts || !HEAVY_MATCH) { /* Choose first neighbor */
    for (i = nvtxs; i; i--) {
      if (mflag[vtx] == 0) { /* Not already matched. */
        /* Select first free edge. */
        matched = FALSE;
        for (j = 1; !matched && j < graph[vtx]->nedges; j++) {
          neighbor = graph[vtx]->edges[j];
          if (mflag[neighbor] == 0) {
            mflag[vtx]      = neighbor;
            mflag[neighbor] = vtx;
            matched         = TRUE;
            nmerged++;
          }
        }
      }

      if (++vtx > nvtxs) {
        vtx = 1;
      }
    }
  }

  else { /* Choose heavy edge neighbor */
    for (i = nvtxs; i; i--) {
      if (mflag[vtx] == 0) { /* Not already matched. */
        /* Select heaviest free edge. */
        jsave    = 0;
        ewgt_max = 0;
        for (j = 1; j < graph[vtx]->nedges; j++) {
          neighbor = graph[vtx]->edges[j];
          if (mflag[neighbor] == 0 && graph[vtx]->ewgts[j] > ewgt_max) {
            ewgt_max = graph[vtx]->ewgts[j];
            jsave    = j;
          }
        }
        if (jsave > 0) {
          neighbor        = graph[vtx]->edges[jsave];
          mflag[vtx]      = neighbor;
          mflag[neighbor] = vtx;
          nmerged++;
        }
      }

      if (++vtx > nvtxs) {
        vtx = 1;
      }
    }
  }

  return (nmerged);
}
