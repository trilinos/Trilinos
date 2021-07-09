/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "smalloc.h" // for sfree, smalloc
#include "structs.h" // for vtx_data

/* Find a maximal matching in a graph using simple greedy algorithm. */
/* Randomly permute vertices, and then have each select an unmatched */
/* neighbor. */

int maxmatch3(struct vtx_data **graph,      /* array of vtx data for graph */
              int               nvtxs,      /* number of vertices in graph */
              int *             mflag,      /* flag indicating vtx selected or not */
              int               using_ewgts /* are edge weights being used? */
)
{
  extern int HEAVY_MATCH; /* pick heavy edges in matching? */
  int *      order;       /* random ordering of vertices */
  int *      iptr, *jptr; /* loops through integer arrays */
  double     prob_sum;    /* sum of probabilities to select from */
  double     val;         /* random value for selecting neighbor */
  float      ewgt;        /* edge weight */
  int        save;        /* neighbor vertex if only one active */
  int        vtx;         /* vertex to process next */
  int        neighbor;    /* neighbor of a vertex */
  int        nmerged;     /* number of edges in matching */
  int        i, j;        /* loop counters */

  double drandom();
  void   randomize();

  /* First, randomly permute the vertices. */
  iptr = order = smalloc((nvtxs + 1) * sizeof(int));
  jptr         = mflag;
  for (i = 1; i <= nvtxs; i++) {
    *(++iptr) = i;
    *(++jptr) = 0;
  }
  randomize(order, nvtxs);

  nmerged = 0;
  if (!using_ewgts || !HEAVY_MATCH) { /* All edges equal. */
    for (i = 1; i <= nvtxs; i++) {
      vtx = order[i];
      if (mflag[vtx] == 0) { /* Not already matched. */
        /* Add up sum of edge weights of neighbors. */
        prob_sum = 0;
        save     = 0;
        for (j = 1; j < graph[vtx]->nedges; j++) {
          neighbor = graph[vtx]->edges[j];
          if (mflag[neighbor] == 0) {
            /* Set flag for single possible neighbor. */
            if (prob_sum == 0) {
              save = neighbor;
            }
            else {
              save = 0;
            }
            prob_sum += 1.0;
          }
        }

        if (prob_sum != 0) { /* Does vertex have contractible edges? */
          nmerged++;
          if (save != 0) { /* Only one neighbor, special case. */
            mflag[vtx]  = save;
            mflag[save] = vtx;
          }
          else { /* Pick randomly neighbor. */
            val      = drandom() * prob_sum * .999999;
            prob_sum = 0;
            for (j = 1; !mflag[vtx]; j++) {
              neighbor = graph[vtx]->edges[j];
              if (mflag[neighbor] == 0) {
                prob_sum += 1.0;
                if (prob_sum >= val) {
                  mflag[vtx]      = neighbor;
                  mflag[neighbor] = vtx;
                }
              }
            }
          }
        }
      }
    }
  }

  else { /* Choose heavy edges preferentially. */
    for (i = 1; i <= nvtxs; i++) {
      vtx = order[i];
      if (mflag[vtx] == 0) { /* Not already matched. */
        /* Add up sum of edge weights of neighbors. */
        prob_sum = 0;
        save     = 0;
        for (j = 1; j < graph[vtx]->nedges; j++) {
          neighbor = graph[vtx]->edges[j];
          if (mflag[neighbor] == 0) {
            /* Set flag for single possible neighbor. */
            if (prob_sum == 0) {
              save = neighbor;
            }
            else {
              save = 0;
            }
            ewgt = graph[vtx]->ewgts[j];
            prob_sum += ewgt;
          }
        }

        if (prob_sum != 0) { /* Does vertex have contractible edges? */
          nmerged++;
          if (save != 0) { /* Only one neighbor, special case. */
            mflag[vtx]  = save;
            mflag[save] = vtx;
          }
          else { /* Pick randomly neighbor, skewed by edge weights. */
            val      = drandom() * prob_sum * .999999;
            prob_sum = 0;
            for (j = 1; !mflag[vtx]; j++) {
              neighbor = graph[vtx]->edges[j];
              if (mflag[neighbor] == 0) {
                ewgt = graph[vtx]->ewgts[j];
                prob_sum += ewgt;
                if (prob_sum >= val) {
                  mflag[vtx]      = neighbor;
                  mflag[neighbor] = vtx;
                }
              }
            }
          }
        }
      }
    }
  }

  sfree(order);
  return (nmerged);
}
