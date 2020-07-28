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
/* neighbor.  Choose the neighbor which results in coarsened vertex of */
/* minimum degree. */

int maxmatch9(struct vtx_data **graph,      /* array of vtx data for graph */
              int               nvtxs,      /* number of vertices in graph */
              int *             mflag,      /* flag indicating vtx selected or not */
              int               using_ewgts /* are edge weights being used? */
)
{
  extern int HEAVY_MATCH;         /* use heavy-edge matching? */
  int *      order;               /* random ordering of vertices */
  int *      neighbors;           /* scatter array for neighbor list */
  int *      iptr, *jptr;         /* loops through integer arrays */
  float      ewgt;                /* edge weight */
  int        save;                /* neighbor vertex if only one active */
  int        vtx;                 /* vertex to process next */
  int        neighbor;            /* neighbor of a vertex */
  int        best;                /* best match found so far */
  float      same = 0, best_same; /* maximum # neighbors in common so far */
  float      best_ewgt;           /* edge weight of possible matching edge */
  int        nmerged;             /* number of edges in matching */
  int        i, j, k;             /* loop counters */

  void randomize();

  /* First, randomly permute the vertices. */
  neighbors = smalloc((nvtxs + 1) * sizeof(int));
  iptr = order = smalloc((nvtxs + 1) * sizeof(int));
  jptr         = mflag;
  for (i = 1; i <= nvtxs; i++) {
    *(++iptr)    = i;
    *(++jptr)    = 0;
    neighbors[i] = 0;
  }
  randomize(order, nvtxs);

  /*    if (!using_ewgts || !HEAVY_MATCH) { */

  nmerged = 0;
  ewgt    = 0;
  for (i = 1; i <= nvtxs; i++) {
    vtx = order[i];
    if (mflag[vtx] == 0) { /* Not already matched. */
      /* Add up sum of edge weights of neighbors. */
      save = -1;
      for (j = 1; j < graph[vtx]->nedges; j++) {
        neighbor            = graph[vtx]->edges[j];
        neighbors[neighbor] = i;
        if (mflag[neighbor] == 0) {
          /* Set flag for single possible neighbor. */
          if (save == -1) {
            save = neighbor;
          }
          else {
            save = 0;
          }
        }
        else {
          neighbors[mflag[neighbor]] = i;
        }
      }

      if (save != -1) { /* Does vertex have contractible edges? */
        nmerged++;
        if (save > 0) { /* Only one neighbor, easy special case. */
          mflag[vtx]  = save;
          mflag[save] = vtx;
        }
        else { /* Merge with best neighbor */
          best      = 0;
          best_same = -1;
          best_ewgt = -1;
          for (j = 1; j < graph[vtx]->nedges; j++) {
            neighbor = graph[vtx]->edges[j];
            if (mflag[neighbor] == 0) {
              if (using_ewgts && HEAVY_MATCH) {
                ewgt = graph[vtx]->ewgts[j];
              }
              if (ewgt > best_ewgt) {
                best      = neighbor;
                best_same = same;
                best_ewgt = ewgt;
              }
              else if (ewgt == best_ewgt) {
                /* break ties by larger same value */
                if (best_same == -1) {
                  /* Compute same value for current best vtx. */
                  best_same = 0;
                  for (k = 1; k < graph[best]->nedges; k++) {
                    if (neighbors[graph[best]->edges[k]] == i) {
                      if (using_ewgts) {
                        best_same += graph[best]->ewgts[k];
                      }
                      else {
                        best_same += 1;
                      }
                    }
                  }
                }
                /* Now compute same value for candidate vtx. */
                same = 0;
                for (k = 1; k < graph[neighbor]->nedges; k++) {
                  if (neighbors[graph[neighbor]->edges[k]] == i) {
                    if (using_ewgts) {
                      same += graph[neighbor]->ewgts[k];
                    }
                    else {
                      same += 1;
                    }
                  }
                }

                if (same > best_same) {
                  best      = neighbor;
                  best_same = same;
                  best_ewgt = ewgt;
                }
              }
            }
          }
          mflag[vtx]  = best;
          mflag[best] = vtx;
        }
      }
    }
  }

  sfree(order);
  sfree(neighbors);
  return (nmerged);
}
