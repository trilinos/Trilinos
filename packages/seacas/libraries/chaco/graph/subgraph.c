/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "structs.h" // for vtx_data
#include <stdio.h>   // for NULL

/* Construct a subgraph of a graph. */
/* This reuses the graph storage, and can be undone. */

void make_subgraph(struct vtx_data **graph,      /* graph data structure */
                   struct vtx_data **subgraph,   /* subgraph data structure */
                   int               subnvtxs,   /* number of vtxs in subgraph */
                   int *             psubnedges, /* ptr to number of edges in subgraph */
                   int *             assignment, /* values designating subgraph inclusion */
                   int               set,        /* assignment value indicating inclusion */
                   int *             glob2loc,   /* mapping from graph to subgraph numbering */
                   int *             loc2glob,   /* mapping from subgraph to graph numbering */
                   int *             degree,     /* degrees of vertices in graph */
                   int               using_ewgts /* are edge weights being used? */
)
{
  struct vtx_data *subgptr = NULL; /* loops through subgraph */
  float *          fptr    = NULL; /* loops through edge weights */
  int *            iptr    = NULL; /* loops through edge list */
  float            tempwgt;        /* weight of vertex being swapped */
  double           ewgtsum;        /* sum of weights of subgraph edges */
  int              subnedges;      /* number of edges in subgraph */
  int              neighbor;       /* neighbor vertex in graph */
  int              newnedges;      /* vertex degree in subgraph */
  int              tempvtx;        /* vertex number being swapped */
  int              i, j;           /* loop counter */

  subnedges = 0;
  for (i = 1; i <= subnvtxs; i++) {
    /* First get the vertices organized properly. */
    subgptr = subgraph[i] = graph[loc2glob[i]];
    newnedges = degree[i] = subgptr->nedges;

    /* Now work on the edges. */
    subgptr->edges[0] = i;
    ewgtsum           = 0;
    /* Move all deactivated edges to the end of the list. */
    iptr = subgptr->edges + 1;
    if (using_ewgts) {
      fptr = subgptr->ewgts + 1;
    }
    for (j = 1; j < newnedges;) {
      neighbor = *iptr;
      if (assignment[neighbor] == set) { /* Keep vertex in edge list. */
        subgptr->edges[j] = glob2loc[neighbor];
        if (using_ewgts) {
          ewgtsum += *fptr++;
        }
        j++;
        iptr++;
      }
      else { /* Move vertex to back of edge list. */
        --newnedges;
        tempvtx                   = subgptr->edges[newnedges];
        subgptr->edges[newnedges] = neighbor;
        *iptr                     = tempvtx;
        if (using_ewgts) {
          tempwgt                   = subgptr->ewgts[newnedges];
          subgptr->ewgts[newnedges] = *fptr;
          *fptr                     = tempwgt;
        }
      }
    }
    subgptr->nedges = newnedges;
    subnedges += newnedges;
    if (using_ewgts) {
      subgptr->ewgts[0] = -ewgtsum;
    }
  }
  *psubnedges = (subnedges - subnvtxs) / 2;
}

/* Undo the construction of the subgraph. */
void remake_graph(struct vtx_data **subgraph,   /* subgraph data structure */
                  int               subnvtxs,   /* number of vtxs in subgraph */
                  int *             loc2glob,   /* mapping from subgraph to graph numbering */
                  int *             degree,     /* degrees of vertices in graph */
                  int               using_ewgts /* are edge weights being used? */
)
{
  struct vtx_data *subgptr; /* loops through subgraph */
  float *          fptr;    /* loops through edge weights */
  int *            iptr;    /* loops through adjacency list */
  double           ewgtsum; /* sum of weights of subgraph edges */
  int              nedges;  /* vertex degree in subgraph */
  int              i, j;    /* loop counter */

  /* For each vertex in subgraph, expand the edge set back out. */
  for (i = 1; i <= subnvtxs; i++) {
    subgptr           = subgraph[i];
    subgptr->edges[0] = loc2glob[i];
    nedges            = subgptr->nedges;
    /* Change edges back to global numbering system. */
    iptr = subgptr->edges;
    for (j = nedges - 1; j; j--) {
      iptr++;
      *iptr = loc2glob[*iptr];
    }
    subgptr->nedges = degree[i];

    /* Now get the diagonal value right. */
    if (using_ewgts) {
      ewgtsum = 0;
      fptr    = subgptr->ewgts;
      for (j = degree[i] - 1; j; j--) {
        ewgtsum += *(++fptr);
      }
      subgptr->ewgts[0] = -ewgtsum;
    }
  }
}
