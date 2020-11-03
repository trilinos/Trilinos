/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "structs.h" // for vtx_data, bilist
#include <stdio.h>   // for NULL

void bucketsortsv(struct vtx_data **graph,      /* graph data structure */
                  int               nvtxs,      /* number of vertices */
                  struct bilist **  lbuckets,   /* array of lists for bucket sort */
                  struct bilist **  rbuckets,   /* array of lists for bucket sort */
                  struct bilist *   llistspace, /* list data structure for each vertex */
                  struct bilist *   rlistspace, /* list data structure for each vertex */
                  int *             ldvals,     /* d-values for each vertex for removing */
                  int *             rdvals,     /* d-values for each vertex for removing */
                  int *             sets,       /* processor each vertex is assigned to */
                  int               maxdval,    /* maximum possible dvalue for a vertex */
                  int               parity,     /* work in forward or backward direction? */
                  int *             bspace,     /* indices for randomly ordering vtxs */
                  int               list_length /* number of values in bspace to work with */
)
{
  extern int      KL_RANDOM;    /* use randomness in KL? */
  extern int      KL_UNDO_LIST; /* only sort vertices who have moved. */
  struct bilist **lbptr;        /* loops through set of buckets */
  struct bilist **rbptr;        /* loops through set of buckets */
  int *           bsptr;        /* loops through bspace */
  int *           edges;        /* edge list for a vertex */
  int             left_weight;  /* my neighbors in 0 set */
  int             right_weight; /* my neighbors in 1 set */
  int             vtx;          /* vertex in graph */
  int             neighbor;     /* neighbor of vertex */
  int             set;          /* set that neighboring vertex belongs to */
  int             i, j;         /* loop counters */
  void            randomize(), add2bilist();

  /* For each vertex, compute d-values and store in buckets. */

  /* Empty all the buckets. */
  rbptr = lbuckets;
  lbptr = rbuckets;
  for (i = (2 * maxdval + 1); i; i--) {
    *lbptr++ = NULL;
    *rbptr++ = NULL;
  }

  /* Randomize the order of the vertices */

  if ((KL_UNDO_LIST && list_length == nvtxs) || (!KL_UNDO_LIST && !KL_RANDOM) || list_length == 0) {
    /* Don't need to reorder if about to randomize. */
    list_length = nvtxs;
    bsptr       = bspace;
    if (parity) {
      for (i = 1; i <= nvtxs; i++) {
        *bsptr++ = i;
      }
    }
    else {
      for (i = nvtxs; i; i--) {
        *bsptr++ = i;
      }
    }
  }
  if (KL_RANDOM) {
    randomize(bspace - 1, list_length);
  }

  /* Now compute d-vals by seeing which sets neighbors belong to. */

  bsptr = bspace;
  for (i = 0; i < list_length; i++) { /* Loop through vertices. */
    vtx = *bsptr++;
    if (sets[vtx] == 2) { /* Only worry about separator vertices. */

      /* Initialize all the preference values. */
      left_weight = right_weight = 0;

      /* First count the neighbors in each set. */
      edges = graph[vtx]->edges;
      for (j = graph[vtx]->nedges - 1; j; j--) {
        neighbor = *(++edges);
        set      = sets[neighbor];
        if (set < 0) {
          set = -set - 1;
        }
        if (set == 0) {
          left_weight += graph[neighbor]->vwgt;
        }
        else if (set == 1) {
          right_weight += graph[neighbor]->vwgt;
        }
      }

      ldvals[vtx] = graph[vtx]->vwgt - right_weight;
      rdvals[vtx] = graph[vtx]->vwgt - left_weight;

      /* Now add to appropriate buckets. */
      add2bilist(&llistspace[vtx], &lbuckets[ldvals[vtx] + maxdval]);
      add2bilist(&rlistspace[vtx], &rbuckets[rdvals[vtx] + maxdval]);
    }
  }
}
