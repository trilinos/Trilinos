/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "params.h"  // for MAXSETS
#include "structs.h" // for vtx_data, bilist
#include <stdio.h>   // for NULL

/* Idea:
   'buckets[i][j]' is a set of buckets to sort moves from i to j.
   listspace[i] is space for lists in buckets[i][j].
   Loop through all nonequal pairs [i][j], taking the first element
   in each list.  Compare them all to find the largest allowed move.
   Make that move, and save it in movelist.
*/

/* Routine slightly streamlined for case with only two sets. */

void bucketsorts_bi(struct vtx_data **graph,       /* graph data structure */
                    int               nvtxs,       /* number of vertices */
                    struct bilist ****buckets,     /* array of lists for bucket sort */
                    struct bilist **  listspace,   /* list data structure for each vertex */
                    int **            dvals,       /* d-values for each vertex for removing */
                    int *             sets,        /* processor each vertex is assigned to */
                    float *           term_wgts[], /* weights for terminal propagation */
                    int               maxdval,     /* maximum possible dvalue for a vertex */
                    int               nsets,       /* number of sets being divided into */
                    int               parity,      /* work in forward or backward direction? */
                    int (*hops)[MAXSETS],          /* hop cost between sets */
                    int *bspace,                   /* indices for randomly ordering vtxs */
                    int  list_length,              /* number of values in bspace to work with */
                    int  npass,                    /* which pass through KL is this? */
                    int  using_ewgts               /* are edge weights being used? */
)
{
  extern int      KL_RANDOM;       /* use randomness in KL? */
  extern int      KL_UNDO_LIST;    /* only sort vertices who have moved. */
  extern double   CUT_TO_HOP_COST; /* ..if so, relative cut/hop importance */
  struct bilist **bptr  = NULL;    /* loops through set of buckets */
  struct bilist * lptr  = NULL;    /* pointer to an element in listspace */
  float *         ewptr = NULL;    /* loops through edge weights */
  float *         twptr = NULL;    /* weights for terminal propagation */
  int *           bsptr = NULL;    /* loops through bspace */
  int *           edges = NULL;    /* edge list for a vertex */
  int             myset;           /* set current vertex belongs to */
  int             other_set;       /* set current vertex doesn't belong to */
  int             set;             /* set that neighboring vertex belongs to */
  int             weight;          /* edge weight for a particular edge */
  int             vtx;             /* vertex in graph */
  int             val;             /* terminal propagation rounded value */
  double          cut_cost;        /* relative cut/hop importance */
  double          hop_cost;        /* relative hop/cut importance */
  int             myhop;           /* hops associated with current vertex */
  int             i, j;            /* loop counters */
  void            randomize(), add2bilist();

  /* For each vertex, compute d-values for each possible transition. */
  /* Then store them in each appropriate bucket. */

  if (npass == 1 || !KL_UNDO_LIST || list_length == nvtxs) {
    /* Empty all the buckets. */
    /* Last clause catches case where lists weren't undone. */
    bptr = buckets[0][1];
    for (i = nsets * (nsets - 1) * (2 * maxdval + 1); i; i--) {
      *bptr++ = NULL;
    }
  }

  /* Randomize the order of the vertices */

  if ((KL_UNDO_LIST && list_length == nvtxs) || !KL_UNDO_LIST) {
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
  cut_cost = hop_cost = 1;
  if (term_wgts[1] != NULL) {
    if (CUT_TO_HOP_COST > 1) {
      cut_cost = CUT_TO_HOP_COST;
    }
    else {
      hop_cost = 1.0 / CUT_TO_HOP_COST;
    }
  }

  weight = cut_cost + .5;

  bsptr = bspace;
  twptr = term_wgts[1];
  for (i = 0; i < list_length; i++) { /* Loop through vertices. */
    vtx       = *bsptr++;
    myset     = sets[vtx];
    other_set = !myset;

    /* Initialize all the preference values. */
    if (twptr != NULL) {
      /* Using terminal propagation.  Round to integer value. */
      if (twptr[vtx] < 0) {
        val = -twptr[vtx] * hop_cost + .5;
        val = -val;
      }
      else {
        val = twptr[vtx] * hop_cost + .5;
      }
      if (myset == 0) {
        dvals[vtx][0] = val;
      }
      else {
        dvals[vtx][0] = -val;
      }
    }
    else {
      dvals[vtx][0] = 0;
    }

    /* First count the neighbors in each set. */
    edges = graph[vtx]->edges;
    if (using_ewgts) {
      ewptr = graph[vtx]->ewgts;
    }
    for (j = graph[vtx]->nedges - 1; j; j--) {
      set = sets[*(++edges)];
      if (set < 0) {
        set = -set - 1;
      }
      if (using_ewgts) {
        weight = *(++ewptr) * cut_cost + .5;
      }
      myhop = hops[myset][set];

      dvals[vtx][0] += weight * (myhop - hops[other_set][set]);
    }

    /* Now add to appropriate buckets. */
    lptr = listspace[0];
    add2bilist(&lptr[vtx], &buckets[myset][other_set][dvals[vtx][0] + maxdval]);
  }
}
