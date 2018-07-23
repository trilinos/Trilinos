/*
 * Copyright (c) 2005-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include "params.h"  // for MAXSETS
#include "structs.h" // for bilist
#include <stdio.h>   // for NULL

/* Perform KL between two sets. */
int kl_refine(struct vtx_data **graph,      /* graph data structure */
              struct vtx_data **subgraph,   /* space for subgraph to refine */
              struct bilist *   set_list,   /* lists of vtxs in each set */
              struct bilist *   vtx_elems,  /* start of storage for lists */
              int *             new_assign, /* set assignments for all vertices */
              int set1, int set2,           /* two sets being refined */
              int *glob2loc,                /* maps vertices to subgraph vertices */
              int *loc2glob,                /* maps subgraph vertices to vertices */
              int *sub_assign,              /* new assignment for subgraphs */
              int *old_sub_assign,          /* current assignment for subgraphs */
              int *degrees,                 /* space for forming subgraphs */
              int  using_ewgts,             /* are edge weights being used? */
              int (*hops)[MAXSETS],         /* KL set preferences */
              double *goal,                 /* desired set sizes */
              int *   sizes,                /* number of vertices in different sets */
              float * term_wgts[],          /* space for terminal propagation weights */
              int     architecture,         /* 0 => hypercube, d => d-dimensional mesh */
              int     mesh_dims[3]          /* if mesh, how big is it? */
)
{
  extern int     TERM_PROP;               /* perform terminal propagation? */
  extern double  KL_IMBALANCE;            /* fractional imbalance allowed in KL */
  struct bilist *ptr;                     /* element in set_list */
  double         subgoal[2];              /* goal within two subgraphs */
  double         weights[2] = {0.0, 0.0}; /* weights for each set */
  double         maxdeg;                  /* largest degree of a vertex */
  double         ratio;                   /* set sizes / goals */
  int *          null_ptr;                /* argument to klspiff */
  int            vwgt_max;                /* largest vertex weight */
  int            max_dev;                 /* largest set deviation allowed in KL */
  int            subnvtxs;                /* number of vtxs in subgraph */
  int            vwgt_sum1;               /* sum of vertex wgts in first set */
  int            vwgt_sum2;               /* sum of vertex wgts in second set */
  int            subnedges;               /* number of edges in subgraph */
  int            setA, setB;              /* two sets being refined */
  int            nsame;                   /* number of vertices not moved */
  int            vtx;                     /* vertex in subgraph */
  int            i;                       /* loop counter */
  double         find_maxdeg();
  void           make_maps_ref(), make_subgraph(), remake_graph();
  void           klspiff(), make_terms_ref(), count_weights();

  /* Compute all the quantities I'll need. */
  null_ptr = NULL;
  make_maps_ref(graph, set_list, vtx_elems, new_assign, sub_assign, set1, set2, glob2loc, loc2glob,
                &subnvtxs, &vwgt_max, &vwgt_sum1, &vwgt_sum2);

  for (i = 1; i <= subnvtxs; i++) {
    old_sub_assign[i] = sub_assign[i];
  }

  /* Set up goals for this KL invocation. */
  ratio      = (vwgt_sum1 + vwgt_sum2) / (goal[set1] + goal[set2]);
  subgoal[0] = ratio * goal[set1];
  subgoal[1] = ratio * goal[set2];

  if (TERM_PROP) {
    make_terms_ref(graph, using_ewgts, subnvtxs, loc2glob, set1, set2, new_assign, architecture,
                   mesh_dims, term_wgts);
  }

  /* New_assign has overwritten set2 with set1. */
  make_subgraph(graph, subgraph, subnvtxs, &subnedges, new_assign, set1, glob2loc, loc2glob,
                degrees, using_ewgts);

  maxdeg = find_maxdeg(subgraph, subnvtxs, using_ewgts, (float *)NULL);

  count_weights(subgraph, subnvtxs, sub_assign, 2, weights, (vwgt_max != 1));

  max_dev = vwgt_max;
  ratio   = (subgoal[0] + subgoal[1]) * KL_IMBALANCE / 2;
  if (ratio > max_dev) {
    max_dev = ratio;
  }

  klspiff(subgraph, subnvtxs, sub_assign, 2, hops, subgoal, term_wgts, max_dev, maxdeg, using_ewgts,
          &null_ptr, weights);

  /* Figure out which modification leaves most vertices intact. */
  nsame = 0;
  for (i = 1; i <= subnvtxs; i++) {
    if (old_sub_assign[i] == sub_assign[i]) {
      nsame++;
    }
  }
  if (2 * nsame > subnvtxs) {
    setA = set1;
    setB = set2;
  }
  else {
    setA = set2;
    setB = set1;
  }

  /* Now update the assignments. */
  sizes[setA] = sizes[setB] = 0;
  for (i = 1; i <= subnvtxs; i++) {
    vtx = loc2glob[i];
    /* Update the set_lists. */
    ptr = &(vtx_elems[vtx]);
    if (ptr->next != NULL) {
      ptr->next->prev = ptr->prev;
    }
    if (ptr->prev != NULL) {
      ptr->prev->next = ptr->next;
    }

    if (sub_assign[i] == 0) {
      new_assign[vtx] = setA;
      ++sizes[setA];
      ptr->next = set_list[setA].next;
      if (ptr->next != NULL) {
        ptr->next->prev = ptr;
      }
      ptr->prev           = &(set_list[setA]);
      set_list[setA].next = ptr;
    }
    else {
      new_assign[vtx] = setB;
      ++sizes[setB];
      ptr->next = set_list[setB].next;
      if (ptr->next != NULL) {
        ptr->next->prev = ptr;
      }
      ptr->prev           = &(set_list[setB]);
      set_list[setB].next = ptr;
    }
  }

  remake_graph(subgraph, subnvtxs, loc2glob, degrees, using_ewgts);

  return (nsame != subnvtxs);
}
