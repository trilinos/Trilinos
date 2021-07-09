/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "defs.h"
#include "internal.h"
#include "smalloc.h"
#include "structs.h"
#include <stdio.h>

/* Greedily increase the number of internal vtxs in each set. */
void force_internal(struct vtx_data **graph,       /* graph data structure */
                    int               nvtxs,       /* number of vertices in graph */
                    int               using_ewgts, /* are edge weights being used? */
                    int *             assign,      /* current assignment */
                    double *          goal,        /* desired set sizes */
                    int               nsets_tot,   /* total number of sets */
                    int               npasses_max  /* number of passes to make */
)
{
  extern int     DEBUG_TRACE;          /* trace main execution path? */
  extern int     DEBUG_INTERNAL;       /* turn on debugging code here? */
  struct bidint *prev;                 /* back pointer for setting up lists */
  struct bidint *int_list      = NULL; /* internal vwgt in each set */
  struct bidint *vtx_elems     = NULL; /* linked lists of vtxs in each set */
  struct bidint *set_list      = NULL; /* headers for vtx_elems lists */
  double *       internal_vwgt = NULL; /* total internal vwgt in each set */
  int *          total_vwgt    = NULL; /* total vertex weight in each set */
  int *          indices       = NULL; /* orders sets by internal vwgt */
  int *          locked        = NULL; /* is vertex allowed to switch sets? */
  int            internal;             /* is a vertex internal or not? */
  int *          space = NULL;         /* space for mergesort */
  int            npasses;              /* number of callse to improve_internal */
  int            nlocked;              /* number of vertices that can't move */
  int            set, set2;            /* sets two vertices belong to */
  int            any_change;           /* did pass improve # internal vtxs? */
  int            niter;                /* counts calls to improve_internal */
  int            vwgt_max;             /* largest vertex weight in graph */
  int            progress;             /* am I improving # internal vertices? */
  int            error;                /* out of space? */
  int            size;                 /* array spacing */
  int            i, j;                 /* loop counters */
  int            improve_internal();
  void           ch_mergesort(double *vals, int nvals, int *indices, int *space);
  void           check_internal(), strout();

  error = 1;

  /* For each set, compute the total weight of internal vertices. */

  if (DEBUG_TRACE > 0) {
    printf("<Entering force_internal>\n");
  }

  indices       = smalloc_ret(nsets_tot * sizeof(int));
  internal_vwgt = smalloc_ret(nsets_tot * sizeof(double));
  total_vwgt    = smalloc_ret(nsets_tot * sizeof(int));
  if (indices == NULL || internal_vwgt == NULL || total_vwgt == NULL) {
    goto skip;
  }

  for (set = 0; set < nsets_tot; set++) {
    total_vwgt[set] = internal_vwgt[set] = 0;
    indices[set]                         = set;
  }

  vwgt_max = 0;
  for (i = 1; i <= nvtxs; i++) {
    internal = TRUE;
    set      = assign[i];
    for (j = 1; j < graph[i]->nedges && internal; j++) {
      set2     = assign[graph[i]->edges[j]];
      internal = (set2 == set);
    }

    total_vwgt[set] += graph[i]->vwgt;
    if (internal) {
      internal_vwgt[set] += graph[i]->vwgt;
    }
    if (graph[i]->vwgt > vwgt_max) {
      vwgt_max = graph[i]->vwgt;
    }
  }

  /* Now sort all the internal_vwgt values. */
  space = smalloc_ret(nsets_tot * sizeof(int));
  if (space == NULL) {
    goto skip;
  }
  ch_mergesort(internal_vwgt, nsets_tot, indices, space);
  sfree(space);
  space = NULL;

  /* Now construct a doubly linked list of sorted, internal_vwgt values. */
  int_list = smalloc_ret((nsets_tot + 1) * sizeof(struct bidint));
  if (int_list == NULL) {
    goto skip;
  }

  prev       = &(int_list[nsets_tot]);
  prev->prev = NULL;
  for (i = 0; i < nsets_tot; i++) {
    set                = indices[i];
    int_list[set].prev = prev;
    int_list[set].val  = internal_vwgt[set];
    prev->next         = &(int_list[set]);
    prev               = &(int_list[set]);
  }
  prev->next              = NULL;
  int_list[nsets_tot].val = -1;

  sfree(internal_vwgt);
  sfree(indices);
  internal_vwgt = NULL;
  indices       = NULL;

  /* Set up convenient data structure for navigating through sets. */
  set_list  = smalloc_ret(nsets_tot * sizeof(struct bidint));
  vtx_elems = smalloc_ret((nvtxs + 1) * sizeof(struct bidint));
  if (set_list == NULL || vtx_elems == NULL) {
    goto skip;
  }

  for (i = 0; i < nsets_tot; i++) {
    set_list[i].next = NULL;
  }

  for (i = 1; i <= nvtxs; i++) {
    set               = assign[i];
    vtx_elems[i].next = set_list[set].next;
    if (vtx_elems[i].next != NULL) {
      vtx_elems[i].next->prev = &(vtx_elems[i]);
    }
    vtx_elems[i].prev  = &(set_list[set]);
    set_list[set].next = &(vtx_elems[i]);
  }

  locked = smalloc_ret((nvtxs + 1) * sizeof(int));
  if (locked == NULL) {
    goto skip;
  }

  nlocked = 0;
  size    = (int)(&(int_list[1]) - &(int_list[0]));

  any_change = TRUE;
  npasses    = 1;
  while (any_change && npasses <= npasses_max) {
    for (i = 1; i <= nvtxs; i++) {
      locked[i] = FALSE;
    }

    /* Now select top guy off the list and improve him. */
    any_change = FALSE;
    progress   = TRUE;
    niter      = 1;
    while (progress) {
      prev = int_list[nsets_tot].next;
      set  = ((int)(prev - int_list)) / size;

      if (DEBUG_INTERNAL > 0) {
        printf("Before iteration %d, nlocked = %d, int[%d] = %d\n", niter, nlocked, set, prev->val);
      }
      if (DEBUG_INTERNAL > 1) {
        check_internal(graph, nvtxs, int_list, set_list, vtx_elems, total_vwgt, assign, nsets_tot);
      }

      progress = improve_internal(graph, nvtxs, assign, goal, int_list, set_list, vtx_elems, set,
                                  locked, &nlocked, using_ewgts, vwgt_max, total_vwgt);
      if (progress) {
        any_change = TRUE;
      }
      niter++;
    }
    npasses++;
  }
  error = 0;

skip:

  if (error) {
    strout("\nWARNING: No space to increase internal vertices.");
    strout("         NO INTERNAL VERTEX INCREASE PERFORMED.\n");
  }

  sfree(internal_vwgt);
  sfree(indices);
  sfree(locked);
  sfree(total_vwgt);
  sfree(vtx_elems);
  sfree(int_list);
  sfree(set_list);
}
