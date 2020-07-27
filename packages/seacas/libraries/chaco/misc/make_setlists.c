/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

/* Make linked list of elements in each set. */
void make_setlists(int *setlists,  /* space for linked lists of vertices in sets */
                   int *list_ptrs, /* headers of each linked list */
                   int  nsets,     /* number of sets being created this step */
                   int *subsets,   /* subsets being created at current step */
                   int *subassign, /* assignment for vertices in (sub)graph */
                   int *loc2glob,  /* mapping from subgraph to graph numbering */
                   int  subnvtxs,  /* number of vertices in (sub)graph */
                   int  first      /* is this full graph or subgraph? */
)
{
  int set;  /* set a vertex belongs to */
  int i, j; /* loop counters */

  if (first) { /* Don't need to indirect set numbers. */
    for (i = 0; i < nsets; i++) {
      list_ptrs[subsets[i]] = 0;
    }

    for (i = subnvtxs; i; i--) {
      set            = subassign[i];
      setlists[i]    = list_ptrs[set];
      list_ptrs[set] = i;
    }
  }

  else {
    for (i = 0; i < nsets; i++) {
      list_ptrs[subsets[i]] = 0;
    }

    for (i = subnvtxs; i; i--) {
      set            = subsets[subassign[i]];
      j              = loc2glob[i];
      setlists[j]    = list_ptrs[set];
      list_ptrs[set] = j;
    }
  }
}
