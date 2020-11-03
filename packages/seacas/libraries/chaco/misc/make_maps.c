/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <stdio.h>

int make_maps(int *setlists,  /* linked list of set vertices */
              int *list_ptrs, /* head of each linked list */
              int  set,       /* set value denoting subgraph */
              int *glob2loc,  /* graph -> subgraph numbering map */
              int *loc2glob   /* subgraph -> graph numbering map */
)
{
  int i, j; /* loop counter */

  j = 0;
  i = list_ptrs[set];

  if (glob2loc != NULL) {
    while (i != 0) {
      loc2glob[++j] = i;
      glob2loc[i]   = j;
      i             = setlists[i];
    }
  }

  else {
    while (i != 0) {
      loc2glob[++j] = i;
      i             = setlists[i];
    }
  }

  return (j);
}

void make_maps2(int *assignment, /* set assignments for graph */
                int  nvtxs,      /* length of assignment */
                int  set,        /* set value denoting subgraph */
                int *glob2loc,   /* graph -> subgraph numbering map */
                int *loc2glob    /* subgraph -> graph numbering map */
)
{
  int i, j; /* loop counter */

  j = 0;
  if (glob2loc != NULL) {
    for (i = 1; i <= nvtxs; i++) {
      if (assignment[i] == set) {
        j++;
        glob2loc[i] = j;
        loc2glob[j] = i;
      }
    }
  }
  else {
    for (i = 1; i <= nvtxs; i++) {
      if (assignment[i] == set) {
        j++;
        loc2glob[j] = i;
      }
    }
  }
}
