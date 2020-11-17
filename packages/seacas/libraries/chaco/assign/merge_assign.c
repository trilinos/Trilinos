/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

/* Combine the old assignment value with the new partition. */
void merge_assignments(int *assignment, /* assignment list for graph */
                       int *subassign,  /* subgraph assignment list */
                       int *subsets,    /* mapping from local to global sets */
                       int  subnvtxs,   /* number of vtxs in subgraph */
                       int *loc2glob    /* subgraph -> graph numbering map */
)
{
  int i; /* loop counter */

  for (i = 1; i <= subnvtxs; i++) {
    assignment[loc2glob[i]] = subsets[subassign[i]];
  }
}
