/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

void makev2cv(
    /* Construct mapping from original graph vtxs to coarsened graph vtxs. */
    int *mflag, /* flag indicating vtx selected or not */
    int  nvtxs, /* number of vtxs in original graph */
    int *v2cv   /* mapping from vtxs to coarsened vtxs */
)
{
  int i, j; /* loop counters */

  j = 1;
  for (i = 1; i <= nvtxs; i++) {
    if (mflag[i] == 0 || mflag[i] > i) {
      v2cv[i] = j++;
    }
    else {
      v2cv[i] = v2cv[mflag[i]];
    }
  }
}
