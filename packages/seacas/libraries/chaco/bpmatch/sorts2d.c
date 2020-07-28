/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "params.h"  // for MAXSETS
#include "smalloc.h" // for smalloc, sfree

void sorts2d(
    /* Sort the lists needed to find the splitter. */
    double *vals[4][MAXSETS],    /* lists of values to sort */
    int *   indices[4][MAXSETS], /* indices of sorted lists */
    int     nvtxs                /* number of vertices */
)
{
  int *space;      /* space for mergesort routine */
  int *temp[4];    /* place holders for indices */
  int  nlists = 4; /* number of directions to sort */
  int  i;          /* loop counter */

  void ch_mergesort(double *vals, int nvals, int *indices, int *space);

  space = smalloc(nvtxs * sizeof(int));

  for (i = 0; i < nlists; i++) {
    temp[i] = smalloc(nvtxs * sizeof(int));
  }

  ch_mergesort(vals[0][1], nvtxs, temp[0], space);
  ch_mergesort(vals[0][2], nvtxs, temp[1], space);
  ch_mergesort(vals[0][3], nvtxs, temp[2], space);
  ch_mergesort(vals[1][2], nvtxs, temp[3], space);

  sfree(space);

  indices[0][1] = indices[1][0] = indices[2][3] = indices[3][2] = temp[0];
  indices[0][2] = indices[2][0] = indices[1][3] = indices[3][1] = temp[1];
  indices[0][3] = indices[3][0] = temp[2];
  indices[1][2] = indices[2][1] = temp[3];
}
