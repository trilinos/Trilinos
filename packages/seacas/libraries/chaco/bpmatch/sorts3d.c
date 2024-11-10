/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "params.h"  // for MAXSETS
#include "smalloc.h" // for smalloc, sfree

void sorts3d(
    /* Sort the lists needed to find the splitter. */
    double *vals[8][MAXSETS],    /* lists of values to sort */
    int    *indices[8][MAXSETS], /* indices of sorted lists */
    int     nvtxs                /* number of vertices */
)
{
  int *space;       /* space for mergesort routine */
  int  nsets  = 8;  /* number of sets */
  int  nlists = 13; /* number of directions to sort */
  int *temp[13];    /* place holders for indices */
  int  i, j;        /* loop counter */

  void ch_mergesort(double *vals, int nvals, int *indices, int *space);

  space = smalloc(nvtxs * sizeof(int));

  for (i = 0; i < nlists; i++) {
    temp[i] = smalloc(nvtxs * sizeof(int));
  }

  ch_mergesort(vals[0][1], nvtxs, temp[0], space);
  ch_mergesort(vals[0][2], nvtxs, temp[1], space);
  ch_mergesort(vals[0][4], nvtxs, temp[2], space);
  ch_mergesort(vals[0][3], nvtxs, temp[3], space);
  ch_mergesort(vals[1][2], nvtxs, temp[4], space);
  ch_mergesort(vals[0][5], nvtxs, temp[5], space);
  ch_mergesort(vals[1][4], nvtxs, temp[6], space);
  ch_mergesort(vals[0][6], nvtxs, temp[7], space);
  ch_mergesort(vals[2][4], nvtxs, temp[8], space);
  ch_mergesort(vals[0][7], nvtxs, temp[9], space);
  ch_mergesort(vals[1][6], nvtxs, temp[10], space);
  ch_mergesort(vals[2][5], nvtxs, temp[11], space);
  ch_mergesort(vals[3][4], nvtxs, temp[12], space);

  sfree(space);

  indices[0][1] = indices[2][3] = indices[4][5] = indices[6][7] = temp[0];
  indices[0][2] = indices[1][3] = indices[4][6] = indices[5][7] = temp[1];
  indices[0][4] = indices[1][5] = indices[2][6] = indices[3][7] = temp[2];
  indices[0][3] = indices[4][7] = temp[3];
  indices[1][2] = indices[5][6] = temp[4];
  indices[0][5] = indices[2][7] = temp[5];
  indices[1][4] = indices[3][6] = temp[6];
  indices[0][6] = indices[1][7] = temp[7];
  indices[2][4] = indices[3][5] = temp[8];
  indices[0][7]                 = temp[9];
  indices[1][6]                 = temp[10];
  indices[2][5]                 = temp[11];
  indices[3][4]                 = temp[12];

  for (i = 0; i < nsets; i++) {
    for (j = i + 1; j < nsets; j++) {
      indices[j][i] = indices[i][j];
    }
  }
}
