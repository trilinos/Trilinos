/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "params.h"  // for MAXSETS
#include "smalloc.h" // for smalloc

void genvals2d(
    /* Create lists of sets of values to be sorted. */
    double **xvecs,            /* vectors to partition */
    double  *vals[4][MAXSETS], /* ptrs to lists of values */
    int      nvtxs             /* number of values */
)
{
  int     nlists = 4; /* number of lists to generate */
  double *temp[4];    /* place holders for vals */
  int     i;          /* loop counter */

  for (i = 0; i < nlists; i++) {
    temp[i] = smalloc(nvtxs * sizeof(double));
  }

  for (i = 1; i <= nvtxs; i++) {
    temp[0][i - 1] = 4 * xvecs[1][i];
    temp[1][i - 1] = 4 * xvecs[2][i];
    temp[2][i - 1] = 4 * (xvecs[1][i] + xvecs[2][i]);
    temp[3][i - 1] = 4 * (xvecs[2][i] - xvecs[1][i]);
  }

  vals[0][1] = vals[1][0] = vals[2][3] = vals[3][2] = temp[0];
  vals[0][2] = vals[2][0] = vals[1][3] = vals[3][1] = temp[1];
  vals[0][3] = vals[3][0] = temp[2];
  vals[1][2] = vals[2][1] = temp[3];
}
