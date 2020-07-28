/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <stdio.h> // for NULL

void y2x(double **xvecs,   /* pointer to list of x-vectors */
         int      ndims,   /* number of divisions to make (# xvecs) */
         int      nmyvtxs, /* number of vertices I own (length xvecs) */
         double * wsqrt    /* sqrt of vertex weights */
)

/* Convert from y to x by dividing by wsqrt. */
{
  double *wptr; /* loops through wsqrt */
  double *xptr; /* loops through elements of a xvec */
  int     i, j; /* loop counters */

  if (wsqrt == NULL) {
    return;
  }

  for (i = 1; i <= ndims; i++) {
    xptr = xvecs[i];
    wptr = wsqrt;
    for (j = nmyvtxs; j; j--) {
      *(++xptr) /= *(++wptr);
    }
  }
}

void x2y(double **yvecs,   /* pointer to list of y-vectors */
         int      ndims,   /* number of divisions to make (# yvecs) */
         int      nmyvtxs, /* number of vertices I own (length yvecs) */
         double * wsqrt    /* sqrt of vertex weights */
)

/* Convert from x to y by multiplying by wsqrt. */
{
  double *wptr; /* loops through wsqrt */
  double *yptr; /* loops through elements of a yvec */
  int     i, j; /* loop counters */

  if (wsqrt == NULL) {
    return;
  }

  for (i = 1; i <= ndims; i++) {
    yptr = yvecs[i];
    wptr = wsqrt;
    for (j = nmyvtxs; j; j--) {
      *(++yptr) *= *(++wptr);
    }
  }
}
