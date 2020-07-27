/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "structs.h"
#include <stdio.h> // for NULL

int aprod(long *lnvtxs, double *x, double *y, double *dA, double *vwsqrt, double *work,
          double *dorthlist /* vectors to orthogonalize against */
)
{
  int               nvtxs; /* int copy of long_nvtxs */
  struct vtx_data **A;
  struct orthlink * orthlist; /* vectors to orthogonalize against */
  void              splarax(), orthog1(), orthogvec(), orthogonalize();

  nvtxs    = (int)*lnvtxs;
  A        = (struct vtx_data **)dA;
  orthlist = (struct orthlink *)dorthlist;

  /* The offset on x and y is because the arrays come originally from Fortran
     declarations which index from 1 */
  splarax(y - 1, A, nvtxs, x - 1, vwsqrt, work - 1);

  /* Now orthogonalize against lower eigenvectors. */
  if (vwsqrt == NULL) {
    orthog1(y - 1, 1, nvtxs);
  }
  else {
    orthogvec(y - 1, 1, nvtxs, vwsqrt);
  }
  orthogonalize(y - 1, nvtxs, orthlist);

  return (0);
}
