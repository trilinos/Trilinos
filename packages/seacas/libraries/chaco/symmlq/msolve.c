/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "structs.h"

int msolve(int nvtxs, double *x, double *y)
{
  int i;

  /* Just do a copy for now. */
  for (i = 0; i < nvtxs; i++) {
    y[i] = x[i];
  }

  return (0);
}
