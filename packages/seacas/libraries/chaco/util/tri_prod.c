/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <stdio.h>

double tri_prod(double *v1, double *v2, double *v3, double *wsqrt, int n)

/* Form inner product of three vectors. */
{
  double dot = 0;
  int    i;

  if (wsqrt == NULL) { /* Unweighted case.  Use weights = 1. */
    for (i = 1; i <= n; i++) {
      dot += v1[i] * v2[i] * v3[i];
    }
  }
  else { /* Unweighted case. */
    for (i = 1; i <= n; i++) {
      dot += v1[i] * v2[i] * v3[i] / wsqrt[i];
    }
  }
  return (dot);
}
