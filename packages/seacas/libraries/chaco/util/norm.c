/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <math.h>

/* Returns 2-norm of a double n-vector over range. */
double ch_norm(double *vec, int beg, int end)
{
  double dot(double *vec1, int beg, int end, double *vec2);

  return (sqrt(dot(vec, beg, end, vec)));
}

/* Returns 2-norm of a float n-vector over range. */
double norm_float(float *vec, int beg, int end)
{
  double temp;
  double dot_float(float *vec1, int beg, int end, float *vec2);

  temp = sqrt(dot_float(vec, beg, end, vec));
  return (temp);
}
