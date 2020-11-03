/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

/* Normalizes a double n-vector over range. */
double ch_normalize(double *vec, int beg, int end)
{
  int    i;
  double scale;
  double ch_norm(double *vec, int beg, int end);

  scale = ch_norm(vec, beg, end);
  vec   = vec + beg;
  for (i = end - beg + 1; i; i--) {
    *vec = *vec / scale;
    vec++;
  }
  return (scale);
}

/* Normalizes such that element k is positive */
double sign_normalize(double *vec, int beg, int end, int k)
{
  int    i;
  double scale, scale2;
  double ch_norm(double *vec, int beg, int end);

  scale = ch_norm(vec, beg, end);
  if (vec[k] < 0) {
    scale2 = -scale;
  }
  else {
    scale2 = scale;
  }
  vec = vec + beg;
  for (i = end - beg + 1; i; i--) {
    *vec = *vec / scale2;
    vec++;
  }
  return (scale);
}

/* Normalizes a float n-vector over range. */
double normalize_float(float *vec, int beg, int end)
{
  int    i;
  float  scale;
  double norm_float(float *vec, int beg, int end);

  scale = norm_float(vec, beg, end);
  vec   = vec + beg;
  for (i = end - beg + 1; i; i--) {
    *vec = *vec / scale;
    vec++;
  }
  return ((double)scale);
}
