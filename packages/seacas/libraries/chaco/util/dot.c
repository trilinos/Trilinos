/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

/* Returns scalar product of two double n-vectors. */
double dot(double *vec1, int beg, int end, double *vec2)
{
  int    i;
  double sum;

  sum  = 0.0;
  vec1 = vec1 + beg;
  vec2 = vec2 + beg;
  for (i = end - beg + 1; i; i--) {
    sum += (*vec1++) * (*vec2++);
  }
  return (sum);
}

/* Returns scalar product of two float n-vectors. */
double dot_float(float *vec1, int beg, int end, float *vec2)
{
  int   i;
  float sum;

  sum  = 0.0;
  vec1 = vec1 + beg;
  vec2 = vec2 + beg;
  for (i = end - beg + 1; i; i--) {
    sum += (*vec1++) * (*vec2++);
  }
  return ((double)sum);
}
