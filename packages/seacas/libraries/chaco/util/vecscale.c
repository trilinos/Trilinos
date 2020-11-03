/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

/* Scale - fills vec1 with alpha*vec2 over range, double version */
void vecscale(double *vec1, int beg, int end, double alpha, double *vec2)
{
  int i;

  vec1 += beg;
  vec2 += beg;
  for (i = end - beg + 1; i; i--) {
    (*vec1++) = alpha * (*vec2++);
  }
}

/* Scale - fills vec1 with alpha*vec2 over range, float version */
void vecscale_float(float *vec1, int beg, int end, float alpha, float *vec2)
{
  int i;

  vec1 += beg;
  vec2 += beg;
  for (i = end - beg + 1; i; i--) {
    (*vec1++) = alpha * (*vec2++);
  }
}
