/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

/* Orthogonalize a double vector to all one's */
void orthog1(double *x, int beg, int end)
{
  int     i;
  double *pntr;
  double  sum;
  int     len;

  len  = end - beg + 1;
  sum  = 0.0;
  pntr = x + beg;
  for (i = len; i; i--) {
    sum += *pntr++;
  }
  sum /= len;
  pntr = x + beg;
  for (i = len; i; i--) {
    *pntr++ -= sum;
  }
}

/* Orthogonalize a float vector to all one's */
void orthog1_float(float *x, int beg, int end)
{
  int    i;
  float *pntr;
  float  sum;
  int    len;

  len  = end - beg + 1;
  sum  = 0.0;
  pntr = x + beg;
  for (i = len; i; i--) {
    sum += *pntr++;
  }
  sum /= len;
  pntr = x + beg;
  for (i = len; i; i--) {
    *pntr++ -= sum;
  }
}
