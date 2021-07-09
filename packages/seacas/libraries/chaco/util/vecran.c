/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

/* Fill double vector with random numbers over a range. */
void vecran(double *vec, int beg, int end)
{
  int     i;
  double *pntr;
  double  ch_normalize(double *vec, int beg, int end);
  double  drandom(void);

  pntr = vec + beg;
  for (i = end - beg + 1; i; i--) {
    (*pntr++) = drandom();
  }
  ch_normalize(vec, beg, end);
}

/* Fill float vector with random numbers over a range. */
void vecran_float(float *vec, int beg, int end)
{
  int    i;
  float *pntr;
  double normalize_float(float *vec, int beg, int end);
  double drandom(void);

  pntr = vec + beg;
  for (i = end - beg + 1; i; i--) {
    (*pntr++) = drandom();
  }
  normalize_float(vec, beg, end);
}
