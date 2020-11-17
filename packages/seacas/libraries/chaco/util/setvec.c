/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

/* Set a double precision vector to constant over range. */
void setvec(double *vec, int beg, int end, double setval)
{
  int i;

  vec = vec + beg;
  for (i = end - beg + 1; i; i--) {
    (*vec++) = setval;
  }
}

/* Set a float precision vector to constant over range. */
void setvec_float(float *vec, int beg, int end, float setval)
{
  int i;

  vec = vec + beg;
  for (i = end - beg + 1; i; i--) {
    (*vec++) = setval;
  }
}
