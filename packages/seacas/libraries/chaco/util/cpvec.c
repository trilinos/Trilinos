/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

/* Copy a range of a double vector to a double vector */
void cpvec(double *copy, int beg, int end, double *vec)
{
  int i;

  copy = copy + beg;
  vec  = vec + beg;
  for (i = end - beg + 1; i; i--) {
    *copy++ = *vec++;
  }
}

/* Copy a range of a float vector to a double vector */
void float_to_double(double *copy, int beg, int end, float *vec)
{
  int i;

  copy = copy + beg;
  vec  = vec + beg;
  for (i = end - beg + 1; i; i--) {
    *copy++ = *vec++;
  }
}

/* Copy a range of a double vector to a float vector */
void double_to_float(float *copy, int beg, int end, double *vec)
{
  int i;

  copy = copy + beg;
  vec  = vec + beg;
  for (i = end - beg + 1; i; i--) {
    *copy++ = *vec++;
  }
}
