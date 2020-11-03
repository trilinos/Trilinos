/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

/* update - fills double vec1 with vec2 + alpha*vec3 over range*/
void update(double *vec1, int beg, int end, double *vec2, double fac, double *vec3)
{
  int i;

  vec1 += beg;
  vec2 += beg;
  vec3 += beg;
  for (i = end - beg + 1; i; i--) {
    (*vec1++) = (*vec2++) + fac * (*vec3++);
  }
}

/* update - fills float vec1 with vec2 + alpha*vec3 over range*/
void update_float(float *vec1, int beg, int end, float *vec2, float fac, float *vec3)
{
  int i;

  vec1 += beg;
  vec2 += beg;
  vec3 += beg;
  for (i = end - beg + 1; i; i--) {
    (*vec1++) = (*vec2++) + fac * (*vec3++);
  }
}
