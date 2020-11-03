/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

/* Scaled add - fills double vec1 with vec1 + alpha*vec2 over range*/
void scadd(double *vec1, int beg, int end, double fac, double *vec2)
{
  int i;

  vec1 = vec1 + beg;
  vec2 = vec2 + beg;
  for (i = end - beg + 1; i; i--) {
    (*vec1++) += fac * (*vec2++);
  }
}

/* Scaled add - fills float vec1 with vec1 + alpha*vec2 over range*/
void scadd_float(float *vec1, int beg, int end, float fac, float *vec2)
{
  int i;

  vec1 = vec1 + beg;
  vec2 = vec2 + beg;
  for (i = end - beg + 1; i; i--) {
    (*vec1++) += fac * (*vec2++);
  }
}

/* Scaled add - fills double vec1 with vec1 + alpha*vec2 where vec2 is float */
void scadd_mixed(double *vec1, int beg, int end, double fac, float *vec2)
{
  int i;

  vec1 = vec1 + beg;
  vec2 = vec2 + beg;
  for (i = end - beg + 1; i; i--) {
    (*vec1++) += fac * (*vec2++);
  }
}
