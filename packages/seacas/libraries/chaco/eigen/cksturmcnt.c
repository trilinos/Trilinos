/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

/* Check sturmcnt */
void cksturmcnt(double *vec, int beg, int end, double x1, double x2, int *x1ck, int *x2ck,
                int *numck)
{

  int i, count;

  count = 0;
  for (i = beg; i <= end; i++) {
    if (vec[i] > x1) {
      count += 1;
    }
  }
  *x1ck = end - count;

  count = 0;
  for (i = beg; i <= end; i++) {
    if (vec[i] > x2) {
      count += 1;
    }
  }
  *x2ck = end - count;

  count = 0;
  for (i = beg; i <= end; i++) {
    if (vec[i] > x1 && vec[i] < x2) {
      count += 1;
    }
  }
  *numck = count;
}
