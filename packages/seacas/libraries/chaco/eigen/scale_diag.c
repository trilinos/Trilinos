/*
 * Copyright(C) 1999-2021 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <stdio.h>

/* Scales vector by diagonal matrix (passed as vector) over range. */
void scale_diag(double *vec,      /* the vector to scale */
                int beg, int end, /* specify the range to norm over */
                double *diag      /* vector to scale by */
)
{
  if (diag != NULL) {
    vec  = vec + beg;
    diag = diag + beg;
    for (int i = end - beg + 1; i; i--) {
      *vec++ *= *diag++;
    }
  }
  /* otherwise return vec unchanged */
}

/* Scales vector by diagonal matrix (passed as vector) over range. */
void scale_diag_float(float *vec,       /* the vector to scale */
                      int beg, int end, /* specify the range to norm over */
                      float *diag       /* vector to scale by */
)
{
  if (diag != NULL) {
    vec  = vec + beg;
    diag = diag + beg;
    for (int i = end - beg + 1; i; i--) {
      *vec++ *= *diag++;
    }
  }
}
