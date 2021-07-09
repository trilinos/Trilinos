/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "defs.h"
#include <math.h>
#include <stdio.h>

/* Print vertically range of double vector. */
void vecout(double *vec, int beg, int end, char *tag, char *file_name)
{
  FILE *file;
  int   i;

  if (file_name != NULL) {
    file = fopen(file_name, "w");
  }
  else {
    file = stdout;
  }

  fprintf(file, "%s:\n", tag);
  for (i = beg; i <= end; i++) {
    if (fabs(vec[i]) >= 1.0e-16) {
      fprintf(file, "%2d.   %24.16f\n", i, vec[i]);
    }
    else {
      fprintf(file, "%2d.         %g \n", i, vec[i]);
    }
  }
  if (file_name != NULL) {
    fclose(file);
  }
}

/* Scale the eigenvector such that the first element is non-negative.
 * This is an attempt to reduce machine-to-machine variance which can
 * result from calculated eigenvectors being mirror-images of each
 * other due to small roundoff..
 */
void vecnorm(double *vec, int beg, int end)
{
  int i;
  if (vec[beg] < 0.0) {
    for (i = beg; i <= end; i++) {
      vec[i] *= -1.0;
    }
  }
}
