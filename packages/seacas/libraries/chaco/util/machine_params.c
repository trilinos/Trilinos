/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __STDC__
#include <float.h>
#endif

/* Returns some machine/compiler specific values. */
/* Note: These values might not be calculated precisely correctly.*/
/*       If you know them for your machine, replace this code! */

void machine_params(double *double_epsilon, double *double_max)
/* returns machine precision */
/* returns largest double value */
{

#ifndef DBL_EPSILON
  double eps; /* machine precision */
#endif

#ifndef DBL_MAX

#ifndef DBL_MIN
  double double_min;    /* smallest double value */
  double min, min_prev; /* values halved to compute double_min */
#endif

  double max; /* largest double precision value */
#endif

#ifndef DBL_EPSILON
  eps = 1.0 / 16.0;
  while (1.0 + eps > 1.0)
    eps /= 2.0;
  *double_epsilon = eps * 2.0;
#else
  *double_epsilon = DBL_EPSILON;
#endif

#ifndef DBL_MAX

#ifndef DBL_MIN
  min_prev = min = 1.0;
  while (min * 1.0 > 0) {
    min_prev = min;
    min /= 32.0;
  }
  min = min_prev;
  while (min * 1.0 > 0) {
    min_prev = min;
    min /= 2.0;
  }
  double_min = min_prev / (*double_epsilon);
#else
  double_min = DBL_MIN;
#endif

  max = 2.0 * (1.0 - *double_epsilon) / (double_min);
  /*
     two_max = max*2.0;
     while (two_max/2.0 == max) {
        max = two_max;
        two_max = max*2.0;
     }
  */
  *double_max = max;
#else
  *double_max     = DBL_MAX;
#endif
}
