/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "smalloc.h"
#include <stdio.h>

/* Allocates a double vector with range [nl..nh]. Dies. */
double *mkvec(int nl, int nh)
{
  double *v;

  v = smalloc((nh - nl + 1) * sizeof(double));
  return (v - nl);
}

/* Allocates a double vector with range [nl..nh]. Returns error code. */
double *mkvec_ret(int nl, int nh)
{
  double *v;

  v = smalloc_ret((nh - nl + 1) * sizeof(double));
  if (v == NULL) {
    return (NULL);
  }

  return (v - nl);
}

/* Free a double vector with range [nl..nh]. */
void frvec(double *v, int nl) { sfree((v + nl)); }

/* Allocates a float vector with range [nl..nh]. Dies. */
float *mkvec_float(int nl, int nh)
{
  float *v;

  v = smalloc((nh - nl + 1) * sizeof(float));
  return (v - nl);
}

/* Allocates a float vector with range [nl..nh]. Returns error code. */
float *mkvec_ret_float(int nl, int nh)
{
  float *v;

  v = smalloc_ret((nh - nl + 1) * sizeof(float));
  if (v == NULL) {
    return (NULL);
  }

  return (v - nl);
}

/* Free a float vector with range [nl..nh]. */
void frvec_float(float *v, int nl) { sfree((v + nl)); }
