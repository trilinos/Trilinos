/*
 * Copyright (c) 2005-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
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
  else {
    return (v - nl);
  }
}

/* Free a double vector with range [nl..nh]. */
void frvec(double *v, int nl)
{

  sfree((v + nl));
  v = NULL;
}

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
  else {
    return (v - nl);
  }
}

/* Free a float vector with range [nl..nh]. */
void frvec_float(float *v, int nl)
{

  sfree((v + nl));
  v = NULL;
}
