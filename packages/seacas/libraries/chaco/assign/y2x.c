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

#include <stdio.h> // for NULL

void y2x(double **xvecs,   /* pointer to list of x-vectors */
         int      ndims,   /* number of divisions to make (# xvecs) */
         int      nmyvtxs, /* number of vertices I own (lenght xvecs) */
         double * wsqrt    /* sqrt of vertex weights */
)

/* Convert from y to x by dividing by wsqrt. */
{
  double *wptr; /* loops through wsqrt */
  double *xptr; /* loops through elements of a xvec */
  int     i, j; /* loop counters */

  if (wsqrt == NULL) {
    return;
  }

  for (i = 1; i <= ndims; i++) {
    xptr = xvecs[i];
    wptr = wsqrt;
    for (j = nmyvtxs; j; j--) {
      *(++xptr) /= *(++wptr);
    }
  }
}

void x2y(double **yvecs,   /* pointer to list of y-vectors */
         int      ndims,   /* number of divisions to make (# yvecs) */
         int      nmyvtxs, /* number of vertices I own (lenght yvecs) */
         double * wsqrt    /* sqrt of vertex weights */
)

/* Convert from x to y by multiplying by wsqrt. */
{
  double *wptr; /* loops through wsqrt */
  double *yptr; /* loops through elements of a yvec */
  int     i, j; /* loop counters */

  if (wsqrt == NULL) {
    return;
  }

  for (i = 1; i <= ndims; i++) {
    yptr = yvecs[i];
    wptr = wsqrt;
    for (j = nmyvtxs; j; j--) {
      *(++yptr) *= *(++wptr);
    }
  }
}
