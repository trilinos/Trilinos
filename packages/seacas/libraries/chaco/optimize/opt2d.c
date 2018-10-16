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

#include "defs.h"
#include "structs.h"
#include <math.h>
#include <stdio.h>

double opt2d(
    /* Compute rotation angle to minimize distance to discrete points. */
    struct vtx_data **graph,  /* data structure with vertex weights */
    double **         yvecs,  /* eigenvectors */
    int               nvtxs,  /* total number of vertices */
    int               nmyvtxs /* number of vertices I own */
)
{
  extern int DEBUG_OPTIMIZE; /* debug flag for optimization */
  double *   aptr, *bptr;    /* loop through yvecs */
  double     coeffs[5];      /* various products of yvecs */
  double     a, b;           /* temporary values */
  double     func = 0.0;     /* value of function to be minimized */
  double     grad, hess;     /* first and 2nd derivatives of function */
  double     grad_min;       /* acceptably small gradient */
  double     theta;          /* angle being optimized */
  double     step;           /* change in angle */
  double     step_max;       /* maximum allowed step */
  double     step_min;       /* minimum step => convergence */
  double     hess_min;       /* value for hessian is < 0 */
  double     hfact;          /* scaling for min tolerated hessian */
  double     w;              /* vertex weight squared */
  double     pdtol;          /* allowed error in hessian pd-ness */
  int        pdflag;         /* is hessian positive semi-definite? */
  int        i;              /* loop counter */
  double     func2d();
  double     grad2d();
  double     hess2d();

  /* Set parameters. */
  step_max = PI / 8;
  step_min = 2.0e-5;
  grad_min = 1.0e-7;
  pdtol    = 1.0e-8;
  hfact    = 2;

  for (i = 0; i < 5; i++) {
    coeffs[i] = 0;
  }
  aptr = yvecs[1] + 1;
  bptr = yvecs[2] + 1;
  for (i = 1; i <= nmyvtxs; i++) {
    a = *aptr++;
    b = *bptr++;
    w = graph[i]->vwgt;
    if (w == 1) {
      coeffs[0] += a * a * a * a;
      coeffs[1] += a * a * a * b;
      coeffs[2] += a * a * b * b;
      coeffs[3] += a * b * b * b;
      coeffs[4] += b * b * b * b;
    }
    else {
      w = 1 / (w * w);
      coeffs[0] += a * a * a * a * w;
      coeffs[1] += a * a * a * b * w;
      coeffs[2] += a * a * b * b * w;
      coeffs[3] += a * b * b * b * w;
      coeffs[4] += b * b * b * b * w;
    }
  }
  /* Adjust for normalization of eigenvectors. */
  /* This should make tolerances independent of vector length */
  for (i = 0; i < 5; i++) {
    coeffs[i] *= nvtxs;
  }

  i      = 0;
  theta  = 0.0;
  step   = step_max;
  pdflag = FALSE;
  grad   = 0;
  while (fabs(step) >= step_min && (!pdflag || fabs(grad) > grad_min)) {
    func = func2d(coeffs, theta);
    grad = grad2d(coeffs, theta);
    hess = hess2d(coeffs);

    if (hess < -pdtol) {
      pdflag = FALSE;
    }
    else {
      pdflag = TRUE;
    }

    hess_min = hfact * fabs(grad) / step_max;
    if (hess < hess_min) {
      hess = hess_min;
    }

    if (fabs(grad) > fabs(hess * step_max)) {
      step = -step_max * sign(grad);
    }
    else {
      step = -grad / hess;
    }

    theta += step;
    if (fabs(step) < step_min && !pdflag) { /* Convergence to non-min. */
      step = step_min;
      theta += step;
    }
    i++;
  }

  if (DEBUG_OPTIMIZE > 0) {
    printf("After %d passes, func=%e, theta = %f\n", i, func, theta);
  }

  return (theta);
}
