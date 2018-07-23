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

#include "defs.h" // for FALSE, TRUE
#include <math.h> // for fabs

/* NOTE: This should only be called if j >= 2. It returns a residual and
         the point where the forward and backward recurrences met. */

/* NOTE: The paper has an error in the residual calculation. And the
         heuristic of increasing the cut-off for the local maximum
         by a factor of 10 when the residual tolerance is not met
         appears to work poorly - letting the recursion run further
         to find a bigger local max often results in <less> accurate
         results. Static cut-off of 1 (recall we set the first element
         to 1) worked better. */

/* Finds eigenvector s of T and returns residual norm. */
double bidir(double *alpha, /* vector of Lanczos scalars */
             double *beta,  /* vector of Lanczos scalars */
             int     j,     /* number of Lanczos iterations taken */
             double  ritz,  /* approximate eigenvalue  of T */
             double *s,     /* approximate eigenvector of T */
             double  hurdle /* hurdle for local maximum in recurrence */
)

{
  int    i;        /* index */
  int    k;        /* meeting point for bidirect. recurrence */
  double sav_val;  /* value at meeting point */
  int    iterate;  /* keep iterating? */
  double scale;    /* scale factor for bidirectional recurrence */
  double residual; /* eigen residual */

  double sign_normalize(); /* normalizes such that given entry positive */

  s[2]    = -(alpha[1] - ritz) / beta[2];
  i       = 3;
  iterate = TRUE;
  while (i <= j && iterate) {
    /* follow the forward recurrence until a local maximum > 1.0 */
    s[i] = -((alpha[i - 1] - ritz) * s[i - 1] + beta[i - 1] * s[i - 2]) / beta[i];
    if (fabs(s[i - 1]) > hurdle && i > 2) {
      if (fabs(s[i]) < fabs(s[i - 1]) && fabs(s[i - 1]) > fabs(s[i - 2])) {
        iterate = FALSE;
        k       = i - 1;
        sav_val = s[k];
      }
    }
    i++;
  }
  if (!iterate) {
    /* we stopped at an intermediate point */
    s[j]     = 1.0;
    s[j - 1] = -(alpha[j] - ritz) / beta[j];
    for (i = j; i >= k + 2; i--) {
      s[i - 2] = -((alpha[i - 1] - ritz) * s[i - 1] + beta[i] * s[i]) / beta[i - 1];
    }
    scale = s[k] / sav_val;
    for (i = 1; i <= k - 1; i++) {
      s[i] = s[i] * scale;
    }
    /* paper gets this expression wrong */
    residual = beta[k] * s[k - 1] + (alpha[k] - ritz) * s[k] + beta[k + 1] * s[k + 1];
  }
  else {
    /* we went all the way forward */
    residual = (alpha[j] - ritz) * s[j] + beta[j] * s[j - 1];
  }
  residual = fabs(residual) / sign_normalize(s, 1, j, j);
  return (residual);
}
