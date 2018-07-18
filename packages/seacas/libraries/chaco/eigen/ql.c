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

/* Eigensolution of real symmetric tridiagonal matrix using the algorithm
   of Numerical Recipies p. 380. Removed eigenvector calculation and added
   return codes: 1 if maximum number of iterations is exceeded, 0 otherwise.
   NOTE CAREFULLY: the vector e is used as workspace, the eigenvals are
   returned in the vector d. */

#include <math.h>

#define SIGN(a, b) ((b) < 0 ? -fabs(a) : fabs(a))

int ql(double d[], double e[], int n)

{
  int    m, l, iter, i;
  double s, r, p, g, f, dd, c, b;

  e[n] = 0.0;

  for (l = 1; l <= n; l++) {
    iter = 0;
    do {
      for (m = l; m <= n - 1; m++) {
        dd = fabs(d[m]) + fabs(d[m + 1]);
        if (fabs(e[m]) + dd == dd) {
          break;
        }
      }
      if (m != l) {
        if (iter++ == 50) {
          return (1);
          /* ... not converging; bail out with error code. */
        }
        g = (d[l + 1] - d[l]) / (2.0 * e[l]);
        r = sqrt((g * g) + 1.0);
        g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
        s = c = 1.0;
        p     = 0.0;
        for (i = m - 1; i >= l; i--) {
          f = s * e[i];
          b = c * e[i];
          if (fabs(f) >= fabs(g)) {
            c        = g / f;
            r        = sqrt((c * c) + 1.0);
            e[i + 1] = f * r;
            c *= (s = 1.0 / r);
          }
          else {
            s        = f / g;
            r        = sqrt((s * s) + 1.0);
            e[i + 1] = g * r;
            s *= (c = 1.0 / r);
          }
          g        = d[i + 1] - p;
          r        = (d[i] - g) * s + 2.0 * c * b;
          p        = s * r;
          d[i + 1] = g + p;
          g        = c * r - b;
        }
        d[l] = d[l] - p;
        e[l] = g;
        e[m] = 0.0;
      }
    } while (m != l);
  }
  return (0); /* ... things seem ok */
}
