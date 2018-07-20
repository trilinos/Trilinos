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

#include <math.h>  // for fabs, sqrt
#include <stdio.h> // for printf

/* Solve the shifted, symmetric tridiagonal linear system (T - lambda*I)v = b
   where b is a multiple of e1 (the first column of the identity matrix).
   T is constructed with diagonal alpha[1], alpha[2] ... alpha[j] and
   off diagonal beta[1], beta[2] ... beta[j-1]. The algorithm is based on
   p. 156 of Golub and Van Loan, 2nd ed. Modified this to (1) incorporate
   shift by lambda, (2) use work vectors rather than overwriting so
   won't have to copy input as this is called repeatedly in the bisection
   loop for the extended eigenproblem, and (3) to squawk and die if a
   pivot below machine precision is encountered and (4) exploit the fact
   that all elements of the rhs except the first are 0. */
void tri_solve(double *alpha,  /* vector of Lanczos scalars */
               double *beta,   /* vector of Lanczos scalars */
               int     j,      /* system size */
               double  lambda, /* diagonal shift */
               double *v,      /* solution vector */
               double  b,      /* scalar multiple of e1 specifying the rhs */
               double *d,      /* work vec. for diagonal of Cholesky factor */
               double *e       /* work vec. for off diagonal of Cholesky factor */
               )
{
  extern int DEBUG_EVECS; /* debug flag for eigen computation */
  int        i;           /* loop index */
  double     tol;         /* cutoff for numerical singularity */
  double     diff;        /* an element of the residual vector */
  double     res;         /* the norm of the residual vector */
  void       bail();      /* our exit routine */

  /* set cut-off for "zero" pivot */
  tol = 1.0e-15;

  /* Cholesky factorization of shifted symmetric tridiagonal */
  d[1] = alpha[1] - lambda;
  if (fabs(d[1]) < tol) {
    bail("ERROR: Zero pivot in tri_solve().", 1);
  }
  if (j == 1) {
    v[1] = b / d[1];
    return;
    /* It's a scalar equation, so we're done. */
  }
  for (i = 2; i <= j; i++) {
    e[i - 1] = beta[i - 1] / d[i - 1];
    d[i]     = alpha[i] - lambda - d[i - 1] * e[i - 1] * e[i - 1];
    if (fabs(d[i]) < tol) {
      bail("ERROR: Zero pivot in tri_solve().", 1);
    }
  }

  /* Substitution solution using Cholesky factors */
  v[1] = b;
  for (i = 2; i <= j; i++) {
    v[i] = -e[i - 1] * v[i - 1];
  }
  v[j] = v[j] / d[j];
  for (i = j - 1; i >= 1; i--) {
    v[i] = v[i] / d[i] - e[i] * v[i + 1];
  }

  /* Check residual */
  if (DEBUG_EVECS > 1) {
    diff = b - ((alpha[1] - lambda) * v[1] + beta[1] * v[2]);
    res  = diff * diff;
    for (i = 2; i <= j - 1; i++) {
      diff = beta[i - 1] * v[i - 1] + (alpha[i] - lambda) * v[i] + beta[i] * v[i + 1];
      res += diff * diff;
    }
    diff = beta[j - 1] * v[j - 1] + (alpha[j] - lambda) * v[j];
    res += diff * diff;
    res = sqrt(res);
    if (res > 1.0e-13) {
      printf("Tridiagonal solve residual %g\n", res);
    }
  }
}
