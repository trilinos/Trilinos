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

#include <math.h>

/* Evaluates principal minor polynomial, returns the index of the eigenvalue just
   left of mu. Based on Wilkinson's algorithm AEP, p.302. Found that this algorithm
   could fail in practice when the p recursion overflowed. Fixed by rescaling the
   recursion if it got too large. Routine returns -1 if this re-scaling doesn't work
   for some reason. Otherwise it returns 0.  Note Wilkinson indexes beta such that
   first off-diagonal entry in T is beta[2]. We index such that it is beta[1]. */

int sturmcnt(double *alpha, /* vector of Lanczos scalars */
             double *beta,  /* vector of Lanczos scalars */
             int     j,     /* index of Lanczos step we're on */
             double  mu,    /* argument to the sequence generating polynomial */
             double *p      /* work vector for sturm sequence */
)
{
  extern double DOUBLE_MAX; /* maximum double precision number to be used */
  int           i;          /* loop index */
  int           cnt;        /* number of sign changes in the sequence */
  int           sign;       /* algebraic sign of current sequence value */
  int           last_sign;  /* algebraic sign of previous sequence value */
  double *      pntr_p;     /* for stepping through sequence array */
  double *      pntr_p1;    /* for stepping through sequence array one behind */
  double *      pntr_p2;    /* for stepping through sequence array two behind */
  double *      pntr_alpha; /* for stepping through alpha array */
  double *      pntr_beta;  /* for stepping through beta array */
  double        limit;      /* the cut-off for re-scaling the recurrence */
  double        scale;      /* magnitude of entry in p recursion */

  if (j == 1) {
    /* have to special case this one */
    if (alpha[1] > mu) {
      cnt = 1;
    }
    else {
      cnt = 0;
    }
  }
  else {
    /* compute the Sturm sequence */
    limit      = sqrt(DOUBLE_MAX);
    p[0]       = 1;
    p[1]       = alpha[1] - mu;
    pntr_p     = &p[2];
    pntr_p1    = &p[1];
    pntr_p2    = &p[0];
    pntr_alpha = &alpha[2];
    pntr_beta  = &beta[1];
    for (i = 2; i <= j; i++) {
      *pntr_p++ = (*pntr_alpha - mu) * (*pntr_p1) - (*pntr_beta) * (*pntr_beta) * (*pntr_p2);
      pntr_alpha++;
      pntr_beta++;
      pntr_p1++;
      pntr_p2++;
      scale = fabs(*pntr_p1);
      if (scale > limit) {
        *pntr_p1 /= scale;
        *pntr_p2 /= scale;
        /* re-scale to avoid overflow in p recursion */
      }
    }

    /* count sign changes */
    cnt       = 0;
    last_sign = 1;
    pntr_p    = &p[1];
    for (i = j; i; i--) {
      if (*pntr_p != *pntr_p || fabs(*pntr_p) > limit) {
        return (-1);
        /* ... re-scaling failed; bail out and return error code.
               Note (x != x) is TRUE iff x is NaN, so this check
               should be ok on non-IEEE floating point systems too. */
      }
      if (*pntr_p > 0) {
        sign = 1;
      }
      else if (*pntr_p < 0) {
        sign = -1;
      }
      else {
        sign = -last_sign;
      }
      if (sign == last_sign) {
        cnt++;
      }
      last_sign = sign;
      pntr_p++;
    }
  }

  /* cnt is number of evals strictly > than mu. Instead send back index of the
     eigenvalue just left of mu. */
  return (j - cnt); /* ... things seem ok. */
}
