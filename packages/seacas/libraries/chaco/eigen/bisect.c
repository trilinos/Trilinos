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

#include <stdio.h> // for printf, NULL

/* Finds selected eigenvalues of T using Sturm sequence bisection. Based
    on Wilkinson's algorithm, AEP, p.302. Returns 1 if sturmcnt() fails and
    2 if it hasn't converged in max_steps of bisection. If neither of these
    errors is detected the return value is 0. */
int bisect(double *alpha,        /* vector of Lanczos scalars */
           double *beta,         /* vector of Lanczos scalars */
           int     j,            /* number of Lanczos iterations taken */
           double  Anorm,        /* Gershgorin estimate */
           double *workj,        /* work vector for Sturm sequence */
           double *ritz,         /* array holding evals */
           int     nevals_left,  /* number of evals on right to find */
           int     nevals_right, /* number of evals on left to find */
           double  tol,          /* tolerance on bracket width */
           double *ritz_sav,     /* space to copy ritzvals for debugging */
           int     max_steps     /* maximum number of bisection steps allowed */
)
{
  extern int    DEBUG_EVECS;  /* debug flag for eigen computation */
  extern double DOUBLE_MAX;   /* largest double value */
  int           index;        /* index of sturm polynomial */
  int           i;            /* loop index */
  double *      pntr;         /* pntr to double array */
  double        x1, x2;       /* the bracketing interval */
  int           x1cnt, x2cnt; /* Sturm counts at x1 and x2 */
  double        x;            /* the inserted point */
  int           xcnt;         /* the Sturm count at x */
  int           steps;        /* number of bisection steps for a Ritzval */
  int           tot_steps;    /* number of bisection steps for all Ritzvals */
  int           numbracketed; /* number of evals between x1 and x2 */
  int           x1ck;         /* debugging check on x1cnt */
  int           x2ck;         /* debugging check on x2cnt */
  int           numck;        /* debugging check on numbracketed */
  int           sturmcnt();   /* counts the Sturm sequence */
  double        diff;         /* debugging register */
  int           ii;           /* debugging loop counter */
  void          cksturmcnt(); /* checks the Sturm sequence count against ql */

  /* If space has been allocated for a copy of the ritz values, assume
     we are to check the Sturm sequence counts directly using ql(). */
  if (ritz_sav != NULL) {
    printf("\nAnorm %g j %d nevals_left %d\n", Anorm, j, nevals_left);
    printf("step              x1                 x2         x1cnt  ck  x2cnt  ck  brack   ck   "
           "x2-x1\n");
    ii = 0;
  }

  /* Initialize portion of ritz we will use (use max double so scanmin will work
     properly when called later on) */
  pntr = &ritz[1];
  for (i = j; i; i--) {
    *pntr++ = DOUBLE_MAX;
  }

  tot_steps = 0;

  /* find evals on left in decreasing index order */
  x2           = Anorm;
  x2cnt        = j;
  numbracketed = j;
  for (index = nevals_left; index >= 1; index--) {
    x1    = 0;
    x1cnt = 0;
    steps = 1; /* ... since started with Anorm bracketing j roots */
    while ((x2 - x1) > tol || numbracketed > 1) {
      x    = 0.5 * (x1 + x2);
      xcnt = sturmcnt(alpha, beta, j, x, workj);
      if (xcnt == -1) {
        return (1);
        /* ... sturmcnt() failed; bail out with error code */
      }
      if (xcnt >= index) {
        x2    = x;
        x2cnt = xcnt;
      }
      else {
        x1    = x;
        x1cnt = xcnt;
      }
      numbracketed = x2cnt - x1cnt;
      steps++;
      if (steps == max_steps) {
        return (2);
        /* ... not converging; bail out with error code */
      }

      if (ritz_sav != NULL) {
        diff = x2 - x1;
        cksturmcnt(ritz_sav, 1, j, x1, x2, &x1ck, &x2ck, &numck);
        printf("%4d %20.16f %20.16f   %3d   %3d  %3d   %3d   %3d   %3d   %g", ii++, x1, x2, x1cnt,
               x1ck, x2cnt, x2ck, numbracketed, numck, diff);
        if (x1cnt != x1ck || x2cnt != x2ck || numbracketed != numck) {
          printf("**\n");
        }
        else {
          printf("\n");
        }
      }
    }
    ritz[index] = 0.5 * (x1 + x2);
    if (ritz_sav != NULL) {
      printf("Ritzval #%d:\n", index);
      printf("            bisection %20.16f\n", ritz[index]);
      printf("                   ql %20.16f\n", ritz_sav[index]);
      printf("           difference %20.16f\n", ritz[index] - ritz_sav[index]);
      printf("---------------------------------------------------\n");
    }
    if (DEBUG_EVECS > 2) {
      printf("    index %d, bisection steps %d, root %20.16f\n", index, steps, ritz[index]);
    }
    tot_steps += steps;
  }

  /* find evals on right in increasing index order */
  x1    = 0;
  x1cnt = 0;
  for (index = j - nevals_right + 1; index <= j; index++) {
    x2    = Anorm;
    x2cnt = j;
    steps = 1; /* ... since started with Anorm bracketing j roots */
    while ((x2 - x1) > tol || numbracketed > 1) {
      x    = 0.5 * (x1 + x2);
      xcnt = sturmcnt(alpha, beta, j, x, workj);
      if (xcnt == -1) {
        return (1);
        /* ... sturmcnt() failed; bail out with error code */
      }
      if (xcnt >= index) {
        x2    = x;
        x2cnt = xcnt;
      }
      else {
        x1    = x;
        x1cnt = xcnt;
      }
      numbracketed = x2cnt - x1cnt;
      steps++;
      if (steps == max_steps) {
        return (2);
        /* ... not converging; bail out with error code */
      }

      if (ritz_sav != NULL) {
        diff = x2 - x1;
        cksturmcnt(ritz_sav, 1, j, x1, x2, &x1ck, &x2ck, &numck);
        printf("%4d %20.16f %20.16f   %3d   %3d  %3d   %3d   %3d   %3d   %g", ii++, x1, x2, x1cnt,
               x1ck, x2cnt, x2ck, numbracketed, numck, diff);
        if (x1cnt != x1ck || x2cnt != x2ck || numbracketed != numck) {
          printf("**\n");
        }
        else {
          printf("\n");
        }
      }
    }
    ritz[index] = 0.5 * (x1 + x2);
    if (ritz_sav != NULL) {
      printf("Ritzval #%d:\n", index);
      printf("            bisection %20.16f\n", ritz[index]);
      printf("                   ql %20.16f\n", ritz_sav[index]);
      printf("           difference %20.16f\n", ritz[index] - ritz_sav[index]);
      printf("---------------------------------------------------\n");
    }
    if (DEBUG_EVECS > 2) {
      printf("    index %d, bisection steps %d, root %20.16f\n", index, steps, ritz[index]);
    }
    tot_steps += steps;
  }
  if (DEBUG_EVECS > 2) {
    printf("  Total number of bisection steps %d.\n", tot_steps);
  }

  return (0); /* ... things seem ok. */
}
