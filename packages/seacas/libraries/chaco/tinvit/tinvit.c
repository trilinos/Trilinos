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
/* tinvit.f -- translated by f2c (version of 16 May 1991  13:06:06).
 */

#include <float.h>
#include <math.h>

int tinvit(long int *nm, long int *n, double *d, double *e, double *e2, long int *m, double *w,
           long int *ind, double *z, long int *ierr, double *rv1, double *rv2, double *rv3,
           double *rv4, double *rv6)
{
  /* System generated locals */
  long int z_dim1, z_offset, i__1, i__2, i__3;
  double   d__1, d__2, d__3, d__4;

  /* Local variables */
  static double   norm;
  static long int i, j, p, q, r, s;
  static double   u, v, order;
  static long int group;
  static double   x0, x1;
  static long int ii, jj, ip;
  static double   uk, xu;
  static long int tag, its;
  static double   eps2, eps3, eps4;

  /*     this subroutine is a translation of the inverse iteration tech- */
  /*     nique in the algol procedure tristurm by peters and wilkinson. */
  /*     handbook for auto. comp., vol.ii-linear algebra, 418-439(1971). */

  /*     this subroutine finds those eigenvectors of a tridiagonal */
  /*     symmetric matrix corresponding to specified eigenvalues, */
  /*     using inverse iteration. */

  /*     on input */

  /*        nm must be set to the row dimension of two-dimensional */
  /*          array parameters as declared in the calling program */
  /*          dimension statement. */

  /*        n is the order of the matrix. */

  /*        d contains the diagonal elements of the input matrix. */

  /*        e contains the subdiagonal elements of the input matrix */
  /*          in its last n-1 positions.  e(1) is arbitrary. */

  /*        e2 contains the squares of the corresponding elements of e, */
  /*          with zeros corresponding to negligible elements of e. */
  /*          e(i) is considered negligible if it is not larger than */
  /*          the product of the relative machine precision and the sum */
  /*          of the magnitudes of d(i) and d(i-1).  e2(1) must contain */
  /*          0.0d0 if the eigenvalues are in ascending order, or 2.0d0 */
  /*          if the eigenvalues are in descending order.  if  bisect, */
  /*          tridib, or  imtqlv  has been used to find the eigenvalues, */
  /*          their output e2 array is exactly what is expected here. */

  /*        m is the number of specified eigenvalues. */

  /*        w contains the m eigenvalues in ascending or descending order.
   */

  /*        ind contains in its first m positions the submatrix indices */
  /*          associated with the corresponding eigenvalues in w -- */
  /*          1 for eigenvalues belonging to the first submatrix from */
  /*          the top, 2 for those belonging to the second submatrix, etc.
   */

  /*     on output */

  /*        all input arrays are unaltered. */

  /*        z contains the associated set of orthonormal eigenvectors. */
  /*          any vector which fails to converge is set to zero. */

  /*        ierr is set to */
  /*          zero       for normal return, */
  /*          -r         if the eigenvector corresponding to the r-th */
  /*                     eigenvalue fails to converge in 5 iterations. */

  /*        rv1, rv2, rv3, rv4, and rv6 are temporary storage arrays. */

  /*     calls hypot for  dsqrt(a*a + b*b) . */

  /*     questions and comments should be directed to burton s. garbow, */
  /*     mathematics and computer science div, argonne national laboratory
   */

  /*     this version dated august 1983. */

  /*     ------------------------------------------------------------------
   */

  /* Parameter adjustments */
  --rv6;
  --rv4;
  --rv3;
  --rv2;
  --rv1;
  z_dim1   = *nm;
  z_offset = z_dim1 + 1;
  z -= z_offset;
  --ind;
  --w;
  --e2;
  --e;
  --d;

  /* Function Body */
  *ierr = 0;
  if (*m == 0) {
    goto L1001;
  }
  tag   = 0;
  order = 1. - e2[1];
  q     = 0;
/*     .......... establish and process next submatrix .......... */
L100:
  p = q + 1;

  i__1 = *n;
  for (q = p; q <= i__1; ++q) {
    if (q == *n) {
      goto L140;
    }
    if (e2[q + 1] == 0.) {
      goto L140;
    }
    /* L120: */
  }
/*     .......... find vectors by inverse iteration .......... */
L140:
  ++tag;
  s = 0;

  i__1 = *m;
  for (r = 1; r <= i__1; ++r) {
    if (ind[r] != tag) {
      goto L920;
    }
    its = 1;
    x1  = w[r];
    if (s != 0) {
      goto L510;
    }
    /*     .......... check for isolated root .......... */
    xu = 1.;
    if (p != q) {
      goto L490;
    }
    rv6[p] = 1.;
    goto L870;
  L490:
    norm = (d__1 = d[p], fabs(d__1));
    ip   = p + 1;

    i__2 = q;
    for (i = ip; i <= i__2; ++i) {
      /* L500: */
      /* Computing MAX */
      d__3 = norm, d__4 = (d__1 = d[i], fabs(d__1)) + (d__2 = e[i], fabs(d__2));
      norm = d__3 > d__4 ? d__3 : d__4;
    }
    /*     .......... eps2 is the criterion for grouping, */
    /*                eps3 replaces zero pivots and equal */
    /*                roots are modified by eps3, */
    /*                eps4 is taken very small to avoid overflow .........
    . */
    eps2 = norm * .001;
    eps3 = DBL_EPSILON * fabs(norm);
    uk   = (double)(q - p + 1);
    eps4 = uk * eps3;
    uk   = eps4 / sqrt(uk);
    s    = p;
  L505:
    group = 0;
    goto L520;
  /*     .......... look for close or coincident roots .......... */
  L510:
    if ((d__1 = x1 - x0, fabs(d__1)) >= eps2) {
      goto L505;
    }
    ++group;
    if (order * (x1 - x0) <= 0.) {
      x1 = x0 + order * eps3;
    }
  /*     .......... elimination with interchanges and */
  /*                initialization of vector .......... */
  L520:
    v = 0.;

    i__2 = q;
    for (i = p; i <= i__2; ++i) {
      rv6[i] = uk;
      if (i == p) {
        goto L560;
      }
      if ((d__1 = e[i], fabs(d__1)) < fabs(u)) {
        goto L540;
      }
      /*     .......... warning -- a divide check may occur here if */
      /*                e2 array has not been specified correctly ......
      .... */
      xu         = u / e[i];
      rv4[i]     = xu;
      rv1[i - 1] = e[i];
      rv2[i - 1] = d[i] - x1;
      rv3[i - 1] = 0.;
      if (i != q) {
        rv3[i - 1] = e[i + 1];
      }
      u = v - xu * rv2[i - 1];
      v = -xu * rv3[i - 1];
      goto L580;
    L540:
      xu         = e[i] / u;
      rv4[i]     = xu;
      rv1[i - 1] = u;
      rv2[i - 1] = v;
      rv3[i - 1] = 0.;
    L560:
      u = d[i] - x1 - xu * v;
      if (i != q) {
        v = e[i + 1];
      }
    L580:;
    }

    if (u == 0.) {
      u = eps3;
    }
    rv1[q] = u;
    rv2[q] = 0.;
    rv3[q] = 0.;
  /*     .......... back substitution */
  /*                for i=q step -1 until p do -- .......... */
  L600:
    i__2 = q;
    for (ii = p; ii <= i__2; ++ii) {
      i      = p + q - ii;
      rv6[i] = (rv6[i] - u * rv2[i] - v * rv3[i]) / rv1[i];
      v      = u;
      u      = rv6[i];
      /* L620: */
    }
    /*     .......... orthogonalize with respect to previous */
    /*                members of group .......... */
    if (group == 0) {
      goto L700;
    }
    j = r;

    i__2 = group;
    for (jj = 1; jj <= i__2; ++jj) {
    L630:
      --j;
      if (ind[j] != tag) {
        goto L630;
      }
      xu = 0.;

      i__3 = q;
      for (i = p; i <= i__3; ++i) {
        /* L640: */
        xu += rv6[i] * z[i + j * z_dim1];
      }

      i__3 = q;
      for (i = p; i <= i__3; ++i) {
        /* L660: */
        rv6[i] -= xu * z[i + j * z_dim1];
      }

      /* L680: */
    }

  L700:
    norm = 0.;

    i__2 = q;
    for (i = p; i <= i__2; ++i) {
      /* L720: */
      norm += (d__1 = rv6[i], fabs(d__1));
    }

    if (norm >= 1.) {
      goto L840;
    }
    /*     .......... forward substitution .......... */
    if (its == 5) {
      goto L830;
    }
    if (norm != 0.) {
      goto L740;
    }
    rv6[s] = eps4;
    ++s;
    if (s > q) {
      s = p;
    }
    goto L780;
  L740:
    xu = eps4 / norm;

    i__2 = q;
    for (i = p; i <= i__2; ++i) {
      /* L760: */
      rv6[i] *= xu;
    }
  /*     .......... elimination operations on next vector */
  /*                iterate .......... */
  L780:
    i__2 = q;
    for (i = ip; i <= i__2; ++i) {
      u = rv6[i];
      /*     .......... if rv1(i-1) .eq. e(i), a row interchange */
      /*                was performed earlier in the */
      /*                triangularization process .......... */
      if (rv1[i - 1] != e[i]) {
        goto L800;
      }
      u          = rv6[i - 1];
      rv6[i - 1] = rv6[i];
    L800:
      rv6[i] = u - rv4[i] * rv6[i - 1];
      /* L820: */
    }

    ++its;
    goto L600;
    /*     .......... set error -- non-converged eigenvector .......... */

  L830:
    *ierr = -r;
    xu    = 0.;
    goto L870;
  /*     .......... normalize so that sum of squares is */
  /*                1 and expand to full order .......... */
  L840:
    u = 0.;

    i__2 = q;
    for (i = p; i <= i__2; ++i) {
      /* L860: */
      u = hypot(u, rv6[i]);
    }

    xu = 1. / u;

  L870:
    i__2 = *n;
    for (i = 1; i <= i__2; ++i) {
      /* L880: */
      z[i + r * z_dim1] = 0.;
    }

    i__2 = q;
    for (i = p; i <= i__2; ++i) {
      /* L900: */
      z[i + r * z_dim1] = rv6[i] * xu;
    }

    x0 = x1;
  L920:;
  }

  if (q < *n) {
    goto L100;
  }
L1001:
  return 0;
} /* tinvit_ */
