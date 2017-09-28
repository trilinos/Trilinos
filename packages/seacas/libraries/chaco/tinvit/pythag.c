/*
 * Copyright (c) 2005 National Technology & Engineering Solutions
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
/* pythag.f -- translated by f2c (version of 16 May 1991  13:06:06).
   You must link the resulting object file with the libraries:
        -link <S|C|M|L>f2c.lib   (in that order)
*/

#include "f2c.h"

doublereal pythag_(doublereal *a, doublereal *b)
{
  /* System generated locals */
  doublereal ret_val, d__1, d__2, d__3;

  /* Local variables */
  static doublereal p, r, s, t, u;

  /*     finds dsqrt(a**2+b**2) without overflow or destructive underflow */

  /* Computing MAX */
  d__1 = abs(*a), d__2 = abs(*b);
  p = max(d__1, d__2);
  if (p == 0.) {
    goto L20;
  }
  /* Computing MIN */
  d__2 = abs(*a), d__3 = abs(*b);
  /* Computing 2nd power */
  d__1 = min(d__2, d__3) / p;
  r    = d__1 * d__1;
L10:
  t = r + 4.;
  if (t == 4.) {
    goto L20;
  }
  s = r / t;
  u = s * 2. + 1.;
  p = u * p;
  /* Computing 2nd power */
  d__1 = s / u;
  r    = d__1 * d__1 * r;
  goto L10;
L20:
  ret_val = p;
  return ret_val;
} /* pythag_ */
