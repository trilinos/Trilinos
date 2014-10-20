/*
 * Copyright (c) 2014, Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
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
 *     * Neither the name of Sandia Corporation nor the names of its
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
/* epslon.f -- translated by f2c (version of 16 May 1991  13:06:06).
   You must link the resulting object file with the libraries:
	-link <S|C|M|L>f2c.lib   (in that order)
*/

#include "f2c.h"

doublereal epslon_(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Local variables */
    static doublereal a, b, c, eps;


/*     estimate unit roundoff in quantities of size x. */


/*     this program should function properly on all systems */
/*     satisfying the following two assumptions, */
/*        1.  the base used in representing floating point */
/*            numbers is not a power of three. */
/*        2.  the quantity  a  in statement 10 is represented to */
/*            the accuracy used in floating point variables */
/*            that are stored in memory. */
/*     the statement number 10 and the go to 10 are intended to */
/*     force optimizing compilers to generate code satisfying */
/*     assumption 2. */
/*     under these assumptions, it should be true that, */
/*            a  is not exactly equal to four-thirds, */
/*            b  has a zero for its last bit or digit, */
/*            c  is not exactly equal to one, */
/*            eps  measures the separation of 1.0 from */
/*                 the next larger floating point number. */
/*     the developers of eispack would appreciate being informed */
/*     about any systems where these assumptions do not hold. */

/*     this version dated 4/6/83. */

    a = 1.3333333333333333;
L10:
    b = a - 1.;
    c = b + b + b;
    eps = (d__1 = c - 1., abs(d__1));
    if (eps == 0.) {
	goto L10;
    }
    ret_val = eps * abs(*x);
    return ret_val;
} /* epslon_ */

