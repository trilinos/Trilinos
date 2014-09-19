/*
 * Copyright (c) 2014, Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Governement retains certain rights in this software.
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
/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include	<stdio.h>
#include	<stdlib.h>
#include	<limits.h>

#ifdef __STDC__
#include	<float.h>
#endif


/* Returns some machine/compiler specific values. */
/* Note: These values might not be calculated precisely correctly.*/
/*       If you know them for your machine, replace this code! */

void      machine_params(double *double_epsilon, double *double_max)
                         	/* returns machine precision */
                     		/* returns largest double value */
{

#ifndef DBL_EPSILON
    double    eps;		/* machine precision */
#endif

#ifndef DBL_MAX

#ifndef DBL_MIN
    double   double_min;	/* smallest double value */
    double    min, min_prev;	/* values halved to compute double_min */
#endif

    double    max;		/* largest double precision value */
#endif

#ifndef DBL_EPSILON
    eps = 1.0 / 16.0;
    while (1.0 + eps > 1.0)
	eps /= 2.0;
    *double_epsilon = eps * 2.0;
#else
    *double_epsilon = DBL_EPSILON;
#endif

#ifndef DBL_MAX

#ifndef DBL_MIN
    min_prev = min = 1.0;
    while (min * 1.0 > 0) {
	min_prev = min;
	min /= 32.0;
    }
    min = min_prev;
    while (min * 1.0 > 0) {
	min_prev = min;
	min /= 2.0;
    }
    double_min = min_prev / (*double_epsilon);
#else
    double_min = DBL_MIN;
#endif

    max = 2.0 * (1.0 - *double_epsilon) / (double_min);
/*
   two_max = max*2.0;
   while (two_max/2.0 == max) {
      max = two_max;
      two_max = max*2.0;
   }
*/
    *double_max = max;
#else
    *double_max = DBL_MAX;
#endif

}
