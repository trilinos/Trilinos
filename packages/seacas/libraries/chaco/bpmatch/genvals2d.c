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

#include "params.h"                     // for MAXSETS
#include "smalloc.h"                    // for smalloc


void 
genvals2d (
/* Create lists of sets of values to be sorted. */
    double **xvecs,		/* vectors to partition */
    double *vals[4][MAXSETS],	/* ptrs to lists of values */
    int nvtxs		/* number of values */
)
{
    int       nlists = 4;	/* number of lists to generate */
    double   *temp[4];		/* place holders for vals */
    int       i;		/* loop counter */

    for (i = 0; i < nlists; i++) {
	temp[i] = smalloc(nvtxs * sizeof(double));
    }

    for (i = 1; i <= nvtxs; i++) {
	temp[0][i - 1] = 4 * xvecs[1][i];
	temp[1][i - 1] = 4 * xvecs[2][i];
	temp[2][i - 1] = 4 * (xvecs[1][i] + xvecs[2][i]);
	temp[3][i - 1] = 4 * (xvecs[2][i] - xvecs[1][i]);
    }

    vals[0][1] = vals[1][0] = vals[2][3] = vals[3][2] = temp[0];
    vals[0][2] = vals[2][0] = vals[1][3] = vals[3][1] = temp[1];
    vals[0][3] = vals[3][0] = temp[2];
    vals[1][2] = vals[2][1] = temp[3];
}
