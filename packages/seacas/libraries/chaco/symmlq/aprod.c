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
/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>                      // for NULL
#include "structs.h"

int 
aprod_ (
    long *lnvtxs,
    double *x,
    double *y,
    double *dA,
    double *vwsqrt,
    double *work,
    double *dorthlist		/* vectors to orthogonalize against */
)
{
    int       nvtxs;		/* int copy of long_nvtxs */
    struct vtx_data **A;
    struct orthlink *orthlist;	/* vectors to orthogonalize against */
    void      splarax(), orthog1(), orthogvec(), orthogonalize();

    nvtxs = (int) *lnvtxs;
    A = (struct vtx_data **) dA;
    orthlist = (struct orthlink *) dorthlist;

    /* The offset on x and y is because the arrays come originally from Fortran
       declarations which index from 1 */
    splarax(y - 1, A, nvtxs, x - 1, vwsqrt, work - 1);

    /* Now orthogonalize against lower eigenvectors. */
    if (vwsqrt == NULL)
	orthog1(y - 1, 1, nvtxs);
    else
	orthogvec(y - 1, 1, nvtxs, vwsqrt);
    orthogonalize(y - 1, nvtxs, orthlist);

    return (0);
}
