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

#include <stdio.h>                      // for fprintf, printf, NULL, FILE
#include "structs.h"

/* Check an extended eigenpair of A by direct multiplication. Uses
   the Ay = extval*y + Dg form of the problem for convenience. */

double 
checkeig_ext (
    double *err,
    double *work,			/* work vector of length n */
    struct vtx_data **A,
    double *y,
    int n,
    double extval,
    double *vwsqrt,
    double *gvec,
    double eigtol,
    int warnings		/* don't want to see warning messages in one of the
				   contexts this is called */
)
{
    extern FILE *Output_File;	/* output file or null */
    extern int DEBUG_EVECS;	/* print debugging output? */
    extern int WARNING_EVECS;	/* print warning messages? */
    double    resid;		/* the extended eigen residual */
    double    ch_norm();		/* vector norm */
    void      splarax();	/* sparse matrix vector mult */
    void      scadd();		/* scaled vector add */
    void      scale_diag();	/* scale vector by another's elements */
    void      cpvec();		/* vector copy */

    splarax(err, A, n, y, vwsqrt, work);
    scadd(err, 1, n, -extval, y);
    cpvec(work, 1, n, gvec);	/* only need if going to re-use gvec */
    scale_diag(work, 1, n, vwsqrt);
    scadd(err, 1, n, -1.0, work);
    resid = ch_norm(err, 1, n);

    if (DEBUG_EVECS > 0) {
	printf("  extended residual: %g\n", resid);
	if (Output_File != NULL) {
	    fprintf(Output_File, "  extended residual: %g\n", resid);
	}
    }
    if (warnings && WARNING_EVECS > 0 && resid > eigtol) {
	printf("WARNING: Extended residual (%g) greater than tolerance (%g).\n",
	       resid, eigtol);
	if (Output_File != NULL) {
	    fprintf(Output_File,
		  "WARNING: Extended residual (%g) greater than tolerance (%g).\n",
		    resid, eigtol);
	}
    }

    return (resid);
}
