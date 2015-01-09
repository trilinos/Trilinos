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

#include <stdio.h>                      // for printf
#include "structs.h"                    // for orthlink, orthlink_float

/* Print out the orthogonalization set, double version */ 
void 
solistout (
    struct orthlink **solist,	/* vector of pntrs to orthlnks */
    int n,			/* length of vecs to orth. against */
    int ngood,		/* number of good vecs on list */
    int j			/* current number of Lanczos steps */
)
{
    int       i;		/* index */
    extern int DEBUG_EVECS;	/* debugging output level for eigen computations */

    /* to placate alint */
    n = n;

    for (i = 1; i <= ngood; i++) {
	if ((solist[i])->index <= (int) (j / 2)) {
	    printf(".");
	}
	else {
	    printf("+");
	}
	/* Really detailed output: printf("\n"); printf("depth
	   %d\n",(solist[i])->depth); printf("index %d\n",(solist[i])->index);
	   printf("ritzval %g\n",(solist[i])->ritzval); printf("betaji
	   %g\n",(solist[i])->betaji); printf("tau %g\n",(solist[i])->tau);
	   printf("prevtau %g\n",(solist[i])->prevtau);
	   vecout((solist[i])->vec,1,n,"vec", NULL); */
    }
    printf("%d\n", ngood);

    if (DEBUG_EVECS > 2) {
	printf("  actual indicies: ");
	for (i = 1; i <= ngood; i++) {
	    printf(" %2d", solist[i]->index);
	}
	printf("\n");
    }
}

/* Print out the orthogonalization set, float version */ 
void 
solistout_float (
    struct orthlink_float **solist,	/* vector of pntrs to orthlnks */
    int n,			/* length of vecs to orth. against */
    int ngood,		/* number of good vecs on list */
    int j			/* current number of Lanczos steps */
)
{
    int       i;		/* index */

    /* to placate alint */
    n = n;

    for (i = 1; i <= ngood; i++) {
	if ((solist[i])->index <= (int) (j / 2)) {
	    printf(".");
	}
	else {
	    printf("+");
	}
    }
    printf("%d\n", ngood);
}
