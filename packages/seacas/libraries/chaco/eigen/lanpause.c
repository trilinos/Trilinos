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
#include "defs.h"                       // for FALSE, TRUE

/* Determine whether to pause in Lanczos */
int 
lanpause (
    int j,			/* current step */
    int lastpause,		/* when last paused */
    int interval,		/* interval between pauses */
    double **q,			/* the Lanczos vectors */
    int n,			/* length of Lanczos vectors */
    int *pausemode,		/* which pausing criterion to use */
    int version,		/* which version of sel. orth. we are using */
    double beta			/* current off-diagonal value */
)
{
    extern int DEBUG_EVECS;	/* debugging level for eigen computation */
    extern double DOUBLE_EPSILON;       /* machine precision */
    double    paige_dot;	/* q[j]^T q[1] */
    double    paigetol;		/* pause if paigedot > paigetol */
    double    dot();		/* standard dot product */
    double    fabs();		/* intrinsic abs. value */
    void      checkorth();


    /* Check orthogonality of last Lanczos vector against previous ones */
    if (DEBUG_EVECS > 3) {
	checkorth(q, n, j);
    }

    /* periodic reorthogonalization */
    if (version == 1 || version == 2) {
	if ((j - lastpause) == interval || beta < 1000 * DOUBLE_EPSILON ) {
	    return (TRUE);
	}
	else {
	    return (FALSE);
	}
    }

    /* Run until orthogonality with first Lanczos vector deteriorates, then switch
       switch to periodic reorthog. */
    if (version == 3) {
	paigetol = 1.0e-3;
	if (*pausemode == 1) {
	    paige_dot = fabs(dot(q[1], 1, n, q[j]));
	    if ((paige_dot > paigetol && j > 1) || beta < 1000 * DOUBLE_EPSILON) {
		if (DEBUG_EVECS > 1) {
		    printf("  Pausing on step %3d with Paige prod. = %g\n", j, paige_dot);
		}
		*pausemode = 2;
		return (TRUE);
	    }
	    else {
		return (FALSE);
	    }
	}
	if (*pausemode == 2) {
	    if ((j - lastpause) == interval || beta < 1000 * DOUBLE_EPSILON) {
		return (TRUE);
	    }
	    else {
		return (FALSE);
	    }
	}
    }

    /* shouldn't ever get this far, but alint really wants a return value */
    return (FALSE);
}


int 
lanpause_float (
    int j,			/* current step */
    int lastpause,		/* when last paused */
    int interval,		/* interval between pauses */
    float **q,			/* the Lanczos vectors */
    int n,			/* length of Lanczos vectors */
    int *pausemode,		/* which pausing criterion to use */
    int version,		/* which version of sel. orth. we are using */
    double beta			/* current off-diagonal value */
)
{
    extern int DEBUG_EVECS;	/* debugging level for eigen computation */
    extern double DOUBLE_EPSILON;       /* machine precision */
    double    paige_dot;	/* q[j]^T q[1] */
    double    paigetol;		/* pause if paigedot > paigetol */
    double    dot_float();	/* standard dot product */
    double    fabs();		/* intrinsic abs. value */
    void      checkorth_float();


    /* Check orthogonality of last Lanczos vector against previous ones */
    if (DEBUG_EVECS > 3) {
	checkorth_float(q, n, j);
    }

    /* periodic reorthogonalization */
    if (version == 1 || version == 2) {
	if ((j - lastpause) == interval || beta < 1000 * DOUBLE_EPSILON ) {
	    return (TRUE);
	}
	else {
	    return (FALSE);
	}
    }

    /* Run until orthogonality with first Lanczos vector deteriorates, then switch
       switch to periodic reorthog. */
    if (version == 3) {
	paigetol = 1.0e-3;
	if (*pausemode == 1) {
	    paige_dot = fabs(dot_float(q[1], 1, n, q[j]));
	    if ((paige_dot > paigetol && j > 1) || beta < 1000 * DOUBLE_EPSILON) {
		if (DEBUG_EVECS > 1) {
		    printf("  Pausing on step %3d with Paige prod. = %g\n", j, paige_dot);
		}
		*pausemode = 2;
		return (TRUE);
	    }
	    else {
		return (FALSE);
	    }
	}
	if (*pausemode == 2) {
	    if ((j - lastpause) == interval || beta < 1000 * DOUBLE_EPSILON ) {
		return (TRUE);
	    }
	    else {
		return (FALSE);
	    }
	}
    }

    /* shouldn't ever get this far, but alint really wants a return value */
    return (FALSE);
}

