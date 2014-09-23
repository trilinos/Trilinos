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

#include <stdio.h>                      // for NULL
#include "smalloc.h"                    // for sfree, smalloc
#include "structs.h"                    // for ipairs

static struct ipairs *pedges;	/* perturbed edges */
static double *pvals;		/* perturbed values */

/* Inititialize the perturbation */
void 
perturb_init (
    int n			/* graph size at this level */
)
{
    extern int NPERTURB;	/* number of edges to perturb */
    extern double PERTURB_MAX;	/* maximum perturbation */
    int       i, j;		/* loop counter */
    double    drandom();

    /* Initialize the diagonal perturbation weights */
    pedges = smalloc(NPERTURB * sizeof(struct ipairs));
    pvals = smalloc(NPERTURB * sizeof(double));

    if (n <= 1) {
	for (i = 0; i < NPERTURB; i++) {
	    pedges[i].val1 = pedges[i].val2 = 0;
	    pvals[i] = 0;
	}
	return;
    }

    for (i = 0; i < NPERTURB; i++) {
	pedges[i].val1 = 1 + (n * drandom());

	/* Find another vertex to define an edge. */
	j = 1 + (n * drandom());
	while (j == i)
	    j = 1 + (n * drandom());
	pedges[i].val2 = 1 + (n * drandom());

	pvals[i] = PERTURB_MAX * drandom();
    }
}

void 
perturb_clear (void)
{


    sfree(pedges);
    sfree(pvals);
    pedges = NULL;
    pvals = NULL;
}


/* Modify the result of splarax to break any graph symmetry */
void 
perturb (
    double *result,		/* result of matrix-vector multiply */
    double *vec			/* vector matrix multiplies */
)
{
    extern int NPERTURB;	/* number of edges to perturb */
    int       i;		/* loop counter */

    for (i = 0; i < NPERTURB; i++) {
	result[pedges[i].val1] +=
	   pvals[i] * (vec[pedges[i].val2] - vec[pedges[i].val1]);
	result[pedges[i].val2] +=
	   pvals[i] * (vec[pedges[i].val1] - vec[pedges[i].val2]);
    }
}


/* Modify the result of splarax to break any graph symmetry, float version */
void 
perturb_float (
    float *result,		/* result of matrix-vector multiply */
    float *vec			/* vector matrix multiplies */
)
{
    extern int NPERTURB;	/* number of edges to perturb */
    int       i;		/* loop counter */

    for (i = 0; i < NPERTURB; i++) {
	result[pedges[i].val1] +=
	   pvals[i] * (vec[pedges[i].val2] - vec[pedges[i].val1]);
	result[pedges[i].val2] +=
	   pvals[i] * (vec[pedges[i].val1] - vec[pedges[i].val2]);
    }
}
