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

#include <stdio.h>                      // for printf, NULL
#include "params.h"                     // for MAXSETS
#include "smalloc.h"                    // for sfree, smalloc
#include "structs.h"                    // for vtx_data


/* Partition vertices into sets in one of several simplistic ways. */
void 
simple_part (
    struct vtx_data **graph,	/* data structure for graph */
    int nvtxs,		/* total number of vtxs in graph */
    int *sets,			/* sets vertices get assigned to */
    int nsets,		/* number of sets at each division */
    int simple_type,		/* type of decomposition */
    double *goal			/* desired set sizes */
)
{
    extern int DEBUG_TRACE;	/* trace the execution of the code */
    double    cutoff;		/* ending weight for a partition */
    double    ratio;		/* weight/goal */
    double    best_ratio;	/* lowest ratio of weight/goal */
    double    sum;              /* sum of vwgts in a set */
    double    vwgt;             /* vertex weight */
    int       using_vwgts;	/* are vertex weights active? */
    int      *order;		/* random ordering of vertices */
    int       weight;		/* sum of vertex weights in a partition */
    int       wgts[MAXSETS];	/* weight assigned to given set so far */
    int       set=-1;		/* set vertex is assigned to */
    int       i, j;		/* loop counters */

    void      randomize();

    using_vwgts = (graph != NULL);

    /* Scattered initial decomposition. */
    if (simple_type == 1) {
	if (DEBUG_TRACE > 0) {
	    printf("Generating scattered partition, nvtxs = %d\n", nvtxs);
	}
	for (j=0; j<nsets; j++) wgts[j] = 0;
	for (i = 1; i <= nvtxs; i++) {
	    best_ratio = 2;
	    for (j=0; j<nsets; j++) {
	        ratio = wgts[j] / goal[j];
	        if (ratio < best_ratio) {
		    best_ratio = ratio;
		    set = (int) j;
		}
	    }
	    if (using_vwgts) {
	        wgts[set] += graph[i]->vwgt;
	    }
	    else {
	        wgts[set]++;
	    }
	    sets[i] = set;
	}
    }

    /* Random initial decomposition. */
    if (simple_type == 2) {
	if (DEBUG_TRACE > 0) {
	    printf("Generating random partition, nvtxs = %d\n", nvtxs);
	}
	/* Construct random order in which to step through graph. */
	order = smalloc((nvtxs + 1) * sizeof(int));
	for (i = 1; i <= nvtxs; i++)
	    order[i] = i;
	randomize(order, nvtxs);

	weight = 0;
	cutoff = goal[0];
	set = 0;
	vwgt = 1;
	sum = 0;
	for (i = 1; i <= nvtxs; i++) {
	    if (using_vwgts)
		vwgt = graph[order[i]]->vwgt;

            if (set < nsets - 1 &&
	        (weight >= cutoff  ||
	         (weight + vwgt >= cutoff &&
		   sum + vwgt - goal[set] > vwgt - goal[set + 1] ))) {
		cutoff += goal[++set];
		sum = 0;
	    }

	    weight += vwgt;
	    sum += vwgt;

	    sets[order[i]] = set;
	}
	sfree(order);
    }

    /* Linearly ordered initial decomposition. */
    if (simple_type == 3) {
	if (DEBUG_TRACE > 0) {
	    printf("Generating linear partition, nvtxs = %d\n", nvtxs);
	}
	weight = 0;
	cutoff = goal[0];
	set = 0;
	vwgt = 1;
	sum = 0;
	for (i = 1; i <= nvtxs; i++) {
	    if (using_vwgts)
		vwgt = graph[i]->vwgt;

            if (set < nsets - 1 &&
	        (weight >= cutoff  ||
	         (weight + vwgt >= cutoff &&
		   sum + vwgt - goal[set] > vwgt - goal[set + 1] ))) {
		cutoff += goal[++set];
		sum = 0;
	    }

	    weight += vwgt;
	    sum += vwgt;

	    sets[i] = set;
	}
    }
}
