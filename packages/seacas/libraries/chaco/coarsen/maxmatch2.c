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

#include "defs.h"                       // for FALSE, TRUE
#include "smalloc.h"                    // for sfree, smalloc
#include "structs.h"                    // for vtx_data


/* Find a maximal matching in a graph using simple greedy algorithm. */
/* Randomly permute vertices, and then have each select first unmatched */
/* neighbor. */

int 
maxmatch2 (
    struct vtx_data **graph,	/* array of vtx data for graph */
    int nvtxs,		/* number of vertices in graph */
    int *mflag,		/* flag indicating vtx selected or not */
    int using_ewgts		/* are edge weights being used? */
)
{
    extern int HEAVY_MATCH;	/* encourage heavy matching edges? */
    float     ewgt_max;		/* heaviest edge seen so far */
    int      *order;		/* random ordering of vertices */
    int      *iptr, *jptr;	/* loops through integer arrays */
    int       matched;		/* has vertex been matched? */
    int       vtx;		/* vertex to process next */
    int       neighbor;		/* neighbor of a vertex */
    int       nmerged;		/* number of edges in matching */
    int       jsave;		/* best edge so far */
    int       i, j;		/* loop counters */

    void      randomize();

    /* First, randomly permute the vertices. */
    iptr = order = smalloc((nvtxs + 1) * sizeof(int));
    jptr = mflag;
    for (i = 1; i <= nvtxs; i++) {
	*(++iptr) = i;
	*(++jptr) = 0;
    }
    randomize(order, nvtxs);

    nmerged = 0;

    if (!using_ewgts || !HEAVY_MATCH) {	/* Simple greedy approach. */
	for (i = 1; i <= nvtxs; i++) {
	    vtx = order[i];
	    if (mflag[vtx] == 0) {	/* Not already matched. */
		/* Find first unmatched neighbor of vtx. */
		matched = FALSE;
		for (j = 1; j < graph[vtx]->nedges && !matched; j++) {
		    neighbor = graph[vtx]->edges[j];
		    if (mflag[neighbor] == 0) {	/* Found one! */
			mflag[vtx] = neighbor;
			mflag[neighbor] = vtx;
			matched = TRUE;
			nmerged++;
		    }
		}
	    }
	}
    }

    else {			/* Find heavy edge to match */
	for (i = 1; i <= nvtxs; i++) {
	    vtx = order[i];
	    if (mflag[vtx] == 0) {	/* Not already matched. */
		/* Find heaviest unmatched neighbor of vtx. */
		jsave = 0;
		ewgt_max = 0;
		for (j = 1; j < graph[vtx]->nedges; j++) {
		    neighbor = graph[vtx]->edges[j];
		    if (mflag[neighbor] == 0 && graph[vtx]->ewgts[j] > ewgt_max) {
			ewgt_max = graph[vtx]->ewgts[j];
			jsave = j;
		    }
		}
		if (jsave > 0) {
		    neighbor = graph[vtx]->edges[jsave];
		    mflag[vtx] = neighbor;
		    mflag[neighbor] = vtx;
		    nmerged++;
		}
	    }
	}
    }

    sfree(order);
    return (nmerged);
}
