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

#include "structs.h"                    // for vtx_data


void 
clear_dvals (
    struct vtx_data **graph,	/* data structure for graph */
    int nvtxs,		/* number of vtxs in graph */
    int *ldvals,		/* d-values for each transition */
    int *rdvals,		/* d-values for each transition */
    int *bspace,		/* list of activated vertices */
    int list_length		/* number of activated vertices */
)
{
    int      *edges;		/* loops through edge list */
    int       vtx;		/* vertex in bspace */
    int       neighbor;		/* neighbor of vtx */
    int       i, j;		/* loop counters */

    if (list_length > .05 * nvtxs) {	/* Do it directly. */
	for (i = 1; i <= nvtxs; i++) {
	    ldvals[i] = rdvals[i] = 0;
	}
    }
    else {			/* Do it more carefully */
	for (i = 0; i < list_length; i++) {
	    vtx = bspace[i];
	    if (vtx < 0)
		vtx = -vtx;
	    ldvals[vtx] = rdvals[vtx] = 0;
	    edges = graph[vtx]->edges;
	    for (j = graph[vtx]->nedges - 1; j; j--) {
		neighbor = *(++edges);
		ldvals[neighbor] = rdvals[neighbor] = 0;
	    }
	}
    }
}
