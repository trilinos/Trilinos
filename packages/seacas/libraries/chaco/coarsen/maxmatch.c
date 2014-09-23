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

#include <stdio.h>                      // for printf
#include "structs.h"


/* Find a maximal matching in a graph using one of several algorithms. */

int 
maxmatch (
    struct vtx_data **graph,	/* array of vtx data for graph */
    int nvtxs,		/* number of vertices in graph */
    int nedges,		/* number of edges in graph */
    int *mflag,		/* flag indicating vtx selected or not */
    int using_ewgts,		/* are edge weights being used? */
    int igeom,		/* geometric dimensionality */
    float **coords		/* coordinates for each vertex */
)
{
    extern int DEBUG_COARSEN;	/* debug output for coarsening? */
    extern int MATCH_TYPE;	/* which matching routine to use */
    int       nmerged=0;		/* number of matching edges found */
    int       maxmatch1(), maxmatch2(), maxmatch3(), maxmatch4();
    int       maxmatch5(), maxmatch9();

    if (MATCH_TYPE == 1) {	/* Dumb, fast routine. */
	nmerged = maxmatch1(graph, nvtxs, mflag, using_ewgts);
    }

    else if (MATCH_TYPE == 2) {	/* More random but somewhat slower. */
	nmerged = maxmatch2(graph, nvtxs, mflag, using_ewgts);
    }

    else if (MATCH_TYPE == 3) {	/* Much more random but slower still. */
	nmerged = maxmatch3(graph, nvtxs, mflag, using_ewgts);
    }

    else if (MATCH_TYPE == 4) {	/* Truly random but very slow. */
	nmerged = maxmatch4(graph, nvtxs, nedges, mflag, using_ewgts);
    }
    else if (MATCH_TYPE == 5) {	/* Geometric nearness. */
	nmerged = maxmatch5(graph, nvtxs, mflag, igeom, coords);
    }

    else if (MATCH_TYPE == 9) {	/* Minimum degree of merged vertex */
	nmerged = maxmatch9(graph, nvtxs, mflag, using_ewgts);
    }

    if (DEBUG_COARSEN > 0) {
	printf("Number of matching edges = %d\n", nmerged);
    }

    return (nmerged);
}
