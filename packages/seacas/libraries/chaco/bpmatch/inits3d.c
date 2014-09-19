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
#include "structs.h"                    // for vtx_data


void 
inits3d (
    struct vtx_data **graph,	/* graph data structure for vertex weights */
    double **xvecs,		/* values to partition with */
    double *vals[8][MAXSETS],	/* values in sorted lists */
    int *indices[8][MAXSETS],	/* indices sorting lists */
    int nvtxs,		/* number of vertices */
    double *dist,			/* trial separation point */
    int startvtx[8][MAXSETS],	/* indices defining separation */
    double *size,			/* size of each set being modified */
    int *sets			/* set each vertex gets assigned to */
)
{
    double    xmid, ymid, zmid;	/* median x, y and z values */
    double    val, bestval;	/* values for determining set preferences */
    int     bestset = 0;	/* set vertex wants to be in */
    int       signx, signy, signz;	/* sign values for different target points */
    int       nsets = 8;	/* number of different sets */
    int       i, j;		/* loop counters */
    int       findindex();

/*
    xmid = .25 * (vals[0][1][indices[0][1][nvtxs / 2]] +
		  vals[0][1][indices[0][1][nvtxs / 2 - 1]]);
    ymid = .25 * (vals[0][2][indices[0][2][nvtxs / 2]] +
		  vals[0][2][indices[0][2][nvtxs / 2 - 1]]);
    zmid = .25 * (vals[0][4][indices[0][4][nvtxs / 2]] +
		  vals[0][4][indices[0][4][nvtxs / 2 - 1]]);
*/
    xmid = .5 * vals[0][1][indices[0][1][nvtxs / 2]];
    ymid = .5 * vals[0][2][indices[0][2][nvtxs / 2]];
    zmid = .5 * vals[0][4][indices[0][4][nvtxs / 2]];

    dist[0] = -xmid - ymid - zmid;
    dist[1] = xmid - ymid - zmid;
    dist[2] = -xmid + ymid - zmid;
    dist[3] = xmid + ymid - zmid;
    dist[4] = -xmid - ymid + zmid;
    dist[5] = xmid - ymid + zmid;
    dist[6] = -xmid + ymid + zmid;
    dist[7] = xmid + ymid + zmid;

    /* Now initialize startvtxs. */
    startvtx[0][1] = startvtx[2][3] = startvtx[4][5] = startvtx[6][7] = nvtxs / 2;
    startvtx[0][2] = startvtx[1][3] = startvtx[4][6] = startvtx[5][7] = nvtxs / 2;
    startvtx[0][4] = startvtx[1][5] = startvtx[2][6] = startvtx[3][7] = nvtxs / 2;
    startvtx[0][3] = startvtx[4][7] =
       findindex(indices[0][3], vals[0][3], dist[3] - dist[0], nvtxs);
    startvtx[1][2] = startvtx[5][6] =
       findindex(indices[1][2], vals[1][2], dist[2] - dist[1], nvtxs);
    startvtx[0][5] = startvtx[2][7] =
       findindex(indices[0][5], vals[0][5], dist[5] - dist[0], nvtxs);
    startvtx[1][4] = startvtx[3][6] =
       findindex(indices[1][4], vals[1][4], dist[4] - dist[1], nvtxs);
    startvtx[0][6] = startvtx[1][7] =
       findindex(indices[0][6], vals[0][6], dist[6] - dist[0], nvtxs);
    startvtx[2][4] = startvtx[3][5] =
       findindex(indices[2][4], vals[2][4], dist[4] - dist[2], nvtxs);
    startvtx[0][7] = findindex(indices[0][7], vals[0][7], dist[7] - dist[0], nvtxs);
    startvtx[1][6] = findindex(indices[1][6], vals[1][6], dist[6] - dist[1], nvtxs);
    startvtx[2][5] = findindex(indices[2][5], vals[2][5], dist[5] - dist[2], nvtxs);
    startvtx[3][4] = findindex(indices[3][4], vals[3][4], dist[4] - dist[3], nvtxs);

    /* Finally, determine the set sizes based on this splitter. */

    for (i = 0; i < nsets; i++)
	size[i] = 0;

    for (i = 1; i <= nvtxs; i++) {
	/* Which set is this vertex in? */
	signx = signy = signz = -1;
	bestval = 0;
	for (j = 0; j < nsets; j++) {
	    val = -dist[j] + 2 * (signx * xvecs[1][i] + signy * xvecs[2][i] + signz * xvecs[3][i]);
	    if (j == 0 || val < bestval) {
		bestval = val;
		bestset = (int) j;
	    }
	    if (signx == 1 && signy == 1)
		signz *= -1;
	    if (signx == 1)
		signy *= -1;
	    signx *= -1;
	}
	sets[i] = bestset;
	size[bestset] += graph[i]->vwgt;
    }
}
