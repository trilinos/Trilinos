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

#include "smalloc.h"                    // for smalloc
#include "structs.h"                    // for vtx_data

/* Make coarse graph vertex coordinates be center-of-mass of their */
/* fine graph constituents. */

void 
makeccoords (
    struct vtx_data **graph,	/* array of vtx data for graph */
    int cnvtxs,		/* number of vertices in coarse graph */
    int *cv2v_ptrs,		/* vtxs corresponding to each cvtx */
    int *cv2v_vals,		/* indices into cv2v_vals */
    int igeom,		/* dimensions of geometric data */
    float **coords,		/* coordinates for vertices */
    float **ccoords		/* coordinates for coarsened vertices */
)
{
    double    mass;		/* total mass of merged vertices */
    float    *cptr;		/* loops through ccoords */
    int       cvtx;		/* coarse graph vertex */
    int       vtx;		/* vertex being merged */
    int       i, j;		/* loop counters */

    for (i = 0; i < igeom; i++) {
	ccoords[i] = cptr = smalloc((cnvtxs + 1) * sizeof(float));
	for (cvtx = cnvtxs; cvtx; cvtx--) {
	    *(++cptr) = 0;
	}
    }
    if (igeom == 1) {
	for (cvtx = 1; cvtx <= cnvtxs; cvtx++) {
	    mass = 0;
	    for (j = cv2v_ptrs[cvtx]; j < cv2v_ptrs[cvtx + 1]; j++) {
		vtx = cv2v_vals[j];
		mass += graph[vtx]->vwgt;
		ccoords[0][cvtx] += graph[vtx]->vwgt * coords[0][vtx];
	    }
	    ccoords[0][cvtx] /= mass;
	}
    }
    else if (igeom == 2) {
	for (cvtx = 1; cvtx <= cnvtxs; cvtx++) {
	    mass = 0;
	    for (j = cv2v_ptrs[cvtx]; j < cv2v_ptrs[cvtx + 1]; j++) {
		vtx = cv2v_vals[j];
		mass += graph[vtx]->vwgt;
		ccoords[0][cvtx] += graph[vtx]->vwgt * coords[0][vtx];
		ccoords[1][cvtx] += graph[vtx]->vwgt * coords[1][vtx];
	    }
	    ccoords[0][cvtx] /= mass;
	    ccoords[1][cvtx] /= mass;
	}
    }
    else if (igeom > 2) {
	for (cvtx = 1; cvtx <= cnvtxs; cvtx++) {
	    mass = 0;
	    for (j = cv2v_ptrs[cvtx]; j < cv2v_ptrs[cvtx + 1]; j++) {
		vtx = cv2v_vals[j];
		mass += graph[vtx]->vwgt;
		ccoords[0][cvtx] += graph[vtx]->vwgt * coords[0][vtx];
		ccoords[1][cvtx] += graph[vtx]->vwgt * coords[1][vtx];
		ccoords[2][cvtx] += graph[vtx]->vwgt * coords[2][vtx];
	    }
	    ccoords[0][cvtx] /= mass;
	    ccoords[1][cvtx] /= mass;
	    ccoords[2][cvtx] /= mass;
	}
    }
}
