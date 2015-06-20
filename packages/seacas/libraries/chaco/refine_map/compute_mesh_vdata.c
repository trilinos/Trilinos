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

#include "refine_map.h"                 // for refine_vdata
#include "structs.h"                    // for vtx_data


void 
compute_mesh_vdata (
    struct refine_vdata *vdata,	/* preference data for a vertex */
    struct vtx_data **comm_graph,	/* communication graph data structure */
    int vtx,			/* current vertex */
    int *vtx2node,		/* maps graph vtxs to mesh nodes */
    int mesh_dims[3],		/* size of mesh */
    int dim			/* dimension we are currently working in */
)
{
    float     above;		/* my preference to move up in each dimension */
    float     below;		/* my preference to move down in each dimension */
    float     same;		/* my preference to stay where I am */
    int     my_loc;		/* my location in mesh */
    int       neighb_loc;	/* neighbor's location in mesh */
    float     ewgt;		/* weight of an edge */
    int       node;		/* set vertex is assigned to */
    int       neighbor;		/* neighboring vtx in comm_graph */
    int       j;		/* loop counter */

    node = vtx2node[vtx];

    neighb_loc = 0;
    my_loc = 0;

    if (dim == 0) {
	my_loc = node % mesh_dims[0];
    }
    else if (dim == 1) {
	my_loc = (node / mesh_dims[0]) % mesh_dims[1];
    }
    else if (dim == 2) {
	my_loc = node / (mesh_dims[0] * mesh_dims[1]);
    }

    below = above = same = 0;
    for (j = 1; j < comm_graph[vtx]->nedges; j++) {
	neighbor = comm_graph[vtx]->edges[j];
	ewgt = comm_graph[vtx]->ewgts[j];
	node = vtx2node[neighbor];

	if (dim == 0) {
	    neighb_loc = node % mesh_dims[0];
	}
	else if (dim == 1) {
	    neighb_loc = (node / mesh_dims[0]) % mesh_dims[1];
	}
	else if (dim == 2) {
	    neighb_loc = node / (mesh_dims[0] * mesh_dims[1]);
	}

	if (neighb_loc < my_loc)
	    below += ewgt;
	else if (neighb_loc > my_loc)
	    above += ewgt;
	else
	    same += ewgt;
    }
    vdata->below = below;
    vdata->above = above;
    vdata->same = same;
}
