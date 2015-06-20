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

#include <stdio.h>                      // for NULL
#include "defs.h"                       // for TRUE
#include "smalloc.h"                    // for sfree, smalloc_ret
#include "structs.h"


/* Given a partition, refine the mapping in a locally greedy fasion. */

void 
refine_map (
    struct vtx_data **graph,	/* graph data structure */
    int nvtxs,		/* number of vertices in graph */
    int using_ewgts,		/* are edge weights being used? */
    int *assign,		/* current assignment */
    int cube_or_mesh,		/* 0 => hypercube, d => d-dimensional mesh */
    int ndims_tot,		/* if hypercube, number of dimensions */
    int mesh_dims[3]		/* if mesh, dimensions of mesh */
)
{
    struct vtx_data **comm_graph;	/* graph for communication requirements */
    int       nsets_tot=0;	/* total number of sets */
    int    *vtx2node = NULL;	/* mapping of comm_graph vtxs to processors */
    int    *node2vtx = NULL;	/* mapping of sets to comm_graph vtxs */
    double    maxdesire=0.0;	/* largest possible desire to flip an edge */
    int       error=0;		/* out of space? */
    int       i;		/* loop counter */

    double    find_maxdeg();
    void      free_graph(), strout();
    int       make_comm_graph(), refine_mesh(), refine_cube();

    if (cube_or_mesh == 0)
	nsets_tot = 1 << ndims_tot;
    else if (cube_or_mesh == 1)
	nsets_tot = mesh_dims[0];
    else if (cube_or_mesh == 2)
	nsets_tot = mesh_dims[0] * mesh_dims[1];
    else if (cube_or_mesh == 3)
	nsets_tot = mesh_dims[0] * mesh_dims[1] * mesh_dims[2];

    node2vtx = vtx2node = NULL;

    /* Construct the weighted quotient graph representing communication. */
    error = make_comm_graph(&comm_graph, graph, nvtxs, using_ewgts, assign, nsets_tot);

    if (!error) {
        maxdesire = 2 * find_maxdeg(comm_graph, nsets_tot, TRUE, (float *) NULL);

        vtx2node = smalloc_ret((nsets_tot + 1) * sizeof(int));
        node2vtx = smalloc_ret(nsets_tot * sizeof(int));
	if (node2vtx == NULL || vtx2node == NULL) {
	    error = 1;
	    goto skip;
	}

        for (i = 1; i <= nsets_tot; i++) {
	    vtx2node[i] = (int) i - 1;
	    node2vtx[i - 1] = (int) i;
        }

        if (cube_or_mesh > 0) {
	    error = refine_mesh(comm_graph, cube_or_mesh, mesh_dims,
		maxdesire, vtx2node, node2vtx);
        }

        else if (cube_or_mesh == 0) {
	    error = refine_cube(comm_graph, ndims_tot, maxdesire,
		 vtx2node, node2vtx);
        }

	if (!error) {
            for (i = 1; i <= nvtxs; i++) {
		assign[i] = vtx2node[assign[i] + 1];
	    }
	}
    }

skip:

    if (error) {
	strout("\nWARNING: No space to refine mapping to processors.");
	strout("         NO MAPPING REFINEMENT PERFORMED.\n");
    }

    sfree(node2vtx);
    sfree(vtx2node);
    free_graph(comm_graph);
}
