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

#include <stdio.h>                      // for FILE
#include "structs.h"

/* Print metrics of partition quality. */

void 
countup (
    struct vtx_data **graph,	/* graph data structure */
    int nvtxs,		/* number of vtxs in graph */
    int *assignment,		/* set number of each vtx (length nvtxs+1) */
    int ndims,		/* number of cuts at each level */
    int architecture,		/* what's the target parallel machine? */
    int ndims_tot,		/* total number of hypercube dimensions */
    int mesh_dims[3],		/* extent of mesh in each dimension */
    int print_lev,		/* level of output */
    FILE *outfile,		/* output file if not NULL */
    int using_ewgts		/* are edge weights being used? */
)
{
    extern int VERTEX_SEPARATOR;/* vertex instead of edge separator? */
    extern int VERTEX_COVER;	/* make/improve vtx separator via matching? */
    void countup_cube(), countup_mesh(), countup_vtx_sep();

    if (VERTEX_SEPARATOR || VERTEX_COVER) {
	countup_vtx_sep(graph, nvtxs, assignment);
    }
    else {
	if (architecture == 0) {
            countup_cube(graph, nvtxs, assignment, ndims, ndims_tot, print_lev, outfile,
		                 using_ewgts);
	}

	else if (architecture > 0) {
            countup_mesh(graph, nvtxs, assignment, mesh_dims, print_lev, outfile,
		                 using_ewgts);
	}
    }
}
