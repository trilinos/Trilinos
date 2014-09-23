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

#include "structs.h"

/* Combine goals of collections of processors for next division. */
void 
merge_goals (
    double *goal,			/* desired set sizes */
    double *merged_goal,		/* sizes of sets at this partition level */
    struct set_info *set_info,	/* information about all sets */
    int *subsets,		/* set numbers of processors to merge */
    int nsets,		/* number of sets created by this division */
    int ndims_tot,		/* total number of dimensions in the hypercube */
    int cube_or_mesh,		/* 0=> hypercube, d=> d-dimensional mesh */
    int mesh_dims[3],		/* shape of mesh */
    double vwgt_sum		/* actual sum of vertex weights */
)
{
    struct set_info *set;	/* set of processors still clumped together */
    double    total_goal;	/* total of desired goals */
    int       index;		/* x, y or z location of a processor */
    int       i, x, y, z;	/* loop counters */

    total_goal = 0;
    for (i = 0; i < nsets; i++) {
	set = &(set_info[subsets[i]]);
	merged_goal[i] = 0;

	if (cube_or_mesh > 0) {	/* Mesh architecture. */
	    for (x = set->low[0]; x < set->low[0] + set->span[0]; x++) {
		for (y = set->low[1]; y < set->low[1] + set->span[1]; y++) {
		    for (z = set->low[2]; z < set->low[2] + set->span[2]; z++) {
			index = z * mesh_dims[0] * mesh_dims[1] + y * mesh_dims[0] + x;
			merged_goal[i] += goal[index];
		    }
		}
	    }
	}

	else if (cube_or_mesh == 0) {	/* Hypercube architecture. */
	    x = 1 << (ndims_tot - set->ndims);
	    y = 1 << ndims_tot;
	    for (z = set->setnum; z < y; z += x) {
		merged_goal[i] += goal[z];
	    }
	}

	total_goal += merged_goal[i];
    }

    /* Now scale goals to reflect actual weight of vertices available. */
    for (i = 0; i < nsets; i++) {
	merged_goal[i] = (merged_goal[i] / total_goal) * vwgt_sum;
    }
}
