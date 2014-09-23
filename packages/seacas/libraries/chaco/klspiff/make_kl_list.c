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

#include <stdio.h>                      // for NULL
#include "structs.h"                    // for bilist, vtx_data


int 
make_kl_list (
    struct vtx_data **graph,	/* data structure for graph */
    struct bilist *movelist,	/* list of vtxs to be moved */
    struct bilist ****buckets,	/* array of lists for bucket sort */
    struct bilist **listspace,	/* list data structure for each vertex */
    int *sets,			/* processor each vertex is assigned to */
    int nsets,		/* number of sets divided into */
    int *bspace,		/* list of active vertices for bucketsort */
    int **dvals,		/* d-values for each transition */
    int maxdval		/* maximum d-value for a vertex */
)
{
    struct bilist **list;	/* bucket to erase element from */
    struct bilist *vptr;	/* loops through movelist */
    int      *bptr;		/* loops through bspace */
    int      *iptr;		/* loops through edge list */
    int       vtx;		/* vertex that was moved */
    int       neighbor;		/* neighbor of a vertex */
    int       myset;		/* set a vertex is in */
    int       newset;		/* loops through other sets */
    int       list_length;	/* number of values put in bspace */
    int       size;		/* array spacing */
    int       i, l;		/* loop counter */
    void      removebilist();

    /* First push all the moved vertices onto list, so they can be flagged. */
    /* They've already been removed from buckets, so want to avoid them. */
    size = (int) (&(listspace[0][1]) - &(listspace[0][0]));
    vptr = movelist;
    bptr = bspace;
    list_length = 0;
    while (vptr != NULL) {
	vtx = ((int) (vptr - listspace[0])) / size;
	*bptr++ = vtx;
	if (sets[vtx] >= 0)
	    sets[vtx] = -sets[vtx] - 1;
	++list_length;
	vptr = vptr->next;
    }

    /* Now look at all the neighbors of moved vertices. */
    vptr = movelist;
    while (vptr != NULL) {
	vtx = ((int) (vptr - listspace[0])) / size;

	iptr = graph[vtx]->edges;
	for (i = graph[vtx]->nedges - 1; i; i--) {
	    neighbor = *(++iptr);
	    if (sets[neighbor] >= 0) {
		*bptr++ = neighbor;
		++list_length;
		myset = sets[neighbor];
		sets[neighbor] = -sets[neighbor] - 1;

		/* Remove neighbor entry from all his buckets. */
		/* Note: vertices in movelist already removed from buckets. */
		l = 0;
		for (newset = 0; newset < nsets; newset++) {
		    if (newset != myset) {
			list = &buckets[myset][newset][dvals[neighbor][l] + maxdval];
			removebilist(&listspace[l][neighbor], list);
			l++;
		    }
		}
	    }
	}
	vptr = vptr->next;
    }

    /* Now that list is constructed, go reconstruct all the set numbers. */
    bptr = bspace;
    for (i = list_length; i; i--) {
	vtx = *bptr++;
	sets[vtx] = -sets[vtx] - 1;
    }

    return (list_length);
}
