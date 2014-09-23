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

#include <stdio.h>                      // for fprintf, printf, NULL, etc
#include "smalloc.h"                    // for smalloc
#include "structs.h"                    // for vtx_data

/* Make a bipartite graph from vertex separator and neighbors. */

void 
make_bpgraph (
    struct vtx_data **graph,	/* list of graph info for each vertex */
    int *sets,			/* local partitioning of vtxs */
    int *bndy_list,		/* list of vertices on boundary (0 ends) */
    int sep_size,		/* length of bndy_list */
    int set_match,		/* side to match against */
    int **ppointers,		/* start/stop of adjacency lists */
    int **pindices,		/* adjacency list for each vertex */
    int **pvweight,		/* weight of each vertex */
    int **ploc2glob,		/* maps bp number to full graph */
    int *pnleft,		/* number of nodes in left half */
    int *pnright,		/* number of nodes in right half */
    int using_vwgts		/* are vertices weighted? */
)
{
    int      *loc2glob=NULL;	/* maps bp number to full graph */
    int      *pointers=NULL;	/* start/stop of adjacency lists */
    int      *indices=NULL;	/* adjacency list for each vertex */
    int      *vwgts=NULL;	/* saves vwgts so they can be overwritten */
    int       nleft, nright;	/* # vertices in halves of bipartite graph */
    int       nedges;		/* number of edges in bipartite graph */
    int       vtx;		/* vertex in graph */
    int       neighbor;		/* neighbor of vertex */
    int       i, j, k;		/* loop counters */

    /* First count everything that needs to be counted. */
    nleft = sep_size;
    nright = 0;
    nedges = 0;
    for (i = 0; i < sep_size; i++) {
	vtx = bndy_list[i];
	for (j = 1; j < graph[vtx]->nedges; j++) {
	    neighbor = graph[vtx]->edges[j];
	    if (sets[neighbor] == set_match) {
		++nedges;
		if (graph[neighbor]->edges[0] > 0) {	/* Not yet seen */
		    ++nright;
		    /* Flag him as seen already. */
		    graph[neighbor]->edges[0] = -1;
		}
	    }
	}
    }

    pointers = smalloc((nleft + nright + 1) * sizeof(int));
    indices = smalloc((2 * nedges + 1) * sizeof(int));

    /* Now set up data structures to make construction easier */
    loc2glob = smalloc((nleft + nright) * sizeof(int));

    if (!using_vwgts) {
        for (i = 0; i < nleft; i++) {
	    vtx = bndy_list[i];
	    loc2glob[i] = vtx;
	    graph[vtx]->edges[0] = i;
        }

        k = nleft;
        for (i = 0; i < nleft; i++) {
	    vtx = bndy_list[i];
	    for (j = 1; j < graph[vtx]->nedges; j++) {
	        neighbor = graph[vtx]->edges[j];
	        if (sets[neighbor] == set_match) {
		    if (graph[neighbor]->edges[0] == -1) {
		        loc2glob[k] = neighbor;
		        /* Reflag him as seen already with glob2loc value. */
		        graph[neighbor]->edges[0] = k;
		        k++;
		    }
		}
	    }
	}
    }
    else {
        vwgts = smalloc((nleft + nright) * sizeof(int));

        for (i = 0; i < nleft; i++) {
	    vtx = bndy_list[i];
	    loc2glob[i] = vtx;
	    vwgts[i] = graph[vtx]->vwgt;
	    /* Use edges[0] as a seen flag and as a glob2loc value. */
	    graph[vtx]->edges[0] = i;
        }

        k = nleft;
        for (i = 0; i < nleft; i++) {
	    vtx = bndy_list[i];
	    for (j = 1; j < graph[vtx]->nedges; j++) {
	        neighbor = graph[vtx]->edges[j];
	        if (sets[neighbor] == set_match) {
		    if (graph[neighbor]->edges[0] == -1) {	/* First occurence. */
		        loc2glob[k] = neighbor;
			vwgts[k] = graph[neighbor]->vwgt;
			/* Use edges[0] as a seen flag and as a glob2loc value. */
		        graph[neighbor]->edges[0] = k;
		        k++;
		    }
		}
	    }
	}
    }

    /* I can now construct graph directly */
    nedges = 0;
    pointers[0] = 0;
    for (i = 0; i < nleft; i++) {
	vtx = loc2glob[i];
	for (j = 1; j < graph[vtx]->nedges; j++) {
	    neighbor = graph[vtx]->edges[j];
	    if (sets[neighbor] == set_match) {
		indices[nedges++] = graph[neighbor]->edges[0];
	    }
	}
	pointers[i+1] = nedges;
    }

    for (i = nleft; i < nleft + nright; i++) {
	vtx = loc2glob[i];
	for (j = 1; j < graph[vtx]->nedges; j++) {
	    neighbor = graph[vtx]->edges[j];
	    if (sets[neighbor] == 2) {
		indices[nedges++] = graph[neighbor]->edges[0];
	    }
	}
	pointers[i+1] = nedges;
    }

    /* Now restore the edges[0] values. */
    for (i = 0; i < nleft + nright; i++) {
	graph[loc2glob[i]]->edges[0] = loc2glob[i];
    }

/*
check_bpgraph(nleft, nright, pointers, indices);
*/

    if (using_vwgts) {
	*pvweight = vwgts;
    }
    else {
	*pvweight = NULL;
    }

    *ploc2glob = loc2glob;
    *ppointers = pointers;
    *pindices = indices;
    *pnleft = nleft;
    *pnright = nright;

}


void 
check_bpgraph (int n_left, int n_right, int *pointers, int *indices)
{
    int i, j, k, neighbor;

    for (i = 0; i < n_left; i++) {
       for (j = pointers[i]; j < pointers[i+1]; j++) {
	   neighbor = indices[j];
	   if (neighbor < n_left || neighbor >= n_left + n_right) {
	       printf("Bad edge (%d, %d)\n", i, neighbor);
	   }
	   /* Check for counter-edge */
	   for (k = pointers[neighbor]; k < pointers[neighbor+1]; k++) {
	       if (indices[k] == i) break;
	   }
	   if (k == pointers[neighbor+1]) {
	       printf("Flip edge (%d, %d) not found\n", k, i);
	   }
        }
    }

    for (i = n_left; i < n_left + n_right; i++) {
       for (j = pointers[i]; j < pointers[i+1]; j++) {
	   neighbor = indices[j];
	   if (neighbor < 0 || neighbor >= n_left) {
	       printf("Bad edge (%d, %d)\n", i, neighbor);
	   }
	   /* Check for counter-edge */
	   for (k = pointers[neighbor]; k < pointers[neighbor+1]; k++) {
	       if (indices[k] == i) break;
	   }
	   if (k == pointers[neighbor+1]) {
	       printf("Flip edge (%d, %d) not found\n", k, i);
	   }
        }
    }
}


void 
print_bpgraph (int nleft, int nright, int *pointers, int *indices, int *vwgts)
{
   FILE *file;
   int i, j, nedges, nvtxs;

   nvtxs = nleft + nright;
   nedges = (pointers[nvtxs] - pointers[0]) / 2;

   file = fopen("BPGRAPH", "w");

   fprintf(file, "%d %d\n", nvtxs, nedges);

   for (i = 0; i < nvtxs; i++) {
       if (vwgts != NULL) fprintf(file, "%d     ", vwgts[i]);
       for (j = pointers[i]; j < pointers[i+1]; j++) {
	   fprintf(file, "%d ", indices[j]);
       }
       fprintf(file, "\n");
    }
    fclose(file);
}
