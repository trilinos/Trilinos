/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include	<stdio.h>
#include	"structs.h"
#include	"params.h"
#include	"defs.h"
#include "smalloc.h"


/* Construct a graph representing the inter-set communication. */
int       refine_part(graph, nvtxs, using_ewgts, assign, architecture,
    ndims_tot, mesh_dims, goal)
struct vtx_data **graph;	/* graph data structure */
int       nvtxs;		/* number of vertices in graph */
int       using_ewgts;		/* are edge weights being used? */
int    *assign;		/* current assignment */
int       architecture;		/* 0 => hypercube, d => d-dimensional mesh */
int       ndims_tot;		/* if hypercube, number of dimensions */
int       mesh_dims[3];		/* if mesh, size in each direction */
double   *goal;			/* desired set sizes */
{
    extern int TERM_PROP;	/* perform terminal propagation? */
    struct bilist *set_list = NULL;	/* lists of vtxs in each set */
    struct bilist *vtx_elems = NULL;	/* space for all vtxs in set_lists */
    struct bilist *ptr = NULL;	/* loops through set_lists */
    struct ipairs *pairs = NULL;/* ordered list of edges in comm graph */
    double   *comm_vals = NULL;	/* edge wgts of comm graph for sorting */
    float    *term_wgts[2];	/* terminal propagation vector */
    int     hops[MAXSETS][MAXSETS];	/* preference weighting */
    double   *temp = NULL;	/* return argument from srealloc_ret() */
    int      *indices = NULL;	/* sorted order for communication edges */
    int      *space = NULL;	/* space for mergesort */
    int      *sizes = NULL;	/* sizes of the different sets */
    int    *sub_assign = NULL;	/* new assignment for subgraph */
    int    *old_sub_assign = NULL;	/* room for current sub assignment */
    int     **edges_list = NULL;/* lists of comm graph edges */
    int     **ewgts_list = NULL;/* lists of comm graph edge wgts */
    int      *ewgts = NULL;	/* loops through ewgts_list */
    int      *edges = NULL;	/* edges in communication graph */
    int      *adj_sets = NULL;	/* weights connecting sets */
    int      *eptr = NULL;	/* loop through edges and edge weights */
    int      *ewptr = NULL;	/* loop through edges and edge weights */
    int       ewgt;		/* weight of an edge */
    struct vtx_data **subgraph = NULL;	/* subgraph data structure */
    int      *nedges = NULL;	/* space for saving graph data */
    int    *degrees = NULL;	/* # neighbors of vertices */
    int      *glob2loc = NULL;	/* maps full to reduced numbering */
    int      *loc2glob = NULL;	/* maps reduced to full numbering */
    int       nmax;		/* largest subgraph I expect to encounter */
    int       set, set1, set2;	/* sets vertices belong to */
    int       vertex;		/* vertex in graph */
    int       ncomm;		/* # edges in communication graph */
    int       dist=-1;		/* architectural distance between two sets */
    int       nsets_tot=0;	/* total number of processors */
    int       change;		/* did change occur in this pass? */
    int       any_change = FALSE;/* has any change occured? */
    int       error;		/* out of space? */
    int       size;		/* array spacing */
    int       i, j, k;		/* loop counters */
    int       abs(), kl_refine();
    void      mergesort(), strout();

    error = 1;
    term_wgts[1] = NULL;

    if (architecture == 0) {
	nsets_tot = 1 << ndims_tot;
    }
    else if (architecture > 0) {
	nsets_tot = mesh_dims[0] * mesh_dims[1] * mesh_dims[2];
    }

    hops[0][0] = hops[1][1] = 0;
    if (!TERM_PROP) {
        hops[0][1] = hops[1][0] = 1;
    }

    /* Set up convenient data structure for navigating through sets. */
    set_list = smalloc_ret(nsets_tot * sizeof(struct bilist));
    vtx_elems = smalloc_ret((nvtxs + 1) * sizeof(struct bilist));
    sizes = smalloc_ret(nsets_tot * sizeof(int));
    if (set_list == NULL || vtx_elems == NULL || sizes == NULL) goto skip;


    for (i = 0; i < nsets_tot; i++) {
	set_list[i].next = NULL;
	sizes[i] = 0;
    }

    for (i = 1; i <= nvtxs; i++) {
	set = assign[i];
	++sizes[set];
	vtx_elems[i].next = set_list[set].next;
	if (vtx_elems[i].next != NULL) {
	    vtx_elems[i].next->prev = &(vtx_elems[i]);
	}
	vtx_elems[i].prev = &(set_list[set]);
	set_list[set].next = &(vtx_elems[i]);
    }


    /* For each set, find connections to all set neighbors. */
    edges_list = smalloc_ret(nsets_tot * sizeof(int *));
    if (edges_list == NULL) goto skip;
    for (set = 0; set < nsets_tot-1; set++) {
	edges_list[set] = NULL;
    }

    ewgts_list = smalloc_ret(nsets_tot * sizeof(int *));
    if (ewgts_list == NULL) goto skip;
    for (set = 0; set < nsets_tot-1; set++) {
	ewgts_list[set] = NULL;
    }

    nedges = smalloc_ret(nsets_tot * sizeof(int));
    adj_sets = smalloc_ret(nsets_tot * sizeof(int));
    if (nedges == NULL || adj_sets == NULL) goto skip;

    size = (int) (&(vtx_elems[1]) - &(vtx_elems[0]));
    ncomm = 0;
    ewgt = 1;
    nmax = 0;
    for (set = 0; set < nsets_tot-1; set++) {
	if (sizes[set] > nmax) nmax = sizes[set];
        for (i = 0; i < nsets_tot; i++) {
	    adj_sets[i] = 0;
        }
	for (ptr = set_list[set].next; ptr != NULL; ptr = ptr->next) {
	    vertex = ((int) (ptr - vtx_elems)) / size;
	    for (j = 1; j < graph[vertex]->nedges; j++) {
		set2 = assign[graph[vertex]->edges[j]];
		if (using_ewgts)
		    ewgt = graph[vertex]->ewgts[j];
		adj_sets[set2] += ewgt;
	    }
	}

	/* Now save adj_sets data to later construct graph. */
	j = 0;
	for (i = set+1; i < nsets_tot; i++) {
	    if (adj_sets[i] != 0)
		j++;
	}
	nedges[set] = j;
	if (j) {
	    edges_list[set] = edges = smalloc_ret(j * sizeof(int));
	    ewgts_list[set] = ewgts = smalloc_ret(j * sizeof(int));
	    if (edges == NULL || ewgts == NULL) goto skip;
	}
	j = 0;
	for (i = set + 1; i < nsets_tot; i++) {
	    if (adj_sets[i] != 0) {
		edges[j] = i;
		ewgts[j] = adj_sets[i];
		j++;
	    }
	}
	ncomm += j;
    }

    sfree(adj_sets);
    adj_sets = NULL;


    /* Now compact all communication weight information into single */
    /* vector for sorting. */

    pairs = smalloc_ret((ncomm + 1) * sizeof(struct ipairs));
    comm_vals = smalloc_ret((ncomm + 1) * sizeof(double));
    if (pairs == NULL || comm_vals == NULL) goto skip;

    j = 0;
    for (set = 0; set < nsets_tot - 1; set++) {
	eptr = edges_list[set];
	ewptr = ewgts_list[set];
	for (k = 0; k < nedges[set]; k++) {
	    set2 = eptr[k];
	    pairs[j].val1 = set;
	    pairs[j].val2 = set2;
	    comm_vals[j] = ewptr[k];
	    j++;
	}
    }

    sfree(nedges);
    nedges = NULL;

    indices = smalloc_ret((ncomm + 1) * sizeof(int));
    space = smalloc_ret((ncomm + 1) * sizeof(int));
    if (indices == NULL || space == NULL) goto skip;

    mergesort(comm_vals, ncomm, indices, space);
    sfree(space);
    sfree(comm_vals);
    space = NULL;
    comm_vals = NULL;

    for (set = 0; set < nsets_tot - 1; set++) {
	if (edges_list[set] != NULL)
	    sfree(edges_list[set]);
	if (ewgts_list[set] != NULL)
	    sfree(ewgts_list[set]);
    }
    sfree(ewgts_list);
    sfree(edges_list);
    ewgts_list = NULL;
    edges_list = NULL;

    /* 2 for 2 subsets, 20 for safety margin. Should check this at run time. */
    nmax = 2 * nmax + 20;

    subgraph = smalloc_ret((nmax + 1) * sizeof(struct vtx_data *));
    degrees = smalloc_ret((nmax + 1) * sizeof(int));
    glob2loc = smalloc_ret((nvtxs + 1) * sizeof(int));
    loc2glob = smalloc_ret((nmax + 1) * sizeof(int));
    sub_assign = smalloc_ret((nmax + 1) * sizeof(int));
    old_sub_assign = smalloc_ret((nmax + 1) * sizeof(int));

    if (subgraph == NULL || degrees == NULL || glob2loc == NULL ||
        loc2glob == NULL || sub_assign == NULL || old_sub_assign == NULL) {
	goto skip;
    }

    if (TERM_PROP) {
        term_wgts[1] = smalloc_ret((nmax + 1) * sizeof(float));
	if (term_wgts[1] == NULL) goto skip;
    }
    else {
	term_wgts[1] = NULL;
    }

    /* Do large boundaries first to encourage convergence. */
    any_change = FALSE;
    for (i = ncomm - 1; i >= 0; i--) {
	j = indices[i];
	set1 = pairs[j].val1;
	set2 = pairs[j].val2;


	/* Make sure subgraphs aren't too big. */
	if (sizes[set1] + sizes[set2] > nmax) {
	    nmax = sizes[set1] + sizes[set2];

	    temp = srealloc_ret(subgraph,
		(nmax + 1) * sizeof(struct vtx_data *));
	    if (temp == NULL) {
		goto skip;
	    }
	    else {
		subgraph = (struct vtx_data **) temp;
	    }

	    temp = srealloc_ret(degrees,
		(nmax + 1) * sizeof(int));
	    if (temp == NULL) {
		goto skip;
	    }
	    else {
		degrees = (int *) temp;
	    }

	    temp = srealloc_ret(loc2glob,
		(nmax + 1) * sizeof(int));
	    if (temp == NULL) {
		goto skip;
	    }
	    else {
		loc2glob = (int *) temp;
	    }

	    temp = srealloc_ret(sub_assign,
		(nmax + 1) * sizeof(int));
	    if (temp == NULL) {
		goto skip;
	    }
	    else {
		sub_assign = (int *) temp;
	    }

	    temp = srealloc_ret(old_sub_assign,
		(nmax + 1) * sizeof(int));
	    if (temp == NULL) {
		goto skip;
	    }
	    else {
		old_sub_assign = (int *) temp;
	    }

	    if (TERM_PROP) {
	        temp = srealloc_ret(term_wgts[1],
		    (nmax + 1) * sizeof(float));
	        if (temp == NULL) {
		    goto skip;
	        }
	        else {
		    term_wgts[1] = (float *) temp;
	        }
	    }
	}


	if (TERM_PROP) {
	    if (architecture == 0) {
		j = set1 ^ set2;
		dist = 0;
		while (j) {
		    if (j & 1) dist++;
		    j >>= 1;
		}
	    }
	    else if (architecture > 0) {
		dist = abs((set1 % mesh_dims[0]) - (set2 % mesh_dims[0]));
		dist += abs(((set1 / mesh_dims[0]) % mesh_dims[1]) -
		            ((set2 / mesh_dims[0]) % mesh_dims[1]));
		dist += abs((set1 / (mesh_dims[0] * mesh_dims[1])) -
		            (set2 / (mesh_dims[0] * mesh_dims[1])));
	    }
	    hops[0][1] = hops[1][0] = (int) dist;
	}

	change = kl_refine(graph, subgraph, set_list, vtx_elems, assign,
	    set1, set2, glob2loc, loc2glob, sub_assign, old_sub_assign,
	    degrees, using_ewgts, hops, goal, sizes, term_wgts, architecture,
	    mesh_dims);
       
	any_change |= change;
    }
    error = 0;

skip:
    if (error) {
	strout("\nWARNING: No space to refine partition.");
	strout("         NO PARTITION REFINEMENT PERFORMED.\n");
    }

    if (edges_list != NULL) {
        for (set = 0; set < nsets_tot - 1; set++) {
	    if (edges_list[set] != NULL)
		sfree(edges_list[set]);
	}
        sfree(edges_list);
    }

    if (ewgts_list != NULL) {
        for (set = 0; set < nsets_tot - 1; set++) {
	    if (ewgts_list[set] != NULL)
	        sfree(ewgts_list[set]);
	}
        sfree(ewgts_list);
    }

    sfree(space);
    sfree(comm_vals);
    sfree(nedges);
    sfree(adj_sets);
    sfree(term_wgts[1]);
    sfree(old_sub_assign);
    sfree(sub_assign);
    sfree(loc2glob);
    sfree(glob2loc);
    sfree(degrees);
    sfree(subgraph);

    sfree(indices);
    sfree(pairs);

    sfree(sizes);
    sfree(vtx_elems);
    sfree(set_list);

    return(any_change);
}
