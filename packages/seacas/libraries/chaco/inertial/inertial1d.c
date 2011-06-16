/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "defs.h"
#include "structs.h"
#include "smalloc.h"


void      inertial1d(graph, nvtxs, cube_or_mesh, nsets, x, sets, goal,
		               using_vwgts)
struct vtx_data **graph;	/* graph data structure */
int       nvtxs;		/* number of vtxs in graph */
int       cube_or_mesh;		/* 0 => hypercube, d => d-dimensional mesh */
int       nsets;		/* number of sets to divide into */
float    *x;			/* x coordinates of vertices */
int    *sets;			/* set each vertex gets assigned to */
double   *goal;			/* desired set sizes */
int       using_vwgts;		/* are vertex weights being used? */
{
    extern double median_time;	/* time to find medians */
    double   *value;		/* values passed to median routine */
    double    time;		/* timing variables */
    int      *space;		/* space required by median routine */
    int       i;		/* loop counter */
    void      rec_median_1();
    double    seconds();


    value = smalloc((nvtxs + 1) * sizeof(double));

    /* Copy values into double precision array. */
    for (i = 1; i <= nvtxs; i++)
	value[i] = x[i];

    /* Now find the median value and partition based upon it. */
    space = smalloc(nvtxs * sizeof(int));

    time = seconds();
    rec_median_1(graph, value, nvtxs, space, cube_or_mesh, nsets, goal,
		 using_vwgts, sets, TRUE);
    median_time += seconds() - time;

    sfree(space);
    sfree(value);
}
