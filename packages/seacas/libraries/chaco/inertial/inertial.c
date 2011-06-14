/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <math.h>
#include <stdio.h>
#include "structs.h"

int       PROJECTION_AXIS = 0;	/* axis to flatten geometry */
				/* => long regions, good for SnRad */

void      inertial(graph, nvtxs, cube_or_mesh, nsets, igeom, coords, sets,
		             goal, using_vwgts)
struct vtx_data **graph;	/* graph data structure */
int       nvtxs;		/* number of vtxs in graph */
int       cube_or_mesh;		/* 0 => hypercube, d => d-dimensional mesh */
int       nsets;		/* number of sets to cut into */
int       igeom;		/* 1, 2 or 3 dimensional geometry? */
float   **coords;		/* x, y and z coordinates of vertices */
int    *sets;			/* set each vertex gets assigned to */
double   *goal;			/* desired set sizes */
int       using_vwgts;		/* are vertex weights being used? */
{
    extern    int DEBUG_TRACE;	/* trace the execution of the code */
    extern    int PROJECTION_AXIS;	/* axis to project out geometry */
    extern double inertial_time;/* time spend in inertial calculations */
    double    time;		/* timing parameter */
    float    *inert_coords[3];	/* coord arrays passed down */
    int       i, j;		/* loop counters */
    double    seconds();
    void      inertial1d(), inertial2d(), inertial3d();

    time = seconds();

    if (DEBUG_TRACE > 0) {
	printf("<Entering inertial, nvtxs = %d>\n", nvtxs);
    }

    if (PROJECTION_AXIS == 0) {
	for (i = 0; i < igeom; i++) {
	    inert_coords[i] = coords[i];
	}
    }

    else {	/* project out an axis to get long regions */
	j = 0;
	for (i = 0; i < igeom; i++) {
	    if (PROJECTION_AXIS != i+1) {
	        inert_coords[j] = coords[i];
		j++;
	    }
	}
	--igeom;
    }


    if (igeom == 1)
	inertial1d(graph, nvtxs, cube_or_mesh, nsets,
		   inert_coords[0], sets, goal, using_vwgts);

    else if (igeom == 2)
	inertial2d(graph, nvtxs, cube_or_mesh, nsets,
		   inert_coords[0], inert_coords[1], sets, goal, using_vwgts);

    else if (igeom == 3)
	inertial3d(graph, nvtxs, cube_or_mesh, nsets,
	       inert_coords[0], inert_coords[1], inert_coords[2],
	       sets, goal, using_vwgts);
    inertial_time += seconds() - time;
}
