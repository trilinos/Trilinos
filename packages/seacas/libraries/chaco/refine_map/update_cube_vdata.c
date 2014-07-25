/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include "refine_map.h"                 // for refine_vdata


void 
update_cube_vdata (
    int old_side,		/* previous side for moved vertex in moved dimension */
    int mask,			/* bit set in current dimension */
    int neighbor_node,	/* node neighbor vertex assigned to */
    double ewgt,			/* weight of edge */
    struct refine_vdata *vdata	/* neighbor connected by that edge */
)
{
    int neighbor_side;		/* side of cube neighbor is on */

    neighbor_side = (neighbor_node & mask);

    if (neighbor_side == old_side) {
	vdata->above += ewgt;
	vdata->same -= ewgt;
    }
    else {
	vdata->above -= ewgt;
	vdata->same += ewgt;
    }
}
