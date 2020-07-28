/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "refine_map.h" // for refine_vdata

void update_cube_vdata(int    old_side,      /* previous side for moved vertex in moved dimension */
                       int    mask,          /* bit set in current dimension */
                       int    neighbor_node, /* node neighbor vertex assigned to */
                       double ewgt,          /* weight of edge */
                       struct refine_vdata *vdata /* neighbor connected by that edge */
)
{
  int neighbor_side; /* side of cube neighbor is on */

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
