/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "refine_map.h" // for refine_edata

void init_cube_edata(struct refine_edata *edata, /* desire data for current edge */
                     int                  node1, /* processor incident to current wire */
                     int                  dim,   /* direction of wire */
                     int                  mask   /* bit set in wire dimension */
)
{

  edata->node1 = node1;
  edata->node2 = node1 ^ mask;
  edata->dim   = dim;
}
