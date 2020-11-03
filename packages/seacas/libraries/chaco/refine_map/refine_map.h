#ifndef CHACO_REFINE_MAP_REFINE_MAP_H
#define CHACO_REFINE_MAP_REFINE_MAP_H

/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
struct refine_vdata
{
  float above; /* sum of edge weights pulling me higher */
  float below; /* sum of edge weights pulling me lower */
  float same;  /* sum of edge weights keeping me here */
};

struct refine_edata
{
  short                node1, node2; /* nodes in mesh connected by this edge */
  short                dim;          /* which dimension of mesh does wire span? */
  float                swap_desire;  /* reduction in hops if edge is flipped */
  struct refine_edata *prev;         /* pointer to previous guy in list */
  struct refine_edata *next;         /* pointer to next guy in list */
};

#endif
