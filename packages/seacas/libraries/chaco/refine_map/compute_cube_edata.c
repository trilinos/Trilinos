/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "refine_map.h" // for refine_vdata, refine_edata
#include "structs.h"

double compute_cube_edata(struct refine_edata *edata,      /* desire data for current edge */
                          struct refine_vdata *vdata,      /* data for all vertices */
                          int                  nsets_tot,  /* total number of processors */
                          struct vtx_data **   comm_graph, /* communication graph */
                          int *                node2vtx    /* maps mesh nodes to graph vertices */
)
{
  double desire;     /* edge's interest in flipping */
  float  ewgt;       /* edge weight */
  int    offset;     /* offset into vdata array */
  int    vtx1, vtx2; /* vertices on either side of wire */
  int    is_an_edge();

  vtx1   = node2vtx[edata->node1];
  vtx2   = node2vtx[edata->node2];
  offset = nsets_tot * edata->dim;

  desire = (vdata[offset + vtx1].above - vdata[offset + vtx1].same) +
           (vdata[offset + vtx2].above - vdata[offset + vtx2].same);

  /* Subtract off potential doubly counted edge. */
  if (is_an_edge(comm_graph[vtx1], vtx2, &ewgt)) {
    desire -= 2 * ewgt;
  }

  return (desire);
}
