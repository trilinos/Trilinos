/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "refine_map.h" // for refine_vdata, refine_edata
#include "structs.h"

double compute_mesh_edata(struct refine_edata *edata,        /* desire data for current edge */
                          struct refine_vdata *vdata,        /* data for all vertices */
                          int                  mesh_dims[3], /* dimensions of processor mesh */
                          struct vtx_data **   comm_graph,   /* communication graph */
                          int *                node2vtx      /* maps mesh nodes to graph vertices */
)
{
  double desire;     /* edge's interest in flipping */
  float  ewgt;       /* edge weight */
  int    vtx1, vtx2; /* vertices on either side of wire */
  int    off;        /* index into vdata */
  int    is_an_edge();

  vtx1 = node2vtx[edata->node1];
  vtx2 = node2vtx[edata->node2];

  off = edata->dim * mesh_dims[0] * mesh_dims[1] * mesh_dims[2];

  desire = (vdata[off + vtx1].above - vdata[off + vtx1].below - vdata[off + vtx1].same) +
           (vdata[off + vtx2].below - vdata[off + vtx2].above - vdata[off + vtx2].same);

  /* Subtract off potential doubly counted edge. */
  if (is_an_edge(comm_graph[vtx1], vtx2, &ewgt)) {
    desire -= 2 * ewgt;
  }

  return (desire);
}
