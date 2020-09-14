/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "structs.h"
#include <stdio.h> // for FILE

/* Print metrics of partition quality. */

void countup(struct vtx_data **graph,        /* graph data structure */
             int               nvtxs,        /* number of vtxs in graph */
             int *             assignment,   /* set number of each vtx (length nvtxs+1) */
             int               ndims,        /* number of cuts at each level */
             int               architecture, /* what's the target parallel machine? */
             int               ndims_tot,    /* total number of hypercube dimensions */
             int               mesh_dims[3], /* extent of mesh in each dimension */
             int               print_lev,    /* level of output */
             FILE *            outfile,      /* output file if not NULL */
             int               using_ewgts   /* are edge weights being used? */
)
{
  extern int VERTEX_SEPARATOR; /* vertex instead of edge separator? */
  extern int VERTEX_COVER;     /* make/improve vtx separator via matching? */
  void       countup_cube(), countup_mesh(), countup_vtx_sep();

  if (VERTEX_SEPARATOR || VERTEX_COVER) {
    countup_vtx_sep(graph, nvtxs, assignment);
  }
  else {
    if (architecture == 0) {
      countup_cube(graph, nvtxs, assignment, ndims, ndims_tot, print_lev, outfile, using_ewgts);
    }

    else if (architecture > 0) {
      countup_mesh(graph, nvtxs, assignment, mesh_dims, print_lev, outfile, using_ewgts);
    }
  }
}
