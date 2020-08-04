/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "defs.h"    // for TRUE
#include "smalloc.h" // for sfree, smalloc_ret
#include "structs.h"
#include <stdio.h> // for NULL

/* Given a partition, refine the mapping in a locally greedy fashion. */

void refine_map(struct vtx_data **graph,        /* graph data structure */
                int               nvtxs,        /* number of vertices in graph */
                int               using_ewgts,  /* are edge weights being used? */
                int *             assign,       /* current assignment */
                int               cube_or_mesh, /* 0 => hypercube, d => d-dimensional mesh */
                int               ndims_tot,    /* if hypercube, number of dimensions */
                int               mesh_dims[3]  /* if mesh, dimensions of mesh */
)
{
  struct vtx_data **comm_graph;       /* graph for communication requirements */
  int               nsets_tot = 0;    /* total number of sets */
  int *             vtx2node  = NULL; /* mapping of comm_graph vtxs to processors */
  int *             node2vtx  = NULL; /* mapping of sets to comm_graph vtxs */
  double            maxdesire = 0.0;  /* largest possible desire to flip an edge */
  int               error     = 0;    /* out of space? */
  int               i;                /* loop counter */

  double find_maxdeg();
  void   free_graph(), strout();
  int    make_comm_graph(), refine_mesh(), refine_cube();

  if (cube_or_mesh == 0) {
    nsets_tot = 1 << ndims_tot;
  }
  else if (cube_or_mesh == 1) {
    nsets_tot = mesh_dims[0];
  }
  else if (cube_or_mesh == 2) {
    nsets_tot = mesh_dims[0] * mesh_dims[1];
  }
  else if (cube_or_mesh == 3) {
    nsets_tot = mesh_dims[0] * mesh_dims[1] * mesh_dims[2];
  }

  node2vtx = vtx2node = NULL;

  /* Construct the weighted quotient graph representing communication. */
  error = make_comm_graph(&comm_graph, graph, nvtxs, using_ewgts, assign, nsets_tot);

  if (!error) {
    maxdesire = 2 * find_maxdeg(comm_graph, nsets_tot, TRUE, (float *)NULL);

    vtx2node = smalloc_ret((nsets_tot + 1) * sizeof(int));
    node2vtx = smalloc_ret(nsets_tot * sizeof(int));
    if (node2vtx == NULL || vtx2node == NULL) {
      error = 1;
      goto skip;
    }

    for (i = 1; i <= nsets_tot; i++) {
      vtx2node[i]     = i - 1;
      node2vtx[i - 1] = i;
    }

    if (cube_or_mesh > 0) {
      error = refine_mesh(comm_graph, cube_or_mesh, mesh_dims, maxdesire, vtx2node, node2vtx);
    }

    else if (cube_or_mesh == 0) {
      error = refine_cube(comm_graph, ndims_tot, maxdesire, vtx2node, node2vtx);
    }

    if (!error) {
      for (i = 1; i <= nvtxs; i++) {
        assign[i] = vtx2node[assign[i] + 1];
      }
    }
  }

skip:

  if (error) {
    strout("\nWARNING: No space to refine mapping to processors.");
    strout("         NO MAPPING REFINEMENT PERFORMED.\n");
  }

  sfree(node2vtx);
  sfree(vtx2node);
  free_graph(comm_graph);
}
