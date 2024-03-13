/*
 * Copyright(C) 1999-2020, 2023 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "defs.h"
#include "structs.h"
#include <stdio.h> // for printf

void assign(struct vtx_data **graph,        /* data structure with vtx weights */
            double          **yvecs,        /* ptr to list of y-vectors (lengths nvtxs+1) */
            int               nvtxs,        /* number of vertices in graph */
            int               ndims,        /* number of vectors for dividing */
            int               cube_or_mesh, /* 0 => hypercube, d => d-dimensional mesh */
            int               nsets,        /* number of sets to divide into */
            double           *wsqrt,        /* sqrt of vertex weights */
            int              *sets,         /* processor assignment for my vtxs */
            int              *active,       /* space for nvtxs integers */
            int               mediantype,   /* which partitioning strategy to use */
            double           *goal,         /* desired set sizes */
            int               vwgt_max      /* largest vertex weight */
)
{
  extern int DEBUG_TRACE;       /* trace execution path of code */
  extern int DEBUG_ASSIGN;      /* turn on debugging in assignment */
  double     theta, phi, gamma; /* angles for optimal rotation */
  double     temp;
  int        using_vwgts; /* are vertex weights active? */

  if (DEBUG_TRACE > 0) {
    printf("<Entering assign, nvtxs = %d, ndims = %d>\n", nvtxs, ndims);
  }

  using_vwgts = (vwgt_max != 1);

  if (ndims == 1) {
    /* Unscale yvecs to get xvecs. */
    y2x(yvecs, ndims, nvtxs, wsqrt);

    mapper(graph, yvecs, nvtxs, active, sets, ndims, cube_or_mesh, nsets, mediantype, goal,
           vwgt_max);
  }

  else if (ndims == 2) {
    theta = opt2d(graph, yvecs, nvtxs, nvtxs);

    rotate2d(yvecs, nvtxs, theta);

    y2x(yvecs, ndims, nvtxs, wsqrt);

    mapper(graph, yvecs, nvtxs, active, sets, ndims, cube_or_mesh, nsets, mediantype, goal,
           vwgt_max);
  }

  else if (ndims == 3) {

    if (DEBUG_ASSIGN > 0) {
      temp = tri_prod(yvecs[1], yvecs[2], yvecs[3], wsqrt, nvtxs);
      printf("Before rotation, 3-way orthogonality = %e\n", temp);
    }

    opt3d(graph, yvecs, nvtxs, nvtxs, wsqrt, &theta, &phi, &gamma, using_vwgts);

    rotate3d(yvecs, nvtxs, theta, phi, gamma);

    if (DEBUG_ASSIGN > 0) {
      temp = tri_prod(yvecs[1], yvecs[2], yvecs[3], wsqrt, nvtxs);
      printf("After rotation (%f,%f,%f), 3-way orthogonality = %e\n", theta, phi, gamma, temp);
    }

    y2x(yvecs, ndims, nvtxs, wsqrt);

    mapper(graph, yvecs, nvtxs, active, sets, ndims, cube_or_mesh, nsets, mediantype, goal,
           vwgt_max);
  }
}
