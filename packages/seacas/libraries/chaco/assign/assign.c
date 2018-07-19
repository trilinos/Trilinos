/*
 * Copyright (c) 2005-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include "structs.h"
#include <stdio.h> // for printf

void assign(struct vtx_data **graph,        /* data structure with vtx weights */
            double **         yvecs,        /* ptr to list of y-vectors (lengths nvtxs+1) */
            int               nvtxs,        /* number of vertices in graph */
            int               ndims,        /* number of vectors for dividing */
            int               cube_or_mesh, /* 0 => hypercube, d => d-dimensional mesh */
            int               nsets,        /* number of sets to divide into */
            double *          wsqrt,        /* sqrt of vertex weights */
            int *             sets,         /* processor assignment for my vtxs */
            int *             active,       /* space for nvtxs integers */
            int               mediantype,   /* which partitioning strategy to use */
            double *          goal,         /* desired set sizes */
            int               vwgt_max      /* largest vertex weight */
            )
{
  extern int DEBUG_TRACE;       /* trace execution path of code */
  extern int DEBUG_ASSIGN;      /* turn on debugging in assignment */
  double     theta, phi, gamma; /* angles for optimal rotation */
  double     temp;
  int        using_vwgts; /* are vertex weights active? */
  double     tri_prod();
  double     opt2d();
  void       y2x(), mapper(), rotate2d(), opt3d(), rotate3d();

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
