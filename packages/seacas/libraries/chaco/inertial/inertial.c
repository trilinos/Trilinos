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

int PROJECTION_AXIS = 0; /* axis to flatten geometry */
                         /* => long regions, good for SnRad */

void inertial(struct vtx_data **graph,        /* graph data structure */
              int               nvtxs,        /* number of vtxs in graph */
              int               cube_or_mesh, /* 0 => hypercube, d => d-dimensional mesh */
              int               nsets,        /* number of sets to cut into */
              int               igeom,        /* 1, 2 or 3 dimensional geometry? */
              float **          coords,       /* x, y and z coordinates of vertices */
              int *             sets,         /* set each vertex gets assigned to */
              double *          goal,         /* desired set sizes */
              int               using_vwgts   /* are vertex weights being used? */
)
{
  extern int    DEBUG_TRACE;     /* trace the execution of the code */
  extern int    PROJECTION_AXIS; /* axis to project out geometry */
  extern double inertial_time;   /* time spend in inertial calculations */
  double        time;            /* timing parameter */
  float *       inert_coords[3]; /* coord arrays passed down */
  int           i, j;            /* loop counters */
  double        seconds();
  void          inertial1d(), inertial2d(), inertial3d();

  time = seconds();

  if (DEBUG_TRACE > 0) {
    printf("<Entering inertial, nvtxs = %d>\n", nvtxs);
  }

  if (PROJECTION_AXIS == 0) {
    for (i = 0; i < igeom; i++) {
      inert_coords[i] = coords[i];
    }
  }

  else { /* project out an axis to get long regions */
    j = 0;
    for (i = 0; i < igeom; i++) {
      if (PROJECTION_AXIS != i + 1) {
        inert_coords[j] = coords[i];
        j++;
      }
    }
    --igeom;
  }

  if (igeom == 1) {
    inertial1d(graph, nvtxs, cube_or_mesh, nsets, inert_coords[0], sets, goal, using_vwgts);
  }
  else if (igeom == 2) {
    inertial2d(graph, nvtxs, cube_or_mesh, nsets, inert_coords[0], inert_coords[1], sets, goal,
               using_vwgts);
  }
  else if (igeom == 3) {
    inertial3d(graph, nvtxs, cube_or_mesh, nsets, inert_coords[0], inert_coords[1], inert_coords[2],
               sets, goal, using_vwgts);
  }
  inertial_time += seconds() - time;
}
