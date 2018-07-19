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

#include "defs.h"    // for TRUE
#include "smalloc.h" // for sfree, smalloc
#include "structs.h"

void inertial1d(struct vtx_data **graph,        /* graph data structure */
                int               nvtxs,        /* number of vtxs in graph */
                int               cube_or_mesh, /* 0 => hypercube, d => d-dimensional mesh */
                int               nsets,        /* number of sets to divide into */
                float *           x,            /* x coordinates of vertices */
                int *             sets,         /* set each vertex gets assigned to */
                double *          goal,         /* desired set sizes */
                int               using_vwgts   /* are vertex weights being used? */
                )
{
  extern double median_time; /* time to find medians */
  double *      value;       /* values passed to median routine */
  double        time;        /* timing variables */
  int *         space;       /* space required by median routine */
  int           i;           /* loop counter */
  void          rec_median_1();
  double        seconds();

  value = smalloc((nvtxs + 1) * sizeof(double));

  /* Copy values into double precision array. */
  for (i = 1; i <= nvtxs; i++) {
    value[i] = x[i];
  }

  /* Now find the median value and partition based upon it. */
  space = smalloc(nvtxs * sizeof(int));

  time = seconds();
  rec_median_1(graph, value, nvtxs, space, cube_or_mesh, nsets, goal, using_vwgts, sets, TRUE);
  median_time += seconds() - time;

  sfree(space);
  sfree(value);
}
