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
#include "smalloc.h" // for smalloc, sfree
#include "structs.h" // for vtx_data

void mapper(struct vtx_data **graph,        /* data structure with vertex weights */
            double **         xvecs,        /* continuous indicator vectors */
            int               nvtxs,        /* number of vtxs in graph */
            int *             active,       /* space for nmyvals ints */
            int *             sets,         /* returned processor assignment for my vtxs */
            int               ndims,        /* number of dimensions being divided into */
            int               cube_or_mesh, /* 0 => hypercube, d => d-dimensional mesh */
            int               nsets,        /* number of sets to divide into */
            int               mediantype,   /* type of eigenvector partitioning to use */
            double *          goal,         /* desired set sizes */
            int               vwgt_max      /* largest vertex weight */
)
{
  double temp_goal[2]; /* combined set goals if using option 1. */
  double wbelow;       /* weight of vertices with negative values */
  double wabove;       /* weight of vertices with positive values */
  int *  temp_sets;    /* sets vertices get assigned to */
  int    vweight;      /* weight of a vertex */
  int    using_vwgts;  /* are vertex weights being used? */
  int    bits;         /* bits for assigning set numbers */
  int    i, j;         /* loop counters */

  void median(), rec_median_k(), rec_median_1(), median_assign();
  void map2d(), map3d();

  /* NOTE: THIS EXPECTS XVECS, NOT YVECS! */

  using_vwgts = (vwgt_max != 1);

  if (ndims == 1 && mediantype == 1) {
    mediantype = 3; /* simpler call than normal option 1. */
  }

  if (mediantype == 0) { /* Divide at zero instead of median. */
    bits      = 1;
    temp_sets = smalloc((nvtxs + 1) * sizeof(int));
    for (j = 1; j <= nvtxs; j++) {
      sets[j] = 0;
    }

    for (i = 1; i <= ndims; i++) {
      temp_goal[0] = temp_goal[1] = 0;
      for (j = 0; j < (1 << ndims); j++) {
        if (bits & j) {
          temp_goal[1] += goal[j];
        }
        else {
          temp_goal[0] += goal[j];
        }
      }
      bits <<= 1;

      wbelow = wabove = 0;
      vweight         = 1;
      for (j = 1; j <= nvtxs; j++) {
        if (using_vwgts) {
          vweight = graph[j]->vwgt;
        }
        if (xvecs[i][j] < 0) {
          wbelow += vweight;
        }
        else if (xvecs[i][j] > 0) {
          wabove += vweight;
        }
      }

      median_assign(graph, xvecs[i], nvtxs, temp_goal, using_vwgts, temp_sets, wbelow, wabove, 0.0);

      for (j = 1; j <= nvtxs; j++) {
        sets[j] = (sets[j] << 1) + temp_sets[j];
      }
    }
    sfree(temp_sets);
  }

  else if (mediantype == 1) { /* Divide using min-cost assignment. */
    if (ndims == 2) {
      map2d(graph, xvecs, nvtxs, sets, goal, vwgt_max);
    }
    else if (ndims == 3) {
      map3d(graph, xvecs, nvtxs, sets, goal, vwgt_max);
    }
  }

  else if (mediantype == 2) { /* Divide recursively using medians. */
    rec_median_k(graph, xvecs, nvtxs, active, ndims, cube_or_mesh, goal, using_vwgts, sets);
  }

  else if (mediantype == 3) { /* Cut with independent medians => unbalanced. */
    bits      = 1;
    temp_sets = smalloc((nvtxs + 1) * sizeof(int));
    for (j = 1; j <= nvtxs; j++) {
      sets[j] = 0;
    }

    for (i = 1; i <= ndims; i++) {
      temp_goal[0] = temp_goal[1] = 0;
      for (j = 0; j < (1 << ndims); j++) {
        if (bits & j) {
          temp_goal[1] += goal[j];
        }
        else {
          temp_goal[0] += goal[j];
        }
      }
      bits <<= 1;

      median(graph, xvecs[i], nvtxs, active, temp_goal, using_vwgts, temp_sets);
      for (j = 1; j <= nvtxs; j++) {
        sets[j] = (sets[j] << 1) + temp_sets[j];
      }
    }
    sfree(temp_sets);
  }

  if (mediantype == 4) { /* Stripe the domain. */
    rec_median_1(graph, xvecs[1], nvtxs, active, cube_or_mesh, nsets, goal, using_vwgts, sets,
                 TRUE);
  }
}
