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

/* Use same ratios as original goals, but adjust based on set sizes. */
/* Note: mesh stuff only works for division into two sets. */

void make_subgoal(double *goal,         /* goals for sets */
                  double *subgoal,      /* goals for subset of sets */
                  int     nsets,        /* number of subsets in one partition */
                  int     cube_or_mesh, /* 0=> hypercube, d=> d-dimensional mesh */
                  int     nsets_tot,    /* total number of sets to divide into */
                  int     mesh_dims[3], /* shape of mesh */
                  int     set,          /* which set am I in? */
                  double  sub_vwgt_sum  /* sum of subgraph vertex weights */
)
{
  double tweight;        /* total weight among all subgoals */
  double ratio;          /* scaling factor */
  int    index;          /* x, y or z location of a processor */
  int    sub_nsets;      /* largest number of processors in submesh */
  int    xstart, xwidth; /* parameters describing submesh */
  int    i, j, x;        /* loop counters */

  if (!cube_or_mesh) { /* First do hypercube case. */
    tweight = 0;
    for (j = 0, i = set; i < nsets_tot; i += nsets, j++) {
      subgoal[j] = goal[i];
      tweight += goal[i];
    }
    sub_nsets = nsets_tot / nsets;
  }

  else {
    if (set == 0) {
      xstart = 0;
      xwidth = mesh_dims[0] - mesh_dims[0] / 2;
    }
    else {
      xwidth = mesh_dims[0] / 2;
      xstart = mesh_dims[0] - mesh_dims[0] / 2;
    }
    i       = 0;
    tweight = 0;
    index   = xstart;
    for (x = xstart; x < xstart + xwidth; x++) {
      subgoal[i] = goal[index++];
      tweight += subgoal[i++];
    }
    sub_nsets = xwidth;
  }

  ratio = sub_vwgt_sum / tweight;
  for (i = 0; i < sub_nsets; i++) {
    subgoal[i] *= ratio;
  }
}
