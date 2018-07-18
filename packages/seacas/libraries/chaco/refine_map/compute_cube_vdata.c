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

#include "refine_map.h" // for refine_vdata
#include "structs.h"    // for vtx_data

void compute_cube_vdata(struct refine_vdata *vdata,      /* preference data for a vertex */
                        struct vtx_data **   comm_graph, /* communication graph data structure */
                        int                  vtx,        /* current vertex */
                        int                  mask,    /* bit set in current hypercube dimension */
                        int *                vtx2node /* maps graph vtxs to mesh nodes */
                        )
{
  float same;        /* my preference to stay where I am */
  float change;      /* my preference to change this bit */
  float ewgt;        /* weight of an edge */
  int   neighbor;    /* neighboring vtx in comm_graph */
  int   my_side;     /* which side of hypercube I'm on */
  int   neighb_side; /* which side of hypercube neighbor's on */
  int   j;           /* loop counter */

  my_side = (vtx2node[vtx] & mask);

  change = 0;
  same   = 0;
  for (j = 1; j < comm_graph[vtx]->nedges; j++) {
    neighbor = comm_graph[vtx]->edges[j];
    ewgt     = comm_graph[vtx]->ewgts[j];

    neighb_side = (vtx2node[neighbor] & mask);

    if (neighb_side != my_side) {
      change += ewgt;
    }
    else {
      same += ewgt;
    }
  }
  vdata->same  = same;
  vdata->above = change;
}
