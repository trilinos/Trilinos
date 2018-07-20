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

#include "refine_map.h" // for refine_vdata, refine_edata
#include "structs.h"

double compute_cube_edata(struct refine_edata *edata,      /* desire data for current edge */
                          struct refine_vdata *vdata,      /* data for all vertices */
                          int                  nsets_tot,  /* total number of processors */
                          struct vtx_data **   comm_graph, /* communication graph */
                          int *                node2vtx    /* maps mesh nodes to graph vertices */
                          )
{
  double desire;     /* edge's interest in flipping */
  float  ewgt;       /* edge weight */
  int    offset;     /* offset into vdata array */
  int    vtx1, vtx2; /* vertices on either side of wire */
  int    is_an_edge();

  vtx1   = node2vtx[edata->node1];
  vtx2   = node2vtx[edata->node2];
  offset = nsets_tot * edata->dim;

  desire = (vdata[offset + vtx1].above - vdata[offset + vtx1].same) +
           (vdata[offset + vtx2].above - vdata[offset + vtx2].same);

  /* Subtract off potential doubly counted edge. */
  if (is_an_edge(comm_graph[vtx1], vtx2, &ewgt)) {
    desire -= 2 * ewgt;
  }

  return (desire);
}
