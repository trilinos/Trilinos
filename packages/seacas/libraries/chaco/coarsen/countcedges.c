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

#include "structs.h" // for vtx_data

void countcedges(
    /* Count edges in coarsened graph and construct start array. */
    struct vtx_data **graph,    /* array of vtx data for graph */
    int               nvtxs,    /* number of vertices in graph */
    int *             start,    /* start of edgevals list for each vertex */
    int *             seenflag, /* flags indicating vtxs already counted */
    int *             mflag,    /* flag indicating vtx matched or not */
    int *             v2cv,     /* mapping from fine to coarse vertices */
    int *             pcnedges  /* number of edges in coarsened graph */
)
{
  int *jptr;       /* loops through edge list */
  int  cnedges;    /* twice number of edges in coarsened graph */
  int  neighbor;   /* neighboring vertex */
  int  cneighbor;  /* neighboring vertex in coarse graph */
  int  nneighbors; /* number of neighboring vertices */
  int  newi;       /* loops over vtxs in coarsened graph */
  int  i, j;       /* loop counters */

  /* Note: seenflag values already set to zero. */

  cnedges  = 0;
  newi     = 1;
  start[1] = 0;
  for (i = 1; i <= nvtxs; i++) {
    if (mflag[i] == 0 || mflag[i] > i) {
      nneighbors = 0;
      jptr       = graph[i]->edges;
      for (j = graph[i]->nedges - 1; j; j--) {
        /* Has edge already been added? */
        neighbor = *(++jptr);
        if (neighbor != mflag[i]) {
          cneighbor = v2cv[neighbor];

          if (seenflag[cneighbor] != i) {
            nneighbors++;
            seenflag[cneighbor] = i;
          }
        }
      }

      if (mflag[i] > i) { /* Take care of matching vertex. */
        jptr = graph[mflag[i]]->edges;
        for (j = graph[mflag[i]]->nedges - 1; j; j--) {
          neighbor = *(++jptr);
          if (neighbor != i) {
            cneighbor = v2cv[neighbor];

            if (seenflag[cneighbor] != i) {
              nneighbors++;
              seenflag[cneighbor] = i;
            }
          }
        }
      }

      start[newi + 1] = start[newi] + nneighbors;
      newi++;
      cnedges += nneighbors;
    }
  }
  *pcnedges = cnedges / 2;
}
