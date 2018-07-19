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
#include <stdio.h>   // for printf

void countup_vtx_sep(struct vtx_data **graph, /* list of graph info for each vertex */
                     int               nvtxs, /* number of vertices in graph */
                     int *             sets   /* local partitioning of vtxs */
                     )
{
  int vtx, set; /* vertex and set in graph */
  int sep_size; /* size of the separator */
  int i, j, k;  /* loop counters */

  sep_size = 0;
  j = k = 0;
  for (i = 1; i <= nvtxs; i++) {
    if (sets[i] == 0) {
      j += graph[i]->vwgt;
    }
    if (sets[i] == 1) {
      k += graph[i]->vwgt;
    }
    if (sets[i] == 2) {
      sep_size += graph[i]->vwgt;
    }
  }
  printf("Set sizes = %d/%d, Separator size = %d\n\n", j, k, sep_size);

  /* Now check that it really is a separator. */
  for (i = 1; i <= nvtxs; i++) {
    set = sets[i];
    if (set != 2) {
      for (j = 1; j < graph[i]->nedges; j++) {
        vtx = graph[i]->edges[j];
        if (sets[vtx] != 2 && sets[vtx] != set) {
          printf("Error: %d (set %d) adjacent to %d (set %d)\n", i, set, vtx, sets[vtx]);
        }
      }
    }
  }
}
