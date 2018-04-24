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
#include "structs.h" // for vtx_data
#include <stdio.h>   // for fprintf, NULL, fclose, etc

/* Print out subgraph in matrix-market format. */
void mm_out(struct vtx_data **graph,       /* graph data structure */
            int               nvtxs,       /* number of vtxs in graph */
            int               using_ewgts, /* Are edges weighted? */
            char *            tag,         /* message to include */
            char *            file_name    /* output file name if not null */
)
{
  FILE *file;   /* output file */
  int   nedges; /* number of edges in graph */
  int   i, j;   /* loop counter */
  int   DIAG;

  DIAG = TRUE;

  if (file_name != NULL) {
    file = fopen(file_name, "w");
  }
  else {
    file = stdout;
  }

  /* Determine all the appropriate parameters. */
  nedges = 0;
  for (i = 1; i <= nvtxs; i++) {
    nedges += graph[i]->nedges - 1;
  }
  if (DIAG) {
    nedges += nvtxs;
  }

  if (tag != NULL) {
    fprintf(file, "%% graph_out: %s\n", tag);
  }
  fprintf(file, " %d %d %d\n", nvtxs, nvtxs, nedges);
  for (i = 1; i <= nvtxs; i++) {
    if (DIAG) {
      if (!using_ewgts) {
        fprintf(file, "%d %d\n", i, i);
      }
      else {
        fprintf(file, "%d %d %.9f\n", i, i, 1.0);
      }
    }
    for (j = 1; j < graph[i]->nedges; j++) {
      if (!using_ewgts) {
        fprintf(file, "%d %d\n", i, graph[i]->edges[j]);
      }
      else {
        fprintf(file, "%d %d %.9f\n", i, graph[i]->edges[j], graph[i]->ewgts[j]);
      }
    }
  }

  if (file_name != NULL) {
    fclose(file);
  }
}
