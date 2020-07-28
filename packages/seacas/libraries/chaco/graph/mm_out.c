/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
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
  nedges += nvtxs;

  if (tag != NULL) {
    fprintf(file, "%% graph_out: %s\n", tag);
  }
  fprintf(file, " %d %d %d\n", nvtxs, nvtxs, nedges);
  for (i = 1; i <= nvtxs; i++) {
    if (!using_ewgts) {
      fprintf(file, "%d %d\n", i, i);
    }
    else {
      fprintf(file, "%d %d %.9f\n", i, i, 1.0);
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
