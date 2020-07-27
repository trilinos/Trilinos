/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "smalloc.h"
#include "structs.h"
#include <stdio.h>

/* Free a graph data structure. */

void free_graph(struct vtx_data **graph)
{

  if (graph != NULL) {
    if (graph[1] != NULL) {
      if (graph[1]->ewgts != NULL) {
        sfree(graph[1]->ewgts);
      }
      if (graph[1]->edges != NULL) {
        sfree(graph[1]->edges);
      }
      sfree(graph[1]);
    }
    sfree(graph);
  }
}
