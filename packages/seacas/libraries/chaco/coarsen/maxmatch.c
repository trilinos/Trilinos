/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "structs.h"
#include <stdio.h> // for printf

/* Find a maximal matching in a graph using one of several algorithms. */

int maxmatch(struct vtx_data **graph,       /* array of vtx data for graph */
             int               nvtxs,       /* number of vertices in graph */
             int               nedges,      /* number of edges in graph */
             int *             mflag,       /* flag indicating vtx selected or not */
             int               using_ewgts, /* are edge weights being used? */
             int               igeom,       /* geometric dimensionality */
             float **          coords       /* coordinates for each vertex */
)
{
  extern int DEBUG_COARSEN; /* debug output for coarsening? */
  extern int MATCH_TYPE;    /* which matching routine to use */
  int        nmerged = 0;   /* number of matching edges found */
  int        maxmatch1(), maxmatch2(), maxmatch3(), maxmatch4();
  int        maxmatch5(), maxmatch9();

  if (MATCH_TYPE == 1 || (MATCH_TYPE == 5 && coords == NULL)) {
    /* Dumb, fast routine. */
    nmerged = maxmatch1(graph, nvtxs, mflag, using_ewgts);
  }

  else if (MATCH_TYPE == 2) { /* More random but somewhat slower. */
    nmerged = maxmatch2(graph, nvtxs, mflag, using_ewgts);
  }

  else if (MATCH_TYPE == 3) { /* Much more random but slower still. */
    nmerged = maxmatch3(graph, nvtxs, mflag, using_ewgts);
  }

  else if (MATCH_TYPE == 4) { /* Truly random but very slow. */
    nmerged = maxmatch4(graph, nvtxs, nedges, mflag, using_ewgts);
  }
  else if (MATCH_TYPE == 5 && coords != NULL) { /* Geometric nearness. */
    nmerged = maxmatch5(graph, nvtxs, mflag, igeom, coords);
  }

  else if (MATCH_TYPE == 9) { /* Minimum degree of merged vertex */
    nmerged = maxmatch9(graph, nvtxs, mflag, using_ewgts);
  }

  if (DEBUG_COARSEN > 0) {
    printf("Number of matching edges = %d\n", nmerged);
  }

  return (nmerged);
}
