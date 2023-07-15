/*
 * Copyright(C) 1999-2020, 2023 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "params.h"  // for MAXSETS
#include "smalloc.h" // for sfree
#include "structs.h"
#include <stdio.h> // for printf

/* Free the space used in the bpmatch routines. */
static void free2d(double *vals[4][MAXSETS], int *indices[4][MAXSETS])
{
  sfree(vals[0][1]);
  sfree(vals[0][2]);
  sfree(vals[0][3]);
  sfree(vals[1][2]);

  sfree(indices[0][1]);
  sfree(indices[0][2]);
  sfree(indices[0][3]);
  sfree(indices[1][2]);
}

void map2d(struct vtx_data **graph,   /* data structure with vertex weights */
           double          **xvecs,   /* vectors to partition */
           int               nvtxs,   /* number of vertices */
           int              *sets,    /* set each vertex gets assigned to */
           double           *goal,    /* desired set sizes */
           int               vwgt_max /* largest vertex weight */
)
{
  extern int DEBUG_BPMATCH;        /* turn on debugging for bipartite matching */
  extern int N_VTX_MOVES;          /* total number of vertex moves */
  extern int N_VTX_CHECKS;         /* total number of moves contemplated */
  double    *vals[4][MAXSETS];     /* values in sorted lists */
  double     dist[4];              /* trial separation point */
  double     size[4];              /* sizes of each set being modified */
  int       *indices[4][MAXSETS];  /* indices sorting lists */
  int        startvtx[4][MAXSETS]; /* indices defining separation */
  int        nsection = 2;         /* number of xvectors */
  int        nsets    = 4;         /* number of sets being divided into */

  N_VTX_CHECKS = N_VTX_MOVES = 0;

  /* Generate all the lists of values that need to be sorted. */
  genvals2d(xvecs, vals, nvtxs);

  /* Sort the lists of values. */
  sorts2d(vals, indices, nvtxs);

  /* Now initialize dists and assign to sets. */
  inits2d(graph, xvecs, vals, indices, nvtxs, dist, startvtx, size, sets);

  /* Determine the largest and smallest allowed set sizes. */
  /* (For now, assume all sets must be same size, but can easily change.) */

  if (DEBUG_BPMATCH > 1) {
    printf(" Calling check before movevtxs\n");
    checkbp(graph, xvecs, sets, dist, nvtxs, nsection);
  }

  movevtxs(graph, nvtxs, nsets, dist, indices, vals, startvtx, sets, size, goal, vwgt_max);

  if (DEBUG_BPMATCH > 0) {
    printf(" N_VTX_CHECKS = %d, N_VTX_MOVES = %d\n", N_VTX_CHECKS, N_VTX_MOVES);
    checkbp(graph, xvecs, sets, dist, nvtxs, nsection);
  }

  free2d(vals, indices);
}
