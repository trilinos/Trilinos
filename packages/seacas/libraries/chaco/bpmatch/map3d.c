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

static void free3d(double *vals[8][MAXSETS], int *indices[8][MAXSETS]);

void map3d(struct vtx_data **graph,   /* graph data structure */
           double          **xvecs,   /* vectors to partition */
           int               nvtxs,   /* number of vertices */
           int              *sets,    /* set each vertex gets assigned to */
           double           *goal,    /* desired set sizes */
           int               vwgt_max /* largest vertex weight */
)
{
  extern int DEBUG_BPMATCH;        /* debug flag for bipartite matching */
  extern int N_VTX_MOVES;          /* number of vertices moved between sets */
  extern int N_VTX_CHECKS;         /* number of vertex moves considered */
  double    *vals[8][MAXSETS];     /* values in sorted lists */
  double     dist[8];              /* trial separation point */
  double     size[8];              /* sizes of each set being modified */
  int       *indices[8][MAXSETS];  /* indices sorting lists */
  int        startvtx[8][MAXSETS]; /* indices defining separation */
  int        nsection = 3;         /* number of xvectors */
  int        nsets    = 8;         /* number of sets being divided into */

  N_VTX_CHECKS = N_VTX_MOVES = 0;

  /* Generate all the lists of values that need to be sorted. */
  genvals3d(xvecs, vals, nvtxs);

  /* Sort the lists of values. */
  sorts3d(vals, indices, nvtxs);

  /* Now initialize distances using median values, and assign to sets. */
  inits3d(graph, xvecs, vals, indices, nvtxs, dist, startvtx, size, sets);

  /* Determine largest and smallest allowed sizes for each set. */
  /* (For now all are the same, but easily changed.) */

  if (DEBUG_BPMATCH > 1) {
    printf(" Calling check before movevtxs\n");
    checkbp(graph, xvecs, sets, dist, nvtxs, nsection);
  }

  movevtxs(graph, nvtxs, nsets, dist, indices, vals, startvtx, sets, size, goal, vwgt_max);

  if (DEBUG_BPMATCH > 0) {
    printf(" N_VTX_CHECKS = %d, N_VTX_MOVES = %d\n", N_VTX_CHECKS, N_VTX_MOVES);
    checkbp(graph, xvecs, sets, dist, nvtxs, nsection);
  }

  free3d(vals, indices);
}

static void free3d(double *vals[8][MAXSETS], int *indices[8][MAXSETS])
{

  sfree(vals[0][1]);
  sfree(vals[0][2]);
  sfree(vals[0][4]);
  sfree(vals[0][3]);
  sfree(vals[1][2]);
  sfree(vals[0][5]);
  sfree(vals[1][4]);
  sfree(vals[0][6]);
  sfree(vals[2][4]);
  sfree(vals[0][7]);
  sfree(vals[1][6]);
  sfree(vals[2][5]);
  sfree(vals[3][4]);

  sfree(indices[0][1]);
  sfree(indices[0][2]);
  sfree(indices[0][4]);
  sfree(indices[0][3]);
  sfree(indices[1][2]);
  sfree(indices[0][5]);
  sfree(indices[1][4]);
  sfree(indices[0][6]);
  sfree(indices[2][4]);
  sfree(indices[0][7]);
  sfree(indices[1][6]);
  sfree(indices[2][5]);
  sfree(indices[3][4]);
}
