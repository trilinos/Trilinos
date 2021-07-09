/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "defs.h"
#include "params.h"
#include "structs.h"
#include <math.h>
#include <stdio.h>

/* Confirm that the bipartite match algorithm did the right thing. */
void checkbp(struct vtx_data **graph, /* graph data structure for vertex weights */
             double **         xvecs, /* values to partition */
             int *             sets,  /* set assignments returned */
             double *          dists, /* distances that separate sets */
             int               nvtxs, /* number of vertices */
             int               ndims  /* number of dimensions for division */
)
{
  int    signs[MAXDIMS];     /* signs for coordinates of target points */
  int    sizes[MAXSETS];     /* size of each set */
  int    weights[MAXSETS];   /* size of each set */
  double setval = 0.0;       /* value from assigned set */
  double val, bestval = 0.0; /* value to decide set assignment */
  double tol   = 1.0e-8;     /* numerical tolerance */
  int    error = FALSE;      /* are errors encountered? */
  int    nsets;              /* number of sets */
  int    bestset = -1;       /* set vtx should be assigned to */
  int    i, j, k;            /* loop counter */
  void   checkpnt();

  nsets = 1 << ndims;

  for (i = 0; i < nsets; i++) {
    sizes[i]   = 0;
    weights[i] = 0;
  }

  for (i = 1; i <= nvtxs; i++) {
    /* Is vertex closest to the set it is assigned to? */
    for (j = 0; j < MAXDIMS; j++) {
      signs[j] = -1;
    }
    bestval = 0;
    for (j = 0; j < nsets; j++) {
      val = -dists[j];
      for (k = 1; k <= ndims; k++) {
        val += 2 * signs[k - 1] * xvecs[k][i];
      }
      if (j == sets[i]) {
        setval = val;
      }
      if (j == 0 || val < bestval) {
        bestval = val;
        bestset = j;
      }
      if (signs[0] == 1 && signs[1] == 1) {
        signs[2] *= -1;
      }
      if (signs[0] == 1) {
        signs[1] *= -1;
      }
      signs[0] *= -1;
    }
    if (fabs(setval - bestval) >= tol * (fabs(setval) + fabs(bestval))) {
      printf(" Vtx %d in set %d (%e), but should be in %d (%e)\n", i, sets[i], setval, bestset,
             bestval);
      error = TRUE;
    }
    ++sizes[sets[i]];
    weights[sets[i]] += graph[i]->vwgt;
  }

  printf(" Sizes:");
  for (i = 0; i < nsets; i++) {
    printf(" %d(%d)", sizes[i], weights[i]);
  }
  printf("\n");

  if (error) {
    checkpnt("ERROR in checkbp");
  }
}
