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
