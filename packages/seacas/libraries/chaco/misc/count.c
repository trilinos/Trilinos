/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "params.h"  // for MAXSETS
#include "smalloc.h" // for sfree, smalloc
#include "structs.h" // for vtx_data
#include <stdio.h>   // for printf

void count(struct vtx_data **graph, /* graph data structure */
           int               nvtxs, /* number of vtxs in graph */
           int              *sets,  /* processor each vertex is assigned to */
           int               nsets, /* number of sets partitioned into */
           int (*hops)[MAXSETS],    /* hops metric between sets */
           int dump,                /* flag for extended output */
           int using_ewgts          /* are edge weights being used? */
)
{
  int *nguys;      /* number of vtxs in each set */
  int  ncross;     /* number of outgoing edges */
  int  nhops;      /* number of hops */
  int  neighbor;   /* neighbor of a vertex */
  int  nmax, nmin; /* largest and smallest set sizes */
  int  i, j;       /* loop counters */

  nguys = smalloc(nsets * sizeof(int));

  for (i = 0; i < nsets; i++) {
    nguys[i] = 0;
  }

  ncross = nhops = 0;
  for (i = 1; i <= nvtxs; i++) {
    nguys[sets[i]] += graph[i]->vwgt;

    for (j = 1; j < graph[i]->nedges; j++) {
      neighbor = graph[i]->edges[j];
      if (sets[neighbor] != sets[i]) {
        if (using_ewgts) {
          ncross += graph[i]->ewgts[j];
          nhops += graph[i]->ewgts[j] * hops[sets[i]][sets[neighbor]];
        }
        else {
          ncross++;
          nhops += hops[sets[i]][sets[neighbor]];
        }
      }
    }
  }

  ncross /= 2;
  nhops /= 2;

  nmax = nguys[0];
  nmin = nguys[0];
  for (i = 1; i < nsets; i++) {
    if (nguys[i] > nmax) {
      nmax = nguys[i];
    }
    if (nguys[i] < nmin) {
      nmin = nguys[i];
    }
  }
  printf("In subgraph: Cuts=%d, Hops=%d; Max=%d, Min=%d (nvtxs=%d).\n", ncross, nhops, nmax, nmin,
         nvtxs);

  if (dump) {
    for (i = 0; i < nsets; i++) {
      printf(" Size of %d = %d\n", i, nguys[i]);
    }

    for (i = 0; i < nvtxs; i++) {
      printf("%d\n", sets[i]);
    }
    printf("\n\n");
  }

  sfree(nguys);
}
