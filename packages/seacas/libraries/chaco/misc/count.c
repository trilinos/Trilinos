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

#include "params.h"  // for MAXSETS
#include "smalloc.h" // for sfree, smalloc
#include "structs.h" // for vtx_data
#include <stdio.h>   // for printf

void count(struct vtx_data **graph, /* graph data structure */
           int               nvtxs, /* number of vtxs in graph */
           int *             sets,  /* processor each vertex is assigned to */
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
