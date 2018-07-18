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

#include "defs.h"    // for FALSE, TRUE
#include "structs.h" // for vtx_data

/* Find a maximal matching in a graph using simple greedy algorithm. */
/* Randomly permute vertices, and then have each select an unmatched */
/* neighbor. */

int maxmatch1(struct vtx_data **graph,      /* array of vtx data for graph */
              int               nvtxs,      /* number of vertices in graph */
              int *             mflag,      /* flag indicating vtx selected or not */
              int               using_ewgts /* are edge weights being used? */
              )
{
  extern int HEAVY_MATCH; /* choose heavy edges in matching? */
  float      ewgt_max;    /* largest edge weight seen so far */
  int *      jptr;        /* loops through integer arrays */
  int        vtx;         /* vertex to process next */
  int        neighbor;    /* neighbor of a vertex */
  int        nmerged;     /* number of edges in matching */
  int        matched;     /* is a vertex matched yet? */
  int        jsave;       /* best matching edge found so far */
  int        i, j;        /* loop counters */
  double     drandom();

  /* Initialize mflag array. */
  jptr = mflag;
  for (i = nvtxs; i; i--) {
    *(++jptr) = 0;
  }

  nmerged = 0;

  /* Select random starting point in list of vertices. */
  vtx = 1 + drandom() * nvtxs;

  if (!using_ewgts || !HEAVY_MATCH) { /* Choose first neighbor */
    for (i = nvtxs; i; i--) {
      if (mflag[vtx] == 0) { /* Not already matched. */
        /* Select first free edge. */
        matched = FALSE;
        for (j = 1; !matched && j < graph[vtx]->nedges; j++) {
          neighbor = graph[vtx]->edges[j];
          if (mflag[neighbor] == 0) {
            mflag[vtx]      = neighbor;
            mflag[neighbor] = vtx;
            matched         = TRUE;
            nmerged++;
          }
        }
      }

      if (++vtx > nvtxs) {
        vtx = 1;
      }
    }
  }

  else { /* Choose heavy edge neighbor */
    for (i = nvtxs; i; i--) {
      if (mflag[vtx] == 0) { /* Not already matched. */
        /* Select heaviest free edge. */
        jsave    = 0;
        ewgt_max = 0;
        for (j = 1; j < graph[vtx]->nedges; j++) {
          neighbor = graph[vtx]->edges[j];
          if (mflag[neighbor] == 0 && graph[vtx]->ewgts[j] > ewgt_max) {
            ewgt_max = graph[vtx]->ewgts[j];
            jsave    = j;
          }
        }
        if (jsave > 0) {
          neighbor        = graph[vtx]->edges[jsave];
          mflag[vtx]      = neighbor;
          mflag[neighbor] = vtx;
          nmerged++;
        }
      }

      if (++vtx > nvtxs) {
        vtx = 1;
      }
    }
  }

  return (nmerged);
}
