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

#include "smalloc.h" // for smalloc_ret
#include "structs.h" // for vtx_data
#include <stdio.h>   // for NULL, fprintf, printf, FILE

/* Change from a FORTRAN graph style to our graph data structure. */

int reformat(int *              start,     /* start of edge list for each vertex */
             int *              adjacency, /* edge list data */
             int                nvtxs,     /* number of vertices in graph */
             int *              pnedges,   /* ptr to number of edges in graph */
             int *              vwgts,     /* weights for all vertices */
             float *            ewgts,     /* weights for all edges */
             struct vtx_data ***pgraph     /* ptr to array of vtx data for graph */
)
{
  extern FILE *     Output_File;      /* output file or null */
  struct vtx_data **graph     = NULL; /* array of vtx data for graph */
  struct vtx_data * links     = NULL; /* space for data for all vtxs */
  int *             edges     = NULL; /* space for all adjacency lists */
  float *           eweights  = NULL; /* space for all edge weights */
  int *             eptr      = NULL; /* steps through adjacency list */
  int *             eptr_save = NULL; /* saved index into adjacency list */
  float *           wptr      = NULL; /* steps through edge weights list */
  int               self_edge;        /* number of self loops detected */
  int               size;             /* length of all edge lists */
  double            sum;              /* sum of edge weights for a vtx */
  int               using_ewgts;      /* are edge weights being used? */
  int               using_vwgts;      /* are vertex weights being used? */
  int               i, j;             /* loop counters */

  using_ewgts = (ewgts != NULL);
  using_vwgts = (vwgts != NULL);

  graph   = smalloc_ret((nvtxs + 1) * sizeof(struct vtx_data *));
  *pgraph = graph;
  if (graph == NULL) {
    return (1);
  }

  graph[1] = NULL;

  /* Set up all the basic data structure for the vertices. */
  /* Replace many small mallocs by a few large ones. */
  links = smalloc_ret((nvtxs) * sizeof(struct vtx_data));
  if (links == NULL) {
    return (1);
  }

  for (i = 1; i <= nvtxs; i++) {
    graph[i] = links++;
  }

  graph[1]->edges = NULL;
  graph[1]->ewgts = NULL;

  /* Now fill in all the data fields. */
  if (start != NULL) {
    *pnedges = start[nvtxs] / 2;
  }
  else {
    *pnedges = 0;
  }
  size  = 2 * (*pnedges) + nvtxs;
  edges = smalloc_ret(size * sizeof(int));
  if (edges == NULL) {
    return (1);
  }

  if (using_ewgts) {
    eweights = smalloc_ret(size * sizeof(float));
    if (eweights == NULL) {
      return (1);
    }
  }

  if (start != NULL) {
    eptr = adjacency + start[0];
    wptr = ewgts;
  }
  self_edge = 0;

  for (i = 1; i <= nvtxs; i++) {
    if (using_vwgts) {
      graph[i]->vwgt = *(vwgts++);
    }
    else {
      graph[i]->vwgt = 1;
    }
    if (start != NULL) {
      size = start[i] - start[i - 1];
    }
    else {
      size = 0;
    }
    graph[i]->nedges = size + 1;
    graph[i]->edges  = edges;
    *edges++         = i;
    eptr_save        = eptr;
    for (j = size; j; j--) {
      if (*eptr != i) {
        *edges++ = *eptr++;
      }
      else { /* Self edge, skip it. */
        if (!self_edge) {
          printf("WARNING: Self edge (%d,%d) being ignored\n", i, i);
          if (Output_File != NULL) {
            fprintf(Output_File, "WARNING: Self edge (%d,%d) being ignored\n", i, i);
          }
        }
        ++self_edge;
        eptr++;
        --(graph[i]->nedges);
        --(*pnedges);
      }
    }
    if (using_ewgts) {
      graph[i]->ewgts = eweights;
      eweights++;
      sum = 0;
      for (j = size; j; j--) {
        if (*eptr_save++ != i) {
          sum += *wptr;
          *eweights++ = *wptr++;
        }
        else {
          wptr++;
        }
      }
      graph[i]->ewgts[0] = -sum;
    }
    else {
      graph[i]->ewgts = NULL;
    }
  }
  if (self_edge > 1) {
    printf("WARNING: %d self edges were detected and ignored\n", self_edge);
    if (Output_File != NULL) {
      fprintf(Output_File, "WARNING: %d self edges were detected and ignored\n", self_edge);
    }
  }

  return (0);
}
