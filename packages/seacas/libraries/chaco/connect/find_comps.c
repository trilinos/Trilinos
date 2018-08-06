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

#include "smalloc.h" // for smalloc
#include "structs.h" // for edgeslist, vtx_data
#include <stdio.h>   // for NULL

static int bfsearch();

/* Breadth first search algorithm to find & mark connected components. */
int find_comps(struct vtx_data **graph,  /* graph data structure */
               int               nvtxs,  /* number of vertices in graph */
               int *             mark,   /* space for nvtxs+1 ints */
               int *             vtxlist /* space for nvtxs ints */
)
{
  int    root;   /* vertex to start the dfs */
  int    count;  /* number of vertices seen so far */
  int    ncomps; /* number of components found */
  int    i;      /* loop counter */
  double drandom();

  for (i = 1; i <= nvtxs; i++) {
    mark[i] = -1;
  }
  count  = 0;
  ncomps = 0;
  root   = nvtxs * drandom() + 1;

  bfsearch(graph, root, &count, mark, vtxlist, ncomps);

  while (count != nvtxs) { /* Are there any remaining vertices? */
    /* Find starting vtx for next BFS. */
    root = nvtxs * drandom() + 1;
    while (mark[root] >= 0) {
      root++;
      if (root > nvtxs) {
        root = 1;
      }
    }
    /* Add new edge to list needed for connectivity. */
    ncomps++;
    bfsearch(graph, root, &count, mark, vtxlist, ncomps);
  }
  return (ncomps + 1);
}

/* Breadth first search algorithm to find & mark connected components. */
/* Returns list of edges to connect them together. */
int find_edges(struct vtx_data ** graph,   /* graph data structure */
               int                nvtxs,   /* number of vertices in graph */
               int *              mark,    /* space for nvtxs+1 ints */
               int *              vtxlist, /* space for nvtxs ints */
               struct edgeslist **edges    /* list of edges connecting graph */
)
{
  struct edgeslist *newedge; /* space to add new edge */
  int               root;    /* vertex to start the dfs */
  int               last;    /* last vertex seen in BFS */
  int               count;   /* number of vertices seen so far */
  int               nadded;  /* number of edges needed to be added */
  int               i;       /* loop counter */
  double            drandom();

  for (i = 1; i <= nvtxs; i++) {
    mark[i] = -1;
  }
  count  = 0;
  nadded = 0;
  *edges = NULL;
  root   = nvtxs * drandom() + 1;

  last = bfsearch(graph, root, &count, mark, vtxlist, nadded);

  while (count != nvtxs) { /* Are there any remaining vertices? */
    /* Find starting vtx for next BFS. */
    root = nvtxs * drandom() + 1;
    while (mark[root] >= 0) {
      root++;
      if (root > nvtxs) {
        root = 1;
      }
    }
    /* Add new edge to list needed for connectivity. */
    newedge       = smalloc(sizeof(struct edgeslist));
    newedge->next = *edges;
    newedge->vtx1 = last;
    newedge->vtx2 = root;
    *edges        = newedge;
    nadded++;
    last = bfsearch(graph, root, &count, mark, vtxlist, nadded);
  }
  return (nadded);
}

/* BFS to find connected component */
static int bfsearch(struct vtx_data **graph,   /* graph data structure */
                    int               root,    /* start vertex for DFS */
                    int *             count,   /* number of vertices in component */
                    int *             mark,    /* has vtx been seen? */
                    int *             vtxlist, /* space for storing vtxs to search */
                    int               comp_num /* current component number */
)
{
  int *iptr;           /* loops through neighbor list */
  int  vtxbeg, vtxend; /* beginning and end of vertices in vtxlist */
  int  vtx;            /* vertex being processed */
  int  neighbor;       /* neighbor of vertex */
  int  i;              /* loop counter */

  vtxbeg = vtxend = 1;
  mark[root]      = comp_num;
  vtxlist[0]      = root;

  /* Copy root's neighbors to vtxlist, incrementing count */
  iptr = graph[root]->edges;
  for (i = graph[root]->nedges - 1; i; i--) {
    neighbor          = *(++iptr);
    vtxlist[vtxend++] = neighbor;
    mark[neighbor]    = comp_num;
  }

  while (vtxbeg < vtxend) {
    vtx = vtxlist[vtxbeg++];
    /* Loop through neighbors, copying to vtxlist if unmarked. */
    iptr = graph[vtx]->edges;
    for (i = graph[vtx]->nedges - 1; i; i--) {
      neighbor = *(++iptr);
      if (mark[neighbor] != comp_num) {
        mark[neighbor]    = comp_num;
        vtxlist[vtxend++] = neighbor;
      }
    }
  }
  *count += vtxend;
  return (vtxlist[vtxend - 1]);
}
