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
#include "structs.h" // for vtx_data, edgeslist, flists, etc
#include <stdio.h>   // for NULL

void add_edges(struct vtx_data **graph,      /* graph data structure */
               struct edgeslist *new_edges,  /* list of edges connecting graph */
               struct ilists **  old_edges,  /* edges data overwritten for connecting */
               struct flists **  old_ewgts,  /* weights of edges overwritten */
               int               using_ewgts /* are edge weights being used? */
)
{
  struct ilists *   save_list;  /* space to save old edge list */
  struct flists *   save_ewgts; /* space to save old edge weights */
  struct edgeslist *edges;      /* loops through new edges */
  float *           new_ewgts;  /* new edge weights */
  int *             new_list;   /* new edge list */
  int               nedges;     /* number of edges a vertex has */
  int               vtx, vtx2;  /* two vertices in edge to be added */
  int               i, j;       /* loop counter */

  *old_edges = NULL;
  *old_ewgts = NULL;
  edges      = new_edges;
  while (edges != NULL) {
    for (j = 0; j < 2; j++) {
      if (j == 0) {
        vtx  = edges->vtx1;
        vtx2 = edges->vtx2;
      }
      else {
        vtx  = edges->vtx2;
        vtx2 = edges->vtx1;
      }

      /* Copy old edge list to new edge list. */
      nedges   = graph[vtx]->nedges;
      new_list = smalloc((nedges + 1) * sizeof(int));
      for (i = 0; i < nedges; i++) {
        new_list[i] = graph[vtx]->edges[i];
      }
      new_list[nedges] = vtx2;

      /* Save old edges. */
      save_list       = smalloc(sizeof(struct ilists));
      save_list->list = graph[vtx]->edges;

      /* Add new list at FRONT of linked list to facilitate uncoarsening. */
      save_list->next = *old_edges;
      *old_edges      = save_list;

      /* Now modify graph to have new edges list. */
      graph[vtx]->nedges++;
      graph[vtx]->edges = new_list;

      /* If using edge weights, I have to modify those too. */
      if (using_ewgts) {
        new_ewgts = smalloc((nedges + 1) * sizeof(float));
        for (i = 1; i < nedges; i++) {
          new_ewgts[i] = graph[vtx]->ewgts[i];
        }
        new_ewgts[nedges] = 1;
        new_ewgts[0]      = graph[vtx]->ewgts[0] - new_ewgts[nedges];

        /* Save old edge weights. */
        save_ewgts       = smalloc(sizeof(struct flists));
        save_ewgts->list = graph[vtx]->ewgts;

        save_ewgts->next = *old_ewgts;
        *old_ewgts       = save_ewgts;

        /* Finally, modify graph to have new edge weights. */
        graph[vtx]->ewgts = new_ewgts;
      }
    }
    edges = edges->next;
  }
}
