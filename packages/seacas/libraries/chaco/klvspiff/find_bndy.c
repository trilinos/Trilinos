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

#include "smalloc.h" // for smalloc, srealloc
#include "structs.h" // for vtx_data
#include <stdio.h>   // for NULL

/* Find vertices on boundary of partition, and change their assignments. */

int find_bndy(struct vtx_data **graph,      /* array of vtx data for graph */
              int               nvtxs,      /* number of vertices in graph */
              int *             assignment, /* processor each vertex gets assigned to */
              int               new_val,    /* assignment value for boundary vtxs */
              int **            pbndy_list  /* returned list, end with zero */
              )
{
  int *bndy_list;   /* returned list, end with zero */
  int *edges;       /* loops through edge list */
  int  list_length; /* returned number of vtxs on boundary */
  int  set, set2;   /* set a vertex is in */
  int  i, j;        /* loop counters */

  bndy_list = smalloc((nvtxs + 1) * sizeof(int));

  list_length = 0;
  for (i = 1; i <= nvtxs; i++) {
    set   = assignment[i];
    edges = graph[i]->edges;
    for (j = graph[i]->nedges - 1; j; j--) {
      set2 = assignment[*(++edges)];
      if (set2 != set) {
        bndy_list[list_length++] = i;
        break;
      }
    }
  }
  bndy_list[list_length] = 0;

  for (i = 0; i < list_length; i++) {
    assignment[bndy_list[i]] = new_val;
  }

  /* Shrink out unnecessary space */
  *pbndy_list = srealloc(bndy_list, (list_length + 1) * sizeof(int));

  return (list_length);
}

/* Find a vertex separator on one side of an edge separator. */

int find_side_bndy(struct vtx_data **graph,      /* array of vtx data for graph */
                   int               nvtxs,      /* number of vertices in graph */
                   int *             assignment, /* processor each vertex gets assigned to */
                   int               side,       /* side to take vertices from */
                   int               new_val,    /* assignment value for boundary vtxs */
                   int **            pbndy_list  /* returned list, end with zero */
                   )

{
  int *edges;       /* loops through edge list */
  int *bndy_list;   /* returned list, end with zero */
  int  list_length; /* returned number of vtxs on boundary */
  int  set, set2;   /* set a vertex is in */
  int  i, j;        /* loop counters */

  if (*pbndy_list != NULL) {
    /* Contains list of all vertices on boundary. */
    bndy_list = *pbndy_list;
    i = list_length = 0;
    while (bndy_list[i] != 0) {
      if (assignment[bndy_list[i]] == side) {
        bndy_list[list_length++] = bndy_list[i];
      }
      ++i;
    }
  }

  else {
    bndy_list = smalloc((nvtxs + 1) * sizeof(int));

    list_length = 0;
    for (i = 1; i <= nvtxs; i++) {
      set = assignment[i];
      if (set == side) {
        edges = graph[i]->edges;
        for (j = graph[i]->nedges - 1; j; j--) {
          set2 = assignment[*(++edges)];
          if (set2 != set) {
            bndy_list[list_length++] = i;
            break;
          }
        }
      }
    }
  }

  bndy_list[list_length] = 0;

  for (i = 0; i < list_length; i++) {
    assignment[bndy_list[i]] = new_val;
  }

  /* Shrink out unnecessary space */
  *pbndy_list = srealloc(bndy_list, (list_length + 1) * sizeof(int));

  return (list_length);
}
