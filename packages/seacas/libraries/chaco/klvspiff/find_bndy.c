/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
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
