/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "smalloc.h" // for sfree, smalloc_ret
#include "structs.h" // for vtx_data
#include <stdio.h>   // for NULL

/* Construct a weighted quotient graph representing the inter-set communication. */
int make_comm_graph(struct vtx_data ***pcomm_graph, /* graph for communication requirements */
                    struct vtx_data ** graph,       /* graph data structure */
                    int                nvtxs,       /* number of vertices in graph */
                    int                using_ewgts, /* are edge weights being used? */
                    int *              assign,      /* current assignment */
                    int                nsets_tot    /* total number of sets */
)
{
  float  ewgt;               /* edge weight in graph */
  int ** edges_list  = NULL; /* lists of edges */
  int ** ewgts_list  = NULL; /* lists of edge weights */
  int *  edges       = NULL; /* edges in communication graph */
  int *  ewgts       = NULL; /* edge weights in communication graph */
  float *float_ewgts = NULL; /* edge weights in floating point */
  int *  adj_sets    = NULL; /* weights connecting sets */
  int *  order       = NULL; /* ordering of vertices by set */
  int *  sizes       = NULL; /* sizes of different sets */
  int *  start       = NULL; /* pointers into adjacency data */
  int *  adjacency   = NULL; /* array with all the edge info */
  int *  eptr        = NULL; /* loops through edges in graph */
  int *  ewptr       = NULL; /* loop through edge weights */
  int    set, set2;          /* sets two vertices belong to */
  int    vertex;             /* vertex in graph */
  int    ncomm_edges;        /* number of edges in communication graph */
  int    error;              /* out of space? */
  int    i, j;               /* loop counters */
  int    reformat();

  error        = 1;
  *pcomm_graph = NULL;

  /* First construct some mappings to ease later manipulations. */
  sizes = smalloc_ret(nsets_tot * sizeof(int));
  if (sizes == NULL) {
    goto skip;
  }

  for (i = 0; i < nsets_tot; i++) {
    sizes[i] = 0;
  }
  for (i = 1; i <= nvtxs; i++) {
    ++(sizes[assign[i]]);
  }

  /* Now make sizes reflect the start index for each set. */
  for (i = 1; i < nsets_tot - 1; i++) {
    sizes[i] += sizes[i - 1];
  }
  for (i = nsets_tot - 1; i; i--) {
    sizes[i] = sizes[i - 1];
  }
  sizes[0] = 0;

  /* Now construct list of all vertices in set 0, all in set 1, etc. */
  order = smalloc_ret(nvtxs * sizeof(int));
  if (order == NULL) {
    goto skip;
  }
  for (i = 1; i <= nvtxs; i++) {
    set               = assign[i];
    order[sizes[set]] = i;
    ++sizes[set];
  }

  /* For each set, find total weight to all neighbors. */
  adj_sets   = smalloc_ret(nsets_tot * sizeof(int));
  edges_list = smalloc_ret(nsets_tot * sizeof(int *));
  ewgts_list = smalloc_ret(nsets_tot * sizeof(int *));
  start      = smalloc_ret((nsets_tot + 1) * sizeof(int));
  if (adj_sets == NULL || edges_list == NULL || ewgts_list == NULL || start == NULL) {
    goto skip;
  }

  start[0]    = 0;
  ewgt        = 1;
  ncomm_edges = 0;

  for (set = 0; set < nsets_tot; set++) {
    edges_list[set] = NULL;
    ewgts_list[set] = NULL;
  }

  for (set = 0; set < nsets_tot; set++) {
    for (i = 0; i < nsets_tot; i++) {
      adj_sets[i] = 0;
    }
    for (i = (set ? sizes[set - 1] : 0); i < sizes[set]; i++) {
      vertex = order[i];
      for (j = 1; j < graph[vertex]->nedges; j++) {
        set2 = assign[graph[vertex]->edges[j]];
        if (set2 != set) {
          if (using_ewgts) {
            ewgt = graph[vertex]->ewgts[j];
          }
          adj_sets[set2] += ewgt;
        }
      }
    }

    /* Now save adj_sets data to later construct graph. */
    j = 0;
    for (i = 0; i < nsets_tot; i++) {
      if (adj_sets[i]) {
        j++;
      }
    }
    ncomm_edges += j;
    start[set + 1] = ncomm_edges;
    if (j) {
      edges_list[set] = edges = smalloc_ret(j * sizeof(int));
      ewgts_list[set] = ewgts = smalloc_ret(j * sizeof(int));
      if (edges == NULL || ewgts == NULL) {
        goto skip;
      }
    }
    j = 0;
    for (i = 0; i < nsets_tot; i++) {
      if (adj_sets[i]) {
        edges[j] = i + 1;
        ewgts[j] = adj_sets[i];
        j++;
      }
    }
  }

  sfree(adj_sets);
  sfree(order);
  sfree(sizes);
  adj_sets = order = sizes = NULL;

  /* I now need to pack the edge and weight data into single arrays. */
  adjacency   = smalloc_ret((ncomm_edges + 1) * sizeof(int));
  float_ewgts = smalloc_ret((ncomm_edges + 1) * sizeof(float));
  if (adjacency == NULL || float_ewgts == NULL) {
    goto skip;
  }

  for (set = 0; set < nsets_tot; set++) {
    j     = start[set];
    eptr  = edges_list[set];
    ewptr = ewgts_list[set];
    for (i = start[set]; i < start[set + 1]; i++) {
      adjacency[i]   = eptr[i - j];
      float_ewgts[i] = ewptr[i - j];
    }
    if (start[set] != start[set + 1]) {
      sfree(edges_list[set]);
      sfree(ewgts_list[set]);
    }
  }
  sfree(edges_list);
  sfree(ewgts_list);
  edges_list = ewgts_list = NULL;

  error =
      reformat(start, adjacency, nsets_tot, &ncomm_edges, (int *)NULL, float_ewgts, pcomm_graph);

skip:
  sfree(adj_sets);
  sfree(order);
  sfree(sizes);
  if (edges_list != NULL) {
    for (set = nsets_tot - 1; set >= 0; set--) {
      if (edges_list[set] != NULL) {
        sfree(edges_list[set]);
      }
    }
    sfree(edges_list);
  }

  if (ewgts_list != NULL) {
    for (set = nsets_tot - 1; set >= 0; set--) {
      if (ewgts_list[set] != NULL) {
        sfree(ewgts_list[set]);
      }
    }
    sfree(ewgts_list);
  }

  sfree(float_ewgts);
  sfree(adjacency);
  sfree(start);

  return (error);
}
