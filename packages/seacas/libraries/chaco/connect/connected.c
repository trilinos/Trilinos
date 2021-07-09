/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "smalloc.h" // for sfree, smalloc
#include "structs.h" // for edgeslist, connect_data, etc
#include <stdio.h>   // for NULL, printf

void make_connected(
    /* Add edges to make graph connected. */
    struct vtx_data **    graph,      /* graph data structure */
    int                   nvtxs,      /* number of vertices in graph */
    int *                 nedges,     /* number of edges in graph */
    int *                 mark,       /* space for nvtxs+1 ints */
    int *                 vtxlist,    /* space for nvtxs ints */
    struct connect_data **cdata,      /* space for connectivity data */
    int                   using_ewgts /* are edges of graph weighted? */
)
{
  struct edgeslist *new_edges; /* list of edges connecting graph */
  struct edgeslist *prev_edge; /* pointer for manipulating edge list */
  struct edgeslist *curr_edge; /* pointer for manipulating edge list */
  struct edgeslist *next_edge; /* pointer for manipulating edge list */
  int               nadded;    /* number of edges being added */
  int               find_edges();
  void              add_edges();

  /* First find edges needed to make graph connected. */
  nadded = find_edges(graph, nvtxs, mark, vtxlist, &new_edges);

  /* Now add these needed edges to graph data structure if needed. */
  if (nadded == 0) {
    *cdata = NULL;
  }
  else {
    *cdata              = smalloc(sizeof(struct connect_data));
    (*cdata)->old_edges = NULL;
    (*cdata)->old_ewgts = NULL;
    add_edges(graph, new_edges, &(*cdata)->old_edges, &(*cdata)->old_ewgts, using_ewgts);
    *nedges += nadded;

    /* Now, reverse the order of the new_edges list for consistency with */
    /* the removal order. */
    curr_edge       = new_edges->next;
    new_edges->next = NULL;
    prev_edge       = new_edges;
    while (curr_edge != NULL) {
      next_edge       = curr_edge->next;
      curr_edge->next = prev_edge;
      prev_edge       = curr_edge;
      curr_edge       = next_edge;
    }
    (*cdata)->new_edges = prev_edge;
  }
}

void make_unconnected(
    /* Restore graph to its pristine state and free space for connectivity. */
    struct vtx_data **    graph,      /* graph data structure */
    int *                 nedges,     /* number of edges in graph */
    struct connect_data **cdata,      /* space for connectivity data */
    int                   using_ewgts /* are edges of graph weighted? */
)
{
  struct ilists *   old_edges = NULL; /* edges overwritten for connecting */
  struct flists *   old_ewgts = NULL; /* weights of edges overwritten */
  struct edgeslist *new_edges;        /* list of edges connecting graph */
  struct ilists *   tempi;            /* used for freeing space */
  struct flists *   tempf;            /* used for freeing space */
  struct edgeslist *tempe;            /* used for freeing edgelist space */
  struct edgeslist *edges;            /* loops through new edges */
  int               vtx;              /* vertex in an added edge */
  int               j;                /* loop counters */

  if (*cdata == NULL) {
    return;
  }

  old_edges = (*cdata)->old_edges;
  old_ewgts = (*cdata)->old_ewgts;
  new_edges = (*cdata)->new_edges;
  sfree(*cdata);
  *cdata = NULL;

  edges = new_edges;
  while (edges != NULL) {
    /* Restore edges and weights to original status. */
    (*nedges)--;
    for (j = 0; j < 2; j++) {
      if (j == 0) {
        vtx = edges->vtx2;
      }
      else {
        vtx = edges->vtx1;
      }

      sfree(graph[vtx]->edges);
      graph[vtx]->edges = old_edges->list;
      graph[vtx]->nedges--;
      tempi     = old_edges;
      old_edges = old_edges->next;
      sfree(tempi);

      if (using_ewgts) {
        sfree(graph[vtx]->ewgts);
        graph[vtx]->ewgts = old_ewgts->list;
        tempf             = old_ewgts;
        old_ewgts         = old_ewgts->next;
        sfree(tempf);
      }
    }
    tempe = edges;
    edges = edges->next;
    sfree(tempe);
  }
}

/* Print out the added edges. */
void print_connected(struct connect_data *cdata /* space for connectivity data */
)
{
  struct edgeslist *edges; /* loops through new edges */

  if (cdata == NULL) {
    printf("No phantom edges\n");
  }
  else {
    printf("Phantom edges: ");
    edges = cdata->new_edges;
    while (edges != NULL) {
      printf("(%d,%d) ", edges->vtx1, edges->vtx2);
      edges = edges->next;
    }
    printf("\n");
  }
}

/* Free the edge list created by find_edges. */
void free_edgeslist(struct edgeslist *edge_list /* list to be freed */
)
{
  struct edgeslist *next_list; /* next guy in list */

  while (edge_list != NULL) {
    next_list = edge_list->next;
    sfree(edge_list);
    edge_list = next_list;
  }
}
