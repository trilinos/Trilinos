/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "smalloc.h" // for smalloc
#include "structs.h" // for vtx_data
#include <stdio.h>   // for fprintf, printf, NULL, etc

/* Make a bipartite graph from vertex separator and neighbors. */

void make_bpgraph(struct vtx_data **graph,      /* list of graph info for each vertex */
                  int              *sets,       /* local partitioning of vtxs */
                  int              *bndy_list,  /* list of vertices on boundary (0 ends) */
                  int               sep_size,   /* length of bndy_list */
                  int               set_match,  /* side to match against */
                  int             **ppointers,  /* start/stop of adjacency lists */
                  int             **pindices,   /* adjacency list for each vertex */
                  int             **pvweight,   /* weight of each vertex */
                  int             **ploc2glob,  /* maps bp number to full graph */
                  int              *pnleft,     /* number of nodes in left half */
                  int              *pnright,    /* number of nodes in right half */
                  int               using_vwgts /* are vertices weighted? */
)
{
  int *loc2glob = NULL; /* maps bp number to full graph */
  int *pointers = NULL; /* start/stop of adjacency lists */
  int *indices  = NULL; /* adjacency list for each vertex */
  int *vwgts    = NULL; /* saves vwgts so they can be overwritten */
  int  nleft, nright;   /* # vertices in halves of bipartite graph */
  int  nedges;          /* number of edges in bipartite graph */
  int  vtx;             /* vertex in graph */
  int  neighbor;        /* neighbor of vertex */
  int  i, j, k;         /* loop counters */

  /* First count everything that needs to be counted. */
  nleft  = sep_size;
  nright = 0;
  nedges = 0;
  for (i = 0; i < sep_size; i++) {
    vtx = bndy_list[i];
    for (j = 1; j < graph[vtx]->nedges; j++) {
      neighbor = graph[vtx]->edges[j];
      if (sets[neighbor] == set_match) {
        ++nedges;
        if (graph[neighbor]->edges[0] > 0) { /* Not yet seen */
          ++nright;
          /* Flag him as seen already. */
          graph[neighbor]->edges[0] = -1;
        }
      }
    }
  }

  pointers = smalloc((nleft + nright + 1) * sizeof(int));
  indices  = smalloc((2 * nedges + 1) * sizeof(int));

  /* Now set up data structures to make construction easier */
  loc2glob = smalloc((nleft + nright) * sizeof(int));

  if (!using_vwgts) {
    for (i = 0; i < nleft; i++) {
      vtx                  = bndy_list[i];
      loc2glob[i]          = vtx;
      graph[vtx]->edges[0] = i;
    }

    k = nleft;
    for (i = 0; i < nleft; i++) {
      vtx = bndy_list[i];
      for (j = 1; j < graph[vtx]->nedges; j++) {
        neighbor = graph[vtx]->edges[j];
        if (sets[neighbor] == set_match) {
          if (graph[neighbor]->edges[0] == -1) {
            loc2glob[k] = neighbor;
            /* Reflag him as seen already with glob2loc value. */
            graph[neighbor]->edges[0] = k;
            k++;
          }
        }
      }
    }
  }
  else {
    vwgts = smalloc((nleft + nright) * sizeof(int));

    for (i = 0; i < nleft; i++) {
      vtx         = bndy_list[i];
      loc2glob[i] = vtx;
      vwgts[i]    = graph[vtx]->vwgt;
      /* Use edges[0] as a seen flag and as a glob2loc value. */
      graph[vtx]->edges[0] = i;
    }

    k = nleft;
    for (i = 0; i < nleft; i++) {
      vtx = bndy_list[i];
      for (j = 1; j < graph[vtx]->nedges; j++) {
        neighbor = graph[vtx]->edges[j];
        if (sets[neighbor] == set_match) {
          if (graph[neighbor]->edges[0] == -1) { /* First occurrence. */
            loc2glob[k] = neighbor;
            vwgts[k]    = graph[neighbor]->vwgt;
            /* Use edges[0] as a seen flag and as a glob2loc value. */
            graph[neighbor]->edges[0] = k;
            k++;
          }
        }
      }
    }
  }

  /* I can now construct graph directly */
  nedges      = 0;
  pointers[0] = 0;
  for (i = 0; i < nleft; i++) {
    vtx = loc2glob[i];
    for (j = 1; j < graph[vtx]->nedges; j++) {
      neighbor = graph[vtx]->edges[j];
      if (sets[neighbor] == set_match) {
        indices[nedges++] = graph[neighbor]->edges[0];
      }
    }
    pointers[i + 1] = nedges;
  }

  for (i = nleft; i < nleft + nright; i++) {
    vtx = loc2glob[i];
    for (j = 1; j < graph[vtx]->nedges; j++) {
      neighbor = graph[vtx]->edges[j];
      if (sets[neighbor] == 2) {
        indices[nedges++] = graph[neighbor]->edges[0];
      }
    }
    pointers[i + 1] = nedges;
  }

  /* Now restore the edges[0] values. */
  for (i = 0; i < nleft + nright; i++) {
    graph[loc2glob[i]]->edges[0] = loc2glob[i];
  }

  /*
  check_bpgraph(nleft, nright, pointers, indices);
  */

  if (using_vwgts) {
    *pvweight = vwgts;
  }
  else {
    *pvweight = NULL;
  }

  *ploc2glob = loc2glob;
  *ppointers = pointers;
  *pindices  = indices;
  *pnleft    = nleft;
  *pnright   = nright;
}

void check_bpgraph(int n_left, int n_right, int *pointers, int *indices)
{
  int i, j, k, neighbor;

  for (i = 0; i < n_left; i++) {
    for (j = pointers[i]; j < pointers[i + 1]; j++) {
      neighbor = indices[j];
      if (neighbor < n_left || neighbor >= n_left + n_right) {
        printf("Bad edge (%d, %d)\n", i, neighbor);
      }
      /* Check for counter-edge */
      for (k = pointers[neighbor]; k < pointers[neighbor + 1]; k++) {
        if (indices[k] == i) {
          break;
        }
      }
      if (k == pointers[neighbor + 1]) {
        printf("Flip edge (%d, %d) not found\n", k, i);
      }
    }
  }

  for (i = n_left; i < n_left + n_right; i++) {
    for (j = pointers[i]; j < pointers[i + 1]; j++) {
      neighbor = indices[j];
      if (neighbor < 0 || neighbor >= n_left) {
        printf("Bad edge (%d, %d)\n", i, neighbor);
      }
      /* Check for counter-edge */
      for (k = pointers[neighbor]; k < pointers[neighbor + 1]; k++) {
        if (indices[k] == i) {
          break;
        }
      }
      if (k == pointers[neighbor + 1]) {
        printf("Flip edge (%d, %d) not found\n", k, i);
      }
    }
  }
}

void print_bpgraph(int nleft, int nright, int *pointers, int *indices, int *vwgts)
{
  FILE *file;
  int   i, j, nedges, nvtxs;

  nvtxs  = nleft + nright;
  nedges = (pointers[nvtxs] - pointers[0]) / 2;

  file = fopen("BPGRAPH", "w");

  fprintf(file, "%d %d\n", nvtxs, nedges);

  for (i = 0; i < nvtxs; i++) {
    if (vwgts != NULL) {
      fprintf(file, "%d     ", vwgts[i]);
    }
    for (j = pointers[i]; j < pointers[i + 1]; j++) {
      fprintf(file, "%d ", indices[j]);
    }
    fprintf(file, "\n");
  }
  fclose(file);
}
