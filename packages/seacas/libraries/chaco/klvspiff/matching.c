/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

/* Despite the formal disclaimer above, this routine is modified from
   code provided by Ed Rothberg at SGI. */

#include "smalloc.h"
#include <stdio.h>

#define TRUE 1
#define FALSE 0

/*
The following code takes a bipartite graph as input (with
n_left+n_right nodes, where node 'i' is adjacent to nodes
'indices[pointers[i]]' through 'indices[pointers[i+1]-1]', all
0-based) and returns a minimum size edge cover (in sep_nodes[]).

*/

static void bpmatching(), reachability(), augment();
static int  touch(), touch2();

void bpcover(int  n_left,   /* number of vertices on left side */
             int  n_right,  /* number of vertices on right side */
             int *pointers, /* start/stop of adjacency lists */
             int *indices,  /* adjacency list for each vertex */
             int *sep_size, /* returned size of separator */
             int *sep_nodes /* list of separator nodes */
)
{
  extern int DEBUG_COVER; /* controls debugging output in this routine */
  int *      matching;    /* array to encode matching */
  int *      touched;     /* flags for each vertex */
  int        i;           /* loop counter */

  void confirm_match();

  if (DEBUG_COVER) {
    printf("-> Entering bpcover, nleft = %d, nright = %d, 2*nedges = %d\n", n_left, n_right,
           pointers[n_left + n_right] - pointers[0]);
  }

  matching = smalloc((n_left + n_right) * sizeof(int));
  touched  = smalloc((n_left + n_right) * sizeof(int));

  bpmatching(n_left, n_right, pointers, indices, matching, touched);

  reachability(n_left, n_right, pointers, indices, matching, touched);

  /* Separator includes untouched nodes on left, touched on right. */
  /* Left separator nodes if unconnected to unmatched left node via */
  /* augmenting path, right separator nodes otherwise. */

  *sep_size = 0;
  for (i = 0; i < n_left; i++) {
    if (!touched[i]) {
      sep_nodes[(*sep_size)++] = i;
    }
  }
  for (i = n_left; i < n_left + n_right; i++) {
    if (touched[i]) {
      sep_nodes[(*sep_size)++] = i;
    }
  }

  sep_nodes[(*sep_size)] = 0;

  if (DEBUG_COVER) {
    confirm_match(n_left, n_right, pointers, indices, matching, *sep_size, sep_nodes);
  }

  sfree(touched);
  sfree(matching);
}

static void bpmatching(int  n_left,   /* number of vertices on left side */
                       int  n_right,  /* number of vertices on right side */
                       int *pointers, /* start/stop of adjacency lists */
                       int *indices,  /* adjacency list for each vertex */
                       int *matching, /* array to encode matching */
                       int *touched   /* flags for each vertex */
)
{
  int *seen; /* space for list of encountered vertices */
  int  i, j; /* loop counters */

  /* First mark all the vertices as unmatched & untouched. */
  for (i = 0; i < n_left + n_right; i++) {
    matching[i] = -1;
    touched[i]  = FALSE;
  }

  /* Now generate a fast, greedy matching to start. */
  for (i = n_left; i < n_left + n_right; i++) {
    for (j = pointers[i]; j < pointers[i + 1]; j++) {
      if (matching[indices[j]] == -1) {
        /* Node not already matched. */
        matching[i]          = indices[j];
        matching[indices[j]] = i;
        break;
      }
    }
  }

  /* Now try to enlarge it via augmenting paths. */

  seen = smalloc((n_left + n_right) * sizeof(int));

  /* Look for an augmenting path. */
  for (i = 0; i < n_left; i++) {
    if (matching[i] == -1) {
      augment(i, pointers, indices, matching, touched, seen);
    }
  }

  sfree(seen);
}

static void augment(int  node,     /* start node in augmenting path */
                    int *pointers, /* start/stop of adjacency lists */
                    int *indices,  /* adjacency list for each vertex */
                    int *matching, /* array to encode matching */
                    int *touched,  /* flags for each vertex */
                    int *seen      /* keeps list of vertices encountered */
)
{
  int nseen;    /* number of vertices encountered */
  int enlarged; /* was matching enlarged? */
  int i;        /* loop counter */

  /* Look for augmenting path in graph. */

  nseen    = 0;
  enlarged = touch(node, pointers, indices, matching, touched, seen, &nseen);

  if (enlarged) { /* Found an augmenting path! */
    /* Free all the vertices encountered in search. */
    /* Otherwise, they can't be involved in augmentation, */
    /* so leave them touched. */
    for (i = 0; i < nseen; i++) {
      touched[*seen++] = FALSE;
    }
  }
}

/* Mark everybody in my alternating path tree, and recursively update */
/* matching if augmenting path found. */
static int touch(int node, int *pointers, /* start/stop of adjacency lists */
                 int *indices,            /* adjacency list for each vertex */
                 int *matching,           /* array to encode matching */
                 int *touched,            /* flags for each vertex */
                 int *seen,               /* list of vertices encountered */
                 int *nseen               /* number of vertices encountered */
)
{
  int neighbor; /* neighbor of a vertex */
  int result;   /* return node number (or -1) */
  int j;        /* loop counter */

  touched[node]    = TRUE;
  seen[(*nseen)++] = node;

  for (j = pointers[node]; j < pointers[node + 1]; j++) {
    neighbor = indices[j];
    if (!touched[neighbor]) {
      touched[neighbor] = TRUE;
      seen[(*nseen)++]  = neighbor;
      if (matching[neighbor] == -1) { /* Found augmenting path! */
        matching[neighbor] = node;
        matching[node]     = neighbor;
        return (TRUE);
      }

      result = touch(matching[neighbor], pointers, indices, matching, touched, seen, nseen);
      if (result) {
        matching[neighbor] = node;
        matching[node]     = neighbor;
        return (TRUE);
      }
    }
  }
  return (FALSE);
}

static void reachability(int  n_left,   /* number of vertices on left side */
                         int  n_right,  /* number of vertices on right side */
                         int *pointers, /* start/stop of adjacency lists */
                         int *indices,  /* adjacency list for each vertex */
                         int *matching, /* array to encode matching */
                         int *touched   /* flags for each vertex */
)
{
  int i; /* loop counter */

  /* Initialize all the vertices to be untouched */
  for (i = 0; i < n_left + n_right; i++) {
    touched[i] = 0;
  }

  for (i = 0; i < n_left; i++) {
    if (!touched[i] && matching[i] == -1) {
      touch2(i, pointers, indices, matching, touched);
    }
  }
}

/* Mark everybody in my alternating path tree, and return vertex at */
/* end of augmenting path if found. */
static int touch2(int node, int *pointers, /* start/stop of adjacency lists */
                  int *indices,            /* adjacency list for each vertex */
                  int *matching,           /* array to encode matching */
                  int *touched             /* flags for each vertex */
)
{
  int neighbor; /* neighbor of a vertex */
  int result;   /* return node number (or -1) */
  int j;        /* loop counter */

  touched[node] = TRUE;
  for (j = pointers[node]; j < pointers[node + 1]; j++) {
    neighbor = indices[j];
    if (!touched[neighbor]) {
      touched[neighbor] = TRUE;
      if (matching[neighbor] == -1) {
        return (TRUE);
      }

      result = touch2(matching[neighbor], pointers, indices, matching, touched);
      if (result) {
        return (TRUE);
      }
    }
  }
  return (FALSE);
}

void confirm_match(int  n_left,   /* number of vertices on left side */
                   int  n_right,  /* number of vertices on right side */
                   int *pointers, /* start/stop of adjacency lists */
                   int *indices,  /* adjacency list for each vertex */
                   int *matching, /* array to encode matching */
                   int  sep_size, /* returned size of separator */
                   int *sep_nodes /* list of separator nodes */
)
{
  int *marked;
  int  neighbor;
  int  i, j; /* loop counter */
  int  match_size();

  marked = smalloc((n_left + n_right) * sizeof(int));

  for (i = 0; i < n_left + n_right; i++) {
    marked[i] = FALSE;
  }

  for (i = 0; i < sep_size; i++) {
    marked[sep_nodes[i]] = TRUE;
  }

  for (i = 0; i < n_left; i++) {
    if (!marked[i]) {
      for (j = pointers[i]; j < pointers[i + 1]; j++) {
        neighbor = indices[j];
        if (!marked[neighbor]) {
          printf("Edge (%d, %d) not covered\n", i, neighbor);
        }
      }
    }
  }

  sfree(marked);

  i = match_size(matching, n_left);
  if (sep_size != i) {
    printf("ERROR: sep_size = %d, but match_size = %d\n", sep_size, i);
  }

  for (i = 0; i < n_left + n_right; i++) {
    if (matching[i] != -1 && matching[matching[i]] != i) {
      printf("ERROR: matching[%d] = %d, but matching[%d] = %d\n", i, matching[i], matching[i],
             matching[matching[i]]);
    }
  }
}

int match_size(int *matching, int nleft)
{
  int i, nmatch;

  nmatch = 0;
  for (i = 0; i < nleft; i++) {
    if (matching[i] != -1) {
      ++nmatch;
    }
  }
  return (nmatch);
}
