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

#include "smalloc.h" // for smalloc, sfree, srealloc
#include "structs.h" // for vtx_data
#include <stdio.h>   // for NULL, printf

static void makecv2v();

void makefgraph(struct vtx_data ** graph,       /* array of vtx data for graph */
                int                nvtxs,       /* number of vertices in graph */
                int                nedges,      /* number of edges in graph */
                struct vtx_data ***pcgraph,     /* coarsened version of graph */
                int                cnvtxs,      /* number of vtxs in coarsened graph */
                int *              pcnedges,    /* number of edges in coarsened graph */
                int *              v2cv,        /* mapping from vtxs to coarsened vtxs */
                int                using_ewgts, /* are edge weights being used? */
                int                igeom,       /* dimensions of geometric data */
                float **           coords,      /* coordinates for vertices */
                float **           ccoords      /* coordinates for coarsened vertices */
)
{
  extern double     make_cgraph_time;
  extern int        DEBUG_COARSEN;    /* debug flag for coarsening output */
  extern int        COARSEN_VWGTS;    /* turn off vertex weights in coarse graph? */
  extern int        COARSEN_EWGTS;    /* turn off edge weights in coarse graph? */
  struct vtx_data **cgraph    = NULL; /* coarsened version of graph */
  struct vtx_data * links     = NULL; /* space for all the vertex data */
  struct vtx_data **gptr      = NULL; /* loops through cgraph */
  struct vtx_data * cgptr     = NULL; /* loops through cgraph */
  int *             iptr      = NULL; /* loops through integer arrays */
  int *             seenflag  = NULL; /* flags for vtxs already put in edge list */
  int *             sptr      = NULL; /* loops through seenflags */
  int *             cv2v_vals = NULL; /* vtxs corresponding to each cvtx */
  int *             cv2v_ptrs = NULL; /* indices into cv2v_vals */
  float *           eweights  = NULL; /* space for edge weights in coarsened graph */
  float *           ewptr     = NULL; /* loops through eweights */
  float *           fptr      = NULL; /* loops through eweights */
  float             ewgt;             /* edge weight */
  double            ewgt_sum;         /* sum of edge weights */
  double            time;             /* timing parameters */
  int               nseen;            /* number of edges of coarse graph seen so far */
  int               vtx;              /* vertex in original graph */
  int               cvtx;             /* vertex in coarse graph */
  int               cnedges;          /* twice number of edges in coarsened graph */
  int               neighbor;         /* neighboring vertex */
  int               size;             /* space needed for coarsened graph */
  int *             edges = NULL;     /* space for edges in coarsened graph */
  int *             eptr  = NULL;     /* loops through edges data structure */
  int               cneighbor;        /* neighboring vertex number in coarsened graph */
  int               i, j;             /* loop counters */
  double            seconds();
  void              makeccoords();

  /* Compute the number of vertices and edges in the coarsened graph, */
  /* and construct start pointers into coarsened edge array. */
  time = seconds();

  /* Construct mapping from original graph vtxs to coarsened graph vtxs. */
  cv2v_vals = smalloc(nvtxs * sizeof(int));
  cv2v_ptrs = smalloc((cnvtxs + 2) * sizeof(int));
  makecv2v(nvtxs, cnvtxs, v2cv, cv2v_vals, cv2v_ptrs);

  /* Compute an upper bound on the number of coarse graph edges. */
  cnedges = nedges - (nvtxs - cnvtxs);

  /* Now allocate space for the new graph.  Overallocate and realloc later. */
  *pcgraph = cgraph = smalloc((cnvtxs + 1) * sizeof(struct vtx_data *));
  links             = smalloc(cnvtxs * sizeof(struct vtx_data));

  size  = 2 * cnedges + cnvtxs;
  edges = smalloc(size * sizeof(int));
  if (COARSEN_EWGTS) {
    ewptr = eweights = smalloc(size * sizeof(float));
  }

  /* Zero all the seen flags. */
  seenflag = smalloc((cnvtxs + 1) * sizeof(int));
  sptr     = seenflag;
  for (i = cnvtxs; i; i--) {
    *(++sptr) = 0;
  }

  /* Use the renumbering to fill in the edge lists for the new graph. */
  cnedges = 0;
  eptr    = edges;
  ewgt    = 1;

  sptr = cv2v_vals;
  for (cvtx = 1; cvtx <= cnvtxs; cvtx++) {
    nseen = 1;

    cgptr = cgraph[cvtx] = links++;

    if (COARSEN_VWGTS) {
      cgptr->vwgt = 0;
    }
    else {
      cgptr->vwgt = 1;
    }

    eptr[0]      = cvtx;
    cgptr->edges = eptr;
    if (COARSEN_EWGTS) {
      cgptr->ewgts = ewptr;
    }
    else {
      cgptr->ewgts = NULL;
    }

    ewgt_sum = 0;
    for (i = cv2v_ptrs[cvtx + 1] - cv2v_ptrs[cvtx]; i; i--) {
      vtx = *sptr++;

      iptr = graph[vtx]->edges;
      if (using_ewgts) {
        fptr = graph[vtx]->ewgts;
      }
      for (j = graph[vtx]->nedges - 1; j; j--) {
        neighbor  = *(++iptr);
        cneighbor = v2cv[neighbor];
        if (cneighbor != cvtx) {
          if (using_ewgts) {
            ewgt = *(++fptr);
          }
          ewgt_sum += ewgt;

          /* Seenflags being used as map from cvtx to index. */
          if (seenflag[cneighbor] == 0) { /* New neighbor. */
            cgptr->edges[nseen] = cneighbor;
            if (COARSEN_EWGTS) {
              cgptr->ewgts[nseen] = ewgt;
            }
            seenflag[cneighbor] = nseen++;
          }
          else { /* Already seen neighbor. */
            if (COARSEN_EWGTS) {
              cgptr->ewgts[seenflag[cneighbor]] += ewgt;
            }
          }
        }
        else if (using_ewgts) {
          ++fptr;
        }
      }
    }

    /* Now clear the seenflag values. */
    iptr = cgptr->edges;
    for (j = nseen - 1; j; j--) {
      seenflag[*(++iptr)] = 0;
    }

    if (COARSEN_EWGTS) {
      cgptr->ewgts[0] = -ewgt_sum;
    }
    /* Increment pointers into edges list. */
    cgptr->nedges = nseen;
    eptr += nseen;
    if (COARSEN_EWGTS) {
      ewptr += nseen;
    }

    cnedges += nseen - 1;
  }

  sfree(seenflag);

  /* Form new vertex weights by adding those from contracted edges. */
  if (COARSEN_VWGTS) {
    gptr = graph;
    for (i = 1; i <= nvtxs; i++) {
      cgraph[v2cv[i]]->vwgt += (*(++gptr))->vwgt;
    }
  }

  /* Reduce arrays to actual sizes */
  cnedges /= 2;
  size  = 2 * cnedges + cnvtxs;
  eptr  = edges;
  edges = srealloc(edges, size * sizeof(int));
  if (eptr != edges) { /* Need to reset pointers in graph. */
    for (i = 1; i <= cnvtxs; i++) {
      cgraph[i]->edges = edges;
      edges += cgraph[i]->nedges;
    }
  }

  if (COARSEN_EWGTS) {
    ewptr    = eweights;
    eweights = srealloc(eweights, size * sizeof(float));
    if (ewptr != eweights) { /* Need to reset pointers in graph. */
      for (i = 1; i <= cnvtxs; i++) {
        cgraph[i]->ewgts = eweights;
        eweights += cgraph[i]->nedges;
      }
    }
  }

  /* If desired, make new vtx coordinates = center-of-mass of their parents. */
  if (coords != NULL && ccoords != NULL && igeom > 0) {
    makeccoords(graph, cnvtxs, cv2v_ptrs, cv2v_vals, igeom, coords, ccoords);
  }

  *pcnedges = cnedges;

  sfree(cv2v_ptrs);
  sfree(cv2v_vals);

  if (DEBUG_COARSEN > 0) {
    printf(" Coarse graph has %d vertices and %d edges\n", cnvtxs, cnedges);
  }

  make_cgraph_time += seconds() - time;
}

static void makecv2v(int  nvtxs,     /* number of vertices in graph */
                     int  cnvtxs,    /* number of vtxs in coarsened graph */
                     int *v2cv,      /* mapping from vtxs to coarsened vtxs */
                     int *cv2v_vals, /* vtxs corresponding to each cvtx */
                     int *cv2v_ptrs  /* indices into cv2c_vals */
)

{
  int sum; /* cumulative offsets into vals array */
  int i;   /* loop counter */

  /* First find number of vtxs associated with each coarse graph vtx. */

  for (i = 1; i <= cnvtxs + 1; i++) {
    cv2v_ptrs[i] = 0;
  }

  for (i = 1; i <= nvtxs; i++) {
    ++cv2v_ptrs[v2cv[i] + 1]; /* +1 offsets and simplifies next loop. */
  }
  cv2v_ptrs[1] = 0;

  /* Now make this a cumulative total to index into cv2v_vals. */
  sum = 0;
  for (i = 2; i <= cnvtxs + 1; i++) {
    cv2v_ptrs[i] += sum;
    sum = cv2v_ptrs[i];
  }

  /* Now go ahead and set the cv2v_vals. */
  for (i = 1; i <= nvtxs; i++) {
    cv2v_vals[cv2v_ptrs[v2cv[i]]] = i;
    ++cv2v_ptrs[v2cv[i]];
  }

  /* Finally, reset the cv2v_ptrs values. */
  for (i = cnvtxs; i; i--) {
    cv2v_ptrs[i] = cv2v_ptrs[i - 1];
  }
  cv2v_ptrs[1] = 0;
}
