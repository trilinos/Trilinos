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

#include "smalloc.h" // for smalloc, sfree
#include "structs.h" // for vtx_data
#include <stdio.h>   // for NULL, printf

extern void makeccoords();

void makecgraph(struct vtx_data ** graph,       /* array of vtx data for graph */
                int                nvtxs,       /* number of vertices in graph */
                struct vtx_data ***pcgraph,     /* coarsened version of graph */
                int *              pcnvtxs,     /* number of vtxs in coarsened graph */
                int *              pcnedges,    /* number of edges in coarsened graph */
                int *              mflag,       /* flag indicating vtx matched or not */
                int *              v2cv,        /* mapping from vtxs to coarsened vtxs */
                int                nmerged,     /* number of merged vertices */
                int                using_ewgts, /* are edge weights being used? */
                int                igeom,       /* dimensions of geometric data */
                float **           coords,      /* coordinates for vertices */
                float **           ccoords      /* coordinates for coarsened vertices */
)
{
  extern double     make_cgraph_time;
  extern int        DEBUG_COARSEN;   /* debug flag for coarsening output */
  extern int        COARSEN_VWGTS;   /* turn off vertex weights in coarse graph? */
  extern int        COARSEN_EWGTS;   /* turn off edge weights in coarse graph? */
  struct vtx_data **cgraph   = NULL; /* coarsened version of graph */
  struct vtx_data * links    = NULL; /* space for all the vertex data */
  struct vtx_data **gptr     = NULL; /* loops through cgraph */
  struct vtx_data * cgptr    = NULL; /* loops through cgraph */
  int *             start    = NULL; /* start of edgevals list for each vertex */
  int *             iptr     = NULL; /* loops through integer arrays */
  int *             seenflag = NULL; /* flags for vtxs already put in edge list */
  int *             sptr     = NULL; /* loops through seenflags */
  float *           eweights = NULL; /* space for edge weights in coarsened graph */
  float *           fptr     = NULL; /* loops through eweights */
  float             ewgt     = 0.0;  /* edge weight */
  double            ewgt_sum;        /* sum of edge weights */
  double            time;            /* timing parameters */
  int               nseen;           /* number of edges of coarse graph seen so far */
  int               cnvtxs;          /* number of vtxs in coarsened graph */
  int               cnedges;         /* twice number of edges in coarsened graph */
  int               neighbor;        /* neighboring vertex */
  int               size;            /* space needed for coarsened graph */
  int *             edges = NULL;    /* space for edges in coarsened graph */
  int               cvtx;            /* vertex number in coarsened graph */
  int               cneighbor;       /* neighboring vertex number in coarsened graph */
  int               i, j;            /* loop counters */
  double            seconds();

  void makev2cv(), countcedges();

  /* Compute the number of vertices and edges in the coarsened graph, */
  /* and construct start pointers into coarsened edge array. */
  time = seconds();

  *pcnvtxs = cnvtxs = nvtxs - nmerged;

  /* Construct mapping from original graph vtxs to coarsened graph vtxs. */
  makev2cv(mflag, nvtxs, v2cv);

  start = smalloc((cnvtxs + 2) * sizeof(int));

  seenflag = smalloc((cnvtxs + 1) * sizeof(int));
  sptr     = seenflag;
  for (i = cnvtxs; i; i--) {
    *(++sptr) = 0;
  }

  countcedges(graph, nvtxs, start, seenflag, mflag, v2cv, pcnedges);
  cnedges = *pcnedges;

  if (DEBUG_COARSEN > 0) {
    printf(" Coarse graph has %d vertices and %d edges\n", cnvtxs, cnedges);
  }

  /* Now allocate space for the new graph. */
  *pcgraph = cgraph = smalloc((cnvtxs + 1) * sizeof(struct vtx_data *));
  links             = smalloc(cnvtxs * sizeof(struct vtx_data));

  size  = 2 * cnedges + cnvtxs;
  edges = smalloc(size * sizeof(int));
  if (COARSEN_EWGTS) {
    eweights = smalloc(size * sizeof(float));
  }

  /* Fill in simple data fields for coarsened graph. */
  /* Edges and weights are put in later. */
  gptr = cgraph;
  sptr = seenflag;
  for (i = 1; i <= cnvtxs; i++) {
    size            = start[i + 1] - start[i] + 1;
    links->nedges   = size;
    links->edges    = edges;
    links->edges[0] = i;
    edges += size;
    if (COARSEN_VWGTS) {
      links->vwgt = 0;
    }
    else {
      links->vwgt = 1;
    }
    if (COARSEN_EWGTS) {
      links->ewgts = eweights;
      eweights += size;
    }
    else {
      links->ewgts = NULL;
    }
    *(++gptr) = links++;
    *(++sptr) = 0;
  }
  sfree(start);

  /* Now form new vertex weights by adding those from contracted edges. */
  if (COARSEN_VWGTS) {
    gptr = graph;
    for (i = 1; i <= nvtxs; i++) {
      cgraph[v2cv[i]]->vwgt += (*(++gptr))->vwgt;
    }
  }

  /* Use the renumbering to fill in the edge lists for the new graph. */
  ewgt = 1;
  for (i = 1; i <= nvtxs; i++) {
    if (mflag[i] > i || mflag[i] == 0) {
      /* Unmatched edge, or first appearance of matched edge. */
      nseen    = 1;
      cvtx     = v2cv[i];
      cgptr    = cgraph[cvtx];
      ewgt_sum = 0;

      iptr = graph[i]->edges;
      if (using_ewgts) {
        fptr = graph[i]->ewgts;
      }
      for (j = graph[i]->nedges - 1; j; j--) {
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

      if (mflag[i] > i) { /* Now handle the matched vertex. */
        iptr = graph[mflag[i]]->edges;
        if (using_ewgts) {
          fptr = graph[mflag[i]]->ewgts;
        }
        for (j = graph[mflag[i]]->nedges - 1; j; j--) {
          neighbor  = *(++iptr);
          cneighbor = v2cv[neighbor];
          if (cneighbor != cvtx) {
            if (using_ewgts) {
              ewgt = *(++fptr);
            }
            ewgt_sum += ewgt;

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
      if (COARSEN_EWGTS) {
        cgptr->ewgts[0] = -ewgt_sum;
      }
      /* Now clear the seenflag values. */
      iptr = cgraph[cvtx]->edges;
      for (j = cgraph[cvtx]->nedges - 1; j; j--) {
        seenflag[*(++iptr)] = 0;
      }
    }
  }

  sfree(seenflag);

  /* If desired, make new vtx coordinates = center-of-mass of their parents. */
  if (coords != NULL && ccoords != NULL && igeom > 0) {
    makeccoords(graph, nvtxs, cnvtxs, mflag, v2cv, igeom, coords, ccoords);
  }

  make_cgraph_time += seconds() - time;
}
