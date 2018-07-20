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

#include "params.h"  // for MAXDIMS
#include "smalloc.h" // for sfree, smalloc
#include "structs.h" // for vtx_data
#include <stdio.h>   // for NULL, fprintf, printf, etc

void sequence(struct vtx_data **graph,       /* graph data structure */
              int               nvtxs,       /* number of vertices in graph */
              int               nedges,      /* number of edges in graph */
              int               using_ewgts, /* are edge weights being used? */
              double *          vwsqrt,      /* sqrt of vertex weights (length nvtxs+1) */
              int               solver_flag, /* which eigensolver should I use? */
              int               rqi_flag,    /* use multilevel eigensolver? */
              int               vmax,        /* if so, how many vtxs to coarsen down to? */
              double            eigtol       /* tolerance on eigenvectors */
              )
{
  extern char       SEQ_FILENAME[];           /* name of sequence file */
  extern int        RQI_CONVERGENCE_MODE;     /* residual or partition mode? */
  extern int        LANCZOS_CONVERGENCE_MODE; /* residual or partition mode? */
  extern double     sequence_time;            /* time spent sequencing */
  struct vtx_data **subgraph  = NULL;         /* subgraph data structure */
  struct edgeslist *edge_list = NULL;         /* edges added for connectivity */
  double *          yvecs[MAXDIMS + 1];       /* space for pointing to eigenvectors */
  double            evals[MAXDIMS + 1];       /* corresponding eigenvalues */
  double *          subvwsqrt = NULL;         /* vwsqrt vector for subgraphs */
  double *          values    = NULL;         /* sorted Fiedler vector values */
  double *          subvals   = NULL;         /* values for one connected component */
  double            goal[2];                  /* needed for eigen convergence mode = 1 */
  double            total_vwgt;               /* sum of all vertex weights */
  double            time;                     /* time spent sequencing */
  float *           term_wgts[2];             /* dummy vector for terminal weights */
  int *             setsize     = NULL;       /* size of each connected component */
  int *             glob2loc    = NULL;       /* maps graph vtxs to subgraph vtxs */
  int *             loc2glob    = NULL;       /* maps subgraph vtxs to graph vtxs */
  int *             permutation = NULL;       /* computed permutation */
  int *             subperm     = NULL;       /* partial permutation */
  int *             compnum     = NULL;       /* component number for each vertex */
  int *             degree      = NULL;       /* degrees of vertices in subgraph */
  double            maxdeg;                   /* maximum weighted degree of a vertex */
  int               using_vwgts;              /* are vertex weights being used? */
  int               subnvtxs;                 /* number of vertices in subgraph */
  int               subnedges;                /* number of edges in subgraph */
  int               maxsize;                  /* size of largest connected component */
  int               subvwgt_max;              /* largest vertex weight in component */
  int               ncomps;                   /* number of connected components */
  int               comp;                     /* loops over connected components */
  int               nused;                    /* number of vertices already handled */
  int               old_rqi_conv_mode;        /* value of RQI_CONVERGENCE_MODE */
  int               old_lan_conv_mode;        /* value of LANCZOS_CONVERGENCE_MODE */
  int               i;                        /* loop counters */
  int *             space     = NULL;
  FILE *            orderfile = NULL;
  void              mergesort(), eigensolve(), free_edgeslist(), y2x();
  void              make_subvector(), make_subgraph(), remake_graph();
  void              make_maps2();
  double            find_maxdeg(), seconds();
  int               find_edges();

  time        = seconds();
  using_vwgts = (vwsqrt != NULL);

  /* Sort each connected component seperately. */
  compnum     = smalloc((nvtxs + 1) * sizeof(int));
  permutation = smalloc(nvtxs * sizeof(int));
  values      = smalloc((nvtxs + 1) * sizeof(double));

  space  = smalloc(nvtxs * sizeof(int));
  ncomps = find_edges(graph, nvtxs, compnum, space, &edge_list);
  ++ncomps;

  free_edgeslist(edge_list);
  yvecs[1] = smalloc((nvtxs + 1) * sizeof(double));

  term_wgts[1] = NULL;
  subvwsqrt    = NULL;

  old_rqi_conv_mode = RQI_CONVERGENCE_MODE;
  old_lan_conv_mode = LANCZOS_CONVERGENCE_MODE;

  RQI_CONVERGENCE_MODE     = 0;
  LANCZOS_CONVERGENCE_MODE = 0;

  maxsize = nvtxs;
  if (ncomps > 1) {

    /* Find size of largest set. */
    setsize = smalloc(ncomps * sizeof(int));
    for (comp = 0; comp < ncomps; comp++) {
      setsize[comp] = 0;
    }
    for (i = 1; i <= nvtxs; i++) {
      ++setsize[compnum[i]];
    }
    maxsize = 0;
    for (comp = 0; comp < ncomps; comp++) {
      if (setsize[comp] > maxsize) {
        maxsize = setsize[comp];
      }
    }

    glob2loc = smalloc((nvtxs + 1) * sizeof(int));
    loc2glob = smalloc((maxsize + 1) * sizeof(int));
    subgraph = smalloc((maxsize + 1) * sizeof(struct vtx_data *));
    degree   = smalloc((maxsize + 1) * sizeof(int));
    if (using_vwgts) {
      subvwsqrt = smalloc((maxsize + 1) * sizeof(double));
    }
  }

  nused = 0;

  for (comp = 0; comp < ncomps; comp++) {
    subperm = &(permutation[nused]);
    subvals = &(values[nused]);

    if (ncomps > 1) {
      make_maps2(compnum, nvtxs, comp, glob2loc, loc2glob);
      subnvtxs = setsize[comp];
      make_subgraph(graph, subgraph, subnvtxs, &subnedges, compnum, comp, glob2loc, loc2glob,
                    degree, using_ewgts);

      if (using_vwgts) {
        make_subvector(vwsqrt, subvwsqrt, subnvtxs, loc2glob);
      }
    }
    else {
      subgraph  = graph;
      subnvtxs  = nvtxs;
      subnedges = nedges;
      subvwsqrt = vwsqrt;
    }

    if (using_vwgts) {
      subvwgt_max = 0;
      total_vwgt  = 0;
      for (i = 1; i <= subnvtxs; i++) {
        if (subgraph[i]->vwgt > subvwgt_max) {
          subvwgt_max = subgraph[i]->vwgt;
        }
        total_vwgt += subgraph[i]->vwgt;
      }
    }
    else {
      subvwgt_max = 1;
      total_vwgt  = subnvtxs;
    }

    goal[0] = goal[1] = total_vwgt / 2;

    if (subnvtxs == 1) {
      yvecs[1][1] = 0;
    }
    else {
      maxdeg = find_maxdeg(subgraph, subnvtxs, using_ewgts, (float *)NULL);

      eigensolve(subgraph, subnvtxs, subnedges, maxdeg, subvwgt_max, subvwsqrt, using_vwgts,
                 using_ewgts, term_wgts, 0, (float **)NULL, yvecs, evals, 0, space, goal,
                 solver_flag, rqi_flag, vmax, 1, 3, eigtol);

      /* Sort values in eigenvector */
      if (using_vwgts) {
        y2x(yvecs, 1, subnvtxs, subvwsqrt);
      }
      mergesort(&(yvecs[1][1]), subnvtxs, subperm, space);
    }

    if (ncomps > 1) {
      remake_graph(subgraph, subnvtxs, loc2glob, degree, using_ewgts);
      for (i = 0; i < subnvtxs; i++) {
        subvals[i] = yvecs[1][subperm[i] + 1];
        subperm[i] = loc2glob[subperm[i] + 1];
      }
    }
    else {
      for (i = 0; i < nvtxs; i++) {
        ++subperm[i];
      }
    }

    nused += subnvtxs;
  }

  /* Turn off timer. */
  sequence_time += seconds() - time;

  /* Print the ordering to a file. */
  /* Note that permutation is 1<->n. */
  orderfile = fopen(SEQ_FILENAME, "w");
  if (ncomps == 1) {
    for (i = 0; i < nvtxs; i++) {
      fprintf(orderfile, "%-7d   %9.6f\n", permutation[i], yvecs[1][permutation[i]]);
    }
  }
  else {
    for (i = 0; i < nvtxs; i++) {
      fprintf(orderfile, "%-7d   %9.6f\n", permutation[i], values[i]);
    }
  }
  fclose(orderfile);

  RQI_CONVERGENCE_MODE     = old_rqi_conv_mode;
  LANCZOS_CONVERGENCE_MODE = old_lan_conv_mode;

  if (ncomps > 1) {
    sfree(degree);
    sfree(subgraph);
    sfree(loc2glob);
    sfree(glob2loc);
    sfree(setsize);
    if (using_vwgts) {
      sfree(subvwsqrt);
    }
  }

  sfree(yvecs[1]);
  sfree(space);
  sfree(values);
  sfree(permutation);
  sfree(compnum);

  printf("%d connected components found\n", ncomps);
  printf("Sequence printed to file `%s'\n", SEQ_FILENAME);
}
