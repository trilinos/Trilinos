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

#include "defs.h"
#include "params.h"
#include "smalloc.h"
#include "structs.h"
#include <math.h>
#include <stdio.h>

void coarsen(
    /* Coarsen until nvtxs <= vmax, compute and uncoarsen. */
    struct vtx_data **graph,       /* array of vtx data for graph */
    int               nvtxs,       /* number of vertices in graph */
    int               nedges,      /* number of edges in graph */
    int               using_vwgts, /* are vertices weights being used? */
    int               using_ewgts, /* are edge weights being used? */
    float *           term_wgts[], /* terminal weights */
    int               igeom,       /* dimension for geometric information */
    float **          coords,      /* coordinates for vertices */
    double **         yvecs,       /* eigenvectors returned */
    int               ndims,       /* number of eigenvectors to calculate */
    int               solver_flag, /* which eigensolver to use */
    int               vmax,        /* largest subgraph to stop coarsening */
    double            eigtol,      /* tolerance in eigen calculation */
    int               nstep,       /* number of coarsenings between RQI steps */
    int               step,        /* current step number */
    int               give_up      /* has coarsening bogged down? */
)
{
  extern FILE *     Output_File;            /* output file or null */
  extern int        DEBUG_COARSEN;          /* debug flag for coarsening */
  extern int        PERTURB;                /* was matrix perturbed in Lanczos? */
  extern double     COARSEN_RATIO_MIN;      /* min vtx reduction for coarsening */
  extern int        COARSEN_VWGTS;          /* use vertex weights while coarsening? */
  extern int        COARSEN_EWGTS;          /* use edge weights while coarsening? */
  extern double     refine_time;            /* time for RQI/Symmlq iterative refinement */
  struct vtx_data **cgraph;                 /* array of vtx data for coarsened graph */
  struct orthlink * orthlist;               /* list of lower evecs to suppress */
  struct orthlink * newlink;                /* lower evec to suppress */
  double *          cyvecs[MAXDIMS + 1];    /* eigenvectors for subgraph */
  double            evals[MAXDIMS + 1];     /* eigenvalues returned */
  double            goal[MAXSETS];          /* needed for convergence mode = 1 */
  double *          r1, *r2, *work;         /* space needed by symmlq/RQI */
  double *          v, *w, *x, *y;          /* space needed by symmlq/RQI */
  double *          gvec;                   /* rhs vector in extended eigenproblem */
  double            evalest;                /* eigenvalue estimate returned by RQI */
  double            maxdeg;                 /* maximum weighted degree of a vertex */
  float **          ccoords;                /* coordinates for coarsened graph */
  float *           cterm_wgts[MAXSETS];    /* coarse graph terminal weights */
  float *           new_term_wgts[MAXSETS]; /* terminal weights for Bui's method*/
  float **          real_term_wgts;         /* one of the above */
  float *           twptr      = NULL;      /* loops through term_wgts */
  float *           twptr_save = NULL;      /* copy of twptr */
  float *           ctwptr;                 /* loops through cterm_wgts */
  double *          vwsqrt = NULL;          /* square root of vertex weights */
  double            norm, alpha;            /* values used for orthogonalization */
  double            initshift;              /* initial shift for RQI */
  double            total_vwgt;             /* sum of all the vertex weights */
  double            w1, w2;                 /* weights of two sets */
  double            term_tot;               /* sum of all terminal weights */
  int *             space;                  /* room for assignment in Lanczos */
  int *             morespace;              /* room for assignment in Lanczos */
  int *             v2cv;                   /* mapping from vertices to coarse vtxs */
  int               vwgt_max;               /* largest vertex weight */
  int               oldperturb;             /* saves PERTURB value */
  int               cnvtxs;                 /* number of vertices in coarsened graph */
  int               cnedges;                /* number of edges in coarsened graph */
  int               nextstep;               /* next step in RQI test */
  int               nsets;                  /* number of sets being created */
  int               i, j;                   /* loop counters */
  double            time;                   /* time marker */

  double           dot(), ch_normalize(), find_maxdeg(), seconds();
  struct orthlink *makeorthlnk();
  void             makevwsqrt(), eigensolve(), coarsen1(), orthogvec(), rqi_ext();
  void             ch_interpolate(), orthog1(), rqi(), scadd(), free_graph();

  if (DEBUG_COARSEN > 0) {
    printf("<Entering coarsen, step=%d, nvtxs=%d, nedges=%d, vmax=%d>\n", step, nvtxs, nedges,
           vmax);
  }

  nsets = 1 << ndims;

  /* Is problem small enough to solve? */
  if (nvtxs <= vmax || give_up) {
    if (using_vwgts) {
      vwsqrt = smalloc((nvtxs + 1) * sizeof(double));
      makevwsqrt(vwsqrt, graph, nvtxs);
    }
    else {
      vwsqrt = NULL;
    }
    maxdeg = find_maxdeg(graph, nvtxs, using_ewgts, (float *)NULL);

    if (using_vwgts) {
      vwgt_max   = 0;
      total_vwgt = 0;
      for (i = 1; i <= nvtxs; i++) {
        if (graph[i]->vwgt > vwgt_max) {
          vwgt_max = graph[i]->vwgt;
        }
        total_vwgt += graph[i]->vwgt;
      }
    }
    else {
      vwgt_max   = 1;
      total_vwgt = nvtxs;
    }
    for (i = 0; i < nsets; i++) {
      goal[i] = total_vwgt / nsets;
    }

    space = smalloc((nvtxs + 1) * sizeof(int));

    /* If not coarsening ewgts, then need care with term_wgts. */
    if (!using_ewgts && term_wgts[1] != NULL && step != 0) {
      twptr      = smalloc((nvtxs + 1) * (nsets - 1) * sizeof(float));
      twptr_save = twptr;
      for (j = 1; j < nsets; j++) {
        new_term_wgts[j] = twptr;
        twptr += nvtxs + 1;
      }

      for (j = 1; j < nsets; j++) {
        twptr  = term_wgts[j];
        ctwptr = new_term_wgts[j];
        for (i = 1; i <= nvtxs; i++) {
          if (twptr[i] > .5) {
            ctwptr[i] = 1;
          }
          else if (twptr[i] < -.5) {
            ctwptr[i] = -1;
          }
          else {
            ctwptr[i] = 0;
          }
        }
      }
      real_term_wgts = new_term_wgts;
    }
    else {
      real_term_wgts   = term_wgts;
      new_term_wgts[1] = NULL;
    }

    eigensolve(graph, nvtxs, nedges, maxdeg, vwgt_max, vwsqrt, using_vwgts, using_ewgts,
               real_term_wgts, igeom, coords, yvecs, evals, 0, space, goal, solver_flag, FALSE, 0,
               ndims, 3, eigtol);

    if (real_term_wgts != term_wgts && new_term_wgts[1] != NULL) {
      sfree(real_term_wgts[1]);
    }
    sfree(space);
    space = NULL;
    sfree(vwsqrt);
    vwsqrt = NULL;
    sfree(twptr_save);
    twptr_save = NULL;
    return;
  }

  /* Otherwise I have to coarsen. */
  if (coords != NULL) {
    ccoords = smalloc(igeom * sizeof(float *));
  }
  else {
    ccoords = NULL;
  }
  coarsen1(graph, nvtxs, nedges, &cgraph, &cnvtxs, &cnedges, &v2cv, igeom, coords, ccoords,
           using_ewgts);

  /* If coarsening isn't working very well, give up and partition. */
  give_up = FALSE;
  if (nvtxs * COARSEN_RATIO_MIN < cnvtxs && cnvtxs > vmax) {
    printf("WARNING: Coarsening not making enough progress, nvtxs = %d, cnvtxs = %d.\n", nvtxs,
           cnvtxs);
    printf("         Recursive coarsening being stopped prematurely.\n");
    if (Output_File != NULL) {
      fprintf(Output_File,
              "WARNING: Coarsening not making enough progress, nvtxs = %d, cnvtxs = %d.\n", nvtxs,
              cnvtxs);
      fprintf(Output_File, "         Recursive coarsening being stopped prematurely.\n");
    }
    give_up = TRUE;
  }

  /* Create space for subgraph yvecs. */
  for (i = 1; i <= ndims; i++) {
    cyvecs[i] = smalloc((cnvtxs + 1) * sizeof(double));
  }

  /* Make coarse version of terminal weights. */
  if (term_wgts[1] != NULL) {
    twptr      = smalloc((cnvtxs + 1) * (nsets - 1) * sizeof(float));
    twptr_save = twptr;
    for (i = (cnvtxs + 1) * (nsets - 1); i; i--) {
      *twptr++ = 0;
    }
    twptr = twptr_save;
    for (j = 1; j < nsets; j++) {
      cterm_wgts[j] = twptr;
      twptr += cnvtxs + 1;
    }
    for (j = 1; j < nsets; j++) {
      ctwptr = cterm_wgts[j];
      twptr  = term_wgts[j];
      for (i = 1; i < nvtxs; i++) {
        ctwptr[v2cv[i]] += twptr[i];
      }
    }
  }
  else {
    cterm_wgts[1] = NULL;
  }

  /* Now recurse on coarse subgraph. */
  nextstep = step + 1;
  coarsen(cgraph, cnvtxs, cnedges, COARSEN_VWGTS, COARSEN_EWGTS, cterm_wgts, igeom, ccoords, cyvecs,
          ndims, solver_flag, vmax, eigtol, nstep, nextstep, give_up);

  ch_interpolate(yvecs, cyvecs, ndims, graph, nvtxs, v2cv, using_ewgts);

  sfree(twptr_save);
  twptr_save = NULL;
  sfree(v2cv);
  v2cv = NULL;

  /* I need to do Rayleigh Quotient Iteration each nstep stages. */
  time = seconds();
  if (!(step % nstep)) {
    oldperturb = PERTURB;
    PERTURB    = FALSE;
    /* Should I do some orthogonalization here against vwsqrt? */
    if (using_vwgts) {
      vwsqrt = smalloc((nvtxs + 1) * sizeof(double));
      makevwsqrt(vwsqrt, graph, nvtxs);

      for (i = 1; i <= ndims; i++) {
        orthogvec(yvecs[i], 1, nvtxs, vwsqrt);
      }
    }
    else {
      for (i = 1; i <= ndims; i++) {
        orthog1(yvecs[i], 1, nvtxs);
      }
    }

    /* Allocate space that will be needed in RQI. */
    r1   = smalloc(7 * (nvtxs + 1) * sizeof(double));
    r2   = &r1[nvtxs + 1];
    v    = &r1[2 * (nvtxs + 1)];
    w    = &r1[3 * (nvtxs + 1)];
    x    = &r1[4 * (nvtxs + 1)];
    y    = &r1[5 * (nvtxs + 1)];
    work = &r1[6 * (nvtxs + 1)];

    if (using_vwgts) {
      vwgt_max   = 0;
      total_vwgt = 0;
      for (i = 1; i <= nvtxs; i++) {
        if (graph[i]->vwgt > vwgt_max) {
          vwgt_max = graph[i]->vwgt;
        }
        total_vwgt += graph[i]->vwgt;
      }
    }
    else {
      vwgt_max   = 1;
      total_vwgt = nvtxs;
    }
    for (i = 0; i < nsets; i++) {
      goal[i] = total_vwgt / nsets;
    }

    space     = smalloc((nvtxs + 1) * sizeof(int));
    morespace = smalloc((nvtxs) * sizeof(int));

    initshift = 0;
    orthlist  = NULL;
    for (i = 1; i < ndims; i++) {
      ch_normalize(yvecs[i], 1, nvtxs);
      rqi(graph, yvecs, i, nvtxs, r1, r2, v, w, x, y, work, eigtol, initshift, &evalest, vwsqrt,
          orthlist, 0, nsets, space, morespace, 3, goal, vwgt_max, ndims);

      /* Now orthogonalize higher yvecs against this one. */
      norm = dot(yvecs[i], 1, nvtxs, yvecs[i]);
      for (j = i + 1; j <= ndims; j++) {
        alpha = -dot(yvecs[j], 1, nvtxs, yvecs[i]) / norm;
        scadd(yvecs[j], 1, nvtxs, alpha, yvecs[i]);
      }

      /* Now prepare for next pass through loop. */
      initshift     = evalest;
      newlink       = makeorthlnk();
      newlink->vec  = yvecs[i];
      newlink->pntr = orthlist;
      orthlist      = newlink;
    }
    ch_normalize(yvecs[ndims], 1, nvtxs);

    if (term_wgts[1] != NULL && ndims == 1) {
      /* Solve extended eigen problem */

      /* If not coarsening ewgts, then need care with term_wgts. */
      if (!using_ewgts && term_wgts[1] != NULL && step != 0) {
        twptr      = smalloc((nvtxs + 1) * (nsets - 1) * sizeof(float));
        twptr_save = twptr;
        for (j = 1; j < nsets; j++) {
          new_term_wgts[j] = twptr;
          twptr += nvtxs + 1;
        }

        for (j = 1; j < nsets; j++) {
          twptr  = term_wgts[j];
          ctwptr = new_term_wgts[j];
          for (i = 1; i <= nvtxs; i++) {
            if (twptr[i] > .5) {
              ctwptr[i] = 1;
            }
            else if (twptr[i] < -.5) {
              ctwptr[i] = -1;
            }
            else {
              ctwptr[i] = 0;
            }
          }
        }
        real_term_wgts = new_term_wgts;
      }
      else {
        real_term_wgts   = term_wgts;
        new_term_wgts[1] = NULL;
      }

      /* Following only works for bisection. */
      w1       = goal[0];
      w2       = goal[1];
      gvec     = smalloc((nvtxs + 1) * sizeof(double));
      term_tot = 0;
      for (j = 1; j <= nvtxs; j++) {
        term_tot += (real_term_wgts[1])[j];
      }
      term_tot /= (w1 + w2);
      if (using_vwgts) {
        for (j = 1; j <= nvtxs; j++) {
          gvec[j] = (real_term_wgts[1])[j] / graph[j]->vwgt - term_tot;
        }
      }
      else {
        for (j = 1; j <= nvtxs; j++) {
          gvec[j] = (real_term_wgts[1])[j] - term_tot;
        }
      }

      rqi_ext();

      sfree(gvec);
      gvec = NULL;
      if (real_term_wgts != term_wgts && new_term_wgts[1] != NULL) {
        sfree(new_term_wgts[1]);
        new_term_wgts[1] = NULL;
      }
    }
    else {
      rqi(graph, yvecs, ndims, nvtxs, r1, r2, v, w, x, y, work, eigtol, initshift, &evalest, vwsqrt,
          orthlist, 0, nsets, space, morespace, 3, goal, vwgt_max, ndims);
    }
    refine_time += seconds() - time;

    /* Free the space allocated for RQI. */
    sfree(morespace);
    sfree(space);
    while (orthlist != NULL) {
      newlink = orthlist->pntr;
      sfree(orthlist);
      orthlist = newlink;
    }
    sfree(r1);
    sfree(vwsqrt);
    vwsqrt  = NULL;
    PERTURB = oldperturb;
  }
  if (DEBUG_COARSEN > 0) {
    printf(" Leaving coarsen, step=%d\n", step);
  }

  sfree(twptr_save);
  twptr_save = NULL;

  /* Free the space that was allocated. */
  if (ccoords != NULL) {
    for (i = 0; i < igeom; i++) {
      sfree(ccoords[i]);
    }
    sfree(ccoords);
  }
  for (i = ndims; i > 0; i--) {
    sfree(cyvecs[i]);
  }
  free_graph(cgraph);
}
