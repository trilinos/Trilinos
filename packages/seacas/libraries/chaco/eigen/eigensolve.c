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

/* Invoke the eigenvector calculation */
void eigensolve(struct vtx_data **graph,        /* graph data structure */
                int               nvtxs,        /* number of vertices in graph */
                int               nedges,       /* number of edges in graph */
                double            maxdeg,       /* largest (weighted) degree of a vertex */
                int               vwgt_max,     /* largest vertex weight */
                double *          vwsqrt,       /* sqrt of vertex weights (length nvtxs+1) */
                int               using_vwgts,  /* are vertex weights being used? */
                int               using_ewgts,  /* are edge weights being used? */
                float *           term_wgts[],  /* terminal propagation weight vector */
                int               igeom,        /* geometric dimensionality if given coords */
                float **          coords,       /* coordinates of vertices */
                double **         yvecs,        /* space for pointing to eigenvectors */
                double *          evals,        /* eigenvalues associated with eigenvectors */
                int               architecture, /* 0 => hypercube, d => d-dimensional mesh */
                int *             assignment,   /* set number of each vtx (length n+1) */
                double *          goal,         /* desired set sizes */
                int               solver_flag,  /* flag indicating which solver to use */
                int               rqi_flag,     /* use multi-level techniques? */
                int               vmax,         /* if so, how many vtxs to coarsen down to? */
                int               ndims,        /* number of eigenvectors (2^d sets) */
                int               mediantype,   /* which partitioning strategy to use */
                double            eigtol        /* tolerance on eigenvectors */
)
{
  extern int    DEBUG_TRACE;              /* trace the execution of the code */
  extern int    DEBUG_EVECS;              /* debug flag for eigenvector generation */
  extern int    DEBUG_PERTURB;            /* debug flag for matrix perturbation */
  extern int    PERTURB;                  /* randomly perturb to break symmetry? */
  extern int    NPERTURB;                 /* number of edges to perturb */
  extern double PERTURB_MAX;              /* maximum size of perturbation */
  extern int    COARSE_NLEVEL_RQI;        /* do RQI this often in uncoarsening */
  extern int    LANCZOS_CONVERGENCE_MODE; /* how to stop Lanczos? */
  extern double lanczos_time;             /* time spent in Lanczos algorithm */
  extern double rqi_symmlq_time;          /* time spent in RQI/Symmlq method */
  extern int    WARNING_EVECS;            /* print warning messages? */
  extern int    LANCZOS_SO_PRECISION;     /* controls precision in eigen calc. */
  extern double SRESTOL;                  /* resid tol for T evec computation */
  extern int    LANCZOS_MAXITNS;          /* maximum Lanczos iterations allowed */
  extern int    EXPERT;                   /* user type */
  double        bound[MAXDIMS + 1];       /* ritz approx bounds to eigenpairs */
  double        time;                     /* time marker */
  float *       dummy_twgt[2];            /* turns off terminal propagation */
  float *       twptr;                    /* terminal propagation weight vector */
  int *         active;                   /* space for nvtxs values */
  int           step;                     /* current step in RQI counting */
  int           nstep;                    /* number of uncoarsening levels between RQIs */
  int           version;                  /* which version of sel. orth. to use */
  int           nsets = 0;                /* number of sets to divide into */
  double *      g;                        /* rhs n-vector in the extended eigenproblem */
  double *      ptr;                      /* loops through yvec */
  double        w1, w2;                   /* desired weights of two sets */
  double        term_tot;                 /* sum of terminal weights */
  double        sigma;                    /* norm constraint on extended eigenvector */
  int           i, j;                     /* loop counter */
  int           normal;                   /* use normal or extended eigensolver? */
  int           autoset_maxitns;          /* set LANCZOS_MAXITNS automatically? */
  int           prev_maxitns = 0;         /* LANCZOS_MAXITNS value above this routine */
  int           autoset_srestol;          /* set SRESTOL automatically? */
  double        prev_srestol = 0;         /* SRESTOL value above this routine */

  double seconds();
  void   coarsen(), lanczos_FO(), lanczos_SO(), vecout(), vecnorm();
  void   lanczos_SO_float(), strout();
  void   perturb_init(), perturb_clear(), x2y(), y2x();
  int    lanczos_ext(), lanczos_ext_float();

  if (DEBUG_TRACE > 0) {
    printf("<Entering eigensolve, nvtxs = %d, nedges = %d>\n", nvtxs, nedges);
  }

  if (nvtxs <= ndims) { /* Pathological special case. */
    for (i = 1; i <= ndims; i++) {
      for (j = 1; j <= nvtxs; j++) {
        yvecs[i][j] = 0;
      }
    }
    return;
  }

  active = NULL;

  /* Autoset (if necessary) some parameters for the eigen calculation */
  autoset_maxitns = FALSE;
  autoset_srestol = FALSE;
  if (LANCZOS_MAXITNS < 0) {
    autoset_maxitns = TRUE;
    prev_maxitns    = LANCZOS_MAXITNS;
    LANCZOS_MAXITNS = 2 * nvtxs;
  }
  if (SRESTOL < 0) {
    autoset_srestol = TRUE;
    prev_srestol    = SRESTOL;
    SRESTOL         = eigtol * eigtol;
  }

  /* Note: When (if ever) rqi_ext is done, change precedence of eigensolvers. */

  if (term_wgts[1] != NULL && ndims == 1) { /* then use lanczos_ext */
    if (PERTURB) {
      if (NPERTURB > 0 && PERTURB_MAX > 0.0) {
        perturb_init(nvtxs);
        if (DEBUG_PERTURB > 0) {
          printf("Matrix being perturbed with scale %e\n", PERTURB_MAX);
        }
      }
      else if (DEBUG_PERTURB > 0) {
        printf("Matrix not being perturbed\n");
      }
    }

    version = 2;
    if (LANCZOS_CONVERGENCE_MODE == 1) {
      active = smalloc(nvtxs * sizeof(int));
    }

    w1    = goal[0];
    w2    = goal[1];
    sigma = sqrt(4 * w1 * w2 / (w1 + w2));
    g     = smalloc((nvtxs + 1) * sizeof(double));

    twptr    = term_wgts[1];
    term_tot = 0;
    for (i = 1; i <= nvtxs; i++) {
      term_tot += twptr[i];
    }
    term_tot /= (w1 + w2);
    if (using_vwgts) {
      for (i = 1; i <= nvtxs; i++) {
        g[i] = twptr[i] / graph[i]->vwgt - term_tot;
      }
    }
    else {
      for (i = 1; i <= nvtxs; i++) {
        g[i] = twptr[i] - term_tot;
      }
    }

    time = seconds();

    if (LANCZOS_SO_PRECISION == 2) { /* double precision */
      normal = lanczos_ext(graph, nvtxs, ndims, yvecs, eigtol, vwsqrt, maxdeg, version, g, sigma);
    }
    else { /* single precision */
      normal =
          lanczos_ext_float(graph, nvtxs, ndims, yvecs, eigtol, vwsqrt, maxdeg, version, g, sigma);
    }

    sfree(g);
    if (active != NULL) {
      sfree(active);
    }
    active = NULL;

    if (normal) {
      if (WARNING_EVECS > 2) {
        strout("WARNING: Not an extended eigenproblem; switching to standard eigensolver.\n");
      }
    }
    else {
      if (w2 != w1) {
        if (using_vwgts) {
          y2x(yvecs, ndims, nvtxs, vwsqrt);
        }
        sigma = (w2 - w1) / (w2 + w1);
        ptr   = yvecs[1];
        for (i = nvtxs; i; i--) {
          *(++ptr) += sigma;
        }
        /* Note: if assign() could skip scaling, next call unnecessary. */
        if (using_vwgts) {
          x2y(yvecs, ndims, nvtxs, vwsqrt);
        }
      }
    }

    if (PERTURB && NPERTURB > 0 && PERTURB_MAX > 0.0) {
      perturb_clear();
    }

    lanczos_time += seconds() - time;
  }
  else {
    normal = TRUE;
  }

  if (normal) {
    if (rqi_flag) {
      /* Solve using multi-level scheme RQI/Symmlq. */
      time          = seconds();
      nstep         = COARSE_NLEVEL_RQI;
      step          = 0;
      dummy_twgt[1] = NULL;
      coarsen(graph, nvtxs, nedges, using_vwgts, using_ewgts, dummy_twgt, igeom, coords, yvecs,
              ndims, solver_flag, vmax, eigtol, nstep, step, FALSE);

      rqi_symmlq_time += seconds() - time;
    }

    else { /* Use standard Lanczos. */
      if (PERTURB) {
        if (NPERTURB > 0 && PERTURB_MAX > 0.0) {
          perturb_init(nvtxs);
          if (DEBUG_PERTURB > 0) {
            printf("Matrix being perturbed with scale %e\n", PERTURB_MAX);
          }
        }
        else if (DEBUG_PERTURB > 0) {
          printf("Matrix not being perturbed\n");
        }
      }

      if (solver_flag == 1) {
        time    = seconds();
        version = 1;
        lanczos_FO(graph, nvtxs, ndims, yvecs, evals, bound, eigtol, vwsqrt, maxdeg, version);
        lanczos_time += seconds() - time;
      }
      if (solver_flag == 2) {
        time    = seconds();
        version = 2;
        lanczos_FO(graph, nvtxs, ndims, yvecs, evals, bound, eigtol, vwsqrt, maxdeg, version);
        lanczos_time += seconds() - time;
      }
      else if (solver_flag == 3) {
        version = 2; /* orthog. against left end only */
        if (LANCZOS_CONVERGENCE_MODE == 1) {
          active = smalloc(nvtxs * sizeof(int));
        }
        nsets = 1 << ndims;
        time  = seconds();
        if (LANCZOS_SO_PRECISION == 2) { /* double precision */
          lanczos_SO(graph, nvtxs, ndims, yvecs, evals, bound, eigtol, vwsqrt, maxdeg, version,
                     architecture, nsets, assignment, active, mediantype, goal, vwgt_max);
        }
        else { /* single precision */
          lanczos_SO_float(graph, nvtxs, ndims, yvecs, evals, bound, eigtol, vwsqrt, maxdeg,
                           version, architecture, nsets, assignment, active, mediantype, goal,
                           vwgt_max);
        }
        lanczos_time += seconds() - time;
      }
      else if (solver_flag == 4) {
        if (EXPERT) {
          version = 1; /* orthog. against both ends */
        }
        else {
          /* this should have been caught earlier ... */
          version = 2;
        }
        if (LANCZOS_CONVERGENCE_MODE == 1) {
          active = smalloc(nvtxs * sizeof(int));
        }
        nsets = 1 << ndims;
        time  = seconds();
        if (LANCZOS_SO_PRECISION == 1) { /* Single precision */
          lanczos_SO_float(graph, nvtxs, ndims, yvecs, evals, bound, eigtol, vwsqrt, maxdeg,
                           version, architecture, nsets, assignment, active, mediantype, goal,
                           vwgt_max);
        }
        else { /* Double precision */
          lanczos_SO(graph, nvtxs, ndims, yvecs, evals, bound, eigtol, vwsqrt, maxdeg, version,
                     architecture, nsets, assignment, active, mediantype, goal, vwgt_max);
        }
        lanczos_time += seconds() - time;
      }
    }

    /*
    file = fopen("CHACO.EVECS", "w");
    for (i = 1; i <= nvtxs; i++) {
     for (j = 1; j <= ndims; j++) {
      fprintf(file, "%g ", (yvecs[j])[i]);
     }
      fprintf(file, "\n");
    }
    fclose(file);
    */

    if (PERTURB && NPERTURB > 0 && PERTURB_MAX > 0.0) {
      perturb_clear();
    }
  }

  /* This is an attempt to reduce some machine-to-machine
   * variance. If the first value in the eigenvector is negative,
   * reflect the eigenvector...  This may not be needed following
   * the addition of the standard random number generator in util/random.c
   */
  for (nstep = 1; nstep <= ndims; nstep++) {
    vecnorm(yvecs[nstep], 1, nvtxs);
  }

  if (DEBUG_EVECS > 4) {
    for (nstep = 1; nstep <= ndims; nstep++) {
      vecout(yvecs[nstep], 1, nvtxs, "Eigenvector", NULL);
    }
  }

  /* Auto-reset (if necessary) some parameters for the eigen calculation */
  if (autoset_maxitns) {
    LANCZOS_MAXITNS = prev_maxitns;
  }
  if (autoset_srestol) {
    SRESTOL = prev_srestol;
  }

  if (active != NULL) {
    sfree(active);
  }

  if (DEBUG_TRACE > 1) {
    printf("<Leaving eigensolve>\n");
  }
}
