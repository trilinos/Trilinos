/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "defs.h"
#include "params.h"
#include "smalloc.h"
#include "structs.h"
#include <stdio.h>

int DEBUG_COVER  = 0;     /* Debug output in min vtx cover routines? */
int FLATTEN      = FALSE; /* Merge indistinguishable vtxs first? */
int COARSE_KLV   = TRUE;  /* Use KLV as multilevel refinement? */
int COARSE_BPM   = FALSE; /* Use vertex cover as ML refinement? */
int VERTEX_COVER = 0;     /* Do vertex cover to refine separator? */

/* Once, or iteratively? */

void coarsen_klv(
    /* Coarsen until nvtxs < vmax, compute and uncoarsen. */
    struct vtx_data **graph,        /* array of vtx data for graph */
    int               nvtxs,        /* number of vertices in graph */
    int               nedges,       /* number of edges in graph */
    int               using_vwgts,  /* are vertices weights being used? */
    int               using_ewgts,  /* are edge weights being used? */
    float *           term_wgts[],  /* weights for terminal propagation */
    int               igeom,        /* dimension for geometric information */
    float **          coords,       /* coordinates for vertices */
    int               vwgt_max,     /* largest vertex weight */
    int *             assignment,   /* processor each vertex gets assigned to */
    double *          goal,         /* desired set sizes */
    int               architecture, /* 0 => hypercube, d => d-dimensional mesh */
    int (*hops)[MAXSETS],           /* cost of edge between sets */
    int     solver_flag,            /* which eigensolver to use */
    int     ndims,                  /* number of eigenvectors to calculate */
    int     nsets,                  /* number of sets being divided into */
    int     vmax,                   /* largest subgraph to stop coarsening */
    int     mediantype,             /* flag for different assignment strategies */
    int     mkconnected,            /* enforce connectivity before eigensolve? */
    double  eigtol,                 /* tolerance in eigen calculation */
    int     nstep,                  /* number of coarsenings between RQI steps */
    int     step,                   /* current step number */
    int **  pbndy_list,             /* pointer to returned boundary list */
    double *weights,                /* weights of vertices in each set */
    int     give_up                 /* has coarsening bogged down? */
)
{
  extern FILE *     Output_File;            /* output file or null */
  extern int        DEBUG_TRACE;            /* trace the execution of the code */
  extern int        DEBUG_COARSEN;          /* debug flag for coarsening */
  extern double     COARSEN_RATIO_MIN;      /* vtx reduction demanded */
  extern int        COARSEN_VWGTS;          /* use vertex weights while coarsening? */
  extern int        COARSEN_EWGTS;          /* use edge weights while coarsening? */
  extern int        LIMIT_KL_EWGTS;         /* limit edges weights in KL? */
  extern int        COARSE_KLV;             /* apply klv as a smoother? */
  extern int        COARSE_BPM;             /* apply bipartite matching/flow as a smoother? */
  extern double     KL_IMBALANCE;           /* fractional imbalance allowed in KL */
  struct vtx_data **cgraph;                 /* array of vtx data for coarsened graph */
  double            new_goal[MAXSETS];      /* new goal if not using vertex weights */
  double *          real_goal;              /* chooses between goal and new_goal */
  double            total_weight;           /* total weight of vertices */
  double            goal_weight;            /* total weight of vertices in goal */
  float *           cterm_wgts[MAXSETS];    /* terminal weights for coarse graph */
  float *           new_term_wgts[MAXSETS]; /* modified for Bui's method */
  float **          real_term_wgts;         /* which of previous two to use */
  float             ewgt_max;               /* largest edge weight in graph */
  float *           twptr      = NULL;      /* loops through term_wgts */
  float *           twptr_save = NULL;      /* copy of twptr */
  float *           ctwptr;                 /* loops through cterm_wgts */
  float **          ccoords;                /* coarse graph coordinates */
  int *             v2cv;                   /* mapping from vtxs to coarse vtxs */
  int *             flag;                   /* scatter array for coarse bndy vtxs */
  int *             bndy_list;              /* list of vertices on boundary */
  int *             cbndy_list;             /* list of vertices of coarse graph on boundary */
  int *             cassignment;            /* set assignments for coarsened vertices */
  int               flattened;              /* was this graph flattened? */
  int               list_length;            /* length of boundary vtx list */
  int               cnvtxs;                 /* number of vertices in coarsened graph */
  int               cnedges;                /* number of edges in coarsened graph */
  int               cvwgt_max;              /* largest vertex weight in coarsened graph */
  int               nextstep;               /* next step in RQI test */
  int               max_dev;                /* largest allowed deviation from balance */
  int               i, j;                   /* loop counters */
  int               find_bndy(), flatten();
  double            find_maxdeg();
  void              makevwsqrt(), make_connected(), print_connected(), eigensolve();
  void              make_unconnected(), assign(), klvspiff(), coarsen1(), free_graph();
  void              compress_ewgts(), restore_ewgts(), count_weights(), bpm_improve();
  void              simple_part();

  if (DEBUG_COARSEN > 0 || DEBUG_TRACE > 0) {
    printf("<Entering coarsen_kl, step=%d, nvtxs=%d, nedges=%d, vmax=%d>\n", step, nvtxs, nedges,
           vmax);
  }

  /* Is problem small enough to solve? */
  if (nvtxs <= vmax || give_up) {
    real_goal = goal;

    simple_part(graph, nvtxs, assignment, nsets, 1, real_goal);
    list_length = find_bndy(graph, nvtxs, assignment, 2, &bndy_list);

    count_weights(graph, nvtxs, assignment, nsets + 1, weights, 1);

    max_dev      = (step == 0) ? vwgt_max : 5 * vwgt_max;
    total_weight = 0;
    for (i = 0; i < nsets; i++) {
      total_weight += real_goal[i];
    }

    if (max_dev > total_weight) {
      max_dev = vwgt_max;
    }
    goal_weight = total_weight * KL_IMBALANCE / nsets;
    if (goal_weight > max_dev) {
      max_dev = goal_weight;
    }

    if (COARSE_KLV) {
      klvspiff(graph, nvtxs, assignment, real_goal, max_dev, &bndy_list, weights);
    }
    if (COARSE_BPM) {
      bpm_improve(graph, assignment, real_goal, max_dev, &bndy_list, weights, using_vwgts);
    }
    *pbndy_list = bndy_list;
    return;
  }

  /* Otherwise I have to coarsen. */
  flattened = FALSE;
  if (coords != NULL) {
    ccoords = smalloc(igeom * sizeof(float *));
  }
  else {
    ccoords = NULL;
  }
  if (FLATTEN && step == 0) {
    flattened = flatten(graph, nvtxs, nedges, &cgraph, &cnvtxs, &cnedges, &v2cv,
                        using_ewgts && COARSEN_EWGTS, igeom, coords, ccoords);
  }
  if (!flattened) {
    coarsen1(graph, nvtxs, nedges, &cgraph, &cnvtxs, &cnedges, &v2cv, igeom, coords, ccoords,
             using_ewgts && COARSEN_EWGTS);
  }

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

  /* If coarsening isn't working very well, give up and partition. */
  give_up = FALSE;
  if (nvtxs * COARSEN_RATIO_MIN < cnvtxs && cnvtxs > vmax && !flattened) {
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

  /* Now recurse on coarse subgraph. */
  if (COARSEN_VWGTS) {
    cvwgt_max = 0;
    for (i = 1; i <= cnvtxs; i++) {
      if (cgraph[i]->vwgt > cvwgt_max) {
        cvwgt_max = cgraph[i]->vwgt;
      }
    }
  }

  else {
    cvwgt_max = 1;
  }

  cassignment = smalloc((cnvtxs + 1) * sizeof(int));
  if (flattened) {
    nextstep = step;
  }
  else {
    nextstep = step + 1;
  }
  coarsen_klv(cgraph, cnvtxs, cnedges, COARSEN_VWGTS, COARSEN_EWGTS, cterm_wgts, igeom, ccoords,
              cvwgt_max, cassignment, goal, architecture, hops, solver_flag, ndims, nsets, vmax,
              mediantype, mkconnected, eigtol, nstep, nextstep, &cbndy_list, weights, give_up);

  /* Interpolate assignment back to fine graph. */
  for (i = 1; i <= nvtxs; i++) {
    assignment[i] = cassignment[v2cv[i]];
  }

  /* Construct boundary list from coarse boundary list. */
  flag = smalloc((cnvtxs + 1) * sizeof(int));
  for (i = 1; i <= cnvtxs; i++) {
    flag[i] = FALSE;
  }
  for (i = 0; cbndy_list[i]; i++) {
    flag[cbndy_list[i]] = TRUE;
  }

  list_length = 0;
  for (i = 1; i <= nvtxs; i++) {
    if (flag[v2cv[i]]) {
      ++list_length;
    }
  }

  bndy_list = smalloc((list_length + 1) * sizeof(int));

  list_length = 0;
  for (i = 1; i <= nvtxs; i++) {
    if (flag[v2cv[i]]) {
      bndy_list[list_length++] = i;
    }
  }
  bndy_list[list_length] = 0;

  sfree(flag);
  sfree(cbndy_list);

  /* Free the space that was allocated. */
  sfree(cassignment);
  if (twptr_save != NULL) {
    sfree(twptr_save);
    twptr_save = NULL;
  }
  free_graph(cgraph);
  sfree(v2cv);

  /* Smooth using KL or BPM every nstep steps. */
  if (!(step % nstep) && !flattened) {
    if (!COARSEN_VWGTS && step != 0) { /* Construct new goal */
      goal_weight = 0;
      for (i = 0; i < nsets; i++) {
        goal_weight += goal[i];
      }
      for (i = 0; i < nsets; i++) {
        new_goal[i] = goal[i] * (nvtxs / goal_weight);
      }
      real_goal = new_goal;
    }
    else {
      real_goal = goal;
    }

    if (LIMIT_KL_EWGTS) {
      find_maxdeg(graph, nvtxs, using_ewgts, &ewgt_max);
      compress_ewgts(graph, nvtxs, nedges, ewgt_max, using_ewgts);
    }

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

    max_dev      = (step == 0) ? vwgt_max : 5 * vwgt_max;
    total_weight = 0;
    for (i = 0; i < nsets; i++) {
      total_weight += real_goal[i];
    }
    if (max_dev > total_weight) {
      max_dev = vwgt_max;
    }
    goal_weight = total_weight * KL_IMBALANCE / nsets;
    if (goal_weight > max_dev) {
      max_dev = goal_weight;
    }

    if (!COARSEN_VWGTS) {
      count_weights(graph, nvtxs, assignment, nsets + 1, weights, (vwgt_max != 1));
    }

    if (COARSE_KLV) {
      klvspiff(graph, nvtxs, assignment, real_goal, max_dev, &bndy_list, weights);
    }
    if (COARSE_BPM) {
      bpm_improve(graph, assignment, real_goal, max_dev, &bndy_list, weights, using_vwgts);
    }

    if (real_term_wgts != term_wgts && new_term_wgts[1] != NULL) {
      sfree(real_term_wgts[1]);
    }

    if (LIMIT_KL_EWGTS) {
      restore_ewgts(graph, nvtxs);
    }
  }

  *pbndy_list = bndy_list;

  if (twptr_save != NULL) {
    sfree(twptr_save);
    twptr_save = NULL;
  }

  /* Free the space that was allocated. */
  if (ccoords != NULL) {
    for (i = 0; i < igeom; i++) {
      sfree(ccoords[i]);
    }
    sfree(ccoords);
  }

  if (DEBUG_COARSEN > 0) {
    printf(" Leaving coarsen_klv, step=%d\n", step);
  }
}

void print_sep_size(int *list, struct vtx_data **graph /* array of vtx data for graph */
)
{
  int sep_size, sep_weight;
  int i;

  sep_size = sep_weight = 0;
  for (i = 0; list[i] != 0; i++) {
    sep_size++;
    sep_weight += graph[list[i]]->vwgt;
  }
  printf(" Sep_size = %d, Sep_weight = %d\n", sep_size, sep_weight);
}
