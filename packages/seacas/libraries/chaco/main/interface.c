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

#include "defs.h"    // for FALSE, TRUE
#include "smalloc.h" // for sfree, smalloc_ret
#include "structs.h"
#include <stdio.h>  // for NULL, printf, fclose, fopen, etc
#include <stdlib.h> // for free

extern int submain(struct vtx_data **graph, int nvtxs, int nedges, int using_vwgts, int using_ewgts,
                   int igeom, float **coords, char *outassignname, char *outfilename,
                   int *assignment, double *goal, int architecture, int ndims_tot, int mesh_dims[3],
                   int global_method, int local_method, int rqi_flag, int vmax, int ndims,
                   double eigtol, long seed);

int Using_Main = FALSE; /* Is main routine being called? */

int interface(int    nvtxs,                 /* number of vertices in full graph */
              int *  start,                 /* start of edge list for each vertex */
              int *  adjacency,             /* edge list data */
              int *  vwgts,                 /* weights for all vertices */
              float *ewgts,                 /* weights for all edges */
              float *x, float *y, float *z, /* coordinates for inertial method */
              char *  outassignname,        /* name of assignment output file */
              char *  outfilename,          /* output file name */
              int *   assignment,           /* set number of each vtx (length n) */
              int     architecture,         /* 0 => hypercube, d => d-dimensional mesh */
              int     ndims_tot,            /* total number of cube dimensions to divide */
              int     mesh_dims[3],         /* dimensions of mesh of processors */
              double *goal,                 /* desired set sizes for each set */
              int     global_method,        /* global partitioning algorithm */
              int     local_method,         /* local partitioning algorithm */
              int     rqi_flag,             /* should I use RQI/Symmlq eigensolver? */
              int     vmax,                 /* how many vertices to coarsen down to? */
              int     ndims,                /* number of eigenvectors (2^d sets) */
              double  eigtol,               /* tolerance on eigenvectors */
              long    seed                  /* for random graph mutations */
)
{
  extern char *     PARAMS_FILENAME;    /* name of file with parameter updates */
  extern int        MAKE_VWGTS;         /* make vertex weights equal to degrees? */
  extern int        MATCH_TYPE;         /* matching routine to use */
  extern int        FREE_GRAPH;         /* free graph data structure after reformat? */
  extern int        DEBUG_PARAMS;       /* debug flag for reading parameters */
  extern int        DEBUG_TRACE;        /* trace main execution path */
  extern double     start_time;         /* time routine is entered */
  extern double     reformat_time;      /* time spent reformatting graph */
  FILE *            params_file = NULL; /* file for reading new parameters */
  struct vtx_data **graph;              /* graph data structure */
  double            vwgt_sum;           /* sum of vertex weights */
  double            time;               /* timing variable */
  float **          coords;             /* coordinates for vertices if used */
  int *             vptr;               /* loops through vertex weights */
  int               flag;               /* return code from balance */
  int               nedges;             /* number of edges in graph */
  int               using_vwgts;        /* are vertex weights being used? */
  int               using_ewgts;        /* are edge weights being used? */
  int               nsets_tot = 0;      /* total number of sets being created */
  int               igeom;              /* geometric dimension for inertial method */
  int               default_goal;       /* using default goals? */
  int               i;                  /* loop counter */
  double            seconds();
  int               reformat();
  void              free_graph(), read_params(), strout();

  if (DEBUG_TRACE > 0) {
    printf("<Entering interface>\n");
  }

  flag   = 0;
  graph  = NULL;
  coords = NULL;

  if (!Using_Main) { /* If not using main, need to read parameters file. */
    start_time  = seconds();
    params_file = fopen(PARAMS_FILENAME, "r");
    if (params_file == NULL && DEBUG_PARAMS > 1) {
      printf("Parameter file `%s' not found; using default parameters.\n", PARAMS_FILENAME);
    }
    read_params(params_file);
  }

  if (goal == NULL) { /* If not passed in, default goals have equal set sizes. */
    default_goal = TRUE;
    if (architecture == 0) {
      nsets_tot = 1 << ndims_tot;
    }
    else if (architecture == 1) {
      nsets_tot = mesh_dims[0];
    }
    else if (architecture == 2) {
      nsets_tot = mesh_dims[0] * mesh_dims[1];
    }
    else if (architecture > 2) {
      nsets_tot = mesh_dims[0] * mesh_dims[1] * mesh_dims[2];
    }

    if (MAKE_VWGTS && start != NULL) {
      vwgt_sum = start[nvtxs] - start[0] + nvtxs;
    }
    else if (vwgts == NULL) {
      vwgt_sum = nvtxs;
    }
    else {
      vwgt_sum = 0;
      vptr     = vwgts;
      for (i = nvtxs; i; i--) {
        vwgt_sum += *(vptr++);
      }
    }

    if (nsets_tot > 0) {
      vwgt_sum /= nsets_tot;
    }
    goal = smalloc_ret(nsets_tot * sizeof(double));
    if (goal == NULL) {
      strout("\nERROR: No room to make goals.\n");
      flag = 1;
      goto skip;
    }
    for (i = 0; i < nsets_tot; i++) {
      goal[i] = vwgt_sum;
    }
  }
  else {
    default_goal = FALSE;
  }

  if (MAKE_VWGTS) {
    /* Generate vertex weights equal to degree of node. */
    if (vwgts != NULL) {
      strout("WARNING: Vertex weights being overwritten by vertex degrees.");
    }
    vwgts = smalloc_ret(nvtxs * sizeof(int));
    if (vwgts == NULL) {
      strout("\nERROR: No room to make vertex weights.\n");
      flag = 1;
      goto skip;
    }
    if (start != NULL) {
      for (i = 0; i < nvtxs; i++) {
        vwgts[i] = 1 + start[i + 1] - start[i];
      }
    }
    else {
      for (i = 0; i < nvtxs; i++) {
        vwgts[i] = 1;
      }
    }
  }

  using_vwgts = (vwgts != NULL);
  using_ewgts = (ewgts != NULL);

  if (start != NULL || vwgts != NULL) { /* Reformat into our data structure. */
    time = seconds();
    flag = reformat(start, adjacency, nvtxs, &nedges, vwgts, ewgts, &graph);

    if (flag) {
      strout("\nERROR: No room to reformat graph.\n");
      goto skip;
    }

    reformat_time += seconds() - time;
  }
  else {
    nedges = 0;
  }

  if (FREE_GRAPH) { /* Free old graph data structures. */
    sfree(start);
    start = NULL;
    sfree(adjacency);
    adjacency = NULL;
    sfree(vwgts);
    vwgts = NULL;
    sfree(ewgts);
    ewgts     = NULL;
    start     = NULL;
    adjacency = NULL;
  }

  if (global_method == 3 ||
      (MATCH_TYPE == 5 && (global_method == 1 || (global_method == 2 && rqi_flag)))) {
    if (x == NULL) {
      igeom = 0;
    }
    else { /* Set up coordinate data structure. */
      coords = smalloc_ret(3 * sizeof(float *));
      if (coords == NULL) {
        strout("\nERROR: No room to make coordinate array.\n");
        flag = 1;
        goto skip;
      }
      /* Minus 1's are to allow remainder of program to index with 1. */
      coords[0] = x - 1;
      igeom     = 1;
      if (y != NULL) {
        coords[1] = y - 1;
        igeom     = 2;
        if (z != NULL) {
          coords[2] = z - 1;
          igeom     = 3;
        }
      }
    }
  }
  else {
    igeom = 0;
  }

  /* Subtract from assignment to allow code to index from 1. */
  assignment = assignment - 1;
  flag       = submain(graph, nvtxs, nedges, using_vwgts, using_ewgts, igeom, coords, outassignname,
                 outfilename, assignment, goal, architecture, ndims_tot, mesh_dims, global_method,
                 local_method, rqi_flag, vmax, ndims, eigtol, seed);

skip:
  sfree(coords);

  if (default_goal) {
    sfree(goal);
  }

  if (graph != NULL) {
    free_graph(graph);
  }

  if (flag && FREE_GRAPH) {
    sfree(start);
    start = NULL;
    sfree(adjacency);
    adjacency = NULL;
    sfree(vwgts);
    vwgts = NULL;
    sfree(ewgts);
    ewgts = NULL;
  }

  if (!Using_Main && params_file != NULL) {
    fclose(params_file);
  }

  sfree(vwgts);
  sfree(ewgts);
  return (flag);
}
