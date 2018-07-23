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
#include <stdio.h>

void balance(struct vtx_data **graph,         /* data structure for graph */
             int               nvtxs,         /* number of vertices in full graph */
             int               nedges,        /* number of edges in graph */
             int               using_vwgts,   /* are vertex weights being used? */
             int               using_ewgts,   /* are edge weights being used? */
             double *          vwsqrt,        /* sqrt of vertex weights (length nvtxs+1) */
             int               igeom,         /* geometric dimension for inertial method */
             float **          coords,        /* coordinates for inertial method */
             int *             assignment,    /* set number of each vtx (length n) */
             double *          goal,          /* desired set sizes */
             int               architecture,  /* 0=> hypercube, d=> d-dimensional mesh */
             int               ndims_tot,     /* number of cuts to make in total */
             int *             mesh_dims,     /* shape of mesh */
             int               global_method, /* global partitioning algorithm */
             int               local_method,  /* local partitioning algorithm */
             int               rqi_flag,      /* should I use multilevel eigensolver? */
             int               vmax,          /* if so, how many vertices to coarsen down to? */
             int               ndims,         /* number of eigenvectors (2^d sets) */
             double            eigtol,        /* tolerance on eigenvectors */
             int (*hops)[MAXSETS]             /* between-set hop cost for KL */
)
{
  extern int        TERM_PROP;                      /* invoking terminal propagation? */
  extern int        DEBUG_TRACE;                    /* trace the execution of the code */
  extern int        MATCH_TYPE;                     /* type of matching to use when coarsening */
  struct vtx_data **subgraph    = NULL;             /* data structure for subgraph */
  struct set_info * set_info    = NULL;             /* information about each processor subset */
  struct set_info **set_buckets = NULL;             /* buckets for sorting processor sets */
  struct set_info * set         = NULL;             /* current processor set information */
  int               hops_special[MAXSETS][MAXSETS]; /* hop mtx for nonstandard cases */
  float *           term_wgts;                      /* net pull of terminal propagation */
  float *           save_term_wgts;                 /* saved location of term_wgts */
  float *           all_term_wgts[MAXSETS];         /* net pull on all sets */
  int *             loc2glob;                       /* mapping from subgraph to graph numbering */
  int *             glob2loc;                       /* mapping from graph to subgraph numbering */
  int *             setlists;                       /* space for linked lists of vertices in sets */
  int *             list_ptrs;                      /* headers of each linked list */
  int *             degree;                         /* degrees of graph vertices from a subgraph */
  int               subsets[MAXSETS];               /* subsets being created at current step */
  int               setsize[MAXSETS];               /* sizes of sets created by division */
  double            merged_goal[MAXSETS];           /* sizes of sets at this partition level */
  double            sub_vwgt_sum;                   /* sum of subgraph vertex weights */
  int               cut_dirs[MAXDIMS];              /* directions of processor cuts if mesh */
  int               hops_flag;                      /* use normal or special hops? */
  int               ndims_real;                     /* actual dimension of partitioning */
  int               nsets_real;                     /* actual # subsets to create */
  int               maxsize;                        /* size of largest subgraph */
  int               nsets_tot;                      /* total sets to divide subgraph into */
  int               subnvtxs;                       /* number of vertices in subgraph */
  int               subnedges;                      /* number of edgess in subgraph */
  double *          subvwsqrt = NULL;               /* vwsqrt array for subgraph */
  int *             subassign = NULL;               /* set assignments for subgraph */
  float **          subcoords = NULL;               /* coordinates for subgraph */
  int               striping;                       /* partition with parallel cuts? */
  int               pass;                           /* counts passes through loop */
  int               max_proc_size;                  /* size of largest processor set */
  int               old_max_proc_size;              /* previous size of largest processor set */
  int               new_dir;                        /* new cut direction? => no terminal prop */
  int               nsets;                          /* typical number of sets at each cut */
  int               done_dir[3];                    /* mesh directions already cut */
  int               i;                              /* loop counter */

  int  make_maps(), divide_procs();
  void merge_goals();
  void divide(), make_subgraph(), make_term_props();
  void make_subvector(), make_subgeom(), remake_graph();
  void merge_assignments(), make_setlists();

  if (DEBUG_TRACE > 0) {
    printf("<Entering balance>\n");
  }

  if (global_method != 7) { /* Not read from file. */
    for (i = 1; i <= nvtxs; i++) {
      assignment[i] = 0;
    }
  }

  if (nvtxs <= 1) {
    return;
  }

  /* Compute some simple parameters. */
  save_term_wgts = NULL;
  maxsize        = 0;
  nsets          = 1 << ndims;
  subsets[2]     = 2; /* Needed for vertex separators */

  nsets_tot = 0;
  if (architecture > 0) {
    nsets_tot = mesh_dims[0] * mesh_dims[1] * mesh_dims[2];
  }
  else if (architecture == 0) {
    nsets_tot = 1 << ndims_tot;
  }

  /* Construct data structures for keeping track of processor sets. */
  set_buckets = smalloc((nsets_tot + 1) * sizeof(struct set_info *));
  set_info    = smalloc(nsets_tot * sizeof(struct set_info));
  for (i = 0; i < nsets_tot; i++) {
    set_info[i].setnum  = i;
    set_info[i].ndims   = -1;
    set_info[i].span[0] = -1;
    set_info[i].next    = NULL;
    set_buckets[i + 1]  = NULL;
  }
  set_buckets[nsets_tot] = &(set_info[0]);

  if (architecture > 0) {
    set_info[0].low[0] = set_info[0].low[1] = set_info[0].low[2] = 0;
    set_info[0].span[0]                                          = mesh_dims[0];
    set_info[0].span[1]                                          = mesh_dims[1];
    set_info[0].span[2]                                          = mesh_dims[2];
  }
  else if (architecture == 0) {
    set_info[0].ndims = ndims_tot;
  }

  done_dir[0] = done_dir[1] = done_dir[2] = FALSE;
  pass                                    = 0;
  loc2glob                                = NULL;
  degree                                  = NULL;
  glob2loc                                = NULL;
  setlists                                = NULL;
  list_ptrs                               = NULL;

  old_max_proc_size = 2 * nsets_tot;
  max_proc_size     = nsets_tot;

  while (max_proc_size > 1) { /* Some set still needs to be divided. */

    pass++;
    set                        = set_buckets[max_proc_size];
    set_buckets[max_proc_size] = set->next;

    /* Divide the processors. */
    hops_flag =
        divide_procs(architecture, ndims, ndims_tot, set_info, set, subsets, (global_method == 3),
                     &ndims_real, &nsets_real, &striping, cut_dirs, mesh_dims, hops_special);

    /* If new cut direction, turn off terminal propagation. */
    new_dir = FALSE;
    if (architecture == 0) {
      if (old_max_proc_size != max_proc_size) {
        new_dir = TRUE;
      }
    }
    else if (architecture > 0) {
      if (!done_dir[cut_dirs[0]]) {
        new_dir = TRUE;
      }
      done_dir[cut_dirs[0]] = TRUE;
    }
    old_max_proc_size = max_proc_size;

    /* Now place the new sets into the set_info data structure. */
    /* This is indexed by number of processors in set. */
    for (i = 0; i < nsets_real; i++) {
      int j = 0;
      if (architecture > 0) {
        j = set_info[subsets[i]].span[0] * set_info[subsets[i]].span[1] *
            set_info[subsets[i]].span[2];
      }
      else if (architecture == 0) {
        j = 1 << set_info[subsets[i]].ndims;
      }
      set_info[subsets[i]].next = set_buckets[j];
      set_buckets[j]            = &(set_info[subsets[i]]);
    }

    /* Construct desired set sizes for this division step. */
    if (pass == 1) { /* First partition. */
      subgraph         = graph;
      subnvtxs         = nvtxs;
      subnedges        = nedges;
      subvwsqrt        = vwsqrt;
      subcoords        = coords;
      subassign        = assignment;
      all_term_wgts[1] = NULL;

      if (!using_vwgts) {
        sub_vwgt_sum = subnvtxs;
      }
      else {
        sub_vwgt_sum = 0;
        for (i = 1; i <= subnvtxs; i++) {
          sub_vwgt_sum += subgraph[i]->vwgt;
        }
      }
      merge_goals(goal, merged_goal, set_info, subsets, nsets_real, ndims_tot, architecture,
                  mesh_dims, sub_vwgt_sum);
    }

    else { /* Not the first cut. */

      /* After first cut, allocate all space we'll need. */
      if (pass == 2) {
        glob2loc = smalloc((nvtxs + 1) * sizeof(int));
        loc2glob = smalloc((maxsize + 1) * sizeof(int));
        if (!using_vwgts) {
          subvwsqrt = NULL;
        }
        else {
          subvwsqrt = smalloc((maxsize + 1) * sizeof(double));
        }

        if (graph != NULL) {
          subgraph = smalloc((maxsize + 1) * sizeof(struct vtx_data *));
          degree   = smalloc((maxsize + 1) * sizeof(int));
        }
        else {
          subgraph = NULL;
        }
        subassign = smalloc((maxsize + 1) * sizeof(int));
        if (global_method == 3 ||
            (MATCH_TYPE == 5 && (global_method == 1 || (global_method == 2 && rqi_flag)))) {
          subcoords    = smalloc(3 * sizeof(float *));
          subcoords[1] = subcoords[2] = NULL;
          subcoords[0]                = smalloc((maxsize + 1) * sizeof(float));
          if (igeom > 1) {
            subcoords[1] = smalloc((maxsize + 1) * sizeof(float));
          }
          if (igeom > 2) {
            subcoords[2] = smalloc((maxsize + 1) * sizeof(float));
          }
        }
        else {
          subcoords = NULL;
        }
        if (TERM_PROP && graph != NULL) {
          term_wgts      = smalloc((nsets - 1) * (maxsize + 1) * sizeof(float));
          save_term_wgts = term_wgts;
          for (i = 1; i < nsets; i++) {
            all_term_wgts[i] = term_wgts;
            term_wgts += maxsize + 1;
          }
        }
      }

      /* Construct mappings between local and global vertex numbering */
      subnvtxs = make_maps(setlists, list_ptrs, set->setnum, glob2loc, loc2glob);

      if (TERM_PROP && !new_dir && graph != NULL) {
        all_term_wgts[1] = save_term_wgts;
        make_term_props(graph, subnvtxs, loc2glob, assignment, architecture, ndims_tot, ndims_real,
                        set_info, set->setnum, nsets_real, nsets_tot, subsets, all_term_wgts,
                        using_ewgts);
      }
      else {
        all_term_wgts[1] = NULL;
      }

      /* Form the subgraph in our graph format. */
      if (graph != NULL) {
        make_subgraph(graph, subgraph, subnvtxs, &subnedges, assignment, set->setnum, glob2loc,
                      loc2glob, degree, using_ewgts);
      }
      else {
        subnedges = 0; /* Otherwise some output is garbage */
      }

      if (!using_vwgts) {
        sub_vwgt_sum = subnvtxs;
      }
      else {
        sub_vwgt_sum = 0;
        for (i = 1; i <= subnvtxs; i++) {
          sub_vwgt_sum += subgraph[i]->vwgt;
        }
      }
      merge_goals(goal, merged_goal, set_info, subsets, nsets_real, ndims_tot, architecture,
                  mesh_dims, sub_vwgt_sum);

      /* Condense the relevant vertex weight array. */
      if (using_vwgts && vwsqrt != NULL) {
        make_subvector(vwsqrt, subvwsqrt, subnvtxs, loc2glob);
      }

      if (global_method == 3 ||
          (MATCH_TYPE == 5 && (global_method == 1 || (global_method == 2 && rqi_flag)))) {
        make_subgeom(igeom, coords, subcoords, subnvtxs, loc2glob);
      }
    }

    if (DEBUG_TRACE > 1) {
      printf("About to call divide with nvtxs = %d, nedges = %d, ", subnvtxs, subnedges);
      if (!architecture) {
        printf("ndims = %d\n", set->ndims);
      }
      else if (architecture == 1) {
        printf("mesh = %d\n", set->span[0]);
      }
      else if (architecture == 2) {
        printf("mesh = %dx%d\n", set->span[0], set->span[1]);
      }
      else if (architecture == 3) {
        printf("mesh = %dx%dx%d\n", set->span[0], set->span[1], set->span[2]);
      }
    }

    /* Perform a single division step. */
    divide(subgraph, subnvtxs, subnedges, using_vwgts, using_ewgts, subvwsqrt, igeom, subcoords,
           subassign, merged_goal, architecture, all_term_wgts, global_method, local_method,
           rqi_flag, vmax, ndims_real, eigtol, (hops_flag ? hops_special : hops), nsets_real,
           striping);

    /* Undo the subgraph construction. */
    if (pass != 1 && graph != NULL) {
      remake_graph(subgraph, subnvtxs, loc2glob, degree, using_ewgts);
    }

    /* Prepare for next division */
    while (max_proc_size > 1 && set_buckets[max_proc_size] == NULL) {
      --max_proc_size;
    }

    /* Merge the subgraph partitioning with the graph partitioning. */
    if (pass == 1) {
      subgraph  = NULL;
      subvwsqrt = NULL;
      subcoords = NULL;
      subassign = NULL;

      /* Find size of largest subgraph for recursing. */
      if (max_proc_size > 1) {
        for (i = 0; i < nsets; i++) {
          setsize[i] = 0;
        }
        for (i = 1; i <= nvtxs; i++) {
          ++setsize[assignment[i]];
        }
        maxsize = 0;
        for (i = 0; i < nsets; i++) {
          if (setsize[i] > maxsize) {
            maxsize = setsize[i];
          }
        }
      }

      for (i = 1; i <= nvtxs; i++) {
        assignment[i] = subsets[assignment[i]];
      }

      /* Construct list of vertices in sets for make_maps. */
      if (max_proc_size > 1) {
        setlists  = smalloc((nvtxs + 1) * sizeof(int));
        list_ptrs = smalloc(nsets_tot * sizeof(int));

        make_setlists(setlists, list_ptrs, nsets_real, subsets, assignment, loc2glob, nvtxs, TRUE);
      }
    }
    else {
      /* Construct list of vertices in sets for make_maps. */
      if (max_proc_size > 1) {
        make_setlists(setlists, list_ptrs, nsets_real, subsets, subassign, loc2glob, subnvtxs,
                      FALSE);
      }

      merge_assignments(assignment, subassign, subsets, subnvtxs, loc2glob);
    }
  }

  /* Free everything allocated for subgraphs. */
  sfree(list_ptrs);
  sfree(setlists);
  sfree(save_term_wgts);
  sfree(subassign);
  sfree(loc2glob);
  sfree(glob2loc);
  if (graph != NULL) {
    sfree(degree);
  }
  if (subgraph != NULL) {
    sfree(subgraph);
  }
  if (subvwsqrt != NULL) {
    sfree(subvwsqrt);
  }
  if (subcoords != NULL) {
    if (subcoords[0] != NULL) {
      sfree(subcoords[0]);
    }
    if (subcoords[1] != NULL) {
      sfree(subcoords[1]);
    }
    if (subcoords[2] != NULL) {
      sfree(subcoords[2]);
    }
    sfree(subcoords);
  }
  sfree(set_info);
  sfree(set_buckets);
}
