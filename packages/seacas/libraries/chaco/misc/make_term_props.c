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

#include "params.h"
#include "smalloc.h"
#include "structs.h"
#include <math.h>

static void   avg_dists_cube(), avg_dists_mesh();
static double avg_dist_mesh(), avg_dist_interval();

/* Compute the terminal constraints for next partition. */
void make_term_props(struct vtx_data **graph,        /* data structure for graph */
                     int               sub_nvtxs,    /* number of vtxs in subgraph */
                     int *             loc2glob,     /* mapping from subgraph to graph */
                     int *             assignment,   /* set for each vertex */
                     int               architecture, /* 0 => hypercube, 1 => mesh */
                     int               ndims_tot,    /* total hypercube dimensions */
                     int               ndims,        /* number of dimensions at this step */
                     struct set_info * set_info,     /* data about all the sets */
                     int               setnum,       /* number of set being divided */
                     int               nsets,        /* number of subsets being created */
                     int               set_max,      /* largest set created so far */
                     int *             subsets,      /* subsets being created */
                     float *           term_wgts[],  /* set of terminal weights for each vertex */
                     int               using_ewgts   /* are edge weights being used? */
)
{
  double term_wgt[MAXSETS]; /* terminal weights */
  float *twptr;             /* one of the term_wgts vectors */
  float *dists[MAXSETS];    /* distances from my subsets to other sets */
  float *dist;              /* one of the dists arrays */
  float  edge_wgt;          /* weight of an edge */
  int    vtx;               /* vertex number */
  int    neighbor;          /* neighboring vertex number */
  int    neighbor_setnum;   /* set neighboring vertex is in */
  int    i, j, k;           /* loop counters */

  /* First compute average distance between my subsets and all other sets. */
  dist = smalloc(nsets * (set_max + 1) * sizeof(float));
  for (i = 0; i < nsets; i++) {
    dists[i] = dist;
    dist += set_max + 1;
  }

  for (k = 0; k < MAXSETS; k++) {
    term_wgt[k] = 0;
  }

  if (architecture == 0) {
    avg_dists_cube(ndims_tot, ndims, set_info, nsets, set_max, subsets, dists);
  }
  else if (architecture > 0) {
    avg_dists_mesh(architecture, set_info, nsets, set_max, subsets, dists);
  }

  edge_wgt = 1;
  for (i = 1; i <= sub_nvtxs; i++) {
    for (k = 1; k < nsets; k++) {
      term_wgt[k] = 0;
    }

    vtx = loc2glob[i];

    for (j = 1; j < graph[vtx]->nedges; j++) {
      neighbor        = graph[vtx]->edges[j];
      neighbor_setnum = assignment[neighbor];
      if (neighbor_setnum != setnum) {
        if (using_ewgts) {
          edge_wgt = graph[vtx]->ewgts[j];
        }
        for (k = 1; k < nsets; k++) {
          dist = dists[k];
          term_wgt[k] += edge_wgt * dist[neighbor_setnum];
        }
      }
    }

    for (k = 1; k < nsets; k++) {
      twptr    = term_wgts[k];
      twptr[i] = term_wgt[k];
    }
  }

  sfree(dists[0]);
}

static void avg_dists_cube(int              ndims_tot, /* total number of hypercube dimensions */
                           int              ndims,     /* number of dimensions created this step */
                           struct set_info *set_info,  /* data about all the sets */
                           int              nsets,     /* number of subsets being created */
                           int              set_max,   /* largest set created so far */
                           int *            subsets,   /* subsets being created */
                           float *dists[MAXSETS]       /* distances from my subsets to other sets */
)
{
  float *dist0;      /* first of dists vectors */
  float *dist;       /* one of dists vectors */
  int    ndims_old;  /* hypercube dimensions not relevant */
  int    ndims_left; /* hypercube dimensions left to do */
  int    myset;      /* subset being analyzed */
  int    start;      /* bit difference between two sets */
  int    val;        /* number of differing bits */
  int    set;        /* loops through all other sets */
  int    i;          /* loop counter */

  /* First compute distances for subset 0. */
  myset      = subsets[0];
  dist0      = dists[0];
  ndims_left = set_info[myset].ndims;
  ndims_old  = ndims_tot - ndims_left - ndims;
  for (set = 0; set < set_max; set++) {
    if (set_info[set].ndims >= 0) {
      val = 0;
      if (ndims_left == set_info[set].ndims) {
        start = (myset ^ set) >> ndims_old;
        while (start) {
          if (start & 1) {
            val++;
          }
          start >>= 1;
        }
      }
      dist0[set] = val;
    }
  }

  /* Now compute all distances relative to subset 0. */
  for (i = 1; i < nsets; i++) {
    myset = subsets[i];
    dist  = dists[i];

    for (set = 0; set < set_max; set++) {
      if (set_info[set].ndims >= 0) {
        val = 0;
        if (ndims_left == set_info[set].ndims) {
          start = (myset ^ set) >> ndims_old;
          while (start) {
            if (start & 1) {
              val++;
            }
            start >>= 1;
          }
        }
        /* Note: this is net preference for set over set 0. */
        dist[set] = dist0[set] - val;
      }
    }
  }
}

static void avg_dists_mesh(int              architecture, /* dimensions of mesh */
                           struct set_info *set_info,     /* data about all the sets */
                           int              nsets,        /* number of subsets being created */
                           int              set_max,      /* largest set created so far */
                           int *            subsets,      /* subsets being created */
                           float *dists[MAXSETS] /* distances from my subsets to other sets */
)
{
  float *dist0; /* first of dists vectors */
  float *dist;  /* one of dists vectors */
  double val;   /* distance from subset to set */
  double sep;   /* distance between two subsets */
  int    set;   /* loops through all other sets */
  int    i;     /* loop counter */

  /* First compute distances for subset 0. */
  dist0 = dists[0];

  for (set = 0; set < set_max; set++) {
    if (set_info[set].span[0] >= 0) {
      val        = avg_dist_mesh(&set_info[subsets[0]], &set_info[set], architecture);
      dist0[set] = val;
    }
  }

  /* Now compute all distances relative to subset 0. */
  for (i = 1; i < nsets; i++) {
    dist = dists[i];
    sep  = avg_dist_mesh(&set_info[subsets[i]], &set_info[subsets[0]], architecture);

    for (set = 0; set < set_max; set++) {
      if (set_info[set].span[0] >= 0) {
        val = avg_dist_mesh(&set_info[subsets[i]], &set_info[set], architecture);
        /* Note: this is net preference for set over 0. */
        dist[set] = (dist0[set] - val) / sep;
      }
    }
  }
}

/* Compute the average distance between two subsets of mesh processors. */
static double avg_dist_mesh(struct set_info *set1,        /* data about all first set */
                            struct set_info *set2,        /* data about all second set */
                            int              architecture /* dimension of mesh */
)
{
  double val; /* distance returned */
  int    i;   /* loop counter */
  double avg_dist_interval();

  val = 0;

  for (i = 0; i < architecture; i++) {
    val += avg_dist_interval(set1->low[i], set1->span[i], set2->low[i], set2->span[i]);
  }

  return (val);
}

/* Compute the average distance between two intervals */
static double avg_dist_interval(int set1_low,  /* lowest point for first interval */
                                int set1_span, /* highest point for first interval */
                                int set2_low,  /* lowest point for second interval */
                                int set2_span  /* highest point for second interval */
)
{
  double set1_high; /* length of first interval */
  double set1_avg;  /* average value in first interval */
  double set2_high; /* length of second interval */
  double set2_avg;  /* average value in second interval */
  double val;       /* average distance between intervals */

  val       = 0;
  set1_high = set1_low + set1_span - 1;
  set1_avg  = .5 * (set1_high + set1_low);
  set2_high = set2_low + set2_span - 1;
  set2_avg  = .5 * (set2_high + set2_low);

  if (set1_low > set2_high || set2_low > set1_high) {
    val = fabs(set1_avg - set2_avg);
  }

  else {
    if (set1_high > set2_high) {
      val += .5 * (set2_high - set2_low + 1) * (set1_high - set2_high) * (set1_high - set2_low + 1);
      set1_high = set2_high;
    }
    else if (set2_high > set1_high) {
      val += .5 * (set1_high - set1_low + 1) * (set2_high - set1_high) * (set2_high - set1_low + 1);
      set2_high = set1_high;
    }
    if (set1_low < set2_low) {
      val += .5 * (set2_high - set2_low + 1) * (set2_low - set1_low) * (set2_high - set1_low + 1);
      set1_low = set2_low;
    }
    else if (set2_low < set1_low) {
      val += .5 * (set1_high - set1_low + 1) * (set1_low - set2_low) * (set1_high - set2_low + 1);
      set2_low = set1_low;
    }
    val += (set1_high - set1_low) * (set1_high - set1_low + 1) * (set1_high - set1_low + 2) / 3.0;

    val /= set1_span * set2_span;
  }

  return (val);
}
