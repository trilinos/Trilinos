/*
 * Copyright(C) 1999-2020, 2023 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "params.h"
#include "smalloc.h"
#include "structs.h"
#include <math.h>

static void   avg_dists_cube(int              ndims_tot, /* total number of hypercube dimensions */
                             int              ndims,    /* number of dimensions created this step */
                             struct set_info *set_info, /* data about all the sets */
                             int              nsets,    /* number of subsets being created */
                             int              set_max,  /* largest set created so far */
                             int             *subsets,  /* subsets being created */
                             float *dists[MAXSETS] /* distances from my subsets to other sets */
  );
static void   avg_dists_mesh(int              architecture, /* dimensions of mesh */
                             struct set_info *set_info,     /* data about all the sets */
                             int              nsets,        /* number of subsets being created */
                             int              set_max,      /* largest set created so far */
                             int             *subsets,      /* subsets being created */
                             float *dists[MAXSETS] /* distances from my subsets to other sets */
  );
static double avg_dist_mesh(struct set_info *set1,        /* data about all first set */
                            struct set_info *set2,        /* data about all second set */
                            int              architecture /* dimension of mesh */
);
static double avg_dist_interval(int set1_low,  /* lowest point for first interval */
                                int set1_span, /* highest point for first interval */
                                int set2_low,  /* lowest point for second interval */
                                int set2_span  /* highest point for second interval */
);

/* Compute the terminal constraints for next partition. */
void make_term_props(struct vtx_data **graph,        /* data structure for graph */
                     int               sub_nvtxs,    /* number of vtxs in subgraph */
                     int              *loc2glob,     /* mapping from subgraph to graph */
                     int              *assignment,   /* set for each vertex */
                     int               architecture, /* 0 => hypercube, 1 => mesh */
                     int               ndims_tot,    /* total hypercube dimensions */
                     int               ndims,        /* number of dimensions at this step */
                     struct set_info  *set_info,     /* data about all the sets */
                     int               setnum,       /* number of set being divided */
                     int               nsets,        /* number of subsets being created */
                     int               set_max,      /* largest set created so far */
                     int              *subsets,      /* subsets being created */
                     float            *term_wgts[],  /* set of terminal weights for each vertex */
                     int               using_ewgts   /* are edge weights being used? */
)
{
  double term_wgt[MAXSETS]; /* terminal weights */
  float *twptr;             /* one of the term_wgts vectors */
  float *dists[MAXSETS];    /* distances from my subsets to other sets */
  float  edge_wgt;          /* weight of an edge */
  int    vtx;               /* vertex number */
  int    neighbor;          /* neighboring vertex number */
  int    neighbor_setnum;   /* set neighboring vertex is in */

  /* First compute average distance between my subsets and all other sets. */
  if (nsets <= 0) {
    return;
  }

  float *dist = smalloc(nsets * (set_max + 1) * sizeof(float));
  for (int i = 0; i < nsets; i++) {
    dists[i] = dist;
    dist += set_max + 1;
  }

  for (int k = 0; k < MAXSETS; k++) {
    term_wgt[k] = 0;
  }

  if (architecture == 0) {
    avg_dists_cube(ndims_tot, ndims, set_info, nsets, set_max, subsets, dists);
  }
  else if (architecture > 0) {
    avg_dists_mesh(architecture, set_info, nsets, set_max, subsets, dists);
  }

  edge_wgt = 1;
  for (int i = 1; i <= sub_nvtxs; i++) {
    for (int k = 1; k < nsets; k++) {
      term_wgt[k] = 0;
    }

    vtx = loc2glob[i];

    for (int j = 1; j < graph[vtx]->nedges; j++) {
      neighbor        = graph[vtx]->edges[j];
      neighbor_setnum = assignment[neighbor];
      if (neighbor_setnum != setnum) {
        if (using_ewgts) {
          edge_wgt = graph[vtx]->ewgts[j];
        }
        for (int k = 1; k < nsets; k++) {
          dist = dists[k];
          term_wgt[k] += edge_wgt * dist[neighbor_setnum];
        }
      }
    }

    for (int k = 1; k < nsets; k++) {
      twptr    = term_wgts[k];
      twptr[i] = term_wgt[k];
    }
  }

  if (dists[0] != NULL) {
    sfree(dists[0]);
  }
}

static void avg_dists_cube(int              ndims_tot, /* total number of hypercube dimensions */
                           int              ndims,     /* number of dimensions created this step */
                           struct set_info *set_info,  /* data about all the sets */
                           int              nsets,     /* number of subsets being created */
                           int              set_max,   /* largest set created so far */
                           int             *subsets,   /* subsets being created */
                           float *dists[MAXSETS]       /* distances from my subsets to other sets */
)
{
  /* First compute distances for subset 0. */
  int    myset      = subsets[0];
  float *dist0      = dists[0];
  int    ndims_left = set_info[myset].ndims;
  int    ndims_old  = ndims_tot - ndims_left - ndims;
  for (int set = 0; set < set_max; set++) {
    if (set_info[set].ndims >= 0) {
      int val = 0;
      if (ndims_left == set_info[set].ndims) {
        int start = (myset ^ set) >> ndims_old;
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
  for (int i = 1; i < nsets; i++) {
    myset       = subsets[i];
    float *dist = dists[i];

    for (int set = 0; set < set_max; set++) {
      if (set_info[set].ndims >= 0) {
        int val = 0;
        if (ndims_left == set_info[set].ndims) {
          int start = (myset ^ set) >> ndims_old;
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
                           int             *subsets,      /* subsets being created */
                           float *dists[MAXSETS] /* distances from my subsets to other sets */
)
{
  /* First compute distances for subset 0. */
  float *dist0 = dists[0];

  for (int set = 0; set < set_max; set++) {
    if (set_info[set].span[0] >= 0) {
      double val = avg_dist_mesh(&set_info[subsets[0]], &set_info[set], architecture);
      dist0[set] = val;
    }
  }

  /* Now compute all distances relative to subset 0. */
  for (int i = 1; i < nsets; i++) {
    float *dist = dists[i];
    double sep  = avg_dist_mesh(&set_info[subsets[i]], &set_info[subsets[0]], architecture);

    for (int set = 0; set < set_max; set++) {
      if (set_info[set].span[0] >= 0) {
        double val = avg_dist_mesh(&set_info[subsets[i]], &set_info[set], architecture);
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
  double val = 0;
  for (int i = 0; i < architecture; i++) {
    val += avg_dist_interval(set1->low[i], set1->span[i], set2->low[i], set2->span[i]);
  }

  return val;
}

/* Compute the average distance between two intervals */
static double avg_dist_interval(int set1_low,  /* lowest point for first interval */
                                int set1_span, /* highest point for first interval */
                                int set2_low,  /* lowest point for second interval */
                                int set2_span  /* highest point for second interval */
)
{
  double val       = 0;
  double set1_high = set1_low + set1_span - 1;
  double set1_avg  = .5 * (set1_high + set1_low);
  double set2_high = set2_low + set2_span - 1;
  double set2_avg  = .5 * (set2_high + set2_low);

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
    }
    val += (set1_high - set1_low) * (set1_high - set1_low + 1) * (set1_high - set1_low + 2) / 3.0;

    val /= set1_span * set2_span;
  }

  return (val);
}
