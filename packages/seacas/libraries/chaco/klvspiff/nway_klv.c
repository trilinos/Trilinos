/*
 * Copyright(C) 1999-2023 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "defs.h"
#include "smalloc.h"
#include "structs.h"
#include <math.h>
#include <stdio.h>
#include <sys/types.h>

/*
   Keep guys moved in and guys moving out of separator.

   To restore, move all undesirable guys back.
     (1) guys moved out of separator get put back in.
     (2) guys moved into separator (bspace) get put back.
         Note: Should be done in this order.
         Note: No neighbors need be considered.

     (3) To clear dvals, I should compute touch all guys that
         were ever in separator (bspace) and their neighbors.
*/

int nway_klv(struct vtx_data **graph,      /* data structure for graph */
             int               nvtxs,      /* number of vtxs in graph */
             struct bilist   **lbuckets,   /* array of lists for bucket sort */
             struct bilist   **rbuckets,   /* array of lists for bucket sort */
             struct bilist    *llistspace, /* list data structure for each vertex */
             struct bilist    *rlistspace, /* list data structure for each vertex */
             int              *ldvals,     /* d-values for each transition */
             int              *rdvals,     /* d-values for each transition */
             int              *sets,       /* processor each vertex is assigned to */
             int               maxdval,    /* maximum d-value for a vertex */
             double           *goal,       /* desired set sizes */
             int               max_dev,    /* largest allowed deviation from balance */
             int             **bndy_list,  /* list of vertices on boundary (0 ends) */
             double           *weightsum   /* sum of vweights in each set (in and out) */
)
{
  extern double kl_bucket_time; /* time spent in KL bucketsort */
  extern int    KL_BAD_MOVES;   /* # bad moves in a row to stop KL */
  extern int    DEBUG_KL;       /* debug flag for KL */
  extern int    KL_NTRIES_BAD;  /* number of unhelpful passes before quitting */
  extern int    KL_MAX_PASS;    /* maximum # outer KL loops */

  int nbadtries = KL_NTRIES_BAD;

  int enforce_balance      = FALSE;
  int enforce_balance_hard = FALSE;

  double total_weight = goal[0] + goal[1];

  int *bspace = smalloc_ret((nvtxs + 1) * sizeof(int));

  if (bspace == NULL) {
    return (1);
  }

  int *bdy_ptr     = *bndy_list;
  int  list_length = 0;
  while (*bdy_ptr != 0) {
    bspace[list_length++] = *bdy_ptr++;
  }

  sfree(*bndy_list);

  clear_dvals(graph, nvtxs, ldvals, rdvals, bspace, list_length);

  int step_cutoff = KL_BAD_MOVES;
  int cost_cutoff = maxdval * step_cutoff / 7;
  if (cost_cutoff < step_cutoff) {
    cost_cutoff = step_cutoff;
  }

  double partial_weight = weightsum[0] + weightsum[1];
  double ratio          = partial_weight / total_weight;
  double delta0         = fabs(weightsum[0] - goal[0] * ratio);
  double delta1         = fabs(weightsum[1] - goal[1] * ratio);
  int    balanced =
      (delta0 + delta1 <= max_dev) && weightsum[0] != total_weight && weightsum[1] != total_weight;

  double bestg_min = -2.0 * nvtxs * maxdval;
  int    parity    = FALSE;
  int    nbad      = 0;
  int    npass     = 0;
  int    improved  = 0;
  int    done      = FALSE;
  while (!done) {
    npass++;
    int    ever_balanced = FALSE;
    double balance_best  = delta0 + delta1;

    /* Initialize various quantities. */
    int ltop = 2 * maxdval;
    int rtop = 2 * maxdval;

    int            gtotal     = 0;
    double         bestg      = bestg_min;
    int            beststep   = -1;
    int            bestlength = list_length;
    struct bilist *out_list   = NULL;

    int neg_steps = 0;

    /* Compute the initial d-values, and bucket-sort them. */
    double time = seconds();
    bucketsortsv(graph, nvtxs, lbuckets, rbuckets, llistspace, rlistspace, ldvals, rdvals, sets,
                 maxdval, parity, bspace, list_length);
    parity = !parity;
    kl_bucket_time += seconds() - time;

    if (DEBUG_KL > 2) {
      printf("After sorting, left buckets:\n");
      p1bucket(lbuckets, llistspace, maxdval);
      printf("              right buckets:\n");
      p1bucket(rbuckets, rlistspace, maxdval);
    }

    /* Now determine the set of vertex moves. */

    int step = 1;
    for (;; step++) {

      /* Find the highest d-value in each set. */
      /* But only consider moves from large to small sets, or moves */
      /* in which balance is preserved. */
      /* Break ties in some nonarbitrary manner. */
      int bestval = -maxdval - 1;

      partial_weight    = weightsum[0] + weightsum[1];
      ratio             = partial_weight / total_weight;
      int left_too_big  = (weightsum[0] > (goal[0] + .5 * max_dev) * ratio);
      int right_too_big = (weightsum[1] > (goal[1] + .5 * max_dev) * ratio);

      while (ltop >= 0 && lbuckets[ltop] == NULL) {
        --ltop;
      }

      int to         = -1; /* sets moving into / out of */
      int bestvtx    = -1; /* best vertex to move */
      int weightfrom = 0;

      double left_imbalance = 0.0; /* imbalance if I move to the left */
      if (ltop >= 0 && !left_too_big) {
        int lvtx       = ((size_t)lbuckets[ltop] - (size_t)llistspace) / sizeof(struct bilist);
        int lweight    = graph[lvtx]->vwgt;
        int rweight    = lweight - (ltop - maxdval);
        weightfrom     = rweight;
        to             = 0;
        bestvtx        = lvtx;
        bestval        = ltop - maxdval;
        partial_weight = weightsum[0] + lweight + weightsum[1] - rweight;
        ratio          = partial_weight / total_weight;
        left_imbalance = max(fabs(weightsum[0] + lweight - goal[0] * ratio),
                             fabs(weightsum[1] - rweight - goal[1] * ratio));
      }

      while (rtop >= 0 && rbuckets[rtop] == NULL) {
        --rtop;
      }
      if (rtop >= 0 && !right_too_big) {
        int rvtx       = ((size_t)rbuckets[rtop] - (size_t)rlistspace) / sizeof(struct bilist);
        int rweight    = graph[rvtx]->vwgt;
        int lweight    = rweight - (rtop - maxdval);
        partial_weight = weightsum[0] - lweight + weightsum[1] + rweight;
        ratio          = partial_weight / total_weight;
        double right_imbalance = max(fabs(weightsum[0] - lweight - goal[0] * ratio),
                                     fabs(weightsum[1] + rweight - goal[1] * ratio));
        if (rtop - maxdval > bestval || (rtop - maxdval == bestval &&
                                         (right_imbalance < left_imbalance ||
                                          (right_imbalance == left_imbalance && drandom() < .5)))) {
          to         = 1;
          weightfrom = lweight;
          bestvtx    = rvtx;
          bestval    = rtop - maxdval;
        }
      }

      if (bestval == -maxdval - 1) { /* No allowed moves */
        if (DEBUG_KL > 0) {
          printf("No KLV moves at step %d.  bestg = %g at step %d.\n", step, bestg, beststep);
        }
        break;
      }

      int            *to_dvals;       /* d-values I'm moving to */
      int            *from_dvals;     /* d-values I'm moving from */
      int            *to_top;         /* ptr to top of set moving to */
      struct bilist  *to_listspace;   /* list structure I'm moving to */
      struct bilist  *from_listspace; /* list structure I'm moving from */
      struct bilist **to_buckets;     /* buckets I'm moving to */
      struct bilist **from_buckets;   /* buckets I'm moving from */
      int             from;
      if (to == 0) {
        from           = 1;
        to_listspace   = llistspace;
        from_listspace = rlistspace;
        to_dvals       = ldvals;
        from_dvals     = rdvals;
        to_buckets     = lbuckets;
        from_buckets   = rbuckets;
        to_top         = &ltop;
      }
      else {
        from           = 0;
        to_listspace   = rlistspace;
        from_listspace = llistspace;
        to_dvals       = rdvals;
        from_dvals     = ldvals;
        to_buckets     = rbuckets;
        from_buckets   = lbuckets;
        to_top         = &rtop;
      }

      int vweight = graph[bestvtx]->vwgt;

      weightsum[to] += vweight;
      weightsum[from] -= weightfrom;

      /* Check if this partition is balanced. */
      partial_weight    = weightsum[0] + weightsum[1];
      ratio             = partial_weight / total_weight;
      delta0            = fabs(weightsum[0] - goal[0] * ratio);
      delta1            = fabs(weightsum[1] - goal[1] * ratio);
      int temp_balanced = (delta0 + delta1 <= max_dev) && weightsum[0] != total_weight &&
                          weightsum[1] != total_weight;
      ever_balanced      = (ever_balanced || temp_balanced);
      double balance_val = delta0 + delta1;

      gtotal += bestval;

      if ((gtotal > bestg && temp_balanced) ||
          (enforce_balance_hard && balance_val < balance_best)) {
        bestg    = gtotal;
        beststep = step;
        if (balance_val < balance_best) {
          balance_best = balance_val;
        }
        if (temp_balanced) {
          enforce_balance_hard = FALSE;
        }
      }

      /* Monitor the stopping criteria. */
      if (bestval < 0) {
        if (!enforce_balance || ever_balanced) {
          neg_steps++;
        }
        int neg_cost; /* decrease in sum of d-values */
        if (bestg != bestg_min) {
          neg_cost = bestg - gtotal;
        }
        else {
          neg_cost = -maxdval - 1;
        }
        if ((neg_steps > step_cutoff || neg_cost > cost_cutoff) &&
            !(enforce_balance && bestg == bestg_min)) {
          if (DEBUG_KL > 0) {
            if (neg_steps > step_cutoff) {
              printf("KLV step cutoff at step %d.  bestg = %g at step %d.\n", step, bestg,
                     beststep);
            }
            else if (neg_cost > cost_cutoff) {
              printf("KLV cost cutoff at step %d.  bestg = %g at step %d.\n", step, bestg,
                     beststep);
            }
          }
          weightsum[to] -= vweight;
          weightsum[from] += weightfrom;
          break;
        }
      }
      else if (bestval > 0) {
        neg_steps = 0;
      }

      /* Remove vertex from its buckets, and flag it as finished. */
      sets[bestvtx] = to;
      removebilist(&to_listspace[bestvtx], &to_buckets[bestval + maxdval]);
      /*
                  printf("After to removebilist\n");
                  p1bucket(to_buckets, to_listspace, maxdval);
      */

      if (from_dvals[bestvtx] != -maxdval - 1) {
        removebilist(&from_listspace[bestvtx], &from_buckets[from_dvals[bestvtx] + maxdval]);
        /*
                        printf("After from removebilist\n");
                        p1bucket(from_buckets, from_listspace, maxdval);
        */
      }
      from_dvals[bestvtx] = -maxdval - 1;

      /* Now keep track of vertices moved out of separator so */
      /* I can restore them as needed. */
      llistspace[bestvtx].next = out_list;
      out_list                 = &(llistspace[bestvtx]);

      /* Now update the d-values of all the neighbors */
      /* And neighbors of neighbors ... */

      /* If left move:
         1. Separator neighbors right gain => infinity
         2. Left neighbors unaffected.
         3. Right neighbors move into separator.
            A. Right gain = infinity.
            B. Left gain = computed.
            C. For any of their neighbors in separator increase left gain.
      */

      int *edges = graph[bestvtx]->edges;
      for (int j = graph[bestvtx]->nedges - 1; j; j--) {
        int neighbor = *(++edges);

        int group = sets[neighbor];

        if (group == 2) { /* In separator. */
          int gain = from_dvals[neighbor] + maxdval;
          /* Gain in the from direction => -infinity */
          if (gain >= 0) {
            removebilist(&from_listspace[neighbor], &from_buckets[gain]);
            /*
                                    printf("\n  After removing %d\n", neighbor);
                                    p1bucket(from_buckets, from_listspace, maxdval);
            */
            from_dvals[neighbor] = -maxdval - 1;
          }
        }
        else if (group == from) {
          /* Gain in the from direction => -infinity */
          sets[neighbor]       = 2;
          from_dvals[neighbor] = -maxdval - 1;

          if (to == 0) {
            bspace[list_length++] = -neighbor;
          }
          else {
            bspace[list_length++] = neighbor;
          }

          int *edges2 = graph[neighbor]->edges;
          int  vwgt   = graph[neighbor]->vwgt;
          int  gain   = graph[neighbor]->vwgt;
          int  flag   = FALSE;
          for (int k = graph[neighbor]->nedges - 1; k; k--) {
            int neighbor2 = *(++edges2);
            int group2    = sets[neighbor2];
            if (group2 == 2) {
              int dval = to_dvals[neighbor2] + maxdval;
              if (dval >= 0) {
                movebilist(&to_listspace[neighbor2], &to_buckets[dval], &to_buckets[dval + vwgt]);
                /*
                                                printf("\n  After moving %d from bucket %d to bucket
                   %d\n", neighbor2, dval, dval + vwgt);
                                                p1bucket(to_buckets, to_listspace, maxdval);
                */
                to_dvals[neighbor2] += vwgt;
                dval += vwgt;
                if (dval > *to_top) {
                  *to_top = dval;
                }
              }
            }
            else if (group2 == from) {
              gain -= graph[neighbor2]->vwgt;
              if (to_dvals[neighbor2] + maxdval < 0) {
                flag = TRUE;
              }
            }
          }

          if (flag) { /* Not allowed to move further. */
            to_dvals[neighbor] = -maxdval - 1;
          }
          else {
            to_dvals[neighbor] = gain;
            /* place in appropriate bucket */

            gain += maxdval;
            add2bilist(&to_listspace[neighbor], &to_buckets[gain]);
            /*
                                    printf("\nAfter adding %d to bucket %d\n", neighbor, gain -
               maxdval);
                                    p1bucket(to_buckets, to_listspace, maxdval);
            */

            if (gain > *to_top) {
              *to_top = gain;
            }
          }
        }
      }
      if (beststep == step) {
        bestlength = list_length;
      }
      if (DEBUG_KL > 2) {
        printf("\n-- After step, left buckets:\n");
        p1bucket(lbuckets, llistspace, maxdval);
        printf("             right buckets:\n");
        p1bucket(rbuckets, rlistspace, maxdval);
      }
    }

    /* Done with a pass; should we actually perform any swaps? */
    if (bestg > 0 || (bestg != bestg_min && !balanced && enforce_balance)) {
      improved += bestg;
    }
    else {
      if (enforce_balance_hard) {
        /* I've done the best I can, give up. */
        done = TRUE;
      }
      if (enforce_balance) {
        enforce_balance_hard = TRUE;
      }
      enforce_balance = TRUE;
      nbad++;
    }

    /* Work backwards, undoing all the undesirable moves. */

    /* First reset vertices moved out of the separator. */
    if (out_list) {
      if (beststep < 0) {
        beststep = 0;
      }
      for (int i = step - 1; i > beststep; i--) {
        int vtx = ((size_t)out_list - (size_t)llistspace) / sizeof(struct bilist);
        if (sets[vtx] != 2) {
          weightsum[sets[vtx]] -= graph[vtx]->vwgt;
        }
        sets[vtx] = 2;
        out_list  = out_list->next;
      }
    }

    for (int i = list_length - 1; i >= bestlength; i--) {
      int vtx = bspace[i];
      if (vtx < 0) {
        if (sets[-vtx] == 2) {
          weightsum[1] += graph[-vtx]->vwgt;
        }
        sets[-vtx] = 1;
      }
      else {
        if (sets[vtx] == 2) {
          weightsum[0] += graph[vtx]->vwgt;
        }
        sets[vtx] = 0;
      }
    }

    partial_weight = weightsum[0] + weightsum[1];
    ratio          = partial_weight / total_weight;
    delta0         = fabs(weightsum[0] - goal[0] * ratio);
    delta1         = fabs(weightsum[1] - goal[1] * ratio);
    balanced       = (delta0 + delta1 <= max_dev) && weightsum[0] != total_weight &&
               weightsum[1] != total_weight;

    done = done || (nbad >= nbadtries && balanced);
    if (KL_MAX_PASS > 0) {
      done = done || (npass == KL_MAX_PASS && balanced);
    }

    if (!done) { /* Rezero dval values. */
      clear_dvals(graph, nvtxs, ldvals, rdvals, bspace, list_length);
    }

    /* Construct list of separator vertices to pass to buckets or return */
    list_length = make_sep_list(bspace, list_length, sets);

    if (done) {
      bspace[list_length] = 0;
      bspace              = srealloc(bspace, (list_length + 1) * sizeof(int));
      *bndy_list          = bspace;
    }

    /*
    int gain = 0;
    int j = 0;
    int k = 0;
    for (int i = 1; i <= nvtxs; i++) {
      if (sets[i] == 0) {
        j += graph[i]->vwgt;
      }
      else if (sets[i] == 1) {
        k += graph[i]->vwgt;
      }
      else if (sets[i] == 2) {
        gain += graph[i]->vwgt;
      }
    }
            printf("\nAfter pass of KLV: sets = %d/%d, sep = %d  (bestg = %g)\n\n\n",
                   j, k, gain, bestg);
    */
  }

  if (DEBUG_KL > 0) {
    printf("   KLV required %d passes to improve by %d.\n", npass, improved);
  }

  return (0);
}
