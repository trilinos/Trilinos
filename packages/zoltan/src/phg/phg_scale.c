/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "phg.h"
#include <float.h>



/****************************************************************************/
/* Scaling the weight of hypergraph edges. Currently there are 5 methods.
   The default should be number 1. */
int Zoltan_PHG_Scale_HGraph_Weight (ZZ *zz, PHGraph *hg, float *new_ewgt,
 int scale)
{
  int i, j;
  static char *yo = "Zoltan_PHG_Scale_HGraph_Weight";

  if (scale == 1) {
    if (hg->vwgt) {
      double sum, factor, weight;
      for (i = 0; i < hg->nEdge; i++) {
        sum = factor = 0.0;
        for (j = hg->hindex[i]; j < hg->hindex[i+1]; j++)
          sum += (double) hg->vwgt[hg->hvertex[j]];
        for (j = hg->hindex[i]; j < hg->hindex[i+1]; j++) {
          weight = (double) hg->vwgt[hg->hvertex[j]];
          factor += (weight * (sum-weight));
        }
        factor /= 2.0;
        if (factor <= 0.0)
          new_ewgt[i] = FLT_MAX/10.0;
        else
          new_ewgt[i] = (hg->ewgt ? hg->ewgt[i] : 1.0) / factor;
      }
    }
    else {
      int size;
      for (i = 0; i < hg->nEdge; i++) {
        size = hg->hindex[i+1] - hg->hindex[i];
        new_ewgt[i] =(hg->ewgt ? hg->ewgt[i] : 1.0)/(double)(size*(size-1)/2);
      }
    }
  }
  else if (scale == 2) {
    for (i = 0; i < hg->nEdge; i++) {
      new_ewgt[i] = (hg->ewgt ? hg->ewgt[i] : 1.0);
      if (hg->vwgt)
        for (j = hg->hindex[i]; j < hg->hindex[i+1]; j++) {
          if (hg->vwgt[hg->hvertex[j]] <= 0.0) {
            new_ewgt[i] = FLT_MAX/10.0;
            break;
          }
          else
            new_ewgt[i] /= (hg->vwgt[hg->hvertex[j]]);
        }
     }
  }
  else if (scale == 3) {
    if (hg->vwgt) {
      double sum ;
      for (i = 0; i < hg->nEdge; i++) {
        sum = 0.0;
        for (j = hg->hindex[i]; j < hg->hindex[i+1]; j++)
          sum += (double)(hg->vwgt[hg->hvertex[j]]);
        if (sum <= 0.0)
          new_ewgt[i] = FLT_MAX/10.0;
        else
          new_ewgt[i] = (hg->ewgt ? hg->ewgt[i] : 1.0) / sum;
      }
    }
    else
      for (i = 0; i < hg->nEdge; i++)
        new_ewgt[i] = (hg->ewgt ? hg->ewgt[i] : 1.0)/
        (double)(hg->hindex[i+1] - hg->hindex[i]);
  }
  else if (scale == 4) {
    if (hg->vwgt) {
      double max_weight;
      for (i = 0; i < hg->nEdge; i++) {
        max_weight = 0.0;
        for (j = hg->hindex[i]; j < hg->hindex[i+1]; j++)
          max_weight = MAX (max_weight, hg->vwgt[hg->hvertex[j]]);
        new_ewgt[i] = (hg->ewgt ? hg->ewgt[i] : 1.0) / max_weight;
      }
     }
     else
       for (i = 0; i < hg->nEdge; i++)
          new_ewgt[i] = hg->ewgt ? hg->ewgt[i] : 1.0;
     }
  else if (scale == 5) {
     if (hg->vwgt) {
        double min_weight;
        for (i = 0; i < hg->nEdge; i++) {
           min_weight = 0.0;
           for (j = hg->hindex[i]; j < hg->hindex[i+1]; j++)
              min_weight = MIN(min_weight, hg->vwgt[hg->hvertex[j]]);
           new_ewgt[i] = (hg->ewgt ? hg->ewgt[i] : 1.0) / min_weight;
           }
        }
     else
        for (i = 0; i < hg->nEdge; i++)
           new_ewgt[i] = hg->ewgt ? hg->ewgt[i] : 1.0;
     }
  else {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid PHG_REDUCTION_SCALING");
     return ZOLTAN_FATAL;
  }

  return ZOLTAN_OK;
}

/****************************************************************************/



/* Scaling the weight of the edges in a graph. There are currently 4
   different methods. Default should be number 1. */
int Zoltan_PHG_Scale_Graph_Weight (ZZ *zz, PGraph *g, float *new_ewgt, int scale)
{
int   i, j;
double vwgt_j;

  if (!g->vwgt)
     return ZOLTAN_FATAL;

  if (scale == 1) {
     for (i = 0; i < g->nVtx; i++)
        for (j = g->nindex[i]; j < g->nindex[i+1]; j++)
           if (g->vwgt[i] <= 0.0 || (vwgt_j = g->vwgt[g->neigh[j]]) <= 0.0)
              new_ewgt[j] = FLT_MAX/10.0;
           else
              new_ewgt[j] = (g->ewgt ? g->ewgt[j] : 1.0)/(g->vwgt[i] * vwgt_j);
     }
  else if (scale == 2) {
     for (i = 0; i < g->nVtx; i++)
        for (j = g->nindex[i]; j < g->nindex[i+1]; j++)
           if (g->vwgt[i] <= 0.0 || (vwgt_j = g->vwgt[g->neigh[j]]) <= 0.0)
              new_ewgt[j] = FLT_MAX / 10.0;
           else
              new_ewgt[j] = (g->ewgt ? g->ewgt[j] : 1.0)/(g->vwgt[i] + vwgt_j);
     }
  else if (scale == 3) {
    for (i = 0; i < g->nVtx; i++)
      for (j = g->nindex[i]; j < g->nindex[i+1]; j++)
        if (g->vwgt[i] <= 0.0 || (vwgt_j = g->vwgt[g->neigh[j]]) <= 0.0)
          new_ewgt[j] = FLT_MAX/10.0;
        else
          new_ewgt[j] = (g->ewgt ? g->ewgt[j] : 1.0)/MAX(g->vwgt[i],vwgt_j);
  }
  else if (scale == 4) {
    for (i = 0; i < g->nVtx; i++)
      for (j = g->nindex[i]; j < g->nindex[i+1]; j++)
        if (g->vwgt[i] <= 0.0 || (vwgt_j = g->vwgt[g->neigh[j]]) <= 0.0)
          new_ewgt[j] = FLT_MAX/10.0;
        else if (scale == 4)
          new_ewgt[j] = (g->ewgt ? g->ewgt[j] : 1.0)/MIN(g->vwgt[i],vwgt_j);
  }
  return ZOLTAN_OK;
}

/****************************************************************************/



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
