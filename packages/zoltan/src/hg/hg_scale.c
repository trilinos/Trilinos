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

#include "hypergraph.h"

/****************************************************************************/

int Zoltan_HG_Scale_Graph_Weight (ZZ *zz, Graph *g, float *new_ewgt, int scale)
{ int   i, j;
  float vwgt_i, vwgt_j;

  if (!g->vwgt)
    return ZOLTAN_FATAL;

  for (i=0; i<g->nVtx; i++)
  { vwgt_i = g->vwgt[i];
    for (j=g->nindex[i]; j<g->nindex[i+1]; j++)
    { vwgt_j = g->vwgt[g->neigh[j]];
      if (vwgt_i<=0.0 || vwgt_j<=0.0)
        new_ewgt[j] = FLT_MAX;
      else if (scale == 1)
        new_ewgt[j] = (g->ewgt?g->ewgt[j]:1.0)/(vwgt_i*vwgt_j);
      else if (scale == 2)
        new_ewgt[j] = (g->ewgt?g->ewgt[j]:1.0)/(vwgt_i+vwgt_j);
      else if (scale == 3)
        new_ewgt[j] = (g->ewgt?g->ewgt[j]:1.0)/MAX(vwgt_i,vwgt_j);
      else if (scale == 4)
        new_ewgt[j] = (g->ewgt?g->ewgt[j]:1.0)/MIN(vwgt_i,vwgt_j);
  } }
  return ZOLTAN_OK;
}

/****************************************************************************/

int Zoltan_HG_Scale_HGraph_Weight (ZZ *zz, HGraph *hg, float *new_ewgt)
{ int    i, j;
  double weight, sum, scale;

  if (!hg->vwgt)
    return ZOLTAN_FATAL;

  for (i=0; i<hg->nEdge; i++)
  { scale = sum = 0.0;
    if (hg->vwgt)
    { for (j=hg->hindex[i]; j<hg->hindex[i+1]; j++)
        sum += (double)(hg->vwgt[hg->hvertex[j]]);
      for (j=hg->hindex[i]; j<hg->hindex[i+1]; j++)
      { weight = (double)(hg->vwgt[hg->hvertex[j]]);
        scale += weight*(sum-weight);
      }
      scale /= 2.0;
    }
    else
      scale = (double)(hg->hindex[i+1]-hg->hindex[i]);

    if (scale == 0.0)
      new_ewgt[i] = FLT_MAX;
    else
      new_ewgt[i] = (hg->ewgt?hg->ewgt[i]:1.0)/scale;
  }

  return ZOLTAN_OK;
}

/****************************************************************************/

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
