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

int Zoltan_HG_Scale_Graph_Weight (ZZ *zz, Graph *g, float *new_ewgt)
{ int   i, j;

  if (!g->vwgt)
    return ZOLTAN_FATAL;

  for (i=0; i<g->nVtx; i++)
    for (j=g->nindex[i]; j<g->nindex[i+1]; j++)
    { if (g->vwgt[i]<=0.0 || g->vwgt[g->neigh[j]]<=0.0)
        new_ewgt[j] = FLT_MAX;
      else
        new_ewgt[j] = (g->ewgt?g->ewgt[j]:1.0)/g->vwgt[i]/g->vwgt[g->neigh[j]];
    }
  return ZOLTAN_OK;
}

/****************************************************************************/

int Zoltan_HG_Scale_HGraph_Weight (ZZ *zz, HGraph *hg, float *new_ewgt)
{ int   i, j;
  float weight, sum1, sum2;

  if (!hg->vwgt)
    return ZOLTAN_FATAL;

  for (i=0; i<hg->nEdge; i++)
  { sum1 = sum2 = 0.0;
    if (hg->vwgt)
    { for (j=hg->hindex[i]; j<hg->hindex[i+1]; j++)
      { weight = hg->vwgt[hg->hvertex[j]];
        sum1 += weight;
        sum2 += weight*weight;
    } }
    else
      sum1 = sum2 = (float)(hg->hindex[i+1]-hg->hindex[i]);
    sum1 = (sum1*sum1-sum2)/2.0;
    if (sum1 == 0.0)
      new_ewgt[i] = FLT_MAX;
    else
      new_ewgt[i] = (hg->ewgt?hg->ewgt[i]:1.0)/sum1;
  }
   
  return ZOLTAN_OK;
}

/****************************************************************************/

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
