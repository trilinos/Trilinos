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

/* Scaling the weight of hypergraph edges. Currently there are 5 methods.
   The default should be number 1.
*/
int Zoltan_HG_Scale_HGraph_Weight (ZZ *zz, HGraph *hg, float *new_ewgt, int scale)
{ int    i, j;

  if (scale == 1)
  { if (hg->vwgt)
    { double sum, factor, weight;
      for (i=0; i<hg->nEdge; i++)
      { sum = factor = 0.0;
        for (j=hg->hindex[i]; j<hg->hindex[i+1]; j++)
          sum += (double)(hg->vwgt[hg->hvertex[j]]);
        for (j=hg->hindex[i]; j<hg->hindex[i+1]; j++)
        { weight = (double)(hg->vwgt[hg->hvertex[j]]);
          factor += weight*(sum-weight);
        }
        factor /= 2.0;
        if (factor <= 0.0)
          new_ewgt[i] = FLT_MAX;
        else
          new_ewgt[i] = (hg->ewgt?hg->ewgt[i]:1.0)/factor;
    } }
    else
    { int size;
      for (i=0; i<hg->nEdge; i++)
      { size = hg->hindex[i+1]-hg->hindex[i];
        new_ewgt[i] = (hg->ewgt?hg->ewgt[i]:1.0)/(float)(size*(size-1)/2);
  } } }
  else if (scale == 2)
  { for (i=0; i<hg->nEdge; i++)
    { new_ewgt[i] = (hg->ewgt?hg->ewgt[i]:1.0);
      if (hg->vwgt)
        for (j=hg->hindex[i]; j<hg->hindex[i+1]; j++)
        { if (hg->vwgt[hg->hvertex[j]] <= 0.0)
          { new_ewgt[i] = FLT_MAX;
            break;
          }
          else
            new_ewgt[i] /= (hg->vwgt[hg->hvertex[j]]);
  } }   }
  else if (scale == 3)
  { if (hg->vwgt)
    { double sum;
      for (i=0; i<hg->nEdge; i++)
      { for (j=hg->hindex[i]; j<hg->hindex[i+1]; j++)
          sum += (double)(hg->vwgt[hg->hvertex[j]]);
        if (sum <= 0.0)
          new_ewgt[i] = FLT_MAX;
        else
          new_ewgt[i] = (hg->ewgt?hg->ewgt[i]:1.0)/sum;
    } }
    else
      for (i=0; i<hg->nEdge; i++) 
        new_ewgt[i] = (hg->ewgt?hg->ewgt[i]:1.0)/
                      (float)(hg->hindex[i+1]-hg->hindex[i]);
  }
  else if (scale == 4)
  { if (hg->vwgt)
    { float max_weight;
      for (i=0; i<hg->nEdge; i++)
      { max_weight = 0.0;
        for (j=hg->hindex[i]; j<hg->hindex[i+1]; j++)
          max_weight = MAX(max_weight,(hg->vwgt[hg->hvertex[j]]));
        new_ewgt[i] = (hg->ewgt?hg->ewgt[i]:1.0)/max_weight;
    } }
    else
      for (i=0; i<hg->nEdge; i++) 
        new_ewgt[i] = (hg->ewgt?hg->ewgt[i]:1.0);
  }
  else if (scale == 5)
  { if (hg->vwgt)
    { float min_weight;
      for (i=0; i<hg->nEdge; i++)
      { min_weight = 0.0;
        for (j=hg->hindex[i]; j<hg->hindex[i+1]; j++)
          min_weight = MIN(min_weight,(hg->vwgt[hg->hvertex[j]]));
        new_ewgt[i] = (hg->ewgt?hg->ewgt[i]:1.0)/min_weight;
    } }
    else
      for (i=0; i<hg->nEdge; i++) 
        new_ewgt[i] = (hg->ewgt?hg->ewgt[i]:1.0);
  }
  return ZOLTAN_OK;
}

/****************************************************************************/

/* Scaling the weight of the edges in a graph. There are currently 4
   different methods. Default should be number 1.
*/
int Zoltan_HG_Scale_Graph_Weight (ZZ *zz, Graph *g, float *new_ewgt, int scale)
{ int   i, j;
  float vwgt_j;

  if (!g->vwgt)
    return ZOLTAN_FATAL;

  if (scale == 1)
  { for (i=0; i<g->nVtx; i++)
      for (j=g->nindex[i]; j<g->nindex[i+1]; j++)
        if (g->vwgt[i]<=0.0 || (vwgt_j=g->vwgt[g->neigh[j]])<=0.0)
          new_ewgt[j] = FLT_MAX;
        else
          new_ewgt[j] = (g->ewgt?g->ewgt[j]:1.0)/(g->vwgt[i]*vwgt_j);
  }
  else if (scale == 2)
  { for (i=0; i<g->nVtx; i++)
      for (j=g->nindex[i]; j<g->nindex[i+1]; j++)
        if (g->vwgt[i]<=0.0 || (vwgt_j=g->vwgt[g->neigh[j]])<=0.0)
          new_ewgt[j] = FLT_MAX;
        else
          new_ewgt[j] = (g->ewgt?g->ewgt[j]:1.0)/(g->vwgt[i]+vwgt_j);
  }
  else if (scale == 3)
  { for (i=0; i<g->nVtx; i++)
      for (j=g->nindex[i]; j<g->nindex[i+1]; j++)
        if (g->vwgt[i]<=0.0 || (vwgt_j=g->vwgt[g->neigh[j]])<=0.0)
          new_ewgt[j] = FLT_MAX;
        else
          new_ewgt[j] = (g->ewgt?g->ewgt[j]:1.0)/MAX(g->vwgt[i],vwgt_j);
  }
  else if (scale == 4)
  { for (i=0; i<g->nVtx; i++)
      for (j=g->nindex[i]; j<g->nindex[i+1]; j++)
        if (g->vwgt[i]<=0.0 || (vwgt_j=g->vwgt[g->neigh[j]])<=0.0)
          new_ewgt[j] = FLT_MAX;
        else if (scale == 4)
          new_ewgt[j] = (g->ewgt?g->ewgt[j]:1.0)/MIN(g->vwgt[i],vwgt_j);
  }
  return ZOLTAN_OK;
}

/****************************************************************************/

/* This is the similarity measure between two vertices in a hypergraph.
   The similarity is equal to the scaled weight of the edge in the
   transformed graph. But with this function we calculate the edge
   weights on the fly without explicitly constructing the graph.
*/

float sim (HGraph *hg, int a, int b)
{ int    i, j, edge, pins;
  float  weight, sim=0.0;

  /* First calculate the edge weight of the graph between a and b */
  for (i=hg->vindex[a]; i<hg->vindex[a+1]; i++)
  { edge = hg->vedge[i];
    j = hg->hindex[edge];
    while (j<hg->hindex[edge+1] && hg->hvertex[j]!=b)
      j++;
    if (j < hg->hindex[edge+1])
    { pins = hg->hindex[edge+1]-hg->hindex[edge];
      weight = 2.0/((pins-1)*pins);
      if (hg->ewgt)
        weight *= hg->ewgt[edge];
      sim += weight;
  } }
  /* Now scale the edge weight...currently only scaling option 1 implemented */
  if (hg->vwgt)
  { if (hg->vwgt[a]<=(float)0.0 || hg->vwgt[b]<=(float)0.0)
      sim = FLT_MAX/10.0;
    else
      sim = sim/(hg->vwgt[a]*hg->vwgt[b]);
  }
  return sim;
}

void sim_check (HGraph *hg, Graph *g)
{ int i, j;
  for (i=0; i<g->nVtx; i++)
  { for (j=g->nindex[i]; j<g->nindex[i+1]; j++)
     if ((fabs(g->ewgt[j] - sim(hg,i,g->neigh[j]))) > EPS*(g->ewgt[j]))
       printf("%d %d %.20f %.20f\n",i,g->neigh[j],g->ewgt[j],sim(hg,i,g->neigh[j]));
  }
}

/****************************************************************************/

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
