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

static int split_hypergraph (HGraph*, HGraph*, Partition, int, ZZ*);



/* recursively divides problem into 2 parts until all p found */
int Zoltan_HG_rdivide (int lo, int hi, Partition final, ZZ *zz, HGraph *hg,
 HGPartParams *hgp, int level)
{
  int i, mid;
  int err;
  Partition part;
  HGraph *new;
  char *yo = "Zoltan_HG_rdivide";

  hg->redl = 2;  /* hg->redl gets changed during execution, reset it */
   
  /* only one part remaining, record results and exit */
  if (lo == hi) {
    for (i = 0; i < hg->nVtx; i++)
      final [hg->vmap[i]] = lo - 1;
    return ZOLTAN_OK;
  }

  part = (Partition) ZOLTAN_MALLOC (hg->nVtx * sizeof (int));
  if (part == NULL) {
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Unable to allocate memory.");
    return ZOLTAN_MEMERR;
  }

  /* bipartition current hypergraph with appropriate split ratio */
  mid = (lo+hi)/2;
  hg->ratio = (double) (mid-lo+1) / (double) (hi-lo+1);
  hg->redl = 0;
  err = Zoltan_HG_HPart_Lib (zz, hg, 2, part, hgp, level);
  if (err != ZOLTAN_OK)
    return err;

  /* if only two parts total, record results and exit */
  if (lo + 1 == hi)  {
    for (i = 0; i < hg->nVtx; i++)
      final [hg->vmap[i]] = ((part[i] == 0) ? (lo-1) : (hi-1));
    ZOLTAN_FREE (&part);
    return ZOLTAN_OK;
  }

  new = (HGraph*) ZOLTAN_MALLOC (sizeof (HGraph));
  if (new == NULL)  {
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Unable to allocate memory.");
    return ZOLTAN_MEMERR;
  }
  Zoltan_HG_HGraph_Init (new);

  /* recursively divide in two parts and repartition hypergraph */
  err = split_hypergraph (hg, new, part, 0, zz);
  if (err != ZOLTAN_OK)
    return err;
    
  err = Zoltan_HG_rdivide (lo, mid, final, zz, new, hgp, level+1);
  Zoltan_HG_HGraph_Free (new);
  if (err != ZOLTAN_OK)
    return err;

  err = split_hypergraph (hg, new, part, 1, zz);
  if (err != ZOLTAN_OK)
    return err;
    
  err = Zoltan_HG_rdivide (mid+1, hi, final, zz, new, hgp, level+1);
  Zoltan_HG_HGraph_Free (new);
  if (err != ZOLTAN_OK)
    return err;

  /* remove alloc'ed structs */
  ZOLTAN_FREE (&part);
  ZOLTAN_FREE (&new);

  return ZOLTAN_OK;
}



static int split_hypergraph (HGraph *old, HGraph *new, Partition part,
 int partid, ZZ *zz)
{
  int *tmap;        /* temporary array mapping from old HGraph info to new */
  int edge, i;                                            /* loop counters */
  int pins;
  char *yo = "split_hypergraph";

  /* allocate memory for dynamic arrays in new HGraph and for tmap array */
  new->vmap = (int*) ZOLTAN_MALLOC (old->nVtx * sizeof (int));
  tmap      = (int*) ZOLTAN_MALLOC (old->nVtx * sizeof (int));
  if (new->vmap == NULL || tmap == NULL)  {
    Zoltan_Multifree (__FILE__, __LINE__, 2, &new->vmap, &tmap);
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Unable to allocate memory 1.");
    return ZOLTAN_MEMERR;
  }

  /* save vertex and edge weights if they exist */
  if (old->vwgt)
    new->vwgt = (float*) ZOLTAN_MALLOC (old->nVtx  * sizeof(float)
     * old->VertexWeightDim);
  if (old->ewgt)
    new->ewgt = (float*) ZOLTAN_MALLOC (old->nEdge * sizeof(float)
     * old->EdgeWeightDim);

  /* fill in tmap array, -1 for ignored vertices, otherwise nonnegative int */
  new->nVtx = 0;
    for (i = 0; i < old->nVtx; i++)
      if (part[i] == partid)  {
        tmap[i] = new->nVtx;
        new->vmap[new->nVtx] = old->vmap[i];
        if (new->vwgt)
          new->vwgt[new->nVtx] = old->vwgt[i];
        new->nVtx++;
      }
      else
        tmap[i] = -1;

  /* continue allocating memory for dynamic arrays in new HGraph */
  new->vmap    = (int*) ZOLTAN_REALLOC (new->vmap, new->nVtx * sizeof (int));
  new->hindex  = (int*) ZOLTAN_MALLOC ((old->nEdge+1) * sizeof (int));
  new->hvertex = (int*) ZOLTAN_MALLOC (old->nInput * sizeof (int));
  if (new->vmap == NULL || new->hindex == NULL || new->hvertex == NULL)  {
    Zoltan_Multifree (__FILE__, __LINE__, 5, &new->vmap, &new->hindex,
     &new->hvertex, &new->vmap, &tmap);
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Unable to allocate memory 2.");
    return ZOLTAN_MEMERR;
  }

  /* fill in hindex and hvertex arrays in new HGraph */
  new->nEdge  = 0;
  new->nInput = 0;
  for (edge = 0; edge < old->nEdge; edge++)  {
    pins = new->nInput;      /* not yet sure if this edge will be included */
    new->hindex[new->nEdge] = new->nInput;
    for (i = old->hindex[edge]; i < old->hindex[edge+1]; i++)
      if (tmap [old->hvertex[i]] >= 0)  {
        new->hvertex[new->nInput] = tmap[old->hvertex[i]];
        new->nInput++;        /* edge has at least one vertex in partition */
      }
    if (pins < new->nInput)  {
      if (new->ewgt)
        new->ewgt[new->nEdge] = old->ewgt[edge];
      new->nEdge++;
    }
  }
  new->hindex[new->nEdge] = new->nInput;
  new->info               = old->info;
  new->VertexWeightDim    = old->VertexWeightDim;
  new->EdgeWeightDim      = old->EdgeWeightDim;

  /* shrink hindex, hvertex arrays to correct size & determine vindex, vedge */
  new->hvertex = (int*) ZOLTAN_REALLOC(new->hvertex, sizeof(int) * new->nInput);
  new->hindex = (int*) ZOLTAN_REALLOC(new->hindex, sizeof(int) * (new->nEdge+1));
  Zoltan_HG_Create_Mirror (zz, new);

  ZOLTAN_FREE (&tmap);
  return ZOLTAN_OK;
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
