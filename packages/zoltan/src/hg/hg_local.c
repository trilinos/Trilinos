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


static ZOLTAN_HG_LOCAL_REF_FN local_hc;
static ZOLTAN_HG_LOCAL_REF_FN local_no;

/****************************************************************************/

ZOLTAN_HG_LOCAL_REF_FN *Zoltan_HG_Set_Local_Ref_Fn(char *str)
{
  if      (!strcasecmp(str, "hc")) return local_hc;
  else if (!strcasecmp(str, "no")) return local_no;
  else                             return NULL;
}

/****************************************************************************/

int Zoltan_HG_Local (ZZ *zz, HGraph *hg, int p, Partition part, HGPartParams *hgp)
{
  return hgp->local_ref(zz,hg,p,part,hgp->bal_tol);
}

/****************************************************************************/

static int local_no (
  ZZ *zz,
  HGraph *hg,
  int p,
  Partition part,
  float bal_tol
)
{ return ZOLTAN_OK;
}

/****************************************************************************/

static int local_hc (
  ZZ *zz, 
  HGraph *hg,
  int p,
  Partition part,
  float bal_tol
)
{ int    i, j, *in_boundary, *cut_edge, in_part,
         best_vertex, best_to;
  float  total_weight, max_weight, best_improvement;
  char   *yo="local_hc";

  if (hg->vwgt)
  { total_weight = 0.0;
    for (i=0; i<hg->nVtx; i++)
      total_weight += hg->vwgt[i];
  }
  else
    total_weight = (float)(hg->nVtx);
  max_weight = (total_weight/(float)p)*(1.0+bal_tol/100.0);
/*
  printf("total_weight:%f max_weight:%f \n",total_weight, max_weight);
*/
/*
  if (!(in_boundary = (int *) ZOLTAN_MALLOC(hg->nVtx  * sizeof(int)))  ||
      !(cut_edge    = (int *) ZOLTAN_MALLOC(hg->nEdge * sizeof(int)))   )
  { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }

  for (i=0; i<hg->nEdge; i++)
  { cut_edge[i] = 0;
    in_part = part[hg->hvertex[hg->hindex[i]]];
    for (j=hg->hindex[i]+1; j<hg->hindex[i+1]; j++)
      if (part[hg->hvertex[j]] != in_part)
      { cut_edge[i] = 1;
        break;
      }
  }
*/

  best_vertex = -1;
  best_to = -1;
  best_improvement = 0.0;
  for (i=0; i<hg->nVtx; i++)
  { 
  }




  ZOLTAN_FREE ((void **) &in_boundary);

  return ZOLTAN_OK;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

