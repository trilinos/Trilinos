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

static ZOLTAN_HG_LOCAL_REF_FN local_no;
static ZOLTAN_HG_LOCAL_REF_FN local_hc;

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
{ int    i, j, side, vertex, edge, *cut[2], *boundary[2], boundary_n[2],
         best_vertex;
  float  total_weight, max_weight, improvement, best_improvement,
         part_weight[2];
  char   *yo="local_hc";

  if (p != 2)
  { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "p!=2 not yet implemented for local_hc.");
    return ZOLTAN_FATAL;
  }

  if (hg->nEdge == 0)
    return ZOLTAN_OK;

  if (hg->vwgt)
  { total_weight = 0.0;
    for (i=0; i<hg->nVtx; i++)
      total_weight += hg->vwgt[i];
  }
  else
    total_weight = (float)(hg->nVtx);
  max_weight = (total_weight/(float)p)*bal_tol;
/*
  printf("total_weight:%f max_weight:%f \n",total_weight, max_weight);
*/
  part_weight[0] = part_weight[1] = 0.0;
  for (i=0; i<hg->nVtx; i++)
    part_weight[part[i]] += (hg->vwgt?hg->vwgt[i]:1.0);
/*
  printf("weights: %f %f\n",part_weight[0],part_weight[1]);
*/
  if (!(boundary[0] = (int *) ZOLTAN_MALLOC(2*hg->nVtx* sizeof(int))) || 
      !(cut[0]      = (int *) ZOLTAN_CALLOC(2*hg->nEdge,sizeof(int)))  )
  { ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  boundary[1] = &(boundary[0][hg->nVtx]);
  cut[1] = &(cut[0][hg->nEdge]);

  for (i=0; i<hg->nEdge; i++)
  { cut[0][i] = cut[1][i] = 0;
    for (j=hg->hindex[i]; j<hg->hindex[i+1]; j++)
      (cut[part[hg->hvertex[j]]][i])++;
  }

  do
  { boundary_n[0] = boundary_n[1] = 0;
    for (i=0; i<hg->nVtx; i++)
    { for (j=hg->vindex[i]; j<hg->vindex[i+1]; j++)
        if (cut[0][hg->vedge[j]]==1 || cut[1][hg->vedge[j]]==1)
        { boundary[part[i]][boundary_n[part[i]]++] = i;
          break;
        }
    }
/*
    printf("boundaries:%d %d\n",boundary_n[0],boundary_n[1]);
*/
    best_vertex = -1;
    best_improvement = 0.0;

    for (side=0; side<2; side++)
    { for (i=0; i<boundary_n[side]; i++)
      { vertex = boundary[side][i];
        if (part_weight[1-part[vertex]]+(hg->vwgt?hg->vwgt[vertex]:1.0) < max_weight)
        { improvement = 0.0;
          for (j=hg->vindex[vertex]; j<hg->vindex[vertex+1]; j++)
          { edge = hg->vedge[j];
            if (cut[part[vertex]][edge] == 1)
              improvement += (hg->ewgt?(hg->ewgt[edge]):1.0);
            else if (cut[1-part[vertex]][edge] == 0)
              improvement -= (hg->ewgt?(hg->ewgt[edge]):1.0); 
          }

          if (improvement > best_improvement)
          { best_vertex = vertex;
            best_improvement = improvement;
          }
        }
      }
    }
/*
    printf("best: %d %f %d\n",best_vertex, best_improvement,part[best_vertex]);
*/
    if (best_improvement > 0.0)
    { if (best_vertex < 0)
        return ZOLTAN_FATAL;
/*
      for (i=hg->vindex[best_vertex]; i<hg->vindex[best_vertex+1]; i++)
      { edge = hg->vedge[i];
        printf("%d %f ",edge,hg->ewgt?hg->ewgt[edge]:1.0);
        for (j=hg->hindex[edge]; j<hg->hindex[edge+1]; j++)
          printf("%d ",part[hg->hvertex[j]]);
        for (j=hg->hindex[edge]; j<hg->hindex[edge+1]; j++)
          printf("%d ",hg->hvertex[j]);
        printf("\n");
      }
*/

      for (i=hg->vindex[best_vertex]; i<hg->vindex[best_vertex+1]; i++)
      { cut[part[best_vertex]][hg->vedge[i]] --;
        cut[1-part[best_vertex]][hg->vedge[i]] ++;
      }
      part_weight[part[best_vertex]] -= (hg->vwgt?hg->vwgt[best_vertex]:1.0);
      part[best_vertex] = 1-part[best_vertex];
      part_weight[part[best_vertex]] += (hg->vwgt?hg->vwgt[best_vertex]:1.0);
    }
  } while (best_improvement > 0.0);

  ZOLTAN_FREE ((void **) &(boundary[0]));
  ZOLTAN_FREE ((void **) &(cut[0]));

  return ZOLTAN_OK;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

