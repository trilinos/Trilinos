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

static ZOLTAN_HG_GLOBAL_PART_FN global_ran;
static ZOLTAN_HG_GLOBAL_PART_FN global_lin;

/****************************************************************************/

ZOLTAN_HG_GLOBAL_PART_FN *Zoltan_HG_Set_Global_Part_Fn(char *str)
{
  if      (!strcasecmp(str, "ran")) return global_ran;
  else if (!strcasecmp(str, "lin")) return global_lin;
  else                              return NULL;
}

/****************************************************************************/

int Zoltan_HG_Global (ZZ *zz, HGraph *hg, int p, Partition part, HGPartParams *hgp)
{
  return hgp->global_part(zz,hg,p,part);
}

/****************************************************************************/

static int global_ran (
  ZZ *zz, 
  HGraph *hg,
  int p,
  Partition part
)
{ 
  int i, j, number, vertex, *order=NULL, new_part;
  float *weight=NULL;
  static int srand_set = 0;
  char *yo = "global_ran" ;

  if (srand_set == 0)
  { srand_set = 1 ;
    srand ((unsigned long) RANDOM_SEED) ;
  }

  if (!(order  = (int *)   ZOLTAN_MALLOC (sizeof (int) * hg->nVtx)) ||
      !(weight = (float *) ZOLTAN_CALLOC (p,sizeof(float)))    )
  { ZOLTAN_FREE ((void **) &order) ;
    ZOLTAN_FREE ((void **) &weight) ;
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  for (i=0; i<hg->nVtx; i++)
    order[i] = i;

  for (i=hg->nVtx; i>0; i--)
  { vertex = order[number=rand()%i];
    order[number] = order[i-1];
    new_part = 0;
    for (j=1; j<p; j++)
      if (weight[j] < weight[new_part])
        new_part = j;
    part[vertex] = new_part;
    weight[new_part] += (hg->vwgt)?hg->vwgt[vertex]:1.0;
  }
  ZOLTAN_FREE ((void **) &order);
  ZOLTAN_FREE ((void **) &weight);
  return ZOLTAN_OK;
}


/****************************************************************************/

static int global_lin (ZZ *zz, HGraph *hg, int p, Partition part)
{ int i, number=1;
  float weight_avg=0.0, weight_sum=0.0;

  if (hg->vwgt)
  { for (i=0; i<hg->nVtx; i++)
      weight_avg += hg->vwgt[i];
  }
  else
    weight_avg = (float)hg->nVtx;
  weight_avg /= (float)p;

  for (i=0; i<hg->nVtx; i++)
  { part[i] = number-1;
    weight_sum += hg->vwgt?hg->vwgt[i]:1.0;
    if (number<p && weight_sum > number*weight_avg)
      number++;
  }
  return ZOLTAN_OK;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
