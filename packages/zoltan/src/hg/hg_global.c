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
  int i, j, number, temp, *order=NULL;
  float weight_avg = 0.0, weight_sum = 0.0;
  static int srand_set = 0;
  char *yo = "global_ran" ;

  if (srand_set == 0)
  { srand_set = 1 ;
    srand ((unsigned long) RANDOM_SEED) ;
  }

  if (!(order  = (int *)   ZOLTAN_MALLOC (sizeof (int) * hg->nVtx)))
  { ZOLTAN_FREE ((void **) &order) ;
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  for (i=0; i<hg->nVtx; i++)
    order[i] = i;

  /* Randomly permute order array */
  for (i=hg->nVtx; i>0; i--)
  { number=rand()%i;
    temp = order[number];
    order[number] = order[i-1];
    order[i-1] = temp;
  }

  /* Follow the algorithm from global_lin to do
     linear partitioning on the permuted vertices. */
  if (hg->vwgt)
  { for (i=0; i<hg->nVtx; i++)
      weight_avg += hg->vwgt[i];
  }
  else
    weight_avg = (float)hg->nVtx;
  weight_avg /= (float)p;

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
    printf("GLOBAL_PART weight_avg:%f\n",weight_avg);
  number = 1; /* Assign next vertex to partition (number-1) */
  for (i=0; i<hg->nVtx; i++)
  {
    j = order[i];
    part[j] = number-1;
    weight_sum += hg->vwgt?hg->vwgt[j]:1.0;
    if (number<p && weight_sum > number*weight_avg)
      number++;
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
      printf("GLOBAL_PART i=%2d, part[%d] = %d weightsum:%f\n",i,j,part[j],weight_sum);
  }

  ZOLTAN_FREE ((void **) &order);
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

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
    printf("GLOBAL_PART weight_avg:%f\n",weight_avg);
  for (i=0; i<hg->nVtx; i++)
  { part[i] = number-1;
    weight_sum += hg->vwgt?hg->vwgt[i]:1.0;
    if (number<p && weight_sum > number*weight_avg)
      number++;
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
      printf("GLOBAL_PART part[%d] = %d weightsum:%f\n",i,part[i],weight_sum);
  }
  return ZOLTAN_OK;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
