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

/* Zoltan_HG_Global computes a global partitioning of a hypergraph.
 * Typically, this routine is called at the bottom level in a 
 * multilevel scheme (V-cycle).
 */ 

int Zoltan_HG_Global (ZZ *zz, HGraph *hg, int p, Partition part, HGPartParams *hgp)
{
  return hgp->global_part(zz,hg,p,part);
}

/****************************************************************************/

/* Sequence partitioning on the vertices of a hypergraph
   in a given order. Currently, even partition sizes
   are assumed. Multi-weights are not yet supported. 

   This function is called by global_lin and global_rand. 

   EB: This is a quick heuristic. We could alternatively 
   use a more expensive but optimal algorithm, see paper by Ali Pinar. */

static int seq_part (ZZ *zz, HGraph *hg, int *order, int p, Partition part)
{ 
  int i, j, number;
  float weight_avg = 0.0, weight_sum = 0.0, old_sum, cutoff;

  /* First sum up all the weights and find average. */
  if (hg->vwgt){
    for (i=0; i<hg->nVtx; i++)
      weight_avg += hg->vwgt[i];
  }
  else
    weight_avg = (float)hg->nVtx;
  weight_avg /= (float)p;

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
    printf("GLOBAL_PART weight_avg:%f\n",weight_avg);
  number = 0; /* Assign next vertex to partition no. number */
  cutoff = weight_avg; /* Cutoff for current partition */
  for (i=0; i<hg->nVtx; i++)
  {
    if (order)
      j = order[i];
    else
      j = i;   /* If order==NULL, then use linear order. */
    part[j] = number;
    old_sum = weight_sum;
    weight_sum += hg->vwgt?hg->vwgt[j]:1.0;
    /* Check if we passed the cutoff and should start a new partition */
    /* cutoff = (number+1)*weight_avg for uniform partition sizes */
    if ((number+1)<p && weight_sum > cutoff){
      /* Check if current vertex should be moved to the next partition */
      if (weight_sum-cutoff > cutoff-old_sum)
        part[j]++;
      number++;
      cutoff += weight_avg;
    }
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
      printf("GLOBAL_PART i=%2d, part[%d] = %d weightsum:%f\n",i,j,part[j],weight_sum);
  }

  return ZOLTAN_OK;
}

/****************************************************************************/

/* Linear partitioning. Sequence partitioning with vertices in linear order. */

static int global_lin (ZZ *zz, HGraph *hg, int p, Partition part)
{ 
  /* Call sequence partitioning with no order array. */
  return seq_part( zz, hg, NULL, p, part);
}

/****************************************************************************/

/* Random partitioning. Sequence partitioning with vertices in random order. */
static int global_ran (
  ZZ *zz, 
  HGraph *hg,
  int p,
  Partition part
)
{ 
  int i, ierr, number, temp, *order=NULL;
  static int srand_set=0;
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

  /* Call sequence partitioning with random order array. */
  ierr = seq_part( zz, hg, order, p, part);

  ZOLTAN_FREE ((void **) &order);
  return (ierr);
}


/****************************************************************************/

#if 0 /* EB: This is work in progress. */

/* BFS partitioning. Sequence partitioning with vertices in breadth-first search 
   order. */
static int global_bfs (
  ZZ *zz, 
  HGraph *hg,
  int p,
  Partition part
)
{ 
  int i, j, number, start, *order=NULL;
  float weight_avg = 0.0, weight_sum = 0.0;
  static int srand_set = 0;
  static int bfs_order();
  char *yo = "global_bfs" ;

  if (srand_set == 0)
  { srand_set = 1 ;
    srand ((unsigned long) RANDOM_SEED) ;
  }

  if (!(order  = (int *)   ZOLTAN_MALLOC (sizeof (int) * hg->nVtx)))
  { ZOLTAN_FREE ((void **) &order) ;
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }

  /* Compute bfs order from random vertex */
  start = rand()%(hg->nVtx);
  bfs_order(zz, hg, start, order);

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


/* Compute a BFS order for hg starting at start_vtx. */

static int bfs_order(
  ZZ *zz, 
  HGraph *hg,
  int start_vtx,
  int *order
)
{
/*
    Breadth first search algorithm:

    label all vertices with -1 (unmark)
    num = 0
    label start_vtx with num++
    list L = start_vtx
    while L nonempty
      choose some vertex v from front of list
      visit v
      remove v from front of list 
      for each hyperedge v shares
      for each unmarked neighbor w
        label w with num++
        add w to end of list
*/
  return ZOLTAN_OK;
}
#endif

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
