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
static ZOLTAN_HG_GLOBAL_PART_FN global_bfs;
static ZOLTAN_HG_GLOBAL_PART_FN global_bfsr;

/****************************************************************************/

ZOLTAN_HG_GLOBAL_PART_FN *Zoltan_HG_Set_Global_Part_Fn(char *str)
{
  if      (!strcasecmp(str, "ran")) return global_ran;
  else if (!strcasecmp(str, "lin")) return global_lin;
  else if (!strcasecmp(str, "bfs")) return global_bfs;
  else if (!strcasecmp(str, "bfsr")) return global_bfsr;
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

/* BFS partitioning. Sequence partitioning with vertices 
   in breadth-first search order. */

static int global_bfs (
  ZZ *zz, 
  HGraph *hg,
  int p,
  Partition part
)
{ 
  int ierr, start, *order=NULL;
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

  /* Pick random start vertex */
  start = rand()%(hg->nVtx);

  /* Compute BFS order */
  bfs_order(zz, hg, start, order);

  /* Call sequence partitioning with random order array. */
  ierr = seq_part( zz, hg, order, p, part);

  ZOLTAN_FREE ((void **) &order);
  return (ierr);
}
 
/****************************************************************************/

/* BFS partitioning with restart. Sequence partitioning with vertices 
   in breadth-first search order, breaking of pieces as we go; 
   that is, the BFS is restarted for each partition. */

static int global_bfsr (
  ZZ *zz, 
  HGraph *hg,
  int p,
  Partition part
)
{ 
  int i, j, vtx, edge, number, start, nbor, nremaining;
  int first, last, *Q=NULL;
  float weight_avg = 0.0, weight_sum = 0.0, cutoff;
  static int srand_set = 0;
  char *yo = "global_bfsr" ;

  if (srand_set == 0)
  { srand_set = 1 ;
    srand ((unsigned long) RANDOM_SEED) ;
  }

  if (!(Q  = (int *)   ZOLTAN_MALLOC (sizeof (int) * hg->nVtx)))
  { ZOLTAN_FREE ((void **) &Q) ;
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }

  /* Compute average weight. */
  if (hg->vwgt)
  { for (i=0; i<hg->nVtx; i++)
      weight_avg += hg->vwgt[i];
  }
  else
    weight_avg = (float)hg->nVtx;
  weight_avg /= (float)p;

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
    printf("GLOBAL_PART weight_avg:%f\n",weight_avg);

  /* Sequence partitioning in BFS order. */
  /* Select random start vertex */
  start = rand()%(hg->nVtx);

/*
    Breadth first search algorithm:
    --------------------------------
    unmark all vertices with -1 
    num = 0
    sum_weight = 0;
    part[start_vtx] = num
    queue Q = start_vtx
    while Q nonempty
      choose some vertex v from front of queue
      part[v] = num
      remove v from front of queue 
      if (sum_weight>cutoff)
        num++
        empty queue
      for each hyperedge v shares
        for each unmarked neighbor w
          add w to end of queue
*/

  for (i=0; i<hg->nVtx; i++)
    part[i] = -1;  /* -1 means this vtx has not yet been assigned to a part */

  /* Initially, the queue contains only the start vertex. */
  Q[0] = start;
  first = 0;
  last = 1;
 
  number = 0; /* Assign next vertex to partition no. number */
  nremaining = hg->nVtx; /* No. of remaining vertices */
  cutoff = weight_avg;
  while (nremaining)
  {
    /* Get next vertex from queue */
    vtx = Q[first];
    first++;
    part[vtx] = number;
    weight_sum += hg->vwgt?hg->vwgt[vtx]:1.0;
    nremaining--;
    if (number+1<p && weight_sum > cutoff){
      number++;
      cutoff += weight_avg;
      for (i=first+1; i<last; i++)
        part[Q[i]] = -1; /* Remove these vertices from the queue */
      last = first+1; /* Only leave one vertex in queue */
    }
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
      printf("GLOBAL_PART i=%2d, part[%d] = %d weightsum:%f\n",i,vtx,part[vtx],weight_sum);

    /* Add nbors to queue */
    for (j=hg->vindex[vtx]; j<hg->vindex[vtx+1]; j++){
      edge = hg->vedge[j];
      for (i=hg->hindex[edge]; i<hg->hindex[edge+1]; i++){
        nbor = hg->hvertex[i];
        if (part[nbor] == -1){
          Q[last++] = nbor;
          if (last > hg->nVtx) ZOLTAN_PRINT_WARN(0, yo, "queue full");
          part[nbor] = -2; /* nbor is now in queue */
        }
      }
    }
  }

  ZOLTAN_FREE ((void **) &Q);

  return ZOLTAN_OK;
}

/****************************************************************************/

/* Compute BFS order on a hypergraph. */

static int bfs_order (
  ZZ *zz, 
  HGraph *hg,
  int start,
  int *order
)
{
  int i, j, vtx, edge, number, nbor;
  int first, last, *Q=NULL;
  static char *yo = "bfs_order";

/*
    Breadth first search algorithm:
    --------------------------------
    unmark all vertices 
    num = 0
    order[start_vtx] = num
    queue Q = start_vtx
    while Q nonempty
      choose a vertex v from front of queue
      order[v] = num++
      remove v from front of queue 
      for each hyperedge v shares
        for each unmarked neighbor w
          add w to end of queue
*/

  if (!(Q  = (int *)   ZOLTAN_MALLOC (sizeof (int) * hg->nVtx)))
  { ZOLTAN_FREE ((void **) &Q) ;
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }

  for (i=0; i<hg->nVtx; i++)
    order[i] = -1;  /* -1 means this vtx has not yet been numbered */

  /* Initially, the queue contains only the start vertex. */
  Q[0] = start;
  first = 0;
  last = 1;
 
  number = 0; /* Assign next vertex this number */
  while (number < hg->nVtx)
  {
    /* Get next vertex from queue */
    vtx = Q[first++];
    order[vtx] = number++;
     
    /* Add nbors to queue */
    for (j=hg->vindex[vtx]; j<hg->vindex[vtx+1]; j++){
      edge = hg->vedge[j];
      for (i=hg->hindex[edge]; i<hg->hindex[edge+1]; i++){
        nbor = hg->hvertex[i];
        if (order[nbor] == -1){
          Q[last++] = nbor;
          if (last > hg->nVtx) ZOLTAN_PRINT_WARN(0, yo, "Queue is full");
          order[nbor] = -2; /* nbor is now in queue */
        }
      }
    }
  }

  ZOLTAN_FREE ((void **) &Q);

  return ZOLTAN_OK;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
