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
static ZOLTAN_HG_GLOBAL_PART_FN global_bfsh;
static ZOLTAN_HG_GLOBAL_PART_FN global_rbfs;
static ZOLTAN_HG_GLOBAL_PART_FN global_rbfsh;

/****************************************************************************/

ZOLTAN_HG_GLOBAL_PART_FN *Zoltan_HG_Set_Global_Part_Fn(char *str)
{
  if      (!strcasecmp(str, "ran")) return global_ran;
  else if (!strcasecmp(str, "lin")) return global_lin;
  else if (!strcasecmp(str, "bfs")) return global_bfs;
  else if (!strcasecmp(str, "rbfs")) return global_rbfs;
  else if (!strcasecmp(str, "bfsh")) return global_bfsh;
  else if (!strcasecmp(str, "rbfsh")) return global_rbfsh;
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

   This function is called by global_lin and global_ran. 

   EB: This is a quick heuristic. We could alternatively use 
   a more expensive but optimal algorithm, see e.g. Ali Pinar's thesis. */

static int seq_part (ZZ *zz, HGraph *hg, int *order, int p, Partition part)
{ 
  int i, j, number;
  float weight_sum = 0.0, part_sum = 0.0, old_sum, cutoff;

  /* First sum up all the weights. */
  if (hg->vwgt){
    for (i=0; i<hg->nVtx; i++)
      weight_sum += hg->vwgt[i];
  }
  else
    weight_sum = (float)hg->nVtx;

  number = 0; /* Assign next vertex to partition no. number */
  cutoff = weight_sum/p;  /* Cutoff for current partition */
  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
    printf("GLOBAL_PART weight_sum=%f, cutoff=%f\n",weight_sum, cutoff);

  for (i=0; i<hg->nVtx; i++)
  { 
    /* If order==NULL, then use linear order. */
    j = order ? order[i] : i;
    part[j] = number;
    old_sum = part_sum;
    part_sum += hg->vwgt?hg->vwgt[j]:1.0;
    /* Check if we passed the cutoff and should start a new partition */
    if ((number+1)<p && part_sum > cutoff){
      number++;
      /* Decide if current vertex should be moved to the next partition */
      if (part_sum-cutoff > cutoff-old_sum){
        part[j]++;
        part_sum = old_sum; 
      }
      weight_sum -= part_sum;
      cutoff = weight_sum/(p-number);
      if (part[j] == number)
        part_sum = hg->vwgt?hg->vwgt[j]:1.0;
      else
        part_sum = 0.0;
    }
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
      printf("GLOBAL_PART i=%2d, part[%2d] = %2d, part_sum=%f\n",i,j,part[j],part_sum);
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
  int i, ierr, *order=NULL;
  char *yo = "global_ran" ;

  if (!(order  = (int *)   ZOLTAN_MALLOC (sizeof (int) * hg->nVtx)))
  { ZOLTAN_FREE ((void **) &order) ;
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  for (i=0; i<hg->nVtx; i++)
    order[i] = i;

  /* Randomly permute order array */
  Zoltan_HG_Rand_Perm_Int(order, hg->nVtx);
 
  /* Call sequence partitioning with random order array. */
  ierr = seq_part( zz, hg, order, p, part);

  ZOLTAN_FREE ((void **) &order);
  return (ierr);
}


/****************************************************************************/

/* Compute BFS order on a hypergraph. 
 * order[0] is the first vertex, order[1] the next, etc. 
 * Optionally, compute a partitioning based on the bfs order.
 * If p>1, restart the bfs after each new partition.
 */

static int bfs_order (
  ZZ *zz, 
  HGraph *hg,		/* Hypergraph. */
  int *order,		/* Order array. On exit, order[i] is the i'th vertex. */
  int start_vtx,	/* Start the BFS from this vertex. */
  int visit_mode,	/* Visit random (0) or heavy (1) hyperedges first? */
  int p,		/* Optional (input):  Number of partitions. */
  Partition part	/* Optional (output): Partition array. */
)
{
  int i, j, vtx, edge, bfsnumber, pnumber, nbor, next_vtx, *rank; 
  int first, last, num_edges, *edges;
  int ierr=ZOLTAN_OK;
  float weight_sum= 0.0, part_sum= 0.0, old_sum, cutoff;
  char msg[128], *mark_edge;
  static char *yo = "bfs_order";

/*
    Breadth first search algorithm:
    --------------------------------
    unmark all vertices 
    num = 0
    order[start_vtx] = num
    queue Q = start_vtx
    while Q nonempty and size < cutoff
      choose a vertex v from front of queue
      order[v] = num++
      remove v from front of queue 
      for each hyperedge v shares
        for each unmarked neighbor w
          add w to end of queue
*/

  bfsnumber = 0;  /* Assign next vertex this bfs number */
  pnumber = 0;    /* Assign next vertex this partition number */

  /* Allocate arrays. */
  if (!(rank  = (int *)   ZOLTAN_MALLOC (sizeof (int) * hg->nVtx))) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    ierr =  ZOLTAN_MEMERR;
    goto error;
  }
  for (i=0; i<hg->nVtx; i++)
    rank[i] = -1;  /* -1 means this vtx has not yet been numbered */
   
  /* array edges only needs to be of size maximum #edges for any vtx */
  num_edges = 0;
  for (i=0; i<hg->nVtx; i++)
    if (hg->vindex[i+1] - hg->vindex[i] > num_edges)
      num_edges = hg->vindex[i+1] - hg->vindex[i];
  
  if (!(mark_edge  = (char *)  ZOLTAN_CALLOC (hg->nEdge, sizeof (char))) ||
      !(edges      = (int *)   ZOLTAN_CALLOC (num_edges, sizeof (int)))) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    ierr =  ZOLTAN_MEMERR;
    goto error;
  }

  if (p){
    /* If partitioning is chosen, sum up all the weights. */
    if (hg->vwgt){
      for (i=0; i<hg->nVtx; i++)
        weight_sum += hg->vwgt[i];
    }
    else
      weight_sum = (float)hg->nVtx;

    cutoff = weight_sum/p;  /* Cutoff for current partition */
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
      printf("GLOBAL_PART weight_sum=%f, cutoff=%f\n",weight_sum, cutoff);
  }

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
    printf("Starting new BFS at vertex %d, part=%2d\n", start_vtx, p); 

  /* Use order array as a queue. Put start_vtx in queue. */
  first = bfsnumber;
  last = first+1;
  order[first] = start_vtx;
  next_vtx = (start_vtx+1 == hg->nVtx ? 0 : start_vtx+1); 

  while (bfsnumber < hg->nVtx ) {
    /* Is queue empty? */
    if (last == first){
      /* ZOLTAN_PRINT_WARN(-1, yo, "Queue is empty; hypergraph must be disconnected"); */
      /* Find an unmarked vertex to put in queue */
      while (next_vtx != start_vtx && rank[next_vtx]>=0) 
        if (++next_vtx == hg->nVtx) next_vtx = 0; /* wrap-around */
      if (next_vtx==start_vtx){
        ZOLTAN_PRINT_ERROR(-1, yo, "All vertices seem to be visited, but that cant be!");
        ierr = ZOLTAN_FATAL;
        goto error;
      }
      order[last++] = next_vtx;
    }
    /* Get next vertex from queue */
    vtx = order[first++];
    if (rank[vtx]<0)
      rank[vtx] = bfsnumber++;
    else{
      sprintf(msg, "Vertex %d in queue already labeled", vtx);
      ZOLTAN_PRINT_ERROR(-1, yo, msg);
      sprintf(msg, "bfsnumber=%d, rank[vtx] = %d", bfsnumber, rank[vtx]);
      ZOLTAN_PRINT_ERROR(-1, yo, msg);
      ierr = ZOLTAN_FATAL;
      goto error;
    }
    if (p){
      old_sum = part_sum;
      part_sum += hg->vwgt?hg->vwgt[vtx]:1.0;
      part[vtx] = pnumber;
      if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
        printf("GLOBAL_PART vtx=%2d, bfsnum=%2d, part[%2d]=%2d, part_sum=%f\n",vtx,bfsnumber-1,vtx,part[vtx],part_sum);
    }

    if (p && (pnumber+1)<p && part_sum > cutoff){
      /* Start new partition. Restart bfs. */
      pnumber++;
      /* Decide if current vertex should be moved to the next partition */
      if (part_sum-cutoff > cutoff-old_sum){
        part[vtx]++;
        part_sum = old_sum; 
      }
      weight_sum -= part_sum;
      cutoff = weight_sum/(p-pnumber);
      if (part[vtx] == pnumber)
        part_sum = hg->vwgt?hg->vwgt[vtx]:1.0;
      else
        part_sum = 0.0;
      if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
        printf("GLOBAL_PART initializing for partition %2d, cutoff = %f\n", 
               pnumber, cutoff);

      /* Clean out queue to restart bfs. */
      last = first;
      for (i=0; i<hg->nVtx; i++)
        if (rank[i] == -2) rank[i] = -1;
      for (i=0; i<hg->nEdge; i++)
        mark_edge[i] = 0;
    }
     
    /* Add nbors to queue. */
    /* Pick edges in random order, or heaviest first. */
    num_edges = hg->vindex[vtx+1] - hg->vindex[vtx];
    for (i=0; i<num_edges; i++)
      edges[i] = hg->vedge[hg->vindex[vtx]+i];
    if (visit_mode==0)
      /* Randomly permute the edges. */
      Zoltan_HG_Rand_Perm_Int(edges, num_edges);
    else if (visit_mode==1)
      /* Sort edges by weight */
      quicksort_pointer_dec_float(edges, hg->ewgt, 0, num_edges-1);

    for (j=0; j<num_edges; j++){
      edge = edges[j];
      if (!mark_edge[edge]){
        mark_edge[edge] = 1;
        for (i=hg->hindex[edge]; i<hg->hindex[edge+1]; i++){
          nbor = hg->hvertex[i];
          if (rank[nbor] == -1){
            if (last >= hg->nVtx) {
              ZOLTAN_PRINT_ERROR(-1, yo, "Queue is full");
              ierr = ZOLTAN_FATAL;
              goto error;
            }
            else{
              order[last++] = nbor;
              rank[nbor] = -2; /* nbor is now in queue */
            }
          }
        }
      }
    }
  }

  /* Sanity check: Order should be the inverse permutation of rank. */
  for (i=0; i<hg->nVtx; i++){
    if (rank[i]>=0) 
      if (order[rank[i]] != i)
         ZOLTAN_PRINT_WARN(-1, yo, "Arrays order and rank are inconsistent.");
  }

error:
  ZOLTAN_FREE ((void **) &rank);
  ZOLTAN_FREE ((void **) &mark_edge);
  ZOLTAN_FREE ((void **) &edges);
  return ierr;
}

/****************************************************************************/

/* BFS partitioning. Sequence partitioning with vertices 
   in breadth-first search order. 
   Random visit order for hyperedges. */

static int global_bfs (
  ZZ *zz, 
  HGraph *hg,
  int p,
  Partition part
)
{ 
  int i, ierr, start, *order=NULL;
  char *yo = "global_bfs" ;

  if (!(order  = (int *)   ZOLTAN_MALLOC (sizeof (int) * hg->nVtx)))
  { ZOLTAN_FREE ((void **) &order) ;
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }

  /* Find pseudo-peripheral start vertex */
  start = Zoltan_HG_Rand()%(hg->nVtx);
  for (i=0; i<2; i++){
    ierr = bfs_order(zz, hg, order, start, 0, 0, NULL);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
      goto error;
    start = order[hg->nVtx -1];
  }

  /* Compute BFS order */
  ierr = bfs_order(zz, hg, order, start, 0, 0, NULL);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
    goto error;

  /* Call sequence partitioning with BFS order array. */
  ierr = seq_part( zz, hg, order, p, part);

error:
  ZOLTAN_FREE ((void **) &order);
  return (ierr);
}
 
/****************************************************************************/

/* BFS partitioning. Sequence partitioning with vertices 
   in breadth-first search order.
   Heavy-first visit order for hyperedges. */

static int global_bfsh (
  ZZ *zz, 
  HGraph *hg,
  int p,
  Partition part
)
{ 
  int i, ierr, start, *order=NULL;
  char *yo = "global_bfsh" ;

  if (!(order  = (int *)   ZOLTAN_MALLOC (sizeof (int) * hg->nVtx)))
  { ZOLTAN_FREE ((void **) &order) ;
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }

  /* Find pseudo-peripheral start vertex */
  start = Zoltan_HG_Rand()%(hg->nVtx);
  for (i=0; i<2; i++){
    ierr = bfs_order(zz, hg, order, start, 0, 0, NULL);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
      goto error;
    start = order[hg->nVtx -1];
  }

  /* Compute BFS order */
  ierr = bfs_order(zz, hg, order, start, 1, 0, NULL);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
    goto error;

  /* Call sequence partitioning with BFS order array. */
  ierr = seq_part( zz, hg, order, p, part);

error:
  ZOLTAN_FREE ((void **) &order);
  return (ierr);
}
 
/****************************************************************************/

/* BFS partitioning with restart. Sequence partitioning with vertices 
   in breadth-first search order, breaking of pieces as we go; 
   that is, the BFS is restarted for each partition. 
   Random visit order for hyperedges. */

static int global_rbfs (
  ZZ *zz, 
  HGraph *hg,
  int p,
  Partition part
)
{ 
  int i, start, *order;
  int ierr = ZOLTAN_OK;
  char *yo = "global_rbfs" ;

  if (!(order  = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nVtx))) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    ierr = ZOLTAN_MEMERR;
    goto error;
  }

  /* Find pseudo-peripheral start vertex */
  start = Zoltan_HG_Rand()%(hg->nVtx);
  for (i=0; i<2; i++){
    ierr = bfs_order(zz, hg, order, start, 0, 0, NULL);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
      goto error;
    start = order[hg->nVtx -1];
  }

  /* Call BFS and partition with restart. */
  ierr = bfs_order(zz, hg, order, start, 0, p, part);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
    goto error;

error:
  /* Free data and return. */
  ZOLTAN_FREE ((void **) &order) ;
  return ierr;
}

/****************************************************************************/

/* BFS partitioning with restart. Sequence partitioning with vertices 
   in breadth-first search order, breaking of pieces as we go; 
   that is, the BFS is restarted for each partition. 
   Heavy-first visit order for hyperedges. */

static int global_rbfsh (
  ZZ *zz, 
  HGraph *hg,
  int p,
  Partition part
)
{ 
  int i, start, *order;
  int ierr = ZOLTAN_OK;
  char *yo = "global_rbfsh" ;

  if (!(order  = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nVtx))) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    ierr = ZOLTAN_MEMERR;
    goto error;
  }

  /* Find pseudo-peripheral start vertex */
  start = Zoltan_HG_Rand()%(hg->nVtx);
  for (i=0; i<2; i++){
    ierr = bfs_order(zz, hg, order, start, 0, 0, NULL);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
      goto error;
    start = order[hg->nVtx -1];
  }

  /* Call BFS and partition with restart. */
  ierr = bfs_order(zz, hg, order, start, 1, p, part);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
    goto error;

error:
  /* Free data and return. */
  ZOLTAN_FREE ((void **) &order) ;
  return ierr;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
