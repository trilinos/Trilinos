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
static ZOLTAN_HG_GLOBAL_PART_FN global_gr0;
static ZOLTAN_HG_GLOBAL_PART_FN global_gr1;
static ZOLTAN_HG_GLOBAL_PART_FN global_gr2;
static ZOLTAN_HG_GLOBAL_PART_FN global_gr3;
static ZOLTAN_HG_GLOBAL_PART_FN global_gr4;


/****************************************************************************/

ZOLTAN_HG_GLOBAL_PART_FN *Zoltan_HG_Set_Global_Part_Fn(char *str)
{
  if      (!strcasecmp(str, "ran"))   return global_ran;
  else if (!strcasecmp(str, "lin"))   return global_lin;
  else if (!strcasecmp(str, "bfs"))   return global_bfs;
  else if (!strcasecmp(str, "rbfs"))  return global_rbfs;
  else if (!strcasecmp(str, "bfsh"))  return global_bfsh;
  else if (!strcasecmp(str, "rbfsh")) return global_rbfsh;
  else if (!strcasecmp(str, "gr0"))   return global_gr0;
  else if (!strcasecmp(str, "gr1"))   return global_gr1;
  else if (!strcasecmp(str, "gr2"))   return global_gr2;
  else if (!strcasecmp(str, "gr3"))   return global_gr3;
  else if (!strcasecmp(str, "gr4"))   return global_gr4;
  else                                return NULL;
}

/****************************************************************************/

/* Zoltan_HG_Global computes a global partitioning of a hypergraph.
 * Typically, this routine is called at the bottom level in a
 * multilevel scheme (V-cycle).
 */

int Zoltan_HG_Global(ZZ *zz,HGraph *hg,int p,Partition part,HGPartParams *hgp)
{
  return hgp->global_part(zz,hg,p,part,hgp);
}

/****************************************************************************/

/* Sequence partitioning on the vertices of a hypergraph
   in a given order. Currently, even partition sizes
   are assumed. Multi-weights are not yet supported.

   This function is called by global_lin and global_ran.

   EBEB: This is a quick heuristic. We could alternatively use
   a more expensive but optimal algorithm, see e.g. Ali Pinar's thesis. */

static int seq_part (
  ZZ *zz, 
  HGraph *hg, 
  int *order, 
  int p, 
  Partition part,
  HGPartParams *hgp
)
{
  int i, j, number;
  double weight_sum = 0.0, part_sum = 0.0, old_sum, cutoff;

  /* First sum up all the weights. */
  if (hg->vwgt){
    for (i=0; i<hg->nVtx; i++)
      weight_sum += hg->vwgt[i];
  }
  else
    weight_sum = (double) hg->nVtx;

  number = 0; /* Assign next vertex to partition no. number */
  cutoff = weight_sum/p;  /* Cutoff for current partition */
  if (hgp->output_level >= HG_DEBUG_ALL)
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
    if (hgp->output_level >= HG_DEBUG_ALL)
      printf("GLOBAL_PART i=%2d, part[%2d] = %2d, part_sum=%f\n",i,j,part[j],
       part_sum);
  }

  return ZOLTAN_OK;
}

/****************************************************************************/

/* Linear partitioning. Sequence partitioning with vertices in linear order. */

static int global_lin (
  ZZ *zz, 
  HGraph *hg, 
  int p, 
  Partition part, 
  HGPartParams *hgp
)
{
  /* Call sequence partitioning with no order array. */
  return seq_part( zz, hg, NULL, p, part, hgp);
}

/****************************************************************************/

/* Random partitioning. Sequence partitioning with vertices in random order. */
static int global_ran (
  ZZ *zz,
  HGraph *hg,
  int p,
  Partition part,
  HGPartParams *hgp
)
{
  int i, ierr, *order=NULL;
  char *yo = "global_ran";

  if (!(order  = (int *)   ZOLTAN_MALLOC (hg->nVtx*sizeof(int))))
  { ZOLTAN_FREE ((void **) &order);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  for (i=0; i<hg->nVtx; i++)
    order[i] = i;

  /* Randomly permute order array */
  Zoltan_HG_Rand_Perm_Int(order, hg->nVtx);

  /* Call sequence partitioning with random order array. */
  ierr = seq_part( zz, hg, order, p, part, hgp);

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
  Partition part,	/* Optional (output): Partition array. */
  HGPartParams *hgp     /* Partitioning parameters. */
)
{
  int i, j, vtx, edge, bfsnumber, pnumber, nbor, next_vtx, *rank = NULL;
  int first, last, num_edges, *edges = NULL;
  int ierr=ZOLTAN_OK;
  double weight_sum= 0.0, part_sum= 0.0, old_sum, cutoff;
  char msg[128], *mark_edge = NULL;
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
  if (!(rank  = (int *)   ZOLTAN_MALLOC (hg->nVtx*sizeof(int)))) {
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


  mark_edge  = (char *)  ZOLTAN_CALLOC (hg->nEdge, sizeof (char));
  edges      = (int *)   ZOLTAN_CALLOC (num_edges, sizeof (int));
  if ((hg->nEdge > 0 && mark_edge == NULL) || (num_edges > 0 && edges == NULL)){
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
      weight_sum = (double) hg->nVtx;

    cutoff = weight_sum/p;  /* Cutoff for current partition */
    if (hgp->output_level >= HG_DEBUG_ALL)
      printf("GLOBAL_PART weight_sum=%f, cutoff=%f\n",weight_sum, cutoff);
  }

  if (hgp->output_level >= HG_DEBUG_ALL)
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
        ZOLTAN_PRINT_ERROR(-1, yo,
         "All vertices seem to be visited, but that cant be!");
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
      if (hgp->output_level >= HG_DEBUG_ALL)
        printf("GLOBAL_PART vtx=%2d, bfsnum=%2d, part[%2d]=%2d, part_sum=%f\n",
         vtx,bfsnumber-1,vtx,part[vtx],part_sum);
    }

    if (p && (pnumber+1)<p && part_sum > cutoff){
      /* Start new partition. Restart bfs. */
      pnumber++;
      /* Decide if current vertex should be moved to the next partition */
      if (part_sum-cutoff > cutoff-old_sum){
        part[vtx]++;
        part_sum = old_sum;
        if (hgp->output_level >= HG_DEBUG_ALL)
          printf("GLOBAL_PART vtx=%2d, bfsnum=%2d, part[%2d]=%2d\n",
           vtx,bfsnumber-1,vtx,part[vtx]);
      }
      weight_sum -= part_sum;
      cutoff = weight_sum/(p-pnumber);
      if (part[vtx] == pnumber)
        part_sum = hg->vwgt?hg->vwgt[vtx]:1.0;
      else
        part_sum = 0.0;
      if (hgp->output_level >= HG_DEBUG_ALL)
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
      Zoltan_quicksort_pointer_dec_float(edges, hg->ewgt, 0, num_edges-1);

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
  Partition part,
  HGPartParams *hgp
)
{
  int i, ierr, start, *order=NULL;
  char *yo = "global_bfs";

  if (!(order  = (int *)   ZOLTAN_MALLOC (hg->nVtx*sizeof(int))))
  { ZOLTAN_FREE ((void **) &order);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }

  /* Find pseudo-peripheral start vertex */
  /* EBEB: Make this a function that can be called
     each time we begin a new connected component. */
  start = Zoltan_HG_Rand()%(hg->nVtx);
  for (i=0; i<2; i++){
    ierr = bfs_order(zz, hg, order, start, 0, 0, NULL, hgp);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
      goto error;
    start = order[hg->nVtx -1];
  }

  /* Compute BFS order */
  ierr = bfs_order(zz, hg, order, start, 0, 0, NULL, hgp);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
    goto error;

  /* Call sequence partitioning with BFS order array. */
  ierr = seq_part( zz, hg, order, p, part, hgp);

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
  Partition part,
  HGPartParams *hgp
)
{
  int i, ierr, start, *order=NULL;
  char *yo = "global_bfsh";

  if (!(order  = (int *)   ZOLTAN_MALLOC (hg->nVtx*sizeof(int))))
  { ZOLTAN_FREE ((void **) &order);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }

  /* Find pseudo-peripheral start vertex */
  start = Zoltan_HG_Rand()%(hg->nVtx);
  for (i=0; i<2; i++){
    ierr = bfs_order(zz, hg, order, start, 0, 0, NULL, hgp);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
      goto error;
    start = order[hg->nVtx -1];
  }

  /* Compute BFS order */
  ierr = bfs_order(zz, hg, order, start, 1, 0, NULL, hgp);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
    goto error;

  /* Call sequence partitioning with BFS order array. */
  ierr = seq_part( zz, hg, order, p, part, hgp);

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
  Partition part,
  HGPartParams *hgp
)
{
  int i, start, *order;
  int ierr = ZOLTAN_OK;
  char *yo = "global_rbfs";

  if (!(order  = (int *) ZOLTAN_MALLOC (hg->nVtx*sizeof(int)))) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    ierr = ZOLTAN_MEMERR;
    goto error;
  }

  /* Find pseudo-peripheral start vertex */
  start = Zoltan_HG_Rand()%(hg->nVtx);
  for (i=0; i<2; i++){
    ierr = bfs_order(zz, hg, order, start, 0, 0, NULL, hgp);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
      goto error;
    start = order[hg->nVtx -1];
  }

  /* Call BFS and partition with restart. */
  ierr = bfs_order(zz, hg, order, start, 0, p, part, hgp);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
    goto error;

error:
  /* Free data and return. */
  ZOLTAN_FREE ((void **) &order);
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
  Partition part,
  HGPartParams *hgp
)
{
  int i, start, *order;
  int ierr = ZOLTAN_OK;
  char *yo = "global_rbfsh";

  if (!(order  = (int *) ZOLTAN_MALLOC (hg->nVtx*sizeof(int)))) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    ierr = ZOLTAN_MEMERR;
    goto error;
  }

  /* Find pseudo-peripheral start vertex */
  start = Zoltan_HG_Rand()%(hg->nVtx);
  for (i=0; i<2; i++){
    ierr = bfs_order(zz, hg, order, start, 0, 0, NULL, hgp);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
      goto error;
    start = order[hg->nVtx -1];
  }

  /* Call BFS and partition with restart. */
  ierr = bfs_order(zz, hg, order, start, 1, p, part, hgp);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
    goto error;

error:
  /* Free data and return. */
  ZOLTAN_FREE ((void **) &order);
  return ierr;
}

/*********************************************************************/

/* Greedy ordering/partitioning based on a priority function
   for selecting vertices. A heap is used as a priority queue. */

static int greedy_order (
  ZZ *zz,
  HGraph *hg,		/* Hypergraph. */
  int *order,		/* Order array. On exit, order[i] is the i'th vertex. */
  int start_vtx,	/* Start the ordering from this vertex. */
  int priority_mode,	/* Priority mode for selecting vertices */
  int p,		/* Optional (input):  Number of partitions. */
  Partition part,	/* Optional (output): Partition array. */
  HGPartParams *hgp     /* Partitioning parameters. */
)
{
  int i, j, vtx, edge, bfsnumber, pnumber, nbor, *rank;
  int esize, *vtx_count=NULL, *visited=NULL, *cut[2];
  int ierr=ZOLTAN_OK;
  double weight_sum= 0.0, part_sum= 0.0, old_sum, cutoff;
  double *gain = NULL, *edge_sum = NULL, delta;
  double damp_factor;
  char msg[128];
  HEAP h[2];
  static char *yo = "greedy_order";

  bfsnumber = 0;  /* Assign next vertex this bfs number */
  pnumber = 0;    /* Assign next vertex this partition number */

  /* Allocate arrays. */
  if (!(rank  = (int *)     ZOLTAN_CALLOC (hg->nVtx, sizeof (int))) ||
      !(gain  = (double *)   ZOLTAN_CALLOC (hg->nVtx, sizeof (double))) ) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    ierr =  ZOLTAN_MEMERR;
    goto error;
  }
  for (i=0; i<hg->nVtx; i++)
    rank[i] = -1;  /* -1 means this vtx has not yet been numbered */

  if (priority_mode && (!(priority_mode&1))){ /* 2,4,6,... */
    if (!(edge_sum = (double *)  ZOLTAN_CALLOC (hg->nVtx, sizeof (double)))){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ierr =  ZOLTAN_MEMERR;
      goto error;
    }
    /* Sum up edge weights incident to each vertex. */
    for (edge=0; edge<hg->nEdge; edge++){
      for (i=hg->hindex[edge]; i<hg->hindex[edge+1]; i++){
        edge_sum[hg->hvertex[i]] += (hg->ewgt ? hg->ewgt[edge] : 1.0);
      }
    }
  }
  if (priority_mode==0){
    cut[0] = (int *)  ZOLTAN_CALLOC (2*hg->nEdge, sizeof (int));
    visited = (int *)  ZOLTAN_CALLOC (hg->nVtx, sizeof (int));
    if ((hg->nEdge > 0 && cut[0] == NULL) || visited == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ierr =  ZOLTAN_MEMERR;
      goto error;
    }
    cut[1] = &(cut[0][hg->nEdge]);
    /* Initialize cut values. */
    for (i=0; i<hg->nEdge; i++)
      for (j=hg->hindex[i]; j<hg->hindex[i+1]; j++)
        (cut[visited[hg->hvertex[j]]][i])++;
    /* Initialize gain values. */
    for (i=0; i<hg->nVtx; i++){
      for (j=hg->vindex[i]; j<hg->vindex[i+1]; j++){
        edge = hg->vedge[j];
        gain[i] -= (hg->ewgt?(hg->ewgt[edge]):1.0);
      }
    }
  }
  else
    cut[0] = cut[1] = NULL;

  if (priority_mode>=3){
    if (!(vtx_count = (int *)  ZOLTAN_CALLOC (hg->nEdge, sizeof (int)))){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ierr =  ZOLTAN_MEMERR;
      goto error;
    }
  }

  if (p){
    /* If partitioning is chosen, sum up all the weights. */
    if (hg->vwgt){
      for (i=0; i<hg->nVtx; i++)
        weight_sum += hg->vwgt[i];
    }
    else
      weight_sum = (double) hg->nVtx;

    cutoff = weight_sum/p;  /* Cutoff for current partition */
    if (hgp->output_level >= HG_DEBUG_ALL)
      printf("GLOBAL_PART weight_sum=%f, cutoff=%f\n",weight_sum, cutoff);
  }

  if (hgp->output_level >= HG_DEBUG_ALL)
    printf("Starting new ordering at vertex %d, part=%2d\n", start_vtx, p);

  /* Initialize heap. */
  gain[start_vtx] = 1.0;           /* Make it highest value in heap. */
  Zoltan_heap_init(zz, &(h[0]), hg->nVtx);
  Zoltan_heap_init(zz, &(h[1]), 0);       /* Dummy heap, not used. */
  for(i=0; i<hg->nVtx; i++)
    Zoltan_heap_input(h, i, gain[i]);
  Zoltan_heap_make(h);

  while (bfsnumber < hg->nVtx ) {

    /* Get next vertex from heap */
    vtx = Zoltan_heap_extract_max(h);

    if (vtx<0){
      /* This should never happen. */
      ZOLTAN_PRINT_ERROR(-1, yo, "Internal error: No vertices in heap.");
      ierr = ZOLTAN_FATAL;
      goto error;
    }
    if (rank[vtx]<0){
      order[bfsnumber] = vtx;
      rank[vtx] = bfsnumber++;
    }
    else{
      sprintf(msg, "Vertex %d in heap already labeled", vtx);
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
      if (hgp->output_level >= HG_DEBUG_ALL)
        printf("GLOBAL_PART vtx=%2d, bfsnum=%2d, part[%2d]=%2d, part_sum=%f\n",
               vtx,bfsnumber-1,vtx,part[vtx],part_sum);
    }

    if (p && (pnumber+1)<p && part_sum > cutoff){
      /* Start new partition. Reset gain values. */
      pnumber++;
      /* Decide if current vertex should be moved to the next partition */
      if (part_sum-cutoff > cutoff-old_sum){
        part[vtx]++;
        part_sum = old_sum;
        if (hgp->output_level >= HG_DEBUG_ALL)
          printf("GLOBAL_PART vtx=%2d, bfsnum=%2d, part[%2d]=%2d\n",
               vtx,bfsnumber-1,vtx,part[vtx]);
      }
      weight_sum -= part_sum;
      cutoff = weight_sum/(p-pnumber);
      if (part[vtx] == pnumber){
        part_sum = hg->vwgt?hg->vwgt[vtx]:1.0;
        j = -1;
      }
      else { /* part[vtx] == pnumber-1 */
        part_sum = 0.0;
        j = Zoltan_heap_peek_max(h); /* j will be the first vertex in the next part. */
      }
      if (hgp->output_level >= HG_DEBUG_ALL)
        printf("GLOBAL_PART initializing for partition %2d, cutoff = %f\n",
               pnumber, cutoff);

      if (priority_mode>0){
        /* Reset all gain values (but one). */
        for (i=0; i<hg->nVtx; i++){
          if (i != j) gain[i] = 0.0;
          if (rank[i] <0) Zoltan_heap_change_value(h, i, gain[i]);
        }
        /* Reset counters. */
        if (vtx_count)
          for (j=0; j<hg->nEdge; j++)
            vtx_count[j] = 0;
      }
    }

    /* Update gain values for nbors. */

    if (priority_mode==0){
      /* Move from visited=0 to visited=1. */
      Zoltan_HG_move_vertex(hg, vtx, 0, 1, visited, cut, gain, h);
    }
    else {
      if (part[vtx]==pnumber){
        /* Don't update if vtx was the last in a partition. */
        for (j=hg->vindex[vtx]; j<hg->vindex[vtx+1]; j++){
          edge = hg->vedge[j];
          esize = hg->hindex[edge+1] - hg->hindex[edge];
          if (vtx_count) vtx_count[edge]++;
          for (i=hg->hindex[edge]; i<hg->hindex[edge+1]; i++){
            nbor = hg->hvertex[i];
            if (rank[nbor] <0){
               switch (priority_mode){
               case 1:
               case 2:
                 /* Absorption metric. */
                 delta = (hg->ewgt?hg->ewgt[edge]:1.0)/(esize-1);
                 break;
               case 3:
               case 4:
                 damp_factor = 0.5; /* Choose a value between 0 and 1. */
                 /* gain contribution from current edge will be
                    hg->ewgt[edge]*pow(damp_factor, esize-vtx_count[edge]-1) */
                 if (vtx_count[edge]==1)
                   delta = (hg->ewgt?hg->ewgt[edge]:1.0)*
                     pow(damp_factor, (double) (esize-2));
                 else
                   delta = (hg->ewgt?hg->ewgt[edge]:1.0)*(1.0-damp_factor)*
                     pow(damp_factor, (double) (esize-vtx_count[edge]-1));
                 break;
               }
               if (priority_mode&1)
                 gain[nbor] += delta;
               else
                 gain[nbor] += delta/edge_sum[nbor];

               Zoltan_heap_change_value(h, nbor, gain[nbor]);
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
  ZOLTAN_FREE ((void **) &gain);
  if (edge_sum)  ZOLTAN_FREE ((void **) &edge_sum);
  if (vtx_count) ZOLTAN_FREE ((void **) &vtx_count);
  if (cut[0])    ZOLTAN_FREE ((void **) &cut[0]);
  if (visited)   ZOLTAN_FREE ((void **) &visited);
  Zoltan_heap_free(&(h[0]));
  Zoltan_heap_free(&(h[1]));
  return ierr;
}

/*****************************************************************/

/* Generic greedy ordering.
 *
 * Priority function 0:
 *    gain = cut size improvement (from FM)
 * Priority function 1:  [absorption]
 *    gain(v,S) = \sum_e wgt(e) * |e \intersect S| / |e|
 * Priority function 2:
 *    gain(v,S) = \sum_e wgt(e)/edge_sum(v) * |e \intersect S| / |e|
 */

static int global_greedy (
  ZZ *zz,
  HGraph *hg,
  int p,
  Partition part,
  int pri_mode,
  HGPartParams *hgp
)
{
  int i, start, *order;
  int ierr = ZOLTAN_OK;
  char *yo = "global_greedy";

  if (!(order  = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nVtx))) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    ierr = ZOLTAN_MEMERR;
    goto error;
  }

  /* Find pseudo-peripheral start vertex */
  start = Zoltan_HG_Rand()%(hg->nVtx);
  for (i=0; i<2; i++){
    ierr = bfs_order(zz, hg, order, start, 0, 0, NULL, hgp);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
      goto error;
    start = order[hg->nVtx -1];
  }

  /* Call greedy_order. */
  ierr = greedy_order(zz, hg, order, start, pri_mode, p, part, hgp);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
    goto error;

error:
  /* Free data and return. */
  ZOLTAN_FREE ((void **) &order);
  return ierr;
}

/*****************************************************************/

/* Entry points for all the greedy methods. */

static int global_gr0 (ZZ *zz, HGraph *hg, int p, Partition part,
  HGPartParams *hgp)
{
  return global_greedy(zz, hg, p, part, 0, hgp);
}

static int global_gr1 (ZZ *zz, HGraph *hg, int p, Partition part,
  HGPartParams *hgp)
{
  return global_greedy(zz, hg, p, part, 1, hgp);
}

static int global_gr2 (ZZ *zz, HGraph *hg, int p, Partition part,
  HGPartParams *hgp)
{
  return global_greedy(zz, hg, p, part, 2, hgp);
}

static int global_gr3 (ZZ *zz, HGraph *hg, int p, Partition part,
  HGPartParams *hgp)
{
  return global_greedy(zz, hg, p, part, 3, hgp);
}

static int global_gr4 (ZZ *zz, HGraph *hg, int p, Partition part,
  HGPartParams *hgp)
{
  return global_greedy(zz, hg, p, part, 4, hgp);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
