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

#include "zz_sort.h"
#include "phg.h"
#include "zz_heap.h"
    
/* If adding a new coarse partitioning fn, add prototype here 
 * AND add entry to CoarsePartitionFns array 
 * AND increment NUM_COARSEPARTITION_FN.
 */
#define NUM_COARSEPARTITION_FNS 3

static ZOLTAN_PHG_COARSEPARTITION_FN coarse_part_greedy;
static ZOLTAN_PHG_COARSEPARTITION_FN coarse_part_random;
static ZOLTAN_PHG_COARSEPARTITION_FN coarse_part_linear;

static ZOLTAN_PHG_COARSEPARTITION_FN* CoarsePartitionFns[] = 
                                      {&coarse_part_greedy,
                                       &coarse_part_random,
                                       &coarse_part_linear,
                                      };

static int local_coarse_partitioner(ZZ *, HGraph *, int, float *, Partition,
  PHGPartParams *, ZOLTAN_PHG_COARSEPARTITION_FN *);

static int pick_best(ZZ*, PHGPartParams*, PHGComm*, HGraph*, int, int, int*, float*);

/****************************************************************************/

ZOLTAN_PHG_COARSEPARTITION_FN *Zoltan_PHG_Set_CoarsePartition_Fn(
  PHGPartParams *hgp,
  int *ierr
)
{
char *str, *str2;

  *ierr = ZOLTAN_OK;

  /* local partitioning option is not effective and
     should only be used for debugging and testing */
  str2 = hgp->coarsepartition_str;
  if (!strncasecmp(str2, "l-", 2)) {
    str = str2+2;
    hgp->LocalCoarsePartition = 1;
  }  
  else {
    str = str2;
    hgp->LocalCoarsePartition = 0;
  }

  if      (!strcasecmp(str, "auto"))   return NULL; /* try all methods */
  else if (!strcasecmp(str, "no"))     return NULL; 
  else if (!strcasecmp(str, "none"))   return NULL; 
  else if (!strcasecmp(str, "greedy")) return coarse_part_greedy;
  else if (!strcasecmp(str, "random")) return coarse_part_random;
  else if (!strcasecmp(str, "linear")) return coarse_part_linear;
  else {                              
    *ierr = ZOLTAN_FATAL; 
    return NULL;
  }
}

/****************************************************************************/

#define NUM_PART_KEEP 1            /* No. of partition vectors to keep; 
                                      must be at least 1! Currently only the
                                      best partition vector is used. */

int Zoltan_PHG_CoarsePartition(
  ZZ *zz, 
  HGraph *phg,         /* Input:  coarse hypergraph -- distributed! */
  int numPart,         /* Input:  number of partitions to generate. */
  float *part_sizes,   /* Input:  array of size numPart listing target sizes
                                  (% of work) for the partitions */
  Partition part,      /* Input:  array of initial partition assignments.
                          Output: array of computed partition assignments.   */
  PHGPartParams *hgp   /* Input:  parameters to use.  */
)
{
/* 
 * Zoltan_PHG_CoarsePartition computes a partitioning of a hypergraph.
 * Typically, this routine is called at the bottom level in a
 * multilevel scheme (V-cycle).
 * It gathers the distributed hypergraph to each processor and computes
 * a decomposition of the serial hypergraph.  
 * It computes a different partition on each processor
 * using different random numbers (and possibly also
 * different algorithms) and selects the best.
 */
char *yo = "Zoltan_PHG_CoarsePartition";
int ierr = ZOLTAN_OK;
int i, si, j;
static PHGComm scomm;          /* Serial communicator info */
static int first_time = 1;
HGraph *shg = NULL;            /* Serial hypergraph gathered from phg */
int *spart = NULL;             /* Partition vectors for shg. */
int *new_part = NULL;          /* Ptr to new partition vector. */
float *bestvals = NULL;        /* Best cut values found so far */
int worst, new_cand;
float bal, cut, worst_cut;
int fine_timing = (hgp->use_timers > 2);
static int timer_cpart=-1, timer_gather=-1, timer_refine=-1; 

/* Number of iterations to try coarse partitioning on each proc. */
/* 10 when p=1, and 1 when p is large. */
const int num_coarse_iter = 1 + 9/zz->Num_Proc; 

  ZOLTAN_TRACE_ENTER(zz, yo);

  if (fine_timing) {
    if (timer_gather < 0)
      timer_gather = Zoltan_Timer_Init(zz->ZTime, 1, "CP Gather");
    if (timer_refine < 0)
      timer_refine = Zoltan_Timer_Init(zz->ZTime, 0, "CP Refine");
    if (timer_cpart < 0)
      timer_cpart = Zoltan_Timer_Init(zz->ZTime, 0, "CP Part");

    ZOLTAN_TIMER_START(zz->ZTime, timer_cpart, phg->comm->Communicator);
  }

  /* Force LocalCoarsePartition if large global graph */
#define LARGE_GRAPH_VTX   32000
#define LARGE_GRAPH_PINS 256000
  if (phg->dist_x[phg->comm->nProc_x] > LARGE_GRAPH_VTX){
    /* TODO: || (global_nPins > LARGE_GRAPH_PINS) */
    hgp->LocalCoarsePartition = 1;
  }

  /* take care of all special cases first */

  if (!strcasecmp(hgp->coarsepartition_str, "no")) {
    /* Do no coarse partitioning. */
    /* Do a sanity test and  mapping to parts [0,...,numPart-1] */
    int first = 1;
    for (i = 0; i < phg->nVtx; i++)
      if (part[i] >= numPart) {
        if (first) {
          ZOLTAN_PRINT_WARN(zz->Proc, yo, "Initial part number > numParts.");
          first = 0;
          ierr = ZOLTAN_WARN;
        }
        part[i] = part[i] % numPart;
      }
  }
  else if (numPart == 1) {            
    /* everything goes in the one partition */
    for (i =  0; i < phg->nVtx; i++)
      part[i] = 0;
  }
  else if (numPart >= phg->dist_x[phg->comm->nProc_x]) { 
    /* more partitions than vertices, trivial answer */
    for (i = 0; i < phg->nVtx; i++)
      part[i] = phg->dist_x[phg->comm->myProc_x]+i;
  }
  else if (hgp->LocalCoarsePartition) {
    /* Apply local partitioner to each column */
    ierr = local_coarse_partitioner(zz, phg, numPart, part_sizes, part, hgp,
                                    hgp->CoarsePartition);
  }
  else {
    /* Normal case:
     * Gather distributed HG to each processor;
     * compute different partitioning on each processor;
     * select the "best" result.
     */
    ZOLTAN_PHG_COARSEPARTITION_FN *CoarsePartition;

    /* Select different coarse partitioners for processors here. */

    CoarsePartition = hgp->CoarsePartition;
    if (CoarsePartition == NULL) { /* auto */
      /* Select a coarse partitioner from the array of coarse partitioners */
      CoarsePartition = CoarsePartitionFns[phg->comm->myProc % 
                                           NUM_COARSEPARTITION_FNS];
    }


    if (phg->comm->nProc == 1) {
      /* Serial and parallel hgraph are the same. */
      shg = phg;
    }
    else {
      /* Set up a serial communication struct for gathered HG */

      if (first_time) {
        scomm.nProc_x = scomm.nProc_y = 1;
        scomm.myProc_x = scomm.myProc_y = 0;
        scomm.Communicator = MPI_COMM_SELF;
        scomm.row_comm = MPI_COMM_SELF;
        scomm.col_comm = MPI_COMM_SELF;
        scomm.myProc = 0;
        scomm.nProc = 1;
        first_time = 0;
      }
      scomm.RNGState = Zoltan_Rand(NULL);
      scomm.RNGState_row = Zoltan_Rand(NULL);
      scomm.RNGState_col = Zoltan_Rand(NULL);
      scomm.zz = zz;

      /* 
       * Gather parallel hypergraph phg to each processor, creating
       * serial hypergraph shg.
       */
      if (fine_timing) {
        ZOLTAN_TIMER_STOP(zz->ZTime, timer_cpart, phg->comm->Communicator);
        ZOLTAN_TIMER_START(zz->ZTime, timer_gather, phg->comm->Communicator);
      }

      ierr = Zoltan_PHG_Gather_To_All_Procs(zz, phg, &scomm, &shg);
      if (ierr < 0) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from gather.");
        goto End;
      }

      if (fine_timing) {
        ZOLTAN_TIMER_STOP(zz->ZTime, timer_gather, phg->comm->Communicator);
        ZOLTAN_TIMER_START(zz->ZTime, timer_cpart, phg->comm->Communicator);
      }

    }

    /* 
     * Allocate partition array spart for the serial hypergraph shg
     * and partition shg.
     */
    spart = (int *) ZOLTAN_CALLOC(shg->nVtx * (NUM_PART_KEEP+1),
                                    sizeof(int));
    bestvals = (float *) ZOLTAN_MALLOC((NUM_PART_KEEP+1)*sizeof(int)); 
    if ((!spart) || (!bestvals)) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Out of memory.");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
    
    /* Compute several coarse partitionings. */
    /* Keep the NUM_PART_KEEP best ones around. */
    /* Currently, only the best one is used. */

    /* Set RNG so different procs compute different parts. */
    Zoltan_Srand(Zoltan_Rand(NULL) + zz->Proc, NULL);

    new_cand = 0;
    new_part = spart;

    for (i=0; i< num_coarse_iter; i++){
      int savefmlooplimit=hgp->fm_loop_limit;
        
      /* Overwrite worst partition with new candidate. */
      ierr = CoarsePartition(zz, shg, numPart, part_sizes, 
               new_part, hgp);
      if (ierr < 0) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                         "Error returned from CoarsePartition.");
        goto End;
      }

      /* time refinement step in coarse partitioner */
      if (fine_timing) {
        ZOLTAN_TIMER_STOP(zz->ZTime, timer_cpart, phg->comm->Communicator);
        ZOLTAN_TIMER_START(zz->ZTime, timer_refine, phg->comm->Communicator);
      }

      /* UVCUVC: Refine new candidate: only one pass is enough. */
      hgp->fm_loop_limit = 1;
      Zoltan_PHG_Refinement(zz, shg, numPart, part_sizes, new_part, hgp);
      hgp->fm_loop_limit = savefmlooplimit;
      
      /* stop refinement timer */
      if (fine_timing) {
        ZOLTAN_TIMER_STOP(zz->ZTime, timer_refine, phg->comm->Communicator);
        ZOLTAN_TIMER_START(zz->ZTime, timer_cpart, phg->comm->Communicator);
      }

      /* Decide if candidate is in the top tier or not. */
      /* Our objective is a combination of cuts and balance */

      bal = Zoltan_PHG_Compute_Balance(zz, shg, part_sizes, numPart, new_part); 
      cut = Zoltan_PHG_Compute_ConCut(shg->comm, shg, new_part, numPart, &ierr);
      
      /* Use ratio-cut as our objective. There are many other options! */
      bestvals[new_cand] = cut/(MAX(2.-bal, 0.0001)); /* avoid divide-by-0 */

      if (ierr < 0) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                         "Error returned from Zoltan_PHG_Compute_ConCut.");
        goto End;
      }
      if (i<NUM_PART_KEEP)
        new_cand = i+1;
      else {
        /* find worst partition vector, to overwrite it */
        /* future optimization: keep bestvals sorted */
        worst = 0;
        worst_cut = bestvals[0];
        for (j=1; j<NUM_PART_KEEP+1; j++){
          if (worst_cut < bestvals[j]){
            worst_cut = bestvals[j];
            worst = j;
          }
        }
        new_cand = worst;
      }
      new_part = spart+new_cand*(shg->nVtx);
    }
    /* Copy last partition vector such that all the best ones
       are contiguous starting at spart.                     */
    for (i=0; i<shg->nVtx; i++){
      new_part[i] = spart[NUM_PART_KEEP*(shg->nVtx)+i];
    }
    /* Also update bestvals */
    bestvals[new_cand] = bestvals[NUM_PART_KEEP];

    /* Evaluate and select the best. */
    /* For now, only pick the best one, in the future we pick the k best. */

    ierr = pick_best(zz, hgp, phg->comm, shg, numPart, 
              MIN(NUM_PART_KEEP, num_coarse_iter), spart,
              bestvals);
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                        "Error returned from pick_best.");
      goto End;
    }
  
    if (phg->comm->nProc > 1) {
      /* Map gathered partition back to 2D distribution */
      for (i = 0; i < phg->nVtx; i++) {
        /* KDDKDD  Assume vertices in serial HG are ordered by GNO of phg */
        si = VTX_LNO_TO_GNO(phg, i);
        part[i] = spart[si];
      }

      Zoltan_HG_HGraph_Free(shg);
      ZOLTAN_FREE(&shg);
    } 
    else { /* single processor */
      for (i = 0; i < phg->nVtx; i++)
        part[i] = spart[i];
    }
    ZOLTAN_FREE(&spart);
    ZOLTAN_FREE(&bestvals);
  }
  
End:
  if (fine_timing) 
    ZOLTAN_TIMER_STOP(zz->ZTime, timer_cpart, phg->comm->Communicator);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

/****************************************************************************/

static int local_coarse_partitioner(
  ZZ *zz,
  HGraph *hg,
  int numPart,
  float *part_sizes,
  Partition part,
  PHGPartParams *hgp,
  ZOLTAN_PHG_COARSEPARTITION_FN *CoarsePartition
)
{
/* 
 * Function that allows any of the coarse partitioning strategies to be applied
 * locally, without a gather operation.
 * Each column group does independent partitioning; each tries to balance
 * local vertex weights; our hope is that this will give somewhat balanced 
 * result.
 */
PHGComm *hgc = hg->comm;
int err=0;
int rootnpins, rootrank;

  /* The column processor with the most pins will be our root.  */
  Zoltan_PHG_Find_Root(hg->nPins, hgc->myProc_y, hgc->col_comm, 
                       &rootnpins, &rootrank);

  if (hgc->myProc_y==rootrank)   /* only root of each column does this */
    err = CoarsePartition(zz, hg, numPart, part_sizes, part, hgp);  

  MPI_Bcast(&err, 1, MPI_INT, rootrank, hgc->col_comm);
  if (!err)
    MPI_Bcast(part, hg->nVtx, MPI_INT, rootrank, hgc->col_comm);
    
  return err;
}


/****************************************************************************/

/* Sequence partitioning on the vertices of a hypergraph
   in a given order. The goal is to assign approximately
   even weights to each partition, or alternatively, proportional
   to the target weights (if such are given), subject to
   a specified linear order (sequence). 

   Multi-weights are not yet supported; only the
   first weight is used in computing the partition.

   This is a quick but effective heuristic. We could alternatively use
   a more expensive but optimal algorithm, see e.g. 
   "Fast Optimal Load Balancing Algorithms for 1D Partitioning"
   by Ali Pinar and C. Aykanat, but for our purpose it is 
   probably not worth the effort.

   Adapted for fixed vertices.

*/

static int seq_part (
  ZZ *zz, 
  HGraph *hg,         /* the hypergraph, containing vertex weights */
  int *order,         /* the ordering of the vertices */
  int p,              /* desired number of partitions */
  float *part_sizes,  /* target partition sizes */
  Partition part,     /* Output: partition numbers */
  PHGPartParams *hgp  /* misc hypergraph parameters */
)
{
  int i, j, pnumber;
  int vwgtdim = hg->VtxWeightDim;
  double weight_sum = 0.0, part_sum, old_sum, cutoff;
  double psize_sum = 0.0;
  double *fixed_wgts = NULL;
  char *yo = "seq_part";

  ZOLTAN_TRACE_ENTER(zz, yo);

  if (part_sizes==NULL){
    /* part_sizes should always exist, even with uniform partitions */
    return ZOLTAN_FATAL;
  }

  if (p<1){
    /* should never happen */
    return ZOLTAN_FATAL;
  }

#ifdef DEBUG_FIXED
  hg->fixed = (int *) ZOLTAN_MALLOC(hg->nVtx*sizeof(int));
  for (i=0; i<hg->nVtx; i++)
    hg->fixed[i] = -1;
  /* TEST: Fix the first two vertices */
  hg->fixed[0] = 0;
  hg->fixed[1] = 1;
#endif

  if (hg->fixed) {
    fixed_wgts = (double *) ZOLTAN_CALLOC(p, sizeof(double));
  }

  /* Sum up all the vertex weights. */
  for (i=0; i<hg->nVtx; i++){
    weight_sum += hg->vwgt[i*vwgtdim];
    if (hg->fixed)
      if (hg->fixed[i] >= 0){
        /* Set partition number for fixed vtx. */
        part[i] = hg->fixed[i];
        /* Add up weights of fixed vertices for each partition */
        fixed_wgts[hg->fixed[i]] += hg->vwgt[i*vwgtdim];
      }
  }

  /* Sum up all the target partition weights. */
  /* Only use first vweight for now. */
  for (i=0; i<p; i++)
    psize_sum += part_sizes[i*vwgtdim];

  pnumber = 0; /* Assign next vertex to partition no. pnumber */
  part_sum = (fixed_wgts ? fixed_wgts[0] : 0.); /* Weight of fixed vertices */

  /* Set cutoff for current partition */
  cutoff = weight_sum*part_sizes[0]/psize_sum - 
           (fixed_wgts ? fixed_wgts[0] : 0.);  

  /* Loop through all vertices in specified order, and assign
     partition numbers.  */                                        
  for (i=0; i<hg->nVtx; i++) {
    /* for non-fixed vertices */
    if ((!hg->fixed) || (hg->fixed[i] == -1)){
      /* If order==NULL, then use linear order. */
      j = order ? order[i] : i;
      part[j] = pnumber;
      old_sum = part_sum;
      part_sum += hg->vwgt[j*vwgtdim];
      /* Check if we passed the cutoff and should start a new partition */
      if ((pnumber+1) < p && part_sum > cutoff) {
        pnumber++; /* Increase current part number */
        /* Decide if current vertex should be moved to the next partition */
        if ((part_sum-cutoff) > (cutoff-old_sum)) { 
          part[j]++;
          part_sum = old_sum;
        }
        weight_sum -= part_sum;
        /* Initialize part_sum for next partition no. */
        part_sum = (fixed_wgts ? fixed_wgts[pnumber] : 0.);
        if (part[j] == pnumber)
          part_sum += hg->vwgt[j*vwgtdim];
        /* Update cutoff. */
        psize_sum -= part_sizes[pnumber-1];
        cutoff = weight_sum*part_sizes[pnumber]/psize_sum
                 - (fixed_wgts ? fixed_wgts[pnumber] : 0.);
      }
    }
    if (hgp->output_level >= PHG_DEBUG_ALL)
      printf("COARSE_PART i=%2d, part[%2d] = %2d, part_sum=%f, cutoff=%f\n", 
       i, j, part[j], part_sum, cutoff);
  }

  if (hg->fixed) {
    ZOLTAN_FREE(&fixed_wgts);
#ifdef DEBUG_FIXED
    ZOLTAN_FREE(&hg->fixed);
#endif
  }

  ZOLTAN_TRACE_EXIT(zz, yo);
  return ZOLTAN_OK;
}

/****************************************************************************/
/* Linear partitioning. Sequence partitioning with vertices in linear order. */

static int coarse_part_linear (
  ZZ *zz, 
  HGraph *hg, 
  int p, 
  float *part_sizes,
  Partition part, 
  PHGPartParams *hgp
)
{
    int i, offset, err=0, *order=NULL;
    static char *yo = "coarse_part_linear";

    if (!(order  = (int*) ZOLTAN_MALLOC (hg->nVtx*sizeof(int)))) {
        ZOLTAN_FREE ((void**) &order);
        ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
        return ZOLTAN_MEMERR;
    }

    /* Make sure all procs do different variations of linear partitioning.
       We achieve diversity by picking different starting vertices.
       The vertex ordering is still linear, but with cyclic "wrap-around".
     */
    if (zz->Proc == 0)
      offset = 0;  /* Special case for proc 0 is not really necessary */
    else
      offset = Zoltan_Rand(NULL) % (hg->nVtx);

    for (i=0; i<hg->nVtx; i++) {
        order[i] = offset + i;
        if (order[i] >= hg->nVtx) 
          order[i] -= hg->nVtx;
    }

    /* Call sequence partitioning with order array. */
    err = seq_part (zz, hg, order, p, part_sizes, part, hgp);

    ZOLTAN_FREE ((void**) &order);
    return err;
}



/****************************************************************************/
/* Random partitioning. Sequence partitioning with vertices in random order. */
static int coarse_part_random (
  ZZ *zz,
  HGraph *hg,
  int p,
  float *part_sizes,
  Partition part,
  PHGPartParams *hgp
)
{
    int i, err=0, *order=NULL;
    char *yo = "coarse_part_random";

    if (!(order  = (int*) ZOLTAN_MALLOC (hg->nVtx*sizeof(int)))) {
        ZOLTAN_FREE ((void**) &order);
        ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
        return ZOLTAN_MEMERR;
    }
    for (i=0; i<hg->nVtx; i++) {
        order[i] = i;
    }

    /* Randomly permute order array */
    Zoltan_Rand_Perm_Int (order, hg->nVtx, NULL);
        
    /* Call sequence partitioning with random order array. */
    err = seq_part (zz, hg, order, p, part_sizes, part, hgp);

    ZOLTAN_FREE ((void**) &order);
    return err;
}


/*********************************************************************/
/* Greedy ordering/partitioning based on a priority function
   for selecting vertices. A heap is used as a priority queue. */
 
int Zoltan_PHG_move_vertex (HGraph *hg, int vertex, int sour, int dest,
                            int *part, int **cut, double *gain, HEAP *heap)
{
    int i, j, edge, v;

    gain[vertex] = 0.0;
    part[vertex] = dest;

    for (i = hg->vindex[vertex]; i < hg->vindex[vertex+1]; i++) {
        edge = hg->vedge[i];
        if (cut[sour][edge] == 1) {
            for (j = hg->hindex[edge]; j < hg->hindex[edge+1]; j++) {
                v = hg->hvertex[j];
                gain[v] -= (hg->ewgt ? hg->ewgt[edge] : 1.0);
                if (heap)
                    Zoltan_Heap_Change_Value(&heap[part[v]], v, gain[v]);
            }
        }
        else if (cut[sour][edge] == 2) {
            for (j = hg->hindex[edge]; j < hg->hindex[edge+1]; j++) {
                v = hg->hvertex[j];
                if (part[v] == sour) {
                    gain[v] += (hg->ewgt ? hg->ewgt[edge] : 1.0);
                    if (heap)
                        Zoltan_Heap_Change_Value(&heap[part[v]], v, gain[v]);
                    break;
                }
            }
        }

        if (cut[dest][edge] == 0) {
            for (j = hg->hindex[edge]; j < hg->hindex[edge+1]; j++) {
                v = hg->hvertex[j];
                gain[v] += (hg->ewgt ? hg->ewgt[edge] : 1.0);
                if (heap)
                    Zoltan_Heap_Change_Value(&heap[part[v]], v, gain[v]);
            }
        }
        else if (cut[dest][edge] == 1) {
            for (j = hg->hindex[edge]; j < hg->hindex[edge+1]; j++) {
                v = hg->hvertex[j];
                if (v != vertex && part[v] == dest) {
                    gain[v] -= (hg->ewgt ? hg->ewgt[edge] : 1.0);
                    if (heap)
                        Zoltan_Heap_Change_Value(&heap[part[v]], v, gain[v]);
                    break;
                }
            }
        }
        cut[sour][edge]--;
        cut[dest][edge]++;
    }
    return ZOLTAN_OK;
}

static int greedy_order (
  ZZ *zz,
  HGraph *hg,		/* Hypergraph. */
  int *order,		/* Order array. On exit, order[i] is the i'th vertex. */
  int start_vtx,	/* Start the ordering from this vertex. */
  int priority_mode,	/* Priority mode for selecting vertices */
  int p,		/* Optional (input):  Number of partitions. */
  float *part_sizes,    /* Array of length p containing the percentages of
                           work to be assigned to each partition. */
  Partition part,	/* Optional (output): Partition array. */
  PHGPartParams *hgp     /* Partitioning parameters. */
)
{
  int i, j, vtx, edge, bfsnumber, pnumber, nbor, *rank;
  int esize, *vtx_count=NULL, *visited=NULL, *cut[2];
  int vwgtdim = hg->VtxWeightDim;
  int err=ZOLTAN_OK;
  double weight_sum= 0.0, part_sum= 0.0, old_sum;
  double delta=0.0, cutoff=0.0;
  double *gain = NULL, *edge_sum = NULL;
  double damp_factor, psize_sum= 0.0;
  char msg[128];
  HEAP h[2];
  static char *yo = "greedy_order";

  bfsnumber = 0;  /* Assign next vertex this bfs number */
  pnumber = 0;    /* Assign next vertex this partition number */

  /* Allocate arrays. */
  if (!(rank  = (int*)    ZOLTAN_CALLOC (hg->nVtx, sizeof (int))) ||
      !(gain  = (double*) ZOLTAN_CALLOC (hg->nVtx, sizeof (double))) ) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    err =  ZOLTAN_MEMERR;
    goto End;
  }
  for (i=0; i<hg->nVtx; i++)
    rank[i] = -1;       /* -1 means this vtx has not yet been numbered */
 
  if (priority_mode && (!(priority_mode&1))) {   /* 2,4,6,... */
    if (!(edge_sum = (double*) ZOLTAN_CALLOC (hg->nVtx, sizeof (double)))){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      err = ZOLTAN_MEMERR;
      goto End;
    }
    /* Sum up edge weights incident to each vertex. */
    for (edge=0; edge<hg->nEdge; edge++) {
      for (i=hg->hindex[edge]; i<hg->hindex[edge+1]; i++){
        edge_sum[hg->hvertex[i]] += (hg->ewgt ? hg->ewgt[edge] : 1.0);
      }
    }
  }
  if (priority_mode == 0) {
    cut[0]  = (int*) ZOLTAN_CALLOC (2*hg->nEdge, sizeof (int));
    visited = (int*) ZOLTAN_CALLOC (hg->nVtx,    sizeof (int));
    if ((hg->nEdge > 0 && cut[0] == NULL) || visited == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      err = ZOLTAN_MEMERR;
      goto End;
    }
    cut[1] = &(cut[0][hg->nEdge]);
    /* Initialize cut values. */
    for (i=0; i<hg->nEdge; i++)
      for (j=hg->hindex[i]; j<hg->hindex[i+1]; j++)
        (cut[visited[hg->hvertex[j]]][i])++;
    /* Initialize gain values. */
    for (i=0; i<hg->nVtx; i++){
      for (j=hg->vindex[i]; j<hg->vindex[i+1]; j++) {
        edge = hg->vedge[j];
        gain[i] -= (hg->ewgt ? (hg->ewgt[edge]) : 1.0);
      }
    }
  }
  else
    cut[0] = cut[1] = NULL;

  if (priority_mode >= 3) {
    if (!(vtx_count = (int*) ZOLTAN_CALLOC (hg->nEdge, sizeof (int)))){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      err = ZOLTAN_MEMERR;
      goto End;
    }
  }

  if (p) {
    /* If partitioning is chosen, sum up all the weights. */
    for (i=0; i<hg->nVtx; i++)
      weight_sum += hg->vwgt[i*vwgtdim];

    /* Sum up all the target partition weights. */
    /* Only use first vweight for now. */
    for (i=0; i<p; i++)
      psize_sum += part_sizes[i*vwgtdim];

    /* Set cutoff for current partition */
    cutoff = weight_sum*part_sizes[0]/psize_sum;
  }

  if (hgp->output_level >= PHG_DEBUG_ALL)
    printf("Starting new ordering at vertex %d, part=%2d\n", start_vtx, p);

  /* Initialize heap. */
  gain[start_vtx] = 1.0;           /* Make it highest value in heap. */
  Zoltan_Heap_Init(zz, &h[0], hg->nVtx);
  Zoltan_Heap_Init(zz, &h[1], 0);       /* Dummy heap, not used. */
  for(i=0; i<hg->nVtx; i++)
    Zoltan_Heap_Input(h, i, gain[i]);
  Zoltan_Heap_Make(h);

  while (bfsnumber < hg->nVtx ) {

    /* Get next vertex from heap */
    vtx = Zoltan_Heap_Extract_Max(h);

    if (vtx < 0) {
      /* This should never happen. */
      ZOLTAN_PRINT_ERROR(-1, yo, "Internal error: No vertices in heap.");
      err = ZOLTAN_FATAL;
      goto End;
    }
    if (rank[vtx] < 0){
      order[bfsnumber] = vtx;
      rank[vtx] = bfsnumber++;
    }
    else{
      sprintf(msg, "Vertex %d in heap already labeled", vtx);
      ZOLTAN_PRINT_ERROR(-1, yo, msg);
      sprintf(msg, "bfsnumber=%d, rank[vtx] = %d", bfsnumber, rank[vtx]);
      ZOLTAN_PRINT_ERROR(-1, yo, msg);
      err = ZOLTAN_FATAL;
      goto End;
    }
    if (p) {
      old_sum = part_sum;
      part_sum += hg->vwgt[vtx*vwgtdim];
      part[vtx] = pnumber;
      if (hgp->output_level >= PHG_DEBUG_ALL)
        printf("COARSE_PART vtx=%2d, bfsnum=%2d, part[%2d]=%2d, part_sum=%f\n",
               vtx,bfsnumber-1,vtx,part[vtx],part_sum);
    }

    if (p && (pnumber+1)<p && part_sum > cutoff) {
      /* Start new partition. Reset gain values. */
      pnumber++;
      /* Decide if current vertex should be moved to the next partition */
      if (part_sum-cutoff > cutoff-old_sum) {
        part[vtx]++;
        part_sum = old_sum;
        if (hgp->output_level >= PHG_DEBUG_ALL)
          printf("COARSE_PART vtx=%2d, bfsnum=%2d, part[%2d]=%2d\n",
           vtx, bfsnumber-1, vtx, part[vtx]);
      }
      weight_sum -= part_sum;
      if (part[vtx] == pnumber){
        part_sum = hg->vwgt[vtx*vwgtdim];
        j = -1;
      }
      else { /* part[vtx] == pnumber-1 */
        part_sum = 0.0;
        j = Zoltan_Heap_Peek_Max(h); /* j will be the first vertex in the next part. */
      }
      /* Update cutoff. */
      psize_sum -= part_sizes[pnumber-1];
      cutoff = weight_sum*part_sizes[pnumber]/psize_sum;

      if (hgp->output_level >= PHG_DEBUG_ALL)
        printf("COARSE_PART vtx=%2d, part[%2d] = %2d, part_sum=%f, cutoff=%f\n",
          vtx, vtx, part[vtx], part_sum, cutoff);

      if (priority_mode > 0) {
        /* Reset all gain values (but one). */
        for (i=0; i<hg->nVtx; i++){
          if (i != j) gain[i] = 0.0;
          if (rank[i] < 0) Zoltan_Heap_Change_Value(h, i, gain[i]);
        }
        /* Reset counters. */
        if (vtx_count)
          for (j=0; j<hg->nEdge; j++)
            vtx_count[j] = 0;
      }
    }

    /* Update gain values for nbors. */
    if (priority_mode == 0) {
      /* Move from visited=0 to visited=1. */
      Zoltan_PHG_move_vertex(hg, vtx, 0, 1, visited, cut, gain, h);
    }
    else {
      if (part[vtx] == pnumber) {
        /* Don't update if vtx was the last in a partition. */
        for (j=hg->vindex[vtx]; j<hg->vindex[vtx+1]; j++) {
          edge = hg->vedge[j];
          esize = hg->hindex[edge+1] - hg->hindex[edge];
          if (vtx_count) vtx_count[edge]++;
          for (i=hg->hindex[edge]; i<hg->hindex[edge+1]; i++) {
            nbor = hg->hvertex[i];
            if (rank[nbor] < 0) {
               switch (priority_mode) {
               case 1:
               case 2:
                 /* Absorption metric. */
                 delta = (hg->ewgt ? hg->ewgt[edge] : 1.0)/(esize-1);
                 break;
               case 3:
               case 4:
                 damp_factor = 0.5; /* Choose a value between 0 and 1. */
                 /* gain contribution from current edge will be
                    hg->ewgt[edge]*pow(damp_factor, esize-vtx_count[edge]-1) */
                 if (vtx_count[edge] == 1)
                   delta = (hg->ewgt ? hg->ewgt[edge] : 1.0)
                    * pow(damp_factor, (double) (esize-2));
                 else
                   delta = (hg->ewgt ? hg->ewgt[edge] : 1.0) * (1.0-damp_factor)
                    * pow(damp_factor, (double) (esize-vtx_count[edge]-1));
                 break;
               }
               if (priority_mode & 1)
                 gain[nbor] += delta;
               else
                 gain[nbor] += delta/edge_sum[nbor];

               Zoltan_Heap_Change_Value(h, nbor, gain[nbor]);
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

End:
  ZOLTAN_FREE ((void**) &rank);
  ZOLTAN_FREE ((void**) &gain);
  if (edge_sum)  ZOLTAN_FREE ((void**) &edge_sum);
  if (vtx_count) ZOLTAN_FREE ((void**) &vtx_count);
  if (cut[0])    ZOLTAN_FREE ((void**) &cut[0]);
  if (visited)   ZOLTAN_FREE ((void**) &visited);
  Zoltan_Heap_Free (&h[0]);
  Zoltan_Heap_Free( &h[1]);
  return err;
}



/*****************************************************************/
/* Generic greedy ordering. 
 * Default priority function (0):
 *    gain = cut size improvement (like FM gain)
 *
 * The options below are no longer supported.
 * Priority function 1:  [absorption]
 *    gain(v,S) = \sum_e wgt(e) * |e \intersect S| / |e|
 * Priority function 2:
 *    gain(v,S) = \sum_e wgt(e)/edge_sum(v) * |e \intersect S| / |e|
 */
static int coarse_part_greedy (
  ZZ *zz,
  HGraph *hg,
  int p,
  float *part_sizes,
  Partition part,
  PHGPartParams *hgp
)
{
  int start, *order;
  int err = ZOLTAN_OK;
  const int pri_mode = 0;  /* set to 0; other options not supported any more  */
  char *yo = "coarse_part_greedy";

  if (!(order  = (int*) ZOLTAN_MALLOC (sizeof(int) * hg->nVtx))) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    err = ZOLTAN_MEMERR;
    goto End;
  }

  /* Start at random vertex */
  start = Zoltan_Rand(NULL) % (hg->nVtx);

  /* Call greedy_order. */
  err = greedy_order(zz, hg, order, start, pri_mode, p, part_sizes, part, hgp);
  if (err != ZOLTAN_OK && err != ZOLTAN_WARN)
    goto End;

End:
  /* Free data and return. */
  ZOLTAN_FREE ((void**) &order);
  return err;
}


/*****************************************************************************/
static int pick_best(
  ZZ *zz, 
  PHGPartParams *hgp,  
  PHGComm *phg_comm,
  HGraph *shg, 
  int numPart,
  int numLocalCandidates,
  int *spart,     /* On input: All candidate partition vectors concatenated.
                     On output: Best partition vectors concatenated. */
  float *cutvals  /* Reuse cut values if available */
)
{
/* Routine to select the best of multiple serial partitions and broadcast
 * it to all processors
 */
struct {
  float val;
  int rank;
} local[2], global[2];

static char *yo = "pick_best";
int i, mybest;
float cut, bal;
int err = ZOLTAN_OK;

  /* find best local partition */
  /* if cut values are given on input use these (may be ratio-cut),
     otherwise, only look at cuts not balance. (EB should never happen) */

  mybest = 0;
  /* local[0].val = Zoltan_PHG_Compute_Balance(zz, shg, part_sizes, numPart, spart); */
  local[0].val = bal = 0.0; /* balance */
  if (cutvals)
    local[1].val = cutvals[0];
  else
    local[1].val = Zoltan_PHG_Compute_ConCut(shg->comm, shg, spart, numPart,
                                             &err);
    if (err < 0) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                         "Error returned from Zoltan_PHG_Compute_ConCut");
      goto End;
    }
  local[0].rank = local[1].rank = phg_comm->myProc;

  /* What do we say is "best"?   For now, say lowest (ratio) cut. */
  for (i=1; i<numLocalCandidates; i++){
    /*bal = Zoltan_PHG_Compute_Balance(zz, shg, part_sizes, numPart, spart+i*(shg->nVtx));*/

    if (cutvals)
      cut = cutvals[i];
    else
      cut = Zoltan_PHG_Compute_ConCut(shg->comm, shg, 
             spart+i*(shg->nVtx), numPart, &err);
    if (err < 0) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                         "Error returned from Zoltan_PHG_Compute_ConCut");
      goto End;
    }

    if (cut < local[1].val){
      mybest = i;
      local[0].val = bal;  /* currently not used */
      local[1].val = cut;
    }
  }
  
  /* copy best partition to beginning of spart */
  for (i=0; i<shg->nVtx; i++)
    spart[i] = spart[mybest*(shg->nVtx)+i];

  /* Pick lowest ratio cut as best. */
  MPI_Allreduce(local, global, 2, MPI_FLOAT_INT, MPI_MINLOC, 
                phg_comm->Communicator);

  if (hgp->output_level)
    uprintf(phg_comm,
            "Local Ratio Cut= %.2lf   Global Ratio Cut= %.2lf\n", 
             local[1].val, global[1].val);

  MPI_Bcast(spart, shg->nVtx, MPI_INT, global[1].rank, 
            phg_comm->Communicator);

End:
  return err;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

