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

#include <math.h>
#include "phg.h"
#include "phg_tree.h"
#include "params_const.h"
#include "all_allo_const.h"
#include "zz_const.h"


/*
#define CHECK_LEFTALONE_VERTICES
*/
    
/*
 *  Main routine for Zoltan interface to hypergraph partitioning. 
 *  Also routines that build input data structures, set parameters, etc.
 */


/******************************************************************************/
/*  Parameters structure for parallel HG method.  */
static PARAM_VARS PHG_params[] = {
  /* Add parameters here. */
  {"HYPERGRAPH_PACKAGE",              NULL,  "STRING", 0},
    /* Software package: PHG (Zoltan) or Patoh */
  {"PHG_MULTILEVEL",                  NULL,  "INT", 0},
    /* Indicate whether or not to use multilevel method (1/0) */
  {"PHG_CUT_OBJECTIVE",               NULL,  "STRING", 0},
    /* connectivity or hyperedges */
  {"PHG_OUTPUT_LEVEL",                NULL,  "INT",    0},
    /* Higher value -> more output */
  {"FINAL_OUTPUT",                    NULL,  "INT",    0},
    /* Output related to final partitioning only */ 
  {"CHECK_GRAPH",                     NULL,  "INT",    0},
    /* Same as CHECK_HYPERGRAPH */
  {"CHECK_HYPERGRAPH",                NULL,  "INT",    0},
    /* Same as CHECK_GRAPH */
  {"PHG_NPROC_VERTEX",                NULL,  "INT",    0},
    /* No. of processors along vertex direction initially in 2D layout */
  {"PHG_NPROC_EDGE",                  NULL,  "INT",    0},
    /* No. of processors along hyperedge direction initially in 2D layout */
  {"PHG_COARSENING_LIMIT",            NULL,  "INT",    0},
    /* When to stop coarsening (global no. of vertices) */
  {"PHG_COARSENING_NCANDIDATE",       NULL,  "INT",    0},  
    /* Max no. of candidate vertices in a round of matching */
  {"PHG_COARSENING_METHOD",           NULL,  "STRING", 0},
    /* Coarsening method (ipm, agg, etc. ) */
  {"PHG_COARSENING_METHOD_FAST",      NULL,  "STRING", 0},
    /* Used by a-ipm to alternate between full ipm and a faster method */
  {"PHG_VERTEX_VISIT_ORDER",          NULL,  "INT",    0},
    /* Vertex ordering for greedy matching (linear, random,  ) */
  {"PHG_EDGE_SCALING",                NULL,  "INT",    0},
    /* Edge scaling schemes to tweak inner product similarity in matching */
  {"PHG_VERTEX_SCALING",              NULL,  "INT",    0},
    /* Vertex scaling schemes to tweak inner product similarity in matching */
  {"PHG_COARSEPARTITION_METHOD",      NULL,  "STRING", 0},
    /* Coarse partitioning method: linear, random, greedy, auto */
  {"PHG_REFINEMENT_METHOD",           NULL,  "STRING", 0},
    /* Only 2-way FM (fm2) for now */
  {"PHG_DIRECT_KWAY",                 NULL,  "INT",    0},
    /* Direct k-way partitioning not yet implemented! */
  {"PHG_REFINEMENT_LOOP_LIMIT",       NULL,  "INT",    0},
    /* Max no. of loops in KL/FM. */
  {"PHG_REFINEMENT_MAX_NEG_MOVE",     NULL,  "INT",    0},    
    /* Max. no. of negative moves allowed before exiting refinement. */
  {"PHG_REFINEMENT_QUALITY",          NULL,  "FLOAT",  0},
    /* 1.0 is default; higher (lower) value gives more (less) refinement. */
  {"PHG_USE_TIMERS",                  NULL,  "INT",    0},    
    /* Same as USE_TIMERS. */
  {"USE_TIMERS",                      NULL,  "INT",    0},    
    /* Same as PHG_USE_TIMERS. */
  {"PHG_EDGE_SIZE_THRESHOLD",         NULL,  "FLOAT",  0},
    /* Ignore hyperedges larger than this threshold times nvertex */
    /* If PHG_EDGE_SIZE_THRESHOLD>1, interpret it as absolute value. */
  {"PHG_MATCH_EDGE_SIZE_THRESHOLD",   NULL,  "INT",    0},
    /* Ignore hyperedges larger than this threshold, in local processor, during matching */
  {"PHG_BAL_TOL_ADJUSTMENT",          NULL,  "FLOAT",  0},  
    /* Adjustment factor for balance in recursive bisection. */
  {"PHG_EDGE_WEIGHT_OPERATION",       NULL,  "STRING",  0},
    /* How to handle inconsistent edge weights across processors */
  {"PARKWAY_SERPART",                 NULL,  "STRING", 0},
    /* Serial partitioner for ParKway (PaToH or hMetis) */
  {"ADD_OBJ_WEIGHT",                  NULL,  "STRING", 0},
    /* Add implicit vertex weight, like no. of pins (nonzeros)? */
  {"PHG_RANDOMIZE_INPUT",             NULL,  "INT",    0},    
    /* Randomizing input often improves load balance within PHG but destroys 
       locality, so may produce lower quality partitions  */
  {"PHG_PROCESSOR_REDUCTION_LIMIT",   NULL,  "FLOAT",  0},
    /* When to move data to fewer processors. */
  {"PHG_REPART_MULTIPLIER",           NULL,  "FLOAT",  0},
    /* Multiplier for communication to migration trade-off in repartitioning. */
  {"HYBRID_REDUCTION_FACTOR",        NULL,  "FLOAT",    0}, /* NEANEA */
    /* Factor by which to reduce the number of parts when using RCB matching. */
  {"PATOH_ALLOC_POOL0",               NULL,  "INT",    0},
    /* Memory allocation parameter for Patoh. */
  {"PATOH_ALLOC_POOL1",               NULL,  "INT",    0},   
  /* Memory allocation parameter for Patoh. */
#ifdef CEDRIC_2D_PARTITIONS
  {"PHG_KEEP_TREE",                   NULL,  "INT",    0},
  /* Keep dissection tree */
#endif /* CEDRIC_2D_PARTITIONS */
  {NULL,                              NULL,  NULL,     0}     
};

/* prototypes for static functions: */

static int Zoltan_PHG_Output_Parts(ZZ*, ZHG*, Partition);
static int Zoltan_PHG_Return_Lists(ZZ*, ZHG*, int*, ZOLTAN_ID_PTR*, 
  ZOLTAN_ID_PTR*, int**, int**);

#ifdef CHECK_LEFTALONE_VERTICES    
static int findAndSaveLeftAloneVertices(ZZ *zz, HGraph *hg, int p, 
                                 Partition parts,
                                 PHGPartParams *hgp) 
{
    char *yo="findAndSaveLeftAloneVertices";
    int *lneigh[2]={NULL, NULL}, *neigh[2]={NULL, NULL}, i, j, ierr=ZOLTAN_OK;
    int *lpins=NULL, *pins=NULL;
    PHGComm *hgc=hg->comm;

    if (hg->nEdge && (!(lpins = (int*) ZOLTAN_CALLOC(p * hg->nEdge, sizeof(int))) ||
                      !(pins  = (int*) ZOLTAN_MALLOC(p * hg->nEdge * sizeof(int)))))
        MEMORY_ERROR;

    for (i = 0; i < hg->nEdge; ++i)
        for (j = hg->hindex[i]; j < hg->hindex[i+1]; ++j)
            ++lpins[i*p+parts[hg->hvertex[j]]];
    if (hg->nEdge)
        MPI_Allreduce(lpins, pins, p*hg->nEdge, MPI_INT, MPI_SUM, 
                      hgc->row_comm);
    
    if (hg->nVtx && !(lneigh[0]  = (int*) ZOLTAN_MALLOC(2 * hg->nVtx * sizeof(int))))
        MEMORY_ERROR;
    if (!hgc->myProc_y) 
        if (hg->nVtx && !(neigh[0]  = (int*) ZOLTAN_MALLOC(2 * hg->nVtx * sizeof(int))))
            MEMORY_ERROR;

    if (hg->nVtx) {
        lneigh[1] = &(lneigh[0][hg->nVtx]);
        if (!hgc->myProc_y)
            neigh[1] = &(neigh[0][hg->nVtx]);
    }

    for (i = 0; i < hg->nVtx; ++i) {
        int pno = parts[i];
        lneigh[0][i] = lneigh[1][i] = 0;
        for (j = hg->vindex[i]; j < hg->vindex[i+1]; j++) {
            int edge = hg->vedge[j], k;
            lneigh[0][i] += (pins[edge*p+pno]-1); /* exclude the vertex's itself */
            for (k=0; k<p; ++k)
                if (k!=pno)
                    lneigh[1][i] += pins[edge*p+k];            
        }
    }
    
    if (hg->nVtx) 
        MPI_Reduce(lneigh[0], neigh[0], 2*hg->nVtx, MPI_INT, MPI_SUM, 0, hgc->col_comm);

    if (!hgc->myProc_y) {
        int alone=0, galone=0;        
        for (i=0; i<hg->nVtx; ++i)
            if (!neigh[0] && neigh[1]) {
                ++alone;
                if (alone<10)
                    uprintf(hgc, "vertex %d is alone in part %d but it has %d neighbours (!overcounted!) on other %d parts\n", i, parts[i], neigh[1]);
            }
        MPI_Reduce(&alone, &galone, 1, MPI_INT, MPI_SUM, 0, hgc->row_comm);
        if (!hgc->myProc)
            uprintf(hgc, "There are %d left-alone vertices\n", galone);        
    }
End:
    Zoltan_Multifree(__FILE__,__LINE__, 4, &lpins, &pins, &lneigh[0], &neigh[0]);
    return ierr;
}
#endif
 

/* UVCUVC DEBUG PRINT
static double detailed_balance_info(
  ZZ *zz,
  HGraph *hg,
  float *part_sizes,
  int p,
  Partition part
)
{
  int i;
  double *lsize_w, *size_w, max_imbal, tot_w;
  char *yo = "detailed_balance_info";
  PHGComm *hgc=NULL;
  int part_dim = (hg->VtxWeightDim ? hg->VtxWeightDim : 1);
  
  if (!hg || !hg->comm || !hg->comm->row_comm)  {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Unable to compute balance");
    return 1.0;
  }  

  hgc=hg->comm;  
  if (!(lsize_w = (double*) ZOLTAN_CALLOC (2*p, sizeof(double)))) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
  }
  size_w = lsize_w + p;
  
  if (hg->vwgt)
    for (i = 0; i < hg->nVtx; i++)
      lsize_w[part[i]] += hg->vwgt[i*hg->VtxWeightDim];
  else
    for (i = 0; i < hg->nVtx; i++)
      lsize_w[part[i]]++;
        
  MPI_Allreduce(lsize_w, size_w, p, MPI_DOUBLE, MPI_SUM, hg->comm->row_comm);
  
  max_imbal = tot_w = 0.0;
  for (i = 0; i < p; i++) 
      tot_w += size_w[i];
  if (tot_w) {
      if (!zz->Proc)
          uprintf(hgc, "PNO\tActual\tTarget\timbal\n");
      for (i = 0; i < p; i++) {
          float this_part_size = part_sizes[i*part_dim];
          if (this_part_size) {
              double ib=(size_w[i]-this_part_size*tot_w)/(this_part_size*tot_w);
              if (!zz->Proc)
                  uprintf(hgc, "%d\t%.1lf\t%.1lf\t%.3lf\n",
                          i, size_w[i], this_part_size*tot_w, 1.0+ib);
              if (ib>max_imbal)
                  max_imbal = ib;
          } else if (!zz->Proc)
              uprintf(hgc, "%d\t%.1lf\t%.1lf\tN/A\n", 
                      i, size_w[i], this_part_size*tot_w);
      }
      if (!zz->Proc)
          uprintf(hgc, "Max Imbal = %.3lf\n", 1.0+max_imbal);
  }

  ZOLTAN_FREE (&lsize_w);

  return  1.0+max_imbal;
}
*/

/******************************************************************************/
/* Main routine for Zoltan interface to hypergraph partitioning. Builds input */
/* data structures, set parameters, calls HG partitioner, builds return lists.*/
/* Type = ZOLTAN_LB_FN.                                                       */

int Zoltan_PHG(
ZZ *zz,                    /* The Zoltan structure  */
float *part_sizes,         /* Input:  Array of size zz->Num_Global_Parts 
                                containing the percentage of work assigned 
                                to each partition. */
int *num_imp,              /* not computed */
ZOLTAN_ID_PTR *imp_gids,   /* not computed */
ZOLTAN_ID_PTR *imp_lids,   /* not computed */
int **imp_procs,           /* not computed */
int **imp_to_part,         /* not computed */
int *num_exp,              /* number of objects to be exported */
ZOLTAN_ID_PTR *exp_gids,   /* global ids of objects to be exported */
ZOLTAN_ID_PTR *exp_lids,   /* local  ids of objects to be exported */
int **exp_procs,           /* list of processors to export to */
int **exp_to_part )         /* list of partitions to which exported objs
                                are assigned. */
{

  char *yo = "Zoltan_PHG";
  ZHG *zoltan_hg = NULL;
  PHGPartParams hgp;               /* Hypergraph parameters. */
  HGraph *hg = NULL;               /* Hypergraph itself */
  Partition parts = NULL;          /* Partition assignments in 
                                      2D distribution. */
  int err = ZOLTAN_OK, p=0;
  struct phg_timer_indices *timer = NULL; 
  int do_timing = 0;

#ifdef CEDRIC_2D_PARTITIONS
  Zoltan_PHG_LB_Data *data;

  int* sizeParts=NULL;
  int  numParts;
  struct Zoltan_DD_Struct *ddPartEdge=NULL;
#endif

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Initialization of return arguments. */
  *num_imp   = *num_exp   = -1;
  *imp_gids  = *exp_gids  = NULL;
  *imp_lids  = *exp_lids  = NULL;
  *imp_procs = *exp_procs = NULL;
  
  /* Initialize HG parameters. */
  err = Zoltan_PHG_Initialize_Params(zz, part_sizes, zz->LB.Method, &hgp);
  if (err != ZOLTAN_OK)
    goto End;

  if (hgp.use_timers) {
    if (!Zoltan_PHG_LB_Data_timers(zz))  {
      Zoltan_PHG_Timers_init(zz);
    }
    timer = Zoltan_PHG_LB_Data_timers(zz);
    if (timer->all < 0) 
      timer->all = Zoltan_Timer_Init(zz->ZTime, 1, "Zoltan_PHG");
  }

  if (hgp.use_timers > 1) {
    do_timing = 1;
    if (timer->build < 0) 
      timer->build = Zoltan_Timer_Init(zz->ZTime, 1, "Build");
    if (timer->setupvmap < 0) 
      timer->setupvmap = Zoltan_Timer_Init(zz->ZTime, 0, "Vmaps");
  }

  if (hgp.use_timers) 
    ZOLTAN_TIMER_START(zz->ZTime, timer->all, zz->Communicator);
    
  if (do_timing)
    ZOLTAN_TIMER_START(zz->ZTime, timer->build, zz->Communicator);
    
  /* build initial Zoltan hypergraph from callback functions. */

  err = Zoltan_PHG_Build_Hypergraph (zz, &zoltan_hg, &parts, &hgp);

  if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error building hypergraph.");
    goto End;
  }

  if (zoltan_hg->globalObj == 0){
    /* degenerate case - no objects to partition */
    hgp.final_output = 0;
    *num_exp = 0;
    goto End;
  }

  hg = &zoltan_hg->HG;
  p = zz->LB.Num_Global_Parts;  
  zoltan_hg->HG.redl = MAX(hgp.redl, p);     /* redl needs to be dynamic */

  /* RTHRTH -- redl may need to be scaled by number of procs */
  /* EBEB -- at least make sure redl > #procs */

  if (hgp.UseFixedVtx)
      hg->bisec_split = 1; /* this will be used only #parts=2
                              otherwise rdivide will set to appropriate
                              value */

  if (hgp.UsePrefPart || hgp.UseFixedVtx) { /* allocate memory for pref_part */
    if (hg->nVtx &&                       
        !(hg->pref_part = (int*) ZOLTAN_MALLOC (sizeof(int) * hg->nVtx))) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error building hypergraph.");
        goto End;
    }
  }
  if (hgp.UsePrefPart) /* copy input parts as pref_part;
                          UVCUVC: TODO: this code (and alloc of pref_part)
                          should go to Build_Hypergraph later */
      memcpy(hg->pref_part, parts, sizeof(int) * hg->nVtx);
  
  if (hgp.UseFixedVtx) {
      int i;
      for (i=0; i<hg->nVtx; ++i)
          if (hg->fixed_part[i]>=0)
              hg->pref_part[i] = hg->fixed_part[i];
          else if (!hgp.UsePrefPart)
              hg->pref_part[i] = -1;      
  }
  hgp.UsePrefPart |= hgp.UseFixedVtx;

  if (do_timing)
    ZOLTAN_TIMER_STOP(zz->ZTime, timer->build, zz->Communicator);


/*
  UVCUVC DEBUG PRINT
  uprintf(hg->comm, "Zoltan_PHG kway=%d #parts=%d\n", hgp.kway, zz->LB.Num_Global_Parts);
*/

  if (!strcasecmp(hgp.hgraph_pkg, "PARKWAY")){
    if (do_timing) {
      if (timer->parkway < 0)
        timer->parkway = Zoltan_Timer_Init(zz->ZTime, 0, "PHG_ParKway");
      ZOLTAN_TIMER_START(zz->ZTime, timer->parkway, zz->Communicator);
    }
    err = Zoltan_PHG_ParKway(zz, hg, p,
                             parts, &hgp);
    if (err != ZOLTAN_OK) 
        goto End;
    if (do_timing)
      ZOLTAN_TIMER_STOP(zz->ZTime, timer->parkway, zz->Communicator);
  } else if (!strcasecmp(hgp.hgraph_pkg, "PATOH")){
    if (hgp.use_timers > 1) {
      if (timer->patoh < 0)
        timer->patoh = Zoltan_Timer_Init(zz->ZTime, 0, "HG_PaToH");
      ZOLTAN_TIMER_START(zz->ZTime, timer->patoh, zz->Communicator);
    }
    err = Zoltan_PHG_PaToH(zz, hg, p,
                           parts, &hgp);
    if (err != ZOLTAN_OK) 
      goto End;
    if (hgp.use_timers > 1)
      ZOLTAN_TIMER_STOP(zz->ZTime, timer->patoh, zz->Communicator);
  }      
  else { /* it must be PHG  */
    /* Create tree structure */
    Zoltan_PHG_Tree_create(p, zz);

    /* UVC: if it is bisection anyways; no need to create vmap etc; 
       rdivide is going to call Zoltan_PHG_Partition anyways... */
    if (hgp.globalcomm.Communicator != MPI_COMM_NULL) {
      /* This processor is part of the 2D data distribution; it should
         participate in partitioning. */

      /* TODO : Construct a fake tree when there are only two parts */
      if (hgp.kway) { /* || zz->LB.Num_Global_Parts == 2)  */
        /* call main V cycle routine */
        err = Zoltan_PHG_Partition(zz, hg, p,
                                   hgp.part_sizes, parts, &hgp);

        if (err != ZOLTAN_OK) {
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error partitioning hypergraph.");
          goto End;
        }

      }
      else {
        int i;
          
        if (do_timing) 
          ZOLTAN_TIMER_START(zz->ZTime, timer->setupvmap, zz->Communicator);
        /* vmap associates original vertices to sub hypergraphs */
        if (hg->nVtx && 
            !(hg->vmap = (int*) ZOLTAN_MALLOC(hg->nVtx*sizeof (int))))  {
          err = ZOLTAN_MEMERR;
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
          goto End;
        }
        for (i = 0; i < hg->nVtx; ++i)
          hg->vmap[i] = i;
  
        if (do_timing) 
          ZOLTAN_TIMER_STOP(zz->ZTime, timer->setupvmap, zz->Communicator);


        /* partition hypergraph */
        err = Zoltan_PHG_rdivide (0, p-1, parts, zz, hg, &hgp, 0, 1);

        if (hgp.output_level >= PHG_DEBUG_LIST)     
          uprintf(hg->comm, "FINAL %3d |V|=%6d |E|=%6d #pins=%6d %s/%s/%s p=%d "
                  "bal=%.2f cutl=%.2f\n", 
                  hg->info, hg->nVtx, hg->nEdge, hg->nPins,
                  hgp.redm_str, hgp.coarsepartition_str, 
                  hgp.refinement_str, p,
                  Zoltan_PHG_Compute_Balance(zz, hg, hgp.part_sizes, 0,
                                             p, parts),
                  Zoltan_PHG_Compute_ConCut(hg->comm, hg, parts, p, &err));
            
        if (err != ZOLTAN_OK)  {
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error partitioning hypergraph.");
          goto End;
        }
        ZOLTAN_FREE (&hg->vmap);

      }

#ifdef CHECK_LEFTALONE_VERTICES
      findAndSaveLeftAloneVertices(zz, hg, p, parts, &hgp);
#endif
    }

#ifdef CEDRIC_2D_PARTITIONS
    if (hgp.keep_tree) {
      /* Build a centralized tree */
      Zoltan_PHG_Tree_centralize(zz);
      if (zoltan_hg->ddHedge != NULL) {
	Zoltan_PHG_2ways_hyperedge_partition (zz, hg, parts, get_tree(zz), zoltan_hg->ddHedge,
					      &ddPartEdge, &numParts, &sizeParts);
	Zoltan_DD_Destroy(&zoltan_hg->ddHedge);
	data = (Zoltan_PHG_LB_Data*)zz->LB.Data_Structure;
	data->ddHedge = ddPartEdge;
	data->numParts = numParts;
	data->sizeParts = sizeParts;
      }
      else
	data = NULL;
    }
    else
#endif /* CEDRIC_2D_PARTITIONS */
      Zoltan_PHG_LB_Data_free_tree(zz);

  }

  if (!strcasecmp(hgp.hgraph_method, "REPARTITION")) {
    Zoltan_PHG_Remove_Repart_Data(zz, zoltan_hg, hg, &hgp);
  }


/* UVC DEBUG PRINT
  if (!strcasecmp(hgp->hgraph_method, "REFINE")){
      uprintf(hg->comm, 
              "UVC ATTHEEND |V|=%6d |E|=%6d #pins=%6d p=%d bal=%.2f cutl=%.2f\n",
              hg->nVtx, hg->nEdge, hg->nPins, p,
              Zoltan_PHG_Compute_Balance(zz, hg, part_sizes, p, parts),
              Zoltan_PHG_Compute_ConCut(hg->comm, hg, parts, p, &err));
      detailed_balance_info(zz, hg, part_sizes, p, parts);
  }
*/

  
  /* Initialize these timers here so their output is near end of printout */
  if (do_timing)
    if (timer->retlist < 0) 
      timer->retlist = Zoltan_Timer_Init(zz->ZTime, 1, "Return_Lists");

  if (hgp.use_timers)
    if (timer->finaloutput < 0) 
      timer->finaloutput = Zoltan_Timer_Init(zz->ZTime, 1, "Final_Output");

  if (do_timing) 
    ZOLTAN_TIMER_START(zz->ZTime, timer->retlist, zz->Communicator);

  /* Build Zoltan's Output_Parts, mapped from 2D distribution 
     to input distribution. */

  Zoltan_PHG_Output_Parts(zz, zoltan_hg, parts);

  /* Build Zoltan's return arguments. */
  Zoltan_PHG_Return_Lists(zz, zoltan_hg, num_exp, exp_gids,
   exp_lids, exp_procs, exp_to_part);
    
  if (do_timing)
    ZOLTAN_TIMER_STOP(zz->ZTime, timer->retlist, zz->Communicator);

End:
  if (err == ZOLTAN_MEMERR)
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Memory error.")
  else if (err != ZOLTAN_OK)
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Error partitioning hypergraph.")
    
  /* KDDKDD The following code prints a final quality result even when
   * KDDKDD phg_output_level is zero.  It is useful for our tests and
   * KDDKDD data collection. */
  if ((err == ZOLTAN_OK) && hgp.final_output) {
    static int nRuns=0;
    static double balsum = 0.0, cutlsum = 0.0, cutnsum = 0.0, movesum = 0.0, repartsum = 0.0;
    static double balmax = 0.0, cutlmax = 0.0, cutnmax = 0.0, movemax = 0.0, repartmax = 0.0;
    static double balmin = 1e100, cutlmin = 1e100, cutnmin = 1e100, movemin = 1e100, repartmin = 1e100;
    double bal = 0.; 
    double cutl = 0.; /* Connnectivity cuts:  sum_over_edges((npart-1)*ewgt) */
    double cutn = 0.; /* Net cuts:  sum_over_edges((nparts>1)*ewgt) */

    double rlocal[2];  /* local cut stats for removed edges */
    double rglobal[2]; /* global cut stats for removed edges */
    int gnremove, i;
    double move=0.0, gmove;  /* local and global migration costs */
    double repart=0.0;   /* total repartitioning cost: comcost x multiplier + migration_cost */

    
    if (hgp.use_timers) {
      /* Do not include final output time in partitioning time */
      ZOLTAN_TIMER_STOP(zz->ZTime, timer->all, zz->Communicator);
      ZOLTAN_TIMER_START(zz->ZTime, timer->finaloutput, zz->Communicator);
    }

    if (hgp.globalcomm.Communicator != MPI_COMM_NULL) {
      /* Processor participated in partitioning */
      bal = Zoltan_PHG_Compute_Balance(zz, hg, hgp.part_sizes, 0,
                                       zz->LB.Num_Global_Parts, parts);
      cutl= Zoltan_PHG_Compute_ConCut(hg->comm, hg, parts,
                                      zz->LB.Num_Global_Parts, &err);
      cutn = Zoltan_PHG_Compute_NetCut(hg->comm, hg, parts,
                                       zz->LB.Num_Global_Parts);
      for (i = 0; i < zoltan_hg->nObj; ++i) {
        /* uprintf(hg->comm, " obj[%d] = %d  in=%d out=%d\n", i, zoltan_hg->AppObjSizes[i], zoltan_hg->Input_Parts[i], zoltan_hg->Output_Parts[i]); */
	if (zoltan_hg->Input_Parts[i] != zoltan_hg->Output_Parts[i])
            move += (double) ((zoltan_hg->AppObjSizes) ? zoltan_hg->AppObjSizes[i] : 1.0);

      }
    }

    if (!err) {
     
      /* Add in cut contributions from removed edges */
      MPI_Allreduce(&(zoltan_hg->nHedges), &gnremove, 1, MPI_INT, MPI_SUM,
                    zz->Communicator);
      if (gnremove) {
        err = Zoltan_PHG_Cuts(zz, zoltan_hg, rlocal);
        MPI_Allreduce(rlocal, rglobal, 2, MPI_DOUBLE,MPI_SUM,zz->Communicator);
        
        cutl += rglobal[0];
        cutn += rglobal[1];
      }

      MPI_Allreduce(&move, &gmove, 1, MPI_DOUBLE, MPI_SUM, zz->Communicator);

      repart = cutl*hgp.RepartMultiplier + gmove;
      repartsum += repart;
      if (repart > repartmax) repartmax = repart;
      if (repart < repartmin) repartmin = repart;
      movesum += gmove;
      if (gmove > movemax) movemax = gmove;
      if (gmove < movemin) movemin = gmove;
      cutlsum += cutl;
      if (cutl > cutlmax) cutlmax = cutl;
      if (cutl < cutlmin) cutlmin = cutl;
      cutnsum += cutn;
      if (cutn > cutnmax) cutnmax = cutn;
      if (cutn < cutnmin) cutnmin = cutn;
      balsum += bal;
      if (bal > balmax) balmax = bal;
      if (bal < balmin) balmin = bal;
      nRuns++;
   
      if (zz->Proc == 0) {
        uprintf(hg->comm, 
                "STATS Runs %d  bal  CURRENT %f  MAX %f  MIN %f  AVG %f\n", 
                nRuns, bal, balmax, balmin, balsum/nRuns);
        uprintf(hg->comm, 
                "STATS Runs %d  cutl CURRENT %f  MAX %f  MIN %f  AVG %f\n", 
                nRuns, cutl, cutlmax, cutlmin, cutlsum/nRuns);
        uprintf(hg->comm, 
                "STATS Runs %d  cutn CURRENT %f  MAX %f  MIN %f  AVG %f\n", 
                nRuns, cutn, cutnmax, cutnmin, cutnsum/nRuns);
	uprintf(hg->comm,
		"STATS Runs %d  %s CURRENT %f  MAX %f  MIN %f  AVG %f\n",
		nRuns, (zoltan_hg->showMoveVol) ? "moveVol" : "moveCnt", gmove, movemax, movemin, movesum/nRuns);
        if (zoltan_hg->showMoveVol) 
            uprintf(hg->comm,
		"STATS Runs %d  repart CURRENT %f  MAX %f  MIN %f  AVG %f\n",
		nRuns, repart, repartmax, repartmin, repartsum/nRuns);        
      }
    }

    if (hgp.use_timers) {
      ZOLTAN_TIMER_STOP(zz->ZTime, timer->finaloutput, zz->Communicator);
      ZOLTAN_TIMER_START(zz->ZTime, timer->all, zz->Communicator);
    }
  }
  /* KDDKDD  End of printing section. */
  
  ZOLTAN_FREE(&parts);
  if (zoltan_hg != NULL) {
    Zoltan_PHG_Free_Hypergraph_Data(zoltan_hg);
    ZOLTAN_FREE(&zoltan_hg);
  }

  if (hgp.use_timers) {
    ZOLTAN_TIMER_STOP(zz->ZTime, timer->all, zz->Communicator);
    if (hgp.globalcomm.Communicator != MPI_COMM_NULL)
      Zoltan_Timer_PrintAll(zz->ZTime, 0, hgp.globalcomm.Communicator, stdout);
  }

  if (hgp.globalcomm.row_comm != MPI_COMM_NULL)
    MPI_Comm_free(&(hgp.globalcomm.row_comm));
  if (hgp.globalcomm.col_comm != MPI_COMM_NULL)
    MPI_Comm_free(&(hgp.globalcomm.col_comm));
  if (hgp.globalcomm.Communicator != MPI_COMM_NULL)
    MPI_Comm_free(&(hgp.globalcomm.Communicator));

  /* Free part_sizes if created new due to ADD_OBJ_WEIGHT */
  if (hgp.part_sizes != part_sizes)
    ZOLTAN_FREE(&hgp.part_sizes);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return err;
}
/*****************************************************************************/

void Zoltan_PHG_Free_Hypergraph_Data(ZHG *zoltan_hg)
{
  if (zoltan_hg != NULL) {
    Zoltan_Input_HG_Free(zoltan_hg);
  }
}
    
/*****************************************************************************/

int Zoltan_PHG_Initialize_Params(
  ZZ *zz,   /* the Zoltan data structure */
  float *part_sizes, /* preallocation assumes object weight dimension is one */
  int hgraph_model,  /* GRAPH or HYPERGRAPH -- some default param values 
                        depend on the model. */
  PHGPartParams *hgp
)
{
  int err = ZOLTAN_OK;
  char *yo = "Zoltan_PHG_Initialize_Params";
  int nProc;
  int usePrimeComm;
  MPI_Comm communicator;
  char add_obj_weight[MAX_PARAM_STRING_LEN];
  char edge_weight_op[MAX_PARAM_STRING_LEN];
  char cut_objective[MAX_PARAM_STRING_LEN];
  char *package = hgp->hgraph_pkg; 
  char *method = hgp->hgraph_method;
  char buf[1024];

  memset(hgp, 0, sizeof(*hgp)); /* in the future if we forget to initialize
                                   another param at least it will be 0 */
  
  Zoltan_Bind_Param(PHG_params, "HYPERGRAPH_PACKAGE", &hgp->hgraph_pkg);
  Zoltan_Bind_Param(PHG_params, "PHG_MULTILEVEL", &hgp->useMultilevel);
  Zoltan_Bind_Param(PHG_params, "PHG_OUTPUT_LEVEL", &hgp->output_level);
  Zoltan_Bind_Param(PHG_params, "FINAL_OUTPUT", &hgp->final_output); 
  Zoltan_Bind_Param(PHG_params, "CHECK_GRAPH", &hgp->check_graph);   
  Zoltan_Bind_Param(PHG_params, "CHECK_HYPERGRAPH", &hgp->check_graph);   
  Zoltan_Bind_Param(PHG_params, "PHG_NPROC_VERTEX", &hgp->nProc_x_req);
  Zoltan_Bind_Param(PHG_params, "PHG_NPROC_EDGE", &hgp->nProc_y_req);
  Zoltan_Bind_Param(PHG_params, "PHG_COARSENING_LIMIT", &hgp->redl);
  Zoltan_Bind_Param(PHG_params, "PHG_COARSENING_NCANDIDATE", &hgp->nCand);
  Zoltan_Bind_Param(PHG_params, "PHG_COARSENING_METHOD", hgp->redm_str);
  Zoltan_Bind_Param(PHG_params, "PHG_COARSENING_METHOD_FAST", hgp->redm_fast);
  Zoltan_Bind_Param(PHG_params, "PHG_VERTEX_VISIT_ORDER", &hgp->visit_order);
  Zoltan_Bind_Param(PHG_params, "PHG_EDGE_SCALING", &hgp->edge_scaling);
  Zoltan_Bind_Param(PHG_params, "PHG_VERTEX_SCALING", &hgp->vtx_scaling);
  Zoltan_Bind_Param(PHG_params, "PHG_REFINEMENT_METHOD", hgp->refinement_str);
  Zoltan_Bind_Param(PHG_params, "PHG_DIRECT_KWAY", &hgp->kway);
  Zoltan_Bind_Param(PHG_params, "HYBRID_REDUCTION_FACTOR", &hgp->rcb_red); /* NEANEA */
#ifdef CEDRIC_2D_PARTITIONS
  Zoltan_Bind_Param(PHG_params, "PHG_KEEP_TREE", &hgp->keep_tree);
#endif /* CEDRIC_2D_PARTITIONS */
  Zoltan_Bind_Param(PHG_params, "PHG_REFINEMENT_LOOP_LIMIT", 
                                &hgp->fm_loop_limit);
  Zoltan_Bind_Param(PHG_params, "PHG_REFINEMENT_MAX_NEG_MOVE", 
                                &hgp->fm_max_neg_move);  
  Zoltan_Bind_Param(PHG_params, "PHG_REFINEMENT_QUALITY", 
                                &hgp->refinement_quality);  
  Zoltan_Bind_Param(PHG_params, "PHG_COARSEPARTITION_METHOD", 
                                 hgp->coarsepartition_str);
  Zoltan_Bind_Param(PHG_params, "PHG_USE_TIMERS",
                                 (void*) &hgp->use_timers);  
  Zoltan_Bind_Param(PHG_params, "USE_TIMERS",
                                 (void*) &hgp->use_timers);  
  Zoltan_Bind_Param(PHG_params, "PHG_EDGE_SIZE_THRESHOLD",
                                 (void*) &hgp->EdgeSizeThreshold);  
  Zoltan_Bind_Param(PHG_params, "PHG_MATCH_EDGE_SIZE_THRESHOLD",
                                 (void*) &hgp->MatchEdgeSizeThreshold);  
  Zoltan_Bind_Param(PHG_params, "PHG_BAL_TOL_ADJUSTMENT",
                                 (void*) &hgp->bal_tol_adjustment);  
  Zoltan_Bind_Param(PHG_params, "PARKWAY_SERPART",
                                 (void *) hgp->parkway_serpart);
  Zoltan_Bind_Param(PHG_params, "PHG_CUT_OBJECTIVE",
                                 (void *) &cut_objective);
  Zoltan_Bind_Param(PHG_params, "ADD_OBJ_WEIGHT",
                                 (void *) add_obj_weight);
  Zoltan_Bind_Param(PHG_params, "PHG_EDGE_WEIGHT_OPERATION",
                                 (void *) edge_weight_op);
  Zoltan_Bind_Param(PHG_params, "PHG_RANDOMIZE_INPUT",
                                 (void*) &hgp->RandomizeInitDist);  
  Zoltan_Bind_Param(PHG_params, "PHG_PROCESSOR_REDUCTION_LIMIT",
		                 (void*) &hgp->ProRedL);
  Zoltan_Bind_Param(PHG_params, "PHG_REPART_MULTIPLIER",
		                 (void*) &hgp->RepartMultiplier);
  Zoltan_Bind_Param(PHG_params, "PATOH_ALLOC_POOL0",
                                 (void*) &hgp->patoh_alloc_pool0);
  Zoltan_Bind_Param(PHG_params, "PATOH_ALLOC_POOL1",
                                 (void*) &hgp->patoh_alloc_pool1);
  
  
  /* Set default values */
  strncpy(hgp->hgraph_pkg,           "phg", MAX_PARAM_STRING_LEN);
  strncpy(hgp->redm_str,             "agg", MAX_PARAM_STRING_LEN);
  hgp->match_array_type = 0;
  strncpy(hgp->redm_fast,          "l-ipm", MAX_PARAM_STRING_LEN);
  strncpy(hgp->coarsepartition_str, "auto", MAX_PARAM_STRING_LEN);
  strncpy(hgp->refinement_str,       "fm2", MAX_PARAM_STRING_LEN);
  strncpy(hgp->parkway_serpart,    "patoh", MAX_PARAM_STRING_LEN);
  strncpy(cut_objective,    "connectivity", MAX_PARAM_STRING_LEN);
  strncpy(add_obj_weight,           "none", MAX_PARAM_STRING_LEN);

  if (hgraph_model == GRAPH)
    strncpy(edge_weight_op,          "sum", MAX_PARAM_STRING_LEN);
  else 
    strncpy(edge_weight_op,          "max", MAX_PARAM_STRING_LEN);

  /* LB.Approach is initialized to "REPARTITION", and set in Set_Key_Params  */
  strncpy(hgp->hgraph_method,  zz->LB.Approach, MAX_PARAM_STRING_LEN);
  if (!strcasecmp(zz->LB.Approach,"REFINE")) 
    hgp->useMultilevel = 0;
  else
    hgp->useMultilevel = 1;

  hgp->use_timers = 0;
  hgp->LocalCoarsePartition = 0;
  hgp->edge_scaling = 0;
  hgp->vtx_scaling = 0;
  hgp->vtx_scal_size = 0;
  hgp->vtx_scal = NULL;  /* Array for storing vertex degree scale vector. 
                            Should perhaps go in hg structure, not the
                            param struct? */
  hgp->connectivity_cut = 1; 
  hgp->visit_order = 0;  /* Random */
  hgp->check_graph = 0;
  hgp->bal_tol = zz->LB.Imbalance_Tol[0]; /* Make vector for multiconstraint */
  hgp->bal_tol_adjustment = 0.7;
  hgp->nCand = 100;
  hgp->redl = MAX(2*zz->LB.Num_Global_Parts, 100);
  hgp->output_level = PHG_DEBUG_NONE;
  hgp->final_output = 0;
  hgp->nProc_x_req = -1;
  hgp->nProc_y_req = -1;
  hgp->kway = 0;
  hgp->rcb_red == 0.5; /* NEANEA probably change */
  hgp->fm_loop_limit = 10;
  hgp->fm_max_neg_move = 250;  
  hgp->refinement_quality = 1;
  hgp->RandomizeInitDist = 0;
  hgp->EdgeSizeThreshold = 0.25;
  hgp->MatchEdgeSizeThreshold = 500;  
  hgp->hybrid_keep_factor = 0.;
  hgp->ProRedL = 0.0; /* UVCUVC: CHECK default set to 0 until we run more experiments */
  hgp->RepartMultiplier = 100.;
  hgp->patoh_alloc_pool0 = 0;
  hgp->patoh_alloc_pool1 = 0;
  hgp->UseFixedVtx = 0;
  hgp->UsePrefPart = 0;
  
  /* Get application values of parameters. */
  err = Zoltan_Assign_Param_Vals(zz->Params, PHG_params, zz->Debug_Level, 
          zz->Proc, zz->Debug_Proc);
  
  nProc = zz->Num_Proc;
  usePrimeComm = 0;

  /* Parse add_obj_weight parameter - PHG only processes one weight */

  hgp->part_sizes = part_sizes;

  if (!strcasecmp(add_obj_weight, "none")) {
    hgp->add_obj_weight = PHG_ADD_NO_WEIGHT;
  }
  else if (zz->Obj_Weight_Dim > 0) {
    /* Do not add_obj_weight until multiconstraint PHG is implemented */
    if (zz->Proc == 0){
      ZOLTAN_PRINT_WARN(zz->Proc, yo,
       "Both application supplied *and* ADD_OBJ_WEIGHT "
       "calculated vertex weights were provided.");
      ZOLTAN_PRINT_WARN(zz->Proc, yo,
        "Only the first application supplied weight per vertex will be used.");
    }
    hgp->add_obj_weight = PHG_ADD_NO_WEIGHT;
  } 
  else {
    if (!strcasecmp(add_obj_weight, "vertices")){
      hgp->add_obj_weight = PHG_ADD_UNIT_WEIGHT;
    } else if (!strcasecmp(add_obj_weight, "unit")){
      hgp->add_obj_weight = PHG_ADD_UNIT_WEIGHT;
    } else if (!strcasecmp(add_obj_weight, "vertex degree")){
      hgp->add_obj_weight = PHG_ADD_PINS_WEIGHT;
    } else if (!strcasecmp(add_obj_weight, "nonzeros")){
      hgp->add_obj_weight = PHG_ADD_PINS_WEIGHT;
    } else if (!strcasecmp(add_obj_weight, "pins")){
      hgp->add_obj_weight = PHG_ADD_PINS_WEIGHT;
    } else{
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid ADD_OBJ_WEIGHT parameter.\n");
      err = ZOLTAN_WARN;
    }
  }

  if ((zz->Obj_Weight_Dim==0) &&      /* no application supplied weights */
      (hgp->add_obj_weight==PHG_ADD_NO_WEIGHT)){ /* no calculated weight */

    hgp->add_obj_weight = PHG_ADD_UNIT_WEIGHT; /* default object weight */
  }

  if (!strcasecmp(cut_objective, "default")
      || !strcasecmp(cut_objective, "connectivity"))
      hgp->connectivity_cut = 1;
  else if (!strcasecmp(cut_objective, "hyperedges"))
      hgp->connectivity_cut = 0;
  else {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid PHG_CUT_OBJECTIVE parameter.\n");
      goto End;
  }

  if (!strcasecmp(edge_weight_op, "max")){
    hgp->edge_weight_op = PHG_MAX_EDGE_WEIGHTS;
  } else if (!strcasecmp(edge_weight_op, "add")){
    hgp->edge_weight_op = PHG_ADD_EDGE_WEIGHTS;
  } else if (!strcasecmp(edge_weight_op, "sum")){
    hgp->edge_weight_op = PHG_ADD_EDGE_WEIGHTS;
  } else if (!strcasecmp(edge_weight_op, "error")){
    hgp->edge_weight_op = PHG_FLAG_ERROR_EDGE_WEIGHTS;
  } else{
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,
      "Invalid PHG_EDGE_WEIGHT_OPERATION parameter.  Zoltan will use \"max\".\n");
    err = ZOLTAN_WARN;
  }

  if ((strcasecmp(method, "PARTITION")) &&
      (strcasecmp(method, "REPARTITION")) &&
      (strcasecmp(method, "REFINE"))) {
    sprintf(buf,"%s is not a valid hypergraph method\n",method);
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, buf);
    err = ZOLTAN_FATAL;
    goto End;
  }

  /* Adjust refinement parameters using hgp->refinement_quality */
  if (hgp->refinement_quality < 0.5/hgp->fm_loop_limit) 
    /* No refinement */
    strncpy(hgp->refinement_str,      "no",   MAX_PARAM_STRING_LEN);
  else {
    /* Scale FM parameters */
    hgp->fm_loop_limit   *= hgp->refinement_quality;
    hgp->fm_max_neg_move *= hgp->refinement_quality;
  }

  if (!strcasecmp(package, "PHG")){
    /* Test to determine whether we should change the number of processors
       used for partitioning to make more efficient 2D decomposition */

    if (hgp->nProc_x_req != 1 && hgp->nProc_y_req != 1)  /* Want 2D decomp */
      if (zz->Num_Proc > SMALL_PRIME && Zoltan_PHG_isPrime(zz->Num_Proc)) 
        /* 2D data decomposition is requested but we have a prime 
         * number of processors. */
        usePrimeComm = 1;

    if ((!strcasecmp(method, "REPARTITION"))){
        zz->LB.Remap_Flag = 0;
    }

    if ((!strcasecmp(method, "REPARTITION")) ||
        (!strcasecmp(method, "REFINE"))) {
        hgp->fm_loop_limit = 4; /* experimental evaluation showed that for
                                repartitioning/refinement small number of passes
                                is "good enough". These are all heuristics hence
                                it is possible to create a pathological cases; 
                                but in general this seems to be sufficient */
    }
    
    if (!hgp->useMultilevel) {
        /* don't do coarsening */
        strncpy(hgp->redm_str, "no", MAX_PARAM_STRING_LEN);

        /* we have modified all coarse partitioners to handle preferred part
           if user wants to choose one she can choose; otherwise default 
           partitioner
           (greedy growing) does work better than previous default partitioning
           for phg_refine ("no"). */        
        hgp->UsePrefPart = 1;

    }
    if (!strcasecmp(method, "REFINE") && hgp->useMultilevel){
        /* UVCUVC: as a heuristic we prefer local matching;
           in our experiments for IPDPS'07 and WileyChapter multilevel_refine
           didn't prove itself useful; it is too costly even with local matching
           hence it will not be be released yet (i.e. not in v3). */
        strncpy(hgp->redm_str, "l-ipm", MAX_PARAM_STRING_LEN);                
        hgp->UsePrefPart = 1;
    }    
  }
  else if (!strcasecmp(package, "PARKWAY")){
    if (hgp->nProc_x_req>1) {
      err = ZOLTAN_FATAL;
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "ParKway requires nProc_x=1 or -1.");
      goto End;
    }
    hgp->nProc_x_req = 1;
  } 
  else if (!strcasecmp(package, "PATOH")){
    if (zz->Num_Proc>1) {
      err = ZOLTAN_FATAL;
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "PaToH only works with Num_Proc=1.");
      goto End;
    }
  }

  if (!usePrimeComm)
    MPI_Comm_dup(zz->Communicator, &communicator);
  else {
    MPI_Group newgrp, zzgrp;
    nProc--;
    MPI_Comm_group(zz->Communicator, &zzgrp);
    MPI_Group_excl(zzgrp, 1, &nProc, &newgrp);
    MPI_Comm_create(zz->Communicator, newgrp, &communicator);
    MPI_Group_free(&newgrp);
    MPI_Group_free(&zzgrp);
  }

  err = Zoltan_PHG_Set_2D_Proc_Distrib(zz, communicator, zz->Proc, 
                                       nProc, hgp->nProc_x_req, 
                                       hgp->nProc_y_req, 
                                       &hgp->globalcomm);
  if (err != ZOLTAN_OK) 
    goto End;

  /* Convert strings to function pointers. */
  err = Zoltan_PHG_Set_Part_Options (zz, hgp);
  
End:
  return err;
}

/*****************************************************************************/

int Zoltan_PHG_Set_Param(
  char *name,                     /* name of variable */
  char *val)                      /* value of variable */
{
  /* associates value to named variable for PHG partitioning parameters */
  PARAM_UTYPE result;         /* value returned from Check_Param */
  int index;                  /* index returned from Check_Param */
  int status;
  int i;

  char *valid_pkg[] = {
        "ZOLTAN", "PHG", "PATOH",
         NULL };

  status = Zoltan_Check_Param(name, val, PHG_params, &result, &index);

  if (status == 0){
    /* OK so far, do sanity check of parameter values */

    if (strcasecmp(name, "HYPERGRAPH_PACKAGE") == 0){
      status = 2;
      for (i=0; valid_pkg[i] != NULL; i++){
        if (strcasecmp(val, valid_pkg[i]) == 0){
          status = 0;
          break;
        }
      }
    }
  }
  return(status);
}

/****************************************************************************/

static int Zoltan_PHG_Output_Parts (
  ZZ *zz,
  ZHG *zhg,
  Partition hg_parts   /* Output partitions relative to the 2D distribution
                          of zhg->HG */
)
{
/* Function to map the computed partition from the distribution in HGraph
 * to the input distribution 
 */

static char *yo = "Zoltan_PHG_Output_Parts";
int i;
int msg_tag = 31000;
int ierr = ZOLTAN_OK;
int nObj = zhg->nObj;
int *outparts = NULL;
int *sendbuf = NULL;  
HGraph *phg = &(zhg->HG);

  zhg->Output_Parts = outparts 
                     = (int*) ZOLTAN_MALLOC (nObj * sizeof(int));
  if (zhg->VtxPlan != NULL) {
    /* Get the partition information from the 2D decomposition back to the
     * original owning processor for each GID.  */
    sendbuf = (int*) ZOLTAN_MALLOC(zhg->nRecv_GNOs * sizeof(int));
    for (i = 0; i < zhg->nRecv_GNOs; i++)
      sendbuf[i] = hg_parts[VTX_GNO_TO_LNO(phg, zhg->Recv_GNOs[i])];
    ierr = Zoltan_Comm_Do_Reverse(zhg->VtxPlan, msg_tag, (char*) sendbuf,
                                  sizeof(int), NULL, (char *) outparts);
    if (ierr) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error from Zoltan_Comm_Do_Reverse");
      goto End;
    }

    ZOLTAN_FREE(&sendbuf);
    Zoltan_Comm_Destroy(&(zhg->VtxPlan));
  }
  else {
    for (i = 0; i < zhg->nRecv_GNOs; i++)
      outparts[i] = hg_parts[zhg->Recv_GNOs[i]];
  }

  if (zz->LB.Remap_Flag) {
    int new_map;
    int *newproc = (int *) ZOLTAN_MALLOC(nObj * sizeof(int));
    int num_gid_entries = zz->Num_GID;

    if (nObj && !newproc) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
    for (i = 0; i < nObj; i++){
      newproc[i] = Zoltan_LB_Part_To_Proc(zz, outparts[i],
                                          &(zhg->objGID[i*num_gid_entries]));
      if (newproc[i]<0){
        ZOLTAN_PRINT_ERROR(zz->Proc, yo,
         "Zoltan_LB_Part_To_Proc returned invalid processor number.");
        ierr = ZOLTAN_FATAL;
        ZOLTAN_FREE(&newproc);
        goto End;
      }
    }
    
    ierr = Zoltan_LB_Remap(zz, &new_map, nObj, newproc, zhg->Input_Parts,
                           outparts, 1);
    if (ierr < 0) 
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_LB_Remap");
    ZOLTAN_FREE(&newproc);
  }

End:
  if (zhg->Recv_GNOs) ZOLTAN_FREE(&(zhg->Recv_GNOs));
  zhg->nRecv_GNOs = 0;
  return ierr;
}

/****************************************************************************/

static int Zoltan_PHG_Return_Lists (
  ZZ *zz,
  ZHG *zhg,
  int *num_exp,
  ZOLTAN_ID_PTR *exp_gids,
  ZOLTAN_ID_PTR *exp_lids,
  int **exp_procs,
  int **exp_to_part)
{
/* Routine to build export lists of ZOLTAN_LB_FN. */
char *yo = "Zoltan_PHG_Return_Lists";
int i, j;
int ierr = ZOLTAN_OK;
int eproc;
int num_gid_entries   = zz->Num_GID;
int num_lid_entries   = zz->Num_LID;
int nObj              = zhg->nObj;
Partition input_parts = zhg->Input_Parts;
ZOLTAN_ID_PTR gids    = zhg->objGID;
ZOLTAN_ID_PTR lids    = zhg->objLID;
int *outparts         = zhg->Output_Parts; 

  if (zz->LB.Return_Lists == ZOLTAN_LB_NO_LISTS) 
    goto End;
  else if (zz->LB.Return_Lists == ZOLTAN_LB_CANDIDATE_LISTS) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Candidate Lists not supported in PHG."
                                     "change RETURN_LISTS parameter");
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  /* Count number of objects with new partitions or new processors. */
  *num_exp = 0;
  for (i = 0; i < nObj; i++) {
    eproc = Zoltan_LB_Part_To_Proc(zz, outparts[i], &gids[i*num_gid_entries]);
    if (outparts[i] != input_parts[i] || zz->Proc != eproc)
      (*num_exp)++;
  }

  /* Allocate memory for return lists. */
  if (*num_exp > 0) {
    if (!Zoltan_Special_Malloc(zz, (void**)exp_gids, *num_exp,
                               ZOLTAN_SPECIAL_MALLOC_GID)
     || !Zoltan_Special_Malloc(zz, (void**)exp_lids, *num_exp,
                               ZOLTAN_SPECIAL_MALLOC_LID)
     || !Zoltan_Special_Malloc(zz, (void**)exp_procs, *num_exp,
                               ZOLTAN_SPECIAL_MALLOC_INT)
     || !Zoltan_Special_Malloc(zz, (void**)exp_to_part, *num_exp,
                               ZOLTAN_SPECIAL_MALLOC_INT)) {
        Zoltan_Special_Free(zz,(void**)exp_gids,   ZOLTAN_SPECIAL_MALLOC_GID);
        Zoltan_Special_Free(zz,(void**)exp_lids,   ZOLTAN_SPECIAL_MALLOC_LID);
        Zoltan_Special_Free(zz,(void**)exp_procs,  ZOLTAN_SPECIAL_MALLOC_INT);
        Zoltan_Special_Free(zz,(void**)exp_to_part,ZOLTAN_SPECIAL_MALLOC_INT);
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
        ierr = ZOLTAN_MEMERR;
        goto End;
     }

    for (j = 0, i = 0; i < nObj; i++) {
      eproc = Zoltan_LB_Part_To_Proc(zz, outparts[i], &gids[i*num_gid_entries]);
      if (outparts[i] != input_parts[i] || eproc != zz->Proc) {
        ZOLTAN_SET_GID(zz, &((*exp_gids)[j*num_gid_entries]),
                           &(gids[i*num_gid_entries]));
        if (num_lid_entries > 0)
          ZOLTAN_SET_LID(zz, &((*exp_lids)[j*num_lid_entries]),
                             &(lids[i*num_lid_entries]));
        (*exp_procs)  [j] = eproc;
        (*exp_to_part)[j] = outparts[i];
        j++;
      }
    }
  }

End:

  return ierr;
}


/*****************************************************************************/

int Zoltan_PHG_Set_2D_Proc_Distrib(
    ZZ *zz,                /* Input:  ZZ struct; for debuging   */
    MPI_Comm Communicator, /* Input:  The MPI Communicator; this communicator
                                      may be MPI_COMM_NULL, as PHG_Redistribute
                                      uses this function with MPI_COMM_NULL
                                      to compute nProc_x and nProc_y.  */
    int proc,              /* Input:  Rank of current processor */
    int nProc,             /* Input:  Total # of processors     */    
    int nProc_x,           /* Input:  Suggested #procs in x-direction */
    int nProc_y,           /* Input:  Suggested #procs in y-direction */
    PHGComm *comm          /* Ouput: filled */
    )    
{
/* Computes the processor distribution for the 2D data distrib.
 * Sets nProc_x, nProc_y.
 * Constraint:  nProc_x * nProc_y == nProc. 
 * For 2D data distrib, default should approximate sqrt(nProc).
 * If nProc_x and nProc_y both equal -1 on input, compute default.
 * Otherwise, compute valid values and/or return error.
 */
char *yo = "Zoltan_PHG_Set_2D_Proc_Distrib";
int tmp;
int ierr = ZOLTAN_OK;
    
  if (nProc_x == -1 && nProc_y == -1) {
    /* Compute default */
    tmp = (int) sqrt((double)nProc+0.1);
    while (nProc % tmp) tmp--;
    comm->nProc_x = tmp;
    comm->nProc_y = nProc / tmp;
  } else if (nProc_x == -1) {
    comm->nProc_y = MIN(nProc_y, nProc);
    comm->nProc_x = nProc / comm->nProc_y;
  } else if (nProc_y == -1) {
    comm->nProc_x = MIN(nProc_x, nProc);
    comm->nProc_y = nProc / comm->nProc_x;
  } else {
    comm->nProc_x = nProc_x;
    comm->nProc_y = nProc_y;    
  }
    
  /* Error check */
  if (comm->nProc_x * comm->nProc_y != nProc) {
    ZOLTAN_PRINT_ERROR(proc, yo,
                       "Values for PHG_NPROC_VERTEX and PHG_NPROC_EDGE "
                       "do not evenly divide the "
                       "total number of processors.");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  comm->nProc = nProc;
  comm->Communicator = Communicator;
  comm->zz = zz;

  if (Communicator==MPI_COMM_NULL) {
    comm->myProc_x = -1;
    comm->myProc_y = -1;
    comm->myProc = -1;
    comm->col_comm = comm->row_comm = MPI_COMM_NULL;
  } else {
    comm->myProc_x = proc % comm->nProc_x;
    comm->myProc_y = proc / comm->nProc_x;
    comm->myProc = proc;
    if ((MPI_Comm_split(Communicator, comm->myProc_x, comm->myProc_y, 
                        &comm->col_comm) != MPI_SUCCESS)
     || (MPI_Comm_split(Communicator, comm->myProc_y, comm->myProc_x, 
                        &comm->row_comm) != MPI_SUCCESS)) {
      ZOLTAN_PRINT_ERROR(proc, yo, "MPI_Comm_Split failed");
      return ZOLTAN_FATAL;
    }
    Zoltan_Srand_Sync(Zoltan_Rand(NULL), &(comm->RNGState_row),
                      comm->row_comm);
    Zoltan_Srand_Sync(Zoltan_Rand(NULL), &(comm->RNGState_col),
                      comm->col_comm);
    Zoltan_Srand_Sync(Zoltan_Rand(NULL), &(comm->RNGState),
                      comm->Communicator);
  } 
/*  printf("(%d, %d) of [%d, %d] -> After Comm_split col_comm=%d  row_comm=%d\n", hgp->myProc_x, hgp->myProc_y, hgp->nProc_x, hgp->nProc_y, (int)hgp->col_comm, (int)hgp->row_comm);  */
  
    
End:

  return ierr;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

