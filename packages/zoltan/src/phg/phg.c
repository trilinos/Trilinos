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
#include "params_const.h"
#include "all_allo_const.h"


/*
 *  Main routine for Zoltan interface to hypergraph partitioning. 
 *  Also routines that build input data structures, set parameters, etc.
 */


/******************************************************************************/
/*  Parameters structure for parallel HG method.  */
static PARAM_VARS PHG_params[] = {
  /* Add parameters here. */
  {"PHG_OUTPUT_LEVEL",                NULL,  "INT",    0},
  {"PHG_FINAL_OUTPUT",                NULL,  "INT",    0},
  {"FINAL_OUTPUT",                    NULL,  "INT",    0},
  {"CHECK_GRAPH",                     NULL,  "INT",    0},
  {"PHG_NPROC_X",                     NULL,  "INT",    0},
  {"PHG_NPROC_Y",                     NULL,  "INT",    0},
  {"PHG_PROC_SPLIT",                  NULL,  "INT",    0},  
  {"PHG_REDUCTION_LIMIT",             NULL,  "INT",    0},
  {"PHG_REDUCTION_METHOD",            NULL,  "STRING", 0},
  {"PHG_REDUCTION_LOCAL_IMPROVEMENT", NULL,  "STRING", 0},
  {"PHG_VERTEX_VISIT_ORDER",          NULL,  "INT",    0},
  {"PHG_EDGE_SCALING",                NULL,  "INT",    0},
  {"PHG_VERTEX_SCALING",              NULL,  "INT",    0},
  {"PHG_COARSE_PARTITIONING",         NULL,  "STRING", 0},
  {"PHG_REFINEMENT",                  NULL,  "STRING", 0},
  {"PHG_DIRECT_KWAY",                 NULL,  "INT",    0},
  {"PHG_FM_LOOP_LIMIT",               NULL,  "INT",    0},
  {"PHG_FM_MAX_NEG_MOVE",             NULL,  "INT",    0},    
  {"PHG_COARSE_ITERATIONS",           NULL,  "INT",    0},    
  {"PHG_USE_TIMERS",                  NULL,  "INT",    0},    
  {"USE_TIMERS",                      NULL,  "INT",    0},    
  {"EDGE_SIZE_THRESHOLD",             NULL,  "FLOAT",    0},    
  {NULL,                              NULL,  NULL,     0}     
};

/* prototypes for static functions: */

static int Zoltan_PHG_Initialize_Params(ZZ*, float *, PHGPartParams*);
static int Zoltan_PHG_Output_Parts(ZZ*, ZHG*, Partition);
static int Zoltan_PHG_Return_Lists(ZZ*, ZHG*, int*, ZOLTAN_ID_PTR*, 
  ZOLTAN_ID_PTR*, int**, int**);

 
 
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
  int nVtx;                        /* Temporary variable for base graph. */
  PHGPartParams hgp;               /* Hypergraph parameters. */
  Partition parts = NULL;          /* Partition assignments in 
                                      2D distribution. */
  int err = ZOLTAN_OK;
  static int timer_all = -1;       /* Note:  this timer includes other
                                      timers and their synchronization time,
                                      so it will be a little high. */
  static int timer_build=-1;       /* timers to be used in this function;
                                      declared static so that, over multiple
                                      runs, can accumulate times.  */
  static int timer_retlist=-1;

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Initialization of return arguments. */
  *num_imp   = *num_exp   = -1;
  *imp_gids  = *exp_gids  = NULL;
  *imp_lids  = *exp_lids  = NULL;
  *imp_procs = *exp_procs = NULL;
  
  /* Initialize HG parameters. */
  err = Zoltan_PHG_Initialize_Params (zz, part_sizes, &hgp);
  if (err != ZOLTAN_OK)
    goto End;

  if (hgp.use_timers) {
    if (timer_all < 0) 
      timer_all = Zoltan_Timer_Init(zz->ZTime, 1, "Zoltan_PHG");
    ZOLTAN_TIMER_START(zz->ZTime, timer_all, zz->Communicator);
  }
    
  if (hgp.use_timers > 1) {
    if (timer_build < 0) 
      timer_build = Zoltan_Timer_Init(zz->ZTime, 1, "Build");
    ZOLTAN_TIMER_START(zz->ZTime, timer_build, zz->Communicator);
  }
    
  /* build initial Zoltan hypergraph from callback functions. */
  err = Zoltan_PHG_Build_Hypergraph (zz, &zoltan_hg, &parts, &hgp);
  if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error building hypergraph.");
    goto End;
  }

  if (hgp.use_timers > 1)
    ZOLTAN_TIMER_STOP(zz->ZTime, timer_build, zz->Communicator);
   
  zz->LB.Data_Structure = zoltan_hg;
  nVtx = zoltan_hg->HG.nVtx;
  zoltan_hg->HG.redl = hgp.redl;     /* redl needs to be dynamic */
  /* RTHRTH -- redl may need to be scaled by number of procs */
  /* EBEB -- at least make sure redl > #procs */
 
/*
  uprintf(zoltan_hg->HG.comm, "Zoltan_PHG kway=%d #parts=%d\n", hgp.kway, zz->LB.Num_Global_Parts);
*/

  /* UVC: if it is bisection anyways; no need to create vmap etc; 
     rdivide is going to call Zoltan_PHG_Partition anyways... */
  if (hgp.kway || zz->LB.Num_Global_Parts == 2) {/* call main V cycle routine */
    err = Zoltan_PHG_Partition(zz, &zoltan_hg->HG, zz->LB.Num_Global_Parts,
     hgp.part_sizes, parts, &hgp, 0);
    if (err != ZOLTAN_OK) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error partitioning hypergraph.");
      goto End;
    }
  }
  else {
    int i, p=zz->LB.Num_Global_Parts;
    HGraph *hg = &zoltan_hg->HG;

    /* vmap associates original vertices to sub hypergraphs */
    if (!(hg->vmap = (int*) ZOLTAN_MALLOC(hg->nVtx*sizeof (int))))  {
      err = ZOLTAN_MEMERR;
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      goto End;
    }
    for (i = 0; i < hg->nVtx; ++i)
      hg->vmap[i] = i;

    /* tighten balance tolerance for recursive bisection process */
    hgp.bal_tol = pow (hgp.bal_tol, 1.0 / ceil (log((double)p) / log(2.0)));

    /* partition hypergraph */
    err = Zoltan_PHG_rdivide (0, p-1, parts, zz, hg, &hgp, 0);
    for (i = 0; i < hg->nVtx; ++i)
        if (parts[i]<0 || parts[i]>=p)
            errexit("invalid partvec[%d]=%d", i, parts[i]);

    if (hgp.output_level >= PHG_DEBUG_LIST)     
      uprintf(hg->comm, "FINAL %3d |V|=%6d |E|=%6d #pins=%6d %s/%s/%s p=%d "
       "bal=%.2f cutl=%.2f\n", hg->info, hg->nVtx, hg->nEdge, hg->nPins,
       hgp.redm_str, hgp.coarsepartition_str, hgp.refinement_str, p,
       Zoltan_PHG_Compute_Balance(zz, hg, p, parts),
       Zoltan_PHG_Compute_ConCut(hg->comm, hg, parts, p, &err));
        
    if (err != ZOLTAN_OK)  {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error partitioning hypergraph.");
      goto End;
    }
    ZOLTAN_FREE (&zoltan_hg->HG.vmap);
  }  
  
        
  if (hgp.use_timers > 1) {
    if (timer_retlist < 0) 
      timer_retlist = Zoltan_Timer_Init(zz->ZTime, 1, "Return Lists");
    ZOLTAN_TIMER_START(zz->ZTime, timer_retlist, zz->Communicator);
  }
    
  /* Build Zoltan's Output_Parts, mapped from 2D distribution 
     to input distribution. */

  Zoltan_PHG_Output_Parts(zz, zoltan_hg, parts);

  /* Build Zoltan's return arguments. */
  Zoltan_PHG_Return_Lists(zz, zoltan_hg, num_exp, exp_gids,
   exp_lids, exp_procs, exp_to_part);
    
  if (hgp.use_timers > 1) 
    ZOLTAN_TIMER_STOP(zz->ZTime, timer_retlist, zz->Communicator);

End:
  if (err == ZOLTAN_MEMERR)
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Memory error.")
  else if (err != ZOLTAN_OK)
    ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Error partitioning hypergraph.")
    
  if (hgp.use_timers)
    ZOLTAN_TIMER_STOP(zz->ZTime, timer_all, zz->Communicator);

  /* KDDKDD The following code prints a final quality result even when
   * KDDKDD phg_output_level is zero.  It is useful for our tests and
   * KDDKDD data collection, but it should NOT be included in the released
   * KDDKDD code.  */
  if ((err == ZOLTAN_OK) && hgp.final_output) {
    HGraph *hg = &zoltan_hg->HG;
    static int nRuns=0;
    static double balsum = 0.0, cutlsum = 0.0, cutnsum = 0.0;
    static double balmax = 0.0, cutlmax = 0.0, cutnmax = 0.0;
    static double balmin = 1e100, cutlmin = 1e100, cutnmin = 1e100;
    double bal = Zoltan_PHG_Compute_Balance(zz, hg, zz->LB.Num_Global_Parts,
                                            parts);
    double cutl;   /* Connnectivity cuts:  sum_over_edges((npart-1)*ewgt) */
    double cutn;   /* Net cuts:  sum_over_edges((nparts>1)*ewgt) */
    double remcutl;   /* Connnectivity cuts of removed edges */
    double remcutn;   /* Net cuts of removed edges */
    int gnremove;

    cutl= Zoltan_PHG_Compute_ConCut(hg->comm, hg, parts,
                                    zz->LB.Num_Global_Parts, &err);
    cutn = Zoltan_PHG_Compute_NetCut(hg->comm, hg, parts,
                                     zz->LB.Num_Global_Parts);

    if (!err) {
     
      /* Add in cut contributions from removed edges */
      MPI_Allreduce(&(zoltan_hg->nRemove), &gnremove, 1, MPI_INT, MPI_SUM,
                    zz->Communicator);
      if (gnremove) {
        err = Zoltan_PHG_Removed_Cuts(zz, zoltan_hg, &remcutl, &remcutn);
        cutl += remcutl;
        cutn += remcutn;
      }
  
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
      }
    }
  }
  /* KDDKDD  End of printing section. */
  
  ZOLTAN_FREE(&parts);
  Zoltan_HG_Free_Structure(zz);

  if (hgp.use_timers && zz->Proc == 0)
    Zoltan_Timer_PrintAll(zz->ZTime, zz->Proc, stdout);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return err;
}

/*****************************************************************************/

static int Zoltan_PHG_Initialize_Params(
  ZZ *zz,   /* the Zoltan data structure */
  float *part_sizes,
  PHGPartParams *hgp
)
{
  int err;
  

  memset(hgp, 0, sizeof(*hgp)); /* in the future if we forget to initialize
                                   another param at least it will be 0 */
  
  Zoltan_Bind_Param(PHG_params, "PHG_OUTPUT_LEVEL", &hgp->output_level);
  Zoltan_Bind_Param(PHG_params, "PHG_FINAL_OUTPUT", &hgp->final_output); 
  Zoltan_Bind_Param(PHG_params, "FINAL_OUTPUT", &hgp->final_output); 
  Zoltan_Bind_Param(PHG_params, "CHECK_GRAPH", &hgp->check_graph);   
  Zoltan_Bind_Param(PHG_params, "PHG_NPROC_X", &hgp->nProc_x_req);
  Zoltan_Bind_Param(PHG_params, "PHG_NPROC_Y", &hgp->nProc_y_req);
  Zoltan_Bind_Param(PHG_params, "PHG_PROC_SPLIT", &hgp->proc_split);  
  Zoltan_Bind_Param(PHG_params, "PHG_REDUCTION_LIMIT", &hgp->redl);
  Zoltan_Bind_Param(PHG_params, "PHG_REDUCTION_METHOD", hgp->redm_str);
  Zoltan_Bind_Param(PHG_params, "PHG_REDUCTION_LOCAL_IMPROVEMENT", 
                                 hgp->redmo_str);  
  Zoltan_Bind_Param(PHG_params, "PHG_VERTEX_VISIT_ORDER", &hgp->visit_order);
  Zoltan_Bind_Param(PHG_params, "PHG_EDGE_SCALING", &hgp->edge_scaling);
  Zoltan_Bind_Param(PHG_params, "PHG_VERTEX_SCALING", &hgp->vtx_scaling);
  Zoltan_Bind_Param(PHG_params, "PHG_REFINEMENT", hgp->refinement_str);
  Zoltan_Bind_Param(PHG_params, "PHG_DIRECT_KWAY", &hgp->kway);
  Zoltan_Bind_Param(PHG_params, "PHG_FM_LOOP_LIMIT", &hgp->fm_loop_limit);
  Zoltan_Bind_Param(PHG_params, "PHG_FM_MAX_NEG_MOVE", &hgp->fm_max_neg_move);  
  Zoltan_Bind_Param(PHG_params, "PHG_COARSE_PARTITIONING", 
                                 hgp->coarsepartition_str);

  Zoltan_Bind_Param(PHG_params, "PHG_COARSE_ITERATIONS",
                                 (void*) &hgp->num_coarse_iter);  
  Zoltan_Bind_Param(PHG_params, "PHG_USE_TIMERS",
                                 (void*) &hgp->use_timers);  
  Zoltan_Bind_Param(PHG_params, "USE_TIMERS",
                                 (void*) &hgp->use_timers);  
  Zoltan_Bind_Param(PHG_params, "EDGE_SIZE_THRESHOLD",
                                 (void*) &hgp->EdgeSizeThreshold);  

  /* Set default values */
  strncpy(hgp->redm_str,            "ipm",   MAX_PARAM_STRING_LEN);
  strncpy(hgp->redmo_str,           "no",    MAX_PARAM_STRING_LEN);  
  strncpy(hgp->coarsepartition_str, "gr0",   MAX_PARAM_STRING_LEN);
  strncpy(hgp->refinement_str,      "fm2",   MAX_PARAM_STRING_LEN);

  hgp->use_timers = 0;
  hgp->proc_split = 1;
  hgp->LocalCoarsePartition = 0;
  hgp->locmatching = NULL;
  hgp->edge_scaling = 0;
  hgp->vtx_scaling = 0;
  hgp->vtx_scal = NULL;  /* Array for storing vertex degree scale vector. 
                            Should perhaps go in hg structure, not the
                            param struct? */
  hgp->visit_order = 1;  /* Random */
  hgp->check_graph = 0;
  hgp->bal_tol = zz->LB.Imbalance_Tol[0];
  hgp->redl = MAX(2*zz->LB.Num_Global_Parts, 100);
  hgp->output_level = PHG_DEBUG_LIST;
  hgp->final_output = 0;
  hgp->nProc_x_req = -1;
  hgp->nProc_y_req = -1;
  hgp->kway = 0;
  hgp->fm_loop_limit = 99;
  hgp->fm_max_neg_move = 250;  
  hgp->num_coarse_iter = 10;
  hgp->part_sizes = part_sizes;
  hgp->EdgeSizeThreshold = 0.5;  

  /* Get application values of parameters. */
  Zoltan_Assign_Param_Vals(zz->Params, PHG_params, zz->Debug_Level, zz->Proc,
                           zz->Debug_Proc);

  err = Zoltan_PHG_Set_2D_Proc_Distrib(zz, zz->Communicator, zz->Proc, 
                                       zz->Num_Proc, hgp->nProc_x_req, 
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

  return Zoltan_Check_Param (name, val, PHG_params, &result, &index);
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
ZOLTAN_ID_PTR gids    = zhg->GIDs;
ZOLTAN_ID_PTR lids    = zhg->LIDs;
int *outparts         = zhg->Output_Parts; 

  if (zz->LB.Return_Lists == ZOLTAN_LB_NO_LISTS) 
    goto End;

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
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient Memory");
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

/****************************************************************************/

void Zoltan_PHG_HGraph_Print(
  ZZ *zz,          /* the Zoltan data structure */
  ZHG *zoltan_hg,
  HGraph *hg,
  Partition parts, 
  FILE *fp
)
{
/* Printing routine. Can be used to print a Zoltan_HGraph or just an HGraph.
 * Set zoltan_hg to NULL if want to print only an HGraph.
 * Lots of output; synchronized across processors, so is a bottleneck.
 */
  int i;
  int num_gid = zz->Num_GID;
  int num_lid = zz->Num_LID;
  char *yo = "Zoltan_PHG_HGraph_Print";

  if (zoltan_hg != NULL  &&  hg != &zoltan_hg->HG) {
    ZOLTAN_PRINT_WARN(zz->Proc, yo, "Input hg != Zoltan HG");
    return;
  }

  Zoltan_Print_Sync_Start (zz->Communicator, 1);

  /* Print Vertex Info */
  fprintf (fp, "%s Proc %d\n", yo, zz->Proc);
  fprintf (fp, "Vertices (GID, LID, index)\n");
  for (i = 0; i < zoltan_hg->nObj; i++) {
    fprintf(fp, "(");
    ZOLTAN_PRINT_GID(zz, &zoltan_hg->GIDs[i * num_gid]);
    fprintf(fp, ", ");
    ZOLTAN_PRINT_LID(zz, &zoltan_hg->LIDs [i * num_lid]);
    fprintf(fp, ", %d)\n", i);
  }
  Zoltan_HG_Print(zz, hg, parts, fp, "Build");  
  Zoltan_Print_Sync_End(zz->Communicator, 1);
}

/*****************************************************************************/

int Zoltan_PHG_Set_2D_Proc_Distrib(
    ZZ *zz,                /* Input:  ZZ struct; for debuging   */
    MPI_Comm Communicator, /* Input:  The MPI Communicator      */
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
                       "Values for PHG_NPROC_X and PHG_NPROC_Y "
                       "do not evenly divide the "
                       "total number of processors.");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  comm->zz = zz;
  comm->myProc_x = proc % comm->nProc_x;
  comm->myProc_y = proc / comm->nProc_x;
  comm->Communicator = Communicator;
  comm->myProc = proc;
  comm->nProc = nProc;

  if (Communicator==MPI_COMM_NULL) {
      comm->col_comm = comm->row_comm = MPI_COMM_NULL;
  } else {
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

