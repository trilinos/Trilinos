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
  {"PHG_REDUCTION_LIMIT",        NULL, "INT",    0},
  {"PHG_EDGE_WEIGHT_SCALING",    NULL, "INT",    0},
  {"PHG_REDUCTION_METHOD",       NULL, "STRING", 0},
  {"PHG_COARSE_PARTITIONING",    NULL, "STRING", 0},
  {"PHG_REFINEMENT",             NULL, "STRING", 0},
  {"PHG_NPROC_X",                NULL, "INT",    0},
  {"PHG_NPROC_Y",                NULL, "INT",    0},
  {"PCHECK_GRAPH",               NULL, "INT",    0},
  {"PHG_OUTPUT_LEVEL",           NULL, "INT",    0},
  {"PHG_DIRECT_KWAY",            NULL, "INT",    0},
  {"PHG_FM_LOOP_LIMIT",          NULL, "INT",    0},
  {"PHG_FM_MAX_NEG_MOVE",        NULL, "INT",    0},    
  {NULL,                         NULL,  NULL,    0} 
};

/* prototypes for static functions: */
static int set_proc_distrib (MPI_Comm, int, int, PHGComm*);
static int Zoltan_PHG_Initialize_Params (ZZ*, PHGPartParams*);
static int Zoltan_PHG_Return_Lists (ZZ*, ZPHG*, Partition, int*,
 ZOLTAN_ID_PTR*, ZOLTAN_ID_PTR*, int**, int**);

 
 
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
  int **exp_to_part          /* list of partitions to which exported objs
                                are assigned. */
)
{
    ZPHG *zoltan_hg = NULL;
    int nVtx;                        /* Temporary variable for base graph. */
    PHGPartParams hgp;               /* Hypergraph parameters. */
    Partition output_parts = NULL;   /* Output partition from HG partitioner. */
    int err = ZOLTAN_OK;
    char *yo = "Zoltan_PHG";
    
    ZOLTAN_TRACE_ENTER(zz, yo);
    
    /* Initialization of return arguments. */
    *num_imp   = *num_exp   = -1;
    *imp_gids  = *exp_gids  = NULL;
    *imp_lids  = *exp_lids  = NULL;
    *imp_procs = *exp_procs = NULL;
    
    /* Initialize HG parameters. */
    err = Zoltan_PHG_Initialize_Params (zz, &hgp);
    if (err != ZOLTAN_OK)
        goto End;
    
    /* build initial Zoltan hypergraph from callback functions. */
    err = Zoltan_PHG_Build_Hypergraph (zz, &zoltan_hg, &hgp);
    if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error building hypergraph.");
        goto End;
    }
   
    zz->LB.Data_Structure = zoltan_hg;
    nVtx = zoltan_hg->PHG.nVtx;
    zoltan_hg->PHG.redl = hgp.redl;     /* redl needs to be dynamic */
    /* RTHRTH -- redl may need to be scaled by number of procs */
 
    /* allocate output partition memory */
    output_parts = (Partition) ZOLTAN_MALLOC (nVtx * sizeof(int));
    if (nVtx && output_parts == NULL) {
        err = ZOLTAN_MEMERR;
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        goto End;
    }

    /*
      uprintf(zoltan_hg->PHG.comm, "Zoltan_PHG kway=%d #parts=%d\n", hgp.kway, zz->LB.Num_Global_Parts);
    */

    /* UVC: if it is bisection anyways; no need to create vmap etc; rdrive is going to call
     Zoltan_PHG_HPart_Lib anyways... */
    if (hgp.kway || zz->LB.Num_Global_Parts==2) { 
        /* call main V cycle routine */
        err = Zoltan_PHG_HPart_Lib(zz, &zoltan_hg->PHG, zz->LB.Num_Global_Parts,
                                   output_parts, &hgp, 0);
        if (err != ZOLTAN_OK) {
            ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error partitioning hypergraph.");
            goto End;
        }
    }
    else {
        int i, p=zz->LB.Num_Global_Parts;
        PHGraph *hg = &zoltan_hg->PHG;

        /* vmap associates original vertices to sub hypergraphs */
        if (!(hg->vmap = (int*) ZOLTAN_MALLOC(hg->nVtx*sizeof (int))))  {
            err = ZOLTAN_MEMERR;
            ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
            goto End;
        }
        for (i = 0; i < hg->nVtx; i++)
            hg->vmap[i] = i;

#if 0
        /* tighten balance tolerance for recursive bisection process */
        hgp.bal_tol = pow (hgp.bal_tol,
                           1.0 / ceil (log((double)p) / log(2.0)));
#endif
        /* partition hypergraph */
        err = Zoltan_PHG_rdivide(1, p, output_parts, zz, 
                                 hg, &hgp, 0);

        if (hgp.output_level >= PHG_DEBUG_LIST)     
            uprintf(hg->comm, "FINAL %3d |V|=%6d |E|=%6d |Z|=%6d %s/%s/%s p=%d bal=%.2f cutl=%.2f\n",
                    hg->info, hg->nVtx, hg->nEdge, hg->nNonZero, hgp.redm_str,
                    hgp.coarsepartition_str, hgp.refinement_str, p,
                    Zoltan_PHG_HPart_balance(zz, hg, p, output_parts),
                    Zoltan_PHG_hcut_size_links(hg->comm, hg, output_parts, p));
        
        if (err != ZOLTAN_OK)  {
            ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error partitioning hypergraph.");
            goto End;
        }
        ZOLTAN_FREE (&zoltan_hg->PHG.vmap);
    }  
  
    /* Build Zoltan's return arguments. */
    Zoltan_PHG_Return_Lists(zz, zoltan_hg, output_parts, num_exp, exp_gids,
                            exp_lids, exp_procs, exp_to_part);
    
End:
    if (err == ZOLTAN_MEMERR) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
    } else if (err != ZOLTAN_OK) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error partitioning hypergraph.");
    }
    ZOLTAN_FREE((void**) &output_parts);
    Zoltan_PHG_Free_Structure(zz);
    ZOLTAN_TRACE_EXIT(zz, yo);
    return err;
}

/*****************************************************************************/

void Zoltan_PHG_Free_Structure(ZZ *zz)
{
  /* frees all data associated with LB.Data_Structure for hypergraphs */
  ZPHG *zoltan_hg = (ZPHG*) zz->LB.Data_Structure;

  if (zoltan_hg != NULL) {
    Zoltan_Multifree(__FILE__, __LINE__, 3, &zoltan_hg->Global_IDs,
     &zoltan_hg->Local_IDs, &zoltan_hg->Parts);
    Zoltan_PHG_HGraph_Free(&zoltan_hg->PHG);
    ZOLTAN_FREE ((void**) &zz->LB.Data_Structure);
  }
}

/*****************************************************************************/

static int Zoltan_PHG_Initialize_Params(
  ZZ *zz,   /* the Zoltan data structure */
  PHGPartParams *hgp
)
{
  int ierr;
  
  Zoltan_Bind_Param(PHG_params, "PHG_OUTPUT_LEVEL",
                    (void*) &hgp->output_level);
  Zoltan_Bind_Param(PHG_params, "PHG_NPROC_X",
                    (void*) &(hgp->comm.nProc_x));
  Zoltan_Bind_Param(PHG_params, "PHG_NPROC_Y",
                    (void*) &(hgp->comm.nProc_y));
  Zoltan_Bind_Param(PHG_params, "PHG_REDUCTION_LIMIT",
                    (void*) &hgp->redl);
  Zoltan_Bind_Param(PHG_params, "PHG_REDUCTION_METHOD",
                    (void*) hgp->redm_str);
  Zoltan_Bind_Param(PHG_params, "PHG_EDGE_WEIGHT_SCALING",
                    (void*) &hgp->ews);
  Zoltan_Bind_Param(PHG_params, "PCHECK_GRAPH",
                    (void*) &hgp->check_graph);   
  Zoltan_Bind_Param(PHG_params, "PHG_REFINEMENT",
                    (void*) hgp->refinement_str);
  Zoltan_Bind_Param(PHG_params, "PHG_DIRECT_KWAY",
                    (void*) &hgp->kway);
  Zoltan_Bind_Param(PHG_params, "PHG_FM_LOOP_LIMIT",
                    (void*) &hgp->fm_loop_limit);
  Zoltan_Bind_Param(PHG_params, "PHG_FM_MAX_NEG_MOVE",
                    (void*) &hgp->fm_max_neg_move);  
  Zoltan_Bind_Param(PHG_params, "PHG_COARSE_PARTITIONING", 
                    (void*) hgp->coarsepartition_str);

  /* Set default values */
  strncpy(hgp->redm_str,            "no",  MAX_PARAM_STRING_LEN);
  strncpy(hgp->coarsepartition_str, "gr0", MAX_PARAM_STRING_LEN);
  strncpy(hgp->refinement_str,      "no",  MAX_PARAM_STRING_LEN);
  
  hgp->ews = 1;
  hgp->check_graph = 1;
  hgp->bal_tol = zz->LB.Imbalance_Tol[0];
  hgp->redl = zz->LB.Num_Global_Parts;
  hgp->output_level = PHG_DEBUG_LIST;
  hgp->comm.nProc_x = -1;
  hgp->comm.nProc_y = -1;
  hgp->kway = 0;
  hgp->fm_loop_limit = 99;
  hgp->fm_max_neg_move = 250;  

  /* Get application values of parameters. */
  Zoltan_Assign_Param_Vals(zz->Params, PHG_params, zz->Debug_Level, zz->Proc,
                           zz->Debug_Proc);

  ierr = set_proc_distrib(zz->Communicator, zz->Proc, zz->Num_Proc, 
                          &hgp->comm);
  if (ierr != ZOLTAN_OK) 
      goto End;

  
  /* Convert strings to function pointers. */
  ierr = Zoltan_PHG_Set_Part_Options (zz, hgp);
  
End:
  return ierr;
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

static int Zoltan_PHG_Return_Lists(
  ZZ *zz,
  ZPHG *zoltan_hg,
  Partition output_parts,
  int *num_exp,
  ZOLTAN_ID_PTR *exp_gids,
  ZOLTAN_ID_PTR *exp_lids,
  int **exp_procs,
  int **exp_to_part
)
{
/* Routine to build export lists of ZOLTAN_LB_FN. */
int i, j;
int msg_tag = 31000;
int ierr = ZOLTAN_OK;
int eproc;
int num_gid_entries   = zz->Num_GID;
int num_lid_entries   = zz->Num_LID;
int nObj              = zoltan_hg->nObj;
Partition input_parts = zoltan_hg->Parts;
ZOLTAN_ID_PTR gids    = zoltan_hg->Global_IDs;
ZOLTAN_ID_PTR lids    = zoltan_hg->Local_IDs;
int *outparts = NULL;     /* Pointers to output partitions. */
char *yo = "Zoltan_PHG_Return_Lists";

  if (zz->LB.Return_Lists == ZOLTAN_LB_NO_LISTS) 
    goto End;

  if (zoltan_hg->VtxPlan != NULL) {
    /* Get the partition information from the 2D decomposition back to the
     * original owning processor for each GID.
     */
    outparts = (int *) ZOLTAN_MALLOC(nObj * sizeof(int));
    ierr = Zoltan_Comm_Do_Reverse(zoltan_hg->VtxPlan, msg_tag, 
                                  (char *) output_parts, 
                                  sizeof(int), NULL, (char *) outparts);
    if (ierr) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error from Zoltan_Comm_Do_Reverse");
      goto End;
    }
  }
  else
    outparts = output_parts;

  /* Count number of objects with new partitions or new processors. */
  *num_exp = 0;
  for (i = 0; i < nObj; i++) {
    eproc = Zoltan_LB_Part_To_Proc(zz, outparts[i], 
                                   &gids[i*num_gid_entries]);
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
      eproc = Zoltan_LB_Part_To_Proc(zz, outparts[i], 
       &gids[i*num_gid_entries]);
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

  if (zoltan_hg->VtxPlan != NULL) {
    ZOLTAN_FREE(&outparts);
    Zoltan_Comm_Destroy(&(zoltan_hg->VtxPlan));
  }

  return ierr;
}

/****************************************************************************/

void Zoltan_PHG_HGraph_Print(
  ZZ *zz,          /* the Zoltan data structure */
  ZPHG *zoltan_hg,
  PHGraph *hg,
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

  if (zoltan_hg != NULL  &&  hg != &zoltan_hg->PHG) {
    ZOLTAN_PRINT_WARN(zz->Proc, yo, "Input hg != Zoltan HG");
    return;
  }

  Zoltan_Print_Sync_Start (zz->Communicator, 1);

  /* Print Vertex Info */
  fprintf (fp, "%s Proc %d\n", yo, zz->Proc);
  fprintf (fp, "Vertices (GID, LID, index)\n");
  for (i = 0; i < hg->nVtx; i++) {
    fprintf(fp, "(");
    ZOLTAN_PRINT_GID(zz, &zoltan_hg->Global_IDs[i * num_gid]);
    fprintf(fp, ", ");
    ZOLTAN_PRINT_LID(zz, &zoltan_hg->Local_IDs [i * num_lid]);
    fprintf(fp, ", %d)\n", i);
  }

  Zoltan_PHG_Print(zz, hg, fp);
  Zoltan_Print_Sync_End(zz->Communicator, 1);
}

/*****************************************************************************/

static int set_proc_distrib(
  MPI_Comm Communicator, /* Input:  The MPI Communicator      */
  int proc,              /* Input:  Rank of current processor */
  int nProc,             /* Input:  Total # of processors     */
  PHGComm *comm          /* Input/Ouput: for nProc_x and nProc_y members;
                            Output: for the rest:
  int *nProc_x,      Input/Output:  # processors in x-direction of 2D 
                        data distrib; if -1 on input, compute. 
  int *nProc_y,      Input/Output:  # processors in y-direction of 2D 
                        data distrib; if -1 on input, compute. 
  int *myProc_x,     Output:  x block of proc in [0,nProc_x-1]. 
  int *myProc_y      Output:  y block of proc in [0,nProc_y-1]. */
)
{
/* Computes the processor distribution for the 2D data distrib.
 * Sets nProc_x, nProc_y.
 * Constraint:  nProc_x * nProc_y == nProc. 
 * For 2D data distrib, default should approximate sqrt(nProc).
 * If nProc_x and nProc_y both equal -1 on input, compute default.
 * Otherwise, compute valid values and/or return error.
 */
char *yo = "set_proc_distrib";
int tmp;
int ierr = ZOLTAN_OK;
    
  if (comm->nProc_x == -1 && comm->nProc_y == -1) {
    /* Compute default */
    tmp = (int) sqrt((double)nProc+0.1);
    while (nProc % tmp) tmp--;
    comm->nProc_y = tmp;
    comm->nProc_x = nProc / tmp;
  } else if (comm->nProc_x == -1) {
    /* nProc_y set by user parameter */
    comm->nProc_x = nProc / comm->nProc_y;
  } else if (comm->nProc_y == -1) {
    /* nProc_x set by user parameter */
    comm->nProc_y = nProc / comm->nProc_x;
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
    
  comm->myProc_x = proc % comm->nProc_x;
  comm->myProc_y = proc / comm->nProc_x;
  comm->Communicator = Communicator;
  comm->Proc = proc;
  comm->Num_Proc = nProc;
    
  if ((MPI_Comm_split(Communicator, comm->myProc_x, comm->myProc_y, 
                      &comm->col_comm) != MPI_SUCCESS)
   || (MPI_Comm_split(Communicator, comm->myProc_y, comm->myProc_x, 
                      &comm->row_comm) != MPI_SUCCESS)) {
    ZOLTAN_PRINT_ERROR(proc, yo, "MPI_Comm_Split failed");
    return ZOLTAN_FATAL;
  }
/*  printf("(%d, %d) of [%d, %d] -> After Comm_split col_comm=%d  row_comm=%d\n", hgp->myProc_x, hgp->myProc_y, hgp->nProc_x, hgp->nProc_y, (int)hgp->col_comm, (int)hgp->row_comm);  */
  
    
End:

  return ierr;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

