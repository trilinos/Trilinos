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

#include "hg.h"
#include "params_const.h"
#include "all_allo_const.h"

/*
 *  Main routine for Zoltan interface to hypergraph partitioning.
 *  Builds input data structures, set parameters, etc.
 */

/*****************************************************************************/
/*  Parameters structure for HG method.  */

static PARAM_VARS HG_params[] = {
 /* Add parameters here. */
 { "HG_REDUCTION_LIMIT", NULL, "INT" },
 { "HG_EDGE_WEIGHT_SCALING", NULL, "INT" },
 { "HG_REDUCTION_METHOD", NULL, "STRING" },
 { "HG_GLOBAL_PARTITIONING", NULL, "STRING" },
 { "HG_LOCAL_REFINEMENT", NULL, "STRING" },
 { "HG_REDUCTION_LOCAL_IMPROVEMENT", NULL, "STRING" },
 { "CHECK_GRAPH", NULL, "INT" },
 { NULL, NULL, NULL } };

static int Zoltan_HG_Initialize_Params(ZZ *, HGPartParams *);
static int Zoltan_HG_Return_Lists(ZZ *, struct Zoltan_HGraph *, Partition,
  int *, ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *, int **, int **);


/*****************************************************************************/

/* Main partitioning routine:  Type = ZOLTAN_LB_FN.  */

int Zoltan_HG(
  ZZ *zz,                    /* The Zoltan structure  */
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
char *yo = "Zoltan_HG";
struct Zoltan_HGraph *zoltan_hg = NULL;
int ierr = ZOLTAN_OK;
int nVtx;                       /* Temporary variable for base graph. */

HGPartParams hgp;               /* Hypergraph parameters. */

Partition output_parts = NULL;  /* Output partition from HG partitioner. */

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Initialization of return arguments. */
  *num_imp = *num_exp = -1;
  *imp_gids = *exp_gids = NULL;
  *imp_lids = *exp_lids = NULL;
  *imp_procs = *exp_procs = NULL;

  /* Initialize HG parameters. */
  ierr = Zoltan_HG_Initialize_Params(zz, &hgp);
  if (ierr != ZOLTAN_OK) 
    goto End;

  /* build initial Zoltan hypergraph. */
  ierr = Zoltan_HG_Build_Hypergraph(zz, &zoltan_hg, hgp.check_graph);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error building hypergraph.");
    goto End;
  }

  zz->LB.Data_Structure = zoltan_hg;
  nVtx = zoltan_hg->HG.nVtx;

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
    Zoltan_HG_HGraph_Print(zz, zoltan_hg, &(zoltan_hg->HG));

  /* Call partitioning routines. */
  output_parts = (Partition) ZOLTAN_MALLOC(nVtx * sizeof(int));
  if (output_parts == NULL) {
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  ierr = Zoltan_HG_HPart_Lib(zz, &(zoltan_hg->HG), zz->LB.Num_Global_Parts,
                             output_parts, &hgp);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error partitioning hypergraph.");
    goto End;
  }

  /* Build Zoltan's return arguments. */
  Zoltan_HG_Return_Lists(zz, zoltan_hg, output_parts, 
                         num_exp, exp_gids, exp_lids, exp_procs, exp_to_part);



End:
  if (ierr == ZOLTAN_MEMERR)
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
  ZOLTAN_FREE((void **) &output_parts);
  Zoltan_HG_Free_Structure(zz);
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}


/*****************************************************************************/

void Zoltan_HG_Free_Structure(ZZ *zz)
{
struct Zoltan_HGraph *zoltan_hg=(struct Zoltan_HGraph*)(zz->LB.Data_Structure);

  if (zoltan_hg != NULL) {
    ZOLTAN_FREE((void **) &(zoltan_hg->Global_IDs));
    ZOLTAN_FREE((void **) &(zoltan_hg->Local_IDs));
    ZOLTAN_FREE((void **) &(zoltan_hg->Parts));
    Zoltan_HG_HGraph_Free(&(zoltan_hg->HG));
    ZOLTAN_FREE((void **) &(zz->LB.Data_Structure));
  }
}


/*****************************************************************************/
static int Zoltan_HG_Initialize_Params(
  ZZ *zz,
  HGPartParams *hgp
)
{
int ierr = ZOLTAN_OK;

  /* Bind storage with parameter. */
  Zoltan_Bind_Param(HG_params, "HG_REDUCTION_LIMIT",
                                (void *) &(hgp->redl));
  Zoltan_Bind_Param(HG_params, "HG_REDUCTION_METHOD", 
                                (void *) hgp->redm_str);
  Zoltan_Bind_Param(HG_params, "HG_REDUCTION_LOCAL_IMPROVEMENT",
                                (void *) hgp->redmo_str);
  Zoltan_Bind_Param(HG_params, "HG_EDGE_WEIGHT_SCALING",
                                (void *) &(hgp->ews));
  Zoltan_Bind_Param(HG_params, "HG_GLOBAL_PARTITIONING",
                                (void *) hgp->global_str);
  Zoltan_Bind_Param(HG_params, "HG_LOCAL_REFINEMENT", 
                                (void *) hgp->local_str);
  Zoltan_Bind_Param(HG_params, "CHECK_GRAPH",
                                (void *) &(hgp->check_graph));

  /* Set default values */
  hgp->redl = zz->LB.Num_Global_Parts;
  strcpy(hgp->redm_str, "grg");
  strcpy(hgp->redmo_str, "aug3");
  hgp->ews = 1;
  strcpy(hgp->global_str, "gr0");
  strcpy(hgp->local_str, "fm");
  hgp->check_graph = 1;
  hgp->bal_tol = zz->LB.Imbalance_Tol;

  /* Get application values of parameters. */
  Zoltan_Assign_Param_Vals(zz->Params, HG_params, zz->Debug_Level, zz->Proc,
                           zz->Debug_Proc);

  /* Convert strings to function pointers. */

  ierr = Zoltan_HG_Set_Part_Options(zz, hgp);

  return ierr;
}


/*****************************************************************************/

int Zoltan_HG_Set_Param(
char *name,                     /* name of variable */
char *val)                      /* value of variable */
{
    int status;
    PARAM_UTYPE result;         /* value returned from Check_Param */
    int index;                  /* index returned from Check_Param */

    status = Zoltan_Check_Param(name, val, HG_params, &result, &index);

    return(status);
}


/*****************************************************************************/

static int Zoltan_HG_Return_Lists(
  ZZ *zz,
  struct Zoltan_HGraph *zoltan_hg,
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
int num_gid_entries = zz->Num_GID;
int num_lid_entries = zz->Num_LID;
int nVtx = zoltan_hg->HG.nVtx;
int eproc;
Partition input_parts = zoltan_hg->Parts;
ZOLTAN_ID_PTR gids = zoltan_hg->Global_IDs;
ZOLTAN_ID_PTR lids = zoltan_hg->Local_IDs;


  if (zz->LB.Return_Lists) {

    /* Count number of objects with new partitions or new processors. */
    *num_exp = 0;
    for (i = 0; i < nVtx; i++)
      if (output_parts[i] != input_parts[i]
       || zz->Proc != Zoltan_LB_Part_To_Proc(zz, output_parts[i]))
         (*num_exp)++;

    /* Allocate memory for return lists. */
    if (*num_exp > 0) {
      if (!Zoltan_Special_Malloc(zz, (void **)exp_gids, *num_exp,
                                 ZOLTAN_SPECIAL_MALLOC_GID) ||
          !Zoltan_Special_Malloc(zz, (void **)exp_lids, *num_exp,
                                 ZOLTAN_SPECIAL_MALLOC_LID) ||
          !Zoltan_Special_Malloc(zz, (void **)exp_procs, *num_exp,
                                 ZOLTAN_SPECIAL_MALLOC_INT) ||
          !Zoltan_Special_Malloc(zz, (void **)exp_to_part, *num_exp,
                                 ZOLTAN_SPECIAL_MALLOC_INT)) {
        Zoltan_Special_Free(zz,(void **)exp_gids, ZOLTAN_SPECIAL_MALLOC_GID);
        Zoltan_Special_Free(zz,(void **)exp_lids, ZOLTAN_SPECIAL_MALLOC_LID);
        Zoltan_Special_Free(zz,(void **)exp_procs, ZOLTAN_SPECIAL_MALLOC_INT);
        Zoltan_Special_Free(zz,(void **)exp_to_part, ZOLTAN_SPECIAL_MALLOC_INT);
        return ZOLTAN_MEMERR;
      }

      for (j = 0, i = 0; i < nVtx; i++) {
        eproc = Zoltan_LB_Part_To_Proc(zz, output_parts[i]);
        if (output_parts[i] != input_parts[i] ||
            eproc != zz->Proc) {
          ZOLTAN_SET_GID(zz, &((*exp_gids)[j*num_gid_entries]),
                             &(gids[i*num_gid_entries]));
          if (num_lid_entries > 0)
            ZOLTAN_SET_LID(zz, &((*exp_lids)[j*num_lid_entries]),
                               &(lids[i*num_lid_entries]));
          (*exp_procs)[j] = eproc;
          (*exp_to_part)[j] = output_parts[i];
          j++;
        }
      }
    }
  }
  return ZOLTAN_OK;
}


/****************************************************************************/
void Zoltan_HG_HGraph_Print(
  ZZ *zz,
  struct Zoltan_HGraph *zoltan_hg,
  HGraph *hg
)
{
/* Printing routine.
 * Can be used to print a Zoltan_HGraph or just an HGraph.
 * Set zoltan_hg to NULL if want to print only an HGraph.
 * Lots of output; synchronized across processors, so is a bottleneck.
 */
char *yo = "Zoltan_HG_HGraph_Print";
int i;
int num_gid = zz->Num_GID;
int num_lid = zz->Num_LID;

  if ((zoltan_hg != NULL) && (hg != &(zoltan_hg->HG))) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Input hg != Zoltan HG");
    return;
  }

  Zoltan_Print_Sync_Start(zz->Communicator, 1);

  /* Print Vertex Info */
  printf("%s Proc %d\n", yo, zz->Proc);
  printf("Vertices (GID, LID, index)\n");

  for (i = 0; i < hg->nVtx; i++) {
    printf("(");
    ZOLTAN_PRINT_GID(zz,  &(zoltan_hg->Global_IDs[i*num_gid]));
    printf(", ");
    ZOLTAN_PRINT_LID(zz,  &(zoltan_hg->Local_IDs[i*num_lid]));
    printf(", %d)\n", i); 
    
  }

  Zoltan_HG_Print(zz, hg);

  Zoltan_Print_Sync_End(zz->Communicator, 1);
}


/*****************************************************************************/
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
