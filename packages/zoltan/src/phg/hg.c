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
#include "hg.h"
#include "params_const.h"
#include "all_allo_const.h"

/* static double hcut_size_links (ZZ *zz, HGraph *hg, int p, Partition part);
   static double hcut_size_total (HGraph *hg, Partition part);
*/

/*
 *  Main routine for Zoltan interface to hypergraph partitioning.
 *  Builds input data structures, set parameters, etc.
 */



/*****************************************************************************/
/*  Parameters structure for HG method.  */
static PARAM_VARS HG_params[] = {
  /* Add parameters here. */
  {"HG_REDUCTION_LIMIT",             NULL, "INT",    0},
  {"HG_EDGE_WEIGHT_SCALING",         NULL, "INT",    0},
  {"HG_REDUCTION_METHOD",            NULL, "STRING", 0},
  {"HG_GLOBAL_PARTITIONING",         NULL, "STRING", 0},
  {"HG_LOCAL_REFINEMENT",            NULL, "STRING", 0},
  {"HG_REDUCTION_LOCAL_IMPROVEMENT", NULL, "STRING", 0},
  {"CHECK_GRAPH",                    NULL, "INT",    0},
  {"HG_OUTPUT_LEVEL",                NULL, "INT",    0},
  {NULL,                             NULL,  NULL,    0} 
};

/* prototypes for static functions: */
static int Zoltan_HG_Initialize_Params (ZZ*, HGPartParams*);
static int Zoltan_HG_Return_Lists (ZZ*, ZHG*, Partition, int*,
 ZOLTAN_ID_PTR*, ZOLTAN_ID_PTR*, int**, int**);

 
 
/******************************************************************************/
/* Main routine for Zoltan interface to hypergraph partitioning. Builds input */
/* data structures, set parameters, calls HG partitioner, builds return lists.*/
/* Type = ZOLTAN_LB_FN.                                                       */

int Zoltan_HG(
  ZZ *zz,                    /* The Zoltan structure  */
  float *part_sizes,   /* Input:  Array of size zz->Num_Global_Parts containing
                       /* the percentage of work assigned to each partition. */
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
  ZHG *zoltan_hg = NULL;
  int nVtx;                       /* Temporary variable for base graph. */
  HGPartParams hgp;               /* Hypergraph parameters. */
  Partition output_parts = NULL;  /* Output partition from HG partitioner. */
  int err = ZOLTAN_OK;
  int i;
  char *yo = "Zoltan_HG";

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Initialization of return arguments. */
  *num_imp   = *num_exp   = -1;
  *imp_gids  = *exp_gids  = NULL;
  *imp_lids  = *exp_lids  = NULL;
  *imp_procs = *exp_procs = NULL;

  /* Initialize HG parameters. */
  err = Zoltan_HG_Initialize_Params (zz, &hgp);
  if (err != ZOLTAN_OK)
    goto End;

  /* build initial Zoltan hypergraph from callback functions. */
  err = Zoltan_HG_Build_Hypergraph(zz, &zoltan_hg, &hgp);
  if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error building hypergraph.");
    goto End;
  }
   
  zz->LB.Data_Structure = zoltan_hg;
  nVtx = zoltan_hg->HG.nVtx;
  zoltan_hg->HG.redl = hgp.redl;     /* redl needs to be dynamic */
 
  /* allocate output partition memory */
  output_parts = (Partition) ZOLTAN_MALLOC (nVtx * sizeof(int));
  if (nVtx && output_parts == NULL) {
    err = ZOLTAN_MEMERR;
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    goto End;
  }

  hgp.kway = ((strstr(hgp.local_str,"kway")) ? 1 : 0);  
  if (hgp.kway) {
    err = Zoltan_HG_HPart_Lib (zz, &zoltan_hg->HG, zz->LB.Num_Global_Parts,
     output_parts, &hgp, 0);
    if (err != ZOLTAN_OK)
      return err;
  }
  else {
    /* vmap associates original vertices to sub hypergraphs */
    zoltan_hg->HG.vmap = (int*)ZOLTAN_MALLOC(zoltan_hg->HG.nVtx * sizeof (int));
    if (zoltan_hg->HG.vmap == NULL)  {
      err = ZOLTAN_MEMERR;
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      goto End;
    }
    for (i = 0; i < zoltan_hg->HG.nVtx; i++)
      zoltan_hg->HG.vmap[i] = i;

    /* tighten balance tolerance for recursive bisection process */
    if (zz->LB.Num_Global_Parts > 2)
      hgp.bal_tol = pow (hgp.bal_tol,
       1.0 / ceil (log((double)zz->LB.Num_Global_Parts) / log(2.0)));

    /* partition hypergraph */
    err = Zoltan_HG_rdivide (1, zz->LB.Num_Global_Parts, output_parts, zz,
     &zoltan_hg->HG, &hgp, 0);
    if (err != ZOLTAN_OK)  {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error partitioning hypergraph.");
      goto End;
    }
    ZOLTAN_FREE (&zoltan_hg->HG.vmap);
  }

  /* Build Zoltan's return arguments. */
  Zoltan_HG_Return_Lists(zz, zoltan_hg, output_parts, num_exp, exp_gids,
   exp_lids, exp_procs, exp_to_part);

End:
  if (err == ZOLTAN_MEMERR)
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
  ZOLTAN_FREE((void**) &output_parts);
  Zoltan_HG_Free_Structure(zz);
  ZOLTAN_TRACE_EXIT(zz, yo);
  return err;
}

/*****************************************************************************/



void Zoltan_HG_Free_Structure(ZZ *zz)
{
  /* frees all data associated with LB.Data_Structure for hypergraphs */
  ZHG *zoltan_hg = (ZHG *)(zz->LB.Data_Structure);

  if (zoltan_hg != NULL) {
    Zoltan_Multifree(__FILE__, __LINE__, 3, &zoltan_hg->Global_IDs,
     &zoltan_hg->Local_IDs, &zoltan_hg->Parts);
    Zoltan_HG_HGraph_Free(&zoltan_hg->HG);
    ZOLTAN_FREE((void**) &zz->LB.Data_Structure);
  }
}

/*****************************************************************************/



static int Zoltan_HG_Initialize_Params(
  ZZ *zz,   /* the Zoltan data structure */
  HGPartParams *hgp
)
{
  Zoltan_Bind_Param(HG_params,"HG_OUTPUT_LEVEL",    (void*) &hgp->output_level);
  Zoltan_Bind_Param(HG_params,"HG_REDUCTION_LIMIT",     (void*) &hgp->redl);
  Zoltan_Bind_Param(HG_params,"HG_REDUCTION_METHOD",    (void*) hgp->redm_str);
  Zoltan_Bind_Param(HG_params,"HG_EDGE_WEIGHT_SCALING", (void*) &hgp->ews);
  Zoltan_Bind_Param(HG_params,"HG_GLOBAL_PARTITIONING", (void*) hgp->global_str);
  Zoltan_Bind_Param(HG_params,"HG_LOCAL_REFINEMENT",    (void*) hgp->local_str);
  Zoltan_Bind_Param(HG_params,"CHECK_GRAPH",          (void*) &hgp->check_graph);
  Zoltan_Bind_Param(HG_params,"HG_REDUCTION_LOCAL_IMPROVEMENT",
                              (void*) hgp->redmo_str);

  /* Set default values */
  strncpy(hgp->redm_str,   "no", MAX_PARAM_STRING_LEN);
  strncpy(hgp->redmo_str,  "no", MAX_PARAM_STRING_LEN);
  strncpy(hgp->global_str, "no", MAX_PARAM_STRING_LEN);
  strncpy(hgp->local_str,  "no", MAX_PARAM_STRING_LEN);
  
  hgp->ews = 1;
  hgp->check_graph = 1;
  hgp->bal_tol = zz->LB.Imbalance_Tol[0];
  hgp->redl = zz->LB.Num_Global_Parts;
  hgp->output_level = HG_DEBUG_LIST;

  /* Get application values of parameters. */
  Zoltan_Assign_Param_Vals(zz->Params, HG_params, zz->Debug_Level, zz->Proc,
   zz->Debug_Proc);

  /* Convert strings to function pointers. */
  return Zoltan_HG_Set_Part_Options(zz, hgp);
}

/*****************************************************************************/



int Zoltan_HG_Set_Param(
  char *name,                     /* name of variable */
  char *val)                      /* value of variable */
{
  /* associates value to named variable for hypergraph partitioning parameters */
  PARAM_UTYPE result;         /* value returned from Check_Param */
  int index;                  /* index returned from Check_Param */

  return Zoltan_Check_Param(name, val, HG_params, &result, &index);
}

/*****************************************************************************/



static int Zoltan_HG_Return_Lists(
  ZZ *zz,
  ZHG *zoltan_hg,
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
  int eproc;
  int num_gid_entries   = zz->Num_GID;
  int num_lid_entries   = zz->Num_LID;
  int nVtx              = zoltan_hg->HG.nVtx;
  Partition input_parts = zoltan_hg->Parts;
  ZOLTAN_ID_PTR gids    = zoltan_hg->Global_IDs;
  ZOLTAN_ID_PTR lids    = zoltan_hg->Local_IDs;
  char *yo = "Zoltan_HG_Return_Lists";

  if (zz->LB.Return_Lists) {
    /* Count number of objects with new partitions or new processors. */
    *num_exp = 0;
    for (i = 0; i < nVtx; i++) {
      eproc = Zoltan_LB_Part_To_Proc(zz, output_parts[i], 
                                     &gids[i*num_gid_entries]);
      if (output_parts[i] != input_parts[i] || zz->Proc != eproc)
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
          return ZOLTAN_MEMERR;
       }

      for (j = 0, i = 0; i < nVtx; i++) {
        eproc = Zoltan_LB_Part_To_Proc(zz, output_parts[i], 
         &gids[i*num_gid_entries]);
        if (output_parts[i] != input_parts[i] || eproc != zz->Proc) {
          ZOLTAN_SET_GID(zz, &((*exp_gids)[j*num_gid_entries]),
           &(gids[i*num_gid_entries]));
          if (num_lid_entries > 0)
            ZOLTAN_SET_LID(zz, &((*exp_lids)[j*num_lid_entries]),
                               &(lids[i*num_lid_entries]));
          (*exp_procs)  [j] = eproc;
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
  ZZ *zz,          /* the Zoltan data structure */
  ZHG *zoltan_hg,
  HGraph *hg,
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
  char *yo = "Zoltan_HG_HGraph_Print";

  if (zoltan_hg != NULL  &&  hg != &zoltan_hg->HG) {
    ZOLTAN_PRINT_WARN(zz->Proc, yo, "Input hg != Zoltan HG");
    return;
  }

  Zoltan_Print_Sync_Start(zz->Communicator, 1);

  /* Print Vertex Info */
  fprintf(fp, "%s Proc %d\n", yo, zz->Proc);
  fprintf(fp, "Vertices (GID, LID, index)\n");
  for (i = 0; i < hg->nVtx; i++) {
    fprintf(fp, "(");
    ZOLTAN_PRINT_GID(zz, &(zoltan_hg->Global_IDs[i * num_gid]));
    fprintf(fp, ", ");
    ZOLTAN_PRINT_LID(zz, &(zoltan_hg->Local_IDs [i * num_lid]));
    fprintf(fp, ", %d)\n", i);
  }

  Zoltan_HG_Print(zz, hg, fp);
  Zoltan_Print_Sync_End(zz->Communicator, 1);
}



/* static double hcut_size_links (ZZ *zz, HGraph *hg, int p, Partition part)
{
  int i, j, *parts, nparts;
  double cut = 0.0;
  char *yo = "hcut_size_links";

  if (!(parts = (int*) ZOLTAN_CALLOC (p, sizeof(int)))) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }

  for (i = 0; i < hg->nEdge; i++) {
    nparts = 0;
    for (j = hg->hindex[i]; j < hg->hindex[i+1]; j++) {
      if (parts[part[hg->hvertex[j]]] < i + 1)
        nparts++;
      parts[part[hg->hvertex[j]]] = i + 1;
    }
    cut += (nparts-1) * (hg->ewgt ? hg->ewgt[i] : 1.0);
  }
  ZOLTAN_FREE ((void**) &parts);
  return cut;
}
*/

/* static double hcut_size_total (HGraph *hg, Partition part)
{
  int i, j, hpart;
  double cut = 0.0;

  for (i = 0; i < hg->nEdge; i++) {
    hpart = part[hg->hvertex[hg->hindex[i]]];
    for (j = hg->hindex[i] + 1; j < hg->hindex[i+1]
     && part[hg->hvertex[j]] == hpart; j++);
      if (j != hg->hindex[i+1])
        cut += (hg->ewgt ? hg->ewgt[i] : 1.0);
  }
  return cut;
}
*/

/*****************************************************************************/
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

