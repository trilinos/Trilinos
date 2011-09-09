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




#include <ctype.h>
#include "zz_const.h"
#include "zz_util_const.h"
#include "all_allo_const.h"
#include "params_const.h"
#include "order_const.h"
#include "third_library.h"


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* Auxiliary function prototypes. */
static int
Zoltan_Postprocess_UnScatter_Graph (ZZ *zz,
				    ZOLTAN_Third_Graph *gr,
				    ZOLTAN_Third_Part *prt,
				    indextype **rank);
static int
Zoltan_Postprocess_Order (ZZ *zz,
			  ZOLTAN_Third_Graph *gr,
			  ZOLTAN_Output_Order *ord);

static int
Zoltan_Postprocess_Partition (ZZ *zz,
			      ZOLTAN_Third_Graph *gr,
			      ZOLTAN_Third_Part  *prt,
			      ZOLTAN_Output_Part *part,
			      ZOLTAN_ID_PTR       global_ids,
			      ZOLTAN_ID_PTR       local_ids);

static int Compute_Bal(ZZ *, int, weighttype *, int, indextype *, double *);
static int Compute_EdgeCut(ZZ *, int, indextype *, float *, indextype *, int *, double *);
static float Compute_NetCut(ZZ *, int, indextype *, float *, indextype *, int *);
static float Compute_ConCut(ZZ *, int, indextype *, float *, indextype *, int *);
static int Compute_Adjpart(ZZ *, int, indextype *, indextype *, indextype *, int *, indextype *, int *);



int Zoltan_Postprocess_Graph(
  ZZ *zz,                               /* Zoltan structure */
  ZOLTAN_ID_PTR      global_ids,
  ZOLTAN_ID_PTR      local_ids,
  ZOLTAN_Third_Graph *gr,              /* Graph for third part libs */
  ZOLTAN_Third_Geom  *geo,
  ZOLTAN_Third_Part  *prt,
  ZOLTAN_Third_Vsize *vsp,
  ZOLTAN_Output_Order *ord,
  ZOLTAN_Output_Part  *part)
{
  int ierr = ZOLTAN_OK;

  static char * yo = "Zoltan_Postprocess_Graph";

  if (gr->scatter > 0) {                    /* Graph has been scattered */
    indextype *rank = NULL;

    if (ord) rank = ord->rank;
    else rank = NULL;
    ierr = Zoltan_Postprocess_UnScatter_Graph (zz, gr, prt, &rank);
    if (ierr) {
      ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL,
			 "Error returned from Zoltan_Postprocess_UnScatter_Graph");
    }
    if (ord) ord->rank = rank;
  }

  if (ord) {                                /* We have done ordering  */
      ierr = Zoltan_Postprocess_Order (zz, gr, ord);
      if (ierr == ZOLTAN_FATAL) {
	ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL,
			      "Error returned from Zoltan_Postprocess_Order");
      }
  }

  if (part) {                               /* We have done partitioning */
      ierr = Zoltan_Postprocess_Partition (zz, gr, prt, part, global_ids, local_ids);
      if (ierr == ZOLTAN_FATAL) {
	ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL,
			      "Error returned from Zoltan_Postprocess_Partition");
      }
  }

  return (ierr);

}


static int
Zoltan_Postprocess_UnScatter_Graph (ZZ *zz,
				    ZOLTAN_Third_Graph *gr,
				    ZOLTAN_Third_Part *prt,
				    indextype **rank)
{
  static char * yo = "Zoltan_Postprocess_UnScatter_Graph";

  int ierr = ZOLTAN_FATAL;
  indextype *src;
  indextype *dst;

  if (gr->scatter >0){
    gr->num_obj = gr->num_obj_orig;
    if (*rank) {                        /* We have to project back rank */
      dst = (indextype*) ZOLTAN_MALLOC(gr->num_obj*sizeof(indextype));
      src = *rank;
    }
    else {
      dst = prt->part_orig;
      src = prt->part;
    }

    ierr = Zoltan_Comm_Do_Reverse(gr->comm_plan, TAG2, (char *) src,
				  sizeof(indextype), NULL, (char *) dst);
    if ((ierr == ZOLTAN_FATAL) || (ierr == ZOLTAN_MEMERR)){
      ZOLTAN_THIRD_ERROR(ierr, "Zoltan_Comm_Do_Reverse returned error.");
    }
    Zoltan_Comm_Destroy(&gr->comm_plan); /* Destroy the comm. plan */
    /* We don't need the partition array with the scattered distribution
     * any more */
    ZOLTAN_FREE(&src);
    if (prt) {
      /* part is now the new partition array under the original distribution */
      prt->part = prt->part_orig;
      prt->part_orig = NULL;
    }
    else {
      *rank = dst;
    }

  }
  return (ierr);
}

static int
Zoltan_Postprocess_Order (ZZ *zz,
			  ZOLTAN_Third_Graph *gr,
			  ZOLTAN_Output_Order *ord)
{
  int ierr = ZOLTAN_OK;
  int i;

  /* Ordering */
  /* ParMetis produces the rank vector in Zoltan lingo */

  if (ord->rank != NULL) {
    /* Check if start_index != 0 */
    if (ord->order_opt && ord->order_opt->start_index) {
      int start_index;

      start_index = ord->order_opt->start_index;
      for (i=0; i<gr->num_obj; i++){
	ord->rank[i] +=  start_index;
      }
    }
  }
  else {
    ZOLTAN_PRINT_WARN(zz->Proc, __func__, "rank is NULL, no data returned");
    ierr = ZOLTAN_WARN;
  }

  /* If we did local ordering via METIS, then we also have the inv. perm. */
  if ((gr->graph_type == LOCAL_GRAPH) && (ord->iperm != NULL)){
      int start_index;

      start_index = ord->order_opt->start_index;
      for (i=0; i<gr->num_obj; i++){
	ord->iperm[i] += start_index;
      }
      /* EBEB: Return parameter that says we have computed both return args? */
  }

  /* Fill in the Zoltan Order Struct */
  /* EBEB: For now, discard separator info */
  if (0){ /* DEBUG: Print separator sizes to file */
    FILE *fp;
    fp = fopen("separators.txt", "w");
    fprintf(fp, "%i\n", ord->num_part);
    for (i=0; i<ord->num_part; ++i){
      fprintf(fp, TPL_IDX_SPEC " ", ord->sep_sizes[i]);
    }
    fprintf(fp, "\n");
    for (i=ord->num_part; i<2*ord->num_part-1; ++i){
      fprintf(fp, TPL_IDX_SPEC " ", ord->sep_sizes[i]);
    }
    fprintf(fp, "\n");
    fclose(fp);
  }

  return (ierr);
}

static int
Zoltan_Postprocess_Partition (ZZ *zz,
			      ZOLTAN_Third_Graph *gr,
			      ZOLTAN_Third_Part  *prt,
			      ZOLTAN_Output_Part *part,
			      ZOLTAN_ID_PTR      global_ids,
			      ZOLTAN_ID_PTR      local_ids)
{
  static char * yo = "Zoltan_Postprocess_Partition";

  int ierr = ZOLTAN_OK;
  int i, j, nsend;
  int *newproc, *tmp_part, *tmp_input_part;

  int num_gid_entries = zz->Num_GID;
  int num_lid_entries = zz->Num_LID;

  /* Partitioning */
  /* Determine new processor and number of objects to export */
  newproc = (int *) ZOLTAN_MALLOC(gr->num_obj * sizeof(int));
  if (gr->num_obj && !newproc){
    /* Not enough memory */
    ZOLTAN_THIRD_ERROR(ZOLTAN_MEMERR, "Out of memory. ");
  }
  for (i=0; i<gr->num_obj; i++){
    newproc[i] = Zoltan_LB_Part_To_Proc(zz, (int)prt->part[i],
					&(global_ids[i*num_gid_entries]));
    if (newproc[i]<0){
      ZOLTAN_FREE(&newproc);
      ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL,
			    "Zoltan_LB_Part_To_Proc returned invalid processor number.");
    }
  }

  if (zz->LB.Remap_Flag) {
    int new_map;

    if (sizeof(indextype) == sizeof(int)){
      ierr = Zoltan_LB_Remap(zz, &new_map, gr->num_obj, newproc, (int *)prt->input_part,
			   (int *)prt->part, 1);
    }
    else{
      tmp_part = (int *)ZOLTAN_MALLOC(sizeof(int) * gr->num_obj);
      tmp_input_part = (int *)ZOLTAN_MALLOC(sizeof(int) * gr->num_obj);

      if (gr->num_obj && (!tmp_part || !tmp_input_part)){
	ZOLTAN_THIRD_ERROR(ZOLTAN_MEMERR, "Not enough memory.");
      }
      for (i=0; i < gr->num_obj; i++){
        tmp_part[i] = (int)prt->part[i];
        tmp_input_part[i] = (int)prt->input_part[i];
      }
      ierr = Zoltan_LB_Remap(zz, &new_map, gr->num_obj, newproc, tmp_input_part, tmp_part, 1);

      for (i=0; i < gr->num_obj; i++){
        prt->part[i] = (indextype)tmp_part[i];
        prt->input_part[i] = (indextype)tmp_input_part[i];
      }

      ZOLTAN_FREE(&tmp_part);
      ZOLTAN_FREE(&tmp_input_part);
    }

    if (ierr < 0) {
      ZOLTAN_FREE(&newproc);
      ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL,
			    "Error returned from Zoltan_LB_Remap");
    }
  }

  nsend = 0;
  for (i=0; i<gr->num_obj; i++){
    if ((prt->part[i] != prt->input_part[i]) || ((!part->compute_only_part_changes) &&
						  (newproc[i] != zz->Proc)))
      nsend++;
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
      printf("[%1d] DEBUG: local object %1d: old part = " TPL_IDX_SPEC ", new part = " TPL_IDX_SPEC "\n",
	     zz->Proc, i, prt->input_part[i], prt->part[i]);
  }

  /* Create export lists */
  if (zz->LB.Return_Lists){
    if (zz->LB.Return_Lists == ZOLTAN_LB_CANDIDATE_LISTS) {
      ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL, "Candidate Lists not supported in GRAPH;"
                                       "change RETURN_LISTS parameter.");
    }
    part->num_exp = nsend;
    if (nsend > 0) {
      if (!Zoltan_Special_Malloc(zz,(void **)part->exp_gids,nsend,ZOLTAN_SPECIAL_MALLOC_GID)) {
	ZOLTAN_THIRD_ERROR(ZOLTAN_MEMERR, "Not enough memory.");
      }
      if (!Zoltan_Special_Malloc(zz,(void **)part->exp_lids,nsend,ZOLTAN_SPECIAL_MALLOC_LID)) {
	Zoltan_Special_Free(zz,(void **)part->exp_gids,ZOLTAN_SPECIAL_MALLOC_GID);
	ZOLTAN_THIRD_ERROR(ZOLTAN_MEMERR, "Not enough memory.");
      }
      if (!Zoltan_Special_Malloc(zz,(void **)part->exp_procs,nsend,ZOLTAN_SPECIAL_MALLOC_INT)) {
	Zoltan_Special_Free(zz,(void **)part->exp_lids,ZOLTAN_SPECIAL_MALLOC_LID);
	Zoltan_Special_Free(zz,(void **)part->exp_gids,ZOLTAN_SPECIAL_MALLOC_GID);
	ZOLTAN_THIRD_ERROR(ZOLTAN_MEMERR, "Not enough memory.");
      }
      if (!Zoltan_Special_Malloc(zz,(void **)part->exp_part,nsend,ZOLTAN_SPECIAL_MALLOC_INT)) {
	Zoltan_Special_Free(zz,(void **)part->exp_lids,ZOLTAN_SPECIAL_MALLOC_LID);
	Zoltan_Special_Free(zz,(void **)part->exp_gids,ZOLTAN_SPECIAL_MALLOC_GID);
	  Zoltan_Special_Free(zz,(void **)part->exp_procs,ZOLTAN_SPECIAL_MALLOC_INT);
	  ZOLTAN_THIRD_ERROR(ZOLTAN_MEMERR, "Not enough memory.");
      }
      j = 0;
      for (i=0; i<gr->num_obj; i++){
	if ((prt->part[i] != prt->input_part[i]) || ((!part->compute_only_part_changes)
						      && (newproc[i] != zz->Proc))){
	  /* Object should move to new partition or processor */
	  ZOLTAN_SET_GID(zz, &((*(part->exp_gids))[j*num_gid_entries]),
			 &(global_ids[i*num_gid_entries]));
	  if (num_lid_entries)
	    ZOLTAN_SET_LID(zz, &((*(part->exp_lids))[j*num_lid_entries]),
			   &(local_ids[i*num_lid_entries]));
	  (*(part->exp_part))[j] = (int)prt->part[i];
	  (*(part->exp_procs))[j] = newproc[i];
/*	  printf("[%1d] Debug: Move object %1d to part %1d, proc %1d\n", */
/*	     zz->Proc, i, prt->part[i], newproc[i]); */
	  j++;
	}
      }
    }
  }

  ZOLTAN_FREE(&newproc);

  return (ZOLTAN_OK);
}


int
Zoltan_Postprocess_FinalOutput (ZZ* zz, ZOLTAN_Third_Graph *gr,
				ZOLTAN_Third_Part *prt, ZOLTAN_Third_Vsize *vsp,
				int use_timers, realtype itr)
{
#define FOMAXDIM 10
  static char * yo = "Zoltan_Postprocess_FinalOutput";
  static int nRuns=0;
  static double balsum[FOMAXDIM], cutesum[FOMAXDIM];
  static double balmax[FOMAXDIM], cutemax[FOMAXDIM];
  static double balmin[FOMAXDIM], cutemin[FOMAXDIM];
  /* following variables are defined double to avoid overflow */
  static double cutlsum = 0.0, cutnsum = 0.0, movesum = 0.0, repartsum = 0.0;
  static float cutlmax = 0, cutnmax = 0;
  static double movemax = 0, repartmax = 0;
  static float cutlmin = INT_MAX, cutnmin = INT_MAX;
  static double movemin = 1e100, repartmin = 1e100;
  static int timer_final_output = -1;
  double bal[FOMAXDIM];    /* Balance:  max / avg */
  double cute[FOMAXDIM];   /* Traditional weighted graph edge cuts */
  float cutl;   /* Connnectivity cuts:  sum_over_edges((npart-1)) */
  float cutn;   /* Net cuts:  sum_over_edges((nparts>1)) */
  double move = 0.0, gmove =0.0;   /* migration cost */
  double repart; /* total repartitioning cost; cutl x multiplier + move */
  int *adjpart = NULL;
  int vdim;
  int edim;
  indextype *vsizeBACKUP = NULL;
  indextype *input_part;
  int i,rc;

/* #define UVC_DORUK_COMP_OBJSIZE */
#ifdef UVC_DORUK_COMP_OBJSIZE
  double minD, maxD, gminD, gmaxD;
#endif
  if (use_timers) {
    if (timer_final_output < 0)
      timer_final_output = Zoltan_Timer_Init(zz->ZTime, 1, "Final_Output");
    ZOLTAN_TIMER_START(zz->ZTime, timer_final_output, zz->Communicator);
  }
  if (nRuns == 0) {
    for (i = 0; i < FOMAXDIM; i++) {
      /* Initialize first time */
      balsum[i] = cutesum[i] = 0.0;
      balmax[i] = cutemax[i] = 0.0;
      balmin[i] = cutemin[i] = 1e100;
    }
  }

  if (vsp && vsp->vsizeBACKUP)
    vsizeBACKUP = vsp->vsizeBACKUP;
  vdim = MAX(gr->obj_wgt_dim,1);
  edim = MAX(zz->Edge_Weight_Dim,1);

  if (gr->obj_wgt_dim < FOMAXDIM && zz->Edge_Weight_Dim < FOMAXDIM) {
    if (gr->xadj[gr->num_obj]){
      adjpart = (int *) ZOLTAN_MALLOC(gr->xadj[gr->num_obj] * sizeof(int));
      if (!adjpart){
        ZOLTAN_THIRD_ERROR(ZOLTAN_MEMERR,
			 "Error 1 returned from Zoltan_Postprocess_FinalOutput");
      }
    }

    Compute_Bal(zz, gr->num_obj, gr->vwgt, gr->obj_wgt_dim, prt->part, bal);

    rc = Compute_Adjpart(zz, gr->num_obj, gr->vtxdist, gr->xadj, gr->adjncy,
                    gr->adjproc, prt->part, adjpart);

    if (rc != ZOLTAN_OK){
      ZOLTAN_THIRD_ERROR(ZOLTAN_MEMERR,
			 "Error 2 returned from Zoltan_Postprocess_FinalOutput");
    }
    Compute_EdgeCut(zz, gr->num_obj, gr->xadj, gr->float_ewgts,
                    prt->part, adjpart, cute);
    cutl = Compute_ConCut(zz, gr->num_obj, gr->xadj,  gr->float_ewgts,
                          prt->part, adjpart);
    cutn = Compute_NetCut(zz, gr->num_obj, gr->xadj, gr->float_ewgts,
                          prt->part, adjpart);

#ifdef UVC_DORUK_COMP_OBJSIZE
    if (vsizeBACKUP) {
      minD = vsizeBACKUP[0];
      maxD = vsizeBACKUP[0];
    }
#endif

  if (gr->scatter != 0) { /* Project input part in the scatter graph */
    input_part = (indextype*) ZOLTAN_MALLOC(gr->num_obj * sizeof(int));
    Zoltan_Comm_Do(gr->comm_plan, TAG1, (char *) prt->input_part, sizeof(indextype),
		   (char *) input_part);
  }
  else
    input_part = prt->input_part;

  for (i=0; i<gr->num_obj; i++) {
    /*printf("obj[%d] = %d\n", i, vsize[i]);*/
    if (prt->part[i] != input_part[i]) {
      move += (double) ((vsizeBACKUP) ? vsizeBACKUP[i] : 1.0);
    }
#ifdef UVC_DORUK_COMP_OBJSIZE
    if (vsizeBACKUP) {
      minD = minD < vsizeBACKUP[i] ? minD : vsizeBACKUP[i];
      maxD = maxD > vsizeBACKUP[i] ? maxD : vsizeBACKUP[i];
    }
#endif
  }

  if (gr->scatter != 0)
    ZOLTAN_FREE(&input_part);

#ifdef UVC_DORUK_COMP_OBJSIZE
    if (gr->showMoveVol) {
      MPI_Allreduce(&minD, &gminD, 1, MPI_DOUBLE, MPI_MIN, zz->Communicator);
      MPI_Allreduce(&maxD, &gmaxD, 1, MPI_DOUBLE, MPI_MAX, zz->Communicator);

      if (zz->Proc == 0)
	printf("minD: %f, maxD: %f, gminD: %f, gmaxD: %f\n", minD, maxD, gminD, gmaxD);
    }
#endif

    MPI_Allreduce(&move, &gmove, 1, MPI_DOUBLE, MPI_SUM, zz->Communicator);

    repart =  (itr) * (double) cutl + gmove;
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
    for (i = 0; i < vdim; i++) {
      balsum[i] += bal[i];
      if (bal[i] > balmax[i]) balmax[i] = bal[i];
      if (bal[i] < balmin[i]) balmin[i] = bal[i];
    }
    for (i = 0; i < edim; i++) {
      cutesum[i] += cute[i];
      if (cute[i] > cutemax[i]) cutemax[i] = cute[i];
      if (cute[i] < cutemin[i]) cutemin[i] = cute[i];
    }
    nRuns++;

    if (zz->Proc == 0) {
      for (i = 0; i < vdim; i++) 
	printf("STATS Runs %d  bal[%d]  CURRENT %f  MAX %f  MIN %f  AVG %f\n",
	       nRuns, i, bal[i], balmax[i], balmin[i], balsum[i]/nRuns);
      printf("STATS Runs %d  cutl CURRENT %f  MAX %f  MIN %f  AVG %f\n",
	     nRuns, cutl, cutlmax, cutlmin, cutlsum/nRuns);
      printf("STATS Runs %d  cutn CURRENT %f  MAX %f  MIN %f  AVG %f\n",
	     nRuns, cutn, cutnmax, cutnmin, cutnsum/nRuns);
      for (i = 0; i < edim; i++) 
	printf("STATS Runs %d  cute[%d] CURRENT %f  MAX %f  MIN %f  AVG %f\n",
	       nRuns, i, cute[i], cutemax[i], cutemin[i], cutesum[i]/nRuns);
      printf("STATS Runs %d  %s CURRENT %f  MAX %f  MIN %f  AVG %f\n",
	     nRuns, gr->showMoveVol ? "moveVol" : "moveCnt", gmove,
             movemax, movemin, movesum/nRuns);
      if (gr->showMoveVol)
	printf("STATS Runs %d  repart CURRENT %f  MAX %f  MIN %f  AVG %f\n",
	       nRuns, repart, repartmax, repartmin, repartsum/nRuns);
    }
    ZOLTAN_FREE(&adjpart);
  }
#undef FOMAXDIM
  if (use_timers)
    ZOLTAN_TIMER_STOP(zz->ZTime, timer_final_output, zz->Communicator);

  return (ZOLTAN_OK);
}


/****************************************************************************/
static int Compute_Bal(
  ZZ *zz,
  int nvtx,
  weighttype *vwgts,
  int obj_wgt_dim,
  indextype *parts,
  double *bal
)
{
/*
 * Compute the load balance of the computed partition.
 */
int i, j;
int dim = MAX(obj_wgt_dim, 1);
int size;
float *sum = NULL, *gsum = NULL;
float *tot = NULL, *max = NULL;

  size = zz->LB.Num_Global_Parts * dim;
  sum = (float *) ZOLTAN_CALLOC(2 * size, sizeof(float));
  gsum = sum + size;
  tot = (float *) ZOLTAN_CALLOC(2 * dim, sizeof(float));
  max = tot + dim;

  for (i = 0; i < nvtx; i++)
    for (j = 0; j < dim; j++)
      sum[parts[i]*dim + j] += (vwgts ? vwgts[i*dim + j] : 1.);

  MPI_Allreduce(sum, gsum, size, MPI_FLOAT, MPI_SUM, zz->Communicator);

  for (i = 0; i < zz->LB.Num_Global_Parts; i++)
    for (j = 0; j < dim; j++) {
      tot[j] += gsum[i*dim+j];
      if (gsum[i*dim+j] > max[j]) max[j] = gsum[i*dim+j];
    }

  for (j = 0; j < dim; j++)
    if (tot[j] > 0.)
      bal[j] = zz->LB.Num_Global_Parts * max[j] / tot[j];

  ZOLTAN_FREE(&sum);
  ZOLTAN_FREE(&tot);

  return ZOLTAN_OK;
}

/****************************************************************************/
static int Compute_EdgeCut(
  ZZ *zz,
  int nvtx,
  indextype *xadj,
  float *ewgts,
  indextype *parts,
  int *nborparts,
  double *cute
)
{
/*
 * Compute total weight of cut graph edges w.r.t. partitions.
 */
indextype i, j, k;
int dim = MAX(zz->Edge_Weight_Dim, 1);
double *cut = NULL;

  cut = (double *) ZOLTAN_CALLOC(dim, sizeof(double));

  for (i = 0; i < nvtx; i++)
    for (j = xadj[i]; j < xadj[i+1]; j++)
      if (parts[i] != nborparts[j])
        for (k = 0; k < dim; k++)
          cut[k] += (ewgts ? ewgts[j*dim+k] : 1.);

  MPI_Allreduce(cut, cute, dim, MPI_DOUBLE, MPI_SUM, zz->Communicator);
  ZOLTAN_FREE(&cut);
  return ZOLTAN_OK;
}

/****************************************************************************/
static float Compute_NetCut(
  ZZ *zz,
  int nvtx,
  indextype *xadj,
  float *ewgts,
  indextype *parts,
  int *nborparts
)
{
/*
 * Compute weight of hyperedges cut w.r.t. partitions.
 * Assume one hyperedge per vertex.
 * Equivalent to number of boundary vertices weighted by edge weights.
 */
indextype i, j;
float cutn = 0., gcutn = 0.;
int dim = zz->Edge_Weight_Dim;

  for (i = 0; i < nvtx; i++) {

    /* Compute the maximum graph edge weight of all edges incident to vertex */
    /* For now, use only first weight per edge. */
    float maxewgt = 0.;
    if (ewgts)
      for (j = xadj[i]; j < xadj[i+1]; j++) 
        if (ewgts[j*dim] > maxewgt) maxewgt = ewgts[j*dim];

    for (j = xadj[i]; j < xadj[i+1]; j++)
      if (parts[i] != nborparts[j]) {
        cutn += (ewgts ? maxewgt : 1.); 
        break;
      }
  }

  MPI_Allreduce(&cutn, &gcutn, 1, MPI_FLOAT, MPI_SUM, zz->Communicator);

  return gcutn;
}

/****************************************************************************/
static float Compute_ConCut(
  ZZ *zz,
  int nvtx,
  indextype *xadj,
  float *ewgts,
  indextype *parts,
  int *nborparts
)
{
/*
 * Compute SUM over hyperedges( (#parts/hyperedge - 1)*ewgt);
 * Assume one hyperedge per vertex.
 */
indextype i, j;
float cutl = 0., gcutl = 0.;
indextype *used = NULL;
int dim = zz->Edge_Weight_Dim;

  used = (indextype *) ZOLTAN_MALLOC(zz->LB.Num_Global_Parts * sizeof(indextype));
  for (i = 0; i < zz->LB.Num_Global_Parts; i++) used[i] = -1;

  for (i = 0; i < nvtx; i++) {
    /* Compute the maximum graph edge weight of all edges incident to vertex */
    /* For now, use only first weight per edge. */
    float maxewgt = 0.;
    if (ewgts)
      for (j = xadj[i]; j < xadj[i+1]; j++) 
        if (ewgts[j*dim] > maxewgt) maxewgt = ewgts[j*dim];

    used[parts[i]] = i;
    for (j = xadj[i]; j < xadj[i+1]; j++)
      if (used[nborparts[j]] < i) {
        used[nborparts[j]] = i;
        cutl += (ewgts ? maxewgt : 1.); 
      }
  }
  ZOLTAN_FREE(&used);

  MPI_Allreduce(&cutl, &gcutl, 1, MPI_FLOAT, MPI_SUM, zz->Communicator);

  return gcutl;
}

/****************************************************************************/
static int Compute_Adjpart(
  ZZ *zz,
  int nvtx,         /* Input:  # vtxs in this processor */
  indextype *vtxdist,     /* Input:  Distribution of vertices across processors */
  indextype *xadj,        /* Input:  Index of adjncy:  adjncy[xadj[i]] to 
                               adjncy[xadj[i]+1] are all edge nbors of vtx i. */
  indextype *adjncy,      /* Input:  Array of nbor vertices. */
  int *adjproc,     /* Input:  adjproc[j] == processor owning adjncy[j]. */
  indextype *part,        /* Input:  Partition assignments of vtxs. */
  int *adjpart      /* Output: adjpart[j] == partition owning adjncy[j] */
)
{
/* Given an adjacency list adjncy, find the partition number of each 
 * vertex in adjncy.  Return it in adjpart.
 */
ZOLTAN_COMM_OBJ *plan;
int i;
indextype start = vtxdist[zz->Proc];  /* First vertex on this processor */
indextype *recv_gno= NULL;
int *send_int=NULL;
int nrecv;
int tag = 24542;

  Zoltan_Comm_Create(&plan, (int)xadj[nvtx], adjproc, zz->Communicator, tag++, &nrecv);

  if (nrecv){
    recv_gno = (indextype *) ZOLTAN_MALLOC(nrecv * sizeof(indextype));
    send_int = (int *) ZOLTAN_MALLOC(nrecv * sizeof(int));
    if (!recv_gno || !send_int){
      return ZOLTAN_MEMERR;
    }
  }

  Zoltan_Comm_Do(plan, tag++, (char *) adjncy, sizeof(indextype), (char *) recv_gno);

  for (i = 0; i < nrecv; i++){
    send_int[i] = part[recv_gno[i] - start];
  }

  ZOLTAN_FREE(&recv_gno);

  Zoltan_Comm_Do_Reverse(plan, tag, (char *)send_int, sizeof(int), NULL, (char *) adjpart);

  ZOLTAN_FREE(&send_int);

  Zoltan_Comm_Destroy(&plan);

  return ZOLTAN_OK;
}


#ifdef __cplusplus
}
#endif
