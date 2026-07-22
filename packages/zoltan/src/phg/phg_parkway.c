// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "zz_const.h"
#include "phg.h"

#ifdef ZOLTAN_PARKWAY

    
void Zoltan_ParaPartKway(int numVertices, int numHedges, const int *vWeights, const int *hEdgeWts, const int *pinList, const int *offsets, int numParts, double constraint, int *k_1cut, const int *options, int *pVector, const char *outFile, MPI_Comm comm);    
    
static int scale_round_weights(float *, int *, int, int, int);

#define ZOLTAN_PARKWAY_ERROR(str, err) \
  {ZOLTAN_PRINT_ERROR(zz->Proc, yo, str); ierr = err; goto End;}

#endif  /* ZOLTAN_PARKWAY */


/*****************************************************************************/
int Zoltan_PHG_ParKway(
  ZZ        *zz,
  HGraph    *hg,
  int       nparts,           /* # of desired partitions */
  Partition partvec,          /* Output:  partition assignment vector */
  PHGPartParams *hgp          /* Input: hypergraph parameters */  
)
{
    char *yo = "Zoltan_HG_ParKway";

#ifndef ZOLTAN_PARKWAY
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "ParKway method selected but Zoltan is not"
                       "built and linked with ParKway.");
    return ZOLTAN_FATAL;
#else

    int ierr = ZOLTAN_OK;
    int options[29];                      /* ParKway options */
    int *ivwgts = NULL, *iewgts = NULL;   /* ParKway expects integer weights. */
    int *pvector=NULL;                    /* partvec for "local" vertices */
    int cut;                              /* Diagnostics from ParKway */
    double constraint;                    /* imbalance ratio */
    int i, anVtx, nVtx;                   /* counter and local vertex cnt (for first k-1 parts)*/
    PHGComm *hgc=hg->comm;
    int *disp=NULL, *recv_size=NULL;      /* for allgatherv */
    static int seed=1;
    
    /* ParKway expects integer weights; convert if weights are provided. */
    ivwgts = (int *) ZOLTAN_MALLOC(hg->nVtx  * sizeof(int));
    iewgts = (int *) ZOLTAN_MALLOC(hg->nEdge * sizeof(int));    
    if (!ivwgts || !iewgts)
        ZOLTAN_PARKWAY_ERROR("Memory error.", ZOLTAN_MEMERR);
    
    
    if (hg->VtxWeightDim > 1) { 
        ZOLTAN_PARKWAY_ERROR("ParKway supports Vtx_Weight_Dim == 0 or 1 only.",
                             ZOLTAN_FATAL);
    } else if (hg->VtxWeightDim == 1) 
        scale_round_weights(hg->vwgt, ivwgts, hg->nVtx, hg->VtxWeightDim, 0);
    else 
        for (i=0; i<hg->nVtx; ++i)
            ivwgts[i] = 1;

        
    if (hg->EdgeWeightDim > 1) {
        ZOLTAN_PARKWAY_ERROR("ParKway supports Edge_Weight_Dim == 0 or 1 only.",
                             ZOLTAN_FATAL);
    } else if (hg->EdgeWeightDim == 1) 
        scale_round_weights(hg->ewgt, iewgts, hg->nEdge, hg->EdgeWeightDim, 0);
    else
        for (i=0; i<hg->nEdge; ++i)
            iewgts[i] = 1;
    
    
    
    anVtx = hg->nVtx / hgc->nProc;
    nVtx = (hgc->myProc==hgc->nProc-1) ? hg->nVtx-(anVtx*(hgc->nProc-1)) : anVtx;

    pvector = (int *) ZOLTAN_MALLOC(nVtx * sizeof(int));
    disp = (int *) ZOLTAN_MALLOC(hgc->nProc * sizeof(int));
    recv_size = (int *) ZOLTAN_MALLOC(hgc->nProc * sizeof(int));
    if ((nVtx && !pvector) || !disp || !recv_size)
        ZOLTAN_PARKWAY_ERROR("Memory error.", ZOLTAN_MEMERR);


    /* ----- Set ParKway's options --------------- */
    options[0] = 1;/*0 -> all options use default, else user define*/
    options[1] = seed++;/*0 -> seed chosen by sprng, else use options[1] as seed*/
    options[2] = 0;/*0 -> no disp info, 1 -> some, 2 -> lots*/
    options[3] = 1;/*0 -> do not write partition to disk, 1 -> do write*/
    options[4] = 1;/*number of parallel runs*/
    options[5] = 0;/*vertex to processor allocation: 0 -> as read in, 1 -> random 2 -> as prescribed in partition file */
    options[6] = 100;/*hyperedge length percentile for approx para coarsening and refinement*/
    options[7] = 1;/*increment in percentile options[6]*/
    options[8] = 200;/*numParts*options[5] -> min number of vertices in coarse hypergraph*/
    options[9] = 7;/*[9] and [10] specify reduction ratio in parallel coarsening*/
    options[10] = 4;/*r = [9]/[10]*/
    options[11] = 3;/*vertex visit order: 3 -> random, 1/2 inc/dec by vertex id, 4/5 inc/dec by vertex wt*/
    options[12] = 3;/*divide connectivity by cluster weight/hyperedge length: 0-neither, 1-only cluster, 2-only hedge len, 3-both   */
    options[13] = 3;/*matching request resolution order: 3 -> random, 2 -> as they arrive */
    options[14] = 1;/*number serial partitioning runs*/
    
    options[15] = 5;/*serial partitioning routine, 1-3 RB, 4 khmetis, 5 patoh, see manual*/
    
    if (!strcasecmp(hgp->parkway_serpart, "patoh"))
        options[15] = 5;
    else if (!strcasecmp(hgp->parkway_serpart, "hmetis"))
        options[15] = 4;
    else if (!strcasecmp(hgp->parkway_serpart, "generic"))
        options[15] = 1;
    else if (!strcasecmp(hgp->parkway_serpart, "genericv"))
        options[15] = 2;
    else if (!strcasecmp(hgp->parkway_serpart, "genericmv"))
        options[15] = 3;
    else {
        ZOLTAN_PARKWAY_ERROR("Invalid ParKway serial partitioner. It should be one of; generic, genericv, genericmv, hmetis, patoh.", ZOLTAN_FATAL);
    }
    
    /* uprintf(hgc, "ParKway serpart='%s'  options[13]=%d\n", hgp->parkway_serpart, options[13]); */
    
    options[16] = 2;/*serial coarsening algorithm (only if [15] = RB, see manual)*/
    options[17] = 2;/*num bisection runs in RB (only if [15] = RB, see manual)*/
    options[18] = 10;/*num initial partitioning runs in RB (only if [13] = RB, see manual)*/
    options[19] = 2;/*hmetis_PartKway coarsening option, vals 1-5, see manual (only if [15] = 4)   */
    options[20] = 2;/*hmetis_PartKway refinement option, vals 0-3, see manual (only if [15] = 4)*/
    options[21] = 3;/*patoh_partition parameter settings, vals 1-3, see manual (only if [15] = 5)*/
    options[22] = 1;/*parallel uncoarsening algorithm, 1 simple, 2 only final V-Cycle, 3 all V-Cycle*/
    options[23] = 5;/*limit on number of V-Cycle iterations (only if [22] = 2/3)*/
    options[24] = 0;/*min allowed gain for V-Cycle (percentage, see manual, only if [21] = 2/3)*/
    options[25] = 0;/*percentage threshold used to reject partitions from a number of runs (see manual)*/
    options[26] = 0;/*reduction in [23] as partitions propagate by factor [24]/100 (see manual)*/
    options[27] = 100;/*early exit criterion in parallel refinement, will exit if see ([25]*num vert)/100 consecutive -ve moves */
    options[28] = 0;/*parallel refinement 0->basic, 1->use approx 2->use early exit 3->use approx and early exit  */
    
    constraint = hgp->bal_tol-1.0;
    
    Zoltan_ParaPartKway(nVtx, hg->nEdge, &ivwgts[hgc->myProc*anVtx], iewgts,
                        hg->hindex, hg->hvertex, nparts,
                        constraint, &cut, options, pvector, NULL, hgc->Communicator);
    
/* KDDKDD
   uprintf(hgc, "ParaPartKway cut=%d\n", cut);
*/
    
    
    /* after partitioning Zoltan needs partvec exist on all procs for nProc_x=1 */       
    disp[0] = 0; 
    for (i = 1; i < hgc->nProc; ++i)
        disp[i] = disp[i-1] + anVtx;
    
    MPI_Allgather (&nVtx, 1, MPI_INT, recv_size, 1, MPI_INT, hgc->Communicator);    
    MPI_Allgatherv(pvector, nVtx, MPI_INT, 
                  partvec, recv_size, disp, MPI_INT, hgc->Communicator);

    
  /* HERE:  Check whether imbalance criteria were met. */

End:

    Zoltan_Multifree(__FILE__,__LINE__, 5, &ivwgts, &iewgts, &pvector, &disp, &recv_size);
    
    return ierr;
#endif
}

    
/*****************************************************************************/

#ifdef ZOLTAN_PARKWAY

    
#define INT_EPSILON (1e-5)

static int scale_round_weights(
  float *fwgts, 
  int *iwgts, 
  int n, 
  int dim,
  int mode
)
{
/* Convert floating point weights to integer weights.
 * This routine is stolen from scale_round_weights in parmetis_jostle.c.
 * Because it needs to run only serially, and because it uses only 
 * integers (not idxtype), it has been largely duplicated here.
 */

  int i, j, tmp, ierr; 
  int max_wgt_sum = INT_MAX/8;
  int *nonint;
  float *scale, *sum_wgt, *max_wgt;
  char msg[256];
  static char *yo = "scale_round_weights";

  ierr = ZOLTAN_OK;

  if (mode == 0) {
    /* No scaling; just convert to int */
    for (i=0; i<n*dim; i++){
      iwgts[i] = (int) ceil((double) fwgts[i]);
    }
  }
  else{
    /* Allocate local arrays */
    nonint = (int *)ZOLTAN_MALLOC(dim*sizeof(int));
    scale = (float *)ZOLTAN_MALLOC(3*dim*sizeof(float));
    sum_wgt = scale + dim;
    max_wgt = sum_wgt + dim;
    if (!(nonint && scale)){
      ZOLTAN_PRINT_ERROR(0, yo, "Out of memory.");
      ZOLTAN_FREE(&nonint);
      ZOLTAN_FREE(&scale);
      return ZOLTAN_MEMERR;
    }

    /* Initialize */
    for (j=0; j<dim; j++){
      nonint[j] = 0;
      sum_wgt[j] = 0;
      max_wgt[j] = 0;
    }

    /* Compute local sums of the weights */
    /* Check if all weights are integers */
    for (i=0; i<n; i++){
      for (j=0; j<dim; j++){
        if (!nonint[j]){ 
          /* tmp = (int) roundf(fwgts[i]);  EB: Valid C99, but not C89 */
          tmp = (int) floor((double) fwgts[i] + .5); /* Nearest int */
          if (fabs((double)tmp-fwgts[i*dim+j]) > INT_EPSILON){
            nonint[j] = 1;
          }
        }
        sum_wgt[j] += fwgts[i*dim+j];
        if (fwgts[i*dim+j] > max_wgt[j])
          max_wgt[j] = fwgts[i*dim+j]; 
      }
    }

    /* Calculate scale factor */
    for (j=0; j<dim; j++){
      scale[j] = 1.;
      /* Scale unless all weights are integers (not all zero) */
      if (nonint[j] || (max_wgt[j] <= INT_EPSILON) 
                    || (sum_wgt[j] > max_wgt_sum)){
        if (sum_wgt[j] == 0){
          ierr = ZOLTAN_WARN;
          sprintf(msg, "All weights are zero in component %1d", j);
          ZOLTAN_PRINT_WARN(0, yo, msg);
        }
        else /* sum_wgt[j] != 0) */
          scale[j] = max_wgt_sum/sum_wgt[j];
      }
    }

    /* If mode==2, let the scale factor be the same for all weights */
    if (mode==2){
      for (j=1; j<dim; j++){
        if (scale[j]<scale[0])
          scale[0] = scale[j];
      }
      for (j=1; j<dim; j++){
        scale[j] = scale[0];
      }
    }

    /* Convert weights to positive integers using the computed scale factor */
    for (i=0; i<n; i++){
      for (j=0; j<dim; j++){
        iwgts[i*dim+j] = (int) ceil((double) fwgts[i*dim+j]*scale[j]);
      }
    }

    ZOLTAN_FREE(&nonint);
    ZOLTAN_FREE(&scale);
  }
  return ierr;
}


#endif /* ZOLTAN_PARKWAY */


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
