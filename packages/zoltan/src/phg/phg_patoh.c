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

#include <limits.h>
#include "zz_const.h"
#include "phg.h"

#ifdef ZOLTAN_PATOH

#include "patoh.h"
static int scale_round_weights(float *, int *, int, int, int);
#define ZOLTAN_PATOH_ERROR(str, err) \
  {ZOLTAN_PRINT_ERROR(zz->Proc, yo, str); ierr = err; goto End;}

#endif  /* ZOLTAN_PATOH */


/*****************************************************************************/
int Zoltan_PHG_PaToH(
  ZZ *zz,
  HGraph *hg,
  int nparts,           /* # of desired partitions */
  int *partvec,         /* Output:  partition assignment vector */
  PHGPartParams *hgp     /* Input: hypergraph parameters */
)
{
char *yo = "Zoltan_HG_PaToH";

#ifndef ZOLTAN_PATOH
  ZOLTAN_PRINT_ERROR(zz->Proc, yo, "PaToH method selected but Zoltan is not"
                     "built and linked with PaToH.");
  return ZOLTAN_FATAL;
}

#else

  int ierr = ZOLTAN_OK;
  PaToH_Parameters pargs;
  int *ivwgts = NULL, *iewgts = NULL;   /* PaToH expects integer weights. */
  int *partweight = NULL;               /* Diagnostics from PaToH */
  int cut;                              /* Diagnostics from PaToH */
  static int cnt = 0;                   /* Counter used to set RNG seed */

  /* PaToH expects integer weights; convert if weights are provided. */
  if (hg->VtxWeightDim) {
    ivwgts = (int *) ZOLTAN_MALLOC(hg->nVtx * hg->VtxWeightDim * sizeof(int));
    if (!ivwgts) ZOLTAN_PATOH_ERROR("Memory error.", ZOLTAN_MEMERR);
    scale_round_weights(hg->vwgt, ivwgts, hg->nVtx, hg->VtxWeightDim, 0);
  }

  if (hg->EdgeWeightDim) {
    if (hg->EdgeWeightDim > 1) 
      ZOLTAN_PATOH_ERROR("PaToH supports Edge_Weight_Dim == 0 or 1 only.",
                         ZOLTAN_FATAL);
    iewgts = (int *) ZOLTAN_MALLOC(hg->nEdge*hg->EdgeWeightDim*sizeof(int));
    if (!iewgts) ZOLTAN_PATOH_ERROR("Memory error.", ZOLTAN_MEMERR);
    scale_round_weights(hg->ewgt, iewgts, hg->nEdge, hg->EdgeWeightDim, 0);
  }


  PaToH_Initialize_Parameters(&pargs, PATOH_CONPART, PATOH_SUGPARAM_DEFAULT);

  pargs._k = nparts;
    /* set the same imbalanace ratio to the required one */
  pargs.init_imbal = pargs.final_imbal = hgp->bal_tol-1.0;
  if (hgp->patoh_alloc_pool0>0)
    pargs.MemMul_CellNet = hgp->patoh_alloc_pool0; 
  if (hgp->patoh_alloc_pool1>0)
    pargs.MemMul_Pins = hgp->patoh_alloc_pool1; 

  PaToH_Alloc(&pargs, hg->nVtx, hg->nEdge, hg->VtxWeightDim, 
              ivwgts, iewgts, hg->hindex, hg->hvertex);


  cnt = (cnt > 12345 ? 0 : cnt + 1);  /* Reset seed for call to PaToH. */
  pargs.seed = cnt;                   /* Differ each call to allow     */
                                      /* randomized testing.           */

  partweight = (int *) ZOLTAN_MALLOC(nparts * MAX(1, hg->VtxWeightDim) 
                                     * sizeof(int));
  if (!partweight)
    ZOLTAN_PATOH_ERROR("Memory error.", ZOLTAN_MEMERR);

  if (hg->VtxWeightDim <= 1){
    if (hgp->UseFixedVtx){
      /* Copy fixed vertices from hg->fixed_part */
      memcpy(partvec, hg->fixed_part, hg->nVtx*sizeof(int) );
      PaToH_Partition_with_FixCells(&pargs, hg->nVtx, hg->nEdge, 
                    ivwgts, iewgts, hg->hindex,
                    hg->hvertex, partvec, partweight, &cut);
      }
    else
      PaToH_Partition(&pargs, hg->nVtx, hg->nEdge, ivwgts, iewgts, hg->hindex,
                    hg->hvertex, partvec, partweight, &cut);
  }
  else 
    PaToH_MultiConst_Partition(&pargs, hg->nVtx, hg->nEdge, hg->VtxWeightDim,
                               ivwgts, hg->hindex, hg->hvertex, partvec,
                               partweight, &cut);
  

  /* HERE:  Check whether imbalance criteria were met. */

  /* HERE:  Allow remapping? */

End:

  ZOLTAN_FREE(&ivwgts);
  ZOLTAN_FREE(&iewgts);
  ZOLTAN_FREE(&partweight);
  PaToH_Free();


  return ierr;
}

/*****************************************************************************/

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


#endif /* ZOLTAN_PATOH */


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
