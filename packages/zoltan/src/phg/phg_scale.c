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
#include "phg.h"
#include "zz_const.h"



/****************************************************************************/

/* Scaling the weight of hypergraph edges. 
   This changes the inner product used in matching,
   usually to the better! 
   Note that the scaled weights are returned in a separate
   array (new_ewgts) and the hypergraph is not changed in this function.

   Currently, three scalings are available:
   1: absorption scaling (1/(size-1))
   2: net size scaling   (1/size)
   3: clique scaling     (2/(size*(size-1)))

   EBEB: Removed Robert's old serial scaling methods. 
         We should look at these later.
 */
int Zoltan_PHG_Scale_Edges (ZZ *zz, HGraph *hg, float *new_ewgt, 
                            int edge_scaling)
{
int    i, err;
int    *lsize = NULL;  /* local edge sizes */
int    *size = NULL;   /* edge sizes */
static char *yo = "Zoltan_PHG_Scale_Edges";

  err = ZOLTAN_OK; 

  /* allocate size arrays */
  if (hg->nEdge == 0)
    return ZOLTAN_OK;
  if (!(size  = (int *) ZOLTAN_MALLOC (sizeof(int) * hg->nEdge)) ||
      !(lsize = (int *) ZOLTAN_MALLOC (sizeof(int) * hg->nEdge)) ){
    ZOLTAN_FREE(&size);
    ZOLTAN_FREE(&lsize);
    ZOLTAN_PRINT_WARN(zz->Proc, yo, "Out of memory");
    return ZOLTAN_MEMERR;
  }

  switch (edge_scaling){

  case 0:
    /* copy current weights; no scaling. */
    for (i = 0; i < hg->nEdge; i++) 
        new_ewgt[i] = (hg->ewgt ? hg->ewgt[i] : 1.0);
    break;

  case 1:
    /* absorption scaling; scale by 1/(size -1) */
    /* intentionally fall through into next case! */
  case 2:
    /* net size scaling; scale by 1/size */
    /* intentionally fall through into next case! */
  case 3:
    /* clique scaling; scale by 2/(size*(size-1)) */

    /* first compute size of all hyperedges */
    for (i = 0; i < hg->nEdge; i++) {
      lsize[i] = hg->hindex[i+1] - hg->hindex[i];
    }

    /* sum up local sizes */
    /* assume SCHEMEA : all procs in a row have same # hyperedges */
    MPI_Allreduce(lsize, size, hg->nEdge, MPI_INT, MPI_SUM, 
                  hg->comm->row_comm);

#ifdef DEBUG_EB
  printf("DEBUG: hyperedge sizes =\n");
#endif

    /* scale edge weights */
    for (i = 0; i < hg->nEdge; i++) {
#ifdef DEBUG_EB
      printf(" %1d,", size[i]);
#endif
      if (size[i]>1) {
        if (edge_scaling==1)
          new_ewgt[i] = (hg->ewgt ? hg->ewgt[i] : 1.0) / (size[i]-1.0);
        else if (edge_scaling==2)
          new_ewgt[i] = (hg->ewgt ? hg->ewgt[i] : 1.0) / size[i];
        else if (edge_scaling==3)
          new_ewgt[i] = (hg->ewgt ? hg->ewgt[i] : 1.0) * 2.0 / 
                        (size[i]*(size[i]-1.0));
      }
      else /* size[i] == 1 */
        new_ewgt[i] = 0.0;
    }
#ifdef DEBUG_EB
    printf("\n");
#endif
    break;

  default:
    /* invalid scaling option */
    err = ZOLTAN_FATAL;
    break;
  }

  ZOLTAN_FREE(&size);
  ZOLTAN_FREE(&lsize);

  return err;
}

/**********************************************************************
  Scaling routine for vertices. This creates an array that is only
  used to modify the inner product in the matching.
***********************************************************************/

int Zoltan_PHG_Scale_Vtx (ZZ *zz, HGraph *hg, PHGPartParams *hgp)
{
  int i;
  int *ldegree, *gdegree; 
  char *yo = "Zoltan_PHG_Scale_Vtx";

  if ((hgp->vtx_scaling==0) || (hg->nVtx==0))
    return ZOLTAN_OK;

  /* See whether vtx_scal is large enough; nVtx may have increased due to
     processor reduction */
  if (hgp->vtx_scal && (hgp->vtx_scal_size < hg->nVtx)) {
    hgp->vtx_scal_size = 0;
    ZOLTAN_FREE(&(hgp->vtx_scal));
  }

  /* Allocate vtx_scal array if necessary */
  if (hgp->vtx_scal==NULL){  
    /* first level in V-cycle or ...*/
    /* nVtx increased due to processor reduction */
    hgp->vtx_scal_size = hg->nVtx;
    if (!(hgp->vtx_scal = (float*) ZOLTAN_MALLOC (hgp->vtx_scal_size *
                           sizeof(float)))) {
       hgp->vtx_scal_size = 0;
       ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
       return ZOLTAN_MEMERR;
    }
  }

  /* Allocate temp arrays */
  if (!(ldegree = (int *) ZOLTAN_MALLOC(2*hg->nVtx * sizeof(int)))){
    ZOLTAN_PRINT_WARN(zz->Proc, yo, "Out of memory");
    return ZOLTAN_MEMERR;
  }
  gdegree = ldegree + hg->nVtx;

  if (hgp->vtx_scaling<=2){

    /* First compute local vertex degrees. */
    for (i=0; i<hg->nVtx; i++){
       ldegree[i] = hg->vindex[i+1] - hg->vindex[i]; /* local degree */
    }
                                                                                
    /* Sum up along columns for global degrees. */
    MPI_Allreduce(ldegree, gdegree, hg->nVtx, MPI_INT, MPI_SUM,
          hg->comm->col_comm);

  }

#ifdef DEBUG_EB
  /* Debug */
  printf("DEBUG: vertex degrees =\n");
  for (i=0; i<hg->nVtx; i++)
    printf(" %1d,", gdegree[i]);
  printf("\n");
#endif

  /* Compute scale factor from degrees. */
  if (hgp->vtx_scaling==1){ /* cosine metric, scale by sqrt of degree */
    for (i=0; i<hg->nVtx; i++)
      if (gdegree[i] == 0)
         hgp->vtx_scal[i] = 1.0;
      else
         hgp->vtx_scal[i] = 1. / sqrt((double)gdegree[i]);
  }
  else if (hgp->vtx_scaling==2){  /* scale by degree */
    for (i=0; i<hg->nVtx; i++)
      if (gdegree[i] == 0)
         hgp->vtx_scal[i] = 1.0;
      else    
         hgp->vtx_scal[i] = 1. / gdegree[i];
  }
  else if (hgp->vtx_scaling==3){  /* scale by sqrt vertex weights */
    if (hg->vwgt)
      for (i=0; i<hg->nVtx; i++)  {
         /* KDD Note:  Scaling by only first weight */
         if (hg->vwgt[i*hg->VtxWeightDim] == 0)
            hgp->vtx_scal[i] = 1.0;
         else      
            /* KDD Note:  Scaling by only first weight */
            hgp->vtx_scal[i] = 1. / sqrt((double)hg->vwgt[i*hg->VtxWeightDim]);
      }
  }
  else if (hgp->vtx_scaling==4){  /* scale by vertex weights */
    if (hg->vwgt)
      for (i=0; i<hg->nVtx; i++)  {
         /* KDD Note:  Scaling by only first weight */
         if (hg->vwgt[i*hg->VtxWeightDim] == 0)
            hgp->vtx_scal[i] = 1.0;
         else            
            /* KDD Note:  Scaling by only first weight */
            hgp->vtx_scal[i] = 1. / hg->vwgt[i*hg->VtxWeightDim];
      }
  }

  ZOLTAN_FREE(&ldegree);

  return ZOLTAN_OK;
}
 
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

