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

#include <stdlib.h>
#include "zz_sort.h"    
#include "phg.h"
#include "zz_const.h"

int Zoltan_PHG_Vertex_Visit_Order(
  ZZ *zz, 
  HGraph *hg, 
  PHGPartParams *hgp, 
  int *order)
{
  int i, j, edge;
  int *ldegree=NULL, *gdegree=NULL; /* local/global degree */
  int *lpins=NULL, *gpins=NULL; /* local/global sum of pins */
  char *yo= "Zoltan_PHG_Vertex_Visit_Order";
  int ierr = ZOLTAN_OK;

  /* Start with linear order. */
  for (i=0; i<hg->nVtx; i++)
    order[i] = i;

  /* Permute order array according to chosen strategy. */
  switch (hgp->visit_order){
    case 0: 
      /* random node visit order (recommended)  */
      /* Synchronize so each proc in column visits in same order */
      Zoltan_Srand_Sync(Zoltan_Rand(NULL), &(hg->comm->RNGState_col),
                        hg->comm->col_comm);
      Zoltan_Rand_Perm_Int (order, hg->nVtx, &(hg->comm->RNGState_col));
      break;

    case 1: 
      /* linear (natural) vertex visit order */
      break;

    case 2:
    {
      /* increasing vertex weight */
      float *tmpvwgt;

      if (hg->VtxWeightDim == 1)
        tmpvwgt = hg->vwgt;
      else {
        /* Sort based on first component of multidimensional weight */
        tmpvwgt = (float *) ZOLTAN_MALLOC(hg->nVtx * sizeof(float));
        if (hg->nVtx && !tmpvwgt) MEMORY_ERROR;

        for (i = 0; i < hg->nVtx; i++)
          tmpvwgt[i] = hg->vwgt[i*hg->VtxWeightDim];
      } 
      
      Zoltan_quicksort_pointer_inc_float (order, tmpvwgt, 0, hg->nVtx-1);
      if (tmpvwgt != hg->vwgt) ZOLTAN_FREE(&tmpvwgt);
      break;
    }

    case 3: 
      /* increasing vertex degree */
      /* intentionally fall through into next case */
    case 4: 
      /* increasing vertex degree, weighted by # pins */

      /* allocate 4 arrays of size hg->nVtx with a single malloc */
      if (!(ldegree = (int *) ZOLTAN_MALLOC (4*sizeof(int) * hg->nVtx))){
        ZOLTAN_PRINT_WARN(zz->Proc, yo, "Out of memory");
        ZOLTAN_FREE (&ldegree);
        return ZOLTAN_MEMERR;
      }
      /* first local data, then global data */
      lpins = ldegree + hg->nVtx;
      gdegree = lpins + hg->nVtx;
      gpins = gdegree + hg->nVtx;

      /* loop over vertices */
      for (i=0; i<hg->nVtx; i++){
         ldegree[i] = hg->vindex[i+1] - hg->vindex[i]; /* local degree */
         lpins[i] = 0;
         /* loop over edges, sum up #pins */
         for (j= hg->vindex[i]; j < hg->vindex[i+1]; j++) {
           edge = hg->vedge[j];
           lpins[i] += hg->hindex[edge+1] - hg->hindex[edge];
         }
      }

      /* sum up local degrees in each column to get global degrees */
      /* also sum up #pins in same communication */
      MPI_Allreduce(ldegree, gdegree, 2*hg->nVtx, MPI_INT, MPI_SUM, 
         hg->comm->col_comm);

      /* sort by global values. same on every processor. */
      if (hgp->visit_order == 3)
        Zoltan_quicksort_pointer_inc_int_int (order, gdegree, gpins,
          0, hg->nVtx-1);
      else /* hgp->visit_order == 4 */
        Zoltan_quicksort_pointer_inc_int_int (order, gpins, gdegree,
          0, hg->nVtx-1);

      ZOLTAN_FREE (&ldegree);
      break;

    /* add more cases here */
  }
End:
  return ierr;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
