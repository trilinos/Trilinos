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
#include "phg.h"
#include "hg_hypergraph.h"



/****************************************************************************/


/* Scaling the weight of hypergraph edges. 
   This changes the inner product used in matching,
   hopefully to the better! 
   Note that the scaled weights are only used for matching,
   and the original weights are restored afterwards (see phg_match.c).

   EBEB: Removed Robert's serial scaling methods. 
         We should look at these later.
 */
int Zoltan_PHG_Scale_Weights (ZZ *zz, HGraph *hg, float *new_ewgt, 
                             PHGPartParams *hgp)
{
int    i;
int    *lsize = NULL;  /* local edge sizes */
int    *size = NULL;   /* edge sizes */
static char *yo = "Zoltan_PHG_Scale_Weights";

  /* allocate size arrays */
  if (!(size  = (int *) ZOLTAN_MALLOC (sizeof(int) * hg->nEdge)) ||
      !(lsize = (int *) ZOLTAN_MALLOC (sizeof(int) * hg->nEdge)) ){
    ZOLTAN_PRINT_WARN(zz->Proc, yo, "Out of memory");
    return ZOLTAN_MEMERR;
  }

  switch (hgp->ews){

  case 0:
    /* do nothing */
    break;

  case 1:
    /* absorption or clique scaling; scale by #pins -1 */
    /* first compute size of all hyperedges */
    for (i = 0; i < hg->nEdge; i++) {
      lsize[i] = hg->hindex[i+1] - hg->hindex[i];
    }

    /* sum up local sizes */
    /* assume SCHEMEA : all procs in a row have same # hyperedges */
    MPI_Allreduce(lsize, size, hg->nEdge, MPI_INT, MPI_SUM, 
                  hg->comm->row_comm);

    /* scale edge weights */
    for (i = 0; i < hg->nEdge; i++) {
      /* EBEB printf("[%1d] Debug: Hyperedge %d has size %d\n", 
         zz->Proc, i, size[i]); */
      if (size[i]>1)
        new_ewgt[i] = (hg->ewgt ? hg->ewgt[i] : 1.0) / (size[i]-1.0);
      else
        new_ewgt[i] = 0.0;
    }
    break;
  }

  ZOLTAN_FREE(&size);
  ZOLTAN_FREE(&lsize);

  return ZOLTAN_OK;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

