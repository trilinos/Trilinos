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

#include <stdlib.h>
#include "phg.h"

int Zoltan_PHG_Vertex_Visit_Order(
  ZZ *zz, 
  HGraph *hg, 
  PHGPartParams *hgp, 
  int *order)
{
  int i;

  /* Start with linear order. */
  for (i=0; i<hg->nVtx; i++)
    order[i] = i;

  /* Permute order array according to chosen strategy. */
  switch (hgp->visit_order){
    case 1: 
      /* random node visit order */
      Zoltan_Rand_Perm_Int (order, hg->nVtx);
      break;
    /* add more cases here */
  }

  return ZOLTAN_OK;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
