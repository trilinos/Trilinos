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



static ZOLTAN_PHG_MATCHING_FN matching_no;   /* template -- matching */


/* static void check_upper_bound_of_matching_weight (Graph*, ZZ*, Matching); */
/* static int graph_connected_components (int, int*, int*, int);             */



/*****************************************************************************/
int Zoltan_PHG_Set_Matching_Fn (PHGPartParams *hgp)
{
    int exist = 1;
    if (!strcasecmp(hgp->redm_str, "no"))   hgp->matching = NULL;
    else {
        exist = 0;
        hgp->matching = NULL;
    }

  return exist;
}



/* Removed sim() from serial version at this point */

/*****************************************************************************/

int Zoltan_PHG_Matching (
  ZZ *zz,
  PHGraph *hg,
  Matching match,
  PHGPartParams *hgp)
{
float *old_ewgt = NULL, *new_ewgt = NULL;
int   err;
char  *yo = "Zoltan_PHG_Matching";

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Scale the weight of the edges */
  if (hg->vwgt && hgp->ews) {
     if (!(new_ewgt = (float*) ZOLTAN_MALLOC (hg->nEdge * sizeof(float)))) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        err = ZOLTAN_MEMERR;
        goto End;
        }
     Zoltan_PHG_Scale_HGraph_Weight (zz, hg, new_ewgt, hgp->ews);
     old_ewgt = hg->ewgt;
     hg->ewgt = new_ewgt;
     }

  /* Do the matching */
  if (hgp->matching) {
     err = hgp->matching (zz, hg, match);
     if (err != ZOLTAN_OK && err != ZOLTAN_WARN)
        goto End;
     }

  /* Removed serial Optimization here */

End:

  /* Restore the old edge weights */
  if (hg->vwgt && hgp->ews)
      hg->ewgt = old_ewgt;

  ZOLTAN_FREE ((void**) &new_ewgt);
  ZOLTAN_TRACE_EXIT(zz, yo);
  return err;
}



/*****************************************************************************/
/* template for matching, hypergraph version */
static int matching_no (ZZ *zz, PHGraph *hg, Matching match)
{
  return ZOLTAN_OK;
}



/*****************************************************************************/
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
