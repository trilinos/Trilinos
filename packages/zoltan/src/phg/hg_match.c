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

#include "hypergraph.h"



static ZOLTAN_HG_MATCHING_FN matching_no;  /* template --- matching */


/*****************************************************************************/

int Zoltan_HG_Set_Matching_Fn(HGPartParams *hgp)
{
  int found = 1;

  if      (!strcasecmp(hgp->redm_str, "null")) hgp->matching = matching_no;
  else if (!strcasecmp(hgp->redm_str, "no"))   hgp->matching = NULL;
  else {
    found = 0;
    hgp->matching = NULL;
  }

  if (hgp->matching) {
    /* If reduction method is a matching, set the improvement and edge weight
       scaling functions accordingly. */
       
    hgp->matching_opt=NULL;
  }
  return found;
}



/*****************************************************************************/

int Zoltan_HG_Matching (
  ZZ *zz,
  HGraph *hg,
  Matching match,
  HGPartParams *hgp,
  int *limit)
{
  float *old_ewgt = NULL, *new_ewgt = NULL;
  int   err;
  char  *yo = "Zoltan_HG_Matching";

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Scale the weight of the edges */
  if (hg->vwgt && hgp->ews) {
    if (!(new_ewgt = (float*) ZOLTAN_MALLOC (hg->nEdge * sizeof(float)))) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      err = ZOLTAN_MEMERR;
      goto End;
    }
    Zoltan_HG_Scale_HGraph_Weight (zz, hg, new_ewgt, hgp->ews);
    old_ewgt = hg->ewgt;
    hg->ewgt = new_ewgt;
  }

  /* Do the matching */
  if (hgp->matching) {
    err = hgp->matching (zz, hg, match, limit);
    if (err != ZOLTAN_OK && err != ZOLTAN_WARN)
      goto End;
  }

  /* Optimization */
  if (hgp->matching_opt)
    err = hgp->matching_opt (zz, hg, match, limit);
     
End:

  /* Restore the old edge weights */
  if (hg->vwgt && hgp->ews)
    hg->ewgt = old_ewgt;

  ZOLTAN_FREE ((void**) &new_ewgt);
  ZOLTAN_TRACE_EXIT(zz, yo);
  return err;
}



/*****************************************************************************/

/* null matching, hypergraph version */
static int matching_no(ZZ *zz, HGraph *hg, Matching match, int *limit)
{
  char *yo = "matching_no";
  return ZOLTAN_OK;
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
