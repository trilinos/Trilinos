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



static ZOLTAN_HG_GROUPING_FN grouping_no; /* template --- no grouping */

/****************************************************************************/



int Zoltan_HG_Set_Grouping_Fn(HGPartParams *hgp)
{
int found = 1;

  
  if (!strcasecmp(hgp->redm_str, "no")) hgp->grouping = NULL;
  else {
    found = 0; hgp->grouping = NULL;
  }

  if (hgp->grouping) {
    /* If reduction method is a grouping, set the improvement and edge weight
     * scaling functions accordingly */

    /* register optimization function */
    hgp->grouping_opt=NULL;
    }
  return found;
}

/****************************************************************************/



int Zoltan_HG_Grouping (ZZ *zz, HGraph *hg, Packing pack, HGPartParams *hgp,
 int *limit)
{
  int   err = ZOLTAN_OK;
  float *old_ewgt = NULL, *new_ewgt = NULL;
  char  *yo = "Zoltan_HG_Grouping";

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Scale the weight of the edges */
  if (hg->vwgt && hgp->ews) {
    if (!(new_ewgt = (float*) ZOLTAN_MALLOC (hg->nEdge * sizeof(float)))) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
    Zoltan_HG_Scale_HGraph_Weight (zz, hg, new_ewgt, hgp->ews);
    old_ewgt = hg->ewgt;
    hg->ewgt = new_ewgt;
  }

  /* Do the grouping */
  if (hgp->grouping) {
    err = hgp->grouping(zz, hg, pack, limit);
    if (err != ZOLTAN_OK && err != ZOLTAN_WARN)
      goto End;
  }

  /* Optimization */
  if (hgp->grouping_opt != NULL)
    err = hgp->grouping_opt (zz, hg, pack, limit);

End:
  /* Restore the old edge weights */
  if (hg->vwgt && hgp->ews)
    hg->ewgt = old_ewgt;

  ZOLTAN_FREE ((void**) &new_ewgt);
  ZOLTAN_TRACE_EXIT(zz, yo);
  return err;
}



/******************************************************************************/
/* null grouping, hypergraph version */
static int grouping_no(ZZ *zz, HGraph *hg, Matching match, int *limit)
{
  char *yo = "grouping_no";
  return ZOLTAN_OK;
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
