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

#include "phypergraph.h"



static ZOLTAN_PHG_MATCHING_FN matching_no;   /* template -- matching */

static double sim (PHGraph*, int, int);

/* static void check_upper_bound_of_matching_weight (Graph*, ZZ*, Matching); */
/* static int graph_connected_components (int, int*, int*, int);             */



/*****************************************************************************/
int Zoltan_PHG_Set_Matching_Fn (PHGPartParams *hgp)
{
  int found = 1;

  if (!strcasecmp(hgp->redm_str, "no"))   hgp->matching = matching_no ;
  else  {
    found = 0;
    hgp->matching = NULL;
  }

  if (hgp->matching) {
     /* If reduction method is a matching, set the improvement and edge weight
        scaling functions accordingly. */

     /* Note: matching_aug1 is identical to matching_mxm -> it was eliminated */

     hgp->matching_opt=NULL;
     }
  return found;
}



/****************************************************************************/
/* This is the similarity measure between two vertices in a hypergraph.
   The similarity is equal to the scaled weight of the edge in the
   transformed graph. But with this function we calculate the edge
   weights on the fly without explicitly constructing the graph. */

static double sim (PHGraph *hg, int a, int b)
{
int    i, j, edge, pins, end;
double  weight, sim = 0.0;

  /* First calculate the edge weight of the graph between a and b */
  for (i = hg->vindex[a]; i < hg->vindex[a+1]; i++) {
     edge = hg->vedge[i];
     end  = hg->hindex[edge+1];
     j    = hg->hindex[edge];
     while (j < end && hg->hvertex[j] != b)
        j++;
     if (j < end) {
        pins = end - hg->hindex[edge];
        weight = 2.0 / ((pins-1) * pins);
        if (hg->ewgt)
           weight *= hg->ewgt[edge];
        sim += weight;
        }
     }
  return sim;
}



/*****************************************************************************/

int Zoltan_PHG_Matching (
  ZZ *zz,
  PHGraph *hg,
  Matching match,
  PHGPartParams *hgp,
  int *limit)
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
/* template for matching, hypergraph version */
static int matching_no (ZZ *zz, PHGraph *hg, Matching match, int *limit)
{
  return ZOLTAN_OK;
}



/*****************************************************************************/
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
