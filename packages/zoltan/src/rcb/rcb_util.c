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


#include "zz_const.h"
#include "rcb.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_RCB_Build_Structure(ZZ *zz, int *num_obj, int *max_obj, int wgtflag,
                           int use_ids)
{
/*
 *  Function to build the geometry-based data structures for 
 *  Steve Plimpton's RCB implementation.
 */
char *yo = "Zoltan_RCB_Build_Structure";
RCB_STRUCT *rcb;                      /* Data structure for RCB.             */
struct rcb_tree *treeptr;
int i, ierr = 0;

  /*
   * Allocate an RCB data structure for this Zoltan structure.
   * If the previous data structure is still there, free the Dots and IDs first;
   * the other fields can be reused.
   */

  if (zz->LB.Data_Structure == NULL) {
    rcb = (RCB_STRUCT *) ZOLTAN_MALLOC(sizeof(RCB_STRUCT));
    if (rcb == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return(ZOLTAN_MEMERR);
    }
    zz->LB.Data_Structure = (void *) rcb;
    rcb->Tree_Ptr = NULL;
    rcb->Box = NULL;
    rcb->Global_IDs = NULL;
    rcb->Local_IDs = NULL;
    rcb->Dots = NULL;

    rcb->Tree_Ptr = (struct rcb_tree *)
      ZOLTAN_MALLOC(zz->LB.Num_Global_Parts * sizeof(struct rcb_tree));
    rcb->Box = (struct rcb_box *) ZOLTAN_MALLOC(sizeof(struct rcb_box));
    if (rcb->Tree_Ptr == NULL || rcb->Box == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      Zoltan_RCB_Free_Structure(zz);
      return(ZOLTAN_MEMERR);
    }
    /* initialize Tree_Ptr */
    for (i = 0; i < zz->LB.Num_Global_Parts; i++) {
       treeptr = &(rcb->Tree_Ptr[i]);
       /* initialize dim to -1 to prevent use of cut */
       treeptr->dim = -1;
       treeptr->cut = 0.0;
       treeptr->parent = treeptr->left_leaf = treeptr->right_leaf = 0;
    }
  }
  else {
    rcb = (RCB_STRUCT *) zz->LB.Data_Structure;
    ZOLTAN_FREE(&(rcb->Global_IDs));
    ZOLTAN_FREE(&(rcb->Local_IDs));
    ZOLTAN_FREE(&(rcb->Dots));
  }

  ierr = Zoltan_RB_Build_Structure(zz, &(rcb->Global_IDs), &(rcb->Local_IDs),
                               &(rcb->Dots), num_obj, max_obj,
                               &(rcb->Num_Dim),
                               wgtflag, use_ids);
  if (ierr) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                       "Error returned from Zoltan_RB_Build_Structure.");
    Zoltan_RCB_Free_Structure(zz);
    return(ierr);
  }

  return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void Zoltan_RCB_Free_Structure(ZZ *zz)
{
/*
 * Deallocate the persistent RCB data structures in zz->Structure.
 */
RCB_STRUCT *rcb;                      /* Data structure for RCB.             */

  rcb = (RCB_STRUCT *) (zz->LB.Data_Structure);

  if (rcb != NULL) {
    ZOLTAN_FREE(&(rcb->Tree_Ptr));
    ZOLTAN_FREE(&(rcb->Box));
    ZOLTAN_FREE(&(rcb->Global_IDs));
    ZOLTAN_FREE(&(rcb->Local_IDs));
    ZOLTAN_FREE(&(rcb->Dots));
    ZOLTAN_FREE(&(zz->LB.Data_Structure));
  }
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
