/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include "lb_const.h"
#include "rcb_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_RCB_Build_Structure(LB *lb, int *num_obj, int *max_obj, int wgtflag,
                           int use_ids)
{
/*
 *  Function to build the geometry-based data structures for 
 *  Steve Plimpton's RCB implementation.
 */
char *yo = "LB_RCB_Build_Structure";
RCB_STRUCT *rcb;                      /* Data structure for RCB.             */
struct rcb_tree *treeptr;
int i, ierr = 0;
int num_geom;

  /*
   * Allocate an RCB data structure for this load balancing structure.
   * If the previous data structure is still there, free the Dots and IDs first;
   * the other fields can be reused.
   */

  if (lb->Data_Structure == NULL) {
    rcb = (RCB_STRUCT *) ZOLTAN_MALLOC(sizeof(RCB_STRUCT));
    if (rcb == NULL) {
      ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      return(ZOLTAN_MEMERR);
    }
    lb->Data_Structure = (void *) rcb;
    rcb->Tree_Ptr = NULL;
    rcb->Box = NULL;
    rcb->Global_IDs = NULL;
    rcb->Local_IDs = NULL;
    rcb->Dots = NULL;

    rcb->Tree_Ptr = (struct rcb_tree *)
      ZOLTAN_MALLOC(lb->Num_Proc* sizeof(struct rcb_tree));
    rcb->Box = (struct rcb_box *) ZOLTAN_MALLOC(sizeof(struct rcb_box));
    if (rcb->Tree_Ptr == NULL || rcb->Box == NULL) {
      ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      LB_RCB_Free_Structure(lb);
      return(ZOLTAN_MEMERR);
    }
    /* initialize Tree_Ptr */
    for (i = 0; i < lb->Num_Proc; i++) {
       treeptr = &(rcb->Tree_Ptr[i]);
       /* initialize dim to -1 to prevent use of cut */
       treeptr->dim = -1;
       treeptr->cut = 0.0;
       treeptr->parent = treeptr->left_leaf = treeptr->right_leaf = 0;
    }
  }
  else {
    rcb = (RCB_STRUCT *) lb->Data_Structure;
    ZOLTAN_FREE(&(rcb->Global_IDs));
    ZOLTAN_FREE(&(rcb->Local_IDs));
    ZOLTAN_FREE(&(rcb->Dots));
  }

  ierr = LB_RB_Build_Structure(lb, &(rcb->Global_IDs), &(rcb->Local_IDs),
                               &(rcb->Dots), num_obj, max_obj, &num_geom,
                               wgtflag, use_ids);
  if (ierr) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Error returned from LB_RB_Build_Structure.");
    LB_RCB_Free_Structure(lb);
    return(ierr);
  }

  return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void LB_RCB_Free_Structure(LB *lb)
{
/*
 * Deallocate the persistent RCB data structures in lb->Structure.
 */
RCB_STRUCT *rcb;                      /* Data structure for RCB.             */

  rcb = (RCB_STRUCT *) lb->Data_Structure;

  if (rcb != NULL) {
    ZOLTAN_FREE(&(rcb->Tree_Ptr));
    ZOLTAN_FREE(&(rcb->Box));
    ZOLTAN_FREE(&(rcb->Global_IDs));
    ZOLTAN_FREE(&(rcb->Local_IDs));
    ZOLTAN_FREE(&(rcb->Dots));
    ZOLTAN_FREE(&(lb->Data_Structure));
  }
}
