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

#include "lbi_const.h"
#include "lb_const.h"
#include "rcb_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_RCB_Build_Structure(LB *lb, int *num_obj, int *max_obj, int wgtflag)
{
/*
 *  Function to build the geometry-based data structures for 
 *  Steve Plimpton's RCB implementation.
 */
char *yo = "LB_RCB_Build_Structure";
RCB_STRUCT *rcb;                      /* Data structure for RCB.             */
int i, ierr = 0;
int num_geom;

  /*
   * Allocate an RCB data structure for this load balancing structure.
   * If the previous data structure is still there, free the Dots and IDs first;
   * the other fields can be reused.
   */

  if (lb->Data_Structure == NULL) {
    rcb = (RCB_STRUCT *) LB_MALLOC(sizeof(RCB_STRUCT));
    if (rcb == NULL) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      return(LB_MEMERR);
    }
    lb->Data_Structure = (void *) rcb;
    rcb->Tree_Ptr = NULL;
    rcb->Box = NULL;
    rcb->Global_IDs = NULL;
    rcb->Local_IDs = NULL;
    rcb->Dots = NULL;

    rcb->Tree_Ptr = (struct rcb_tree *)
      LB_MALLOC(lb->Num_Proc* sizeof(struct rcb_tree));
    rcb->Box = (struct rcb_box *) LB_MALLOC(sizeof(struct rcb_box));
    if (rcb->Tree_Ptr == NULL || rcb->Box == NULL) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      LB_RCB_Free_Structure(lb);
      return(LB_MEMERR);
    }
    /* initialize dim to -1 to prevent use of cut */
    for (i = 0; i < lb->Num_Proc; i++)
       rcb->Tree_Ptr[i].dim = -1;
  }
  else {
    rcb = (RCB_STRUCT *) lb->Data_Structure;
    LB_FREE(&(rcb->Global_IDs));
    LB_FREE(&(rcb->Local_IDs));
    LB_FREE(&(rcb->Dots));
  }

  ierr = LB_RB_Build_Structure(lb, &(rcb->Global_IDs), &(rcb->Local_IDs),
                               &(rcb->Dots), num_obj, max_obj, &num_geom,
                               wgtflag);
  if (ierr) {
    LB_PRINT_ERROR(lb->Proc, yo, "Error returned from LB_RB_Build_Structure.");
    LB_RCB_Free_Structure(lb);
    return(ierr);
  }

  return(LB_OK);
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
    LB_FREE(&(rcb->Tree_Ptr));
    LB_FREE(&(rcb->Box));
    LB_FREE(&(rcb->Global_IDs));
    LB_FREE(&(rcb->Local_IDs));
    LB_FREE(&(rcb->Dots));
    LB_FREE(&(lb->Data_Structure));
  }
}
