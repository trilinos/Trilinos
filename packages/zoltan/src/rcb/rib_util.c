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
#include "rib_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_RIB_Build_Structure(LB *lb, int *num_obj, int *max_obj, int wgtflag,
                           int use_ids)
{
/* Function to build the geometry-based data structures for RIB method. */
char           *yo = "LB_RIB_Build_Structure";
RIB_STRUCT     *rib;                  /* Data structure for RIB.             */
struct rib_tree *treeptr;
int            i, ierr = 0;

  /* Allocate an RIB data structure for this load balancing structure.
     If the previous data structure is still there, free the Dots and IDs first;
     the other fields can be reused. */

  if (lb->Data_Structure == NULL) {
    rib = (RIB_STRUCT *) LB_MALLOC(sizeof(RIB_STRUCT));
    if (rib == NULL) {
      ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      return(ZOLTAN_MEMERR);
    }
    lb->Data_Structure = (void *) rib;
    rib->Tree_Ptr = NULL;
    rib->Global_IDs = NULL;
    rib->Local_IDs = NULL;
    rib->Dots = NULL;

    rib->Tree_Ptr = (struct rib_tree *)
                    LB_MALLOC(lb->Num_Proc* sizeof(struct rib_tree));
    if (rib->Tree_Ptr == NULL) {
      ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      LB_RIB_Free_Structure(lb);
      return(ZOLTAN_MEMERR);
    }
    /* initialize Tree_Ptr */
    for (i = 0; i < lb->Num_Proc; i++) {
      treeptr = &(rib->Tree_Ptr[i]);
      treeptr->cm[0] = treeptr->cm[1] = treeptr->cm[2] = 0.0;
      treeptr->ev[0] = treeptr->ev[1] = treeptr->ev[2] = 0.0;
      treeptr->cut = 0.0;
      treeptr->parent = treeptr->left_leaf = treeptr->right_leaf = 0;
    }
  }
  else {
    rib = (RIB_STRUCT *) lb->Data_Structure;
    LB_FREE(&(rib->Global_IDs));
    LB_FREE(&(rib->Local_IDs));
    LB_FREE(&(rib->Dots));
  }

  ierr = LB_RB_Build_Structure(lb, &(rib->Global_IDs), &(rib->Local_IDs),
                               &(rib->Dots), num_obj, max_obj, &(rib->Num_Geom),
                               wgtflag, use_ids);
  if (ierr) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Error returned from LB_RB_Build_Structure.");
    LB_RIB_Free_Structure(lb);
    return(ierr);
  }

  return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void LB_RIB_Free_Structure(LB *lb)
{
/* Deallocate the persistent RIB data structures in lb->Structure.  */
RIB_STRUCT    *rib;                   /* Data structure for RIB. */

  rib = (RIB_STRUCT *) lb->Data_Structure;

  if (rib != NULL) {
    LB_FREE(&(rib->Tree_Ptr));
    LB_FREE(&(rib->Global_IDs));
    LB_FREE(&(rib->Local_IDs));
    LB_FREE(&(rib->Dots));
    LB_FREE(&(lb->Data_Structure));
  }
}
