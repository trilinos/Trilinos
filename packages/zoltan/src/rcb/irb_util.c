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
#include "irb_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_IRB_Build_Structure(LB *lb, int *num_obj, int *max_obj, int wgtflag)
{
/* Function to build the geometry-based data structures for IRB method. */
char           *yo = "LB_IRB_Build_Structure";
IRB_STRUCT     *irb;                  /* Data structure for IRB.             */
int            ierr = 0;

  /* Allocate an IRB data structure for this load balancing structure.
     If the previous data structure is still there, free the Dots and IDs first;
     the other fields can be reused. */

  if (lb->Data_Structure == NULL) {
    irb = (IRB_STRUCT *) LB_MALLOC(sizeof(IRB_STRUCT));
    if (irb == NULL) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      return(LB_MEMERR);
    }
    lb->Data_Structure = (void *) irb;
    irb->Tree_Ptr = NULL;
    irb->Global_IDs = NULL;
    irb->Local_IDs = NULL;
    irb->Dots = NULL;

    irb->Tree_Ptr = (struct irb_tree *)
                    LB_MALLOC(lb->Num_Proc* sizeof(struct irb_tree));
    if (irb->Tree_Ptr == NULL) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      LB_IRB_Free_Structure(lb);
      return(LB_MEMERR);
    }
  }
  else {
    irb = (IRB_STRUCT *) lb->Data_Structure;
    LB_FREE(&(irb->Global_IDs));
    LB_FREE(&(irb->Local_IDs));
    LB_FREE(&(irb->Dots));
  }

  ierr = LB_RB_Build_Structure(lb, &(irb->Global_IDs), &(irb->Local_IDs),
                               &(irb->Dots), num_obj, max_obj, &(irb->Num_Geom),
                               wgtflag);
  if (ierr) {
    LB_PRINT_ERROR(lb->Proc, yo, "Error returned from LB_RB_Build_Structure.");
    LB_IRB_Free_Structure(lb);
    return(ierr);
  }

  return(LB_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void LB_IRB_Free_Structure(LB *lb)
{
/* Deallocate the persistent IRB data structures in lb->Structure.  */
IRB_STRUCT    *irb;                   /* Data structure for IRB. */

  irb = (IRB_STRUCT *) lb->Data_Structure;

  if (irb != NULL) {
    LB_FREE(&(irb->Tree_Ptr));
    LB_FREE(&(irb->Global_IDs));
    LB_FREE(&(irb->Local_IDs));
    LB_FREE(&(irb->Dots));
    LB_FREE(&(lb->Data_Structure));
  }
}
