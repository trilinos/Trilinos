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
#include "rib_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_RIB_Build_Structure(LB *lb, int *num_obj, int *max_obj, int wgtflag)
{
/* Function to build the geometry-based data structures for RIB method. */
char           *yo = "LB_RIB_Build_Structure";
RIB_STRUCT     *rib;                  /* Data structure for RIB.             */
int            ierr = 0;

  /* Allocate an RIB data structure for this load balancing structure.
     If the previous data structure is still there, free the Dots and IDs first;
     the other fields can be reused. */

  if (lb->Data_Structure == NULL) {
    rib = (RIB_STRUCT *) LB_MALLOC(sizeof(RIB_STRUCT));
    if (rib == NULL) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      return(LB_MEMERR);
    }
    lb->Data_Structure = (void *) rib;
    rib->Tree_Ptr = NULL;
    rib->Global_IDs = NULL;
    rib->Local_IDs = NULL;
    rib->Dots = NULL;

    rib->Tree_Ptr = (struct rib_tree *)
                    LB_MALLOC(lb->Num_Proc* sizeof(struct rib_tree));
    if (rib->Tree_Ptr == NULL) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      LB_RIB_Free_Structure(lb);
      return(LB_MEMERR);
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
                               wgtflag);
  if (ierr) {
    LB_PRINT_ERROR(lb->Proc, yo, "Error returned from LB_RB_Build_Structure.");
    LB_RIB_Free_Structure(lb);
    return(ierr);
  }

  return(LB_OK);
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
