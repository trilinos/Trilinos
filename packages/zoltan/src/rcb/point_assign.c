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

#include <stdio.h>
#include "lb_const.h"
#include "rcb_const.h"
#include "irb_const.h"


int LB_Point_Assign(
LB       *lb,                   /* The load-balancing structure */
double   *coords,               /* vector of point coordinates */
int      *proc)                 /* processor that point lands in */
{
/* Locate which processor a point is inside within the tree defined
   by the recursive bisection algorithm chosen. */

     char             *yo = "LB_Point_Assign";
     int               procmid; /* 1st processor in upper half */
     RCB_STRUCT        *rcb;    /* Pointer to data structures for RCB.  */
     struct rcb_tree   *treept; /* tree of RCB cuts */
     IRB_STRUCT        *irb;    /* Pointer to data structures for IRB. */
     struct irb_tree   *itree;  /* tree of IRB cuts */

     if (lb->Data_Structure == NULL) {
        LB_PRINT_ERROR(-1, yo, 
                       "No LB_Data_Structure available for LB_Point_Assign");
        *proc = -1;
        return(LB_FATAL);
     }

     if (lb->Method == RCB) {
        rcb = (RCB_STRUCT *) (lb->Data_Structure);
        treept = rcb->Tree_Ptr;
        if (treept[0].dim < 0) { /* RCB tree was never created. */
           LB_PRINT_ERROR(lb->Proc, yo, "No RCB tree saved for Point_Assign; "
                                        "Must set parameter KEEP_CUTS to 1.");
           *proc = -1;
           return(LB_FATAL);
        }

        procmid = treept[0].right_leaf;

        while (procmid > 0)
           if (coords[treept[procmid].dim] <= treept[procmid].cut)
              procmid = treept[procmid].left_leaf;
           else
              procmid = treept[procmid].right_leaf;

        *proc = -procmid;

        return(LB_OK);
     }
     else if (lb->Method == IRB) {
        irb = (IRB_STRUCT *) (lb->Data_Structure);
        itree = irb->Tree_Ptr;
        if ((procmid = itree[0].right_leaf) < 0) { /* IRB tree never created */
           LB_PRINT_ERROR(lb->Proc, yo, "No IRB tree saved for Point_Assign; "
                                     "Must set parameter IRB_KEEP_CUTS to 1.");
           *proc = -1;
           return(LB_FATAL);
        }

        switch (irb->Num_Geom) {
           case 3:
              while (procmid > 0)
                 if(((coords[0] - itree[procmid].cm[0])*itree[procmid].ev[0] +
                     (coords[1] - itree[procmid].cm[1])*itree[procmid].ev[1] +
                     (coords[2] - itree[procmid].cm[2])*itree[procmid].ev[2])
                       <= itree[procmid].cut)
                    procmid = itree[procmid].left_leaf;
                 else
                    procmid = itree[procmid].right_leaf;
              break;
           case 2:
              while (procmid > 0)
                 if(((coords[0] - itree[procmid].cm[0])*itree[procmid].ev[0] +
                     (coords[1] - itree[procmid].cm[1])*itree[procmid].ev[1])
                       <= itree[procmid].cut)
                    procmid = itree[procmid].left_leaf;
                 else
                    procmid = itree[procmid].right_leaf;
              break;
           case 1:
              while (procmid > 0)
                 if (coords[0] <= itree[procmid].cut)
                    procmid = itree[procmid].left_leaf;
                 else
                    procmid = itree[procmid].right_leaf;
              break;
        }

        *proc = -procmid;

        return(LB_OK);
     }
     else {
        LB_PRINT_ERROR(lb->Proc, yo, "LB_Point_Assign valid only when method "
                                     "is RCB or IRB.");
        *proc = -1;
        return(LB_FATAL);
     }
}
