/*==================================================================== 
 * ------------------------
 * | CVS File Information |
 * ------------------------
 * 
 * $RCSfile$ 
 * 
 * $Author$ 
 * 
 * $Date$ 
 * 
 * $Revision$ 
 * 
 *====================================================================*/ 

/* Locate which rcb_box a point is inside within the rcb tree. */

#include "lb_const.h"
#include "rcb_const.h"

int LB_point_drop(
LB       *lb,                   /* The load-balancing structure with info for
                                   the RCB balancer. */
double   *coord,                /* array of node coordinates */
int       proclower,            /* lower bound of partition */
int       procupper)            /* upper bound of partition */
{
   int       procmid;           /* 1st processor in upper half */
   RCB_STRUCT *rcb;             /* Pointer to data structures for RCB.  */
   struct rcb_tree *treept;     /* tree of RCB cuts */

   rcb = (RCB_STRUCT *) (lb->Data_Structure);
   treept = rcb->Tree_Ptr;

   while (proclower < procupper) {
      procmid = proclower + (procupper - proclower) / 2 + 1;
      if (coord[treept[procmid].dim] <= treept[procmid].cut)
         procupper = procmid-1;
      else
         proclower = procmid;
   }

   return proclower;
}
