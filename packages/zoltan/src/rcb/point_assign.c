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

#include <stdio.h>
#include "lb_const.h"
#include "rcb_const.h"


int LB_Point_Assign(
LB       *lb,                   /* The load-balancing structure */
double   *coords,               /* vector of point coordinates */
int      *proc)			/* processor that point lands in */
{
   int       proclower;         /* lower bound of partition */
   int       procupper;         /* upper bound of partition */
   int       procmid;           /* 1st processor in upper half */
   RCB_STRUCT *rcb;             /* Pointer to data structures for RCB.  */
   struct rcb_tree *treept;     /* tree of RCB cuts */

   if (lb->Method != RCB) {
     fprintf(stderr, "ERROR: Point_Assign only valid when method is RCB\n");
     *proc = -1;
     return(LB_FATAL);
   }
   if (lb->Data_Structure == NULL) {
     fprintf(stderr, "ERROR: No LB_Data_Structure available for Point_Assign\n"); 
     *proc = -1;
     return(LB_FATAL);
   }


   rcb = (RCB_STRUCT *) (lb->Data_Structure);
   treept = rcb->Tree_Ptr;
   if (treept[0].dim < 0) {	/* RCB tree was never created. */
     fprintf(stderr, "ERROR: No RCB tree not saved for Point_Assign.\n"); 
     fprintf(stderr, "       Must set parameter KEEP_CUTS to 1.\n"); 
     *proc = -1;
     return(LB_FATAL);
   }

   proclower = 0;
   procupper = lb->Num_Proc - 1;

   while (proclower < procupper) {
      procmid = proclower + (procupper - proclower) / 2 + 1;
      if (coords[treept[procmid].dim] <= treept[procmid].cut)
         procupper = procmid-1;
      else
         proclower = procmid;
   }

   *proc = proclower;

   return(LB_OK);
}
