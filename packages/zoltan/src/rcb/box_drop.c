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

/* recursive drop of an extended box thru an RCB tree */

#include <stdio.h>
#include <math.h>
#include "lb_const.h"
#include "rcb_const.h"

static void box_drop(struct rcb_tree *, struct rcb_box *, int *, int *,
   int, int);

void LB_box_drop(
LB             *lb,             /* The load-balancing structure with info for
                                   the RCB balancer. */
struct rcb_box *boxpt,          /* extended box */
int            *proc,           /* processors that box is in */
int            *procnum,        /* current number of processors on list */
int             proclower,      /* lower bound of partition */
int             procupper)      /* upper bound of partition */
{
   RCB_STRUCT *rcb;             /* Pointer to data structures for RCB.  */
   struct rcb_tree *treept;     /* tree of RCB cuts */

   rcb = (RCB_STRUCT *) (lb->Data_Structure);
   treept = rcb->Tree_Ptr;
   *procnum = 0;

   box_drop(treept, boxpt, proc, procnum, proclower, procupper);
}

static void box_drop(
struct rcb_tree *treept,        /* RCB tree */
struct rcb_box  *boxpt,         /* extended box */
int             *proc,          /* processors that box is in */
int             *procnum,       /* current number of processors on list */
int              proclower,     /* lower bound of partition */
int              procupper)     /* upper bound of partition */
{
  int       procmid;            /* 1st processor in upper half */
  int       dim;
  double    cut;

/* end recursion when partition size is a single proc */
/* add processor to list of processors */

  if (proclower == procupper) {
    proc[*procnum] = proclower;
    (*procnum)++;
    return;
  }

/* drop box on each side of cut if it extends beyond it */
/* important to use >= and <= criteria since either proc may own boundary */
/* procmid = 1st processor in upper half of partition, loc that stores cut
             for this partition in treept */
/* dim = dimension 0,1,2 of cut for this partition */
/* cut = position of cut for this partition */

  procmid = proclower + (procupper - proclower) / 2 + 1;
  dim = treept[procmid].dim;
  cut = treept[procmid].cut;

  if (boxpt->lo[dim] <= cut)
    box_drop(treept, boxpt, proc, procnum, proclower, procmid-1);
  if (boxpt->hi[dim] >= cut)
    box_drop(treept, boxpt, proc, procnum, procmid, procupper);
}
