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
#ifndef lint
static char *cvs_box_assign_c = "$Id$";
#endif

/* Determine which processors a box intersects. */
/* Currently assumes that partitioning has used RCB, */
/* but should be modified to support IRB and to return */
/* an error message if other method was used */

#include <stdio.h>
#include <math.h>
#include "lb_const.h"
#include "rcb_const.h"

static void Box_Assign(struct rcb_tree *, struct rcb_box *, int *, int *,
   int, int);

int LB_Box_Assign(
LB             *lb,             /* The load-balancing structure */
double          xmin,		/* lower x extent of box */
double          ymin,		/* lower y extent of box */
double          zmin,		/* lower z extent of box */
double          xmax,		/* upper x extent of box */
double          ymax,		/* upper y extent of box */
double          zmax,		/* upper z extent of box */
int            *procs,          /* list of procs that box intersects */
int            *numprocs)       /* number of processors in proc list */
{
   RCB_STRUCT *rcb;             /* Pointer to data structures for RCB.  */
   struct rcb_tree *treept;     /* tree of RCB cuts */
   struct rcb_box  box;         /* box data structure */

   if (lb->Method != RCB) {
     fprintf(stderr, "ERROR: Box_Assign only valid when method is RCB\n");
     *procs = -1;
     *numprocs = 0;
     return(LB_FATAL);
   }
   if (lb->Data_Structure == NULL) {
     fprintf(stderr, "ERROR: No LB_Data_Structure available for Box_Assign\n"); 
     *procs = -1;
     *numprocs = 0;
     return(LB_FATAL);
   }

   rcb = (RCB_STRUCT *) (lb->Data_Structure);
   treept = rcb->Tree_Ptr;
   if (treept[0].dim < 0) {     /* RCB tree was never created. */
     fprintf(stderr, "ERROR: No RCB tree not saved for Box_Assign.\n");
     fprintf(stderr, "       Must set parameter KEEP_CUTS to 1.\n");
     *procs = -1;
     *numprocs = 0;
     return(LB_FATAL);
   }

   *numprocs = 0;
   box.lo[0] = xmin;
   box.lo[1] = ymin;
   box.lo[2] = zmin;
   box.hi[0] = xmax;
   box.hi[1] = ymax;
   box.hi[2] = zmax;

   Box_Assign(treept, &box, procs, numprocs, 0, lb->Num_Proc - 1);

   return(LB_OK);
}

static void Box_Assign(
struct rcb_tree *treept,        /* RCB tree */
struct rcb_box  *boxpt,         /* extended box */
int             *procs,         /* processors that box is in */
int             *numprocs,       /* current number of processors on list */
int              proclower,     /* lower bound of partition */
int              procupper)     /* upper bound of partition */
{
  int       procmid;            /* 1st processor in upper half */
  int       dim;
  double    cut;

/* end recursion when partition size is a single proc */
/* add processor to list of processors */

  if (proclower == procupper) {
    procs[*numprocs] = proclower;
    (*numprocs)++;
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
    Box_Assign(treept, boxpt, procs, numprocs, proclower, procmid-1);
  if (boxpt->hi[dim] >= cut)
    Box_Assign(treept, boxpt, procs, numprocs, procmid, procupper);
}
