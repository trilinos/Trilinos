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
#include <math.h>
#include "zz_const.h"
#include "rcb.h"
#include "rib.h"

static void Box_Assign(struct rcb_tree *, struct rcb_box *, int *, int *, int);
static void Box_Assign3(struct rib_tree *,struct rcb_box *, int *, int *, int);
static void Box_Assign2(struct rib_tree *,struct rcb_box *, int *, int *, int);
static void Box_Assign1(struct rib_tree *,struct rcb_box *, int *, int *, int);

int Zoltan_LB_Box_Assign(
ZZ             *zz,             /* The Zoltan structure */
double          xmin,           /* lower x extent of box */
double          ymin,           /* lower y extent of box */
double          zmin,           /* lower z extent of box */
double          xmax,           /* upper x extent of box */
double          ymax,           /* upper y extent of box */
double          zmax,           /* upper z extent of box */
int            *procs,          /* list of procs that box intersects */
int            *numprocs)       /* number of processors in proc list */
{
/* Determine which processors a box intersects.
   Currently assumes that partitioning has used RCB or RIB, but should be
   modified to return an error message if other method was used */

     static char       *yo = "Zoltan_LB_Box_Assign";
     RCB_STRUCT        *rcb;    /* Pointer to data structures for RCB. */
     struct rcb_tree   *treept; /* tree of RCB cuts */
     RIB_STRUCT        *rib;    /* Pointer to data structures for RIB. */
     struct rib_tree   *itree;  /* tree of RIB cuts */
     struct rcb_box    box;     /* box data structure */

     if (zz->LB.Data_Structure == NULL) {
        ZOLTAN_PRINT_ERROR(-1, yo, 
          "No Decomposition Data available; use KEEP_CUTS parameter.");
        *procs = -1;
        *numprocs = 0;
        return(ZOLTAN_FATAL);
     }

     if (zz->LB.Method == RCB) {
        rcb = (RCB_STRUCT *) (zz->LB.Data_Structure);
        treept = rcb->Tree_Ptr;
        if (treept[0].dim < 0) {     /* RCB tree was never created. */
           ZOLTAN_PRINT_ERROR(zz->Proc, yo, "No RCB tree saved; "
             " Must set parameter KEEP_CUTS to 1.");
           *procs = -1;
           *numprocs = 0;
           return(ZOLTAN_FATAL);
        }

        *numprocs = 0;
        box.lo[0] = xmin;
        box.lo[1] = ymin;
        box.lo[2] = zmin;
        box.hi[0] = xmax;
        box.hi[1] = ymax;
        box.hi[2] = zmax;

        Box_Assign(treept, &box, procs, numprocs, treept[0].right_leaf);
     }
     else if (zz->LB.Method == RIB) {
        rib = (RIB_STRUCT *) (zz->LB.Data_Structure);
        itree = rib->Tree_Ptr;
        if (itree[0].right_leaf < 0) { /* RIB tree was never created. */
           ZOLTAN_PRINT_ERROR(zz->Proc, yo, "No RIB tree saved;"
             " Must set parameter KEEP_CUTS to 1.");
           *procs = -1;
           return(ZOLTAN_FATAL);
        }

        *numprocs = 0;

        switch (rib->Num_Geom) {
           case 3:
              box.lo[0] = xmin;
              box.lo[1] = ymin;
              box.lo[2] = zmin;
              box.hi[0] = xmax;
              box.hi[1] = ymax;
              box.hi[2] = zmax;

              Box_Assign3(itree, &box, procs, numprocs, itree[0].right_leaf);

              break;
           case 2:
              box.lo[0] = xmin;
              box.lo[1] = ymin;
              box.hi[0] = xmax;
              box.hi[1] = ymax;

              Box_Assign2(itree, &box, procs, numprocs, itree[0].right_leaf);

              break;
           case 1:
              box.lo[0] = xmin;
              box.hi[0] = xmax;

              Box_Assign1(itree, &box, procs, numprocs, itree[0].right_leaf);

              break;
        }
     }
     else {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                       "Box_Assign valid only when method is RCB or RIB.");
        *procs = -1;
        *numprocs = 0;
        return(ZOLTAN_FATAL);
     }

   return(ZOLTAN_OK);
}

static void Box_Assign(
struct rcb_tree *treept,        /* RCB tree */
struct rcb_box  *boxpt,         /* extended box */
int             *procs,         /* processors that box is in */
int             *numprocs,      /* current number of processors on list */
int              procmid)       /* 1st processor in upper half */
{
     int       dim;
     double    cut;

     /* end recursion when partition size is a single proc */
     /* add processor to list of processors */

     if (procmid <= 0) {
        procs[*numprocs] = -procmid;
        (*numprocs)++;
        return;
     }

     /* drop box on each side of cut if it extends beyond it */
     /* important to use >= and <= criteria since either proc may own
        boundary */
     /* procmid = 1st processor in upper half of partition, loc that stores cut
                  for this partition in treept */
     /* dim = dimension 0,1,2 of cut for this partition */
     /* cut = position of cut for this partition */

     dim = treept[procmid].dim;
     cut = treept[procmid].cut;

     if (boxpt->lo[dim] <= cut)
        Box_Assign(treept, boxpt, procs, numprocs, treept[procmid].left_leaf);
     if (boxpt->hi[dim] >= cut)
        Box_Assign(treept, boxpt, procs, numprocs, treept[procmid].right_leaf);
}

static void Box_Assign3(
struct rib_tree   *itree,       /* RIB tree */
struct rcb_box    *box,         /* extended box */
int               *procs,       /* processors that box is in */
int               *numprocs,    /* current number of processors on list */
int               procmid)      /* 1st processor in upper half */
{
     double p1[3], p2[3];       /* two points of the box used to test */
     double min, max;           /* values for two points */
     double cut;                /* current cut */

     /* end recursion when partition size is a single proc */
     /* add processor to list of processors */
     if (procmid <= 0) {
        procs[*numprocs] = -procmid;
        (*numprocs)++;
        return;
     }

     /* determine which two points to check based on the direction of
        the eigenvector */
     if (itree[procmid].ev[0] >= 0.0) {
        p1[0] = box->lo[0];
        p2[0] = box->hi[0];
     }
     else {
        p1[0] = box->hi[0];
        p2[0] = box->lo[0];
     }
     if (itree[procmid].ev[1] >= 0.0) {
        p1[1] = box->lo[1];
        p2[1] = box->hi[1];
     }
     else {
        p2[1] = box->lo[1];
        p1[1] = box->hi[1];
     }
     if (itree[procmid].ev[2] >= 0.0) {
        p1[2] = box->lo[2];
        p2[2] = box->hi[2];
     }
     else {
        p2[2] = box->lo[2];
        p1[2] = box->hi[2];
     }

     /* determine the distance from the center of mass point for each point */
     /* this distance is compared to cut which is the distance of the plane */
     min = (p1[0] - itree[procmid].cm[0])*itree[procmid].ev[0] +
           (p1[1] - itree[procmid].cm[1])*itree[procmid].ev[1] +
           (p1[2] - itree[procmid].cm[2])*itree[procmid].ev[2];
     max = (p2[0] - itree[procmid].cm[0])*itree[procmid].ev[0] +
           (p2[1] - itree[procmid].cm[1])*itree[procmid].ev[1] +
           (p2[2] - itree[procmid].cm[2])*itree[procmid].ev[2];

     /* order the points */
     if (min > max) {
        cut = min;
        min = max;
        max = cut;
     }

     cut = itree[procmid].cut;

     if (min <= cut)
        Box_Assign3(itree, box, procs, numprocs, itree[procmid].left_leaf);
     if (max >= cut)
        Box_Assign3(itree, box, procs, numprocs, itree[procmid].right_leaf);
}

static void Box_Assign2(
struct rib_tree   *itree,       /* RIB tree */
struct rcb_box    *box,         /* extended box */
int               *procs,       /* processors that box is in */
int               *numprocs,    /* current number of processors on list */
int               procmid)      /* 1st processor in upper half */
{
     double p1[2], p2[2];       /* two points of the box used to test */
     double min, max;           /* values for two points */
     double cut;                /* current cut */

     /* end recursion when partition size is a single proc */
     /* add processor to list of processors */
     if (procmid <= 0) {
        procs[*numprocs] = -procmid;
        (*numprocs)++;
        return;
     }

     /* determine which two points to check based on the direction of
        the eigenvector */
     if (itree[procmid].ev[0] >= 0.0) {
        p1[0] = box->lo[0];
        p2[0] = box->hi[0];
     }
     else {
        p1[0] = box->hi[0];
        p2[0] = box->lo[0];
     }
     if (itree[procmid].ev[1] >= 0.0) {
        p1[1] = box->lo[1];
        p2[1] = box->hi[1];
     }
     else {
        p2[1] = box->lo[1];
        p1[1] = box->hi[1];
     }

     /* determine the distance from the center of mass point for each point */
     /* this distance is compared to cut which is the distance of the plane */
     min = (p1[0] - itree[procmid].cm[0])*itree[procmid].ev[0] +
           (p1[1] - itree[procmid].cm[1])*itree[procmid].ev[1];
     max = (p2[0] - itree[procmid].cm[0])*itree[procmid].ev[0] +
           (p2[1] - itree[procmid].cm[1])*itree[procmid].ev[1];

     /* order the points using cut as a temporary */
     if (min > max) {
        cut = min;
        min = max;
        max = cut;
     }

     cut = itree[procmid].cut;

     if (min <= cut)
        Box_Assign2(itree, box, procs, numprocs, itree[procmid].left_leaf);
     if (max >= cut)
        Box_Assign2(itree, box, procs, numprocs, itree[procmid].right_leaf);
}

static void Box_Assign1(
struct rib_tree   *itree,       /* RIB tree */
struct rcb_box    *box,         /* extended box */
int               *procs,       /* processors that box is in */
int               *numprocs,    /* current number of processors on list */
int               procmid)      /* 1st processor in upper half */
{
     double cut;                /* current cut */

     /* end recursion when partition size is a single proc */
     /* add processor to list of processors */
     if (procmid <= 0) {
        procs[*numprocs] = -procmid;
        (*numprocs)++;
        return;
     }

     cut = itree[procmid].cut;

     if (box->lo[0] <= cut)
        Box_Assign1(itree, box, procs, numprocs, itree[procmid].left_leaf);
     if (box->hi[0] >= cut)
        Box_Assign1(itree, box, procs, numprocs, itree[procmid].right_leaf);
}
