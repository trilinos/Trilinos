/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdio.h>
#include <math.h>
#include "zz_const.h"
#include "rcb.h"
#include "rib.h"

static void Box_Assign(ZZ *, struct rcb_tree *, struct rcb_box *,
  int, int, int *, int *, int *, int);
static void Box_Assign3(ZZ *, struct rib_tree *,struct rcb_box *,
  int, int, int *, int *, int *, int);
static void Box_Assign2(ZZ *, struct rib_tree *,struct rcb_box *,
  int, int, int *, int *, int *, int);
static void Box_Assign1(ZZ *, struct rib_tree *,struct rcb_box *,
  int, int, int *, int *, int *, int);
static void add_to_list(ZZ *, int, int, int *, int *, int *, int);

/****************************************************************************/
int Zoltan_RB_Box_Assign(
ZZ             *zz,             /* The Zoltan structure */
double          xmin,           /* lower x extent of box */
double          ymin,           /* lower y extent of box */
double          zmin,           /* lower z extent of box */
double          xmax,           /* upper x extent of box */
double          ymax,           /* upper y extent of box */
double          zmax,           /* upper z extent of box */
int            *procs,          /* list of procs that box intersects */
int            *numprocs,       /* number of processors in proc list */
int            *parts,          /* list of parts that box intersects */
int            *numparts)       /* number of partitions in part list */
{
/* Determine which partitions and processors a box intersects.
   Currently assumes that partitioning has used RCB or RIB, but should be
   modified to return an error message if other method was used */

     static char       *yo = "Zoltan_RB_Box_Assign";
     RCB_STRUCT        *rcb;    /* Pointer to data structures for RCB. */
     struct rcb_tree   *treept; /* tree of RCB cuts */
     RIB_STRUCT        *rib;    /* Pointer to data structures for RIB. */
     struct rib_tree   *itree;  /* tree of RIB cuts */
     struct rcb_box    box;     /* box data structure */
     int               *proc_array = NULL;  
                                /* Array of size zz->Num_Proc; initialized
                                   to 0; entry i incremented each time 
                                   a found partition is on processor i. 
                                   Added to support 
                                   !zz->LB.Single_Proc_Per_Part. */
     int               include_procs = (procs != NULL);
     int               include_parts = (parts != NULL);
     int               ierr = ZOLTAN_OK;
     int               i;

     if (zz->LB.Data_Structure == NULL) {
        ZOLTAN_PRINT_ERROR(-1, yo, 
          "No Decomposition Data available; use KEEP_CUTS parameter.");
        ierr = ZOLTAN_FATAL;
        goto End;
     }

     if (include_procs) {
        proc_array = (int *) ZOLTAN_CALLOC(zz->Num_Proc, sizeof(int));
        if (!proc_array) {
           ierr = ZOLTAN_MEMERR;
           goto End;
        }
     }

     *numprocs = *numparts = 0;

     if (zz->LB.Method == RCB) {
        rcb = (RCB_STRUCT *) (zz->LB.Data_Structure);
        treept = rcb->Tree_Ptr;
        if (treept[0].dim < 0) {     /* RCB tree was never created. */
           ZOLTAN_PRINT_ERROR(zz->Proc, yo, "No RCB tree saved; "
             " Must set parameter KEEP_CUTS to 1.");
           ierr = ZOLTAN_FATAL;
           goto End;
        }

        box.lo[0] = xmin;
        box.lo[1] = ymin;
        box.lo[2] = zmin;
        box.hi[0] = xmax;
        box.hi[1] = ymax;
        box.hi[2] = zmax;

        Box_Assign(zz, treept, &box, include_procs, include_parts, 
                   proc_array, parts, numparts, treept[0].right_leaf);
     }
     else if (zz->LB.Method == RIB) {
        rib = (RIB_STRUCT *) (zz->LB.Data_Structure);
        itree = rib->Tree_Ptr;
        if (itree[0].right_leaf < 0) { /* RIB tree was never created. */
           ZOLTAN_PRINT_ERROR(zz->Proc, yo, "No RIB tree saved;"
             " Must set parameter KEEP_CUTS to 1.");
           ierr = ZOLTAN_FATAL;
           goto End;
        }


        switch (rib->Num_Geom) {
           case 3:
              box.lo[0] = xmin;
              box.lo[1] = ymin;
              box.lo[2] = zmin;
              box.hi[0] = xmax;
              box.hi[1] = ymax;
              box.hi[2] = zmax;

              Box_Assign3(zz, itree, &box, include_procs, include_parts, 
                          proc_array, parts, numparts, itree[0].right_leaf);

              break;
           case 2:
              box.lo[0] = xmin;
              box.lo[1] = ymin;
              box.hi[0] = xmax;
              box.hi[1] = ymax;

              Box_Assign2(zz, itree, &box, include_procs, include_parts,
                          proc_array, parts, numparts, itree[0].right_leaf);

              break;
           case 1:
              box.lo[0] = xmin;
              box.hi[0] = xmax;

              Box_Assign1(zz, itree, &box, include_procs, include_parts,
                          proc_array, parts, numparts, itree[0].right_leaf);

              break;
        }
     }
     else {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
              "Valid only when load-balancing method is RCB and RIB.");
        ierr = ZOLTAN_FATAL;
        goto End;
     }

     if (include_procs) {
        for (i = 0; i < zz->Num_Proc; i++)
           if (proc_array[i] > 0)
              procs[(*numprocs)++] = i;
     }

End:
     ZOLTAN_FREE(&proc_array);

     if (ierr != ZOLTAN_OK) {
        if (include_procs) *procs = -1;
        *numprocs = 0;
        if (include_parts) *parts = -1;
        *numparts = 0;
     }
     
     return ierr;
}

/****************************************************************************/
static void Box_Assign(
ZZ              *zz,
struct rcb_tree *treept,        /* RCB tree */
struct rcb_box  *boxpt,         /* extended box */
int              include_procs, /* Flag:  Compute proc lists. */
int              include_parts, /* Flag:  Compute part lists. */
int             *proc_array,    /* Array of size Num_Proc; entry i is
                                   incremented if a found partition is on 
                                   proc i. */
int             *parts,         /* partitions that box is in */
int             *numparts,      /* current number of partitions on list */
int              partmid)       /* 1st partition in upper half */
{
     int       dim;
     double    cut;

     /* end recursion when partition size is a single partition */
     /* add partition to list of partitions */

     if (partmid <= 0) {
        add_to_list(zz, include_procs, include_parts, proc_array, 
                    parts, numparts, -partmid);
        return;
     }

     /* drop box on each side of cut if it extends beyond it */
     /* important to use >= and <= criteria since either partition may own
        boundary */
     /* partmid = 1st partition in upper half, loc that stores cut
                  for this partition in treept */
     /* dim = dimension 0,1,2 of cut for this partition */
     /* cut = position of cut for this partition */

     dim = treept[partmid].dim;
     cut = treept[partmid].cut;

     if (boxpt->lo[dim] <= cut)
        Box_Assign(zz, treept, boxpt, include_procs, include_parts, proc_array,
                    parts, numparts, treept[partmid].left_leaf);
     if (boxpt->hi[dim] >= cut)
        Box_Assign(zz, treept, boxpt, include_procs, include_parts, proc_array,
                    parts, numparts, treept[partmid].right_leaf);
}

/****************************************************************************/
static void Box_Assign3(
ZZ              *zz,
struct rib_tree *itree,         /* RIB tree */
struct rcb_box  *box,           /* extended box */
int              include_procs, /* Flag:  Compute proc lists. */
int              include_parts, /* Flag:  Compute part lists. */
int             *proc_array,    /* Array of size Num_Proc; entry i is
                                   incremented if a found partition is on 
                                   proc i. */
int             *parts,         /* partitions that box is in */
int             *numparts,      /* current number of partitions on list */
int              partmid)       /* 1st partition in upper half */
{
     double p1[3], p2[3];       /* two points of the box used to test */
     volatile double min, max;  /* values for two points */
     double cut;                /* current cut */

     /* end recursion when partition size is a single partition */
     /* add partition to list of partitions */
     if (partmid <= 0) {
        add_to_list(zz, include_procs, include_parts, proc_array, 
                    parts, numparts, -partmid);
        return;
     }

     /* determine which two points to check based on the direction of
        the eigenvector */
     if (itree[partmid].ev[0] >= 0.0) {
        p1[0] = box->lo[0];
        p2[0] = box->hi[0];
     }
     else {
        p1[0] = box->hi[0];
        p2[0] = box->lo[0];
     }
     if (itree[partmid].ev[1] >= 0.0) {
        p1[1] = box->lo[1];
        p2[1] = box->hi[1];
     }
     else {
        p2[1] = box->lo[1];
        p1[1] = box->hi[1];
     }
     if (itree[partmid].ev[2] >= 0.0) {
        p1[2] = box->lo[2];
        p2[2] = box->hi[2];
     }
     else {
        p2[2] = box->lo[2];
        p1[2] = box->hi[2];
     }

     /* determine the distance from the center of mass point for each point */
     /* this distance is compared to cut which is the distance of the plane */
     min = (p1[0] - itree[partmid].cm[0])*itree[partmid].ev[0] +
           (p1[1] - itree[partmid].cm[1])*itree[partmid].ev[1] +
           (p1[2] - itree[partmid].cm[2])*itree[partmid].ev[2];
     max = (p2[0] - itree[partmid].cm[0])*itree[partmid].ev[0] +
           (p2[1] - itree[partmid].cm[1])*itree[partmid].ev[1] +
           (p2[2] - itree[partmid].cm[2])*itree[partmid].ev[2];

     /* order the points */
     if (min > max) {
        cut = min;
        min = max;
        max = cut;
     }

     cut = itree[partmid].cut;

     if (min <= cut)
        Box_Assign3(zz, itree, box, include_procs, include_parts, proc_array,
                    parts, numparts, itree[partmid].left_leaf);
     if (max >= cut)
        Box_Assign3(zz, itree, box, include_procs, include_parts, proc_array,
                    parts, numparts, itree[partmid].right_leaf);
}

/****************************************************************************/
static void Box_Assign2(
ZZ              *zz,
struct rib_tree *itree,         /* RIB tree */
struct rcb_box  *box,           /* extended box */
int              include_procs, /* Flag:  Compute proc lists. */
int              include_parts, /* Flag:  Compute part lists. */
int             *proc_array,    /* Array of size Num_Proc; entry i is
                                   incremented if a found partition is on 
                                   proc i. */
int             *parts,         /* partitions that box is in */
int             *numparts,      /* current number of partitions on list */
int              partmid)       /* 1st partition in upper half */
{
     double p1[2], p2[2];       /* two points of the box used to test */
     volatile double min, max;  /* values for two points */
     double cut;                /* current cut */

     /* end recursion when partition size is a single partition */
     /* add partition to list of partitions */
     if (partmid <= 0) {
        add_to_list(zz, include_procs, include_parts, proc_array, 
                    parts, numparts, -partmid);
        return;
     }

     /* determine which two points to check based on the direction of
        the eigenvector */
     if (itree[partmid].ev[0] >= 0.0) {
        p1[0] = box->lo[0];
        p2[0] = box->hi[0];
     }
     else {
        p1[0] = box->hi[0];
        p2[0] = box->lo[0];
     }
     if (itree[partmid].ev[1] >= 0.0) {
        p1[1] = box->lo[1];
        p2[1] = box->hi[1];
     }
     else {
        p2[1] = box->lo[1];
        p1[1] = box->hi[1];
     }

     /* determine the distance from the center of mass point for each point */
     /* this distance is compared to cut which is the distance of the plane */
     min = (p1[0] - itree[partmid].cm[0])*itree[partmid].ev[0] +
           (p1[1] - itree[partmid].cm[1])*itree[partmid].ev[1];
     max = (p2[0] - itree[partmid].cm[0])*itree[partmid].ev[0] +
           (p2[1] - itree[partmid].cm[1])*itree[partmid].ev[1];

     /* order the points using cut as a temporary */
     if (min > max) {
        cut = min;
        min = max;
        max = cut;
     }

     cut = itree[partmid].cut;

     if (min <= cut)
        Box_Assign2(zz, itree, box, include_procs, include_parts, proc_array,
                    parts, numparts, itree[partmid].left_leaf);
     if (max >= cut)
        Box_Assign2(zz, itree, box, include_procs, include_parts, proc_array,
                    parts, numparts, itree[partmid].right_leaf);
}

/****************************************************************************/
static void Box_Assign1(
ZZ              *zz,
struct rib_tree *itree,         /* RIB tree */
struct rcb_box  *box,           /* extended box */
int              include_procs, /* Flag:  Compute proc lists. */
int              include_parts, /* Flag:  Compute part lists. */
int             *proc_array,    /* Array of size Num_Proc; entry i is
                                   incremented if a found partition is on 
                                   proc i. */
int             *parts,         /* partitions that box is in */
int             *numparts,      /* current number of partitions on list */
int              partmid)       /* 1st partition in upper half */
{
     double cut;                /* current cut */

     /* end recursion when partition size is a single partition */
     /* add partition to list of partitions */
     if (partmid <= 0) {
        add_to_list(zz, include_procs, include_parts, proc_array, 
                    parts, numparts, -partmid);
        return;
     }

     cut = itree[partmid].cut;

     if (box->lo[0] <= cut)
        Box_Assign1(zz, itree, box, include_procs, include_parts, proc_array,
                    parts, numparts, itree[partmid].left_leaf);
     if (box->hi[0] >= cut)
        Box_Assign1(zz, itree, box, include_procs, include_parts, proc_array,
                    parts, numparts, itree[partmid].right_leaf);
}

/****************************************************************************/
void add_to_list(
  ZZ *zz,
  int include_procs,  /* Flag:  Compute proc lists. */
  int include_parts,  /* Flag:  Compute part lists. */
  int *proc_array,    /* Array of size Num_Proc; entry i is
                         incremented if a found partition is on proc i. */
  int *parts,         /* partitions that box is in */
  int *numparts,      /* current number of partitions on list */
  int  add_part       /* partition to be added to list */
)
{
/* Adds partitions and processors to arrays that for the box */
int last_proc;        /* First processor for partition add_part+1. */
int add_proc;         /* First processor for partition add_part.   */
int i;

     if (include_parts) {
        /* Add partition to partition list */
        parts[*numparts] = add_part;
        (*numparts)++;
     }

     if (include_procs) {
        /* Increment appropriate entry of proc_array for partition add_part. */
        add_proc = Zoltan_LB_Part_To_Proc(zz, add_part, NULL);
        proc_array[add_proc]++;
        if (!zz->LB.Single_Proc_Per_Part) {
           /* Partition may be spread across multiple procs.
              Include them all. */
           if (add_part < zz->LB.Num_Global_Parts - 1)
              last_proc = Zoltan_LB_Part_To_Proc(zz, add_part+1, NULL);
           else
              last_proc = zz->Num_Proc;
  
           for (i = add_proc+1; i < last_proc; i++)
              proc_array[i]++;
        }
     }
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
