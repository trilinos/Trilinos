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
#include "zz_const.h"
#include "rcb.h"
#include "rib.h"


int Zoltan_RB_Point_Assign(
ZZ       *zz,                   /* The Zoltan structure */
double   *coords,               /* vector of point coordinates */
int      *proc,                 /* processor that point lands in;
                                   if NULL, processor info is not returned. */
int      *part                  /* partition that point lands in; 
                                   if NULL, partition info is not returned. */
)
{
/* Locate which processor a point is inside within the tree defined
   by the recursive bisection algorithm chosen. */

     char             *yo = "Zoltan_RB_Point_Assign";
     int               partmid; /* 1st partition in upper half */
     RCB_STRUCT        *rcb;    /* Pointer to data structures for RCB.  */
     struct rcb_tree   *treept; /* tree of RCB cuts */
     RIB_STRUCT        *rib;    /* Pointer to data structures for RIB. */
     struct rib_tree   *itree;  /* tree of RIB cuts */
     int ierr = ZOLTAN_OK;
     volatile double t;         /* Temporary variable; volatile to get matching
                                   results with and without optimization. */

     if (zz->LB.Data_Structure == NULL) {
        ZOLTAN_PRINT_ERROR(-1, yo, 
                   "No Decomposition Data available; use KEEP_CUTS parameter.");
        ierr = ZOLTAN_FATAL;
        goto End;
     }

     if (zz->LB.Method == RCB) {
        rcb = (RCB_STRUCT *) (zz->LB.Data_Structure);
        treept = rcb->Tree_Ptr;
        if (treept[0].dim < 0) { /* RCB tree was never created. */
           ZOLTAN_PRINT_ERROR(zz->Proc, yo, "No RCB tree saved; "
                                        "Must set parameter KEEP_CUTS to 1.");
           ierr = ZOLTAN_FATAL;
           goto End;
        }

        partmid = treept[0].right_leaf;

        while (partmid > 0)
           if (coords[treept[partmid].dim] <= treept[partmid].cut)
              partmid = treept[partmid].left_leaf;
           else
              partmid = treept[partmid].right_leaf;

        if (part != NULL)
           *part = -partmid;
        if (proc != NULL)
           *proc = Zoltan_LB_Part_To_Proc(zz, -partmid, NULL);
     }
     else if (zz->LB.Method == RIB) {
        rib = (RIB_STRUCT *) (zz->LB.Data_Structure);
        itree = rib->Tree_Ptr;
        if ((partmid = itree[0].right_leaf) < 0) { /* RIB tree never created */
           ZOLTAN_PRINT_ERROR(zz->Proc, yo, "No RIB tree saved; "
                                     "Must set parameter KEEP_CUTS to 1.");
           ierr = ZOLTAN_FATAL;
           goto End;
        }

        switch (rib->Num_Geom) {
           case 3:
              while (partmid > 0) {
                 t = ((coords[0] - itree[partmid].cm[0])*itree[partmid].ev[0]) +
                     ((coords[1] - itree[partmid].cm[1])*itree[partmid].ev[1]) +
                     ((coords[2] - itree[partmid].cm[2])*itree[partmid].ev[2]);
                 if (t <= itree[partmid].cut)
                    partmid = itree[partmid].left_leaf;
                 else
                    partmid = itree[partmid].right_leaf;
              }
              break;
           case 2:
              while (partmid > 0) {
                 t = ((coords[0] - itree[partmid].cm[0])*itree[partmid].ev[0]) +
                     ((coords[1] - itree[partmid].cm[1])*itree[partmid].ev[1]);
                 if (t <= itree[partmid].cut)
                    partmid = itree[partmid].left_leaf;
                 else
                    partmid = itree[partmid].right_leaf;
              }
              break;
           case 1:
              while (partmid > 0)
                 if (coords[0] <= itree[partmid].cut)
                    partmid = itree[partmid].left_leaf;
                 else
                    partmid = itree[partmid].right_leaf;
              break;
        }

        if (part != NULL)
           *part = -partmid;
        if (proc != NULL)
           *proc = Zoltan_LB_Part_To_Proc(zz, -partmid, NULL);
     }
     else {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
          "Valid only when load-balancing method is RCB or RIB.");
        ierr = ZOLTAN_FATAL;
        goto End;
     }
End:
     if (ierr == ZOLTAN_FATAL) {
        if (part != NULL)
           *part = -1;
        if (proc != NULL)
           *proc = -1;
     }
     return ierr;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
