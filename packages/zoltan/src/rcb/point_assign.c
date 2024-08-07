// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdio.h>
#include "zz_const.h"
#include "zz_util_const.h"
#include "rcb.h"
#include "rib.h"

int Zoltan_RB_Point_Assign(
ZZ       *zz,                   /* The Zoltan structure */
double   *coords,
int      *proc,                 /* processor that point lands in;
                                   if NULL, processor info is not returned. */
int      *part                  /* part that point lands in; 
                                   if NULL, part info is not returned. */
)
{
/* Locate which processor a point is inside within the tree defined
   by the recursive bisection algorithm chosen. */

     char             *yo = "Zoltan_RB_Point_Assign";
     int           partmid = 0; /* 1st partition in upper half */
     RCB_STRUCT        *rcb;    /* Pointer to data structures for RCB.  */
     struct rcb_tree   *treept; /* tree of RCB cuts */
     RIB_STRUCT        *rib;    /* Pointer to data structures for RIB. */
     struct rib_tree   *itree;  /* tree of RIB cuts */
     int ierr = ZOLTAN_OK;
     int num_geom;
     volatile double t;         /* Temporary variable; volatile to get matching
                                   results with and without optimization. */
     double cnew[3];
     double *c = coords;

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

        if (rcb->Tran.Target_Dim > 0){  /* degenerate geometry */
          Zoltan_Transform_Point(c, rcb->Tran.Transformation, 
            rcb->Tran.Permutation, rcb->Num_Dim, rcb->Tran.Target_Dim, cnew);
          c = cnew;
        }

        partmid = treept[0].right_leaf;

        while (partmid > 0)
           if (c[treept[partmid].dim] <= treept[partmid].cut)
              partmid = treept[partmid].left_leaf;
           else
              partmid = treept[partmid].right_leaf;
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

        if (rib->Tran.Target_Dim > 0){ /* degenerate geometry */
          Zoltan_Transform_Point(c, rib->Tran.Transformation, 
            rib->Tran.Permutation, rib->Num_Geom, rib->Tran.Target_Dim, cnew);
          c = cnew;
          num_geom = rib->Tran.Target_Dim;
        }
        else{
          num_geom = rib->Num_Geom;
        }

        switch (num_geom) {
           case 3:
              while (partmid > 0) {
                 t = ((c[0] - itree[partmid].cm[0])*itree[partmid].ev[0]) +
                     ((c[1] - itree[partmid].cm[1])*itree[partmid].ev[1]) +
                     ((c[2] - itree[partmid].cm[2])*itree[partmid].ev[2]);
                 if (t <= itree[partmid].cut)
                    partmid = itree[partmid].left_leaf;
                 else
                    partmid = itree[partmid].right_leaf;
              }
              break;
           case 2:
              while (partmid > 0) {
                 t = ((c[0] - itree[partmid].cm[0])*itree[partmid].ev[0]) +
                     ((c[1] - itree[partmid].cm[1])*itree[partmid].ev[1]);
                 if (t <= itree[partmid].cut)
                    partmid = itree[partmid].left_leaf;
                 else
                    partmid = itree[partmid].right_leaf;
              }
              break;
           case 1:
              while (partmid > 0){
                 if (c[0] <= itree[partmid].cut){
                    partmid = itree[partmid].left_leaf;
                 }
                 else{
                    partmid = itree[partmid].right_leaf;
                 }
              }
              break;
        }
     }

     if (part != NULL) {
        if (zz->LB.Remap)
           *part = zz->LB.Remap[-partmid];
        else
           *part = -partmid;
     }

     if (proc != NULL) {
        if (zz->LB.Remap)
           *proc = Zoltan_LB_Part_To_Proc(zz, zz->LB.Remap[-partmid], NULL);
        else 
           *proc = Zoltan_LB_Part_To_Proc(zz, -partmid, NULL);
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
