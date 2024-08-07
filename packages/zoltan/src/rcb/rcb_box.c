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
#include <math.h>
#include <float.h>
#include "zz_const.h"
#include "rcb.h"

static void compute_box(ZZ *, int, struct rcb_tree *, int, struct rcb_box *);

/****************************************************************************/
int Zoltan_RCB_Box(
ZZ     *zz,             /* The Zoltan structure */
int     part,           /* Partition whose box should be returned */
int    *ndim,           /* Number of dimensions in the geometry partitioned
                           (and, thus, in the RCB box for the part) */
double *xmin,           /* lower x extent of box */
double *ymin,           /* lower y extent of box */
double *zmin,           /* lower z extent of box */
double *xmax,           /* upper x extent of box */
double *ymax,           /* upper y extent of box */
double *zmax            /* upper z extent of box */
)
{
/* Return the bounding box for a processor's subdomain.
 */

static char       *yo = "Zoltan_RCB_Box";
RCB_STRUCT        *rcb;    /* Pointer to data structures for RCB. */
struct rcb_tree   *treept; /* tree of RCB cuts */
struct rcb_box     box;     /* box data structure */
int                i, ierr = ZOLTAN_OK;

  box.lo[0] = -DBL_MAX;
  box.lo[1] = -DBL_MAX;
  box.lo[2] = -DBL_MAX;
  box.hi[0] = DBL_MAX;
  box.hi[1] = DBL_MAX;
  box.hi[2] = DBL_MAX;
  *ndim = -1;

  if (zz->LB.Data_Structure == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
      "No Decomposition Data available; use KEEP_CUTS parameter.");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  if (zz->LB.Method != RCB) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
      "Function can be used only with LB_METHOD == RCB.");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  if (part < 0 || part >= zz->LB.Num_Global_Parts) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
      "Invalid part number.");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  if (zz->LB.Remap) /* Partitions are re-mapped; need to find old part number
                       before going through tree of cuts. */
    for (i = 0; i < zz->LB.Num_Global_Parts; i++)
      if (zz->LB.Remap[i] == part) {
        part = i;
        break;
      }

  rcb = (RCB_STRUCT *) (zz->LB.Data_Structure);
  treept = rcb->Tree_Ptr;
  if (treept[0].dim < 0) {     /* RCB tree was never created. */
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "No RCB tree saved; "
       " Must set parameter KEEP_CUTS to 1.");
     ierr = ZOLTAN_FATAL;
     goto End;
  }

  if (treept[0].right_leaf > 0)
    compute_box(zz, part, treept, treept[0].right_leaf, &box);
  else{
    box.lo[0] = box.lo[1] = box.lo[2] = 0;
    box.hi[0] = box.hi[1] = box.hi[2] = 0;
  }

  *ndim = rcb->Num_Dim;

End:

  *xmin = box.lo[0];
  *ymin = box.lo[1];
  *zmin = box.lo[2];
  *xmax = box.hi[0];
  *ymax = box.hi[1];
  *zmax = box.hi[2];

  return ierr;
}

/****************************************************************************/
static void compute_box(
ZZ              *zz,
int              part,          /* Partition whose box is being computed */
struct rcb_tree *treept,        /* RCB tree */
int              node,          /* node of the RCB tree to be examined */
struct rcb_box  *boxpt          /* extended box */
)
{
/* Routine that traverses RCB tree and assigns lo and hi coordinates of box */
struct rcb_tree *nodept = &(treept[node]);

  if (part >= node) {
    boxpt->lo[nodept->dim] = nodept->cut;
    if (nodept->right_leaf > 0)
      compute_box(zz, part, treept, nodept->right_leaf, boxpt);
  }
  else {
    boxpt->hi[nodept->dim] = nodept->cut;
    if (nodept->left_leaf > 0)
      compute_box(zz, part, treept, nodept->left_leaf, boxpt);
  }
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
