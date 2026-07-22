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

#include "zz_const.h"
#include "rcb.h"
#include "zoltan_partition_tree.h"

/****************************************************************************/

int Zoltan_RCB_Partition_Tree(
struct Zoltan_Struct     *zz,              /* The Zoltan structure */
int    treeNodeIndex,    /* tree node index in zoltan rcb */
int    *parent,          /* parent partition number */
int    *left_leaf,       /* left leaf partition number */
int    *right_leaf       /* right leaf partition number */
)
{
/* Return the rcb tree information.
 */

static char       *yo = "Zoltan_RCB_Partition_Tree";
struct rcb_tree   *treept; /* tree of RCB cuts */
int                ierr = ZOLTAN_OK;

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

  RCB_STRUCT * rcb = (RCB_STRUCT *) (zz->LB.Data_Structure);
  treept = rcb->Tree_Ptr;

  if (treept[0].dim < 0) {     /* RCB tree was never created. */
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "No RCB tree saved; "
      " Must set parameter KEEP_CUTS to 1.");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  *parent = treept[treeNodeIndex].parent;
  *left_leaf = treept[treeNodeIndex].left_leaf;
  *right_leaf = treept[treeNodeIndex].right_leaf;

End:
  return ierr;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
