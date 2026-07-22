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
#include "phg_tree.h"
#include "zoltan_partition_tree.h"

/****************************************************************************/

int
Zoltan_PHG_Partition_Tree_Size(
struct Zoltan_Struct     *zz,
int * tree_size
)
{
/* Return the phg tree size.
   See zoltan_partition_tree.h for a more precise definition of tree size.
 */

static char  *yo = "Zoltan_PHG_Partition_Tree";
int          ierr = ZOLTAN_OK;

  if (zz->LB.Data_Structure == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,
      "No Decomposition Data available. use PHG_KEEP_TREE parameter.");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  if (zz->LB.Method != HYPERGRAPH) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,
      "Function can be used only with LB_METHOD == RCB.");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  Zoltan_PHG_LB_Data* ptr = (Zoltan_PHG_LB_Data*)zz->LB.Data_Structure;
  *tree_size = ptr->tree->size;

End:
  return ierr;
}

int
Zoltan_PHG_Partition_Tree(
struct Zoltan_Struct     *zz,
int    treeNodeIndex,    /* tree node index in zoltan PHG */
int    *lo_index,        /* low index */
int    *hi_index         /* high index */
)
{
/* Return the phg tree node information.
   See zoltan_partition_tree.h for a more precise definition of parameters.
 */

static char              *yo = "Zoltan_PHG_Partition_Tree";
struct Zoltan_PHG_Tree_  *treept; /* tree of RCB cuts */
int                      ierr = ZOLTAN_OK;
int                      loArrayIndex;

  if (zz->LB.Data_Structure == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,
      "No Decomposition Data available. use PHG_KEEP_TREE parameter.");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  if (zz->LB.Method != HYPERGRAPH) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,
      "Function can be used only with LB_METHOD == HYPERGRAPH.");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  Zoltan_PHG_LB_Data* ptr = (Zoltan_PHG_LB_Data*)zz->LB.Data_Structure;
  treept = ptr->tree;

  if (treept->size <= 0) {     /* PHG tree was never created. */
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "No PHG tree saved; "
      " Must set parameter PHG_KEEP_TREE to 1.");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  /* hi will be lo index + 1 */
  /* treept is shifted -2 from true mem ptr */
  loArrayIndex = treeNodeIndex * 2;
  if (treept->size*2 + 2 <= loArrayIndex + 1) {     /* PHG tree not big enough. */
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "PHG tree data requested but"
      " tree was not allocated with this size.");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  *lo_index = -treept->array[loArrayIndex]; /* neg convention */
  *hi_index = treept->array[loArrayIndex+1]; /* hi is low + 1 */

End:
  return ierr;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
