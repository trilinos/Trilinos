/*
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */

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
