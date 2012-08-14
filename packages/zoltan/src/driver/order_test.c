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

#include <stdio.h>
#include "zoltan.h"

#undef DEBUG_PRINT

int Zoltan_Order_Test(
  struct Zoltan_Struct *zz,
  int *num_gid_entries,
  int *num_lid_entries,
  int num_obj,
  ZOLTAN_ID_PTR global_ids,
  ZOLTAN_ID_PTR local_ids,
  int *rank,
  int *iperm
)
{
  int i;
  int nbr;
  int first, last;
  int tmp;
#ifdef DEBUG_PRINT
  int ancestor;
#endif /* DEBUG_PRINT */

  nbr = Zoltan_Order_Get_Num_Blocks(zz);
#ifdef DEBUG_PRINT
  fprintf(stderr, "Nbr block : %d\n", nbr);
  fprintf(stderr,"details : \n");
#endif /* DEBUG_PRINT */

  tmp = 0;
  for (i = 0 ; i<nbr ; ++i) {
    int size;

    Zoltan_Order_Get_Block_Bounds(zz, i, &first, &last);
#ifdef DEBUG_PRINT
    fprintf(stderr, "block %i : %d to %d\n", i, first, last);
#endif /* DEBUG_PRINT */
    if (first != tmp) {
      fprintf(stderr, "Error, non consecutive numbering\n");
      return (ZOLTAN_FATAL);
    }
    tmp = last;
    size = Zoltan_Order_Get_Block_Size(zz, i);
    if (size != last - first) {
      fprintf(stderr, "Error, size doesn't match\n");
      return (ZOLTAN_FATAL);
    }

#ifdef DEBUG_PRINT
    ancestor = Zoltan_Order_Get_Block_Parent(zz, i);
    fprintf(stderr, "Father of %d : %d\n", i, ancestor);
#endif /* DEBUG_PRINT */
  }

  {
    int nbrleaves;
    int *leaves;
    int *blocks;

    nbrleaves = Zoltan_Order_Get_Num_Leaves(zz);
#ifdef DEBUG_PRINT
    fprintf(stderr, "Number of leaves : %d\n", nbrleaves);
#endif /* DEBUG_PRINT */

    blocks = (int*)ZOLTAN_MALLOC(nbr*sizeof(int));
    for (i = 0 ; i < nbr ; ++i)
      blocks[i] = 0;

    leaves = (int*)ZOLTAN_MALLOC((nbrleaves+1)*sizeof(int));
    Zoltan_Order_Get_Block_Leaves(zz, leaves);
    for (i = 0 ; i<nbrleaves ; ++i) {
      int parent;
#ifdef DEBUG_PRINT
      fprintf(stderr, "Leaf %d : %d\n", i, leaves[i]);
#endif /* DEBUG_PRINT */
      if (leaves[i] == -1) {
	fprintf(stderr, "Error, leavessize doesn't match\n");
	return (ZOLTAN_FATAL);
      }

      parent = leaves[i];
      do {
	blocks[parent] = 1;
	parent = Zoltan_Order_Get_Block_Parent(zz, parent);
      } while (parent != -1);
    }
    if (leaves[nbrleaves] != -1) {
      fprintf(stderr, "Error, leaves array is not valid\n");
      return (ZOLTAN_FATAL);
    }

    for (i = 0 ; i <nbr ; ++i) {
      if (blocks[i] != 1) {
	fprintf(stderr, "Error, not all nodes in array are in the tree\n");
	return (ZOLTAN_FATAL);
      }
    }
    ZOLTAN_FREE(&leaves);
    ZOLTAN_FREE(&blocks);
  }

 return (ZOLTAN_OK);
}


#ifdef __cplusplus
}
#endif
