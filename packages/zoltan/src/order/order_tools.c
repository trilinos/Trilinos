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


#define CHECK_ZOLTAN_FREE(ptr) do { if ((ptr) != NULL) ZOLTAN_FREE(&(ptr)); } while (0)

#include "zz_const.h"
#include "order_const.h"


/****************************************************************************
 *  Function which initializes elimination tree structure.
 ****************************************************************************/

int  Zoltan_Order_Init_Tree (struct Zoltan_Order_Struct *order, int blocknbr, int leavesnbr)
{
  Zoltan_Order_Free_Struct(order);

  order->ancestor = (int *) ZOLTAN_MALLOC(blocknbr*sizeof(int));
  order->start = (int *) ZOLTAN_MALLOC((blocknbr+1)*sizeof(int));
  order->leaves = (int *) ZOLTAN_MALLOC((leavesnbr+1)*sizeof(int));

  if ((order->ancestor == NULL) || (order->start == NULL) || (order->leaves == NULL)) {
    Zoltan_Order_Free_Struct(order);
    return (ZOLTAN_MEMERR);
  }
  order->needfree = 1;
  return (ZOLTAN_OK);
}


/****************************************************************************
 *  Function which frees the memory used by the ordering structure
 ****************************************************************************/
void Zoltan_Order_Free_Struct(struct Zoltan_Order_Struct *order)
{
  if (order->needfree == 0)
    return;

  CHECK_ZOLTAN_FREE(order->start);
  CHECK_ZOLTAN_FREE(order->ancestor);
  CHECK_ZOLTAN_FREE(order->leaves);

  order->needfree = 0;
}

#ifdef __cplusplus
}
#endif
