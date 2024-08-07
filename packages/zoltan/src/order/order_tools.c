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
