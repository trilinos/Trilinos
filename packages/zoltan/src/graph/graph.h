/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
=======
 *    $RCSfile$
 *    $Author$
 *    $Date$
>>>>>>> d38681a1388a71b0ee4b652621f3b2f81a79e97f:packages/zoltan/src/graph/graph.h
 *    $Revision$
 ****************************************************************************/

#ifndef __GRAPH_H
#define __GRAPH_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "matrix.h"

typedef struct ZG_ {
  Zoltan_matrix_2d mtx;
  int         *fixed_vertices;
  int          bipartite;
  int          fixObj;
} ZG;

int
ZG_Build (ZZ* zz, ZG* graph, int bipartite, int fixObj);

int
ZG_Register(ZZ* zz, ZG* graph, int* properties);

int
ZG_Query (ZZ* zz, ZG* graph, ZOLTAN_ID_PTR GID, int GID_length, int* properties);

void
ZG_Free(ZZ *zz, ZG *m);


#ifdef __cplusplus
}
#endif

#endif
