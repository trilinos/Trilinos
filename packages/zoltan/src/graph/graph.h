// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __GRAPH_H
#define __GRAPH_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "zoltan_matrix.h"

typedef struct ZG_ {
  Zoltan_matrix_2d mtx;
  int         *fixed_vertices; 
  int          bipartite;
  int          fixObj;
} ZG;

int
Zoltan_ZG_Build (ZZ* zz, ZG* graph, int local, int, int, ZOLTAN_ID_PTR,
                 ZOLTAN_GNO_TYPE *);

int
Zoltan_ZG_Export (ZZ* zz, const ZG* const graph, ZOLTAN_GNO_TYPE *gvtx, int *nvtx, int *obj_wgt_dim, int *edge_wgt_dim,
	   ZOLTAN_GNO_TYPE **vtxdist, int **xadj, ZOLTAN_GNO_TYPE **adjncy, int **adjproc,
	   /* float **xwgt, */ float **ewgt, int **partialD2);

int
Zoltan_ZG_Vertex_Info(ZZ* zz, const ZG *const graph,
		      ZOLTAN_ID_PTR *pgid, ZOLTAN_ID_PTR *plid, float **pwwgt, int **pinput_part);

int
Zoltan_ZG_Register(ZZ* zz, ZG* graph, int* properties);

int
Zoltan_ZG_Query (ZZ* zz, const ZG *graph, const ZOLTAN_ID_PTR GID,
	  int GID_length, int* properties);

void
Zoltan_ZG_Free(ZG *m);


#ifdef __cplusplus
}
#endif

#endif
