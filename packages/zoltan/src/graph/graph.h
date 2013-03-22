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
