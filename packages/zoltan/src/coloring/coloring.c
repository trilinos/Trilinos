/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <ctype.h>
#include "zz_const.h"
#include "zz_util_const.h"
#include "parmetis_jostle.h"
    
/**********************************************************/
/* Interface routine for Graph Coloring                   */
/**********************************************************/

int Zoltan_Color(
  ZZ *zz,               /* Zoltan structure */
  int distance,         /* Input: which coloring to perform;
                           currently only supports D1 and D2 coloring */
  int nvtx,             /* Input: number of vertices on this proc */
  ZOLTAN_ID_PTR *GIDs,  /* Input: global ids of the vertices */
  ZOLTAN_ID_PTR *LIDs,  /* Input: local  ids of the vertices */
  int *color            /* Output: Colors assigned to local vertices */
)
{
  idxtype *vtxdist, *xadj, *adjncy, *vwgt, *adjwgt;
  int *adjproc;

  vtxdist = xadj = adjncy = vwgt = adjwgt = adjproc = NULL;

#if 0
  /* Build ParMetis data structures, or just get vtxdist. */
  ierr = Zoltan_Build_Graph(zz, graph_type, check_graph, num_obj,
         global_ids, local_ids, obj_wgt_dim, edge_wgt_dim,
         &vtxdist, &xadj, &adjncy, &ewgts, &adjproc);

#endif
    return ZOLTAN_OK;
}


#ifdef __cplusplus
}
#endif
