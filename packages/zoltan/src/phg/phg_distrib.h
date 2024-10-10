// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef __ZOLTAN_PHG_DISTRIB_H
#define __ZOLTAN_PHG_DISTRIB_H

#include "phg.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

int Zoltan_PHG_Gno_To_Proc_Block(ZOLTAN_GNO_TYPE gno, ZOLTAN_GNO_TYPE  *dist_dim, int nProc_dim);

    
int Zoltan_PHG_Redistribute(
    ZZ *zz,
    PHGPartParams *hgp,     /* Input: parameters; used only for user's
                               request of nProc_x and nProc_y */
    HGraph  *ohg,           /* Input: Local part of distributed hypergraph */
    int     lo, int hi,     /* Input: range of proc ranks (inclusive)
                               to be included in new communicator: ncomm */
    PHGComm *ncomm,         /* Output: Communicators of new distribution */
    HGraph  *nhg,           /* Output: Newly redistributed hypergraph */
    int     **vmap,         /* Output: allocated with the size nhg->nVtx and
                               vertex map from nhg to ohg's local vertex number*/
    int     **vdest         /* Output: allocated with the size nhg->nVtx and
                               stores dest proc in ocomm */    
    );
    

    

    
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* __ZOLTAN_PHG_DISTRIB_H */
