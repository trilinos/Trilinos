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


/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */
/* It was modified by Courtenay T. Vaughan for use in Zoltan */


#include <stdio.h>
#include "rib.h"
#include "zz_const.h"

int Zoltan_RIB_inertial1d(
     struct Dot_Struct *dotpt,  /* graph data structure for weights */
     int             *dindx,    /* index array into dotpt; if NULL, access dotpt
                                   directly */
     int              dotnum,   /* number of vtxs in graph */
     int              wgtflag,  /* are vertex weights being used? */
     double           cm[3],    /* center of mass in each direction */
     double           evec[3],  /* eigenvector */
     double           *value    /* array for value to sort on */
)
{
     int i, j;               /* loop counter */

     /* Copy values into double precision array. */
     for (j = 0; j < dotnum; j++) {
        i = (dindx ? dindx[j] : j);
        value[j] = dotpt->X[i];
     }

     /* zero unused center of mass and eigenvector */
     cm[0] = cm[1] = cm[2] = 0.0;
     evec[0] = evec[1] = evec[2] = 0.0;

     return(ZOLTAN_OK);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
