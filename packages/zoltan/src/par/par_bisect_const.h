// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef __PAR_BISECT_CONST_H
#define __PAR_BISECT_CONST_H

#include <mpi.h>
#include "zz_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

extern int Zoltan_RB_find_bisector(ZZ *, int, double *, double *, double,
  int *, int, int, int, double *, MPI_Comm,
  double *, int, int, int, int, double, double, 
  double *, double *, double *, double *, int *, int, int);

/* Note: MAX_BISECT_WGTS should be >= RB_MAX_WEIGHTS in RCB. */
/* EBEB: Should we rather include rcb_const.h and use RB_MAX_WEIGHTS? */
#define MAX_BISECT_WGTS 8

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
