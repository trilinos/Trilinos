// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef __COMMON_H
#define __COMMON_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "zoltan_util.h"
#include "third_library_const.h"

extern int Zoltan_Verify_Graph(MPI_Comm, indextype *, indextype *,
       indextype *, weighttype *, weighttype *,
       int, int, int, int, int);

extern int Zoltan_Scatter_Graph(indextype **, indextype **,
       indextype **, weighttype **, indextype **, weighttype **,
       realtype **, int, int, ZZ *, ZOLTAN_COMM_OBJ **);


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
