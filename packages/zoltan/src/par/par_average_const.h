// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef __PAR_AVERAGE_CONST_H
#define __PAR_AVERAGE_CONST_H

#include <mpi.h>

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

double Zoltan_RB_Average_Cut(int, double *, int *, int, int, int, int,
  MPI_Comm, double);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
