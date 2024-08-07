// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __HA_CONST_H
#define __HA_CONST_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


/* Function prototypes */

extern int Zoltan_Get_Processor_Name(ZZ *, char *);

extern int Zoltan_Divide_Machine(ZZ *, int, float *, int, MPI_Comm, int *, 
                                 int *, int *, int *, int *, int *, int *,
                                 double *);
extern int Zoltan_Divide_Parts(ZZ *, int, float *, int, int *, int *, double *);


/* Misc. constants */
#define MAX_PROC_NAME_LEN (MPI_MAX_PROCESSOR_NAME+1)

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
