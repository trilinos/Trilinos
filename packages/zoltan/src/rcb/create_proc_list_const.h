// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef __CREATE_PROC_LIST_CONST_H
#define __CREATE_PROC_LIST_CONST_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


/* function prototype */

extern int Zoltan_RB_Create_Proc_List(ZZ *, int, int, int, int *, MPI_Comm, int, int);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
