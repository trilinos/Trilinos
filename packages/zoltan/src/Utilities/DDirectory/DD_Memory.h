// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ZOLTAN_DD_MEMORY_H_
#define ZOLTAN_DD_MEMORY_H_

#include <stdio.h>
#include <stdlib.h>

#include "zoltan_dd_const.h"


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

extern int *DD_Memory_Alloc_Nodelist(Zoltan_DD_Directory *, DD_NodeIdx , float);
extern DD_NodeIdx DD_Memory_Alloc_Node(Zoltan_DD_Directory *);
extern void DD_Memory_Free_Node(Zoltan_DD_Directory*, DD_NodeIdx);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
