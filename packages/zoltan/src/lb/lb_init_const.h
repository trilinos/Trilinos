// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef __LB_INIT_CONST_H
#define __LB_INIT_CONST_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


extern void Zoltan_Migrate_Init(struct Zoltan_Migrate_Struct *);

extern void Zoltan_LB_Init(struct Zoltan_LB_Struct *, int);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
