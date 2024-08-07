// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __PHG_CONST_H
#define __PHG_CONST_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

int Zoltan_PHG_Set_Param (char*, char*);

#define PHG_ADD_UNIT_WEIGHT 1
#define PHG_ADD_PINS_WEIGHT 2
#define PHG_ADD_NO_WEIGHT 3

#define PHG_MAX_EDGE_WEIGHTS 1
#define PHG_ADD_EDGE_WEIGHTS 2
#define PHG_FLAG_ERROR_EDGE_WEIGHTS 3

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* __PHG_CONST_H */
