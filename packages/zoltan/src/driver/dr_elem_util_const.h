// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef _DR_ELEM_UTIL_CONST_H_
#define _DR_ELEM_UTIL_CONST_H_


#include "dr_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* Function prototypes */
extern void initialize_element(ELEM_INFO *elem);
extern void free_mesh_arrays(MESH_INFO_PTR mesh);
extern void free_element_arrays(ELEM_INFO *elem, MESH_INFO_PTR mesh);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
