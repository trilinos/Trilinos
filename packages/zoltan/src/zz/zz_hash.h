// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef __ZOLTAN_HASH_H
#define __ZOLTAN_HASH_H

#include <mpi.h>
#include "zoltan_types.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
extern unsigned int Zoltan_Hash(ZOLTAN_ID_PTR, int, unsigned int);
extern unsigned int Zoltan_Recommended_Hash_Size (unsigned int n);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
