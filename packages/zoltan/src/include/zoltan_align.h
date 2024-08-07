// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __ZOLTAN_ALIGN_H
#define __ZOLTAN_ALIGN_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* 
 * Plauger alignment algorithm, The Standard C Library.
 * Forces malloc'ed variable size struct alignment.
 * ZOLTAN_ALIGN_VAL is defined in Zoltan/include/zoltan_align.h;
 * values are 0,1,3,7U depending upon machine.
 * (E.g., 7U == 8 byte alignment.)
 */

#ifndef ZOLTAN_ALIGN_VAL
#define ZOLTAN_ALIGN_VAL 7U
#endif

extern int Zoltan_Align(int);
extern size_t Zoltan_Align_size_t(size_t);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
