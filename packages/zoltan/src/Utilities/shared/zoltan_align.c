// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stddef.h>
#include "zoltan_align.h"
#include "zoltan_util.h"

/*****************************************************************************/
/*
 *  Routines for properly aligning data.
 */
/*****************************************************************************/

/* 
 * Plauger alignment algorithm, The Standard C Library.
 * Forces malloc'ed variable size struct alignment.
 * ZOLTAN_ALIGN_VAL is defined in Zoltan/include/zoltan_align.h;
 * values are 0,1,3,7U depending upon machine.
 */

int Zoltan_Align(int a)
{
return((ZOLTAN_ALIGN_VAL + a) & ~ZOLTAN_ALIGN_VAL);
}

size_t Zoltan_Align_size_t(size_t a)
{
return((ZOLTAN_ALIGN_VAL + a) & ~ZOLTAN_ALIGN_VAL);
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
