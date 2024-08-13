// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __COLORING_H
#define __COLORING_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include <ctype.h>
#include "zz_const.h"
#include "zz_util_const.h"
#include "coloring_const.h"

/* Metis also has a swap */
#ifdef SWAP
#undef SWAP
#endif

#define SWAP(a,b) tmp=(a);(a)=(b);(b)=tmp;

/* Macros for error handling */
#define ZOLTAN_COLOR_ERROR(error,str) {ierr = error ; \
 ZOLTAN_PRINT_ERROR(zz->Proc, yo, str) ; goto End ;}

#ifdef __cplusplus
}
#endif
#endif
