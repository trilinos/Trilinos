// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __ZOLTAN_ID_H
#define __ZOLTAN_ID_H

#include "zoltan_types.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/*
 *  This file contains the data types and comparison functions
 *  for IDs used by Zoltan and its Utilities.  The basic data type
 *  is ZOLTAN_ID.
 */

/****************************************************************************/
/*
 * Default value of ZOLTAN_ID_TYPE; 
 * IDs allocated with ZOLTAN_MALLOC_ID are initialized with this value.
 */
#define ZOLTAN_ID_DEFAULT 0

/*
 * Macros for initializing single IDs.
 */

#define ZOLTAN_INIT_ID(n,id) \
  {int ZOLTAN_ID_LOOP;       \
  for (ZOLTAN_ID_LOOP = 0; ZOLTAN_ID_LOOP < (n); ZOLTAN_ID_LOOP++)  \
    (id)[ZOLTAN_ID_LOOP] = ZOLTAN_ID_DEFAULT;                   \
  }

/****************************************************************************/
/*
 *  Macros to copy IDs.
 */
#define ZOLTAN_SET_ID(n,a,b)                                            \
   {int ZOLTAN_ID_LOOP;                                                 \
    for (ZOLTAN_ID_LOOP = 0; ZOLTAN_ID_LOOP < (n); ZOLTAN_ID_LOOP++)    \
      (a)[ZOLTAN_ID_LOOP] = (b)[ZOLTAN_ID_LOOP];                        \
   }

/****************************************************************************/
/*
 *  Prototypes for ID functions in id.c
 */

extern ZOLTAN_ID_PTR ZOLTAN_Malloc_ID(int n, char *file, int line);
extern void ZOLTAN_PRINT_ID(int n, ZOLTAN_ID_PTR a);
extern int ZOLTAN_EQ_ID(int n, ZOLTAN_ID_PTR a, ZOLTAN_ID_PTR b);

#ifdef ZOLTAN_NEEDED
/* Commented out since never used */
extern int ZOLTAN_LT_ID(int n, ZOLTAN_ID_PTR a, ZOLTAN_ID_PTR b);
extern int ZOLTAN_GT_ID(int n, ZOLTAN_ID_PTR a, ZOLTAN_ID_PTR b);
#endif  /* ZOLTAN_NEEDED */

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
