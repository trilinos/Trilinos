/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */

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
