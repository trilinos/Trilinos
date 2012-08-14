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

#ifndef __ZOLTAN_ID_CONST_H
#define __ZOLTAN_ID_CONST_H

#include "zoltan_id.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/*
 *  Macros that access the Zoltan ID manipulation functions.
 *  These macros assume a Zoltan_Struct is available, and access the
 *  appropriate members of Zoltan_Struct in calling the manipulation functions.
 *  Macros are provided for both global IDs (GIDs) and local IDs (LIDs).
 *
 *  The manipulation functions can be called directly; these macros just
 *  make it easier when a Zoltan_Struct is available.
 */

/****************************************************************************/

#define ZOLTAN_INIT_GID(zz,gid) ZOLTAN_INIT_ID((zz)->Num_GID,gid);
#define ZOLTAN_INIT_LID(zz,lid) ZOLTAN_INIT_ID((zz)->Num_LID,lid);

/*
 * Macros for allocating single IDs; IDs are also initialized.
 * If zz->Num_LID is zero, the macro returns NULL.
 */
#define ZOLTAN_MALLOC_GID(zz) ZOLTAN_Malloc_ID((zz)->Num_GID, __FILE__, __LINE__)
#define ZOLTAN_MALLOC_LID(zz) \
    ((zz)->Num_LID \
      ? ZOLTAN_Malloc_ID((zz)->Num_LID, __FILE__, __LINE__) \
      : NULL)

/*
 * Macros for allocating arrays of IDs; arrays are also initialized.
 * If zz->Num_LID is zero, the macro returns NULL.
 */
#define ZOLTAN_MALLOC_GID_ARRAY(zz,num_obj) \
    ZOLTAN_Malloc_ID((num_obj) * (zz)->Num_GID, __FILE__, __LINE__)
#define ZOLTAN_MALLOC_LID_ARRAY(zz,num_obj) \
    ((zz)->Num_LID \
       ? ZOLTAN_Malloc_ID((num_obj) * (zz)->Num_LID, __FILE__, __LINE__) \
       : NULL)

#define ZOLTAN_COPY_GID_ARRAY(to, from, zz,num_obj) \
  memcpy(to, from, (num_obj) * (zz)->Num_GID * sizeof(ZOLTAN_ID_TYPE));

#define ZOLTAN_COPY_LID_ARRAY(to, from, zz, num_obj) \
  if ((zz)->Num_LID) { \
     memcpy(to, from, (num_obj) * (zz)->Num_LID * sizeof(ZOLTAN_ID_TYPE)); \
  }


/*
 * Macros for reallocating arrays of IDs.
 */
#define ZOLTAN_REALLOC_GID_ARRAY(zz,ptr,num_new_obj) \
  (ZOLTAN_ID_PTR) ZOLTAN_REALLOC(ptr,\
               (num_new_obj)*(zz)->Num_GID*sizeof(ZOLTAN_ID_TYPE))

#define ZOLTAN_REALLOC_LID_ARRAY(zz,ptr,num_new_obj) \
  ((zz)->Num_LID \
    ? (ZOLTAN_ID_PTR)ZOLTAN_REALLOC(ptr, \
                               (num_new_obj)*(zz)->Num_LID*sizeof(ZOLTAN_ID_TYPE)) \
    : NULL)

/****************************************************************************/
/*
 *  Macros to copy IDs.
 */
#define ZOLTAN_SET_GID(zz,a,b) ZOLTAN_SET_ID((zz)->Num_GID,a,b)
#define ZOLTAN_SET_LID(zz,a,b) ZOLTAN_SET_ID((zz)->Num_LID,a,b)


/****************************************************************************/
/*
 * Macros to print IDs.
 */
#define ZOLTAN_PRINT_GID(zz,a) ZOLTAN_PRINT_ID((zz)->Num_GID,a)
#define ZOLTAN_PRINT_LID(zz,a) ZOLTAN_PRINT_ID((zz)->Num_LID,a)

/****************************************************************************/
/*
 * Macros to compare global IDs. (Comparisons of local IDs are not
 * needed as Zoltan only copies these IDs; it does not use them
 * in its computations.)
 */
#define ZOLTAN_EQ_GID(zz,a,b) ZOLTAN_EQ_ID((zz)->Num_GID,a,b)

#ifdef ZOLTAN_NEEDED
/* Commented out since never used */
#define ZOLTAN_LT_GID(zz,a,b) ZOLTAN_LT_ID((zz)->Num_GID,a,b)
#define ZOLTAN_GT_GID(zz,a,b) ZOLTAN_GT_ID((zz)->Num_GID,a,b)
#endif

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
