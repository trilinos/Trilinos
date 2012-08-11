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


#ifndef __ALL_ALLO_H
#define __ALL_ALLO_H

#include "zz_const.h"
#include "zoltan_mem.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


extern int Zoltan_Set_Malloc_Param(char *, char *);

/* function prototypes for Fortran allocation functions */

#ifdef PGI
typedef void ZOLTAN_FORT_MALLOC_INT_FN(int *arg, int *size, int **ret, int *hidden);
typedef void ZOLTAN_FORT_FREE_INT_FN(int *arg, int *hidden);
#else
#ifdef FUJITSU
typedef void ZOLTAN_FORT_MALLOC_INT_FN(int *arg, int *size, int **ret, int *hidden, int hidden2, int hidden3);
typedef void ZOLTAN_FORT_FREE_INT_FN(int *arg, int *hidden);
#else
typedef void ZOLTAN_FORT_MALLOC_INT_FN(int *arg, int *size, int **ret);
typedef void ZOLTAN_FORT_FREE_INT_FN(int *arg);
#endif
#endif

/* these are not as complicated - just a simple typedef for all
   compilers, since we are not passing any Fortran structures */
typedef void ZOLTAN_FORT_MALLOC_SET_STRUCT_FN(int *arg, int **ret);

/* type selector for Zoltan_Special_Malloc */

enum Zoltan_Special_Malloc_Type {
  ZOLTAN_SPECIAL_MALLOC_INT,
  ZOLTAN_SPECIAL_MALLOC_GID,
  ZOLTAN_SPECIAL_MALLOC_LID
};

typedef enum Zoltan_Special_Malloc_Type ZOLTAN_SPECIAL_MALLOC_TYPE;

/* function declarations for special malloc */

extern int Zoltan_Special_Malloc(ZZ *zz, void **array, int size,
                      ZOLTAN_SPECIAL_MALLOC_TYPE type);
extern int Zoltan_Special_Free(ZZ *zz, void **array,
                      ZOLTAN_SPECIAL_MALLOC_TYPE type);
extern void Zoltan_Register_Fort_Malloc(ZOLTAN_FORT_MALLOC_INT_FN *,
                                        ZOLTAN_FORT_FREE_INT_FN *,
					ZOLTAN_FORT_MALLOC_SET_STRUCT_FN *);
extern int Zoltan_Special_Fort_Malloc_Set_Struct(int *zz_addr_bytes, 
						 int **fort_zz);
extern int Zoltan_Special_Fort_Free_Struct(int *fort_zz);


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
