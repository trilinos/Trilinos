// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


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
