/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __ALL_ALLO_H
#define __ALL_ALLO_H

#include "lb_const.h"
#include "mem_const.h"

extern int Zoltan_Set_Malloc_Param(char *, char *);
extern void Zoltan_Free_Structure(ZZ *);

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
                                        ZOLTAN_FORT_FREE_INT_FN *);

#endif
