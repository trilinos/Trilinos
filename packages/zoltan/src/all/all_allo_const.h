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

extern int LB_Set_Malloc_Param(char *, char *);
extern void LB_Free_Structure(LB *);

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

/* type selector for LB_Special_Malloc */

enum LB_Special_Malloc_Type {
  LB_SPECIAL_MALLOC_INT,
  LB_SPECIAL_MALLOC_GID,
  LB_SPECIAL_MALLOC_LID
};

typedef enum LB_Special_Malloc_Type LB_SPECIAL_MALLOC_TYPE;

/* function declarations for special malloc */

extern int LB_Special_Malloc(LB *lb, void **array, int size,
                      LB_SPECIAL_MALLOC_TYPE type);
extern int LB_Special_Free(LB *lb, void **array,
                      LB_SPECIAL_MALLOC_TYPE type);
extern void Zoltan_Register_Fort_Malloc(ZOLTAN_FORT_MALLOC_INT_FN *,
                                        ZOLTAN_FORT_FREE_INT_FN *);

#endif
