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

#ifndef __MEM_CONST_H
#define __MEM_CONST_H

#ifdef __STDC__
#include <string.h>
#else
#include <strings.h>
#endif  /* __STDC__ */

#ifndef HAVE_PROTOTYPES
#   if defined(__STDC__) || defined(__GNUC__) || defined(__cplusplus) || defined(c_plusplus)
#       define	HAVE_PROTOTYPES
#   endif
#endif

#undef PROTO
#ifdef HAVE_PROTOTYPES
#   define	PROTO(x)	x
#else
#   define	PROTO(x)	()
#endif

#define LB_MALLOC(a) LB_Malloc((a), __FILE__, __LINE__)
#define LB_REALLOC(a, b) LB_Realloc((a), (b), __FILE__, __LINE__)
#define LB_FREE(a) LB_Free((void **) (a), __FILE__, __LINE__)

/* function declarations for dynamic array allocation */

#ifdef __STDC__
extern double *LB_Array_Alloc(char *file, int lineno, int numdim, ...);
#else
extern double *LB_Array_Alloc();
#endif

extern void LB_Set_Memory_Debug(int);
extern void LB_Free(void **ptr, char *file, int lineno);
extern double *LB_Malloc(int n, char *file, int lineno);
extern double *LB_Realloc(void *ptr, int n, char *filename, int lineno);
extern void LB_Memory_Stats();
extern int LB_Memory_Num();

#ifdef __STDC__
extern void LB_Multifree(int n, ...);
#else
extern void LB_Multifree();
#endif

#endif
