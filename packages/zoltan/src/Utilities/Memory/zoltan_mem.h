/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
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

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


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

#define ZOLTAN_MALLOC(a) Zoltan_Malloc((a), __FILE__, __LINE__)
#define ZOLTAN_REALLOC(a, b) Zoltan_Realloc((a), (b), __FILE__, __LINE__)
#define ZOLTAN_FREE(a) Zoltan_Free((void **) (a), __FILE__, __LINE__)

/* function declarations for dynamic array allocation */

#ifdef __STDC__
extern double *Zoltan_Array_Alloc(char *file, int lineno, int numdim, ...);
#else
extern double *Zoltan_Array_Alloc();
#endif

extern void Zoltan_Memory_Debug(int);
extern void Zoltan_Free(void **ptr, char *file, int lineno);
extern double *Zoltan_Malloc(int n, char *file, int lineno);
extern double *Zoltan_Realloc(void *ptr, int n, char *filename, int lineno);
extern void Zoltan_Memory_Stats();
extern int Zoltan_Memory_Num();

#ifdef __STDC__
extern void Zoltan_Multifree(int n, ...);
#else
extern void Zoltan_Multifree();
#endif

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
