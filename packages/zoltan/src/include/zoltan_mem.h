// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __MEM_CONST_H
#define __MEM_CONST_H


#include <string.h>

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


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

#define ZOLTAN_MALLOC(a)     Zoltan_Malloc((a), __FILE__, __LINE__)
#define ZOLTAN_CALLOC(a, b)  Zoltan_Calloc((a), (b), __FILE__, __LINE__)
#define ZOLTAN_FREE(a)       Zoltan_Free((void**)(void*) (a), __FILE__, __LINE__)

#define ZOLTAN_REALLOC(a, b) Zoltan_Realloc((a), (b),  __FILE__, __LINE__)

#define ZOLTAN_MEM_STAT_TOTAL   0
#define ZOLTAN_MEM_STAT_MAXIMUM 1

/* function declarations for dynamic array allocation */

#ifdef __STDC__
extern double *Zoltan_Array_Alloc(char *file, int lineno, int numdim, ...);
#else
extern double *Zoltan_Array_Alloc();
#endif

extern void    Zoltan_Memory_Debug(int);
extern int Zoltan_Memory_Get_Debug();
extern void    Zoltan_Free(void **, char *, int);
extern double *Zoltan_Calloc(size_t, size_t, char *, int);
extern double *Zoltan_Malloc(size_t, char *, int);
extern double *Zoltan_Realloc(void *, size_t, char *, int);
extern void    Zoltan_Memory_Stats(void);
extern size_t  Zoltan_Memory_Usage(int);
extern void    Zoltan_Memory_Reset(int);

#ifdef __STDC__
extern void Zoltan_Multifree(char *, int, int n, ...);
#else
extern void Zoltan_Multifree();
#endif

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
