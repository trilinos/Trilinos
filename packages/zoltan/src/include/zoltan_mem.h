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
