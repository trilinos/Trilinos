/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 * $Name$
 *====================================================================*/

#ifndef __ALL_ALLO_H
#define __ALL_ALLO_H

#ifndef lint
static char *cvs_all_allo_h =
  "$Id$";
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

#define LB_MALLOC(a) LB_Malloc((a), __FILE__, __LINE__)
#define LB_REALLOC(a, b) LB_Realloc((a), (b), __FILE__, __LINE__)

/* function declarations for dynamic array allocation */

#ifdef __STDC__
extern double *LB_Array_Alloc(char *file, int lineno, int numdim, ...);
#else
extern double *LB_Array_Alloc();
#endif
extern int LB_Set_Malloc_Param(char *, char *);
extern void LB_Free(void **ptr);
extern double *LB_Malloc(int n, char *file, int lineno);
extern double *LB_Realloc(void *ptr, int n, char *filename, int lineno);
extern void LB_Memory_Stats();
extern int LB_Memory_Num();

#endif
