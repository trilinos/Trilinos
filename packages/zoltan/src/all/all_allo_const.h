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

#define LB_SMALLOC(a) LB_smalloc((a), __FILE__, __LINE__)
#define LB_SREALLOC(a, b) LB_srealloc((a), (b), __FILE__, __LINE__)

/* function declarations for dynamic array allocation */

#ifdef __STDC__
extern double *LB_array_alloc(char *file, int lineno, int numdim, ...);
#else
extern double *LB_array_alloc();
#endif
extern int LB_Malloc_Set_Param(char *, char *);
extern void LB_safe_free(void **ptr);
extern double *LB_smalloc(int n, char *file, int lineno);
extern double *LB_srealloc(void *ptr, int n, char *filename, int lineno);

#endif
