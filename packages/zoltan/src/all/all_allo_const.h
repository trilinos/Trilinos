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


/* function declarations for dynamic array allocation */

extern double *LB_array_alloc(char *file, int lineno, int numdim, ...);
extern void LB_safe_free(void **ptr);
extern double *LB_smalloc(int n);

#endif
