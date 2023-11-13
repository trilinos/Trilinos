/* ========================================================================= */
/* === trilinos_amd_internal.h ====================================================== */
/* ========================================================================= */

/* ------------------------------------------------------------------------- */
/* AMD, Copyright (c) Timothy A. Davis,					     */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README.txt for License.     */
/* email: davis at cise.ufl.edu    CISE Department, Univ. of Florida.        */
/* web: http://www.cise.ufl.edu/research/sparse/amd                          */
/* ------------------------------------------------------------------------- */

/* This file is for internal use in AMD itself, and does not normally need to
 * be included in user code (it is included in UMFPACK, however).   All others
 * should use amd.h instead.
 *
 * The following compile-time definitions affect how AMD is compiled.
 *
 *	-DNPRINT
 *
 *	    Disable all printing.  stdio.h will not be included.  Printing can
 *	    be re-enabled at run-time by setting the global pointer trilinos_amd_printf
 *	    to printf (or mexPrintf for a MATLAB mexFunction).
 *
 *	-DNMALLOC
 *
 *	    No memory manager is defined at compile-time.  You MUST define the
 *	    function pointers trilinos_amd_malloc, trilinos_amd_free, trilinos_amd_realloc, and
 *	    trilinos_amd_calloc at run-time for AMD to work properly.
 */

/* ========================================================================= */
/* === NDEBUG ============================================================== */
/* ========================================================================= */

/*
 * Turning on debugging takes some work (see below).   If you do not edit this
 * file, then debugging is always turned off, regardless of whether or not
 * -DNDEBUG is specified in your compiler options.
 *
 * If AMD is being compiled as a mexFunction, then MATLAB_MEX_FILE is defined,
 * and mxAssert is used instead of assert.  If debugging is not enabled, no
 * MATLAB include files or functions are used.  Thus, the AMD library libamd.a
 * can be safely used in either a stand-alone C program or in another
 * mexFunction, without any change.
 */

/*
    AMD will be exceedingly slow when running in debug mode.  The next three
    lines ensure that debugging is turned off.
*/
#ifndef NDEBUG
#define NDEBUG
#endif

/*
    To enable debugging, uncomment the following line:
#undef NDEBUG
*/

/* ------------------------------------------------------------------------- */
/* ANSI include files */
/* ------------------------------------------------------------------------- */

/* from stdlib.h:  size_t, malloc, free, realloc, and calloc */
#include <stdlib.h>

#if !defined(NPRINT) || !defined(NDEBUG)
/* from stdio.h:  printf.  Not included if NPRINT is defined at compile time.
 * fopen and fscanf are used when debugging. */
#include <stdio.h>
#endif

/* from limits.h:  INT_MAX and LONG_MAX */
#include <limits.h>

/* from math.h: sqrt */
#include <math.h>

/* ------------------------------------------------------------------------- */
/* MATLAB include files (only if being used in or via MATLAB) */
/* ------------------------------------------------------------------------- */

#ifdef MATLAB_MEX_FILE
#include "matrix.h"
#include "mex.h"
#endif

/* ------------------------------------------------------------------------- */
/* basic definitions */
/* ------------------------------------------------------------------------- */

#ifdef FLIP
#undef FLIP
#endif

#ifdef MAX
#undef MAX
#endif

#ifdef MIN
#undef MIN
#endif

#ifdef TRILINOS_AMD_EMPTY
#undef TRILINOS_AMD_EMPTY
#endif

#ifdef GLOBAL
#undef GLOBAL
#endif

#ifdef PRIVATE
#undef PRIVATE
#endif

/* FLIP is a "negation about -1", and is used to mark an integer i that is
 * normally non-negative.  FLIP (TRILINOS_AMD_EMPTY) is TRILINOS_AMD_EMPTY.  FLIP of a number > TRILINOS_AMD_EMPTY
 * is negative, and FLIP of a number < EMTPY is positive.  FLIP (FLIP (i)) = i
 * for all integers i.  UNFLIP (i) is >= TRILINOS_AMD_EMPTY. */
#define TRILINOS_AMD_EMPTY (-1)
#define FLIP(i) (-(i)-2)
#define UNFLIP(i) ((i < TRILINOS_AMD_EMPTY) ? FLIP (i) : (i))

/* for integer MAX/MIN, or for doubles when we don't care how NaN's behave: */
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

/* logical expression of p implies q: */
#define IMPLIES(p,q) (!(p) || (q))

/* Note that the IBM RS 6000 xlc predefines TRUE and FALSE in <types.h>. */
/* The Compaq Alpha also predefines TRUE and FALSE. */
#ifdef TRUE
#undef TRUE
#endif
#ifdef FALSE
#undef FALSE
#endif

#define TRUE (1)
#define FALSE (0)
#define PRIVATE static
#define GLOBAL
#define TRILINOS_AMD_EMPTY (-1)

/* Note that Linux's gcc 2.96 defines NULL as ((void *) 0), but other */
/* compilers (even gcc 2.95.2 on Solaris) define NULL as 0 or (0).  We */
/* need to use the ANSI standard value of 0. */
#ifdef NULL
#undef NULL
#endif

#define NULL 0

/* largest value of size_t */
#ifndef SIZE_T_MAX
#define SIZE_T_MAX ((size_t) (-1))
#endif

/* ------------------------------------------------------------------------- */
/* integer type for AMD: int or UF_long */
/* ------------------------------------------------------------------------- */

/* define UF_long */
#include "trilinos_UFconfig.h"

#if defined (DLONG) || defined (ZLONG)

#define Int UF_long
#define ID  UF_long_id
#define Int_MAX UF_long_max

#define TRILINOS_AMD_order trilinos_amd_l_order
#define TRILINOS_AMD_defaults trilinos_amd_l_defaults
#define TRILINOS_AMD_control trilinos_amd_l_control
#define TRILINOS_AMD_info trilinos_amd_l_info
#define TRILINOS_AMD_1 trilinos_amd_l1
#define TRILINOS_AMD_2 trilinos_amd_l2
#define TRILINOS_AMD_valid trilinos_amd_l_valid
#define TRILINOS_AMD_aat trilinos_amd_l_aat
#define TRILINOS_AMD_postorder trilinos_amd_l_postorder
#define TRILINOS_AMD_post_tree trilinos_amd_l_post_tree
#define TRILINOS_AMD_dump trilinos_amd_l_dump
#define TRILINOS_AMD_debug trilinos_amd_l_debug
#define TRILINOS_AMD_debug_init trilinos_amd_l_debug_init
#define TRILINOS_AMD_preprocess trilinos_amd_l_preprocess

#else

#define Int int
#define ID "%d"
#define Int_MAX INT_MAX

#define TRILINOS_AMD_order trilinos_amd_order
#define TRILINOS_AMD_defaults trilinos_amd_defaults
#define TRILINOS_AMD_control trilinos_amd_control
#define TRILINOS_AMD_info trilinos_amd_info
#define TRILINOS_AMD_1 trilinos_amd_1
#define TRILINOS_AMD_2 trilinos_amd_2
#define TRILINOS_AMD_valid trilinos_amd_valid
#define TRILINOS_AMD_aat trilinos_amd_aat
#define TRILINOS_AMD_postorder trilinos_amd_postorder
#define TRILINOS_AMD_post_tree trilinos_amd_post_tree
#define TRILINOS_AMD_dump trilinos_amd_dump
#define TRILINOS_AMD_debug trilinos_amd_debug
#define TRILINOS_AMD_debug_init trilinos_amd_debug_init
#define TRILINOS_AMD_preprocess trilinos_amd_preprocess

#endif

/* ========================================================================= */
/* === PRINTF macro ======================================================== */
/* ========================================================================= */

/* All output goes through the PRINTF macro.  */
#define PRINTF(params) { if (trilinos_amd_printf != NULL) (void) trilinos_amd_printf params ; }

/* ------------------------------------------------------------------------- */
/* AMD routine definitions (user-callable) */
/* ------------------------------------------------------------------------- */

#include "trilinos_amd.h"

/* ------------------------------------------------------------------------- */
/* AMD routine definitions (not user-callable) */
/* ------------------------------------------------------------------------- */

GLOBAL size_t TRILINOS_AMD_aat
(
    Int n,
    const Int Ap [ ],
    const Int Ai [ ],
    Int Len [ ],
    Int Tp [ ],
    double Info [ ]
) ;

GLOBAL void TRILINOS_AMD_1
(
    Int n,
    const Int Ap [ ],
    const Int Ai [ ],
    Int P [ ],
    Int Pinv [ ],
    Int Len [ ],
    Int slen,
    Int S [ ],
    double Control [ ],
    double Info [ ]
) ;

GLOBAL void TRILINOS_AMD_postorder
(
    Int nn,
    Int Parent [ ],
    Int Npiv [ ],
    Int Fsize [ ],
    Int Order [ ],
    Int Child [ ],
    Int Sibling [ ],
    Int Stack [ ]
) ;

GLOBAL Int TRILINOS_AMD_post_tree
(
    Int root,
    Int k,
    Int Child [ ],
    const Int Sibling [ ],
    Int Order [ ],
    Int Stack [ ]
#ifndef NDEBUG
    , Int nn
#endif
) ;

GLOBAL void TRILINOS_AMD_preprocess
(
    Int n,
    const Int Ap [ ],
    const Int Ai [ ],
    Int Rp [ ],
    Int Ri [ ],
    Int W [ ],
    Int Flag [ ]
) ;

/* ------------------------------------------------------------------------- */
/* debugging definitions */
/* ------------------------------------------------------------------------- */

#ifndef NDEBUG

/* from assert.h:  assert macro */
#include <assert.h>

#ifndef EXTERN
#define EXTERN extern
#endif

EXTERN Int TRILINOS_AMD_debug ;

GLOBAL void TRILINOS_AMD_debug_init ( char *s ) ;

GLOBAL void TRILINOS_AMD_dump
(
    Int n,
    Int Pe [ ],
    Int Iw [ ],
    Int Len [ ],
    Int iwlen,
    Int pfree,
    Int Nv [ ],
    Int Next [ ],
    Int Last [ ],
    Int Head [ ],
    Int Elen [ ],
    Int Degree [ ],
    Int W [ ],
    Int nel
) ;

#ifdef ASSERT
#undef ASSERT
#endif

/* Use mxAssert if AMD is compiled into a mexFunction */
#ifdef MATLAB_MEX_FILE
#define ASSERT(expression) (mxAssert ((expression), ""))
#else
#define ASSERT(expression) (assert (expression))
#endif

#define TRILINOS_AMD_DEBUG0(params) { PRINTF (params) ; }
#define TRILINOS_AMD_DEBUG1(params) { if (TRILINOS_AMD_debug >= 1) PRINTF (params) ; }
#define TRILINOS_AMD_DEBUG2(params) { if (TRILINOS_AMD_debug >= 2) PRINTF (params) ; }
#define TRILINOS_AMD_DEBUG3(params) { if (TRILINOS_AMD_debug >= 3) PRINTF (params) ; }
#define TRILINOS_AMD_DEBUG4(params) { if (TRILINOS_AMD_debug >= 4) PRINTF (params) ; }

#else

/* no debugging */
#define ASSERT(expression)
#define TRILINOS_AMD_DEBUG0(params)
#define TRILINOS_AMD_DEBUG1(params)
#define TRILINOS_AMD_DEBUG2(params)
#define TRILINOS_AMD_DEBUG3(params)
#define TRILINOS_AMD_DEBUG4(params)

#endif
