/* ========================================================================== */
/* === Include/cholmod_internal.h =========================================== */
/* ========================================================================== */

/*
 * CHOLMOD version 0.1.  May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Authors: Timothy A. Davis,
 * William W. Hager, and Sivasankaran Rajamanickam.  Note that each of
 * CHOLMOD's modules are licensed separately.  The GNU LGPL license applies to
 * the Core, Check, and Partition modules.  The Cholesky, MatrixOps, Modify,
 * and Supernodal modules are not yet released; the GNU LGPL license will NOT
 * directly apply to them.  Since this file is required by all modules,
 * the GNU LGPL does apply to this file.
 * http://www.cise.ufl.edu/research/sparse
 */

/* CHOLMOD internal include file.
 *
 * Internal definitions for CHOLMOD, not meant to be included in user code,
 * since they define macros that are not prefixed with CHOLMOD_.  This file
 * can safely #include'd in user code if you want to make use of the macros
 * defined here, and don't mind the possible name conflicts with your code,
 * however.
 *
 * Required by all CHOLMOD routines.  Not required by any user routine that
 * uses CHOLMOMD.  Normally, does not require any CHOLMOD module (not even the
 * Core module).
 *
 * If debugging is enabled, all CHOLMOD modules require the Check module.
 * Enabling debugging requires that this file be editted.  Debugging cannot be
 * enabled with a compiler flag.  This is because CHOLMOD is exceedingly slow
 * when debugging is enabled.  Debugging is meant for development of CHOLMOD
 * itself, not by users of CHOLMOD.
 */

#ifndef CHOLMOD_INTERNAL_H
#define CHOLMOD_INTERNAL_H

/* turn off debugging */
#ifndef NDEBUG
#define NDEBUG
#endif

/* Uncomment this line to enable debugging.  CHOLMOD will be very slow.
#undef NDEBUG
 */

/* TODO: remove when done: */
#define GOTCHA(s,done) exit (((int *) NULL) [0] = 0)

/* Uncomment this line to enable METIS memory testing. (debugging only)
#define TEST_MEMORY
 */

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#include "matrix.h"
#else
#include <assert.h>
#endif

#include <stddef.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <stdlib.h>

/* Some non-conforming compilers insist on defining TRUE and FALSE. */
#undef TRUE
#undef FALSE
#define TRUE 1
#define FALSE 0

/* NULL should already be defined, but ensure it is here. */
#ifndef NULL
#define NULL ((void *) 0)
#endif

/* FLIP is a "negation about -1", and is used to mark an integer i that is
 * normally non-negative.  FLIP (EMPTY) is EMPTY.  FLIP of a number > EMPTY
 * is negative, and FLIP of a number < EMTPY is positive.  FLIP (FLIP (i)) = i
 * for all integers i.  UNFLIP (i) is >= EMPTY. */
#define EMPTY (-1)
#define FLIP(i) (-(i)-2)
#define UNFLIP(i) (((i) < EMPTY) ? FLIP (i) : (i))

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MAX3(a,b,c) (((a) > (b)) ? (MAX (a,c)) : (MAX (b,c)))
#define MAX4(a,b,c,d) (((a) > (b)) ? (MAX3 (a,c,d)) : (MAX3 (b,c,d)))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define IMPLIES(p,q) (!(p) || (q))

/* round up an integer x to a multiple of s */
#define ROUNDUP(x,s) ((s) * (((x) + ((s) - 1)) / (s)))

/* find the sign: -1 if x < 0, 1 if x > 0, zero otherwise.  Note that
 * IEEE -0 and NaN are treated as +0 */
#define SIGN(x) (((x) < 0) ? (-1) : (((x) > 0) ? 1 : 0))

/* This test for NaN's conforms to the IEEE 754, but at least one version of
 * the Microsoft Visual C/C++ breaks test. */
#define ISNAN(x) ((x) != (x))

/* Check a pointer and return if null.  Set status to invalid, unless the
 * status is already "out of memory" */
#define RETURN_IF_NULL(A,result) \
{ \
    if ((A) == NULL) \
    { \
	if (Common->status != CHOLMOD_OUT_OF_MEMORY) \
	{ \
	    cholmod_error (CHOLMOD_INVALID, "argument missing", Common) ; \
	} \
	return (result) ; \
    } \
}

/* Return if Common is NULL. No error code can be returned in Common->status. */
#define RETURN_IF_NULL_COMMON(result) \
{ \
    if (Common == NULL) \
    { \
	return (result) ; \
    } \
}

/* index to one entry past the end of the space available for column j of L */
#define COLEND(j) (L->ftype == CHOLMOD_LDL_DYNAMIC ? \
	(Lp [Lnext [j]]) : (((j) < n-1) ? Lp [(j)+1] : ((int) L->nzmax)))

/* ========================================================================== */
/* === debugging definitions ================================================ */
/* ========================================================================== */

#ifndef NDEBUG

#include <stdio.h>
#include "cholmod.h"

/* The cholmod_dump routines are in the Check module.  No CHOLMOD routine
 * calls the cholmod_check_* or cholmod_print_* routines in the Check module,
 * since they use Common workspace that may already be in use.  Instead, they
 * use the cholmod_dump_* routines defined there, which allocate their own
 * workspace if they need it. */
extern int cholmod_dump, cholmod_dump_malloc ;
long cholmod_dump_sparse    (cholmod_sparse  *, char *, cholmod_common *) ;
int  cholmod_dump_factor    (cholmod_factor  *, char *, cholmod_common *) ;
int  cholmod_dump_triplet   (cholmod_triplet *, char *, cholmod_common *) ;
int  cholmod_dump_dense     (cholmod_dense   *, char *, cholmod_common *) ;
int  cholmod_dump_subset    (void *, size_t, size_t, char *, cholmod_common *) ;
int  cholmod_dump_perm      (void *, size_t, size_t, char *, cholmod_common *) ;
int  cholmod_dump_parent    (void *, size_t, char *, cholmod_common *) ;
void cholmod_dump_init      (char *) ;
void cholmod_dump_real      (char *, double *, int, int) ;
void cholmod_dump_col       (char *, int, int, int, int *, double *, int) ;
void cholmod_dump_super     (int, int *, int *, int *, int *, double *) ;
int  cholmod_dump_mem       (char *, int, cholmod_common *) ;
int  cholmod_dump_mdelta    (int, int) ;
int  cholmod_dump_partition (int, int *, int *, int *, int *, int) ;
int  cholmod_dump_work	    (int, int, int, cholmod_common *) ;
void cholmod_dump_getwork   (size_t) ;
void cholmod_dump_freework  (void) ;

#define DEBUG_INIT(s)  { cholmod_dump_init(s) ; }

#ifdef MATLAB_MEX_FILE
#define ASSERT(expression) (mxAssert ((expression), ""))
#define PRINTF mexPrintf
#else
#define ASSERT(expression) { assert (expression) ; }
#define PRINTF printf
#endif

#define PRINT0(params) { (void) PRINTF params ; fflush (stdout) ; }
#define PRINT1(params) { if (cholmod_dump >= 1) (void) PRINTF params ; fflush (stdout) ; }
#define PRINT2(params) { if (cholmod_dump >= 2) (void) PRINTF params ; fflush (stdout) ; }
#define PRINT3(params) { if (cholmod_dump >= 3) (void) PRINTF params ; fflush (stdout) ; }
#define PRINTM(params) { if (cholmod_dump_malloc) (void) PRINTF params ; fflush (stdout) ; }
#define DEBUG(statement) statement

#else

/* Debugging disabled (the normal case) */
#define DEBUG_INIT(s)
#define PRINT0(params)
#define PRINT1(params)
#define PRINT2(params)
#define PRINT3(params)
#define PRINTM(params)
#define ASSERT(expression)
#define DEBUG(statement)
#endif

#endif
