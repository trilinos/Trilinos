/* ========================================================================== */
/* === klu_kernel.h ========================================================= */
/* ========================================================================== */

/* This file should not be included in any user routine. */

int klu_kernel
(
    /* input, not modified */
    int n,	    /* A is n-by-n */
    int Ap [ ],	    /* size n+1, column pointers for A */
    int Ai [ ],	    /* size nz = Ap [n], row indices for A */
    double Ax [ ],  /* size nz, values of A */
    int Q [ ],	    /* size n, optional input permutation */
    double tol,	    /* partial pivoting tolerance parameter */
    double growth,  /* memory growth factor */
    int lsize,	    /* initial size of Li and Lx */
    int usize,	    /* initial size of Ui and Ux */

    /* output, not defined on input */
    int Lp [ ],	    /* size n+1 */
    int **p_Li,	    /* size lsize */
    double **p_Lx,  /* size lsize */
    int Up [ ],	    /* size n+1 */
    int **p_Ui,	    /* size usize */
    double **p_Ux,  /* size usize */
    int Pinv [ ],   /* size n */
    int P [ ],	    /* size n */
    int *p_noffdiag,	/* # of off-diagonal pivots chosen */
    double *p_umin, /* smallest entry on the diagonal of U */
    double *p_umax, /* largest entry on the diagonal of U */

    /* workspace, not defined on input */
    double X [ ],   /* size n, zero on output */

    /* workspace, not defined on input or output */
    int Stack [ ],  /* size n */
    int Flag [ ],   /* size n */
    int adj_pos [ ],	/* size n */
    
    /* workspace for pruning only */
    int Lpend [ ]	/* size n workspace */
) ;

/* ========================================================================== */
/* make sure debugging and printing is turned off */
#ifndef NDEBUG
#define NDEBUG
#endif
#ifndef NPRINT
#define NPRINT
#endif
/* To enable debugging and assertions, uncomment this line: 
#undef NDEBUG
*/
/* To enable diagnostic printing, uncomment this line: 
#undef NPRINT
*/
/* ========================================================================== */

#ifdef MATLAB_MEX_FILE
#include "matrix.h"
#include "mex.h"
#define ASSERT(a) mxAssert(a, "")
#define ALLOCATE mxMalloc
#define _REALLOC mxRealloc
#define _FREE mxFree
#else
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#define ASSERT(a) assert(a)
#define ALLOCATE malloc
#define _REALLOC realloc
#define _FREE free
#endif

#include <stdlib.h>
#include <math.h>

#define REALLOCATE(p,type,size,ok) \
    { \
	type *pnew ; \
	size_t s ; \
	s = (size_t) ((sizeof (type)) * size) ; \
	pnew = (type *) _REALLOC ((void *) p, s) ; \
	ok = (pnew != (type *) NULL) ; \
	if (ok) \
	{ \
	    p = pnew ; \
	} \
    }

#define FREE(p,type) \
    { \
	if (p != (type *) NULL) \
	{ \
	    _FREE (p) ; \
	    p = (type *) NULL ; \
	} \
    }

#define SCALAR_IS_NAN(x) ((x) != (x))

/* true if an integer (stored in double x) would overflow (or if x is NaN) */
#define INT_OVERFLOW(x) ((!((x) * (1.0+1e-8) <= (double) INT_MAX)) \
			|| SCALAR_IS_NAN (x))

#undef TRUE
#undef FALSE
#undef MAX
#undef MIN
#undef ABS
#undef PRINTF
#undef FLIP

#ifndef NPRINT
#define PRINTF(s) { printf s ; } ;
#else
#define PRINTF(s)
#endif

#define TRUE 1
#define FALSE 0
#define MAX(a,b) (((a) > (b)) ?  (a) : (b))
#define MIN(a,b) (((a) < (b)) ?  (a) : (b))
#define ABS(a)   (((a) <  0 ) ? -(a) : (a))

/* FLIP is a "negation about -1", and is used to mark an integer i that is
 * normally non-negative.  FLIP (EMPTY) is EMPTY.  FLIP of a number > EMPTY
 * is negative, and FLIP of a number < EMTPY is positive.  FLIP (FLIP (i)) = i
 * for all integers i.  UNFLIP (i) is >= EMPTY. */
#define EMPTY (-1)
#define FLIP(i) (-(i)-2)
#define UNFLIP(i) ((i < EMPTY) ? FLIP (i) : (i))

/* get a parameter from the Control array */ 
#define GET_CONTROL(i,default) \
    ((Control != (double *) NULL) ? \
	(SCALAR_IS_NAN (Control [i]) ? default : Control [i]) \
	: default)

