/* ========================================================================== */
/* === klu_kernel.h ========================================================= */
/* ========================================================================== */
#include "limits.h"

/* This file should not be included in any user routine. */

typedef struct
{
    int j ;	    /* column j of L */
    int p ;	    /* Li [p..pend-1] currently being traversed */
    int pend ;

} WorkStackType ;


int klu_prune
(
    /* input, not modified */
    int n,	    /* A is n-by-n */
    int Ap [ ],	    /* size n+1, column pointers for A */
    int Ai [ ],	    /* size nz = Ap [n], row indices for A */
    double Ax [ ],  /* size nz, values of A */
    double tol,	    /* partial pivoting tolerance parameter */

    /* input and output */
    int *p_lsize,
    int *p_usize,

    /* output, not defined on input */
    int Lp [ ],	    /* size n+1 */
    int **p_Li,	    /* size lsize */
    double **p_Lx,  /* size lsize */
    int Up [ ],	    /* size n+1 */
    int **p_Ui,	    /* size usize */
    double **p_Ux,  /* size usize */
    int Pinv [ ],   /* size n */
    int P [ ],	    /* size n */

    /* workspace, not defined on input or output */
    double X [ ],   /* size n */
    int Stack [ ],  /* size n */

    /* workspace for non-recursive version only */
    WorkStackType WorkStack [ ],    /* size n */

    /* workspace for pruning only */
    int Lpend [ ],	/* size n workspace */
    int Lpruned [ ]	/* size n workspace */
) ;

int klu_noprune
(
    /* input, not modified */
    int n,	    /* A is n-by-n */
    int Ap [ ],	    /* size n+1, column pointers for A */
    int Ai [ ],	    /* size nz = Ap [n], row indices for A */
    double Ax [ ],  /* size nz, values of A */
    double tol,	    /* partial pivoting tolerance parameter */

    /* input and output */
    int *p_lsize,
    int *p_usize,

    /* output, not defined on input */
    int Lp [ ],	    /* size n+1 */
    int **p_Li,	    /* size lsize */
    double **p_Lx,  /* size lsize */
    int Up [ ],	    /* size n+1 */
    int **p_Ui,	    /* size usize */
    double **p_Ux,  /* size usize */
    int Pinv [ ],   /* size n */
    int P [ ],	    /* size n */

    /* workspace, not defined on input or output */
    double X [ ],   /* size n */
    int Stack [ ],  /* size n */

    /* workspace for non-recursive version only */
    WorkStackType WorkStack [ ],    /* size n */

    /* workspace for pruning only */
    int Lpend [ ],	/* size n workspace */
    int Lpruned [ ]	/* size n workspace */
) ;

int klu_prune_nonrecursive
(
    /* input, not modified */
    int n,	    /* A is n-by-n */
    int Ap [ ],	    /* size n+1, column pointers for A */
    int Ai [ ],	    /* size nz = Ap [n], row indices for A */
    double Ax [ ],  /* size nz, values of A */
    double tol,	    /* partial pivoting tolerance parameter */

    /* input and output */
    int *p_lsize,
    int *p_usize,

    /* output, not defined on input */
    int Lp [ ],	    /* size n+1 */
    int **p_Li,	    /* size lsize */
    double **p_Lx,  /* size lsize */
    int Up [ ],	    /* size n+1 */
    int **p_Ui,	    /* size usize */
    double **p_Ux,  /* size usize */
    int Pinv [ ],   /* size n */
    int P [ ],	    /* size n */

    /* workspace, not defined on input or output */
    double X [ ],   /* size n */
    int Stack [ ],  /* size n */

    /* workspace for non-recursive version only */
    WorkStackType WorkStack [ ],    /* size n */

    /* workspace for pruning only */
    int Lpend [ ],	/* size n workspace */
    int Lpruned [ ]	/* size n workspace */
) ;

int klu_noprune_nonrecursive
(
    /* input, not modified */
    int n,	    /* A is n-by-n */
    int Ap [ ],	    /* size n+1, column pointers for A */
    int Ai [ ],	    /* size nz = Ap [n], row indices for A */
    double Ax [ ],  /* size nz, values of A */
    double tol,	    /* partial pivoting tolerance parameter */

    /* input and output */
    int *p_lsize,
    int *p_usize,

    /* output, not defined on input */
    int Lp [ ],	    /* size n+1 */
    int **p_Li,	    /* size lsize */
    double **p_Lx,  /* size lsize */
    int Up [ ],	    /* size n+1 */
    int **p_Ui,	    /* size usize */
    double **p_Ux,  /* size usize */
    int Pinv [ ],   /* size n */
    int P [ ],	    /* size n */

    /* workspace, not defined on input or output */
    double X [ ],   /* size n */
    int Stack [ ],  /* size n */

    /* workspace for non-recursive version only */
    WorkStackType WorkStack [ ],    /* size n */

    /* workspace for pruning only */
    int Lpend [ ],	/* size n workspace */
    int Lpruned [ ]	/* size n workspace */
) ;

#define NDEBUG
/* To enable debugging, uncomment this line:
#undef NDEBUG
*/

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
#define ASSERT(a) assert(a)
#define ALLOCATE malloc
#define _REALLOC realloc
#define _FREE free
#endif

#include <stdlib.h>

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
	else \
	{ \
	    printf ("failed\n") ; \
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

#ifndef NDEBUG
#define PRINTF(s) { printf s ; } ;
#else
#define PRINTF(s)
#endif

#define TRUE 1
#define FALSE 0
#define MAX(a,b) (((a) > (b)) ?  (a) : (b))
#define MIN(a,b) (((a) < (b)) ?  (a) : (b))
#define ABS(a)   (((a) <  0 ) ? -(a) : (a))
