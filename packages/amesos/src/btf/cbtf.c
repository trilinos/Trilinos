/* ========================================================================== */
/* === cbtf ================================================================= */
/* ========================================================================== */

/* cbtf is a user-callable C wrapper for permuting a matrix to block triangular
 * form. */

/* ========================================================================== */

/* make sure debugging is turned off */
#ifndef NDEBUG
#define NDEBUG
#endif
/* To enable debugging, uncomment this line: 
#undef NDEBUG
*/

#ifdef MATLAB_MEX_FILE
#include "matrix.h"
#include "mex.h"
#define ASSERT(a) mxAssert(a, "")
#define ALLOCATE mxMalloc
#define _FREE mxFree
#else
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#define ASSERT(a) assert(a)
#define ALLOCATE malloc
#define _FREE free
#endif

#include <stdlib.h>

#define FREE(p,type) \
    { \
	if (p != (type *) NULL) \
	{ \
	    _FREE (p) ; \
	    p = (type *) NULL ; \
	} \
    }

#include "cbtf.h"





/* the btf fortran interface to mc13 and mc21 */
int btf_ (int *Ai, int *Ap, int *nz, int *n, int *A1p,
    int *A1i, int *Perm, int *Lenc, int *Cp, int *Bp,
    int *Cperm, int *Rperm, int *Work, int *nblocks) ;

int cbtf_single_block
(
    int n,
    int P [ ],
    int Q [ ],
    int R [ ]
)
{
    /* return the BTF form as a single block (no BTF) */
    int k ;
    R [0] = 0 ;
    R [1] = n ;
    for (k = 0 ; k < n ; k++)
    {
	P [k] = k ;
	Q [k] = k ;
    }
    return (1) ;
}

int cbtf	/* return nblocks, or -1 if out of memory */
(
    /* inputs, not modified */
    int n,
    int Ap [ ],
    int Ai [ ],
    /* outputs, not defined on input */
    int P [ ],	/* size n, row permutation */
    int Q [ ],	/* size n, column permutation */
    int R [ ]	/* size n+1, R [b] = k if row/col k is the start
		   of block b.  R [nblocks] = n.  */
)
{
    int nz, *A1p, *A1i, *Perm, *Cperm, *Rperm, *Lenc, *Work, *Cp, *Bp,
	nblocks, k ;

    nblocks = 0 ;
    nz = Ap [n] ;

    if (nz == 0)
    {
	/* all zero matrix */
	return (cbtf_single_block (n, P, Q, R)) ;
    }

    /* allocate workspace */
    A1p = (int *) ALLOCATE ((n+1) * sizeof (int)) ;
    A1i = (int *) ALLOCATE ((nz+1) * sizeof (int)) ;
    Lenc  = (int *) ALLOCATE (n * sizeof (int)) ;
    Perm  = (int *) ALLOCATE (n * sizeof (int)) ;
    Cp    = (int *) ALLOCATE ((n+1) * sizeof (int)) ;
    Work  = (int *) ALLOCATE ((4*n) * sizeof (int)) ;

    /* allocate 1-based outputs */
    Cperm = (int *) ALLOCATE (n * sizeof (int)) ;
    Rperm = (int *) ALLOCATE (n * sizeof (int)) ;
    Bp    = (int *) ALLOCATE (n * sizeof (int)) ;

    if ((A1p == (int *) NULL) || (A1i == (int *) NULL) ||
	(Lenc == (int *) NULL) || (Perm == (int *) NULL) ||
	(Cp == (int *) NULL) || (Work == (int *) NULL) ||
	(Cperm == (int *) NULL) || (Rperm == (int *) NULL) ||
	(Bp == (int *) NULL))
    {
	/* out of memory */
	nblocks = -1 ;
    }

    if (nblocks != -1)
    {
	(void) btf_ (Ai, Ap, &nz, &n, A1p, A1i, Perm, Lenc, Cp, Bp,
	    Cperm, Rperm, Work, &nblocks) ;
    }

    /* free workspace */
    FREE (A1p, int) ;
    FREE (A1i, int) ;
    FREE (Lenc, int) ;
    FREE (Work, int) ;
    FREE (Cp, int) ;
    FREE (Perm, int) ;

    /* convert 1-based outputs to 0-based outputs (if successful) */
    if (nblocks == 0)
    {
	/* matrix is structurally singular */
	nblocks = cbtf_single_block (n, P, Q, R) ;
    }
    else if (nblocks > 0)
    {
	for (k = 0 ; k < n ; k++)
	{
	    P [k] = Rperm [k] - 1 ;
	    Q [k] = Cperm [k] - 1 ;
	}
	for (k = 0 ; k < nblocks ; k++)
	{
	    R [k] = Bp [k] - 1 ;
	}
	R [nblocks] = n ;
    }

    /* free the 1-based outputs */
    FREE (Cperm, int) ;
    FREE (Rperm, int) ;
    FREE (Bp, int) ;

    return (nblocks) ;
}
