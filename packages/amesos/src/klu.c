/* ========================================================================== */
/* === klu ================================================================== */
/* ========================================================================== */

/* KLU: factorizes P*A into L*U, using the Gilbert-Peierls algorithm, with
 * optional symmetric pruning by Eisenstat and Liu.  The code is by Tim Davis.
 * This algorithm is what appears as the default sparse LU routine in MATLAB
 * version 6.0, and still appears in MATLAB 6.5 as [L,U,P] = lu (A).  Note
 * that no column ordering is provided (see COLAMD or AMD for suitable
 * orderings).  SuperLU is based on this algorithm, except that it adds the
 * use of dense matrix operations on "supernodes" (adjacent columns with
 * identical).  This code doesn't use supernodes, thus its name (Kent LU,
 * as in Clark Kent, in contrast with Super-LU...).  This algorithm is slower
 * than SuperLU and UMFPACK for large matrices with lots of nonzeros in their
 * factors (such as for most finite-element problems).  However, for matrices
 * with very sparse LU factors, this algorithm is typically faster than both
 * SuperLU and UMFPACK, since in this case there is little chance to exploit
 * dense matrix kernels (the BLAS).
 *
 * NOTE: no error checking is done on the inputs.  This version is not yet in
 * pure library-quality shape.  It is a reasonable implementation, but a robust
 * code would provide a mechanism for verifying and printing its inputs and
 * outputs, use fill-reducing orderings, provide scaling options, and so on.
 * The 7 output arrays
 * should be allocated in one call to malloc, and then grown in size via
 * realloc.  This typically allows realloc to perform without copying.
 * UMFPACK includes all of the above features.
 *
 * No fill-reducing ordering is provided!  The ordering quality of klu is your
 * responsibility.  You must pre-permute A to reduce fill-in.  Two typical
 * scenarios (in MATLAB notation) are:
 *
 *	q = colamd (A) ;
 *	[L,U,P] = klu (A (:,q)) ;
 *
 * which returns L, U, and P such that L*U = P*A(:,q) or
 *
 *	q = symamd (A) ;	    % or q = amd (A)
 *	[L,U,P] = klu (A (q,q)) ;
 *
 * which returns L, U, and P such that L*U = P*A(q,q).  In the latter case, a
 * pivot tolerance of 1e-3 would be appropriate, to preserve symmetry.  See the
 * MATLAB interface klu.m for more details.
 *
 * The input matrix A must be in compressed-column form, with sorted row
 * indices.  Row indices for column j of A is in Ai [Ap [j] ... Ap [j+1]-1]
 * and the same range of indices in Ax holds the numerical values.
 * Control is ignored if it is (double *) NULL.  To set the
 * Control array to default values before calling klu, you may use:
 *
 *	klu_defaults (Control) ;
 *
 * You can then modify the control parameters before calling klu.
 *
 * Seven arrays are allocated on output of klu, and returned as (int **) or
 * (double **) pointers.  The call to klu should look like the following:
 *
 *	int n, Ap [n+1], Ai [anz] ;	Ap, Ai, and Ax allocated and defined
 *	double Ax [anz] ;
 *	int result ;
 *	double Control [KLU_CONTROL] ;
 *
 *	int *Lp, *Li, *Up, *Ui, *P ;	not allocated or defined
 *	double *Lx, *Ux ;
 *
 *	result = klu (n, Ap, Ai, Ax, Control,
 *		&Lp, &Li, &Lx, &Up, &Ui, &Ux, &P) ;
 *
 * At this point, the matrices P, L, and U are allocated and defined if
 * klu returns KLU_OK (zero).  P [k] = i if row i is the kth pivot row.
 * Row indices of column j of L are stored in Li [Lp [j] ... Lp [j+1]-1],
 * and the values are in the same area of Lx.  U is stored similarly.
 * The columns of L and U are not sorted, except the unit diagonal of L
 * always appears first in each column, and the diagonal of U always appears
 * as the last entry in each column.  You must free the 7 arrays Lp, Li, Lx,
 * Up, Ui, Ux, and P after you are done with them.  You can simply use the
 * ANSI free routine (or mxFree if you are using klu within MATLAB) or call
 * klu_free, as in:
 *
 *	klu_free (&Lp, &Li, &Lx, &Up, &Ui, &Ux, &P) ;
 *
 * To solve Ax=b after calling klu, do the following:
 *
 *	in C notation:			    equivalent MATLAB notation:
 *	klu_permute (n, P, X, B) ;	    x = P*b
 *	klu_lsolve (n, Lp, Li, Lx, X) ;	    x = L\x
 *	klu_usolve (n, Up, Ui, Ux, X) ;	    x = U\x
 *
 * Copyright August 2003, Tim Davis.  All rights reserved.  See the README
 * file for details on permitted use.  Note that no code from The MathWorks,
 * Inc, or from SuperLU, or from any other source appears here.  The code is
 * written from scratch, from the algorithmic description in Gilbert & Peierls'
 * and Eisenstat & Liu's journal papers.
 */

/* ========================================================================== */

#include "klu.h"
#include "klu_kernel.h"

int klu	/* returns 0 if OK, negative if error */
(
    /* inputs, not modified */
    int n,	    /* A is n-by-n. n must be > 0. */
    int Ap [ ],	    /* size n+1, column pointers for A */
    int Ai [ ],	    /* size nz = Ap [n], row indices for A */
    double Ax [ ],  /* size nz, values of A */
    double User_Control [ ],	    /* Control parameters (optional) */

    /* outputs, allocated on output (or NULL if error occurs) */
    int **p_Lp,	    /* Column pointers for L, of size n+1 */
    int **p_Li,	    /* row indices for L */
    double **p_Lx,  /* values of L */
    int **p_Up,	    /* Column pointers for U, of size n+1 */
    int **p_Ui,	    /* row indices for U */
    double **p_Ux,  /* values of U */
    int **p_P	    /* row permutation */
)
{
    int *Pinv, lsize, usize, result, *P, k, anz, *Lp, *Li, *Up, *Ui, prune,
	*Lpend, *Lpruned, *Stack, recursive ;
    double tol, s, *Lx, *Ux, Control [KLU_CONTROL], *X ;
    WorkStackType *WorkStack ;

    /* ---------------------------------------------------------------------- */
    /* get control parameters, or use defaults */
    /* ---------------------------------------------------------------------- */

    n = MAX (1, n) ;
    anz = Ap [n] ;

    if (User_Control == (double *) NULL)
    {
	klu_defaults (Control) ;
    }
    else
    {
	for (k = 0 ; k < KLU_CONTROL ; k++)
	{
	    Control [k] = User_Control [k] ;
	}
    }

    s = Control [KLU_LSIZE] ;
    if (s < 0)
    {
	lsize = (int) ((-s) * anz) ;
    }
    else
    {
	lsize = (int) s ;
    }
    s = Control [KLU_USIZE] ;
    if (s < 0)
    {
	usize = (int) ((-s) * anz) ;
    }
    else
    {
	usize = (int) s ;
    }
    tol = Control [KLU_TOL] ;
    prune = Control [KLU_PRUNE] != 0 ;
    recursive = Control [KLU_RECURSIVE] != 0 ;

    /* make sure control parameters are in legal range */
    lsize = MAX (1, lsize) ;
    usize = MAX (1, usize) ;
    tol   = MIN (tol, 1.0) ;
    tol   = MAX (0.0, tol) ;

    /* return arguments L, U, and P are not yet assigned */
    *p_Lp = (int *) NULL ;
    *p_Li = (int *) NULL ;
    *p_Lx = (double *) NULL ;
    *p_Up = (int *) NULL ;
    *p_Ui = (int *) NULL ;
    *p_Ux = (double *) NULL ;
    *p_P = (int *) NULL ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace and outputs */
    /* ---------------------------------------------------------------------- */

    X  = (double *) ALLOCATE (n * sizeof (double)) ;
    Pinv  = (int *) ALLOCATE (n * sizeof (int)) ;
    Stack = (int *) ALLOCATE (n * sizeof (int)) ;
    if (prune)
    {
	Lpend   = (int *) ALLOCATE (n * sizeof (int)) ;
	Lpruned = (int *) ALLOCATE (n * sizeof (int)) ;
    }
    else
    {
	Lpend = (int *) NULL ;
	Lpruned = (int *) NULL ;
    }
    if (!recursive)
    {
	WorkStack = (WorkStackType *) ALLOCATE (n * sizeof (WorkStackType)) ;
    }
    else
    {
	WorkStack = (WorkStackType *) NULL ;
    }

    /* create sparse matrix for P, L, and U */
    P  = (int *) ALLOCATE (n * sizeof (int)) ;
    Lp = (int *) ALLOCATE ((n+1) * sizeof (int)) ;
    Up = (int *) ALLOCATE ((n+1) * sizeof (int)) ;
    Li = (int *) ALLOCATE (lsize * sizeof (int)) ;
    Ui = (int *) ALLOCATE (usize * sizeof (int)) ;
    Lx = (double *) ALLOCATE (lsize * sizeof (double)) ;
    Ux = (double *) ALLOCATE (usize * sizeof (double)) ;

    if ((X == (double *) NULL) || (Pinv == (int *) NULL) ||
	(Stack == (int *) NULL) || (Lp == (int *) NULL) ||
	(Li == (int *) NULL) || (Lx == (double *) NULL) ||
	(Up == (int *) NULL) || (Ui == (int *) NULL) ||
	(Ux == (double *) NULL) || (P == (int *) NULL) ||
	(prune && ((Lpend == (int *) NULL) || (Lpruned == (int *) NULL)))
	)
    {
	FREE (X, double) ;
	FREE (Pinv, int) ;
	FREE (Stack, int) ;
	FREE (Lpend, int) ;
	FREE (Lpruned, int) ;
	klu_free (&Lp, &Li, &Lx, &Up, &Ui, &Ux, &P) ;
	return (KLU_OUT_OF_MEMORY) ;
    }

    /* ---------------------------------------------------------------------- */
    /* factorize */
    /* ---------------------------------------------------------------------- */

    if (recursive)
    {
#if 1
#include "assert.h"      
      assert( 0 ) ; 
#else
        Ken Stanley commented this out so 
	if (prune)
	{
	    result = klu_prune
		(n, Ap, Ai, Ax, tol, &lsize, &usize,
		Lp, &Li, &Lx, Up, &Ui, &Ux, Pinv, P, X, Stack,
		WorkStack, Lpend, Lpruned) ;
	}
	else
	{
	    result = klu_noprune
		(n, Ap, Ai, Ax, tol, &lsize, &usize,
		Lp, &Li, &Lx, Up, &Ui, &Ux, Pinv, P, X, Stack,
		WorkStack, Lpend, Lpruned) ;
	}
#endif
    }
    else
    {
	if (prune)
	{
	    result = klu_prune_nonrecursive
		(n, Ap, Ai, Ax, tol, &lsize, &usize,
		Lp, &Li, &Lx, Up, &Ui, &Ux, Pinv, P, X, Stack,
		WorkStack, Lpend, Lpruned) ;
	}
	else
	{
	    result = klu_noprune_nonrecursive
		(n, Ap, Ai, Ax, tol, &lsize, &usize,
		Lp, &Li, &Lx, Up, &Ui, &Ux, Pinv, P, X, Stack,
		WorkStack, Lpend, Lpruned) ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* free workspace */
    /* ---------------------------------------------------------------------- */

    FREE (X, double) ;
    FREE (Pinv, int) ;
    FREE (Stack, int) ;
    FREE (Lpend, int) ;
    FREE (Lpruned, int) ;
    FREE (WorkStack, WorkStackType) ;

    /* ---------------------------------------------------------------------- */
    /* return P, L, and U, or return nothing if an error occurred */
    /* ---------------------------------------------------------------------- */

    if (result != KLU_OK)
    {
	klu_free (&Lp, &Li, &Lx, &Up, &Ui, &Ux, &P) ;
    }
    else
    {
	*p_P = P ;
	*p_Lp = Lp ;
	*p_Li = Li ;
	*p_Lx = Lx ;
	*p_Up = Up ;
	*p_Ui = Ui ;
	*p_Ux = Ux ;
    }

    return (result) ;
}


/* ========================================================================== */
/* === klu_free ============================================================= */
/* ========================================================================== */

/* Frees the P, L, and U matrices constructed by klu.  Note that if any or all
 * of the arrays are already free'd, they are not free'd again. */

void klu_free
(
    int **p_Lp,
    int **p_Li,
    double **p_Lx,
    int **p_Up,
    int **p_Ui,
    double **p_Ux,
    int **p_P
)
{
    FREE (*p_Lp, int) ;
    FREE (*p_Li, int) ;
    FREE (*p_Lx, double) ;
    FREE (*p_Up, int) ;
    FREE (*p_Ui, int) ;
    FREE (*p_Ux, double) ;
    FREE (*p_P, int) ;
}


/* ========================================================================== */
/* === klu_lsolve =========================================================== */
/* ========================================================================== */

/* Solve Lx=b.  Assumes L is unit lower triangular and where the unit diagonal
 * entry is stored (and appears first in each column of L).  Overwrites B
 * with the solution X. */

void klu_lsolve
(
    /* inputs, not modified: */
    int n,
    int Lp [ ],
    int Li [ ],
    double Lx [ ],
    /* right-hand-side on input, solution to Lx=b on output */
    double X [ ]
)
{
    int k, p, pend ;
    double xk ;
    for (k = 0 ; k < n ; k++)
    {
	xk = X [k] ;
	if (xk != 0.0)
	{
	    pend = Lp [k+1] ;
	    for (p = Lp [k] + 1 ; p < pend ; p++)
	    {
		X [Li [p]] -= Lx [p] * xk ;
	    }
	}
    }
}


/* ========================================================================== */
/* === klu_usolve =========================================================== */
/* ========================================================================== */

/* Solve Ux=b.  Assumes U is non-unit upper triangular and where the diagonal
 * entry is stored (and appears last in each column of U).  Overwrites B
 * with the solution X. */

void klu_usolve
(
    /* inputs, not modified: */
    int n,
    int Up [ ],
    int Ui [ ],
    double Ux [ ],
    /* right-hand-side on input, solution to Ux=b on output */
    double X [ ]
)
{
    int k, p, pend ;
    double xk ;
    for (k = n-1 ; k >= 0 ; k--)
    {
	pend = Up [k+1] - 1 ;
	xk = X [k] / Ux [pend] ;
	X [k] = xk ;
	if (xk != 0.0)
	{
	    for (p = Up [k] ; p < pend ; p++)
	    {
		X [Ui [p]] -= Ux [p] * xk ;
	    }
	}
    }
}


/* ========================================================================== */
/* === klu_permute ========================================================== */
/* ========================================================================== */

/* Permute a vector with the permutation matrix P, x = P*b. */

void klu_permute
(
    /* inputs, not modified: */
    int n,
    int P [ ],
    double B [ ],
    /* output */
    double X [ ]
)
{
    int k ;
    for (k = 0 ; k < n ; k++)
    {
	X [k] = B [P [k]] ;
    }
}


/* ========================================================================== */
/* === klu_defaults ========================================================= */
/* ========================================================================== */

/* Set the default parameters:
 *
 * Control [KLU_TOL] is the partial pivoting tolerance.  The diagonal entry is
 * used if its absolute value is >= Control [KLU_TOL] times the largest absolute
 * entry in the column.  A value of 1.0 is convential partial pivoting.
 *
 * Control [KLU_LSIZE] controls the initial size of L.  If negative, then L
 * starts out with space enough for (-Control [KLU_LSIZE]) times the number
 * of nonzeros in A.  If positive, then L starts would with a size equal to
 * Control [KLU_LSIZE].
 *
 * Control [KLU_USIZE] is the same, except that it controls the initial size
 * of U.
 *
 * If either L or U are found not to be large enough during factorization, they
 * are expanded in size (via the ANSI realloc, or the MATLAB mxRealloc routine).
 * When klu is done, L and U are shrunk in size, via realloc or mxRealloc, to be
 * just large enough to hold all of the nonzeros in L and U.
 *
 * Control [KLU_PRUNE] determines whether or not to use symmetric pruning.
 * This has no affect on the output factorization.  Symmetric pruning has the
 * potential to reduce the run time for matrices with nonzero pattern that is
 * at least moderately symmetric.  Default is to do pruning.
 *
 * Control [KLU_RECURSIVE] determines whether the recursive or non-recursive
 * version of the depth-first-search is to be used.  The non-recursive version
 * is immune to possible stack overflow.  Default is to use the non-recursive
 * version.
 */

void klu_defaults
(
    double Control [ ]
)
{
    Control [KLU_TOL] = 1.0 ;	    /* partial pivoting tolerance */
    Control [KLU_LSIZE] = -10 ;	    /* L starts out as size 10*nnz(A) */
    Control [KLU_USIZE] = -10 ;	    /* U starts out as size 10*nnz(A) */
    Control [KLU_PRUNE] = TRUE ;    /* use symmetric pruning */
    Control [KLU_RECURSIVE] = FALSE ;	/* use non-recursive DFS */
}
