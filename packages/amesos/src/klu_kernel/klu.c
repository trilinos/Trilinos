/* ========================================================================== */
/* === klu_user ============================================================= */
/* ========================================================================== */

/* KLU: factorizes P*A into L*U, using the Gilbert-Peierls algorithm, with
 * optional symmetric pruning by Eisenstat and Liu.  The code is by Tim Davis.
 * This algorithm is what appears as the default sparse LU routine in MATLAB
 * version 6.0, and still appears in MATLAB 6.5 as [L,U,P] = lu (A).  Note
 * that no column ordering is provided (see COLAMD or AMD for suitable
 * orderings).  SuperLU is based on this algorithm, except that it adds the
 * use of dense matrix operations on "supernodes" (adjacent columns with
 * identical).  This code doesn't use supernodes, thus its name ("Kent" LU,
 * as in "Clark Kent", in contrast with Super-LU...).  This algorithm is slower
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
 *
 * See klu_btf for the BTF+AMD ordering, and scaling.
 *
 * The 7 output arrays
 * should be allocated in one call to malloc, and then grown in size via
 * realloc.  This typically allows realloc to perform without copying.
 * UMFPACK includes all of the above features (except it doesn't do BTF).
 *
 * No fill-reducing ordering is provided!  The ordering quality of klu is your
 * responsibility.  You must pre-permute A to reduce fill-in, or provide a
 * fill-reducing input permutation Q.  Two typical scenarios (in MATLAB
 * notation) are:
 *
 *	Q = colamd (A) ;
 *	[L,U,P] = klu (A, Q) ;
 *
 * or
 *
 *	Q = symamd (A) ;	    % or Q = amd (A)
 *	[L,U,P] = klu (A, Q) ;
 *
 * Both return L, U, and P such that L*U = A(P,Q).  In the latter case, a
 * pivot tolerance of 1e-3 would be appropriate, to preserve symmetry.  In
 * the former case, use a pivot tolerance of 1.0.  See the MATLAB interface
 * klu.m for more details.
 *
 * The input matrix A must be in compressed-column form, with either sorted 
 * or unsorted row indices.  Row indices for column j of A is in
 * Ai [Ap [j] ... Ap [j+1]-1] and the same range of indices in Ax holds the
 * numerical values.  No duplicate entries are allowed.
 *
 * Control is ignored if it is (double *) NULL.  To set the
 * Control array to default values before calling klu, you may use:
 *
 *	klu_defaults (Control) ;
 *
 * You can then modify the control parameters before calling klu.
 *
 * Four arrays are allocated on output of klu, and returned as (int **) or
 * (double **) pointers.  The call to klu should look like the following:
 *
 *	int n, Ap [n+1], Ai [anz] ;	Ap, Ai, and Ax allocated and defined
 *	double Ax [anz] ;
 *	int result ;
 *	int Q [n] ;			Q is optional, it may be (int *) NULL
 *	double Control [KLU_CONTROL] ;
 *
 *	int Lp [n+1], Up [n+1], P [n] ;
 *	int *Li, *Ui, ;			not allocated or defined
 *	double *Lx, *Ux ;
 *
 *	TODO: add discussion of Xwork and Iwork, BTF arguments
 *
 *	result = klu_factor (n, Ap, Ai, Ax, Q, Control,
 *		Lp, &Li, &Lx, Up, &Ui, &Ux, P,
 *		...) ;
 *
 * At this point, the arrays Li, Lx, Ui, and Ux are allocated and defined if
 * klu returns KLU_OK (zero).  P [k] = i if row i is the kth pivot row.
 * Row indices of column j of L are stored in Li [Lp [j] ... Lp [j+1]-1],
 * and the values are in the same area of Lx.  U is stored similarly.
 * The columns of L and U are not sorted, except the unit diagonal of L
 * always appears first in each column, and the diagonal of U always appears
 * as the last entry in each column.  You must free the 4 arrays Li, Lx,
 * Ui, and Ux after you are done with them.  You can simply use the
 * ANSI free routine (or mxFree if you are using klu within MATLAB) or call
 * klu_free, as in:
 *
 *	klu_free (&Li, &Lx, &Ui, &Ux) ;
 *
 * To solve Ax=b after calling klu, do the following:
 *
 *	in C notation:			    equivalent MATLAB notation:
 *	klu_permute (n, P, B, n, nrhs, X) ;		x = P*b
 *	klu_lsolve (n, Lp, Li, Lx, n, nrhs, X) ;	x = L\x
 *	klu_usolve (n, Up, Ui, Ux, n, nrhs, X) ;	x = U\x
 *	klu_permute (n, Q, B, n, nrhs, X) ;		x = Q*x
 *
 * Copyright August 2004, Tim Davis.  All rights reserved.  See the README
 * file for details on permitted use.  Note that no code from The MathWorks,
 * Inc, or from SuperLU, or from any other source appears here.  The code is
 * written from scratch, from the algorithmic description in Gilbert & Peierls'
 * and Eisenstat & Liu's journal papers.
 *
 * Feb, 2004: DFS code modified.
 *
 * Aug 20, 2003:  the input permutation Q added.  In MATLAB
 * notation, klu now computes the factorization L*U = A (P,Q), where P is
 * determined by partial pivoting, and Q is the input ordering.  If the
 * pivot tolerance is less than 1, the "diagonal" entry that klu attempts to
 * choose is the diagonal of A (Q,Q).  In other words, the input permutation
 * is applied symmetrically to the input matrix.  The output permutation P
 * includes both the partial pivoting ordering and the input permutation.
 * If Q is (int *) NULL, then it is assumed to be the identity permutation.
 * Q is not modified.
 *
 * TODO: Cite Gilbert, Eisenstat here.  Write README file.
 * TODO: do not allocate Lp, Up, or P.  Require them as input.
 * TODO: require Xwork and Iwork to be allocated on input.
 */

/* ========================================================================== */

#include "klu.h"
#include "klu_kernel.h"

int klu_factor	/* returns 0 if OK, negative if error */
(
    /* inputs, not modified */
    int n,	    /* A is n-by-n. n must be > 0. */
    int Ap [ ],	    /* size n+1, column pointers for A */
    int Ai [ ],	    /* size nz = Ap [n], row indices for A */
    double Ax [ ],  /* size nz, values of A */
    int Q [ ],	    /* size n, optional column permutation */
    double Control [ ],	    /* Control parameters (optional) */

    /* outputs, not defined on input */
    int Lp [ ],	    /* Column pointers for L, of size n+1 */
    int **p_Li,	    /* row indices for L */
    double **p_Lx,  /* values of L */
    int Up [ ],	    /* Column pointers for U, of size n+1 */
    int **p_Ui,	    /* row indices for U */
    double **p_Ux,  /* values of U */
    int P [ ],	    /* row permutation, size n */
    int *p_noffdiag,	/* # of off-diagonal pivots chosen */
    double *p_umin,
    double *p_umax,

    /* workspace, undefined on input */
    double *X,	    /* size n double's, zero on output */
    int *Work,	    /* size 5n int's */

    /* ---- the following are only used in the BTF case --- */

    /* inputs, not modified on output */
    int k1,	    /* the block of A is from k1 to k2-1 */
    int PSinv [ ],  /* inverse of P from symbolic factorization */
    double Rs [ ],  /* scale factors for A */

    /* inputs, modified on output */
    int Offp [ ],   /* off-diagonal matrix (modified by this routine) */
    int Offi [ ],
    double Offx [ ]
)
{
    double tol, s, *Lx, *Ux, umin, umax, growth, maxlnz ;
    int *Pinv, lsize, usize, result, anz, *Li, *Ui, *Lpend, *Stack, *Flag,
	noffdiag, *Ap_pos, get_work, get_X, *W, no_btf ;

    /* ---------------------------------------------------------------------- */
    /* get control parameters, or use defaults */
    /* ---------------------------------------------------------------------- */

    no_btf = (Offp == (int *) NULL) ;
    n = MAX (1, n) ;

    if (no_btf)
    {
	anz = Ap [n] ;
    }
    else
    {
	anz = Ap [n+k1] - Ap [k1] ;
    }

    s = GET_CONTROL (KLU_LSIZE, -10.0) ;
    if (s <= 0)
    {
	s = -s ;
	s = MAX (s, 1.0) ;
	lsize = s * anz + n ;
    }
    else
    {
	lsize = s ;
    }

    s = GET_CONTROL (KLU_USIZE, -10.0) ;
    if (s <= 0)
    {
	s = -s ;
	s = MAX (s, 1.0) ;
	usize = s * anz + n ;
    }
    else
    {
	usize = s ;
    }

    tol = GET_CONTROL (KLU_TOL, 1.0) ;
    growth = GET_CONTROL (KLU_GROWTH, 1.2) ;

    /* make sure control parameters are in legal range */
    lsize  = MAX (n+1, lsize) ;
    usize  = MAX (n+1, usize) ;

    maxlnz = (((double) n) * ((double) n) + ((double) n)) / 2. ;
    maxlnz = MIN (maxlnz, ((double) INT_MAX)) ;

    lsize  = MIN (maxlnz, lsize) ;
    usize  = MIN (maxlnz, usize) ;

    tol    = MIN (tol, 1.0) ;
    tol    = MAX (0.0, tol) ;
    growth = MAX (1.0, growth) ;

    PRINTF (("Welcome to klu: n %d anz %d k1 %d lsize %d usize %d maxlnz %g\n",
	n, anz, k1, lsize, usize, maxlnz)) ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace and outputs */
    /* ---------------------------------------------------------------------- */

    /* return arguments L, U, and P are not yet assigned */
    *p_Li = (int *) NULL ;
    *p_Lx = (double *) NULL ;
    *p_Ui = (int *) NULL ;
    *p_Ux = (double *) NULL ;

    get_work = (Work == (int *) NULL) ;
    if (get_work)
    {
	/* allocate workspace of size 5*n*int */
	Work = (int *) ALLOCATE (5*n * sizeof (int)) ;
	if (Work == (int *) NULL)
	{
	    return (KLU_OUT_OF_MEMORY) ;
	}
    }

    W = Work ;
    Pinv = (int *) W ;	    W += n ;
    Stack = (int *) W ;	    W += n ;
    Flag = (int *) W ;	    W += n ;
    Lpend = (int *) W ;	    W += n ;
    Ap_pos = (int *) W ;    W += n ;

    get_X = (X == (double *) NULL) ;
    if (get_X)
    {
	/* allocate workspace of size n*double */
	X = (double *) ALLOCATE (n * sizeof (double)) ;
	if (X == (double *) NULL)
	{
	    return (KLU_OUT_OF_MEMORY) ;
	}
    }

    /* create sparse matrix for P, L, and U */
    Li = (int *) ALLOCATE (lsize * sizeof (int)) ;
    Lx = (double *) ALLOCATE (lsize * sizeof (double)) ;
    Ui = (int *) ALLOCATE (usize * sizeof (int)) ;
    Ux = (double *) ALLOCATE (usize * sizeof (double)) ;

    if ((Li == (int *) NULL) || (Lx == (double *) NULL) ||
        (Ui == (int *) NULL) || (Ux == (double *) NULL))
    {
	if (get_work) FREE (Work, int) ;
	if (get_X)    FREE (X, double) ;
	klu_free (&Li, &Lx, &Ui, &Ux) ;
	return (KLU_OUT_OF_MEMORY) ;
    }

    /* ---------------------------------------------------------------------- */
    /* factorize */
    /* ---------------------------------------------------------------------- */

    /* with pruning, and non-recursive depth-first-search */
    result = klu_kernel (n, Ap, Ai, Ax, Q, tol, growth, lsize, usize,
	    Lp, &Li, &Lx, Up, &Ui, &Ux, Pinv, P, &noffdiag, &umin, &umax,
	    X, Stack, Flag, Ap_pos, Lpend,
	    /* BTF case: */
	    no_btf, k1, PSinv, Rs, Offp, Offi, Offx) ;

    /* ---------------------------------------------------------------------- */
    /* free workspace */
    /* ---------------------------------------------------------------------- */

    if (get_work) FREE (Work, int) ;
    if (get_X)    FREE (X, double) ;

    /* ---------------------------------------------------------------------- */
    /* return P, L, and U, or return nothing if an error occurred */
    /* ---------------------------------------------------------------------- */

    if (result != KLU_OK)
    {
	klu_free (&Li, &Lx, &Ui, &Ux) ;
    }
    else
    {
	*p_Li = Li ;
	*p_Lx = Lx ;
	*p_Ui = Ui ;
	*p_Ux = Ux ;
    }
    *p_noffdiag = noffdiag ;
    *p_umin = umin ;
    *p_umax = umax ;
    PRINTF ((" in klu noffdiag %d\n", noffdiag)) ;

    return (result) ;
}


/* ========================================================================== */
/* === klu_free ============================================================= */
/* ========================================================================== */

/* Frees the L and U matrices constructed by klu.  Note that if any or all
 * of the arrays are already free'd, they are not free'd again. */

void klu_free
(
    int **p_Li,
    double **p_Lx,
    int **p_Ui,
    double **p_Ux
)
{
    FREE (*p_Li, int) ;
    FREE (*p_Lx, double) ;
    FREE (*p_Ui, int) ;
    FREE (*p_Ux, double) ;
}


/* ========================================================================== */
/* === klu_lsolve =========================================================== */
/* ========================================================================== */

/* Solve Lx=b.  Assumes L is unit lower triangular and where the unit diagonal
 * entry is stored (and appears first in each column of L).  Overwrites B
 * with the solution X.  B is n-by-nrhs and is stored in column form with
 * leading dimension d. */

void klu_lsolve
(
    /* inputs, not modified: */
    int n,
    int Lp [ ],
    int Li [ ],
    double Lx [ ],
    int d,
    int nrhs,
    /* right-hand-side on input, solution to Lx=b on output */
    double X [ ]
)
{
    double y [4], lik ;
    double *Y ;
    int k, p, pend, nblocks, block, i, d2, d3, d4 ;

    if (nrhs == 1)
    {
	for (k = 0 ; k < n ; k++)
	{
	    y [0] = X [k] ;
	    if (y [0] != 0.0)
	    {
		pend = Lp [k+1] ;
		for (p = Lp [k] + 1 ; p < pend ; p++)
		{
		    X [Li [p]] -= Lx [p] * y [0] ;
		}
	    }
	}
    }
    else
    {

	/* determine how many blocks of 4 to do */
	nblocks = nrhs / 4 ;

	/* do all the blocks of size 4 */
	Y = X ;
	d2 = d*2 ;
	d3 = d*3 ;
	d4 = d*4 ;
	for (block = 0 ; block < nblocks ; block++)
	{
	    for (k = 0 ; k < n ; k++)
	    {
		y [0] = Y [k     ] ;
		y [1] = Y [k + d ] ;
		y [2] = Y [k + d2] ;
		y [3] = Y [k + d3] ;
		pend = Lp [k+1] ;
		for (p = Lp [k] + 1 ; p < pend ; p++)
		{
		    i = Li [p] ;
		    lik = Lx [p] ;
		    Y [i     ] -= lik * y [0] ;
		    Y [i + d ] -= lik * y [1] ;
		    Y [i + d2] -= lik * y [2] ;
		    Y [i + d3] -= lik * y [3] ;
		}
	    }
	    Y += d4 ;
	}

	/* clean-up for last block of size 1,2,3 */
	switch (nrhs % 4)
	{
	case 0:
	    break ;

	case 1:

	    for (k = 0 ; k < n ; k++)
	    {
		y [0] = Y [k] ;
		pend = Lp [k+1] ;
		for (p = Lp [k] + 1 ; p < pend ; p++)
		{
		    Y [Li [p]] -= Lx [p] * y [0] ;
		}
	    }

	    break ;

	case 2:

	    for (k = 0 ; k < n ; k++)
	    {
		y [0] = Y [k    ] ;
		y [1] = Y [k + d] ;
		pend = Lp [k+1] ;
		for (p = Lp [k] + 1 ; p < pend ; p++)
		{
		    i = Li [p] ;
		    lik = Lx [p] ;
		    Y [i    ] -= lik * y [0] ;
		    Y [i + d] -= lik * y [1] ;
		}
	    }

	    break ;

	case 3:

	    for (k = 0 ; k < n ; k++)
	    {
		y [0] = Y [k     ] ;
		y [1] = Y [k + d ] ;
		y [2] = Y [k + d2] ;
		pend = Lp [k+1] ;
		for (p = Lp [k] + 1 ; p < pend ; p++)
		{
		    i = Li [p] ;
		    lik = Lx [p] ;
		    Y [i     ] -= lik * y [0] ;
		    Y [i + d ] -= lik * y [1] ;
		    Y [i + d2] -= lik * y [2] ;
		}
	    }

	    break ;
	}
    }
}


/* ========================================================================== */
/* === klu_usolve =========================================================== */
/* ========================================================================== */

/* Solve Ux=b.  Assumes U is non-unit upper triangular and where the diagonal
 * entry is stored (and appears last in each column of U).  Overwrites B
 * with the solution X.  B is n-by-nrhs and is stored in column form. */

void klu_usolve
(
    /* inputs, not modified: */
    int n,
    int Up [ ],
    int Ui [ ],
    double Ux [ ],
    int d,
    int nrhs,
    /* right-hand-side on input, solution to Ux=b on output */
    double X [ ]
)
{
    double y [4], uik, ukk ;
    double *Y ;
    int k, p, pend, nblocks, block, i, d2, d3, d4 ;

    if (nrhs == 1)
    {

	for (k = n-1 ; k >= 0 ; k--)
	{
	    pend = Up [k+1] - 1 ;
	    y [0] = X [k] / Ux [pend] ;
	    X [k] = y [0] ;
	    if (y [0] != 0.0)
	    {
		for (p = Up [k] ; p < pend ; p++)
		{
		    X [Ui [p]] -= Ux [p] * y [0] ;
		}
	    }
	}

    }
    else
    {

	/* determine how many blocks of 4 to do */
	nblocks = nrhs / 4 ;

	/* do all the blocks of size 4 */
	Y = X ;
	d2 = d*2 ;
	d3 = d*3 ;
	d4 = d*4 ;
	for (block = 0 ; block < nblocks ; block++)
	{
	    for (k = n-1 ; k >= 0 ; k--)
	    {
		pend = Up [k+1] - 1 ;
		ukk = Ux [pend] ;
		y [0] = Y [k     ] / ukk ;
		y [1] = Y [k + d ] / ukk ;
		y [2] = Y [k + d2] / ukk ;
		y [3] = Y [k + d3] / ukk ;
		Y [k     ] = y [0] ;
		Y [k + d ] = y [1] ;
		Y [k + d2] = y [2] ;
		Y [k + d3] = y [3] ;
		for (p = Up [k] ; p < pend ; p++)
		{
		    i = Ui [p] ;
		    uik = Ux [p] ;
		    Y [i     ] -= uik * y [0] ;
		    Y [i + d ] -= uik * y [1] ;
		    Y [i + d2] -= uik * y [2] ;
		    Y [i + d3] -= uik * y [3] ;
		}
	    }
	    Y += d4 ;
	}

	/* clean-up for last block of size 1,2,3 */
	switch (nrhs % 4)
	{
	case 0:
	    break ;

	case 1:

	    for (k = n-1 ; k >= 0 ; k--)
	    {
		pend = Up [k+1] - 1 ;
		y [0] = Y [k] / Ux [pend] ;
		Y [k] = y [0] ;
		for (p = Up [k] ; p < pend ; p++)
		{
		    Y [Ui [p]] -= Ux [p] * y [0] ;
		}
	    }

	    break ;

	case 2:

	    for (k = n-1 ; k >= 0 ; k--)
	    {
		pend = Up [k+1] - 1 ;
		ukk = Ux [pend] ;
		y [0] = Y [k    ] / ukk ;
		y [1] = Y [k + d] / ukk ;
		Y [k    ] = y [0] ;
		Y [k + d] = y [1] ;
		for (p = Up [k] ; p < pend ; p++)
		{
		    i = Ui [p] ;
		    uik = Ux [p] ;
		    Y [i    ] -= uik * y [0] ;
		    Y [i + d] -= uik * y [1] ;
		}
	    }

	    break ;

	case 3:

	    for (k = n-1 ; k >= 0 ; k--)
	    {
		pend = Up [k+1] - 1 ;
		ukk = Ux [pend] ;
		y [0] = Y [k     ] / ukk ;
		y [1] = Y [k + d ] / ukk ;
		y [2] = Y [k + d2] / ukk ;
		Y [k     ] = y [0] ;
		Y [k + d ] = y [1] ;
		Y [k + d2] = y [2] ;
		for (p = Up [k] ; p < pend ; p++)
		{
		    i = Ui [p] ;
		    uik = Ux [p] ;
		    Y [i     ] -= uik * y [0] ;
		    Y [i + d ] -= uik * y [1] ;
		    Y [i + d2] -= uik * y [2] ;
		}
	    }

	    break ;
	}
    }
}


/* ========================================================================== */
/* === klu_permute ========================================================== */
/* ========================================================================== */

/* Permute a dense matrix with the permutation matrix P, X = P*B. */

void klu_permute
(
    /* inputs, not modified: */
    int n,
    int P [ ],
    int d,
    int nrhs,
    double B [ ],
    /* output */
    double X [ ]
)
{
    double *Y, *Z ;
    int k, s ;
    Y = X ;
    Z = B ;
    for (s = 0 ; s < nrhs ; s++)
    {
	for (k = 0 ; k < n ; k++)
	{
	    Y [k] = Z [P [k]] ;
	}
	Y += d ;
	Z += d ;
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
 */

void klu_defaults
(
    double Control [ ]
)
{
    int i ;
    if (Control == (double *) NULL)
    {
	/* nothing to do - return silently (this is not an error) */
	return ;
    }
    for (i = 0 ; i < KLU_CONTROL ; i++)
    {
	Control [i] = 0 ;
    }
    Control [KLU_TOL] = 1.0 ;	    /* partial pivoting tolerance */
    Control [KLU_LSIZE] = -10 ;	    /* L starts out as size 10*nnz(A) */
    Control [KLU_USIZE] = -10 ;	    /* U starts out as size 10*nnz(A) */
    Control [KLU_GROWTH] = 1.2 ;    /* memory growth factor */
}
