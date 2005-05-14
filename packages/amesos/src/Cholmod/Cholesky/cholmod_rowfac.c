/* ========================================================================== */
/* === Cholesky/cholmod_rowfac ============================================== */
/* ========================================================================== */

/*
 * CHOLMOD/Cholesky version 0.1, May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Authors: Timothy A. Davis
 * and Sivasankaran Rajamanickam
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* Full or incremental numerical LDL' factorization (simplicial, not supernodal)
 * cholmod_factorize is the "easy" wrapper for this code, but it does not
 * provide access to incremental factorization.
 *
 * cholmod_rowfac computes the full or incremental LDL' factorization of
 * A+beta*I (where A is symmetric) or A*F+beta*I (where A and F are unsymmetric
 * and only the upper triangular part of A*F+beta*I is used).  It computes
 * L and D one row at a time.
 *
 * A is nrow-by-ncol or nrow-by-nrow.  In "packed" form it is a conventional
 * column-oriented sparse matrix.  Row indices of column j are in
 * Ai [Ap [j] ... Ap [j+1]-1] and values in the same locations of Ax.
 * will be faster if A has sorted columns.  In "unpacked" form the column
 * of A ends at Ap [j] + Anz [j] - 1 instead of Ap [j+1] - 1.
 *
 * Row indices in each column of A can be sorted or unsorted, but the routine
 * routine works fastest if A is sorted, or if only triu(A) is provided
 * for the symmetric case.
 *
 * The unit-diagonal nrow-by-nrow output matrix L is returned in "unpacked"
 * column form, with row indices of column j in Li [Lp [j] ...
 * Lp [j] + Lnz [j] - 1] and values in the same location in Lx.  The row
 * indices in each column of L are in sorted order.  The unit diagonal of L
 * is not stored.
 *
 * L can be a simplicial symbolic (CHOLMOD_SYMBOLIC), unpacked LDL', or
 * dynamic LDL'.  A symbolic factor is converted immediately into a numeric
 * factor containing the identity matrix.
 *
 * For a full factorization, kstart = 0 and kend = nrow.  To compute an
 * incremental factorization, select kstart and kend as the range of rows of
 * L you wish to compute.  A correct factorization will be computed only if all
 * descendants of all nodes k = kstart to kend-1 in the etree have been
 * factorized by a prior call to this routine.
 *
 * ---------------
 * Symmetric case:
 * ---------------
 *
 *	The factorization (in MATLAB notation) is:
 *
 *	S = beta*I + A
 *	S = triu (S) + triu (S,1)'
 *	L*D*L' = S
 *
 *	A is a conventional sparse matrix in compressed column form.  Only the
 *	diagonal and upper triangular part of A is accessed; the lower
 *	triangular part is ignored and assumed to be equal to the upper
 *	triangular part.  For an incremental factorization, only columns kstart
 *	to kend-1 of A are accessed.  F is not used.
 *
 * ---------------
 * Unsymmetric case:
 * ---------------
 *
 *	The factorization (in MATLAB notation) is:
 *
 *	S = beta*I + A*F
 *	S = triu (S) + triu (S,1)'
 *	L*D*L' = S
 *
 *	The typical case is F=A'.  Alternatively, if F=A(:,f)', then this
 *	routine factorizes S = beta*I + A(:.f)*A(:,f)'.
 *
 *	All of A and F are accessed, but only the upper triangular part of A*F
 *	is used.  F must be of size A->ncol by A->nrow.  F is used for the
 *	unsymmetric case only.  F can be packed or unpacked and it need not be
 *	sorted.
 *
 *	For a complete factorization of beta*I + A*A',
 *	this routine performs a number of flops exactly equal to:
 *
 *	sum (for each column j of A) of (Anz (j)^2 + Anz (j)), to form S
 *	+
 *	sum (for each column j of L) of (Lnz (j)^2 + 3*Lnz (j)), to factorize S
 *
 *	where Anz (j) is the number of nonzeros in column j of A, and Lnz (j)
 *	is the number of nonzero in column j of L below the diagonal.
 *
 *
 * workspace: Flag (nrow), W (nrow), Iwork (nrow)
 */

#include "cholmod_cholesky.h"
#include "cholmod_internal.h"

/* ========================================================================== */
/* === subtree ============================================================== */
/* ========================================================================== */

/* Compute the nonzero pattern of the sparse triangular solve Lx=b, where L in
 * this case is L(0:k-1,0:k-1), and b is a column of A.  This is done by
 * traversing the kth row-subtree of the elimination tree of L, starting from
 * each nonzero entry in b.  The pattern is returned postordered, and is valid
 * for a subsequent numerical triangular solve of Lx=b.  The elimination tree
 * can be provided in a Parent array, or extracted from the pattern of L itself.
 *
 * The pattern of x = inv(L)*b is returned in Stack [top...].
 * Also scatters b, or a multiple of b, into the work vector W. */

/*
static int subtree
(
    int symmetric, double ftk,
    int top, int k, int p, int pend, int *Ai, double *Ax, int sorted,
    int *Flag, int mark, int *Stack, double *W,
    int *Parent, int *Lp, int *Li, int *Lnz
)
{
*/

#define SUBTREE(symmetric,values,etree) \
    for ( ; p < pend ; p++) \
    { \
	i = Ai [p] ; \
	if (i > k) \
	{ \
	    if (sorted) \
	    { \
		break ; \
	    } \
	    else \
	    { \
		continue ; \
	    } \
	} \
	if (values) \
	{ \
	    /* scatter the column, or a multiple of it, into W */ \
	    if (symmetric) \
	    { \
		W [i] = Ax [p] ; \
	    } \
	    else \
	    { \
		W [i] += Ax [p] * ftk ; \
	    } \
	} \
	/* start at node i and traverse up the subtree, stop at node k */ \
	for (len = 0 ; i < k && i != EMPTY && Flag [i] < mark ; i = parent) \
	{ \
	    /* L(k,i) is nonzero, and seen for the first time */ \
	    Stack [len++] = i ;	/* place i on the stack */ \
	    Flag [i] = mark ;	/* mark i as visited */ \
	    if (etree) \
	    { \
		/* parent of i is given in the Parent array */ \
		parent = Parent [i] ; \
	    } \
	    else \
	    { \
		/* parent of i is first off-diagonal entry in column i of L */ \
		parent = ((Lnz [i] > 1) ? (Li [Lp [i] + 1]) : EMPTY) ; \
	    } \
	} \
	/* move the path down to the bottom of the stack */ \
	while (len > 0) \
	{ \
	    Stack [--top] = Stack [--len] ; \
	} \
    }


/* ========================================================================== */
/* === cholmod_row_subtree ================================================== */
/* ========================================================================== */

/* Compute the nonzero pattern of the solution to the lower triangular system
 * L(0:k-1,0:k-1) * x = A (0:k-1,k) if A is symmetric, or
 * L(0:k-1,0:k-1) * x = A (0:k-1,:) * A (:,k)' if A is unsymmetric.
 * This gives the nonzero pattern of row k of L (excluding the diagonal).
 * The pattern is returned postordered.
 *
 * The result is returned in R, a pre-allocated sparse matrix of size nrow-by-1,
 * with R->nzmax >= nrow.  R is assumed to be packed (Rnz [0] is not updated);
 * the number of entries in R is given by Rp [0].
 *
 * FUTURE WORK:  a very minor change to this routine could allow it to compute
 * the nonzero pattern of x for any system Lx=b.
 *
 * workspace: Flag (nrow)
 */

int cholmod_row_subtree
(
    /* inputs, not modified */
    cholmod_sparse *A,
    cholmod_sparse *F,	    /* F = A', for unsymmetric case only */
    size_t krow,
    void *Parent_p,
    
    /* output, must be allocated on input */
    cholmod_sparse *R,	    /* a 1-by-(A->nrow) matrix, with R->nzmax >= nrow.
			     * Contains the pattern of L (k,0:k-1) on output. */
    cholmod_common *Common
)
{
    double ftk = 0. ;
    double *Ax = NULL, *W = NULL ;
    int *Rp, *Stack, *Flag, *Li = NULL, *Lp = NULL, *Lnz = NULL,
	*Ap, *Ai, *Anz, *Fp, *Fi, *Fnz, *Parent ;
    int p, pend, parent, t, symmetric, nrow, k, pf, pfend, Fpacked, packed,
	sorted, top, len, i, mark ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (R, FALSE) ;
    Parent = Parent_p ;
    RETURN_IF_NULL (Parent, FALSE) ;
    symmetric = A->stype ;
    if (!symmetric)
    {
	RETURN_IF_NULL (F, FALSE) ;
    }
    if (krow >= A->nrow)
    {
	cholmod_error (CHOLMOD_INVALID, "cholmod_row: k invalid", Common) ;
	return (FALSE) ;
    }
    if (R->ncol != 1 || A->nrow != R->nrow || A->nrow > R->nzmax)
    {
	cholmod_error (CHOLMOD_INVALID, "cholmod_row: R invalid", Common) ;
	return (FALSE) ;
    }

    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    nrow = A->nrow ;
    cholmod_allocate_work (nrow, 0, 0, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;
    }
    ASSERT (cholmod_dump_work (TRUE, TRUE, 0, Common)) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    if (symmetric > 0)
    {
	/* symmetric upper case: F is not needed.  It may be NULL */
	Fp = NULL ;
	Fi = NULL ;
	Fnz = NULL ;
	Fpacked = TRUE ;
    }
    else if (!symmetric)
    {
	/* unsymmetric case: F is required. */
	Fp = F->p ;
	Fi = F->i ;
	Fnz = F->nz ;
	Fpacked = F->packed ;
    }
    else
    {
	/* symmetric lower triangular form not supported */
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_row: symmetric lower not supported", Common) ;
	return (FALSE) ;
    }

    Ap = A->p ;
    Ai = A->i ;
    Anz = A->nz ;
    packed = A->packed ;
    sorted = A->sorted ;

    k = krow ;
    Stack = R->i ;

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    Flag = Common->Flag ;	/* size nrow, Flag [i] < mark must hold */
    mark = cholmod_clear_flag (Common) ;

    /* ---------------------------------------------------------------------- */
    /* compute the pattern of L(k,:) */
    /* ---------------------------------------------------------------------- */

    top = nrow ;		/* Stack is empty */
    Flag [k] = mark ;		/* do not include diagonal entry in Stack */

    if (symmetric)
    {
	/* scatter kth col of triu (beta*I+AA'), get pattern L(k,:) */
	p = Ap [k] ;
	pend = (packed) ? (Ap [k+1]) : (p + Anz [k]) ;
	SUBTREE (TRUE, FALSE, TRUE) ;
    }
    else
    {
	/* scatter kth col of triu (beta*I+AA'), get pattern L(k,:) */
	pf = Fp [k] ;
	pfend = (Fpacked) ? (Fp [k+1]) : (pf + Fnz [k]) ;
	for ( ; pf < pfend ; pf++)
	{
	    /* get nonzero entry F (t,k) */
	    t = Fi [pf] ;
	    p = Ap [t] ;
	    pend = (packed) ? (Ap [t+1]) : (p + Anz [t]) ;
	    SUBTREE (FALSE, FALSE, TRUE) ;
	}
    }

    /* shift the stack upwards, to the first part of R */
    len = nrow - top ;
    for (i = 0 ; i < len ; i++)
    {
	Stack [i] = Stack [top + i] ;
    }

    Rp = R->p ;
    Rp [0] = 0 ;
    Rp [1] = len ;
    R->sorted = FALSE ;

    cholmod_clear_flag (Common) ;
    ASSERT (cholmod_dump_work (TRUE, TRUE, 0, Common)) ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_rowfac ======================================================= */
/* ========================================================================== */

int cholmod_rowfac
(
    /* inputs, not modified */
    cholmod_sparse *A,	/* nrow-by-ncol, stored in column form */
    cholmod_sparse *F,	/* used for LDL'=AA' case only, stored in row form */
    cholmod_scalar beta,/* factorize beta*I+A or beta*I+AA' */
    size_t kstart,	/* first row to factorize */
    size_t kend,	/* last row to factorize is kend-1 */

    /* input, modified on output: */
    cholmod_factor *L,

    cholmod_common *Common
)
{
    double dk, di, yi, l_ki, ftk = 0. ;
    double *Ax, *Lx, *W, *Fx ;
    int *Ap, *Anz, *Ai, *Lp, *Lnz, *Li, *Lnext, *Flag, *Stack, *Fp, *Fi, *Fnz,
	*Iwork, *Parent ;
    int i, p, k, n, t, pf, pfend, top, s, mark, pend, nrow, lnz,
	use_dmin, packed, symmetric, Fpacked, sorted, nzmax, len, parent ;

    DEBUG (int orig) ;
    DEBUG (int new = 0) ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    symmetric = A->stype ;
    if (!symmetric)
    {
	RETURN_IF_NULL (F, FALSE) ;
    }
    nrow = A->nrow ;
    n = L->n ;
    if (kend > L->n)
    {
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_rowfac: kend invalid", Common) ;
	return (FALSE) ;
    }
    if (nrow != n)
    {
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_rowfac: dimensions of A and L do not match", Common) ;
	return (FALSE) ;
    }

    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    cholmod_allocate_work (nrow, nrow, nrow, sizeof (double), Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;
    }
    DEBUG (orig = Common->malloc_count) ;
    ASSERT (cholmod_dump_work (TRUE, TRUE, nrow, Common)) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    PRINT1 (("\nin cholmod_rowfac, kstart %d kend %d stype %d\n",
		kstart, kend, A->stype)) ;
    DEBUG (cholmod_dump_factor (L, "Initial L\n", Common)) ;

    if (!(L->ftype == CHOLMOD_LDL_UNPACKED
	 || L->ftype == CHOLMOD_LDL_DYNAMIC
	 || L->ftype == CHOLMOD_SYMBOLIC))
    {
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_rowfac: can only do LDL' factorization", Common) ;
	return (FALSE) ;
    }

    if (symmetric > 0)
    {
	/* symmetric upper case: F is not needed.  It may be NULL */
	Fp = NULL ;
	Fi = NULL ;
	Fx = NULL ;
	Fnz = NULL ;
	Fpacked = TRUE ;
    }
    else if (!symmetric)
    {
	/* unsymmetric case: F is required. */
	Fp = F->p ;
	Fi = F->i ;
	Fx = F->x ;
	Fnz = F->nz ;
	Fpacked = F->packed ;
    }
    else
    {
	/* symmetric lower triangular form not supported */
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_rowfac: symmetric lower not supported", Common) ;
	return (FALSE) ;
    }

    Ap = A->p ;		/* size A->ncol+1, column pointers of A */
    Ai = A->i ;		/* size nz = Ap [A->ncol], row indices of A */
    Ax = A->x ;		/* size nz, numeric values of A */
    Anz = A->nz ;
    packed = A->packed ;
    sorted = A->sorted ;

    use_dmin = (Common->dmin > 0) ;	/* fl.pt. compare; false if NaN */

    /* get the current factors L and D, and allocate space if needed */
    if (L->ftype == CHOLMOD_SYMBOLIC)
    {
	/* ------------------------------------------------------------------ */
	/* L is symbolic only; allocate and initialize L and D */
	/* ------------------------------------------------------------------ */

	/* workspace: none */
	if (!cholmod_change_ftype (L, CHOLMOD_LDL_DYNAMIC, Common))
	{
	    /* out of memory */
	    ASSERT (Common->malloc_count == orig) ;
	    return (FALSE) ;
	}
	DEBUG (new = cholmod_dump_mdelta
		(CHOLMOD_SYMBOLIC, CHOLMOD_LDL_DYNAMIC)) ;
    }
    ASSERT (Common->malloc_count == orig + new) ;
    ASSERT (L->ftype == CHOLMOD_LDL_UNPACKED
	 || L->ftype == CHOLMOD_LDL_DYNAMIC) ;
    ASSERT (L->xtype != CHOLMOD_PATTERN) ;

    /* inputs, can be modified on output: */
    Lp = L->p ;		/* size nrow+1 */
    ASSERT (Lp != NULL) ;

    /* outputs, contents defined on input for incremental case only: */
    Lnz = L->nz ;	/* size nrow */
    Lnext = L->next ;
    Li = L->i ;		/* size L->nzmax, can change in size */
    Lx = L->x ;		/* size L->nzmax, can change in size */
    nzmax = L->nzmax ;
    ASSERT (Lnz != NULL && Li != NULL && Lx != NULL) ;

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    Iwork = Common->Iwork ;
    Stack = Iwork ;		/* size nrow (i/i/l) */
    Flag = Common->Flag ;	/* size nrow, Flag [i] < mark must hold */
    W = Common->Xwork ;		/* size nrow, W [i] == 0 must hold */
    mark = Common->mark ;

    /* ---------------------------------------------------------------------- */
    /* compute LDL' factorization by rows */
    /* ---------------------------------------------------------------------- */

    for (k = kstart ; k < ((int) kend) ; k++)
    {
	PRINT1 (("\n===============K %d\n", k)) ;

	/* ------------------------------------------------------------------ */
	/* compute pattern of kth row of L and scatter kth input column */
	/* ------------------------------------------------------------------ */

	/* column k of L is currently empty */
	ASSERT (Lnz [k] == 1) ;
	ASSERT (cholmod_dump_work (TRUE, TRUE, nrow, Common)) ;

	top = nrow ;		/* Stack is empty */
	Flag [k] = mark ;	/* do not include diagonal entry in Stack */

	if (symmetric)
	{
	    /* scatter kth col of triu (beta*I+AA'), get pattern L(k,:) */
	    p = Ap [k] ;
	    pend = (packed) ? (Ap [k+1]) : (p + Anz [k]) ;
	    SUBTREE (TRUE, TRUE, FALSE) ;
	}
	else
	{
	    /* scatter kth col of triu (beta*I+AA'), get pattern L(k,:) */
	    pf = Fp [k] ;
	    pfend = (Fpacked) ? (Fp [k+1]) : (pf + Fnz [k]) ;
	    for ( ; pf < pfend ; pf++)
	    {
		/* get nonzero entry F (t,k) */
		t = Fi [pf] ;
		ftk = Fx [pf] ;
		p = Ap [t] ;
		pend = (packed) ? (Ap [t+1]) : (p + Anz [t]) ;
		SUBTREE (FALSE, TRUE, FALSE) ;
	    }
	}

	/* nonzero pattern of kth row of L is now in Stack [top..nrow-1].
	 * Flag [Stack [top..nrow-1]] is equal to mark, but no longer needed */

	mark = cholmod_clear_flag (Common) ;

	/* ------------------------------------------------------------------ */
	/* compute kth row of L and store in column form */
	/* ------------------------------------------------------------------ */

	/* Solve L (0:k-1, 0:k-1) * y (0:k-1) = b (0:k-1) where
	 * b (0:k) = A (0:k,k) or A(0:k,:) * F(:,k) is in W and Stack.
	 * scale kth row of L, compute diagonal D [k], and store in col form.
	 * L (k, 0:k-1) = y (0:k-1) ./ D (0:k-1)
	 * D (k) = b (k) - L (k, 0:k-1) * y (0:k-1)
	 */

	dk = W [k] + beta.x ;
	W [k] = 0.0 ;
	for (s = top ; s < nrow ; s++)
	{
	    /* get i for each nonzero entry L(k,i) */
	    i = Stack [s] ;

	    /* forward solve using L ((i+1):(k-1),i) */
	    yi = W [i] ;
	    W [i] = 0.0 ;
	    lnz = Lnz [i] ;
	    p = Lp [i] ;
	    ASSERT (lnz > 0 && Li [p] == i) ;
	    pend = p + lnz ;
	    di = Lx [p++] ;
	    for ( ; p < pend ; p++)
	    {
		W [Li [p]] -= Lx [p] * yi ;
	    }
	    /* Scale L (k,0:k-1) and compute dot product for D(k) */
	    l_ki = yi / di ;
	    dk -= l_ki * yi ;

	    /* determine if column i of L can hold the new L(k,i) entry */
	    if (p >= COLEND (i))
	    {
		/* column i needs to grow */
		PRINT1 (("Factor Colrealloc %d, old Lnz %d\n", i, Lnz [i])) ;
		if (!cholmod_reallocate_column (L, i, lnz + 1, Common))
		{
		    /* out of memory, L is now simplicial symbolic */
		    for (i = 0 ; i < nrow ; i++)
		    {
			W [i] = 0 ;
		    }
		    ASSERT (cholmod_dump_work (TRUE, TRUE, nrow, Common)) ;
		    return (FALSE) ;
		}
		Li = L->i ;		/* L->i and L->x may have moved */
		Lx = L->x ;
		p = Lp [i] + lnz ;	/* contents of L->p changed */
		Lnext = L->next ;	/* L->ftype may have changed */
		ASSERT (p < COLEND (i)) ;
	    }

	    /* store L (k,i) in the column form matrix of L */
	    Li [p] = k ;
	    Lx [p] = l_ki ;
	    Lnz [i]++ ;
	}

	/* ------------------------------------------------------------------ */
	/* ensure ABS (dk) >= dmin if dmin is given, and store it in L */
	/* ------------------------------------------------------------------ */

	p = Lp [k] ;
	Lx [p] = (use_dmin ? cholmod_dmin (dk, Common) : dk) ;
	Li [p] = k ;
    }

    DEBUG (cholmod_dump_factor (L, "final cholmod_rowfac", Common)) ;
    ASSERT (cholmod_dump_work (TRUE, TRUE, nrow, Common)) ;
    return (TRUE) ;
}
