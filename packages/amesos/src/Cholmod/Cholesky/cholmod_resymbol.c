/* ========================================================================== */
/* === Cholesky/cholmod_resymbol ============================================ */
/* ========================================================================== */

/*
 * CHOLMOD/Cholesky version 0.1, May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Authors: Timothy A. Davis
 * and Sivasankaran Rajamanickam
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* Recompute the symbolic pattern of L.  Entries not in the symbolic pattern
 * are dropped.  L->Perm is used to permute the input matrix A.
 *
 * This routine is used after a supernodal factorization is converted into
 * a simplicial one, to remove zero entries that were added due to relaxed
 * supernode amalgamation.  It can also be used after a series of downdates
 * to remove entries that would no longer be present if the matrix were
 * factorized from scratch.  A downdate (cholmod_updown) does not remove any
 * entries from L.
 *
 * workspace: Flag (nrow), Head (nrow+1),
 *	if symmetric:   Iwork (2*nrow)
 *	if unsymmetric: Iwork (2*nrow+ncol).
 *	Allocates up to 2 copies of its input matrix A (pattern only).
 */

#include "cholmod_cholesky.h"
#include "cholmod_internal.h"

int cholmod_resymbol
(
    cholmod_sparse *A,
    void *fset,
    size_t fsize,
    cholmod_factor *L,
    int pack,
    cholmod_common *Common
)
{
    cholmod_sparse *H, *F, *G ;
    int ok, symmetric, nrow, ncol ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    Common->status = CHOLMOD_OK ;
    if (L->ftype <= CHOLMOD_SYMBOLIC_SUPER || L->ftype >= CHOLMOD_LL_SUPER)
    {
	/* cannot operate on a supernodal factorization */
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_resymbol: cannot operate on supernodal L", Common) ;
	return (FALSE) ;
    }
    if (L->n != A->nrow)
    {
	/* dimensions must agree */
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_resymbol: A and L dimensions do not match", Common) ;
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    symmetric = A->stype ;
    nrow = A->nrow ;
    ncol = A->ncol ;
    cholmod_allocate_work (nrow, 2*nrow + (symmetric ? 0 : ncol), 0, 0, Common);
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* permute the input matrix if necessary */
    /* ---------------------------------------------------------------------- */

    H = NULL ;
    G = NULL ;

    if (symmetric > 0)
    {
	if (L->ordering == CHOLMOD_NATURAL)
	{
	    /* F = triu(A)' */
	    /* workspace: Iwork (nrow) */
	    G = cholmod_transpose (A, FALSE, NULL, NULL, 0, Common) ;
	}
	else
	{
	    /* F = triu(A(p,p))' */
	    /* workspace: Iwork (2*nrow) */
	    G = cholmod_transpose (A, FALSE, L->Perm, NULL, 0, Common) ;
	}
	F = G ;
    }
    else if (symmetric < 0)
    {
	if (L->ordering == CHOLMOD_NATURAL)
	{
	    F = A ;
	}
	else
	{
	    /* G = triu(A(p,p))' */
	    /* workspace: Iwork (2*nrow) */
	    G = cholmod_transpose (A, FALSE, L->Perm, NULL, 0, Common) ;
	    /* H = G' */
	    /* workspace: Iwork (nrow) */
	    H = cholmod_transpose (G, FALSE, NULL, NULL, 0, Common) ;
	    F = H ;
	}
    }
    else
    {
	if (L->ordering == CHOLMOD_NATURAL)
	{
	    F = A ;
	}
	else
	{
	    /* G = A(p,f)' */
	    /* workspace: Iwork (nrow if no fset; MAX (nrow,ncol) if fset)*/
	    G = cholmod_transpose (A, FALSE, L->Perm, fset, fsize, Common) ;
	    /* H = G' */
	    /* workspace: Iwork (ncol) */
	    H = cholmod_transpose (G, FALSE, NULL, NULL, 0, Common) ;
	    F = H ;
	}
    }

    /* No need to check for failure here.  cholmod_resymbol_noperm will return
     * FALSE if F is NULL. */

    /* ---------------------------------------------------------------------- */
    /* resymbol */
    /* ---------------------------------------------------------------------- */

    ok = cholmod_resymbol_noperm (F, fset, fsize, L, pack, Common) ;

    /* ---------------------------------------------------------------------- */
    /* free the temporary matrices, if they exist */
    /* ---------------------------------------------------------------------- */

    cholmod_free_sparse (&H, Common) ;
    cholmod_free_sparse (&G, Common) ;
    return (ok) ;
}


/* ========================================================================== */
/* === cholmod_resymbol_noperm ============================================== */
/* ========================================================================== */

/* Redo symbolic LDL' or LL' factorization of I + F*F' or I+A, where F=A(:,f).
 *
 * L already exists, but is a superset of the true dynamic pattern (simple
 * column downdates and row deletions haven't pruned anything).  Just redo the
 * symbolic factorization and drop entries that are no longer there.  The
 * diagonal is not modified.  The number of nonzeros in column j of L
 * (L->nz[j]) can decrease, but the column pointers (L->p[j]) remain unchanged.
 * L is in unpacked form.
 *
 * For the symmetric case, the columns of the lower triangular part of A
 * are accessed by column.  NOTE that this the transpose of the general case.
 *
 * For the unsymmetric case, F=A(:,f) is accessed by column.
 *
 * A need not be sorted, and can be packed or unpacked.  If L->Perm is not
 * identity, then A must already be permuted according to the permutation used
 * to factorize L.  The advantage of using this routine is that it does not
 * need to create permuted copies of A first.
 *
 * This routine can be called if L is only partially factored via cholmod_rowfac
 * since all it does is prune.  If an entry is in F*F' or A, but not in L, it
 * isn't added to L.
 *
 * L can be packed, unpacked, or dynamic.  If L is unpacked
 * (but not dynamic) then the output L can be packed or unpacked, depending
 * on the "pack" input parameter.  If L is packed or dynamic, it remains in
 * that state on output, regardless of the "pack" parameter.
 *
 * L must be simplicial LDL' or LL'; it cannot be supernodal.
 *
 * The set f is held in fset and fsize.
 *	fset = NULL means ":" in MATLAB. fset is ignored.
 *	fset != NULL means f = fset [0..fset-1].
 *	fset != NULL and fsize = 0 means f is the empty set.
 *	There can be no duplicates in fset.
 *	Common->status is set to CHOLMOD_INVALID if fset is invalid.
 *
 * workspace: Flag (nrow), Head (nrow+1),
 *	if symmetric:   Iwork (2*nrow)
 *	if unsymmetric: Iwork (2*nrow+ncol).
 *	Unlike cholmod_resymbol, this routine does not allocate any temporary
 *	copies of its input matrix.
 */

int cholmod_resymbol_noperm	/* TRUE if successful, FALSE otherwise */
(
    /* inputs, not modified */
    cholmod_sparse *A,	/* nrow-by-ncol, stored in column form */
    void *fset_p,	/* if present and A->stype=0, then LDL'=F*F' is
			 * analyzed, where F = A (:, fset [0..fsize-1]).
			 * Otherwise, F=A. */
    size_t fsize,	/* number of columns in the set f */

    /* The factors L and D.  Some terms are inputs, some both input & output */
    cholmod_factor *L,	/* nrow-by-nrow */
    int pack,		/* if TRUE and L not dynamic, resulting L is packed */

    /* workspace: */
    cholmod_common *Common	/* uses Flag, Head, Iwork (nrow+ncol) */
)
{
    double *Lx ;
    int i, j, k, row, parent, p, pend, pdest, ncol, apacked, sorted, nrow, nf,
	use_fset, mark, jj, symmetric, lpacked ;
    int *Ap, *Ai, *Anz, *Li, *Lp, *Lnz, *Flag, *Head, *Link, *Anext, *fset,
	*Iwork ;

    fset = fset_p ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    ncol = A->ncol ;
    nrow = A->nrow ;
    symmetric = A->stype ;
    ASSERT (IMPLIES (symmetric, nrow == ncol)) ;
    if (symmetric > 0)
    {
	/* symmetric, with upper triangular part, not supported */
	/* TODO: transpose A */
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_resymbol: symmetric upper not supported ", Common) ;
	return (FALSE) ;
    }

    if (L->ftype <= CHOLMOD_SYMBOLIC_SUPER || L->ftype >= CHOLMOD_LL_SUPER)
    {
	/* cannot operate on a supernodal factorization */
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_resymbol: cannot operate on supernodal L", Common) ;
	return (FALSE) ;
    }
    if (L->n != A->nrow)
    {
	/* dimensions must agree */
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_resymbol: A and L dimensions do not match", Common) ;
	return (FALSE) ;
    }

    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    cholmod_allocate_work (nrow, 2*nrow + (symmetric ? 0 : ncol), 0, 0, Common);
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;	/* out of memory */
    }
    ASSERT (cholmod_dump_work (TRUE, TRUE, 0, Common)) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    Ai = A->i ;
    Ap = A->p ;
    Anz = A->nz ;
    apacked = A->packed ;
    sorted = A->sorted ;

    Li = L->i ;
    Lx = L->x ;
    Lp = L->p ;

    lpacked = (L->ftype == CHOLMOD_LDL_PACKED ||
	       L->ftype == CHOLMOD_LL_PACKED) ;

    if (lpacked)
    {
	/* if L is already packed on input, it remains packed on output */
	pack = TRUE ;
    }
    else if (L->ftype == CHOLMOD_LDL_DYNAMIC)
    {
	/* cannot pack a dynamic matrix */
	pack = FALSE ;
    }

    ASSERT (L->nzmax >= (size_t) (Lp [L->n])) ;

    /* If L is of type CHOLMOD_LDL_UNPACKED on input, then it can be packed or
     * unpacked on output, depending on the pack input parameter. */

    pdest = 0 ;

    PRINT1 (("\n\n===================== Resymbol pack %d Lwhat %d Apacked %d\n",
	pack, L->ftype, A->packed)) ;
    ASSERT (cholmod_dump_sparse (A, "ReSymbol A:", Common) >= 0) ;
    DEBUG (cholmod_dump_factor (L, "ReSymbol initial L (i, x):", Common)) ;

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    Flag  = Common->Flag ;	/* size nrow */
    Head  = Common->Head ;	/* size nrow+1 */
    Iwork = Common->Iwork ;
    Link  = Iwork ;		/* size nrow (i/i/l) [ */
    Lnz   = Iwork + nrow ;	/* size nrow (i/i/l), if L not packed */
    Anext = Iwork + 2*nrow ;	/* size ncol (i/i/l), unsym. only */
    for (j = 0 ; j < nrow ; j++)
    {
	Link [j] = EMPTY ;
    }
    if (lpacked)
    {
	for (j = 0 ; j < nrow ; j++)
	{
	    Lnz [j] = Lp [j+1] - Lp [j] ;
	    ASSERT (Lnz [j] > 0) ;
	}
    }
    else
    {
	/* use Lnz in L itself */
	Lnz = L->nz ;
	ASSERT (Lnz != NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* for the unsymmetric case, queue each column of A (:,f) */
    /* ---------------------------------------------------------------------- */

    /* place each column of the basis set on the link list corresponding to */
    /* the smallest row index in that column */

    if (!symmetric)
    {
	use_fset = (fset != NULL) ;
	if (use_fset)
	{
	    nf = fsize ;
	    /* This is the only O(ncol) loop in cholmod_resymbol.
	     * It is required only to check the fset. */
	    for (j = 0 ; j < ncol ; j++)
	    {
		Anext [j] = -2 ;
	    }
	    for (jj = 0 ; jj < nf ; jj++)
	    {
		j = fset [jj] ;
		if (j < 0 || j > ncol || Anext [j] != -2)
		{
		    /* out-of-range or duplicate entry in fset */
		    cholmod_error (CHOLMOD_INVALID,
			    "cholmod_resymbol: fset invalid", Common) ;
		    ASSERT (cholmod_dump_work (TRUE, TRUE, 0, Common)) ;
		    return (FALSE) ;
		}
		/* flag column j as having been seen */
		Anext [j] = EMPTY ;
	    }
	    /* the fset is now valid */
	    ASSERT (cholmod_dump_perm (fset, nf, ncol, "fset", Common)) ;
	}
	else
	{
	    nf = ncol ;
	}
	for (jj = 0 ; jj < nf ; jj++)
	{
	    j = (use_fset) ? (fset [jj]) : jj ;
	    /* column j is the fset; find the smallest row (if any) */
	    p = Ap [j] ;
	    pend = (apacked) ? (Ap [j+1]) : (p + Anz [j]) ;
	    if (pend > p)
	    {
		k = Ai [p] ;
		if (!sorted)
		{
		    for ( ; p < pend ; p++)
		    {
			k = MIN (k, Ai [p]) ;
		    }
		}
		/* place column j on link list k */
		ASSERT (k >= 0 && k < nrow) ;
		Anext [j] = Head [k] ;
		Head [k] = j ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* recompute symbolic LDL' factorization */
    /* ---------------------------------------------------------------------- */

    for (k = 0 ; k < nrow ; k++)
    {

#ifndef NDEBUG
	PRINT1 (("\n\n================== Initial column k = %d\n", k)) ;
	for (p = Lp [k] ; p < Lp [k] + Lnz [k] ; p++)
	{
	    PRINT1 ((" row: %d  value: %e\n", Li [p], Lx [p])) ;
	}
	PRINT1 (("Recomputing LDL, column k = %d\n", k)) ;
#endif

	/* ------------------------------------------------------------------ */
	/* compute column k of I+F*F' or I+A */
	/* ------------------------------------------------------------------ */

	/* flag the diagonal entry */
	mark = cholmod_clear_flag (Common) ;
	Flag [k] = mark ;
	PRINT1 (("	row: %d (diagonal)\n", k)) ;

	if (symmetric)
	{
	    /* merge column k of A into Flag (lower triangular part only) */
	    p = Ap [k] ;
	    pend = (apacked) ? (Ap [k+1]) : (p + Anz [k]) ;
	    for ( ; p < pend ; p++)
	    {
		i = Ai [p] ;
		if (i > k)
		{
		    Flag [i] = mark ;
		}
	    }
	}
	else
	{
	    /* for each column j whos first row index is in row k */
	    for (j = Head [k] ; j != EMPTY ; j = Anext [j])
	    {
		/* merge column j of A into Flag */
		PRINT1 (("	---- A column %d\n", j)) ;
		p = Ap [j] ;
		pend = (apacked) ? (Ap [j+1]) : (p + Anz [j]) ;
		PRINT1 (("  length %d  adding\n", pend-p)) ;
		for ( ; p < pend ; p++)
		{
#ifndef NDEBUG
		    ASSERT (Ai [p] >= k && Ai [p] < nrow) ;
		    if (Flag [Ai [p]] < mark) PRINT1 (("  row: %d\n", Ai [p])) ;
#endif
		    Flag [Ai [p]] = mark ;
		}
	    }
	    /* clear the kth link list */
	    Head [k] = EMPTY ;
	}

	/* ------------------------------------------------------------------ */
	/* compute pruned pattern of kth column of L = union of children */
	/* ------------------------------------------------------------------ */

	/* for each column j of L whose parent is k */
	for (j = Link [k] ; j != EMPTY ; j = Link [j])
	{
	    /* merge column j of L into Flag */
	    PRINT1 (("	---- L column %d\n", k)) ;
	    ASSERT (j < k) ;
	    ASSERT (Lnz [j] > 0) ;
	    p = Lp [j] ;
	    pend = p + Lnz [j] ;
	    ASSERT (Li [p] == j && Li [p+1] == k) ;
	    p++ ;	    /* skip past the diagonal entry */
	    for ( ; p < pend ; p++)
	    {
		/* add to pattern */
#ifndef NDEBUG
		ASSERT (Li [p] >= k && Li [p] < nrow) ;
		if (Flag [Li [p]] < mark) PRINT1 (("  row: %d\n", Li [p])) ;
#endif
		Flag [Li [p]] = mark ;
	    }
	}

	/* ------------------------------------------------------------------ */
	/* prune the kth column of L */
	/* ------------------------------------------------------------------ */

	PRINT1 (("Final column of L:\n")) ;
	p = Lp [k] ;
	pend = p + Lnz [k] ;

	if (pack)
	{
	    /* shift column k upwards */
	    Lp [k] = pdest ;
	}
	else
	{
	    /* leave column k in place, just reduce Lnz [k] */
	    pdest = p ;
	}

	for ( ; p < pend ; p++)
	{
	    ASSERT (pdest < pend) ;
	    ASSERT (pdest <= p) ;
	    row = Li [p] ;
	    ASSERT (row >= k && row < nrow) ;
	    if (Flag [row] == mark)
	    {
		/* keep this entry */
		PRINT1 ((" row: %d  value: %e\n", row, Lx [p])) ;
		Li [pdest] = row ;
		Lx [pdest] = Lx [p] ;
		pdest++ ;
	    }
#ifndef NDEBUG
	    else
	    {
		PRINT1 ((" row: %d  value: %e  DROPPED ", row, Lx [p])) ;
		PRINT1 (("\n")) ;
	    }
#endif
	}

	/* ------------------------------------------------------------------ */
	/* prepare this column for its parent */
	/* ------------------------------------------------------------------ */

	Lnz [k] = pdest - Lp [k] ;

	PRINT1 ((" L(%d) length %d\n", k, Lnz [k])) ;
	ASSERT (Lnz [k] > 0) ;

	/* parent is the first entry in the column after the diagonal */
	parent = (Lnz [k] > 1) ? (Li [Lp [k] + 1]) : EMPTY ;

	PRINT1 (("parent (%d) = %d\n", k, parent)) ;
	ASSERT ((parent > k && parent < nrow) || (parent == EMPTY)) ;

	if (parent != EMPTY)
	{
	    Link [k] = Link [parent] ;
	    Link [parent] = k ;
	}
    }

    /* done using Iwork for Link, Lnz (if needed), and Anext ] */

    /* ---------------------------------------------------------------------- */
    /* convert L to packed, if requested */
    /* ---------------------------------------------------------------------- */

    if (pack)
    {
	/* finalize Lp */
	Lp [nrow] = pdest ;
	ASSERT (L->ftype == CHOLMOD_LDL_PACKED ||
		L->ftype == CHOLMOD_LDL_UNPACKED ||
		L->ftype == CHOLMOD_LL_PACKED) ;

	if (L->ftype == CHOLMOD_LDL_UNPACKED)
	{
	    /* This conversion takes O(n) time since the columns of L are
	     * already packed and monotonic.  It cannot fail.  It resizes L
	     * to be just large enough. */
	    /* workspace: none */
	    cholmod_change_ftype (L, CHOLMOD_LDL_PACKED, Common) ;
	}
	else
	{
	    /* Shrink L to be just large enough.  It cannot fail. */
	    /* workspace: none */
	    ASSERT ((size_t) (Lp [nrow]) <= L->nzmax) ;
	    cholmod_reallocate_factor (L, Lp [nrow], Common) ;
	}
	ASSERT (L->ftype == CHOLMOD_LDL_PACKED
	     || L->ftype == CHOLMOD_LL_PACKED) ;
	ASSERT (Common->status >= CHOLMOD_OK) ;
    }

    /* ---------------------------------------------------------------------- */
    /* clear workspace */
    /* ---------------------------------------------------------------------- */

    cholmod_clear_flag (Common) ;
    DEBUG (cholmod_dump_factor (L, "ReSymbol final L (i, x):", Common)) ;
    ASSERT (cholmod_dump_work (TRUE, TRUE, 0, Common)) ;
    return (TRUE) ;
}
