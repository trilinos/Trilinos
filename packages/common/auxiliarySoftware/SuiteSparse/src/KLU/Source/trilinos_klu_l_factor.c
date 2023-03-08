/* ========================================================================== */
/* === TRILINOS_KLU_factor =========================================================== */
/* ========================================================================== */

/* Factor the matrix, after ordering and analyzing it with TRILINOS_KLU_analyze
 * or TRILINOS_KLU_analyze_given.
 */

/* This file should make the long int version of KLU */
#define DLONG 1

#include "trilinos_klu_internal.h"

/* ========================================================================== */
/* === KLU_factor2 ========================================================== */
/* ========================================================================== */

static void factor2
(
    /* inputs, not modified */
    Int Ap [ ],		/* size n+1, column pointers */
    Int Ai [ ],		/* size nz, row indices */
    Entry Ax [ ],
    TRILINOS_KLU_symbolic *Symbolic,

    /* inputs, modified on output: */
    TRILINOS_KLU_numeric *Numeric,
    TRILINOS_KLU_common *Common
)
{
    double lsize ;
    double *Lnz, *Rs ;
    Int *P, *Q, *R, *Pnum, *Offp, *Offi, *Pblock, *Pinv, *Iwork,
	*Lip, *Uip, *Llen, *Ulen ;
    Entry *Offx, *X, s, *Udiag ;
    Unit **LUbx ;
    Int k1, k2, nk, k, block, oldcol, pend, oldrow, n, lnz, unz, p, newrow,
	nblocks, poff, nzoff, lnz_block, unz_block, scale, max_lnz_block,
	max_unz_block ;

    /* ---------------------------------------------------------------------- */
    /* initializations */
    /* ---------------------------------------------------------------------- */

    /* get the contents of the Symbolic object */
    n = Symbolic->n ;
    P = Symbolic->P ;
    Q = Symbolic->Q ;
    R = Symbolic->R ;
    Lnz = Symbolic->Lnz ;
    nblocks = Symbolic->nblocks ;
    nzoff = Symbolic->nzoff ;

    Pnum = Numeric->Pnum ;
    Offp = Numeric->Offp ;
    Offi = Numeric->Offi ;
    Offx = (Entry *) Numeric->Offx ;

    Lip = Numeric->Lip ;
    Uip = Numeric->Uip ;
    Llen = Numeric->Llen ;
    Ulen = Numeric->Ulen ;
    LUbx = (Unit **) Numeric->LUbx ;
    Udiag = (double*) Numeric->Udiag ;

    Rs = Numeric->Rs ;
    Pinv = Numeric->Pinv ;
    X = (Entry *) Numeric->Xwork ;		/* X is of size n */
    Iwork = Numeric->Iwork ;			/* 5*maxblock for TRILINOS_KLU_factor */
						/* 1*maxblock for Pblock */
    Pblock = Iwork + 5*((size_t) Symbolic->maxblock) ;
    Common->nrealloc = 0 ;
    scale = Common->scale ;
    max_lnz_block = 1 ;
    max_unz_block = 1 ;

    /* compute the inverse of P from symbolic analysis.  Will be updated to
     * become the inverse of the numerical factorization when the factorization
     * is done, for use in TRILINOS_KLU_refactor */
#ifndef NDEBUG
    for (k = 0 ; k < n ; k++)
    {
	Pinv [k] = TRILINOS_KLU_EMPTY ;
    }
#endif
    for (k = 0 ; k < n ; k++)
    {
	ASSERT (P [k] >= 0 && P [k] < n) ;
	Pinv [P [k]] = k ;
    }
#ifndef NDEBUG
    for (k = 0 ; k < n ; k++) ASSERT (Pinv [k] != TRILINOS_KLU_EMPTY) ;
#endif

    lnz = 0 ;
    unz = 0 ;
    Common->noffdiag = 0 ;
    Offp [0] = 0 ;

    /* ---------------------------------------------------------------------- */
    /* optionally check input matrix and compute scale factors */
    /* ---------------------------------------------------------------------- */

    if (scale >= 0)
    {
	/* use Pnum as workspace. NOTE: scale factors are not yet permuted
	 * according to the final pivot row ordering, so Rs [oldrow] is the
	 * scale factor for A (oldrow,:), for the user's matrix A.  Pnum is
	 * used as workspace in TRILINOS_KLU_scale.  When the factorization is done,
	 * the scale factors are permuted according to the final pivot row
	 * permutation, so that Rs [k] is the scale factor for the kth row of
	 * A(p,q) where p and q are the final row and column permutations. */
	TRILINOS_KLU_scale (scale, n, Ap, Ai, (double *) Ax, Rs, Pnum, Common) ;
	if (Common->status < TRILINOS_KLU_OK)
	{
	    /* matrix is invalid */
	    return ;
	}
    }

#ifndef NDEBUG
    if (scale > 0)
    {
	for (k = 0 ; k < n ; k++) PRINTF (("Rs [%d] %g\n", k, Rs [k])) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* factor each block using klu */
    /* ---------------------------------------------------------------------- */

    for (block = 0 ; block < nblocks ; block++)
    {

	/* ------------------------------------------------------------------ */
	/* the block is from rows/columns k1 to k2-1 */
	/* ------------------------------------------------------------------ */

	k1 = R [block] ;
	k2 = R [block+1] ;
	nk = k2 - k1 ;
	PRINTF (("FACTOR BLOCK %d, k1 %d k2-1 %d nk %d\n", block, k1,k2-1,nk)) ;

	if (nk == 1)
	{

	    /* -------------------------------------------------------------- */
	    /* singleton case */
	    /* -------------------------------------------------------------- */

	    poff = Offp [k1] ;
	    oldcol = Q [k1] ;
	    pend = Ap [oldcol+1] ;
	    CLEAR (s) ;

	    if (scale <= 0)
	    {
		/* no scaling */
		for (p = Ap [oldcol] ; p < pend ; p++)
		{
		    oldrow = Ai [p] ;
		    newrow = Pinv [oldrow] ;
		    if (newrow < k1)
		    {
			Offi [poff] = oldrow ;
			Offx [poff] = Ax [p] ;
			poff++ ;
		    }
		    else
		    {
			ASSERT (newrow == k1) ;
			PRINTF (("singleton block %d", block)) ;
			PRINT_ENTRY (Ax [p]) ;
			s = Ax [p] ;
		    }
		}
	    }
	    else
	    {
		/* row scaling.  NOTE: scale factors are not yet permuted
		 * according to the pivot row permutation, so Rs [oldrow] is
		 * used below.  When the factorization is done, the scale
		 * factors are permuted, so that Rs [newrow] will be used in
		 * klu_solve, klu_tsolve, and klu_rgrowth */
		for (p = Ap [oldcol] ; p < pend ; p++)
		{
		    oldrow = Ai [p] ;
		    newrow = Pinv [oldrow] ;
		    if (newrow < k1)
		    {
			Offi [poff] = oldrow ;
			/* Offx [poff] = Ax [p] / Rs [oldrow] ; */
			SCALE_DIV_ASSIGN (Offx [poff], Ax [p], Rs [oldrow]) ;
			poff++ ;
		    }
		    else
		    {
			ASSERT (newrow == k1) ;
			PRINTF (("singleton block %d ", block)) ;
			PRINT_ENTRY (Ax[p]) ;
			SCALE_DIV_ASSIGN (s, Ax [p], Rs [oldrow]) ;
		    }
		}
	    }

	    Udiag [k1] = s ;

	    if (IS_ZERO (s))
	    {
		/* singular singleton */
		Common->status = TRILINOS_KLU_SINGULAR ;
		Common->numerical_rank = k1 ;
		Common->singular_col = oldcol ;
		if (Common->halt_if_singular)
		{
		    return ;
		}
	    }

	    Offp [k1+1] = poff ;
	    Pnum [k1] = P [k1] ;
	    lnz++ ;
	    unz++ ;

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* construct and factorize the kth block */
	    /* -------------------------------------------------------------- */

	    if (Lnz [block] < 0)
	    {
		/* TRILINOS_COLAMD was used - no estimate of fill-in */
		/* use 10 times the nnz in A, plus n */
		lsize = -(Common->initmem) ;
	    }
	    else
	    {
		lsize = Common->initmem_amd * Lnz [block] + nk ;
	    }

	    /* allocates 1 arrays: LUbx [block] */
	    Numeric->LUsize [block] = TRILINOS_KLU_kernel_factor (nk, Ap, Ai, Ax, Q,
		    lsize, &LUbx [block], Udiag + k1, Llen + k1, Ulen + k1,
		    Lip + k1, Uip + k1, Pblock, &lnz_block, &unz_block,
		    X, Iwork, k1, Pinv, Rs, Offp, Offi, Offx, Common) ;

	    if (Common->status < TRILINOS_KLU_OK ||
	       (Common->status == TRILINOS_KLU_SINGULAR && Common->halt_if_singular))
	    {
		/* out of memory, invalid inputs, or singular */
		return ;
	    }

	    PRINTF (("\n----------------------- L %d:\n", block)) ;
	    ASSERT (TRILINOS_KLU_valid_LU (nk, TRUE, Lip+k1, Llen+k1, LUbx [block])) ;
	    PRINTF (("\n----------------------- U %d:\n", block)) ;
	    ASSERT (TRILINOS_KLU_valid_LU (nk, FALSE, Uip+k1, Ulen+k1, LUbx [block])) ;

	    /* -------------------------------------------------------------- */
	    /* get statistics */
	    /* -------------------------------------------------------------- */

	    lnz += lnz_block ;
	    unz += unz_block ;
	    max_lnz_block = MAX (max_lnz_block, lnz_block) ;
	    max_unz_block = MAX (max_unz_block, unz_block) ;

	    if (Lnz [block] == TRILINOS_KLU_EMPTY)
	    {
		/* revise estimate for subsequent factorization */
		Lnz [block] = MAX (lnz_block, unz_block) ;
	    }

	    /* -------------------------------------------------------------- */
	    /* combine the klu row ordering with the symbolic pre-ordering */
	    /* -------------------------------------------------------------- */

	    PRINTF (("Pnum, 1-based:\n")) ;
	    for (k = 0 ; k < nk ; k++)
	    {
		ASSERT (k + k1 < n) ;
		ASSERT (Pblock [k] + k1 < n) ;
		Pnum [k + k1] = P [Pblock [k] + k1] ;
		PRINTF (("Pnum (%d + %d + 1 = %d) = %d + 1 = %d\n",
		    k, k1, k+k1+1, Pnum [k+k1], Pnum [k+k1]+1)) ;
	    }

	    /* the local pivot row permutation Pblock is no longer needed */
	}
    }
    ASSERT (nzoff == Offp [n]) ;
    PRINTF (("\n------------------- Off diagonal entries:\n")) ;
    ASSERT (TRILINOS_KLU_valid (n, Offp, Offi, Offx)) ;

    Numeric->lnz = lnz ;
    Numeric->unz = unz ;
    Numeric->max_lnz_block = max_lnz_block ;
    Numeric->max_unz_block = max_unz_block ;

    /* compute the inverse of Pnum */
#ifndef NDEBUG
    for (k = 0 ; k < n ; k++)
    {
	Pinv [k] = TRILINOS_KLU_EMPTY ;
    }
#endif
    for (k = 0 ; k < n ; k++)
    {
	ASSERT (Pnum [k] >= 0 && Pnum [k] < n) ;
	Pinv [Pnum [k]] = k ;
    }
#ifndef NDEBUG
    for (k = 0 ; k < n ; k++) ASSERT (Pinv [k] != TRILINOS_KLU_EMPTY) ;
#endif

    /* permute scale factors Rs according to pivotal row order */
    if (scale > 0)
    {
	for (k = 0 ; k < n ; k++)
	{
	    REAL (X [k]) = Rs [Pnum [k]] ;
	}
	for (k = 0 ; k < n ; k++)
	{
	    Rs [k] = REAL (X [k]) ;
	}
    }

    PRINTF (("\n------------------- Off diagonal entries, old:\n")) ;
    ASSERT (TRILINOS_KLU_valid (n, Offp, Offi, Offx)) ;

    /* apply the pivot row permutations to the off-diagonal entries */
    for (p = 0 ; p < nzoff ; p++)
    {
	ASSERT (Offi [p] >= 0 && Offi [p] < n) ;
	Offi [p] = Pinv [Offi [p]] ;
    }

    PRINTF (("\n------------------- Off diagonal entries, new:\n")) ;
    ASSERT (TRILINOS_KLU_valid (n, Offp, Offi, Offx)) ;

#ifndef NDEBUG
    {
	PRINTF (("\n ############# KLU_BTF_FACTOR done, nblocks %d\n",nblocks));
	Entry ss, *Udiag = Numeric->Udiag ;
	for (block = 0 ; block < nblocks && Common->status == TRILINOS_KLU_OK ; block++)
	{
	    k1 = R [block] ;
	    k2 = R [block+1] ;
	    nk = k2 - k1 ;
	    PRINTF (("\n======================TRILINOS_KLU_factor output: k1 %d k2 %d nk %d\n",k1,k2,nk)) ;
	    if (nk == 1)
	    {
		PRINTF (("singleton  ")) ;
		/* ENTRY_PRINT (singleton [block]) ; */
		ss = Udiag [k1] ;
		PRINT_ENTRY (ss) ;
	    }
	    else
	    {
		Int *Lip, *Uip, *Llen, *Ulen ;
		Unit *LU ;
		Lip = Numeric->Lip + k1 ;
		Llen = Numeric->Llen + k1 ;
		LU = (Unit *) Numeric->LUbx [block] ;
		PRINTF (("\n---- L block %d\n", block));
		ASSERT (TRILINOS_KLU_valid_LU (nk, TRUE, Lip, Llen, LU)) ;
		Uip = Numeric->Uip + k1 ;
		Ulen = Numeric->Ulen + k1 ;
		PRINTF (("\n---- U block %d\n", block)) ;
		ASSERT (TRILINOS_KLU_valid_LU (nk, FALSE, Uip, Ulen, LU)) ;
	    }
	}
    }
#endif
}



/* ========================================================================== */
/* === TRILINOS_KLU_factor =========================================================== */
/* ========================================================================== */

TRILINOS_KLU_numeric *TRILINOS_KLU_factor		/* returns NULL if error, or a valid
				   TRILINOS_KLU_numeric object if successful */
(
    /* --- inputs --- */
    Int Ap [ ],		/* size n+1, column pointers */
    Int Ai [ ],		/* size nz, row indices */
    double Ax [ ],
    TRILINOS_KLU_symbolic *Symbolic,
    /* -------------- */
    TRILINOS_KLU_common *Common
)
{
    Int n, nzoff, nblocks, maxblock, k, ok = TRUE ;
    Int *R ;
    TRILINOS_KLU_numeric *Numeric ;
    size_t n1, nzoff1, s, b6, n3 ;

    if (Common == NULL)
    {
	return (NULL) ;
    }
    Common->status = TRILINOS_KLU_OK ;
    Common->numerical_rank = TRILINOS_KLU_EMPTY ;
    Common->singular_col = TRILINOS_KLU_EMPTY ;

    /* ---------------------------------------------------------------------- */
    /* get the contents of the Symbolic object */
    /* ---------------------------------------------------------------------- */

    /* check for a valid Symbolic object */
    if (Symbolic == NULL)
    {
	Common->status = TRILINOS_KLU_INVALID ;
	return (NULL) ;
    }

    n = Symbolic->n ;
    nzoff = Symbolic->nzoff ;
    nblocks = Symbolic->nblocks ;
    maxblock = Symbolic->maxblock ;
    R = Symbolic->R ;
    PRINTF (("TRILINOS_KLU_factor:  n %d nzoff %d nblocks %d maxblock %d\n",
	n, nzoff, nblocks, maxblock)) ;

    /* ---------------------------------------------------------------------- */
    /* get control parameters and make sure they are in the proper range */
    /* ---------------------------------------------------------------------- */

    Common->initmem_amd = MAX (1.0, Common->initmem_amd) ;
    Common->initmem = MAX (1.0, Common->initmem) ;
    Common->tol = MIN (Common->tol, 1.0) ;
    Common->tol = MAX (0.0, Common->tol) ;
    Common->memgrow = MAX (1.0, Common->memgrow) ;

    /* ---------------------------------------------------------------------- */
    /* allocate the Numeric object  */
    /* ---------------------------------------------------------------------- */

    /* this will not cause size_t overflow (already checked by TRILINOS_KLU_symbolic) */
    n1 = ((size_t) n) + 1 ;
    nzoff1 = ((size_t) nzoff) + 1 ;

    Numeric = (TRILINOS_KLU_numeric*) TRILINOS_KLU_malloc (sizeof (TRILINOS_KLU_numeric), 1, Common) ;
    if (Common->status < TRILINOS_KLU_OK)
    {
	/* out of memory */
	Common->status = TRILINOS_KLU_OUT_OF_MEMORY ;
	return (NULL) ;
    }
    Numeric->n = n ;
    Numeric->nblocks = nblocks ;
    Numeric->nzoff = nzoff ;
    Numeric->Pnum = (Int*) TRILINOS_KLU_malloc (n, sizeof (Int), Common) ;
    Numeric->Offp = (Int*) TRILINOS_KLU_malloc (n1, sizeof (Int), Common) ;
    Numeric->Offi = (Int*) TRILINOS_KLU_malloc (nzoff1, sizeof (Int), Common) ;
    Numeric->Offx = (Entry*) TRILINOS_KLU_malloc (nzoff1, sizeof (Entry), Common) ;

    Numeric->Lip  = (Int*) TRILINOS_KLU_malloc (n, sizeof (Int), Common) ;
    Numeric->Uip  = (Int*) TRILINOS_KLU_malloc (n, sizeof (Int), Common) ;
    Numeric->Llen = (Int*) TRILINOS_KLU_malloc (n, sizeof (Int), Common) ;
    Numeric->Ulen = (Int*) TRILINOS_KLU_malloc (n, sizeof (Int), Common) ;

    Numeric->LUsize = (size_t*) TRILINOS_KLU_malloc (nblocks, sizeof (size_t), Common) ;

    Numeric->LUbx = (void**) TRILINOS_KLU_malloc (nblocks, sizeof (Unit *), Common) ;
    if (Numeric->LUbx != NULL)
    {
	for (k = 0 ; k < nblocks ; k++)
	{
	    Numeric->LUbx [k] = NULL ;
	}
    }

    Numeric->Udiag = TRILINOS_KLU_malloc (n, sizeof (Entry), Common) ;

    if (Common->scale > 0)
    {
	Numeric->Rs = (double*) TRILINOS_KLU_malloc (n, sizeof (double), Common) ;
    }
    else
    {
	/* no scaling */
	Numeric->Rs = NULL ;
    }

    Numeric->Pinv = (Int*) TRILINOS_KLU_malloc (n, sizeof (Int), Common) ;

    /* allocate permanent workspace for factorization and solve.  Note that the
     * solver will use an Xwork of size 4n, whereas the factorization codes use
     * an Xwork of size n and integer space (Iwork) of size 6n. TRILINOS_KLU_condest
     * uses an Xwork of size 2n.  Total size is:
     *
     *    n*sizeof(Entry) + max (6*maxblock*sizeof(Int), 3*n*sizeof(Entry))
     */
    s = TRILINOS_KLU_mult_size_t (n, sizeof (Entry), &ok) ;
    n3 = TRILINOS_KLU_mult_size_t (n, 3 * sizeof (Entry), &ok) ;
    b6 = TRILINOS_KLU_mult_size_t (maxblock, 6 * sizeof (Int), &ok) ;
    Numeric->worksize = TRILINOS_KLU_add_size_t (s, MAX (n3, b6), &ok) ;
    Numeric->Work = TRILINOS_KLU_malloc (Numeric->worksize, 1, Common) ;
    Numeric->Xwork = Numeric->Work ;
    Numeric->Iwork = (Int *) ((Entry *) Numeric->Xwork + n) ;
    if (!ok || Common->status < TRILINOS_KLU_OK)
    {
	/* out of memory or problem too large */
	Common->status = ok ? TRILINOS_KLU_OUT_OF_MEMORY : TRILINOS_KLU_TOO_LARGE ;
	TRILINOS_KLU_free_numeric (&Numeric, Common) ;
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* factorize the blocks */
    /* ---------------------------------------------------------------------- */

    factor2 (Ap, Ai, (Entry *) Ax, Symbolic, Numeric, Common) ;

    /* ---------------------------------------------------------------------- */
    /* return or free the Numeric object */
    /* ---------------------------------------------------------------------- */

    if (Common->status < TRILINOS_KLU_OK)
    {
	/* out of memory or inputs invalid */
	TRILINOS_KLU_free_numeric (&Numeric, Common) ;
    }
    else if (Common->status == TRILINOS_KLU_SINGULAR)
    {
	if (Common->halt_if_singular)
	{
	    /* Matrix is singular, and the Numeric object is only partially
	     * defined because we halted early.  This is the default case for
	     * a singular matrix. */
	    TRILINOS_KLU_free_numeric (&Numeric, Common) ;
	}
    }
    else if (Common->status == TRILINOS_KLU_OK)
    {
	/* successful non-singular factorization */
	Common->numerical_rank = n ;
	Common->singular_col = n ;
    }
    return (Numeric) ;
}
