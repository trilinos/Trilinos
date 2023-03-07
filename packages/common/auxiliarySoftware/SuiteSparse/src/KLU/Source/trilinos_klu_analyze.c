/* ========================================================================== */
/* === klu_analyze ========================================================== */
/* ========================================================================== */

/* Order the matrix using TRILINOS_BTF (or not), and then AMD, TRILINOS_COLAMD, the natural
 * ordering, or the user-provided-function on the blocks.  Does not support
 * using a given ordering (use klu_analyze_given for that case). */

#include "trilinos_klu_internal.h"

/* ========================================================================== */
/* === analyze_worker ======================================================= */
/* ========================================================================== */

static Int analyze_worker	/* returns TRILINOS_KLU_OK or < 0 if error */
(
    /* inputs, not modified */
    Int n,		/* A is n-by-n */
    Int Ap [ ],		/* size n+1, column pointers */
    Int Ai [ ],		/* size nz, row indices */
    Int nblocks,	/* # of blocks */
    Int Pbtf [ ],	/* TRILINOS_BTF row permutation */
    Int Qbtf [ ],	/* TRILINOS_BTF col permutation */
    Int R [ ],		/* size n+1, but only Rbtf [0..nblocks] is used */
    Int ordering,	/* what ordering to use (0, 1, or 3 for this routine) */

    /* output only, not defined on input */
    Int P [ ],		/* size n */
    Int Q [ ],		/* size n */
    double Lnz [ ],	/* size n, but only Lnz [0..nblocks-1] is used */

    /* workspace, not defined on input or output */
    Int Pblk [ ],	/* size maxblock */
    Int Cp [ ],		/* size maxblock+1 */
    Int Ci [ ],		/* size MAX (nz+1, Cilen) */
    Int Cilen,		/* nz+1, or COLAMD_recommend(nz,n,n) for TRILINOS_COLAMD */
    Int Pinv [ ],	/* size maxblock */

    /* input/output */
    TRILINOS_KLU_symbolic *Symbolic,
    TRILINOS_KLU_common *Common
)
{
    double amd_Info [TRILINOS_AMD_INFO], lnz, lnz1, flops, flops1 ;
    Int k1, k2, nk, k, block, oldcol, pend, newcol, result, pc, p, newrow,
	maxnz, nzoff, cstats [TRILINOS_COLAMD_STATS], ok, err = TRILINOS_KLU_INVALID ;

    /* ---------------------------------------------------------------------- */
    /* initializations */
    /* ---------------------------------------------------------------------- */

    /* compute the inverse of Pbtf */
#ifndef NDEBUG
    for (k = 0 ; k < n ; k++)
    {
	P [k] = TRILINOS_KLU_EMPTY ;
	Q [k] = TRILINOS_KLU_EMPTY ;
	Pinv [k] = TRILINOS_KLU_EMPTY ;
    }
#endif
    for (k = 0 ; k < n ; k++)
    {
	ASSERT (Pbtf [k] >= 0 && Pbtf [k] < n) ;
	Pinv [Pbtf [k]] = k ;
    }
#ifndef NDEBUG
    for (k = 0 ; k < n ; k++) ASSERT (Pinv [k] != TRILINOS_KLU_EMPTY) ;
#endif
    nzoff = 0 ;
    lnz = 0 ;
    maxnz = 0 ;
    flops = 0 ;
    Symbolic->symmetry = TRILINOS_KLU_EMPTY ;	/* only computed by AMD */

    /* ---------------------------------------------------------------------- */
    /* order each block */
    /* ---------------------------------------------------------------------- */

    for (block = 0 ; block < nblocks ; block++)
    {

	/* ------------------------------------------------------------------ */
	/* the block is from rows/columns k1 to k2-1 */
	/* ------------------------------------------------------------------ */

	k1 = R [block] ;
	k2 = R [block+1] ;
	nk = k2 - k1 ;
	PRINTF (("BLOCK %d, k1 %d k2-1 %d nk %d\n", block, k1, k2-1, nk)) ;

	/* ------------------------------------------------------------------ */
	/* construct the kth block, C */
	/* ------------------------------------------------------------------ */

	Lnz [block] = TRILINOS_KLU_EMPTY ;
	pc = 0 ;
	for (k = k1 ; k < k2 ; k++)
	{
	    newcol = k-k1 ;
	    Cp [newcol] = pc ;
	    oldcol = Qbtf [k] ;
	    pend = Ap [oldcol+1] ;
	    for (p = Ap [oldcol] ; p < pend ; p++)
	    {
		newrow = Pinv [Ai [p]] ;
		if (newrow < k1)
		{
		    nzoff++ ;
		}
		else
		{
		    /* (newrow,newcol) is an entry in the block */
		    ASSERT (newrow < k2) ;
		    newrow -= k1 ;
		    Ci [pc++] = newrow ;
		}
	    }
	}
	Cp [nk] = pc ;
	maxnz = MAX (maxnz, pc) ;
	ASSERT (TRILINOS_KLU_valid (nk, Cp, Ci, NULL)) ;

	/* ------------------------------------------------------------------ */
	/* order the block C */
	/* ------------------------------------------------------------------ */

	if (nk <= 3)
	{

	    /* -------------------------------------------------------------- */
	    /* use natural ordering for tiny blocks (3-by-3 or less) */
	    /* -------------------------------------------------------------- */

	    for (k = 0 ; k < nk ; k++)
	    {
		Pblk [k] = k ;
	    }
	    lnz1 = nk * (nk + 1) / 2 ;
	    flops1 = nk * (nk - 1) / 2 + (nk-1)*nk*(2*nk-1) / 6 ;
	    ok = TRUE ;

	}
	else if (ordering == 0)
	{

	    /* -------------------------------------------------------------- */
	    /* order the block with AMD (C+C') */
	    /* -------------------------------------------------------------- */

	    result = TRILINOS_AMD_order (nk, Cp, Ci, Pblk, NULL, amd_Info) ;
	    ok = (result >= TRILINOS_AMD_OK) ;
	    if (result == TRILINOS_AMD_OUT_OF_MEMORY)
	    {
		err = TRILINOS_KLU_OUT_OF_MEMORY ;
	    }

	    /* account for memory usage in AMD */
	    Common->mempeak = ( size_t) ( MAX (Common->mempeak,
		Common->memusage + amd_Info [TRILINOS_AMD_MEMORY]) ) ;

	    /* get the ordering statistics from AMD */
	    lnz1 = (Int) (amd_Info [TRILINOS_AMD_LNZ]) + nk ;
	    flops1 = 2 * amd_Info [TRILINOS_AMD_NMULTSUBS_LU] + amd_Info [TRILINOS_AMD_NDIV] ;
	    if (pc == maxnz)
	    {
		/* get the symmetry of the biggest block */
		Symbolic->symmetry = amd_Info [TRILINOS_AMD_SYMMETRY] ;
	    }

	}
	else if (ordering == 1)
	{

	    /* -------------------------------------------------------------- */
	    /* order the block with TRILINOS_COLAMD (C) */
	    /* -------------------------------------------------------------- */

	    /* order (and destroy) Ci, returning column permutation in Cp.
	     * TRILINOS_COLAMD "cannot" fail since the matrix has already been checked,
	     * and Ci allocated. */

	    ok = TRILINOS_COLAMD (nk, nk, Cilen, Ci, Cp, NULL, cstats) ;
	    lnz1 = TRILINOS_KLU_EMPTY ;
	    flops1 = TRILINOS_KLU_EMPTY ;

	    /* copy the permutation from Cp to Pblk */
	    for (k = 0 ; k < nk ; k++)
	    {
		Pblk [k] = Cp [k] ;
	    }

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* pass the block to the user-provided ordering function */
	    /* -------------------------------------------------------------- */

	    lnz1 = (Common->user_order) (nk, Cp, Ci, Pblk, Common) ;
	    flops1 = TRILINOS_KLU_EMPTY ;
	    ok = (lnz1 != 0) ;
	}

	if (!ok)
	{
	    return (err) ;  /* ordering method failed */
	}

	/* ------------------------------------------------------------------ */
	/* keep track of nnz(L) and flops statistics */
	/* ------------------------------------------------------------------ */

	Lnz [block] = lnz1 ;
	lnz = (lnz == TRILINOS_KLU_EMPTY || lnz1 == TRILINOS_KLU_EMPTY) ? TRILINOS_KLU_EMPTY : (lnz + lnz1) ;
	flops = (flops == TRILINOS_KLU_EMPTY || flops1 == TRILINOS_KLU_EMPTY) ? TRILINOS_KLU_EMPTY : (flops + flops1) ;

	/* ------------------------------------------------------------------ */
	/* combine the preordering with the TRILINOS_BTF ordering */
	/* ------------------------------------------------------------------ */

	PRINTF (("Pblk, 1-based:\n")) ;
	for (k = 0 ; k < nk ; k++)
	{
	    ASSERT (k + k1 < n) ;
	    ASSERT (Pblk [k] + k1 < n) ;
	    Q [k + k1] = Qbtf [Pblk [k] + k1] ;
	}
	for (k = 0 ; k < nk ; k++)
	{
	    ASSERT (k + k1 < n) ;
	    ASSERT (Pblk [k] + k1 < n) ;
	    P [k + k1] = Pbtf [Pblk [k] + k1] ;
	}
    }

    PRINTF (("nzoff %d  Ap[n] %d\n", nzoff, Ap [n])) ;
    ASSERT (nzoff >= 0 && nzoff <= Ap [n]) ;

    /* return estimates of # of nonzeros in L including diagonal */
    Symbolic->lnz = lnz ;	    /* TRILINOS_KLU_EMPTY if TRILINOS_COLAMD used */
    Symbolic->unz = lnz ;
    Symbolic->nzoff = nzoff ;
    Symbolic->est_flops = flops ;   /* TRILINOS_KLU_EMPTY if TRILINOS_COLAMD or user-ordering used */
    return (TRILINOS_KLU_OK) ;
}


/* ========================================================================== */
/* === order_and_analyze ==================================================== */
/* ========================================================================== */

/* Orders the matrix with or with TRILINOS_BTF, then orders each block with AMD, TRILINOS_COLAMD,
 * or the user ordering function.  Does not handle the natural or given
 * ordering cases. */

static TRILINOS_KLU_symbolic *order_and_analyze	/* returns NULL if error, or a valid
					   TRILINOS_KLU_symbolic object if successful */
(
    /* inputs, not modified */
    Int n,		/* A is n-by-n */
    Int Ap [ ],		/* size n+1, column pointers */
    Int Ai [ ],		/* size nz, row indices */
    /* --------------------- */
    TRILINOS_KLU_common *Common
)
{
    double work ;
    TRILINOS_KLU_symbolic *Symbolic ;
    double *Lnz ;
    Int *Qbtf, *Cp, *Ci, *Pinv, *Pblk, *Pbtf, *P, *Q, *R ;
    Int nblocks, nz, block, maxblock, k1, k2, nk, do_btf, ordering, k, Cilen,
	*Work ;

    /* ---------------------------------------------------------------------- */
    /* allocate the Symbolic object, and check input matrix */
    /* ---------------------------------------------------------------------- */

    Symbolic = TRILINOS_KLU_alloc_symbolic (n, Ap, Ai, Common) ;
    if (Symbolic == NULL)
    {
	return (NULL) ;
    }
    P = Symbolic->P ;
    Q = Symbolic->Q ;
    R = Symbolic->R ;
    Lnz = Symbolic->Lnz ;
    nz = Symbolic->nz ;

    ordering = Common->ordering ;
    if (ordering == 1)
    {
	/* TRILINOS_COLAMD */
	Cilen = TRILINOS_COLAMD_recommended (nz, n, n) ;
    }
    else if (ordering == 0 || (ordering == 3 && Common->user_order != NULL))
    {
	/* AMD or user ordering function */
	Cilen = nz+1 ;
    }
    else
    {
	/* invalid ordering */
	Common->status = TRILINOS_KLU_INVALID ;
	TRILINOS_KLU_free_symbolic (&Symbolic, Common) ;
	return (NULL) ;
    }

    /* AMD memory management routines */
    trilinos_amd_malloc  = Common->malloc_memory ;
    trilinos_amd_free    = Common->free_memory ;
    trilinos_amd_calloc  = Common->calloc_memory ;
    trilinos_amd_realloc = Common->realloc_memory ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace for TRILINOS_BTF permutation */
    /* ---------------------------------------------------------------------- */

    Pbtf = (Int*) TRILINOS_KLU_malloc (n, sizeof (Int), Common) ;
    Qbtf = (Int*) TRILINOS_KLU_malloc (n, sizeof (Int), Common) ;
    if (Common->status < TRILINOS_KLU_OK)
    {
	TRILINOS_KLU_free (Pbtf, n, sizeof (Int), Common) ;
	TRILINOS_KLU_free (Qbtf, n, sizeof (Int), Common) ;
	TRILINOS_KLU_free_symbolic (&Symbolic, Common) ;
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get the common parameters for TRILINOS_BTF and ordering method */
    /* ---------------------------------------------------------------------- */

    do_btf = Common->btf ;
    do_btf = (do_btf) ? TRUE : FALSE ;
    Symbolic->ordering = ordering ;
    Symbolic->do_btf = do_btf ;
    Symbolic->structural_rank = TRILINOS_KLU_EMPTY ;

    /* ---------------------------------------------------------------------- */
    /* find the block triangular form (if requested) */
    /* ---------------------------------------------------------------------- */

    Common->work = 0 ;

    if (do_btf)
    {
	Work = (Int*) TRILINOS_KLU_malloc (5*n, sizeof (Int), Common) ;
	if (Common->status < TRILINOS_KLU_OK)
	{
	    /* out of memory */
	    TRILINOS_KLU_free (Pbtf, n, sizeof (Int), Common) ;
	    TRILINOS_KLU_free (Qbtf, n, sizeof (Int), Common) ;
	    TRILINOS_KLU_free_symbolic (&Symbolic, Common) ;
	    return (NULL) ;
	}

	nblocks = BTF_order (n, Ap, Ai, Common->maxwork, &work, Pbtf, Qbtf, R,
		&(Symbolic->structural_rank), Work) ;
	Common->structural_rank = Symbolic->structural_rank ;
	Common->work += work ;

	TRILINOS_KLU_free (Work, 5*n, sizeof (Int), Common) ;

	/* unflip Qbtf if the matrix does not have full structural rank */
	if (Symbolic->structural_rank < n)
	{
	    for (k = 0 ; k < n ; k++)
	    {
		Qbtf [k] = TRILINOS_BTF_UNFLIP (Qbtf [k]) ;
	    }
	}

	/* find the size of the largest block */
	maxblock = 1 ;
	for (block = 0 ; block < nblocks ; block++)
	{
	    k1 = R [block] ;
	    k2 = R [block+1] ;
	    nk = k2 - k1 ;
	    PRINTF (("block %d size %d\n", block, nk)) ;
	    maxblock = MAX (maxblock, nk) ;
	}
    }
    else
    {
	/* TRILINOS_BTF not requested */
	nblocks = 1 ;
	maxblock = n ;
	R [0] = 0 ;
	R [1] = n ;
	for (k = 0 ; k < n ; k++)
	{
	    Pbtf [k] = k ;
	    Qbtf [k] = k ;
	}
    }

    Symbolic->nblocks = nblocks ;

    PRINTF (("maxblock size %d\n", maxblock)) ;
    Symbolic->maxblock = maxblock ;

    /* ---------------------------------------------------------------------- */
    /* allocate more workspace, for analyze_worker */
    /* ---------------------------------------------------------------------- */

    Pblk = (Int*) TRILINOS_KLU_malloc (maxblock, sizeof (Int), Common) ;
    Cp   = (Int*) TRILINOS_KLU_malloc (maxblock + 1, sizeof (Int), Common) ;
    Ci   = (Int*) TRILINOS_KLU_malloc (MAX (Cilen, nz+1), sizeof (Int), Common) ;
    Pinv = (Int*) TRILINOS_KLU_malloc (n, sizeof (Int), Common) ;

    /* ---------------------------------------------------------------------- */
    /* order each block of the TRILINOS_BTF ordering, and a fill-reducing ordering */
    /* ---------------------------------------------------------------------- */

    if (Common->status == TRILINOS_KLU_OK)
    {
	PRINTF (("calling analyze_worker\n")) ;
	Common->status = analyze_worker (n, Ap, Ai, nblocks, Pbtf, Qbtf, R,
	    ordering, P, Q, Lnz, Pblk, Cp, Ci, Cilen, Pinv, Symbolic, Common) ;
	PRINTF (("analyze_worker done\n")) ;
    }

    /* ---------------------------------------------------------------------- */
    /* free all workspace */
    /* ---------------------------------------------------------------------- */

    TRILINOS_KLU_free (Pblk, maxblock, sizeof (Int), Common) ;
    TRILINOS_KLU_free (Cp, maxblock+1, sizeof (Int), Common) ;
    TRILINOS_KLU_free (Ci, MAX (Cilen, nz+1), sizeof (Int), Common) ;
    TRILINOS_KLU_free (Pinv, n, sizeof (Int), Common) ;
    TRILINOS_KLU_free (Pbtf, n, sizeof (Int), Common) ;
    TRILINOS_KLU_free (Qbtf, n, sizeof (Int), Common) ;

    /* ---------------------------------------------------------------------- */
    /* return the symbolic object */
    /* ---------------------------------------------------------------------- */

    if (Common->status < TRILINOS_KLU_OK)
    {
	TRILINOS_KLU_free_symbolic (&Symbolic, Common) ;
    }
    return (Symbolic) ;
}


/* ========================================================================== */
/* === TRILINOS_KLU_analyze ========================================================== */
/* ========================================================================== */

TRILINOS_KLU_symbolic *TRILINOS_KLU_analyze	/* returns NULL if error, or a valid
				   TRILINOS_KLU_symbolic object if successful */
(
    /* inputs, not modified */
    Int n,		/* A is n-by-n */
    Int Ap [ ],		/* size n+1, column pointers */
    Int Ai [ ],		/* size nz, row indices */
    /* -------------------- */
    TRILINOS_KLU_common *Common
)
{

    /* ---------------------------------------------------------------------- */
    /* get the control parameters for TRILINOS_BTF and ordering method */
    /* ---------------------------------------------------------------------- */

    if (Common == NULL)
    {
	return (NULL) ;
    }
    Common->status = TRILINOS_KLU_OK ;
    Common->structural_rank = TRILINOS_KLU_EMPTY ;

    /* ---------------------------------------------------------------------- */
    /* order and analyze */
    /* ---------------------------------------------------------------------- */

    if (Common->ordering == 2)
    {
	/* natural ordering */
	return (TRILINOS_KLU_analyze_given (n, Ap, Ai, NULL, NULL, Common)) ;
    }
    else
    {
	/* order with P and Q */
	return (order_and_analyze (n, Ap, Ai, Common)) ;
    }
}
