/* ========================================================================== */
/* === klu_btf_analyze ====================================================== */
/* ========================================================================== */

/* Order the matrix using BTF (or not), and then AMD or COLAMD on the blocks.
 * Does not support using the given ordering (use klu_btf_analyze_given
 * for that case). */

#include "klu_btf_internal.h"

/* ========================================================================== */
/* === klu_btf_analyze_worker =============================================== */
/* ========================================================================== */

/* klu_btf_analyze_worker is not-user callable.  See klu_btf_analyze below */

static int klu_btf_analyze_worker	/* returns KLU_OK or < 0 if error */
(
    /* inputs, not modified */
    int n,		/* A is n-by-n */
    int Ap [ ],		/* size n+1, column pointers */
    int Ai [ ],		/* size nz, row indices */
    int nblocks,	/* # of blocks */
    int Pbtf [ ],	/* BTF row permutation */
    int Qbtf [ ],	/* BTF col permutation */
    int R [ ],		/* size n+1, but only Rbtf [0..nblocks] is used */
    int ordering,	/* what ordering to use */

    /* output only, not defined on input */
    int P [ ],		/* size n */
    int Q [ ],		/* size n */
    double Lnz [ ],	/* size n, but only Lnz [0..nblocks-1] is used */
    int *p_nzoff,
    double *p_lnz,

    /* workspace, not defined on input or output */
    int Pamd [ ],	/* size maxblock */
    int Cp [ ],		/* size maxblock+1 */
    int Ci [ ],		/* size nz */
    int Ep [ ],		/* size maxblock+1 */
    int Ei [ ],		/* size nz */
    int Pinv [ ],	/* size maxblock */

    double Info [ ]	/* output statistics */ 
)
{
    int *RowCount, *W, k1, k2, nk, k, block, oldcol, pend, row, oldrow, newcol,
	result, pc, p, newrow, col, maxnz, nzoff ;
    double amd_Info [AMD_INFO], lnz, lnz1, flops ;

    /* ---------------------------------------------------------------------- */
    /* initializations */
    /* ---------------------------------------------------------------------- */

    /* use Pamd as workspace for RowCount and W */
    RowCount = Pamd ;
    W = Pamd ;

    /* compute the inverse of Pbtf */
#ifndef NDEBUG
    for (k = 0 ; k < n ; k++)
    {
	P [k] = EMPTY ;
	Q [k] = EMPTY ;
	Pinv [k] = EMPTY ;
    }
#endif
    for (k = 0 ; k < n ; k++)
    {
	ASSERT (Pbtf [k] >= 0 && Pbtf [k] < n) ;
	Pinv [Pbtf [k]] = k ;
    }
#ifndef NDEBUG
    for (k = 0 ; k < n ; k++) ASSERT (Pinv [k] != EMPTY) ;
#endif
    nzoff = 0 ;
    lnz = 0 ;
    maxnz = 0 ;
    flops = 0 ;

    /* ---------------------------------------------------------------------- */
    /* order each block using AMD or COLAMD */
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

	if (nk == 1)
	{

	    /* -------------------------------------------------------------- */
	    /* singleton case */
	    /* -------------------------------------------------------------- */

	    Lnz [block] = 1 ;
	    P [k1] = Pbtf [k1] ;
	    Q [k1] = Qbtf [k1] ;
	    oldcol = Q [k1] ;
	    pend = Ap [oldcol+1] ;

	    /* nzoff += ((Ap [oldcol+1] - Ap [oldcol]) - 1) ; */

	    for (p = Ap [oldcol] ; p < pend ; p++)
	    {
		oldrow = Ai [p] ;
		newrow = Pinv [oldrow] ;
		if (newrow < k1)
		{
		    nzoff++ ;
		}
		else
		{
		    lnz++ ;
		}
	    }

	    PRINTF (("nzoff so far %d\n", nzoff)) ;
	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* construct the kth block, C */
	    /* -------------------------------------------------------------- */

	    for (row = 0 ; row < nk ; row++)
	    {
		RowCount [row] = 0 ;
	    }
	    pc = 0 ;
	    for (k = k1 ; k < k2 ; k++)
	    {
		newcol = k-k1 ;
		Cp [newcol] = pc ;
		oldcol = Qbtf [k] ;
		pend = Ap [oldcol+1] ;
		for (p = Ap [oldcol] ; p < pend ; p++)
		{
		    oldrow = Ai [p] ;
		    newrow = Pinv [oldrow] ;
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
			RowCount [newrow]++ ;
		    }
		}
	    }
	    Cp [nk] = pc ;
	    maxnz = MAX (maxnz, pc) ;
	    ASSERT (klu_valid (nk, Cp, Ci, (double *) NULL)) ;

	    /* -------------------------------------------------------------- */
	    /* order the block C */
	    /* -------------------------------------------------------------- */

	    if (ordering == KLU_BTF_CONTROL_USE_AMD)
	    {

		/* ---------------------------------------------------------- */
		/* order the block with AMD (C+C') */
		/* ---------------------------------------------------------- */

		/* since AMD requires a sorted input, compute E = C' */
		Ep [0] = 0 ;
		for (row = 0 ; row < nk ; row++)
		{
		    Ep [row+1] = Ep [row] + RowCount [row] ;
		}
		ASSERT (Ep [nk] == pc) ;
		for (row = 0 ; row < nk ; row++)
		{
		    W [row] = Ep [row] ;
		}
		for (col = 0 ; col < nk ; col++)
		{
		    pend = Cp [col+1] ;
		    for (p = Cp [col] ; p < pend ; p++)
		    {
			ASSERT (W [Ci [p]] >= 0 && W [Ci [p]] < pc) ;
			Ei [W [Ci [p]]++] = col ;
		    }
		}
		ASSERT (klu_valid (nk, Ep, Ei, (double *) NULL)) ;

		PRINTF (("calling AMD\n")) ;
		result = amd_order (nk, Ep, Ei, Pamd, (double *)NULL, amd_Info);
		PRINTF (("AMD done\n")) ;
		if (result == AMD_OUT_OF_MEMORY)
		{
		    return (KLU_OUT_OF_MEMORY) ;
		}
		else if (result == AMD_INVALID)
		{
		    /* fatal error - something is corrupted in the matrix E */
		    PRINTF (("AMD invalid!\n")) ;
		    return (KLU_INVALID) ;
		}

		/* get the ordering statistics from AMD */
		lnz1 = (int) (amd_Info [AMD_LNZ]) + nk ;
		Lnz [block] = lnz1 ;
		lnz += lnz1 ;
		if (pc == maxnz)
		{
		    /* get the symmetry of the biggest block */
		    Info [KLU_BTF_INFO_SYMMETRY] = amd_Info [AMD_SYMMETRY] ;
		}
		flops += 2 * amd_Info [AMD_NMULTSUBS_LU] + amd_Info [AMD_NDIV] ;

	    }
	    else
	    {

		/* ---------------------------------------------------------- */
		/* order the block with COLAMD (C) */
		/* ---------------------------------------------------------- */

		int Alen, *AA, ok, cstats [COLAMD_STATS] ;

		PRINTF (("calling COLAMD\n")) ;

		Alen = colamd_recommended (pc, n, n) ;
		AA = (int *) ALLOCATE (Alen * sizeof (int)) ;
		if (AA == (int *) NULL)
		{
		    PRINTF (("out of memory for colamd\n")) ;
		    return (KLU_OUT_OF_MEMORY) ;
		}

		/* copy the matrix C into AA and Ep */
		for (k = 0 ; k <= nk ; k++)
		{
		    Ep [k] = Cp [k] ;
		}
		for (p = 0 ; p < pc ; p++)
		{
		    /* TODO: use Ei instead of AA */
		    AA [p] = Ci [p] ;
		}

		/* order (and destroy) AA, returning col permutation in Ep */
		ok = colamd (nk, nk, Alen, AA, Ep, (double *) NULL, cstats) ;
		PRINTF (("COLAMD done\n")) ;

		/* free the workspace */
		FREE (AA, int) ;

		if (!ok)
		{
		    /* fatal error - something is corrupted in the matrix */
		    PRINTF (("COLAMD invalid!\n")) ;
		    return (KLU_INVALID) ;
		}

		/* copy the permutation from Ep to Pamd */
		for (k = 0 ; k < nk ; k++)
		{
		    Pamd [k] = Ep [k] ;
		}

		/* COLAMD does not return an estimate on fill-in or flops */
		flops = EMPTY ;
		lnz = EMPTY ;
		Lnz [block] = EMPTY ;

		PRINTF (("done COLAMD\n")) ;
	    }

	    /* -------------------------------------------------------------- */
	    /* combine the preordering with the BTF ordering */
	    /* -------------------------------------------------------------- */

	    PRINTF (("Pamd, 1-based:\n")) ;
	    for (k = 0 ; k < nk ; k++)
	    {
		ASSERT (k + k1 < n) ;
		ASSERT (Pamd [k] + k1 < n) ;
		Q [k + k1] = Qbtf [Pamd [k] + k1] ;
	    }
	    for (k = 0 ; k < nk ; k++)
	    {
		ASSERT (k + k1 < n) ;
		ASSERT (Pamd [k] + k1 < n) ;
		P [k + k1] = Pbtf [Pamd [k] + k1] ;
	    }
	}
    }

    PRINTF (("nzoff %d  Ap[n] %d\n", nzoff, Ap [n])) ;
    ASSERT (nzoff >= 0 && nzoff <= Ap [n]) ;

    /* return estimates of # of nonzeros in L including diagonal */
    *p_lnz = lnz ;				/* EMPTY if COLAMD used */
    *p_nzoff = nzoff ;
    Info [KLU_BTF_INFO_EST_FLOPS] = flops ;	/* EMPTY if COLAMD used */
    return (KLU_OK) ;
}

/* ========================================================================== */
/* === klu_btf_analyze ====================================================== */
/* ========================================================================== */

klu_symbolic *klu_btf_analyze	/* returns NULL if error, or a valid
				   klu_symbolic object if successful */
(
    /* inputs, not modified */
    int n,		/* A is n-by-n */
    int Ap [ ],		/* size n+1, column pointers */
    int Ai [ ],		/* size nz, row indices */

    double Control [KLU_BTF_CONTROL],	/* optional; may be NULL */
    double User_Info [KLU_BTF_INFO]	/* optional; may be NULL */
)
{
    int nblocks, *Qbtf, nz, *Cp, *Ci, *Pinv, *Pamd, *Ep, *Ei, *Pbtf,
	block, maxblock, k1, k2, nk, *P, *Q, *R, nzoff, result, j, i, p,
	pend, do_btf, nzdiag, ordering, k, nfound ;
    klu_symbolic *Symbolic ;
    double *Lnz, lnz, Info2 [KLU_BTF_INFO], *Info ;

    /* ---------------------------------------------------------------------- */
    /* get Info array for statistics, and clear it */
    /* ---------------------------------------------------------------------- */

    Info = (User_Info != (double *) NULL) ? User_Info : Info2 ;
    for (i = 0 ; i < KLU_BTF_INFO ; i++)
    {
	Info [i] = EMPTY ;
    }
    Info [KLU_BTF_INFO_N] = n ;

    /* ---------------------------------------------------------------------- */
    /* determine if input matrix is valid, and get # of nonzeros */
    /* ---------------------------------------------------------------------- */

    /* A is n-by-n, with n > 0.  Ap [0] = 0 and nz = Ap [n] >= 0 required.
     * Ap [j] <= Ap [j+1] must hold for all j = 0 to n-1.  Row indices in Ai
     * must be in the range 0 to n-1, and no duplicate entries can be present.
     * The list of row indices in each column of A need not be sorted.
     */

    if (n <= 0 || (Ap == (int *) NULL) || (Ai == (int *) NULL))
    {
	/* Ap and Ai must be present, and n must be > 0 */
	Info [KLU_BTF_INFO_STATUS] = KLU_INVALID ;
	return ((klu_symbolic *) NULL) ;
    }
    nz = Ap [n] ;
    Info [KLU_BTF_INFO_NZ] = nz ;
    if (Ap [0] != 0 || nz < 0)
    {
	/* nz must be >= 0 and Ap [0] must equal zero */
	Info [KLU_BTF_INFO_STATUS] = KLU_INVALID ;
	return ((klu_symbolic *) NULL) ;
    }
    for (j = 0 ; j < n ; j++)
    {
	if (Ap [j] > Ap [j+1])
	{
	    /* column pointers must be non-decreasing */
	    Info [KLU_BTF_INFO_STATUS] = KLU_INVALID ;
	    return ((klu_symbolic *) NULL) ;
	}
    }
    P = (int *) ALLOCATE (n * sizeof (int)) ;
    if (P == (int *) NULL)
    {
	/* out of memory */
	Info [KLU_BTF_INFO_STATUS] = KLU_OUT_OF_MEMORY ;
	return ((klu_symbolic *) NULL) ;
    }
    for (i = 0 ; i < n ; i++)
    {
	P [i] = EMPTY ;
    }
    nzdiag = 0 ;
    for (j = 0 ; j < n ; j++)
    {
	pend = Ap [j+1] ;
	for (p = Ap [j] ; p < pend ; p++)
	{
	    i = Ai [p] ;
	    if (i < 0 || i >= n || P [i] == j)
	    {
		/* row index out of range, or duplicate entry */
		FREE (P, int) ;
		Info [KLU_BTF_INFO_STATUS] = KLU_INVALID ;
		return ((klu_symbolic *) NULL) ;
	    }
	    if (i == j)
	    {
		/* count the number of diagonal entries */
		/* TODO: move this code to maxtrans */
		nzdiag++ ;
	    }
	    /* flag row i as appearing in column j */
	    P [i] = j ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* allocate the Symbolic object */
    /* ---------------------------------------------------------------------- */

    Symbolic = (klu_symbolic *) ALLOCATE (sizeof (klu_symbolic)) ;
    if (Symbolic == (klu_symbolic *) NULL)
    {
	FREE (P, int) ;
	Info [KLU_BTF_INFO_STATUS] = KLU_OUT_OF_MEMORY ;
	return ((klu_symbolic *) NULL) ;
    }

    Q = (int *) ALLOCATE (n * sizeof (int)) ;
    R = (int *) ALLOCATE ((n+1) * sizeof (int)) ;
    Lnz = (double *) ALLOCATE (n * sizeof (double)) ;

    Symbolic->n = n ;
    Symbolic->P = P ;
    Symbolic->Q = Q ;
    Symbolic->R = R ;
    Symbolic->Lnz = Lnz ;

    if ((P == (int *) NULL) || (Q == (int *) NULL) || (R == (int *) NULL) ||
	(Lnz == (double *) NULL))
    {
	Info [KLU_BTF_INFO_STATUS] = KLU_OUT_OF_MEMORY ;
	klu_btf_free_symbolic (&Symbolic) ;
	return ((klu_symbolic *) NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* allocate workspace for BTF permutation */
    /* ---------------------------------------------------------------------- */

    Pbtf = (int *) ALLOCATE (n * sizeof (int)) ;
    Qbtf = (int *) ALLOCATE (n * sizeof (int)) ;
    if ((Pbtf == (int *) NULL) || (Qbtf == (int *) NULL))
    {
	FREE (Pbtf, int) ; 
	FREE (Qbtf, int) ; 
	klu_btf_free_symbolic (&Symbolic) ;
	Info [KLU_BTF_INFO_STATUS] = KLU_OUT_OF_MEMORY ;
	return ((klu_symbolic *) NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get the Control parameters for BTF and ordering method */
    /* ---------------------------------------------------------------------- */

    do_btf   = (int) GET_CONTROL (KLU_BTF_CONTROL_BTF, TRUE) ;
    ordering = (int) GET_CONTROL (KLU_BTF_CONTROL_ORDERING, KLU_BTF_CONTROL_USE_AMD) ;
    do_btf   = (do_btf) ? TRUE : FALSE ;
    ordering = MAX (ordering, KLU_BTF_CONTROL_USE_AMD) ;	/* 0 */
    ordering = MIN (ordering, KLU_BTF_CONTROL_USE_COLAMD) ;	/* 1 */
    Symbolic->ordering = ordering ;
    Symbolic->do_btf = do_btf ;

    /* ---------------------------------------------------------------------- */
    /* find the block triangular form (if requested) */
    /* ---------------------------------------------------------------------- */

    if (do_btf)
    {

#ifndef HARWELL

	int *Work ;
	Work = (int *) ALLOCATE (5*n * sizeof (int)) ;

	if (Work == (int *) NULL)
	{
	    /* out of memory */
	    FREE (Pbtf, int) ; 
	    FREE (Qbtf, int) ; 
	    klu_btf_free_symbolic (&Symbolic) ;
	    Info [KLU_BTF_INFO_STATUS] = KLU_OUT_OF_MEMORY ;
	    return ((klu_symbolic *) NULL) ;
	}

	nblocks = btf_order (n, Ap, Ai, Pbtf, Qbtf, R, &nfound, Work) ;

	/* TODO: add nfound to Symbolic, and to Info */

	FREE (Work, int) ; 

	/* unflip Qbtf if singular */
	if (nfound < n)
	{
	    for (k = 0 ; k < n ; k++)
	    {
		Qbtf [k] = MAXTRANS_UNFLIP (Qbtf [k]) ;
	    }
	}

#else

	/* call mc21 + mc13 */
	nblocks = charwell (n, Ap, Ai, Pbtf, Qbtf, R) ;
	if (nblocks < 0)
	{
	    /* out of memory */
	    FREE (Pbtf, int) ; 
	    FREE (Qbtf, int) ; 
	    klu_btf_free_symbolic (&Symbolic) ;
	    Info [KLU_BTF_INFO_STATUS] = KLU_OUT_OF_MEMORY ;
	    return ((klu_symbolic *) NULL) ;
	}

#endif

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
	/* BTF not requested */
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
    Info [KLU_BTF_INFO_NBLOCKS] = nblocks ;

    PRINTF (("maxblock size %d\n", maxblock)) ;
    Symbolic->maxblock = maxblock ;
    Info [KLU_BTF_INFO_MAXBLOCK] = maxblock ;

    /* TODO: merge adjacent 1-by-1 blocks into an upper triangular block */

    /* ---------------------------------------------------------------------- */
    /* allocate more workspace, for klu_btf_analyze_worker */
    /* ---------------------------------------------------------------------- */

    Pamd = (int *) ALLOCATE (maxblock * sizeof (int)) ;
    Cp = (int *) ALLOCATE ((maxblock+1) * sizeof (int)) ;
    Ep = (int *) ALLOCATE ((maxblock+1) * sizeof (int)) ;
    Ci = (int *) ALLOCATE ((nz+1) * sizeof (int)) ;
    Ei = (int *) ALLOCATE ((nz+1) * sizeof (int)) ;
    Pinv = (int *) ALLOCATE (n * sizeof (int)) ;

    result = ((Pamd != (int *) NULL) &&
	(Cp != (int *) NULL) &&
	(Ci != (int *) NULL) &&
	(Ep != (int *) NULL) &&
	(Ei != (int *) NULL) &&
	(Pinv != (int *) NULL)) ? KLU_OK : KLU_OUT_OF_MEMORY ;

    /* ---------------------------------------------------------------------- */
    /* order each block of the BTF matrix using AMD or COLAMD */
    /* ---------------------------------------------------------------------- */

    if (result == KLU_OK)
    {
	PRINTF (("calling klu_btf_analyze_worker\n")) ;
	result = klu_btf_analyze_worker (n, Ap, Ai, nblocks, Pbtf, Qbtf, R,
	    ordering, P, Q, Lnz, &nzoff, &lnz, Pamd, Cp, Ci, Ep, Ei,
	    Pinv, Info) ;
	PRINTF (("klu_btf_analyze_worker done\n")) ;
	if (!do_btf) ASSERT (nzoff == 0) ;
    }
    Info [KLU_BTF_INFO_STATUS] = result ;

    /* ---------------------------------------------------------------------- */
    /* free all workspace */
    /* ---------------------------------------------------------------------- */

    FREE (Pamd, int) ;
    FREE (Cp, int) ;
    FREE (Ci, int) ;
    FREE (Ep, int) ;
    FREE (Ei, int) ;
    FREE (Pinv, int) ;
    FREE (Pbtf, int) ;
    FREE (Qbtf, int) ;

    PRINTF (("klu_btf_analyze done\n")) ;

    /* ---------------------------------------------------------------------- */
    /* return the symbolic object */
    /* ---------------------------------------------------------------------- */

    if (result != KLU_OK)
    {
	klu_btf_free_symbolic (&Symbolic) ;
	return ((klu_symbolic *) NULL) ;
    }
    else
    {
	Symbolic->lnz = lnz ;
	Symbolic->unz = lnz ;	/* estimated fill-in is symmetric, currently */
	Symbolic->nzoff = nzoff ;
	Info [KLU_BTF_INFO_EST_LNZ] = lnz ;
	Info [KLU_BTF_INFO_NZOFF] = nzoff ;
	return (Symbolic) ;
    }
}
