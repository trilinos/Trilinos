/* ========================================================================== */
/* === klu_btf_analyze_given ================================================ */
/* ========================================================================== */

/* Given an input permutation P and Q, create the Symbolic object.  BTF can
 * be done to modify the user's P and Q. */

#include "klu_btf_internal.h"

klu_symbolic *klu_btf_analyze_given	/* returns NULL if error, or a valid
					   klu_symbolic object if successful */
(
    /* inputs, not modified */
    int n,		/* A is n-by-n */
    int Ap [ ],		/* size n+1, column pointers */
    int Ai [ ],		/* size nz, row indices */
    int Puser [ ],	/* size n, user's row permutation (may be NULL) */
    int Quser [ ],	/* size n, user's column permutation (may be NULL) */

    double Control [KLU_BTF_CONTROL],	/* optional; may be NULL */
    double User_Info [KLU_BTF_INFO]	/* optional; may be NULL */
)
{
    int nblocks, nz, block, maxblock, *P, *Q, *R, nzoff, j,
	i, p, pend, do_btf, nzdiag, k ;
    klu_symbolic *Symbolic ;
    double *Lnz, Info2 [KLU_BTF_INFO], *Info ;

    PRINTF (("in klu_btf_analyze_given\n")) ;

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

    /* TODO: make this a separate routine.  Uses size-n int workspace. */

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
    PRINTF (("matrix OK\n")) ;

    /* ---------------------------------------------------------------------- */
    /* allocate the Symbolic object */
    /* ---------------------------------------------------------------------- */

    Symbolic = (klu_symbolic *) ALLOCATE (sizeof (klu_symbolic)) ;
    if (Symbolic == (klu_symbolic *) NULL)
    {
	/* out of memory */
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

    if ((Q == (int *) NULL) || (R == (int *) NULL) || (Lnz == (double *) NULL))
    {
	/* out of memory */
	Info [KLU_BTF_INFO_STATUS] = KLU_OUT_OF_MEMORY ;
	klu_btf_free_symbolic (&Symbolic) ;
	return ((klu_symbolic *) NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* Q = Quser, or identity if Quser is NULL */
    /* ---------------------------------------------------------------------- */

    if (Quser == (int *) NULL)
    {
	for (k = 0 ; k < n ; k++)
	{
	    Q [k] = k ;
	}
    }
    else
    {
	for (k = 0 ; k < n ; k++)
	{
	    Q [k] = Quser [k] ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* get the Control parameters for BTF and ordering method */
    /* ---------------------------------------------------------------------- */

    do_btf = (int) GET_CONTROL (KLU_BTF_CONTROL_BTF, TRUE) ;
    do_btf = (do_btf) ? TRUE : FALSE ;
    Symbolic->ordering = KLU_BTF_CONTROL_USE_GIVEN ;
    Symbolic->do_btf = do_btf ;

    /* ---------------------------------------------------------------------- */
    /* find the block triangular form, if requested */
    /* ---------------------------------------------------------------------- */

    if (do_btf)
    {

	/* ------------------------------------------------------------------ */
	/* get workspace for strongcomp */
	/* ------------------------------------------------------------------ */

	int *Pinv, *Work, *Bi, k1, k2, nk, oldcol ;

	Work = (int *) ALLOCATE (4*n * sizeof (int)) ;
	Pinv = (int *) ALLOCATE (n * sizeof (int)) ;
	if (Puser != (int *) NULL)
	{
	    Bi = (int *) ALLOCATE ((nz+1) * sizeof (int)) ;
	}
	else
	{
	    Bi = Ai ;
	}

	if ((Work == (int *) NULL) || (Pinv == (int *) NULL) ||
	    (Bi == (int *) NULL))
	{
	    /* out of memory */
	    FREE (Work, int) ;
	    FREE (Pinv, int) ;
	    if (Puser != (int *) NULL)
	    {
		FREE (Bi, int) ;
	    }
	    Info [KLU_BTF_INFO_STATUS] = KLU_OUT_OF_MEMORY ;
	    klu_btf_free_symbolic (&Symbolic) ;
	    return ((klu_symbolic *) NULL) ;
	}

	/* ------------------------------------------------------------------ */
	/* B = Puser * A */
	/* ------------------------------------------------------------------ */

	if (Puser != (int *) NULL)
	{
	    for (k = 0 ; k < n ; k++)
	    {
		Pinv [Puser [k]] = k ;
	    }
	    for (p = 0 ; p < nz ; p++)
	    {
		Bi [p] = Pinv [Ai [p]] ;
	    }
	}

	/* TODO: strongcomp could take a Pinv input argument, avoid making copy
	 * of the input matrix */

	/* ------------------------------------------------------------------ */
	/* find the strongly-connected components */
	/* ------------------------------------------------------------------ */

	/* modifies Q, and determines P and R */
	PRINTF (("strongcomp start\n")) ;
	nblocks = strongcomp (n, Ap, Bi, Q, P, R, Work) ;
	PRINTF (("strongcomp done\n")) ;

	/* ------------------------------------------------------------------ */
	/* P = P * Puser */
	/* ------------------------------------------------------------------ */

	if (Puser != (int *) NULL)
	{
	    for (k = 0 ; k < n ; k++)
	    {
		Work [k] = Puser [P [k]] ;
	    }
	    for (k = 0 ; k < n ; k++)
	    {
		P [k] = Work [k] ;
	    }
	}

	/* ------------------------------------------------------------------ */
	/* Pinv = inverse of P */
	/* ------------------------------------------------------------------ */

	for (k = 0 ; k < n ; k++)
	{
	    Pinv [P [k]] = k ;
	}

	/* ------------------------------------------------------------------ */
	/* analyze each block */
	/* ------------------------------------------------------------------ */

	nzoff = 0 ;	    /* nz in off-diagonal part */
	maxblock = 1 ;	    /* size of the largest block */

	for (block = 0 ; block < nblocks ; block++)
	{

	    /* -------------------------------------------------------------- */
	    /* the block is from rows/columns k1 to k2-1 */
	    /* -------------------------------------------------------------- */

	    k1 = R [block] ;
	    k2 = R [block+1] ;
	    nk = k2 - k1 ;
	    PRINTF (("BLOCK %d, k1 %d k2-1 %d nk %d\n", block, k1, k2-1, nk)) ;
	    maxblock = MAX (maxblock, nk) ;

	    /* -------------------------------------------------------------- */
	    /* scan the kth block, C */
	    /* -------------------------------------------------------------- */

	    for (k = k1 ; k < k2 ; k++)
	    {
		oldcol = Q [k] ;
		pend = Ap [oldcol+1] ;
		for (p = Ap [oldcol] ; p < pend ; p++)
		{
		    if (Pinv [Ai [p]] < k1)
		    {
			nzoff++ ;
		    }
		}
	    }

	    /* fill-in not estimated (yet...) */
	    Lnz [block] = EMPTY ;
	}

	/* ------------------------------------------------------------------ */
	/* free all workspace */
	/* ------------------------------------------------------------------ */

	FREE (Work, int) ;
	FREE (Pinv, int) ;
	if (Puser != (int *) NULL)
	{
	    FREE (Bi, int) ;
	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* BTF not requested */
	/* ------------------------------------------------------------------ */

	nzoff = 0 ;
	nblocks = 1 ;
	maxblock = n ;
	R [0] = 0 ;
	R [1] = n ;

	/* ------------------------------------------------------------------ */
	/* P = Puser, or identity if Puser is NULL */
	/* ------------------------------------------------------------------ */

	if (Puser == (int *) NULL)
	{
	    for (k = 0 ; k < n ; k++)
	    {
		P [k] = k ;
	    }
	}
	else
	{
	    for (k = 0 ; k < n ; k++)
	    {
		P [k] = Puser [k] ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* return the symbolic object */
    /* ---------------------------------------------------------------------- */

    Info [KLU_BTF_INFO_STATUS] = KLU_OK ;

    Symbolic->nblocks = nblocks ;
    Info [KLU_BTF_INFO_NBLOCKS] = nblocks ;

    Symbolic->maxblock = maxblock ;
    Info [KLU_BTF_INFO_MAXBLOCK] = maxblock ;

    Symbolic->lnz = EMPTY ;
    Symbolic->unz = EMPTY ;
    Symbolic->nzoff = nzoff ;
    Info [KLU_BTF_INFO_EST_LNZ] = EMPTY ;
    Info [KLU_BTF_INFO_NZOFF] = nzoff ;

    return (Symbolic) ;
}
