/* ========================================================================== */
/* === klu_btf_refactor ===================================================== */
/* ========================================================================== */

/* Factor the matrix, after ordering and analyzing it with klu_btf_analyze,
 * and factoring it once with klu_btf_factor.  This routine cannot do any
 * numerical pivoting.
 */

#include "klu_btf.h"
#include "klu_kernel.h"
#include "klu_dump.h"

/* ========================================================================== */

/* Sparse left-looking LU factorization, with no pivoting.  Assumes pattern of
 * L and U are already computed, and in topologically sorted order. */

/* TODO: NaN handling */

static int lufact
(
    /* input, not modified */
    int n,	    /* A is n-by-n */
    int Ap [ ],	    /* size n+1, column pointers for A */
    int Ai [ ],	    /* size nz = Ap [n], row indices for A */
    double Ax [ ],  /* size nz, values of A */
    int Lp [ ],	    /* size n+1 */
    int Li [ ],	    /* size lsize = Lp [n] */
    int Up [ ],	    /* size n+1 */
    int Ui [ ],	    /* size usize = Up [n] */

    /* output */
    double Lx [ ],  /* size lsize */
    double Ux [ ],  /* size usize */

    /* workspace, not defined on input or output */
    double X [ ]   /* size n */
)
{
    int i, j, k, p, pend, up, upend ;
    double ukk, ujk ;

    for (k = 0 ; k < n ; k++)
    {
	X [k] = 0 ;
    }

    for (k = 0 ; k < n ; k++)
    {

	/* scatter kth column of A into workspace X */
	pend = Ap [k+1] ;
	for (p = Ap [k] ; p < pend ; p++)
	{
	    X [Ai [p]] = Ax [p] ;
	}

	/* compute the kth column of U, and update the kth column of A */
	upend = Up [k+1] - 1 ;
	for (up = Up [k] ; up < upend ; up++)	/* skip U_kk */
	{
	    j = Ui [up] ;
	    ujk = X [j] ;
	    X [j] = 0 ;
	    Ux [up] = ujk ;
	    pend = Lp [j+1] ;
	    for (p = Lp [j] + 1 ; p < pend ; p++)
	    {
		X [Li [p]] -= Lx [p] * ujk ;
	    }
	}

	/* get the diagonal entry of U */
	ukk = X [k] ;
	X [k] = 0 ;
	Ux [upend] = ukk ;

	/* singular case */
	if (ukk == 0.0)
	{
	    return (KLU_SINGULAR) ;
	}

	/* create the unit diagonal of L */
	p = Lp [k] ;
	Lx [p++] = 1 ;

	/* gather and scale to get the kth column of L */
	pend = Lp [k+1] ;
	for ( ; p < pend ; p++)
	{
	    i = Li [p] ;
	    Lx [p] = X [i] / ukk ;
	    X [i] = 0 ;
	}

    }

    return (KLU_OK) ;
}
/* ========================================================================== */

static int klu_btf_refactor2
(
    /* inputs, not modified */
    int Ap [ ],		/* size n+1, column pointers */
    int Ai [ ],		/* size nz, row indices */
    double Ax [ ],
    klu_symbolic *Symbolic,

    /* defined on input, numerical values modified on output */
    klu_numeric *Numeric,

    /* workspace */
    int Cp [ ],
    int Ci [ ],
    double Cx [ ],
    int Pinv [ ],
    double W [ ]
)
{
    int k1, k2, nk, k, block, oldcol, pend, row, oldrow, newcol, n,
	result, pc, p, newrow, col, *P, *Q, *R, nblocks, poff, *Pnum, *Lp, *Up,
	*Offp, *Offi, **Lbp, **Lbi, **Ubp, **Ubi, *Pblock ;
    double Control [KLU_CONTROL], *Lnz, *Singleton, **Lbx, **Ubx, *Offx, *Rs ;

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

    /* get the contents of the Numeric object */
    Pnum = Numeric->Pnum ;
    Offp = Numeric->Offp ;
    Offi = Numeric->Offi ;
    Offx = Numeric->Offx ;
    Singleton = Numeric->Singleton ;
    Lbp = Numeric->Lbp ;
    Lbi = Numeric->Lbi ;
    Lbx = Numeric->Lbx ;
    Ubp = Numeric->Ubp ;
    Ubi = Numeric->Ubi ;
    Ubx = Numeric->Ubx ;
    Rs = Numeric->Rs ;

    /* compute the inverse of Pnum (TODO could be done in klu_btf_analyze) */
#ifndef NDEBUG
    for (k = 0 ; k < n ; k++)
    {
	Pinv [k] = EMPTY ;
    }
#endif
    for (k = 0 ; k < n ; k++)
    {
	ASSERT (Pnum [k] >= 0 && Pnum [k] < n) ;
	Pinv [Pnum [k]] = k ;
    }
#ifndef NDEBUG
    for (k = 0 ; k < n ; k++) ASSERT (Pinv [k] != EMPTY) ;
#endif

    for (block = 0 ; block < nblocks ; block++)
    {
	Singleton [block] = 0 ;
    }
    poff = 0 ;

    /* ---------------------------------------------------------------------- */
    /* for testing only */
    /* ---------------------------------------------------------------------- */

#ifdef TESTING
    /* randomnly mangle the numerical values to ensure we really recompute the
     * factors.  This is for testing only! */
    for (p = 0 ; p < Symbolic->nzoff ; p++) Offx [p] = p ;
    for (k = 0 ; k < n ; k++) Rs [k] = 1 + k/1000. ;
#endif

    /* ---------------------------------------------------------------------- */
    /* compute the row scale factors, Rs */
    /* ---------------------------------------------------------------------- */

    if (!klu_btf_scale (n, Ap, Ai, Ax, Rs))
    {
	/* matrix is singular */
	PRINTF (("Row scale: matrix is singular\n")) ;
	return (FALSE) ;
    }

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

	    ASSERT (Offp [k1] == poff) ;
	    oldcol = Q [k1] ;
	    pend = Ap [oldcol+1] ;
	    for (p = Ap [oldcol] ; p < pend ; p++)
	    {
		oldrow = Ai [p] ;
		newrow = Pinv [oldrow] ;
		if (newrow < k1)
		{
		    ASSERT (poff < Symbolic->nzoff) ; 
		    ASSERT (Offi [poff] == newrow) ;
		    Offx [poff++] = Ax [p] / Rs [oldrow] ;
		}
		else
		{
		    ASSERT (newrow == k1) ;
		    PRINTF (("Singleton block %d value %g\n", block, Ax [p])) ;
		    Singleton [block] = Ax [p] / Rs [oldrow] ;
		}
	    }
	    if (Singleton [block] == 0.0)
	    {
		PRINTF (("Matrix is singular (zero 1-by-1 block)\n")) ;
		return (FALSE) ;
	    }

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* construct the kth block, C */
	    /* -------------------------------------------------------------- */

	    pc = 0 ;
	    for (k = k1 ; k < k2 ; k++)
	    {
		newcol = k-k1 ;
		Cp [newcol] = pc ;
		ASSERT (Offp [k] == poff) ;
		oldcol = Q [k] ;
		pend = Ap [oldcol+1] ;
		for (p = Ap [oldcol] ; p < pend ; p++)
		{
		    oldrow = Ai [p] ;
		    newrow = Pinv [oldrow] ;
		    if (newrow < k1)
		    {
			/* this is an entry in the off-diagonal part */
			ASSERT (poff < Symbolic->nzoff) ; 
			ASSERT (Offi [poff] == newrow) ;
			Offx [poff++] = Ax [p] / Rs [oldrow] ;
		    }
		    else
		    {
			/* (newrow,newcol) is an entry in the block */
			ASSERT (newrow < k2) ;
			ASSERT (pc < Symbolic->maxnz) ;
			newrow -= k1 ;
			Ci [pc] = newrow ;
			Cx [pc] = Ax [p] / Rs [oldrow] ;
			pc++ ;
		    }
		}
	    }
	    ASSERT (pc <= Symbolic->maxnz) ;
	    Cp [nk] = pc ;
	    PRINTF (("\n----------------------- block %d, C:\n", block)) ;
	    ASSERT (klu_valid (nk, Cp, Ci, Cx)) ;

	    /* -------------------------------------------------------------- */
	    /* factor the block */
	    /* -------------------------------------------------------------- */

	    PRINTF (("calling lufact\n")) ;
	    result = lufact (nk, Cp, Ci, Cx,
		    Lbp [block], Lbi [block],
		    Ubp [block], Ubi [block],
		    Lbx [block], Ubx [block], W) ;
	    PRINTF (("lufact done\n")) ;
	    if (result != KLU_OK)
	    {
		return (result) ;
	    }
	    PRINTF (("\n----------------------- L %d:\n", block)) ;
	    ASSERT (klu_valid (nk, Lbp [block], Lbi [block], Lbx [block])) ;
	    PRINTF (("\n----------------------- U %d:\n", block)) ;
	    ASSERT (klu_valid (nk, Ubp [block], Ubi [block], Ubx [block])) ;

	}
    }

    ASSERT (Offp [n] == poff) ;
    ASSERT (Symbolic->nzoff == poff) ;
    PRINTF (("\n------------------- Off diagonal entries, new:\n")) ;
    ASSERT (klu_valid (n, Offp, Offi, Offx)) ;

    return (KLU_OK) ;
}

/* ========================================================================== */
/* === klu_btf_refactor ===================================================== */
/* ========================================================================== */

int klu_btf_refactor	/* returns KLU_OK if OK, < 0 if error */
(
    /* inputs, not modified */
    int Ap [ ],		/* size n+1, column pointers */
    int Ai [ ],		/* size nz, row indices */
    double Ax [ ],
    klu_symbolic *Symbolic,
    klu_numeric *Numeric
)
{
    int n, nzoff, nblocks, maxblock, maxnz, *Cp, *Ci, *Pinv, result ;
    double *Cx, *W ;

    /* ---------------------------------------------------------------------- */
    /* get the contents of the Symbolic object */
    /* ---------------------------------------------------------------------- */

    n = Symbolic->n ;
    nzoff = Symbolic->nzoff ;
    nblocks = Symbolic->nblocks ;
    maxblock = Symbolic->maxblock ;
    maxnz = Symbolic->maxnz ;
    PRINTF (("klu_btf_factor:  n %d nzoff %d nblocks %d maxblock %d maxnz %d\n",
	n, nzoff, nblocks, maxblock, maxnz)) ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    Cp = (int *) ALLOCATE ((maxblock + 1) * sizeof (int)) ;
    Ci = (int *) ALLOCATE ((maxnz + 1) * sizeof (int)) ;
    Cx = (double *) ALLOCATE ((maxnz + 1) * sizeof (double)) ;
    Pinv = (int *) ALLOCATE (n * sizeof (int)) ;
    W = (double *) ALLOCATE (n * sizeof (double)) ;

    if ((Cp == (int *) NULL) || (Ci == (int *) NULL) ||
	(Cx == (double *) NULL) || (Pinv == (int *) NULL) ||
	(W == (double *) NULL))
    {
	FREE (Cp, int) ;
	FREE (Ci, int) ;
	FREE (Cx, double) ;
	FREE (Pinv, int) ;
	FREE (W, double) ;
	return (KLU_OUT_OF_MEMORY) ;
    }

#ifdef TESTING
    /* randomly mangle the numerical values to test the refactor code */
    {
	int block, k1, k2, nk, p ;
	for (block = 0 ; block < nblocks ; block++)
	{
	    k1 = Symbolic->R [block] ;
	    k2 = Symbolic->R [block+1] ;
	    nk = k2 - k1 ;
	    if (nk == 1)
	    {
		Numeric->Singleton [block] = 100 * block ;
	    }
	    else
	    {
		int *Lp, *Up ;
		double *Ux, *Lx ;
		Lp = Numeric->Lbp [block] ;
		Lx = Numeric->Lbx [block] ;
		for (p = 0 ; p < Lp [nk] ; p++) Lx [p] = 42 + p / 10000. ;
		Up = Numeric->Ubp [block] ;
		Ux = Numeric->Ubx [block] ;
		for (p = 0 ; p < Up [nk] ; p++) Ux [p] = 99 + p / 1000. ;
	    }
	}
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* factorize the blocks */
    /* ---------------------------------------------------------------------- */

    PRINTF (("calling klu_btf_refactor2\n")) ;
    result = klu_btf_refactor2 (Ap, Ai, Ax, Symbolic, Numeric,
	    Cp, Ci, Cx, Pinv, W) ;
    PRINTF (("klu_btf_refactor2 done\n")) ;

    /* ---------------------------------------------------------------------- */
    /* free workspace */
    /* ---------------------------------------------------------------------- */

    FREE (Cp, int) ;
    FREE (Ci, int) ;
    FREE (Cx, double) ;
    FREE (Pinv, int) ;
    FREE (W, double) ;

    /* ---------------------------------------------------------------------- */
    /* return the error code */
    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    {
	int block, k1, k2, nk ;
	PRINTF (("\n ########### KLU_BTF_REFACTOR done, nblocks %d\n",nblocks));
	for (block = 0 ; block < nblocks ; block++)
	{
	    k1 = Symbolic->R [block] ;
	    k2 = Symbolic->R [block+1] ;
	    nk = k2 - k1 ;
	    PRINTF (("\n======================klu_btf_refactor output: k1 %d k2 %d nk %d\n",k1,k2,nk)) ;
	    if (nk == 1)
	    {
		PRINTF (("singleton %g\n", Numeric->Singleton [block])) ;
	    }
	    else
	    {
		int *Lp, *Li, *Up, *Ui ;
		double *Ux, *Lx ;
		Lp = Numeric->Lbp [block] ;
		Li = Numeric->Lbi [block] ;
		Lx = Numeric->Lbx [block] ;
		PRINTF (("\n---- L block %d\n", block)); 
		ASSERT (klu_valid (nk, Lp, Li, Lx)) ;
		Up = Numeric->Ubp [block] ;
		Ui = Numeric->Ubi [block] ;
		Ux = Numeric->Ubx [block] ;
		PRINTF (("\n---- U block %d\n", block)) ; 
		ASSERT (klu_valid (nk, Up, Ui, Ux)) ;
	    }
	}
    }
#endif

    return (result) ;
}
