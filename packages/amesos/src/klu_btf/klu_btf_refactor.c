/* ========================================================================== */
/* === klu_btf_refactor ===================================================== */
/* ========================================================================== */

/* Factor the matrix, after ordering and analyzing it with klu_btf_analyze,
 * and factoring it once with klu_btf_factor.  This routine cannot do any
 * numerical pivoting.  The pattern of the input matrix (Ap, Ai) must be
 * identical to the pattern given to klu_btf_factor.
 */

#include "klu_btf_internal.h"

/* ========================================================================== */

int klu_btf_refactor	/* returns KLU_OK or KLU_INVALID */
(
    /* inputs, not modified */
    int Ap [ ],		/* size n+1, column pointers */
    int Ai [ ],		/* size nz, row indices */
    double Ax [ ],
    klu_symbolic *Symbolic,

    /* defined on input, numerical values modified on output */
    klu_numeric *Numeric
)
{
    double Control [KLU_CONTROL], *Lnz, *Singleton, **Lbx, **Ubx, *Offx, *Rs,
	ukk, ujk, *Lx, *Ux, *X ;

    int k1, k2, nk, k, block, oldcol, pend, row, oldrow, newcol, n,
	pc, p, newrow, col, *P, *Q, *R, nblocks, poff, *Pnum, *Lp, *Up,
	*Offp, *Offi, **Lbp, **Lbi, **Ubp, **Ubi, *Pblock, maxnz, result,
	i, j, up, upend, upstart, *Ui, *Li, *Pinv, maxblock ;

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
    maxnz = Symbolic->maxnz ;
    maxblock = Symbolic->maxblock ;

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
    Pinv = Numeric->Pinv ;
    X = Numeric->X ;

    PRINTF (("klu_btf_factor:  n %d nzoff %d nblocks %d maxblock %d maxnz %d\n",
	n, Symbolic->nzoff, nblocks, maxblock, maxnz)) ;

#ifdef TESTING
    /* randomly mangle the numerical values to test the refactor code */
    {
	int block, k1, k2, nk, p ;
	for (block = 0 ; block < nblocks ; block++)
	{
	    k1 = R [block] ;
	    k2 = R [block+1] ;
	    nk = k2 - k1 ;
	    if (nk == 1)
	    {
		Singleton [block] = 100 * block ;
	    }
	    else
	    {
		int *Lp, *Up ;
		double *Ux, *Lx ;
		Lp = Lbp [block] ;
		Lx = Lbx [block] ;
		for (p = 0 ; p < Lp [nk] ; p++) Lx [p] = 42 + p / 10000. ;
		Up = Ubp [block] ;
		Ux = Ubx [block] ;
		for (p = 0 ; p < Up [nk] ; p++) Ux [p] = 99 + p / 1000. ;
	    }
	}
	for (p = 0 ; p < Symbolic->nzoff ; p++) Offx [p] = p ;
	for (k = 0 ; k < n ; k++) Rs [k] = 1 + k/1000. ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* check the input matrix compute the row scale factors, Rs */
    /* ---------------------------------------------------------------------- */

    /* check for out-of-range indices, but do not check for duplicates */
    result = klu_btf_scale (n, Ap, Ai, Ax, Rs, (int *) NULL) ;
    if (result != KLU_OK)
    {
	return (result) ;
    }

    /* TODO: no need to set Singleton to zero */
    for (block = 0 ; block < nblocks ; block++)
    {
	Singleton [block] = 0 ;
    }
    poff = 0 ;

    /* note that klu_btf_refactor only uses X [0..maxblock-1] */
    for (k = 0 ; k < maxblock ; k++)
    {
	X [k] = 0 ;
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
	ASSERT (nk <= maxblock) ;

	if (nk == 1)
	{

	    /* -------------------------------------------------------------- */
	    /* singleton case */
	    /* -------------------------------------------------------------- */

	    if (Offp [k1] != poff)
	    {
		return (KLU_INVALID) ;
	    }
	    oldcol = Q [k1] ;
	    pend = Ap [oldcol+1] ;
	    for (p = Ap [oldcol] ; p < pend ; p++)
	    {
		oldrow = Ai [p] ;
		newrow = Pinv [oldrow] ;
		if (newrow < k1)
		{
		    ASSERT (poff < Symbolic->nzoff) ; 
		    if (Offi [poff] != newrow)
		    {
			return (KLU_INVALID) ;
		    }
		    Offx [poff++] = Ax [p] / Rs [oldrow] ;
		}
		else
		{
		    if (newrow != k1)
		    {
			return (KLU_INVALID) ;
		    }
		    PRINTF (("Singleton block %d value %g\n", block, Ax [p])) ;
		    Singleton [block] = Ax [p] / Rs [oldrow] ;
		}
	    }

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* construct and factor the kth block */
	    /* -------------------------------------------------------------- */

	    PRINTF (("\n----------------------- block %d\n", block)) ;
	    Lp = Lbp [block] ;
	    Li = Lbi [block] ;
	    Up = Ubp [block] ;
	    Ui = Ubi [block] ;
	    Ux = Ubx [block] ;
	    Lx = Lbx [block] ;

	    for (k = 0 ; k < nk ; k++)
	    {

		/* scatter kth column of the block into workspace X */
		if (Offp [k+k1] != poff)
		{
		    return (KLU_INVALID) ;
		}
		oldcol = Q [k+k1] ;
		pend = Ap [oldcol+1] ;
		for (p = Ap [oldcol] ; p < pend ; p++)
		{
		    oldrow = Ai [p] ;
		    newrow = Pinv [oldrow] ;
		    if (newrow < k1)
		    {
			/* this is an entry in the off-diagonal part */
			ASSERT (poff < Symbolic->nzoff) ; 
			if (Offi [poff] != newrow)
			{
			    return (KLU_INVALID) ;
			}
			Offx [poff++] = Ax [p] / Rs [oldrow] ;
		    }
		    else
		    {
			/* (newrow,k) is an entry in the block */
			ASSERT (newrow < k2) ;
			newrow -= k1 ;
			if (newrow >= nk || pc >= maxnz)
			{
			    return (KLU_INVALID) ;
			}
			X [newrow] = Ax [p] / Rs [oldrow] ;
			/*
			Ci [pc] = newrow ;
			Cx [pc] = Ax [p] / Rs [oldrow] ;
			pc++ ;
			*/
		    }
		}

		/* compute the kth column of U, and update the kth column of A.
		 * Do in reverse order, since the column of U is stored in
		 * reverse topological order. */
		upend = Up [k+1] - 1 ;
		upstart = Up [k] ;
		for (up = upend-1 ; up >= upstart ; up--)	/* skip U_kk */
		{
		    j = Ui [up] ;
		    ASSERT (j >= 0 && j < nk) ;
		    ujk = X [j] ;
		    X [j] = 0 ;
		    Ux [up] = ujk ;
		    pend = Lp [j+1] ;
		    for (p = Lp [j] + 1 ; p < pend ; p++)
		    {
			ASSERT (Li [p] >= 0 && Li [p] < nk) ;
			X [Li [p]] -= Lx [p] * ujk ;
		    }
		}

		/* get the diagonal entry of U */
		ukk = X [k] ;
		X [k] = 0 ;
		Ux [upend] = ukk ;

		/* create the unit diagonal of L */
		p = Lp [k] ;
		Lx [p++] = 1 ;

		/* gather and scale to get the kth column of L */
		pend = Lp [k+1] ;
		for ( ; p < pend ; p++)
		{
		    i = Li [p] ;
		    ASSERT (i >= 0 && i < nk) ;
		    Lx [p] = X [i] / ukk ;
		    X [i] = 0 ;
		}
#ifndef NDEBUG
		for (i = 0 ; i < nk ; i++) ASSERT (X [i] == 0) ;
#endif
	    }

	    PRINTF (("lufact done\n")) ;
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

#ifndef NDEBUG
    if (result == KLU_OK)
    {
	int block, k1, k2, nk ;
	PRINTF (("\n ########### KLU_BTF_REFACTOR done, nblocks %d\n",nblocks));
	for (block = 0 ; block < nblocks ; block++)
	{
	    k1 = R [block] ;
	    k2 = R [block+1] ;
	    nk = k2 - k1 ;
	    PRINTF (("\n======================klu_btf_refactor output: k1 %d k2 %d nk %d\n",k1,k2,nk)) ;
	    if (nk == 1)
	    {
		PRINTF (("singleton %g\n", Singleton [block])) ;
	    }
	    else
	    {
		int *Lp, *Li, *Up, *Ui ;
		double *Ux, *Lx ;
		Lp = Lbp [block] ;
		Li = Lbi [block] ;
		Lx = Lbx [block] ;
		PRINTF (("\n---- L block %d\n", block)); 
		ASSERT (klu_valid (nk, Lp, Li, Lx)) ;
		Up = Ubp [block] ;
		Ui = Ubi [block] ;
		Ux = Ubx [block] ;
		PRINTF (("\n---- U block %d\n", block)) ; 
		ASSERT (klu_valid (nk, Up, Ui, Ux)) ;
	    }
	}
    }
#endif

    return (KLU_OK) ;
}
