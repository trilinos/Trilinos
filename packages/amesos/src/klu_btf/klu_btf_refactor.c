/* ========================================================================== */
/* === klu_btf_refactor ===================================================== */
/* ========================================================================== */

/* Factor the matrix, after ordering and analyzing it with klu_btf_analyze,
 * and factoring it once with klu_btf_factor.  This routine cannot do any
 * numerical pivoting.  The pattern of the input matrix (Ap, Ai) must be
 * identical to the pattern given to klu_btf_factor.
 *
 * TODO: compute Numeric->umin and Numeric->umax in refactor!
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
    klu_control *user_control,

    /* input, values (Rs, Lx, Ux, Offx, umin, umax) modified on output */
    klu_numeric *Numeric
)
{
    double ukk, ujk, umin, umax, s ;
#ifndef NRECIPROCAL
    double ukk_inverse ; 
#endif
    double *Lnz, *Singleton, **Lbx, **Ubx, *Offx, *Rs, *Lx, *Ux, *X ;
    int k1, k2, nk, k, block, oldcol, pend, oldrow, n, p, newrow, scale,
	nblocks, poff, result, i, j, up, upend, upstart, maxblock ;
    int *P, *Q, *R, *Pnum, *Lp, *Up, *Offp, *Offi, *Ui, *Li, *Pinv ;
    int **Lbp, **Lbi, **Ubp, **Ubi ;
    klu_control *control, default_control ;

    /* ---------------------------------------------------------------------- */
    /* get the contents of the Symbolic object */
    /* ---------------------------------------------------------------------- */

    n = Symbolic->n ;
    P = Symbolic->P ;
    Q = Symbolic->Q ;
    R = Symbolic->R ;
    Lnz = Symbolic->Lnz ;
    nblocks = Symbolic->nblocks ;
    maxblock = Symbolic->maxblock ;

    /* ---------------------------------------------------------------------- */
    /* get the contents of the Numeric object */
    /* ---------------------------------------------------------------------- */

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
    X = Numeric->Xwork ;
    Numeric->nlrealloc = 0 ;
    Numeric->nurealloc = 0 ;

    /* ---------------------------------------------------------------------- */
    /* get control parameters and make sure they are in the proper range */
    /* ---------------------------------------------------------------------- */

    if (user_control == (klu_control *) NULL)
    {
	control = &default_control ;
	klu_btf_defaults (control) ;
    }
    else
    {
	control = user_control ;
    }

    scale = control->scale ;
    scale = MAX (0, scale) ;
    scale = MIN (2, scale) ;
    Numeric->scale = scale ;

    PRINTF (("klu_btf_factor scale %d:  n %d nzoff %d nblocks %d maxblock %d\n",
	scale, n, Symbolic->nzoff, nblocks, maxblock)) ;

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
		/* note that this will cause a failure for singular matrices */
		Singleton [block] = 100 * block ;
	    }
	    else
	    {
		int *Lp, *Up ;
		double *Ux, *Lx ;
		Lp = Lbp [block] ;
		Lx = Lbx [block] ;
		for (p = 0 ; p < Lp [nk] ; p++) Lx [p] = 42 + p * 1e-5 ;
		Up = Ubp [block] ;
		Ux = Ubx [block] ;
		for (p = 0 ; p < Up [nk] ; p++) Ux [p] = 99 + p * 1e-4 ;
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
    result = klu_btf_scale (scale, n, Ap, Ai, Ax, Rs, (int *) NULL) ;
    if (result != KLU_OK)
    {
	return (result) ;
    }

    poff = 0 ;

    /* clear X for use in lufact */
    for (k = 0 ; k < maxblock ; k++)
    {
	X [k] = 0 ;
    }

    umin = 0 ;
    umax = 0 ;

    /* ---------------------------------------------------------------------- */
    /* factor each block */
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
	    s = 0 ;

	    if (scale == 0)
	    {
		/* no scaling */
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
			Offx [poff++] = Ax [p] ;
		    }
		    else
		    {
			if (newrow != k1)
			{
			    return (KLU_INVALID) ;
			}
			PRINTF (("Singleton block %d %g\n", block, Ax [p])) ;
			s = Ax [p] ;
		    }
		}
	    }
	    else
	    {
		/* row scaling */
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
#ifndef NRECIPROCAL
			Offx [poff++] = Ax [p] * Rs [oldrow] ;
#else
			Offx [poff++] = Ax [p] / Rs [oldrow] ;
#endif
		    }
		    else
		    {
			if (newrow != k1)
			{
			    return (KLU_INVALID) ;
			}
			PRINTF (("Singleton block %d %g\n", block, Ax [p])) ;
#ifndef NRECIPROCAL
			s = Ax [p] * Rs [oldrow] ;
#else
			s = Ax [p] / Rs [oldrow] ;
#endif
		    }
		}
	    }

#ifndef NRECIPROCAL
	    Singleton [block] = 1.0 / s ;
#else
	    Singleton [block] = s ;
#endif

	    /* keep track of the smallest and largest diagonal entry */
	    ukk = fabs (s) ;
	    if (block == 0)
	    {
		umin = ukk ;
		umax = ukk ;
	    }
	    else if (ukk < umin)
	    {
		umin = ukk ;
	    }
	    else if (ukk > umax)
	    {
		umax = ukk ;
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

		if (scale == 0)
		{
		    /* no scaling */
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
			    Offx [poff++] = Ax [p] ;
			}
			else
			{
			    /* (newrow,k) is an entry in the block */
			    ASSERT (newrow < k2) ;
			    newrow -= k1 ;
			    if (newrow >= nk)
			    {
				return (KLU_INVALID) ;
			    }
			    X [newrow] = Ax [p] ;
			}
		    }
		}
		else
		{
		    /* row scaling */
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
#ifndef NRECIPROCAL
			    Offx [poff++] = Ax [p] * Rs [oldrow] ;
#else
			    Offx [poff++] = Ax [p] / Rs [oldrow] ;
#endif
			}
			else
			{
			    /* (newrow,k) is an entry in the block */
			    ASSERT (newrow < k2) ;
			    newrow -= k1 ;
			    if (newrow >= nk)
			    {
				return (KLU_INVALID) ;
			    }
#ifndef NRECIPROCAL
			    X [newrow] = Ax [p] * Rs [oldrow] ;
#else
			    X [newrow] = Ax [p] / Rs [oldrow] ;
#endif
			}
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
#ifndef NRECIPROCAL
		ukk_inverse = 1.0 / ukk ;
		Ux [upend] = ukk_inverse ;
#else
		Ux [upend] = ukk ;
#endif

		/* create the unit diagonal of L */
		p = Lp [k] ;
		Lx [p++] = 1 ;

		/* gather and scale to get the kth column of L */
		pend = Lp [k+1] ;
		if (ukk == 0.0)
		{
		    for ( ; p < pend ; p++)
		    {
			i = Li [p] ;
			if (Lx [p] != 0)
			{
			    Lx [p] = X [i] / ukk ;
			}
			X [i] = 0 ;
		    }
		}
		else
		{
		    for ( ; p < pend ; p++)
		    {
			i = Li [p] ;
#ifndef NRECIPROCAL
			Lx [p] = X [i] * ukk_inverse ;
#else
			Lx [p] = X [i] / ukk ;
#endif
			X [i] = 0 ;
		    }
		}

#ifndef NDEBUG
		for (i = 0 ; i < nk ; i++) ASSERT (X [i] == 0) ;
#endif

		/* keep track of the smallest and largest diagonal entry */
		ukk = fabs (ukk) ;
		if (k == 0 && block == 0)
		{
		    umin = ukk ;
		    umax = ukk ;
		}
		else if (ukk < umin)
		{
		    umin = ukk ;
		}
		else if (ukk > umax)
		{
		    umax = ukk ;
		}

	    }

	    PRINTF (("lufact done\n")) ;
	    PRINTF (("\n----------------------- L %d:\n", block)) ;
	    ASSERT (klu_valid (nk, Lbp [block], Lbi [block], Lbx [block])) ;
	    PRINTF (("\n----------------------- U %d:\n", block)) ;
	    ASSERT (klu_valid (nk, Ubp [block], Ubi [block], Ubx [block])) ;

	}
    }

    Numeric->umin = umin ;
    Numeric->umax = umax ;

    /* ---------------------------------------------------------------------- */
    /* permute scale factors Rs according to pivotal row order */
    /* ---------------------------------------------------------------------- */

    if (scale != 0)
    {
	for (k = 0 ; k < n ; k++)
	{
	    X [k] = Rs [Pnum [k]] ;
	}
	for (k = 0 ; k < n ; k++)
	{
	    Rs [k] = X [k] ;
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
