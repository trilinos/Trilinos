/* ========================================================================== */
/* === klu_btf_solve ======================================================== */
/* ========================================================================== */

/* Solve Ax=b using the symbolic and numeric objects from klu_btf_analyze and
 * klu_btf_factor.  Note that no iterative refinement is performed.
 * Uses Numeric->Xwork as workspace (undefined on input and output),
 * of size 4n double's (note that columns 2 to 4 of Xwork overlap with
 * Numeric->Iwork).
 *
 * TODO: add iterative refinement?
 * TODO: solve A'x=b
 * TODO: sparse right-hand-sides
 * TODO: add error checking of inputs (and allow it to be turned off)
 */  

#include "klu_btf_internal.h"

void klu_btf_solve
(
    /* inputs, not modified */
    klu_symbolic *Symbolic,
    klu_numeric *Numeric,
    int d,		    /* leading dimension of B */
    int nrhs,		    /* number of right-hand-sides */

    /* right-hand-side on input, overwritten with solution to Ax=b on output */
    double B [ ]	    /* size n*nrhs, in column-oriented form, with
			     * leading dimension d. */
)
{
    double x [4], offik, r ;
    double *Singleton, **Lbx, **Ubx, *Offx, *Rs, *X, *Y, *C ;
    int k1, k2, nk, k, block, pend, row, n, p, *Q, *R, nblocks, poff, *Pnum,
	*Offp, *Offi, **Lbp, **Lbi, **Ubp, **Ubi, s, chunk, nr, n2, n3, i ;

    /* ---------------------------------------------------------------------- */
    /* get the contents of the Symbolic object */
    /* ---------------------------------------------------------------------- */

    n = Symbolic->n ;
    nblocks = Symbolic->nblocks ;
    Q = Symbolic->Q ;
    R = Symbolic->R ;

    PRINTF (("klu_btf_solve:  n %d nblocks %d \n", n, nblocks)) ;

    /* ---------------------------------------------------------------------- */
    /* get the contents of the Numeric object */
    /* ---------------------------------------------------------------------- */

    ASSERT (nblocks == Numeric->nblocks) ;
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
    X = Numeric->Xwork ;
    ASSERT (klu_valid (n, Offp, Offi, Offx)) ;

    n2 = n*2 ;
    n3 = n*3 ;

    /* ---------------------------------------------------------------------- */
    /* solve in chunks of 4 columns at a time */
    /* ---------------------------------------------------------------------- */

    for (chunk = 0 ; chunk < nrhs ; chunk += 4)
    {

	/* ------------------------------------------------------------------ */
	/* get the size of the current chunk */
	/* ------------------------------------------------------------------ */

	nr = MIN (nrhs - chunk, 4) ;

	/* ------------------------------------------------------------------ */
	/* scale and permute the right hand side, X = P*(R\B) */
	/* ------------------------------------------------------------------ */

	Y = X ;
	C = B ;
	for (s = 0 ; s < nr ; s++)
	{
	    for (k = 0 ; k < n ; k++)
	    {
		Y [k] = C [Pnum [k]] / Rs [k] ;
	    }
	    Y += n ;
	    C += d ;
	}

	/* ------------------------------------------------------------------ */
	/* solve X = (L*U + Off)\X */
	/* ------------------------------------------------------------------ */

	for (block = nblocks-1 ; block >= 0 ; block--)
	{

	    /* -------------------------------------------------------------- */
	    /* the block of size nk is from rows/columns k1 to k2-1 */
	    /* -------------------------------------------------------------- */

	    k1 = R [block] ;
	    k2 = R [block+1] ;
	    nk = k2 - k1 ;
	    PRINTF (("solve %d, k1 %d k2-1 %d nk %d\n", block, k1,k2-1,nk)) ;

	    /* solve the block system */
	    if (nk == 1)
	    {
		switch (nr)
		{
		    case 4: X [k1 + n3] /= Singleton [block] ;
		    case 3: X [k1 + n2] /= Singleton [block] ;
		    case 2: X [k1 + n ] /= Singleton [block] ;
		    case 1: X [k1     ] /= Singleton [block] ;
		}
	    }
	    else
	    {
		/* solve the block system */
		klu_lsolve (nk, Lbp [block], Lbi [block], Lbx [block],
			n, nr, X + k1) ;
		klu_usolve (nk, Ubp [block], Ubi [block], Ubx [block],
			n, nr, X + k1) ;
	    }

	    /* -------------------------------------------------------------- */
	    /* block back-substitution for the off-diagonal-block entries */
	    /* -------------------------------------------------------------- */

	    if (block > 0)
	    {
		switch (nr)
		{
		case 4:

		    for (k = k1 ; k < k2 ; k++)	/* order is not relevant */
		    {
			pend = Offp [k+1] ;
			x [0] = X [k     ] ;
			x [1] = X [k + n ] ;
			x [2] = X [k + n2] ;
			x [3] = X [k + n3] ;
			for (p = Offp [k] ; p < pend ; p++)
			{
			    i = Offi [p] ;
			    offik = Offx [p] ;
			    X [i     ] -= offik * x [0] ;
			    X [i + n ] -= offik * x [1] ;
			    X [i + n2] -= offik * x [2] ;
			    X [i + n3] -= offik * x [3] ;
			}
		    }
		    break ;

		case 3:

		    for (k = k1 ; k < k2 ; k++)
		    {
			pend = Offp [k+1] ;
			x [0] = X [k     ] ;
			x [1] = X [k + n ] ;
			x [2] = X [k + n2] ;
			for (p = Offp [k] ; p < pend ; p++)
			{
			    i = Offi [p] ;
			    offik = Offx [p] ;
			    X [i     ] -= offik * x [0] ;
			    X [i + n ] -= offik * x [1] ;
			    X [i + n2] -= offik * x [2] ;
			}
		    }
		    break ;

		case 2:

		    for (k = k1 ; k < k2 ; k++)
		    {
			pend = Offp [k+1] ;
			x [0] = X [k    ] ;
			x [1] = X [k + n] ;
			for (p = Offp [k] ; p < pend ; p++)
			{
			    i = Offi [p] ;
			    offik = Offx [p] ;
			    X [i    ] -= offik * x [0] ;
			    X [i + n] -= offik * x [1] ;
			}
		    }
		    break ;

		case 1:

		    for (k = k1 ; k < k2 ; k++)
		    {
			pend = Offp [k+1] ;
			x [0] = X [k] ;
			for (p = Offp [k] ; p < pend ; p++)
			{
			    X [Offi [p]] -= Offx [p] * x [0] ;
			}
		    }
		    break ;
		}
	    }
	}

	/* ------------------------------------------------------------------ */
	/* permute the result, B = Q*X */
	/* ------------------------------------------------------------------ */

	Y = X ;
	C = B ;
	for (s = 0 ; s < nr ; s++)
	{
	    for (k = 0 ; k < n ; k++)
	    {
		C [Q [k]] = Y [k] ;
	    }
	    Y += n ;
	    C += d ;
	}

	/* ------------------------------------------------------------------ */
	/* go to the next chunk of B */
	/* ------------------------------------------------------------------ */

	B += 4*d ;
    }
}
