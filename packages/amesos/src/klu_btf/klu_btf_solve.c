/* ========================================================================== */
/* === klu_btf_solve ======================================================== */
/* ========================================================================== */

/* Solve Ax=b using the symbolic and numeric objects from klu_btf_analyze and
 * klu_btf_factor.  Note that no iterative refinement is performed.
 * Uses Numeric->X as workspace (undefined on input and output).
 *
 * TODO:
 *	* add iterative refinement?  
 *	* solve A'x=b
 *	* multiple and/or sparse right-hand-sides
 *	* add error checking of inputs
 *	* exploit sparsity better.  Note that if b is all zero, this routine
 *	    still takes O(n) time.  Total time is O (n + flop count).
 */  

#include "klu_btf_internal.h"

void klu_btf_solve
(
    /* inputs, not modified */
    klu_symbolic *Symbolic,
    klu_numeric *Numeric,

    /* right-hand-side on input, overwritten with solution to Ax=b on output */
    double B [ ]
)
{
    double *Singleton, **Lbx, **Ubx, *Offx, *Rs, xk, *X ;
    int k1, k2, nk, k, block, pend, row, n, p, *Q, *R,
	nblocks, poff, *Pnum, *Offp, *Offi, **Lbp, **Lbi, **Ubp, **Ubi ;

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
    X = Numeric->X ;
    ASSERT (klu_valid (n, Offp, Offi, Offx)) ;

    /* ---------------------------------------------------------------------- */
    /* scale and permute the right hand side, x = P*(R\b) */
    /* ---------------------------------------------------------------------- */

    for (k = 0 ; k < n ; k++)
    {
	row = Pnum [k] ;
	X [k] = B [row] / Rs [row] ;
	PRINTF (("X (%d) = %g scaled and permuted\n", k, X [k])) ;
    }

    /* ---------------------------------------------------------------------- */
    /* solve (L*U + Off) x = c, where x starts out equal to c = P*(R\b) */
    /* ---------------------------------------------------------------------- */

    for (block = nblocks-1 ; block >= 0 ; block--)
    {

	/* ------------------------------------------------------------------ */
	/* the block is from rows/columns k1 to k2-1 */
	/* ------------------------------------------------------------------ */

	k1 = R [block] ;
	k2 = R [block+1] ;
	nk = k2 - k1 ;
	PRINTF (("SOLVE BLOCK %d, k1 %d k2-1 %d nk %d\n", block, k1,k2-1,nk)) ;

	if (nk == 1)
	{

	    /* solve the singleton system */
	    xk = X [k1] / Singleton [block] ;
	    X [k1] = xk ;

	    /* back-substitution for the off-diagonal entries */
	    if (xk != 0.0)
	    {
		pend = Offp [k1+1] ;
		for (p = Offp [k1] ; p < pend ; p++)
		{
		    ASSERT (Offi [p] < k1) ;
		    X [Offi [p]] -= Offx [p] * xk ;
		}
	    }

	}
	else
	{

	    /* solve the block system */
	    ASSERT (klu_valid (nk, Lbp [block], Lbi [block], Lbx [block])) ;
	    ASSERT (klu_valid (nk, Ubp [block], Ubi [block], Ubx [block])) ;
	    klu_lsolve (nk, Lbp [block], Lbi [block], Lbx [block], X + k1) ;
#ifndef NDEBUG
	    for (k = k1 ; k < k2 ; k++) PRINTF (("Lsol X[%d] = %g\n",k,X[k])) ;
#endif
	    klu_usolve (nk, Ubp [block], Ubi [block], Ubx [block], X + k1) ;
#ifndef NDEBUG
	    for (k = k1 ; k < k2 ; k++) PRINTF (("Usol X[%d] = %g\n",k,X[k])) ;
#endif

	    /* block back-substitution for the off-diagonal entries */
	    if (block > 0)
	    {
		for (k = k1 ; k < k2 ; k++) /* order is not relevant */
		{
		    xk = X [k] ;
		    if (xk != 0.0)
		    {
			pend = Offp [k+1] ;
			for (p = Offp [k] ; p < pend ; p++)
			{
			    ASSERT (Offi [p] < k1) ;
			    X [Offi [p]] -= Offx [p] * xk ;
			}
		    }
		}
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* permute the result, B = Q*X */
    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    for (k = 0 ; k < n ; k++) PRINTF (("X[%d] = %g before Q\n",k,X[k])) ;
    for (k = 0 ; k < n ; k++) PRINTF (("Q[%d] = %d\n",k,Q[k])) ;
#endif

    for (k = 0 ; k < n ; k++)
    {
	B [Q [k]] = X [k] ;
    }

#ifndef NDEBUG
    for (k = 0 ; k < n ; k++) PRINTF (("Soln[%d] = %g after Q\n",k,B[k])) ;
#endif
}
