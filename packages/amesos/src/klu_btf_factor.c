/* ========================================================================== */
/* === klu_btf_factor ======================================================= */
/* ========================================================================== */

/* Factor the matrix, after ordering and analyzing it with klu_btf_analyze
 *
 * WARNING: this is NOT a general purpose solver.  It will not provide good
 * orderings for matrices with very unsymmetric matrices.  For those matrices,
 * COLAMD must be used on the diagonal blocks.
 *
 * TODO:  
 *	* provide a user switch to control scaling (none, row-sum, max row).
 *	* Error checking of inputs.
 *	* test error cases / make sure there are no memory leaks
 *	* merge adjacent 1-by-1 blocks into a single upper triangular block,
 *	    for faster forward/backsolves
 *	* provide other orderings (natural, COLAMD, nested dissection, ...)
 *	* allow BTF to be turned off
 *	* return statistics, including error codes
 *	* allow for refactorization using existing LU pattern
 *	* handle small or dense diagonal blocks using dense LU
 */

#include "klu_btf.h"
#include "klu_kernel.h"
#include "klu_dump.h"

/* ========================================================================== */

static int klu_btf_factor2
(
    /* inputs, not modified */
    int Ap [ ],		/* size n+1, column pointers */
    int Ai [ ],		/* size nz, row indices */
    double Ax [ ],
    double tol,
    klu_symbolic *Symbolic,

    /* outputs, not defined on input */
    klu_numeric *Numeric,

    /* workspace */
    int Cp [ ],
    int Ci [ ],
    double Cx [ ],
    int Pinv [ ]
)
{
    int k1, k2, nk, k, block, oldcol, pend, row, oldrow, newcol, n, lnz, unz,
	result, pc, p, newrow, col, *P, *Q, *R, nblocks, poff, *Pnum, *Lp, *Up,
	*Offp, *Offi, **Lbp, **Lbi, **Ubp, **Ubi, *Pblock, noffdiag2, noffdiag ;
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

    /* compute the inverse of P (TODO could be done in klu_btf_analyze) */
#ifndef NDEBUG
    for (k = 0 ; k < n ; k++)
    {
	Pinv [k] = EMPTY ;
    }
#endif
    for (k = 0 ; k < n ; k++)
    {
	ASSERT (P [k] >= 0 && P [k] < n) ;
	Pinv [P [k]] = k ;
    }
#ifndef NDEBUG
    for (k = 0 ; k < n ; k++) ASSERT (Pinv [k] != EMPTY) ;
#endif

    for (block = 0 ; block < nblocks ; block++)
    {
	Singleton [block] = 0 ;
    }
    klu_defaults (Control) ;
    Control [KLU_TOL] = tol ;
    poff = 0 ;
    lnz = 0 ;
    unz = 0 ;
    noffdiag = 0 ;

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

	    Offp [k1] = poff ;
	    oldcol = Q [k1] ;
	    pend = Ap [oldcol+1] ;
	    for (p = Ap [oldcol] ; p < pend ; p++)
	    {
		oldrow = Ai [p] ;
		newrow = Pinv [oldrow] ;
		if (newrow < k1)
		{
		    Offi [poff] = oldrow ;
		    Offx [poff] = Ax [p] / Rs [oldrow] ;
		    poff++ ;
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
	    Pnum [k1] = P [k1] ;
	    lnz++ ;
	    unz++ ;

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
		Offp [k] = poff ;
		oldcol = Q [k] ;
		pend = Ap [oldcol+1] ;
		for (p = Ap [oldcol] ; p < pend ; p++)
		{
		    oldrow = Ai [p] ;
		    newrow = Pinv [oldrow] ;
		    if (newrow < k1)
		    {
			/* this is an entry in the off-diagonal part */
			Offi [poff] = oldrow ;
			Offx [poff] = Ax [p] / Rs [oldrow] ;
			poff++ ;
		    }
		    else
		    {
			/* (newrow,newcol) is an entry in the block */
			ASSERT (newrow < k2) ;
			newrow -= k1 ;
			Ci [pc] = newrow ;
			Cx [pc] = Ax [p] / Rs [oldrow] ;
			pc++ ;
		    }
		}
	    }
	    Cp [nk] = pc ;
	    PRINTF (("\n----------------------- block %d, C:\n", block)) ;
	    ASSERT (klu_valid (nk, Cp, Ci, Cx)) ;

	    /* -------------------------------------------------------------- */
	    /* factor the block */
	    /* -------------------------------------------------------------- */

	    PRINTF (("calling klu\n")) ;
	    Control [KLU_LSIZE] = 1.2 * Lnz [block] + nk ;
	    Control [KLU_USIZE] = Control [KLU_LSIZE] ;
	    result = klu (nk, Cp, Ci, Cx, (int *) NULL, Control,
		    &Lbp [block], &Lbi [block], &Lbx [block],
		    &Ubp [block], &Ubi [block], &Ubx [block], &Pblock,
		    &noffdiag2) ;
	    PRINTF (("klu done\n")) ;
	    if (result != KLU_OK)
	    {
		return (FALSE) ;
	    }
	    PRINTF (("\n----------------------- L %d:\n", block)) ;
	    ASSERT (klu_valid (nk, Lbp [block], Lbi [block], Lbx [block])) ;
	    PRINTF (("\n----------------------- U %d:\n", block)) ;
	    ASSERT (klu_valid (nk, Ubp [block], Ubi [block], Ubx [block])) ;

	    /* -------------------------------------------------------------- */
	    /* get statistics */
	    /* -------------------------------------------------------------- */

	    Lp = Lbp [block] ;
	    lnz += Lp [nk] ;
	    Up = Ubp [block] ;
	    unz += Up [nk] ;
	    noffdiag += noffdiag2 ;

	    /* -------------------------------------------------------------- */
	    /* combine the klu ordering with the symbolic pre-ordering */
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

	    /* the local pivot row permutation is no longer needed */
	    FREE (Pblock, int) ;
	}
    }

    /* finish the off-diagonal part */
    Offp [n] = poff ;
    ASSERT (Symbolic->nzoff == poff) ;
    PRINTF (("\n------------------- Off diagonal entries:\n")) ;
    ASSERT (klu_valid (n, Offp, Offi, Offx)) ;

    Numeric->lnz = lnz ;
    Numeric->unz = unz ;
    Numeric->noffdiag = noffdiag ;

    /* compute the inverse of Pnum */
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

    /* apply the pivot row permutations to the off-diagonal entries */
    for (p = 0 ; p < poff ; p++)
    {
	oldrow = Offi [p] ;
	newrow = Pinv [oldrow] ;
	Offi [p] = newrow ;
    }

    PRINTF (("\n------------------- Off diagonal entries, new:\n")) ;
    ASSERT (klu_valid (n, Offp, Offi, Offx)) ;

    return (TRUE) ;
}

/* ========================================================================== */
/* === klu_btf_factor ======================================================= */
/* ========================================================================== */

#define CLEAR(X,size,type) \
{ \
    int i ; \
    if (X != (type **) NULL) \
    { \
	for (i = 0 ; i < size ; i++) \
	{ \
	    X [i] = (type *) NULL ; \
	} \
    } \
}

klu_numeric *klu_btf_factor	/* returns NULL if error, or a valid
				   klu_numeric object if successful */
(
    /* inputs, not modified */
    int Ap [ ],		/* size n+1, column pointers */
    int Ai [ ],		/* size nz, row indices */
    double Ax [ ],
    double tol,
    klu_symbolic *Symbolic
)
{
    int n, nzoff, nblocks, maxblock, maxnz, *Cp, *Ci, *Pinv, ok, lnz, unz ;
    klu_numeric *Numeric ;
    double *Cx ;

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
    /* allocate the Numeric object */
    /* ---------------------------------------------------------------------- */

    Numeric = (klu_numeric *) ALLOCATE (sizeof (klu_numeric)) ;
    if (Numeric == (klu_numeric *) NULL)
    {
	return ((klu_numeric *) NULL) ;
    }
    Numeric->nblocks = nblocks ;
    Numeric->Pnum = (int *) ALLOCATE (n * sizeof (int)) ;
    Numeric->Offp = (int *) ALLOCATE ((n+1) * sizeof (int)) ;
    Numeric->Offi = (int *) ALLOCATE ((nzoff+1) * sizeof (int)) ;
    Numeric->Offx = (double *) ALLOCATE ((nzoff+1) * sizeof (double)) ;
    Numeric->Singleton = (double *) ALLOCATE (nblocks * sizeof (double)) ;
    Numeric->Lbp = (int **) ALLOCATE (nblocks * sizeof (int *)) ; 
    Numeric->Lbi = (int **) ALLOCATE (nblocks * sizeof (int *)) ; 
    Numeric->Lbx = (double **) ALLOCATE (nblocks * sizeof (double *)) ; 
    Numeric->Ubp = (int **) ALLOCATE (nblocks * sizeof (int *)) ; 
    Numeric->Ubi = (int **) ALLOCATE (nblocks * sizeof (int *)) ; 
    Numeric->Ubx = (double **) ALLOCATE (nblocks * sizeof (double *)) ; 
    Numeric->Rs = (double *) ALLOCATE (n * sizeof (double)) ; 

    /* clear the pointer arrays, so that klu_btf_free_numeric works OK */
    CLEAR (Numeric->Lbp, nblocks, int) ;
    CLEAR (Numeric->Lbi, nblocks, int) ;
    CLEAR (Numeric->Lbx, nblocks, double) ;
    CLEAR (Numeric->Ubp, nblocks, int) ;
    CLEAR (Numeric->Ubi, nblocks, int) ;
    CLEAR (Numeric->Ubx, nblocks, double) ;


    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    Cp = (int *) ALLOCATE ((maxblock + 1) * sizeof (int)) ;
    Ci = (int *) ALLOCATE ((maxnz + 1) * sizeof (int)) ;
    Cx = (double *) ALLOCATE ((maxnz + 1) * sizeof (double)) ;
    Pinv = (int *) ALLOCATE (n * sizeof (int)) ;

    if ((Numeric->Pnum == (int *) NULL) || (Numeric->Offp == (int *) NULL) ||
	(Numeric->Offi == (int *) NULL) || (Numeric->Offx == (double *) NULL) ||
	(Numeric->Singleton == (double *) NULL) ||
	(Numeric->Lbp == (int **) NULL) || (Numeric->Lbi == (int **) NULL) ||
	(Numeric->Lbx == (double **) NULL) || (Numeric->Ubp == (int **) NULL) ||
	(Numeric->Ubi == (int **) NULL) || (Numeric->Ubx == (double **) NULL) ||
	(Cp == (int *) NULL) || (Ci == (int *) NULL) ||
	(Cx == (double *) NULL) || (Pinv == (int *) NULL))
    {
	FREE (Cp, int) ;
	FREE (Ci, int) ;
	FREE (Cx, double) ;
	FREE (Pinv, int) ;
	klu_btf_free_numeric (&Numeric) ;
	return ((klu_numeric *) NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* factorize the blocks */
    /* ---------------------------------------------------------------------- */

    PRINTF (("calling klu_btf_factor2\n")) ;
    ok = klu_btf_factor2 (Ap, Ai, Ax, tol, Symbolic, Numeric,
	    Cp, Ci, Cx, Pinv) ;
    PRINTF (("klu_btf_factor2 done\n")) ;

    /* ---------------------------------------------------------------------- */
    /* free workspace */
    /* ---------------------------------------------------------------------- */

    FREE (Cp, int) ;
    FREE (Ci, int) ;
    FREE (Cx, double) ;
    FREE (Pinv, int) ;

    /* ---------------------------------------------------------------------- */
    /* return the numeric object */
    /* ---------------------------------------------------------------------- */

    if (!ok)
    {
	PRINTF (("klu_btf_factor failed!\n")) ;
	klu_btf_free_numeric (&Numeric) ;
	return ((klu_numeric *) NULL) ;
    }
    else
    {

#ifndef NDEBUG
	int block, k1, k2, nk ;
	PRINTF (("\n ############# KLU_BTF_FACTOR done, nblocks %d\n",nblocks));
	for (block = 0 ; block < nblocks ; block++)
	{
	    k1 = Symbolic->R [block] ;
	    k2 = Symbolic->R [block+1] ;
	    nk = k2 - k1 ;
	    PRINTF (("\n======================klu_btf_factor output: k1 %d k2 %d nk %d\n",k1,k2,nk)) ;
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
#endif

	return (Numeric) ;
    }
}
