/* ========================================================================== */
/* === klu_btf_factor ======================================================= */
/* ========================================================================== */

/* Factor the matrix, after ordering and analyzing it with klu_btf_analyze
 * or klu_btf_analyze_given
 *
 * TODO: provide a user switch to control scaling (none, row-sum, max row).
 * TODO: error checking of inputs.
 * TODO: test error cases / make sure there are no memory leaks
 * TODO: merge adjacent 1-by-1 blocks into a single upper triangular block,
 *	    for faster forward/backsolves
 * TODO: provide other orderings (natural, nested dissection, ...)
 * TODO: handle small or dense diagonal blocks using dense LU
 */

#include "klu_btf_internal.h"

/* ========================================================================== */

static int klu_btf_factor2
(
    /* inputs, not modified */
    int Ap [ ],		/* size n+1, column pointers */
    int Ai [ ],		/* size nz, row indices */
    double Ax [ ],
    double tol,
    double growth,
    double initmem_amd,
    double initmem_colamd,
    klu_symbolic *Symbolic,

    /* outputs, not defined on input */
    klu_numeric *Numeric,
    double Info [ ]
)
{
    double klu_Control [KLU_CONTROL], umin, umax, umin2, umax2,
	*Lnz, *Singleton, **Lbx, **Ubx, *Offx, *Rs, *X ;
    int k1, k2, nk, k, block, oldcol, pend, oldrow, n, lnz, unz,
	result, p, newrow, *P, *Q, *R, nblocks, poff, *Pnum, *Lp, *Up,
	*Offp, *Offi, **Lbp, **Lbi, **Ubp, **Ubi, *Pblock, noffdiag2, noffdiag,
	*Pinv, *Iwork, nzoff ;

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
    X = Numeric->Xwork ;			/* X is of size n */
    Iwork = Numeric->Iwork ;			/* 5*maxblock for klu_factor */
    Pblock = Iwork + 5*(Symbolic->maxblock) ;	/* 1*maxblock for Pblock */

    /* compute the inverse of P from symbolic analysis.  Will be updated to
     * become the inverse of the numerical factorization when the factorization
     * is done, for use in klu_btf_refactor */
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
    klu_defaults (klu_Control) ;
    klu_Control [KLU_TOL] = tol ;
    klu_Control [KLU_GROWTH] = growth ;
    lnz = 0 ;
    unz = 0 ;
    noffdiag = 0 ;
    umin = 0 ;
    umax = 0 ;
    Offp [0] = 0 ;

    /* ---------------------------------------------------------------------- */
    /* check the input matrix and compute the row scale factors, Rs */
    /* ---------------------------------------------------------------------- */

    /* use Pnum as workspace */
    result = klu_btf_scale (n, Ap, Ai, Ax, Rs, Pnum) ;
    Info [KLU_BTF_INFO_STATUS] = result ;
    if (result != KLU_OK)
    {
	/* matrix is invalid */
	return (result) ;
    }

#ifndef NDEBUG
    for (k = 0 ; k < n ; k++) PRINTF (("Rs [%d] = %g\n", k, Rs [k])) ;
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
	    Offp [k1+1] = poff ;
	    Pnum [k1] = P [k1] ;
	    lnz++ ;
	    unz++ ;

	    umin2 = Singleton [block] ;
	    umax2 = Singleton [block] ;

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* construct and factorize the kth block */
	    /* -------------------------------------------------------------- */

	    PRINTF (("calling klu\n")) ;
	    if (Lnz [block] == EMPTY)
	    {
		/* COLAMD was used - no estimate of fill-in */
		/* use 10 times the nnz in A, plus n */
		klu_Control [KLU_LSIZE] = -initmem_colamd ;
		klu_Control [KLU_USIZE] = klu_Control [KLU_LSIZE] ;
	    }
	    else
	    {
		klu_Control [KLU_LSIZE] = initmem_amd * Lnz [block] + nk ;
		klu_Control [KLU_USIZE] = klu_Control [KLU_LSIZE] ;
	    }

	    /* allocates 4 arrays:
	     * Lbi [block], Lbx [block], Ubi [block], Ubx [block] */
	    result = klu_factor (nk, Ap, Ai, Ax, Q, klu_Control,
		    Lbp [block], &Lbi [block], &Lbx [block],
		    Ubp [block], &Ubi [block], &Ubx [block], Pblock,
		    &noffdiag2, &umin2, &umax2, X, Iwork,
		    /* BTF-related arguments: */
		    k1, Pinv, Rs, Offp, Offi, Offx) ;

	    PRINTF (("klu done\n")) ;
	    if (result != KLU_OK)
	    {
		/* The 4 arrays (Li,Lx,Ui,Ux) are in the Numeric object and
		 * will be free'd later in klu_btf_free_numeric. */
		Info [KLU_BTF_INFO_STATUS] = result ;
		return (result) ;
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

	    if (Lnz [block] == EMPTY)
	    {
		/* revise estimate for subsequent factorization */
		Lnz [block] = Lp [nk] ;
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

	/* keep track of the largest and smallest diagonal entry of U */
	umin = (block == 0) ? umin2 : MIN (umin, umin2) ;
	umax = (block == 0) ? umax2 : MAX (umax, umax2) ;

    }

    ASSERT (nzoff == Offp [n]) ;
    PRINTF (("\n------------------- Off diagonal entries:\n")) ;
    ASSERT (klu_valid (n, Offp, Offi, Offx)) ;

    Numeric->lnz = lnz ;
    Numeric->unz = unz ;
    Numeric->noffdiag = noffdiag ;

    Info [KLU_BTF_INFO_UMIN] = umin ;
    Info [KLU_BTF_INFO_UMAX] = umax ;
    Info [KLU_BTF_INFO_LNZ] = lnz ;
    Info [KLU_BTF_INFO_UNZ] = unz ;
    Info [KLU_BTF_INFO_FLOPS] = EMPTY ;		/* TODO not yet computed */
    Info [KLU_BTF_INFO_REALLOC] = EMPTY ;	/* TODO not yet computed */
    Info [KLU_BTF_INFO_NOFFDIAG] = noffdiag ;

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

    /* permute scale factors Rs according to pivotal row order */
    for (k = 0 ; k < n ; k++)
    {
	oldrow = Pnum [k] ;
	X [k] = Rs [oldrow] ;
    }
    for (k = 0 ; k < n ; k++)
    {
	Rs [k] = X [k] ;
    }

    PRINTF (("\n------------------- Off diagonal entries, old:\n")) ;
    ASSERT (klu_valid (n, Offp, Offi, Offx)) ;

    /* apply the pivot row permutations to the off-diagonal entries */
    for (p = 0 ; p < nzoff ; p++)
    {
	ASSERT (Offi [p] >= 0 && Offi [p] < n) ;
	Offi [p] = Pinv [Offi [p]] ;
    }

    PRINTF (("\n------------------- Off diagonal entries, new:\n")) ;
    ASSERT (klu_valid (n, Offp, Offi, Offx)) ;

    return (KLU_OK) ;
}

/* ========================================================================== */
/* === klu_btf_factor ======================================================= */
/* ========================================================================== */

#define CLEAR(Ptr,size,type) \
{ \
    int ii ; \
    if (Ptr != (type **) NULL) \
    { \
	for (ii = 0 ; ii < size ; ii++) \
	{ \
	    Ptr [ii] = (type *) NULL ; \
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
    klu_symbolic *Symbolic,

    double Control [KLU_BTF_CONTROL],	/* optional; may be NULL */
    double User_Info [KLU_BTF_INFO]	/* optional; may be NULL */
)
{
    double *Info, Info2 [KLU_BTF_INFO], tol, growth, initmem_amd,
	initmem_colamd ;
    int n, nzoff, nblocks, maxblock, result, i, *R, block, k1, k2, nk,
	**Lbp, **Ubp ;
    klu_numeric *Numeric ;

    /* ---------------------------------------------------------------------- */
    /* get Info array, and clear the parts used by klu_btf_factor */
    /* ---------------------------------------------------------------------- */

    Info = (User_Info != (double *) NULL) ? User_Info : Info2 ;
    for (i = KLU_BTF_INFO_LNZ ; i <= KLU_BTF_INFO_REALLOC ; i++)
    {
	Info [i] = EMPTY ;
    }

    /* ---------------------------------------------------------------------- */
    /* get the contents of the Symbolic object */
    /* ---------------------------------------------------------------------- */

    /* TODO:  check for a valid Symbolic object */

    n = Symbolic->n ;
    nzoff = Symbolic->nzoff ;
    nblocks = Symbolic->nblocks ;
    maxblock = Symbolic->maxblock ;
    R = Symbolic->R ;
    PRINTF (("klu_btf_factor:  n %d nzoff %d nblocks %d maxblock %d\n",
	n, nzoff, nblocks, maxblock)) ;

    /* ---------------------------------------------------------------------- */
    /* get control parameters and make sure they are in the proper range */
    /* ---------------------------------------------------------------------- */

    tol = GET_CONTROL (
	    KLU_BTF_CONTROL_TOL,
	    KLU_BTF_CONTROL_TOL_DEFAULT) ;

    growth = GET_CONTROL (
	    KLU_BTF_CONTROL_GROWTH,
	    KLU_BTF_CONTROL_GROWTH_DEFAULT) ;

    initmem_amd = GET_CONTROL (
	    KLU_BTF_CONTROL_INITMEM_AMD,
	    KLU_BTF_CONTROL_INITMEM_AMD_DEFAULT) ;

    initmem_colamd = GET_CONTROL (
	    KLU_BTF_CONTROL_INITMEM_COLAMD,
	    KLU_BTF_CONTROL_INITMEM_COLAMD_DEFAULT) ;

    initmem_amd    = MAX (1.0, initmem_amd) ;
    initmem_colamd = MAX (1.0, initmem_colamd) ;
    tol    = MIN (tol, 1.0) ;
    tol    = MAX (0.0, tol) ;
    growth = MAX (1.0, growth) ;

    /* ---------------------------------------------------------------------- */
    /* allocate the Numeric object (except Li,Lx,Ui,Ux) */
    /* ---------------------------------------------------------------------- */

    Numeric = (klu_numeric *) ALLOCATE (sizeof (klu_numeric)) ;
    if (Numeric == (klu_numeric *) NULL)
    {
	Info [KLU_BTF_INFO_STATUS] = KLU_OUT_OF_MEMORY ;
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
    Numeric->Pinv = (int *) ALLOCATE (n * sizeof (int)) ; 

    /* allocate permanent workspace for factorization and solve.  On typical
     * computers (sizeof (double) = 2*sizeof (int), this is exactly
     * 4n double's.  Note that the solver will use an X of size 4n, whereas
     * the factorization codes use an X of size n and integer space (Iwork)
     * of size 6n. */
    Numeric->worksize = n * sizeof (double) +
	MAX (3*n * sizeof (double), 6*maxblock * sizeof (int)) ;
    Numeric->Work = (void *) ALLOCATE (Numeric->worksize) ; 
    Numeric->Xwork = (double *) Numeric->Work ;
    Numeric->Iwork = (int *) (Numeric->Xwork + n) ;

    /* clear the pointer arrays, so that klu_btf_free_numeric works OK */
    CLEAR (Numeric->Lbp, nblocks, int) ;
    CLEAR (Numeric->Lbi, nblocks, int) ;
    CLEAR (Numeric->Lbx, nblocks, double) ;
    CLEAR (Numeric->Ubp, nblocks, int) ;
    CLEAR (Numeric->Ubi, nblocks, int) ;
    CLEAR (Numeric->Ubx, nblocks, double) ;

    if ((Numeric->Pnum == (int *) NULL) || (Numeric->Offp == (int *) NULL) ||
	(Numeric->Offi == (int *) NULL) || (Numeric->Offx == (double *) NULL) ||
	(Numeric->Singleton == (double *) NULL) ||
	(Numeric->Lbp == (int **) NULL) || (Numeric->Lbi == (int **) NULL) ||
	(Numeric->Lbx == (double **) NULL) || (Numeric->Ubp == (int **) NULL) ||
	(Numeric->Ubi == (int **) NULL) || (Numeric->Ubx == (double **) NULL) ||
	(Numeric->Pinv == (int *) NULL))
    {
	klu_btf_free_numeric (&Numeric) ;
	Info [KLU_BTF_INFO_STATUS] = KLU_OUT_OF_MEMORY ;
	return ((klu_numeric *) NULL) ;
    }

    /* allocate the column pointer arrays for each block */
    Lbp = Numeric->Lbp ;
    Ubp = Numeric->Ubp ;
    for (block = 0 ; block < nblocks ; block++)
    {
	k1 = R [block] ;
	k2 = R [block+1] ;
	nk = k2 - k1 ;
	if (nk > 1)
	{
	    Lbp [block] = (int *) ALLOCATE ((nk+1) * sizeof (int)) ;
	    Ubp [block] = (int *) ALLOCATE ((nk+1) * sizeof (int)) ;
	    if ((Lbp [block] == (int *) NULL) || (Ubp [block] == (int *) NULL))
	    {
		klu_btf_free_numeric (&Numeric) ;
		Info [KLU_BTF_INFO_STATUS] = KLU_OUT_OF_MEMORY ;
		return ((klu_numeric *) NULL) ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* factorize the blocks */
    /* ---------------------------------------------------------------------- */

    PRINTF (("calling klu_btf_factor2\n")) ;
    result = klu_btf_factor2 (Ap, Ai, Ax, tol, growth,
	    initmem_amd, initmem_colamd, Symbolic, Numeric, Info) ;
    PRINTF (("klu_btf_factor2 done\n")) ;
    Info [KLU_BTF_INFO_STATUS] = result ;

    /* ---------------------------------------------------------------------- */
    /* return the numeric object */
    /* ---------------------------------------------------------------------- */

    if (result != KLU_OK)
    {
	PRINTF (("klu_btf_factor failed!\n")) ;
	klu_btf_free_numeric (&Numeric) ;
	return ((klu_numeric *) NULL) ;
    }

#ifndef NDEBUG
    {
	PRINTF (("\n ############# KLU_BTF_FACTOR done, nblocks %d\n",nblocks));
	for (block = 0 ; block < nblocks ; block++)
	{
	    k1 = R [block] ;
	    k2 = R [block+1] ;
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
		Lp = Lbp [block] ;
		Li = Numeric->Lbi [block] ;
		Lx = Numeric->Lbx [block] ;
		PRINTF (("\n---- L block %d\n", block)); 
		ASSERT (klu_valid (nk, Lp, Li, Lx)) ;
		Up = Ubp [block] ;
		Ui = Numeric->Ubi [block] ;
		Ux = Numeric->Ubx [block] ;
		PRINTF (("\n---- U block %d\n", block)) ; 
		ASSERT (klu_valid (nk, Up, Ui, Ux)) ;
	    }
	}
    }
#endif

    return (Numeric) ;
}
