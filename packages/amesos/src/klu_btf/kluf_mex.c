/* ========================================================================== */
/* === kluf mexFunction ===================================================== */
/* ========================================================================== */

/* Factor A using symbolic info from klua.
 *
 * [P,Q,R,Lnz,Info1] = klua (A) ;
 * [L,U,Off,Pnum,Rs,Info2] = kluf (A, P,Q,R,Lnz,Info1, Control) ;
 *
 * The factorization is L*U + Off = Rs (Pnum,Pnum) \ (A (Pnum,Q)), where Rs is
 * a diagonal matrix of row scale factors.  If Pnum and Q are converted to
 * permutation matrices, then L*U + Off = Pnum * (Rs\A) * Q. 
 */


/* ========================================================================== */

#include "klu_btf_internal.h"

void mexFunction
(
    int	nargout,
    mxArray *pargout [ ],
    int	nargin,
    const mxArray *pargin [ ]
)
{
    int n, *Ap, *Ai, k, p, col, block, nblocks, *Off2i, *Off2p, *L2i, *L2p, nk,
	*U2p, *U2i, *Lp, *Li, *Up, *Ui, lnz, unz, k1, k2, nz, pl, pu, nzoff,
	*Rsi, *Rsp, i ;
    double *Ax, *Px, *Qx, *Rx, *Lnzx, *Info_in, *Pnumx, *Lx, *Ux, *Off2x, *L2x,
	*U2x, *Rsx, *Info, *Control ;
    klu_symbolic *Symbolic ;
    klu_numeric *Numeric ;
    klu_control control ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    if (nargin < 6 || nargin > 7 || nargout != 6)
    {
	mexErrMsgTxt ("Usage: [L,U,Off,Pnum,Rs,Info] = "
		      "kluf (A,P,Q,R,Lnz,Info_in,Control)") ;
    }
    n = mxGetM (pargin [0]) ;
    if (!mxIsSparse (pargin [0]) || n != mxGetN (pargin [0]))
    {
    	mexErrMsgTxt ("klua: A must be sparse, square, real, and non-empty") ;
    }

    /* get sparse matrix A */
    Ap = mxGetJc (pargin [0]) ;
    Ai = mxGetIr (pargin [0]) ;
    Ax = mxGetPr (pargin [0]) ;
    nz = Ap [n] ;

    /* get control parameters */
    klu_btf_defaults (&control) ;
    if (nargin > 6)
    {
	int s ;
	if (!mxIsDouble (pargin [6]))
	{
	    mexErrMsgTxt ("klu: control must be real") ;
	}
	Control = mxGetPr (pargin [6]) ;
	s = mxGetNumberOfElements (pargin [6]) ;
	/* if (s > 0) control.prl      = Control [0] ; */
	if (s > 1) control.btf         = Control [1] ;
	if (s > 2) control.scale       = Control [2] ;
	if (s > 3) control.ordering    = Control [3] ;
	if (s > 4) control.tol         = Control [4] ;
	if (s > 5) control.growth      = Control [5] ;
	if (s > 6) control.initmem_amd = Control [6] ;
	if (s > 7) control.initmem     = Control [7] ;
    }
    PRINTF (("control: btf %d ord %d tol %g gro %g inita %g init %g\n",
	 control.btf, control.ordering, control.tol, control.growth,
	 control.initmem_amd, control.initmem)) ;

    /* ---------------------------------------------------------------------- */
    /* reconstruct the symbolic object */
    /* ---------------------------------------------------------------------- */

    Symbolic = (klu_symbolic *) ALLOCATE (sizeof (klu_symbolic)) ;

    /* get Info */
    Info_in = mxGetPr (pargin [5]) ;
    pargout [5] = mxCreateDoubleMatrix (1, 90, mxREAL) ;
    Info = mxGetPr (pargout [5]) ;
    for (i = 0 ; i < 90 ; i++)
    {
	Info [i] = Info_in [i] ;
    }

    Symbolic->n        = n ;		/* dimension of A */
    Symbolic->nz       = nz ;		/* # entries in input matrix */
    Symbolic->nblocks  = Info_in [ 3] ; /* # of blocks in BTF form */
    Symbolic->maxblock = Info_in [ 4] ; /* dimension of largest block */
    Symbolic->nzoff    = Info_in [ 7] ; /* nz in off-diagonal blocks of A */
    Symbolic->symmetry = Info_in [ 8] ; /* symmetry of largest block */
    Symbolic->lnz      = Info_in [10] ; /* nz in L, estimated (incl diagonal) */
    Symbolic->unz      = Info_in [11] ; /* nz in U, estimated (incl diagonal) */
    Symbolic->est_flops= Info_in [12] ; /* est. factorization flop count */

    PRINTF (("kluf: n %d nzoff %d nblocks %d maxblock %d\n",
	Symbolic->n, Symbolic->nzoff, Symbolic->nblocks,
	Symbolic->maxblock)) ; 

    nblocks = Symbolic->nblocks ;
    ASSERT (nblocks > 0) ;

    Symbolic->P = (int *) ALLOCATE (n * sizeof (int)) ;
    Symbolic->Q = (int *) ALLOCATE (n * sizeof (int)) ;
    Symbolic->R = (int *) ALLOCATE ((nblocks+1) * sizeof (int)) ;
    Symbolic->Lnz = (double *) ALLOCATE (nblocks * sizeof (double)) ;

    ASSERT (Symbolic->nzoff >= 0 && Symbolic->nzoff <= nz) ;
    ASSERT (Symbolic->maxblock > 0 && Symbolic->maxblock <= n - nblocks + 1) ;

    /* get P */
    PRINTF (("get P\n")) ;
    Px = mxGetPr (pargin [1]) ;
    for (k = 0 ; k < n ; k++)
    {
	Symbolic->P [k] = Px [k] - 1 ;	    /* convert to 0-based */
    }

    /* get Q */
    PRINTF (("get P\n")) ;
    Qx = mxGetPr (pargin [2]) ;
    for (k = 0 ; k < n ; k++)
    {
	Symbolic->Q [k] = Qx [k] - 1 ;	    /* convert to 0-based */
    }

    /* get R */
    PRINTF (("get P\n")) ;
    Rx = mxGetPr (pargin [3]) ;
    for (k = 0 ; k <= nblocks ; k++)
    {
	Symbolic->R [k] = Rx [k] - 1 ;	    /* convert to 0-based */
    }
    ASSERT (Symbolic->R [nblocks] == n) ;

    /* get Lnz */
    PRINTF (("get P\n")) ;
    Lnzx = mxGetPr (pargin [4]) ;
    for (k = 0 ; k < nblocks ; k++)
    {
	Symbolic->Lnz [k] = Lnzx [k] ;
    }

    /* ---------------------------------------------------------------------- */
    /* factorize */
    /* ---------------------------------------------------------------------- */

    PRINTF (("calling klu_btf_factor\n")) ;
    Numeric = klu_btf_factor (Ap, Ai, Ax, Symbolic, &control) ;
    PRINTF (("klu_btf_factor done\n")) ;

    if (Numeric == (klu_numeric *) NULL)
    {
	mexErrMsgTxt ("kluf: error") ;
    }

#ifndef NDEBUG
    /* dump */
    for (block = 0 ; block < nblocks ; block++)
    {
	k1 = Symbolic->R [block] ;
	k2 = Symbolic->R [block+1] ;
	nk = k2 - k1 ;
	PRINTF (("block %d k1 %d k2 %d nk %d\n", block, k1, k2, nk)) ;
	if (nk > 1)
	{
	    Lp = Numeric->Lbp [block] ;
	    Li = Numeric->Lbi [block] ;
	    Lx = Numeric->Lbx [block] ;
	    Up = Numeric->Ubp [block] ;
	    Ui = Numeric->Ubi [block] ;
	    Ux = Numeric->Ubx [block] ;
	    PRINTF (("\n---- L block %d\n", block)) ; 
	    ASSERT (klu_valid (nk, Lp, Li, Lx)) ;
	    PRINTF (("\n---- U block %d\n", block)) ; 
	    ASSERT (klu_valid (nk, Up, Ui, Ux)) ;
	}
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* return the results to MATLAB */
    /* ---------------------------------------------------------------------- */

    /* create Info output */
    Info [30] = Numeric->lnz ;		/* nz in L, actual (incl. diagonal) */
    Info [31] = Numeric->unz ;		/* nz in U, actual (incl. diagonal) */
    Info [36] = Numeric->noffdiag ;	/* number of off-diagonal pivots */
    Info [33] = Numeric->umin ;		/* min abs diagonal entry in U */
    Info [34] = Numeric->umax ;		/* max abs diagonal entry in U */

    /* create permutation vector for Pnum */
    pargout [3] = mxCreateDoubleMatrix (1, n, mxREAL) ;
    Pnumx = mxGetPr (pargout [3]) ;
    for (k = 0 ; k < n ; k++)
    {
	Pnumx [k] = Numeric->Pnum [k] + 1 ;	/* convert to 1-based */
	PRINTF (("Pnum (%d) = %g\n", k+1, Pnumx [k])) ;
    }

#ifndef NDEBUG
    /* dump again */
    for (block = 0 ; block < nblocks ; block++)
    {
	k1 = Symbolic->R [block] ;
	k2 = Symbolic->R [block+1] ;
	nk = k2 - k1 ;
	PRINTF (("again, block %d k1 %d k2 %d nk %d\n", block, k1, k2, nk)) ;
	if (nk > 1)
	{
	    Lp = Numeric->Lbp [block] ;
	    Li = Numeric->Lbi [block] ;
	    Lx = Numeric->Lbx [block] ;
	    Up = Numeric->Ubp [block] ;
	    Ui = Numeric->Ubi [block] ;
	    Ux = Numeric->Ubx [block] ;
	    PRINTF (("\n---- L block %d\n", block)) ; 
	    ASSERT (klu_valid (nk, Lp, Li, Lx)) ;
	    PRINTF (("\n---- U block %d\n", block)) ; 
	    ASSERT (klu_valid (nk, Up, Ui, Ux)) ;
	}
    }
#endif

    /* create Off */
    PRINTF (("\n----------------------- Off input:\n")) ;
    ASSERT (klu_valid (n, Numeric->Offp, Numeric->Offi, Numeric->Offx)) ;
    nzoff = Symbolic->nzoff ;
    pargout [2] = mxCreateSparse (n, n, nzoff+1, mxREAL) ;
    Off2p = mxGetJc (pargout [2]) ;
    Off2i = mxGetIr (pargout [2]) ;
    Off2x = mxGetPr (pargout [2]) ;
    for (col = 0 ; col <= n ; col++) Off2p [col] = Numeric->Offp [col] ;
    for (p = 0 ; p < nzoff ; p++)    Off2i [p]   = Numeric->Offi [p] ;
    for (p = 0 ; p < nzoff ; p++)    Off2x [p]   = Numeric->Offx [p] ;
    PRINTF (("\n----------------------- Off output:\n")) ;
    ASSERT (klu_valid (n, Off2p, Off2i, Off2x)) ;

#ifndef NDEBUG
    /* determine # of nonzeros in L and U */
    lnz = 0 ;
    unz = 0 ;
    for (block = 0 ; block < nblocks ; block++)
    {
	k1 = Symbolic->R [block] ;
	k2 = Symbolic->R [block+1] ;
	nk = k2 - k1 ;
	PRINTF (("block %d k1 %d k2 %d nk %d\n", block, k1, k2, nk)) ;
	if (nk == 1)
	{
	    lnz++ ;
	    unz++ ;
	}
	else
	{
	    Lp = Numeric->Lbp [block] ;
	    Li = Numeric->Lbi [block] ;
	    Lx = Numeric->Lbx [block] ;
	    Up = Numeric->Ubp [block] ;
	    Ui = Numeric->Ubi [block] ;
	    Ux = Numeric->Ubx [block] ;
	    PRINTF (("\n---- L block %d\n", block)) ; 
	    ASSERT (klu_valid (nk, Lp, Li, Lx)) ;
	    PRINTF (("\n---- U block %d\n", block)) ; 
	    ASSERT (klu_valid (nk, Up, Ui, Ux)) ;
	    lnz += Lp [nk] ;
	    unz += Up [nk] ;
	}
    }
    PRINTF (("Total lnz %d unz %d for all blocks\n", lnz, unz)) ;
    ASSERT (lnz == Numeric->lnz) ;
    ASSERT (unz == Numeric->unz) ;
#endif
    lnz = Numeric->lnz ;
    unz = Numeric->unz ;

    /* create L */
    pargout [0] = mxCreateSparse (n, n, lnz+1, mxREAL) ;
    L2p = mxGetJc (pargout [0]) ;
    L2i = mxGetIr (pargout [0]) ;
    L2x = mxGetPr (pargout [0]) ;

    /* create U */
    pargout [1] = mxCreateSparse (n, n, unz+1, mxREAL) ;
    U2p = mxGetJc (pargout [1]) ;
    U2i = mxGetIr (pargout [1]) ;
    U2x = mxGetPr (pargout [1]) ;

#ifndef NDEBUG
    /* dump again */
    for (block = 0 ; block < nblocks ; block++)
    {
	k1 = Symbolic->R [block] ;
	k2 = Symbolic->R [block+1] ;
	nk = k2 - k1 ;
	PRINTF (("yet again block %d k1 %d k2 %d nk %d\n", block, k1, k2, nk)) ;
	if (nk > 1)
	{
	    Lp = Numeric->Lbp [block] ;
	    Li = Numeric->Lbi [block] ;
	    Lx = Numeric->Lbx [block] ;
	    Up = Numeric->Ubp [block] ;
	    Ui = Numeric->Ubi [block] ;
	    Ux = Numeric->Ubx [block] ;
	    PRINTF (("\n---- L block %d\n", block)) ; 
	    ASSERT (klu_valid (nk, Lp, Li, Lx)) ;
	    PRINTF (("\n---- U block %d\n", block)) ; 
	    ASSERT (klu_valid (nk, Up, Ui, Ux)) ;
	}
    }
#endif

    /* fill L and U */
    pl = 0 ;
    pu = 0 ;
    for (block = 0 ; block < nblocks ; block++)
    {
	k1 = Symbolic->R [block] ;
	k2 = Symbolic->R [block+1] ;
	nk = k2 - k1 ;
	if (nk == 1)
	{
	    L2p [k1] = pl ;
	    L2i [pl] = k1 ;
	    L2x [pl] = 1 ;
	    pl++ ;

	    U2p [k1] = pu ;
	    U2i [pu] = k1 ;
	    U2x [pu] = Numeric->Singleton [block] ;
	    pu++ ;
	}
	else
	{
	    Lp = Numeric->Lbp [block] ;
	    Li = Numeric->Lbi [block] ;
	    Lx = Numeric->Lbx [block] ;
	    Up = Numeric->Ubp [block] ;
	    Ui = Numeric->Ubi [block] ;
	    Ux = Numeric->Ubx [block] ;

	    for (k = 0 ; k < nk ; k++)
	    {
		L2p [k+k1] = pl ;
		for (p = Lp [k] ; p < Lp [k+1] ; p++)
		{
		    L2i [pl] = Li [p] + k1 ;
		    L2x [pl] = Lx [p] ;
		    pl++ ;
		}

		U2p [k+k1] = pu ;
		for (p = Up [k] ; p < Up [k+1] ; p++)
		{
		    U2i [pu] = Ui [p] + k1 ;
		    U2x [pu] = Ux [p] ;
		    pu++ ;
		}
	    }
	}
    }

    L2p [n] = pl ;
    U2p [n] = pu ;
    ASSERT (pl == lnz) ;
    ASSERT (pu == unz) ;

    /* create Rs */
    pargout [4] = mxCreateSparse (n, n, n, mxREAL) ;
    Rsp = mxGetJc (pargout [4]) ;
    Rsi = mxGetIr (pargout [4]) ;
    Rsx = mxGetPr (pargout [4]) ;
    for (k = 0 ; k < n ; k++)
    {
	Rsp [k] = k ;
	Rsi [k] = k ;
	Rsx [k] = Numeric->Rs [k] ;
	PRINTF (("Rsx [k %d] %g\n", k, Rsx [k])) ;
    }
    Rsp [n] = n ;

    PRINTF (("\n------------------ All of output L:\n")) ;
    ASSERT (klu_valid (n, L2p, L2i, L2x)) ;
    PRINTF (("\n------------------ All of output U:\n")) ;
    ASSERT (klu_valid (n, U2p, U2i, U2x)) ;

    /* destroy the symbolic object */
    klu_btf_free_symbolic (&Symbolic) ;

    /* destroy the numeric object */
    klu_btf_free_numeric (&Numeric) ;
}
