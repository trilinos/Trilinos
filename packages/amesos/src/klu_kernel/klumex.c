/* ========================================================================== */
/* === klu mexFunction ====================================================== */
/* ========================================================================== */

/* Factorizes A(p,q) into L*U.  Usage:
 *
 * [L,U,p,t] = klu (A)			% Q assumed to be identity
 * [L,U,p,t] = klu (A, q)
 * [L,U,p,t] = klu (A, q, Control)
 *
 * Control (1): threshold partial pivoting tolerance (0 means force diagonal
 *		pivoting, 1 means conventional partial pivoting.  A value of
 *		0.5 means the diagonal will be picked if its absolute value
 *		is >= 0.5 times the largest value in the column; otherwise
 *		the largest value is selected.  Default: 1.0.
 * Control (2): if positive, this is the initial size of the L matrix, in # of
 *		nonzero entries.  If negative, the initial size of L is
 *		(-Control (2) * nnz (A)).  Default: -10.
 * Control (3): if positive, this is the initial size of the U matrix, in # of
 *		nonzero entries.  If negative, the initial size of U is
 *		(-Control (2) * nnz (A)).  Default: -10.
 * Control (4): memory growth factor
 *
 * t = [cputime noffdiag umin umax], output statistics.
 *
 * If Control is not present, or present but not of the correct size,
 * defaults are used.
 */

/* ========================================================================== */

#include "klu.h"
#include "mex.h"
#include "tictoc.h"

void mexFunction
(
    int	nargout,
    mxArray *pargout [ ],
    int	nargin,
    const mxArray *pargin [ ]
)
{
    int n, *Ap, *Ai, *Lp, *Li, *Up, *Ui, result, *P, *Pp, *Pi, k, anz, s,
	*Lp2, *Li2, *Up2, *Ui2, lnz, unz, j, p, *Q, col, noffdiag ;
    double *Lx, *Ux, *Ax, *Px, Control [KLU_CONTROL], *User_Control,
	*Lx2, *Ux2, tt [2], *T, *Qx, umin, umax ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    if (nargin < 1 || nargin > 3 || !(nargout == 3 || nargout == 4))
    {
	mexErrMsgTxt ("Usage: [L,U,P,t] = klu (A,Q,Control), where t, Q, and Control are optional.") ;
    }
    n = mxGetM (pargin [0]) ;
    if (!mxIsSparse (pargin [0]) || n != mxGetN (pargin [0])
	|| mxIsComplex (pargin [0]) || n == 0)
    {
    	mexErrMsgTxt ("klu: A must be sparse, square, real, and non-empty") ;
    }

    /* get sparse matrix A */
    Ap = mxGetJc (pargin [0]) ;
    Ai = mxGetIr (pargin [0]) ;
    Ax = mxGetPr (pargin [0]) ;
    anz = Ap [n] ;

    /* get input column permutation Q */
    if (nargin > 1)
    {
	if (!mxIsDouble (pargin [1]) || n != mxGetNumberOfElements (pargin [1]))
	{
	    mexErrMsgTxt ("klu: Q must be a dense 1-by-n vector") ; 
	}
	Qx = mxGetPr (pargin [1]) ;
	Q = mxMalloc (n * sizeof (int)) ;
	for (k = 0 ; k < n ; k++)
	{
	    col = (int) (Qx [k]) - 1 ;	/* convert from 1-based to 0-based */
	    if (col < 0 || col >= n)
	    {
		mexErrMsgTxt ("klu: Q not a valid permutation\n") ;
	    }
	    Q [k] = col ;
	}
    }
    else
    {
	/* klu will assume that Q is the identity permutation */
	Q = (int *) NULL ;
    }

    /* get control parameters */
    klu_defaults (Control) ;
    if (nargin > 2)
    {
	if (!mxIsDouble (pargin [2]))
	{
	    mexErrMsgTxt ("klu: Control must be real") ;
	}
	User_Control = mxGetPr (pargin [2]) ;
	s = mxGetNumberOfElements (pargin [2]) ;
	for (k = 0 ; k < s ; k++)
	{
	    Control [k] = User_Control [k] ;
	}
    }

    Lp = mxMalloc ((n+1) * sizeof (int)) ;
    Up = mxMalloc ((n+1) * sizeof (int)) ;
    P  = mxMalloc (n * sizeof (int)) ;

    /* ---------------------------------------------------------------------- */
    /* factorize */
    /* ---------------------------------------------------------------------- */

    my_tic (tt) ;

    result = klu_factor (n, Ap, Ai, Ax, Q, Control, Lp, &Li, &Lx, Up, &Ui,
	    &Ux, P, &noffdiag, &umin, &umax, (double *) NULL, (int *) NULL,
	    /* no BTF or scaling here */
	    0, (int *) NULL, (double *) NULL, 0, (int *) NULL, (int *) NULL,
	    (double *) NULL) ;

    my_toc (tt) ;

    if (nargout == 4)
    {
	pargout [3] = mxCreateDoubleMatrix (1, 4, mxREAL) ;
	T = mxGetPr (pargout [3]) ;
	T [0] = tt [1] ;
	T [1] = noffdiag ;
	T [2] = umin ;
	T [3] = umax ;
    }

    if (result == KLU_OUT_OF_MEMORY)
    {
	mexErrMsgTxt ("klu: out of memory") ;
    }

    /* NOTE: this should be done in a better way, without copying, but when
     * I tried to do it I got segmentation faults, and gave up ... */

    /* create sparse matrix for L */
    lnz = Lp [n] ;
    pargout [0] = mxCreateSparse (n, n, lnz, mxREAL) ;
    Lp2 = mxGetJc (pargout [0]) ;
    Li2 = mxGetIr (pargout [0]) ;
    Lx2 = mxGetPr (pargout [0]) ;
    for (j = 0 ; j <= n ; j++)
    {
	Lp2 [j] = Lp [j] ;
    }
    for (p = 0 ; p < lnz ; p++)
    {
	Li2 [p] = Li [p] ;
	Lx2 [p] = Lx [p] ;
    }
    mxFree (Lp) ;
    mxFree (Li) ;
    mxFree (Lx) ;

    /* create sparse matrix for U */
    unz = Up [n] ;
    pargout [1] = mxCreateSparse (n, n, unz, mxREAL) ;
    Up2 = mxGetJc (pargout [1]) ;
    Ui2 = mxGetIr (pargout [1]) ;
    Ux2 = mxGetPr (pargout [1]) ;
    for (j = 0 ; j <= n ; j++)
    {
	Up2 [j] = Up [j] ;
    }
    for (p = 0 ; p < unz ; p++)
    {
	Ui2 [p] = Ui [p] ;
	Ux2 [p] = Ux [p] ;
    }
    mxFree (Up) ;
    mxFree (Ui) ;
    mxFree (Ux) ;

    /* create permutation vector for P */
    pargout [2] = mxCreateDoubleMatrix (1, n, mxREAL) ;
    Px = mxGetPr (pargout [2]) ;
    for (k = 0 ; k < n ; k++)
    {
	Px [k] = P [k] + 1 ;	/* convert to 1-based */
    }

    mxFree (P) ;

}
