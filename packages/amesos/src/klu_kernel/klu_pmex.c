/* ========================================================================== */
/* === klu_p mexFunction ==================================================== */
/* ========================================================================== */

/* computes X=B(p,:) using klu_permute, and where p is a permutation vector
 *
 * X = klu_p (p,B)
 */

/* ========================================================================== */

#include "klu.h"
#include "mex.h"

void mexFunction
(
    int	nargout,
    mxArray *pargout[],
    int	nargin,
    const mxArray *pargin[]
)
{
    double *Px, *X, *B ;
    int k, n, *P, nrhs ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    if (nargin != 2 || nargout != 1)
    {
	mexErrMsgTxt ("Usage: x = klu_p (p,b)") ;
    }
    n = mxGetM (pargin [1]) ;

    /* get the permutation P */
    Px = mxGetPr (pargin [0]) ;
    P = (int *) mxMalloc (n * sizeof (int)) ;
    for (k = 0 ; k < n ; k++)
    {
	P [k] = (int) (Px [k]) - 1 ;
    }

    /* get the right-hand-side b */
    B = mxGetPr (pargin [1]) ;
    nrhs = mxGetN (pargin [1]) ;

    /* create the solution, x */
    pargout [0] = mxCreateDoubleMatrix (n, nrhs, mxREAL) ;
    X = mxGetPr (pargout [0]) ;

    /* ---------------------------------------------------------------------- */
    /* X = P*B */
    /* ---------------------------------------------------------------------- */

    klu_permute (n, P, n, nrhs, B, X) ;

    mxFree (P) ;
}
