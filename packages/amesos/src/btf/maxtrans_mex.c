/* ========================================================================== */
/* === maxtrans mexFunction ================================================= */
/* ========================================================================== */

/* MAXTRANS: Find a column permutation for a zero-free diagonal.
 *
 * Usage:
 *
 * p = maxtrans (A) ;
 *
 * A (:,p) has a zero-free diagonal unless A is structurally singular.
 * If the matrix is structurally singular, the p will contain zeros.  Similar
 * to p = dmperm (A), except that dmperm returns a row permutation.
 */

/* ========================================================================== */

#include "mex.h"
#include "btf.h"

void mexFunction
(
    int	nargout,
    mxArray *pargout [ ],
    int	nargin,
    const mxArray *pargin [ ]
)
{
    int n, i, *Ap, *Ai, *Match, nfound, *Work ;
    double *Matchx ;

    /* ---------------------------------------------------------------------- */
    /* get inputs and allocate workspace */
    /* ---------------------------------------------------------------------- */

    if (nargin != 1 || nargout > 1)
    {
	mexErrMsgTxt ("Usage: p = maxtrans (A)") ;
    }
    n = mxGetM (pargin [0]) ;
    if (!mxIsSparse (pargin [0]) || n != mxGetN (pargin [0]))
    {
    	mexErrMsgTxt ("maxtrans: A must be sparse, square, and non-empty") ;
    }

    /* get sparse matrix A */
    Ap = mxGetJc (pargin [0]) ;
    Ai = mxGetIr (pargin [0]) ;

    /* get output array */
    Match = mxMalloc (n * sizeof (int)) ;

    /* get workspace of size 5n (recursive version needs only 2n) */
    Work = mxMalloc (5*n * sizeof (int)) ;

    /* ---------------------------------------------------------------------- */
    /* perform the maximum transversal */
    /* ---------------------------------------------------------------------- */

    nfound = maxtrans (n, Ap, Ai, Match, Work) ;
    if (nfound < n)
    {
	printf ("maxtrans: A is singular\n") ;
    }

    /* ---------------------------------------------------------------------- */
    /* create outputs and free workspace */
    /* ---------------------------------------------------------------------- */

    pargout [0] = mxCreateDoubleMatrix (1, n, mxREAL) ;
    Matchx = mxGetPr (pargout [0]) ;
    for (i = 0 ; i < n ; i++)
    {
	Matchx [i] = Match [i] + 1 ;	/* convert to 1-based */
    }

    mxFree (Work) ;
    mxFree (Match) ;
}
