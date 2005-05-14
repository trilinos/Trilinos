/* ========================================================================== */
/* === ccolamd mexFunction ================================================== */
/* ========================================================================== */

/*
 * CCOLAMD version 0.1, May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Authors: Timothy A. Davis
 * and Sivasankaran Rajamanickam
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* TODO: update comments, here and throughout the code (still refers to colamd).

    Usage:

	P = colamd (A) ;

	P = colamd (A, knobs) ;

	[ P, stats ] = colamd (A) ;

	[ P, stats ] = colamd (A, knobs) ;

    Returns a permutation vector P such that the LU factorization of A (:,P)
    tends to be sparser than that of A.  The Cholesky factorization of
    (A (:,P))'*(A (:,P)) will also tend to be sparser than that of A'*A.
    This routine provides the same functionality as COLMMD, but is much faster
    and returns a better permutation vector.  Note that the COLMMD m-file in
    MATLAB 5.2 also performs a column elimination tree post-ordering.  This
    mexFunction does not do this post-ordering.  This mexFunction is a
    replacement for the p = sparsfun ('colmmd', A) operation.

    The knobs and stats vectors are optional:

	knobs (1)	rows with more than (knobs (1))*n_col entries
			are removed prior to ordering.  If knobs is not present,
			then the default is used (0.5).

	knobs (2)	columns with more than (knobs (2))*n_row entries
			are removed prior to ordering, and placed last in the
			column permutation.  If knobs is not present,
			then the default is used (0.5).

	knobs (3)	print level, similar to spparms ('spumoni')

	stats (1)	the number of dense (or empty) rows ignored

	stats (2)	the number of dense (or empty) columms.  These
			are ordered last, in their natural order.

	stats (3)	the number of garbage collections performed.

	stats (4)	return status:

			0:  matrix is a valid MATLAB matrix.

			1:  matrix has duplicate entries or unsorted columns.
			    This should not occur in a valid MATLAB matrix,
			    but the ordering proceeded anyway by sorting the
			    row indices in each column and by ignoring the
			    duplicate row indices in each column.  See
			    stats (5:7) for more information.

	stats (5)	highest numbered column that is unsorted or has
			duplicate entries (zero if none)

	stats (6)	last seen duplicate or unsorted row index
			(zero if none)

	stats (7)	number of duplicate or unsorted row indices

*/

/* ========================================================================== */
/* === Include files ======================================================== */
/* ========================================================================== */

#include "ccolamd.h"
#include "mex.h"
#include "matrix.h"
#include <stdlib.h>
#include <string.h>

/* ========================================================================== */
/* === colamd mexFunction =================================================== */
/* ========================================================================== */

void mexFunction
(
    /* === Parameters ======================================================= */

    int nlhs,			/* number of left-hand sides */
    mxArray *plhs [],		/* left-hand side matrices */
    int nrhs,			/* number of right--hand sides */
    const mxArray *prhs []	/* right-hand side matrices */
)
{
    /* === Local variables ================================================== */

    int *A ;			/* colamd's copy of the matrix, and workspace */
    int *cset;			/* colamd's copy of the constraint set */
    double *dcs;		/* colamd's copy of the constraint set */
    int *p ;			/* colamd's copy of the column pointers */
    int Alen ;			/* size of A */
    int cslen ;			/* size of CS  */
    int n_col ;			/* number of columns of A */
    int n_row ;			/* number of rows of A */
    int nnz ;			/* number of entries in A */
    int full ;			/* TRUE if input matrix full, FALSE if sparse */
    double knobs [CCOLAMD_KNOBS] ; /* colamd user-controllable parameters */
    double *out_perm ;		/* output permutation vector */
    double *out_stats ;		/* output stats vector */
    double *in_knobs ;		/* input knobs vector */
    int i ;			/* loop counter */
    mxArray *Ainput ;		/* input matrix handle */
    int spumoni ;		/* verbosity variable */
    int stats [CCOLAMD_STATS] ;	/* stats for colamd */

    /* === Check inputs ===================================================== */

    if (nrhs < 1 || nrhs > 3 || nlhs < 0 || nlhs > 2)
    {
	mexErrMsgTxt (
	"colamd: incorrect number of input and/or output arguments") ;
    }

    if (nrhs >= 2)
    {
	dcs =  mxGetPr (prhs[1]);
	cslen = mxGetNumberOfElements( prhs[1] ) ;

	cset = (int *) mxCalloc (cslen, sizeof (int)) ;

	for (i=0 ; i < cslen ; i++ )
	     cset[i] = ((int)dcs[i]-1) ; /* Make the cset zero based */
    }

    /* === Get knobs ======================================================== */

    ccolamd_set_defaults (knobs) ;
    spumoni = 0 ;

    /* check for user-passed knobs */
    if (nrhs == 3)
    {
	in_knobs = mxGetPr (prhs [2]) ;
	i = mxGetNumberOfElements (prhs [2]) ;
	if (i > 0) knobs [CCOLAMD_DENSE_ROW] = in_knobs [CCOLAMD_DENSE_ROW] ;
	if (i > 1) knobs [CCOLAMD_DENSE_COL] = in_knobs [CCOLAMD_DENSE_COL] ;
	if (i > 2) spumoni = (int) in_knobs [2] ;
    }

    /* print knob settings if spumoni is set */
    if (spumoni > 0)
    {
	mexPrintf ("colamd: dense row fraction: %f\n",
	    knobs [CCOLAMD_DENSE_ROW]) ;
	mexPrintf ("colamd: dense col fraction: %f\n",
	    knobs [CCOLAMD_DENSE_COL]) ;
    }

    /* === If A is full, convert to a sparse matrix ========================= */

    Ainput = (mxArray *) prhs [0] ;
    if (mxGetNumberOfDimensions (Ainput) != 2)
    {
	mexErrMsgTxt ("colamd: input matrix must be 2-dimensional") ;
    }
    full = !mxIsSparse (Ainput) ;
    if (full)
    {
	mexCallMATLAB (1, &Ainput, 1, (mxArray **) prhs, "sparse") ;
    }

    /* === Allocate workspace for colamd ==================================== */

    /* get size of matrix */
    n_row = mxGetM (Ainput) ;
    n_col = mxGetN (Ainput) ;

    /* get column pointer vector so we can find nnz */
    p = (int *) mxCalloc (n_col+1, sizeof (int)) ;
    (void) memcpy (p, mxGetJc (Ainput), (n_col+1)*sizeof (int)) ;
    nnz = p [n_col] ;
    Alen = ccolamd_recommended (nnz, n_row, n_col) ;

    /* === Copy input matrix into workspace ================================= */

    A = (int *) mxCalloc (Alen, sizeof (int)) ;
    (void) memcpy (A, mxGetIr (Ainput), nnz*sizeof (int)) ;

    if (full)
    {
	mxDestroyArray (Ainput) ;
    }

    /* Check constraint set size */
    if ( cslen != n_col )
    {
    	mexErrMsgTxt ("colamd : Constraint set length != no of columns !") ;
    }

    /* === Order the columns (destroys A) =================================== */

    if (!ccolamd (n_row, n_col, Alen, A, p, knobs, stats, cset))
    {
	ccolamd_report (stats) ;
	mexErrMsgTxt ("colamd error!") ;
    }
    mxFree (A) ;
    mxFree (cset) ;

    /* === Return the permutation vector ==================================== */

    plhs [0] = mxCreateDoubleMatrix (1, n_col, mxREAL) ;
    out_perm = mxGetPr (plhs [0]) ;
    for (i = 0 ; i < n_col ; i++)
    {
	/* colamd is 0-based, but MATLAB expects this to be 1-based */
	out_perm [i] = p [i] + 1 ;
    }
    mxFree (p) ;

    /* === Return the stats vector ========================================== */

    /* print stats if spumoni > 0 */
    if (spumoni > 0)
    {
	ccolamd_report (stats) ;
    }

    if (nlhs == 2)
    {
	plhs [1] = mxCreateDoubleMatrix (1, CCOLAMD_STATS, mxREAL) ;
	out_stats = mxGetPr (plhs [1]) ;
	for (i = 0 ; i < CCOLAMD_STATS ; i++)
	{
	    out_stats [i] = stats [i] ;
	}

	/* fix stats (5) and (6), for 1-based information on jumbled matrix. */
	/* note that this correction doesn't occur if symamd returns FALSE */
	out_stats [CCOLAMD_INFO1] ++ ; 
	out_stats [CCOLAMD_INFO2] ++ ; 
    }
}
