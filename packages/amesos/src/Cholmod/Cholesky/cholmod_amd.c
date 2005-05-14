/* ========================================================================== */
/* === Cholesky/cholmod_amd ================================================= */
/* ========================================================================== */

/*
 * CHOLMOD/Cholesky version 0.1, May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Authors: Timothy A. Davis
 * and Sivasankaran Rajamanickam
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* Cholmod interface to the AMD ordering routine.  The input matrix A must
 * be symmetric.  On output, Perm [k] = i if row/column i of A is the kth
 * row/column of P*A*P'.  This corresponds to A(p,p) in MATLAB notation.
 *
 * workspace: Iwork (nrow).  Allocates a temporary copy of A+A' (with both
 * upper and lower triangular parts) as input to AMD.
 *
 * Requires AMD v1.2 or later.  It can use v1.1, but the prototype for amd_2
 * is not in the amd.h include file in that version, but in amd_internal.h
 * instead.
 */

#include "amd.h"
#include "cholmod_cholesky.h"
#include "cholmod_internal.h"

int cholmod_amd
(
    /* input, not modified on output: */
    cholmod_sparse *A,

    /* output, not defined on input: */
    void *Perm,	    /* size n, output permutation */

    cholmod_common *Common
)
{
    double *W ;
    int *Cp, *Len, *Nv, *Head, *Elen, *Degree, *Wi, *Iwork, *Work4n, *Next ;
    cholmod_sparse *C ;
    int anz, j, n, wsize ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (Perm, FALSE) ;
    Common->status = CHOLMOD_OK ;

    if (A->stype == 0 || A->nrow != A->ncol)
    {
	/* AMD can only operate on symmetric matrices */
	cholmod_error (CHOLMOD_INVALID, "cholmod_amd: matrix must be symmetric",
		Common) ;
	return (FALSE) ;
    }
    PRINT1 (("cholmod_amd sorted %d\n", A->sorted)) ;
    n = A->nrow ;
    if (n == 0)
    {
	/* nothing to do */
	return (TRUE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    /* Note: this is less than the space used in cholmod_analyze, so if
     * cholmod_amd is being called by that routine, no space will be
     * allocated.
     */

    wsize = cholmod_allocate_work (n, 2*n, 4*n, sizeof (int), Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;
    }

    W = Common->Xwork ;
    Work4n = Common->Xwork ;
    Iwork = Common->Iwork ;

    Head = Common->Head ;   /* size n+1, but only n is used */

    Len  = Work4n ;	    /* size n */
    Nv   = Work4n + n ;	    /* size n */
    Next = Work4n + 2*n ;   /* size n */
    Elen = Work4n + 3*n ;   /* size n */

    Degree = Iwork ;	    /* size n */
    Wi = Iwork + n ;	    /* size n */

    /* ---------------------------------------------------------------------- */
    /* construct the input matrix (A+A') for AMD */
    /* ---------------------------------------------------------------------- */

    anz = cholmod_nnz (A, Common) ;

    /* B = A+A', but use only the upper triangular part of A if A->stype = 1
     * and only the lower part of A if A->stype = -1 */
    C = cholmod_copy_sym_to_unsym (A, anz/2 + n, -1, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;
    }

    Cp = C->p ;
    for (j = 0 ; j < n ; j++)
    {
	Len [j] = Cp [j+1] - Cp [j] ;
    }

    /* ---------------------------------------------------------------------- */
    /* order C using AMD with default parameters */
    /* ---------------------------------------------------------------------- */

    amd_2 (n, Cp, C->i, Len, C->nzmax, Cp [n], Nv, Next, Perm, Head, Elen,
	    Degree, Wi, NULL, NULL) ;

    /* ---------------------------------------------------------------------- */
    /* free the AMD workspace and clear the persistent workspace in Common */
    /* ---------------------------------------------------------------------- */

    cholmod_free_sparse (&C, Common) ;
    ASSERT (IMPLIES (Common->status == CHOLMOD_OK,
		cholmod_dump_perm (Perm, n, n, "AMD2 perm", Common))) ;
    for (j = 0 ; j <= n ; j++)
    {
	Head [j] = EMPTY ;
    }
    for (j = 0 ; j < wsize ; j++)
    {
	W [j] = 0. ;
    }
    return (TRUE) ;
}
