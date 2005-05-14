/* ========================================================================== */
/* === Partition/cholmod_ccolamd ============================================ */
/* ========================================================================== */

/*
 * CHOLMOD/Partition version 0.1.  May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* Cholmod interface to the CCOLAMD ordering routine.  Finds a permutation
 * p such that the Cholesky factorization of PAA'P' is sparser than AA'.
 * The column etree is found and postordered, and the ccolamd ordering is then
 * combined with its postordering.  A must be unsymmetric.
 *
 * workspace: Iwork (MAX (nrow,ncol))
 *	Allocates a copy of its input matrix, which is
 *	then used as CCOLAMD's workspace.
 *
 * TODO: make sure ccolamd handles int overflow (alen, etc)
 */

#include "cholmod_partition.h"
#include "cholmod_internal.h"
#include "ccolamd.h"

/* ========================================================================== */
/* === ccolamd_interface ==================================================== */
/* ========================================================================== */

/* Order with ccolamd */

static int ccolamd_interface
(
    cholmod_sparse *A,
    int alen,
    int *Perm,
    void *Cmember,
    int *fset,
    int fsize,
    cholmod_sparse *C,
    cholmod_common *Common
)
{
    int *Cp = NULL, *Ci = NULL ;
    int ok, k, nrow, ncol, stats [CCOLAMD_STATS] ;

    nrow = A->nrow ;
    ncol = A->ncol ;

    /* ---------------------------------------------------------------------- */
    /* copy (and transpose) the input matrix A into the ccolamd workspace */
    /* ---------------------------------------------------------------------- */

    /* C = A (:,f)', which also packs A if needed. */
    /* workspace: Iwork (nrow if no fset; MAX (nrow,ncol) if fset non-NULL) */
    ok = cholmod_transpose_unsym (A, FALSE, NULL, fset, fsize, C, Common) ;

    /* ---------------------------------------------------------------------- */
    /* order the matrix (destroys the contents of C->i and C->p) */
    /* ---------------------------------------------------------------------- */

    if (ok)
    {
	Cp = C->p ;
	Ci = C->i ;
	ccolamd (ncol, nrow, alen, Ci, Cp, NULL, stats, Cmember) ;
	ok = stats [CCOLAMD_STATUS] ;
	ok = (ok == CCOLAMD_OK || ok == CCOLAMD_OK_BUT_JUMBLED) ;
	/* permutation returned in C->p, if the ordering succeeded */
	for (k = 0 ; k < nrow ; k++)
	{
	    Perm [k] = Cp [k] ;
	}
    }

    return (ok) ;
}


/* ========================================================================== */
/* === cholmod_ccolamd ====================================================== */
/* ========================================================================== */

int cholmod_ccolamd
(
    /* inputs, not modified */
    cholmod_sparse *A,	    /* user's input matrix, nrow-by-ncol */
    void *fset,
    size_t fsize,

    void *Cmember,	    /* size nrow.  Cmember [i] = c if row i is in the
			     * constraint set c.  c must be >= 0.  The # of
			     * constraint sets is max (Cmember) + 1.  If
			     * Cmember is NULL, then it is interpretted as
			     * Cmember [i] = 0 for all i */

    /* output, not defined on input */
    void *Perm,		    /* size nrow, output permutation */

    cholmod_common *Common
)
{
    cholmod_sparse *C ;
    int ok, alen, nrow, ncol ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (Perm, FALSE) ;
    if (A->stype != 0)
    {
	cholmod_error (CHOLMOD_INVALID, "cholmod_ccolamd:"
		"matrix must be unsymmetric", Common) ;
	return (FALSE) ;
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    nrow = A->nrow ;
    ncol = A->ncol ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    cholmod_allocate_work (0, MAX (nrow,ncol), 0, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* allocate workspace for ccolamd */
    /* ---------------------------------------------------------------------- */

    alen = ccolamd_recommended (A->nzmax, ncol, nrow) ;
    if (alen < 0)
    {
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_ccolamd: matrix invalid or too large", Common) ;
	return (FALSE) ;
    }
    C = cholmod_allocate_sparse (ncol, nrow, alen, TRUE, TRUE, 0, FALSE,
	    Common) ;

    /* ---------------------------------------------------------------------- */
    /* order with ccolamd */
    /* ---------------------------------------------------------------------- */

    ok = ccolamd_interface (A, alen, Perm, Cmember, fset, fsize, C, Common) ;

    /* ---------------------------------------------------------------------- */
    /* free the workspace and return result */
    /* ---------------------------------------------------------------------- */

    cholmod_free_sparse (&C, Common) ;
    return (ok) ;
}
