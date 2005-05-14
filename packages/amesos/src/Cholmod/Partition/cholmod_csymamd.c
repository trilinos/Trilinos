/* ========================================================================== */
/* === Partition/cholmod_csymamd ============================================ */
/* ========================================================================== */

/*
 * CHOLMOD/Partition version 0.1.  May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* Cholmod interface to the csymamd ordering routine.  Finds a permutation
 * p such that the Cholesky factorization of PAP' is sparser than A.
 * The column etree is found and postordered, and the csymamd
 * ordering is then combined with its postordering.  If A is unsymmetric,
 * A+A' is ordered (A must be square).
 *
 * workspace: Head (nrow+1)
 *
 * TODO check int overflow case in csymamd
 */

#include "cholmod_partition.h"
#include "cholmod_internal.h"
#include "ccolamd.h"

int cholmod_csymamd
(
    /* inputs, not modified */
    cholmod_sparse *A,	    /* user's input matrix, nrow-by-ncol */
    void *Cmember,	    /* size nrow.  Cmember [i] = c if row i is in the
			     * constraint set c.  c must be >= 0.  The # of
	* constraint sets is max (Cmember) + 1.  If Cmember is NULL, then it
	* is interpretted as Cmember [i] = 0 for all i */

    /* output, not defined on input */
    void *Perm_p,	    /* size nrow, output permutation */

    cholmod_common *Common
)
{
    double knobs [CCOLAMD_KNOBS] ;
    int *perm, *Perm, *Ap, *Ai, *Head ;
    int ok, i, nrow, stats [CCOLAMD_STATS] ;

    ASSERT (cholmod_dump_sparse (A, "csymamd", Common) >= 0) ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    Perm = Perm_p ;
    RETURN_IF_NULL (Perm, FALSE) ;
    Common->status = CHOLMOD_OK ;

    if (A->nrow != A->ncol || !(A->packed))
    {
	cholmod_error (CHOLMOD_INVALID,
	    "cholmod_csymamd: matrix must be square and packed", Common) ;
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    nrow = A->nrow ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    cholmod_allocate_work (nrow, 0, 0, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* order the matrix (does not affect A->p or A->i) */
    /* ---------------------------------------------------------------------- */

    perm = Common->Head ;	/* size nrow+1 (i/l/l) */

    ccolamd_set_defaults (knobs) ;
    /* TODO: fix this knob, 0: Cholesky, otherwise LU */
    knobs [CCOLAMD_FACT_TYPE] = 2 ; /* do permutation for Cholesky */

    Ap = A->p ;
    Ai = A->i ;
    csymamd (nrow, Ai, Ap, perm, knobs, stats,
		Common->calloc_memory, Common->free_memory, Cmember,
		A->stype) ;

    ok = stats [CCOLAMD_STATUS] ;
    if (ok == CCOLAMD_ERROR_out_of_memory)
    {
	cholmod_error (CHOLMOD_OUT_OF_MEMORY, "out of memory", Common) ; 
    }
    ok = (ok == CCOLAMD_OK || ok == CCOLAMD_OK_BUT_JUMBLED) ;

    /* ---------------------------------------------------------------------- */
    /* free the workspace and return result */
    /* ---------------------------------------------------------------------- */

    /* permutation returned in perm [0..n-1] */
    for (i = 0 ; i < nrow ; i++)
    {
	Perm [i] = perm [i] ;
    }

    /* clear Head workspace (used for perm, in csymamd): */
    Head = Common->Head ;
    for (i = 0 ; i <= nrow ; i++)
    {
	Head [i] = EMPTY ;
    }

    return (ok) ;
}
