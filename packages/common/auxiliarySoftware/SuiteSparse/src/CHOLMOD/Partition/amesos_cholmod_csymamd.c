/* ========================================================================== */
/* === Partition/cholmod_csymamd ============================================ */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Partition Module.
 * Copyright (C) 2005-2006, Univ. of Florida.  Author: Timothy A. Davis
 * The CHOLMOD/Partition Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.cise.ufl.edu/research/sparse
 * -------------------------------------------------------------------------- */

/* CHOLMOD interface to the CSYMAMD ordering routine.  Finds a permutation
 * p such that the Cholesky factorization of PAP' is sparser than A.
 * The column etree is found and postordered, and the CSYMAMD
 * ordering is then combined with its postordering.  If A is unsymmetric,
 * A+A' is ordered (A must be square).
 *
 * workspace: Head (nrow+1)
 *
 * Supports any xtype (pattern, real, complex, or zomplex).
 */

#ifndef NPARTITION

#include "amesos_cholmod_internal.h"
#include "trilinos_ccolamd.h"
#include "amesos_cholmod_partition.h"

#if (TRILINOS_CCOLAMD_VERSION < TRILINOS_CCOLAMD_VERSION_CODE (2,5))
#error "CCOLAMD v2.0 or later is required"
#endif

/* ========================================================================== */
/* === cholmod_csymamd ====================================================== */
/* ========================================================================== */

int CHOLMOD(csymamd)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to order */
    /* ---- output --- */
    Int *Cmember,	/* size nrow.  see cholmod_ccolamd.c for description */
    Int *Perm,		/* size A->nrow, output permutation */
    /* --------------- */
    cholmod_common *Common
)
{
    double knobs [TRILINOS_CCOLAMD_KNOBS] ;
    Int *perm, *Head ;
    Int ok, i, nrow, stats [TRILINOS_CCOLAMD_STATS] ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (Perm, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    Common->status = CHOLMOD_OK ;

    if (A->nrow != A->ncol || !(A->packed))
    {
	ERROR (CHOLMOD_INVALID, "matrix must be square and packed") ;
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    nrow = A->nrow ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    CHOLMOD(allocate_work) (nrow, 0, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* order the matrix (does not affect A->p or A->i) */
    /* ---------------------------------------------------------------------- */

    perm = Common->Head ;	/* size nrow+1 (i/l/l) */

    /* get parameters */
#ifdef LONG
    trilinos_ccolamd_l_set_defaults (knobs) ;
#else
    trilinos_ccolamd_set_defaults (knobs) ;
#endif
    if (Common->current >= 0 && Common->current < CHOLMOD_MAXMETHODS)
    {
	/* get the knobs from the Common parameters */
	knobs [TRILINOS_CCOLAMD_DENSE_ROW] =Common->method[Common->current].prune_dense ;
	knobs [TRILINOS_CCOLAMD_AGGRESSIVE]=Common->method[Common->current].aggressive ;
    }

    {
#ifdef LONG
	trilinos_csymamd_l (nrow, A->i, A->p, perm, knobs, stats, Common->calloc_memory,
		Common->free_memory, Cmember, A->stype) ;
#else
	trilinos_csymamd (nrow, A->i, A->p, perm, knobs, stats, Common->calloc_memory,
		Common->free_memory, Cmember, A->stype) ;
#endif
	ok = stats [TRILINOS_CCOLAMD_STATUS] ;
    }

    if (ok == TRILINOS_CCOLAMD_ERROR_out_of_memory)
    {
	ERROR (CHOLMOD_OUT_OF_MEMORY, "out of memory") ; 
    }
    ok = (ok == TRILINOS_CCOLAMD_OK || ok == TRILINOS_CCOLAMD_OK_BUT_JUMBLED) ;

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
	Head [i] = TRILINOS_CHOLMOD_EMPTY ;
    }

    return (ok) ;
}
#endif
