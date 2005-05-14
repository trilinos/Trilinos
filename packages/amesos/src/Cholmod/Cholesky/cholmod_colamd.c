/* ========================================================================== */
/* === Cholesky/cholmod_colamd ============================================== */
/* ========================================================================== */

/*
 * CHOLMOD/Cholesky version 0.1, May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Authors: Timothy A. Davis
 * and Sivasankaran Rajamanickam
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* Cholmod interface to the COLAMD ordering routine (version 2.3 or later).
 * Finds a permutation p such that the Cholesky factorization of PAA'P' is
 * sparser than AA' using colamd.  The column etree is found and postordered,
 * and the colamd ordering is then combined with its postordering.  A must be
 * unsymmetric.
 *
 * There can be no duplicate entries in f.
 * f can be length 0 to n if A is m-by-n.
 *
 * workspace: Iwork (nrow+MAX(nrow,ncol)).
 *	Allocates a copy of its input matrix, which
 *	is then used as CCOLAMD's workspace.
 *
 * TODO: check int overflow case in colamd
 *
 * Note that v2.4 corrected a bug in v2.3, but the bug only affected symamd,
 * not colamd.  This routine can safely use colamd v2.3 or later.
 */

#include "colamd.h"
#include "cholmod_cholesky.h"
#include "cholmod_internal.h"

/* ========================================================================== */
/* === colamd_interface ===================================================== */
/* ========================================================================== */

/* Order with colamd, and then postorder the matrix using the etree of
 * A(p,f)*A(p,f)' */

static int colamd_interface
(
    cholmod_sparse *A,
    int alen,
    int *Perm,
    int *fset,
    int fsize,
    cholmod_sparse *C,
    cholmod_common *Common
)
{
    int *Cp = NULL, *Ci = NULL ;
    int *NewPerm, *Parent, *Post ;
    int ok, k, nrow, ncol, stats [COLAMD_STATS] ;

    nrow = A->nrow ;
    ncol = A->ncol ;

    /* ---------------------------------------------------------------------- */
    /* copy (and transpose) the input matrix A into the colamd workspace */
    /* ---------------------------------------------------------------------- */

    /* C = A (:,f)', which also packs A if needed. */
    /* workspace: Iwork (nrow if no fset; MAX (nrow,ncol) if fset) */
    ok = cholmod_transpose_unsym (A, FALSE, NULL, fset, fsize, C, Common) ;

    /* ---------------------------------------------------------------------- */
    /* order the matrix (destroys the contents of C->i and C->p) */
    /* ---------------------------------------------------------------------- */

    if (ok)
    {
	Cp = C->p ;
	Ci = C->i ;
	colamd (ncol, nrow, alen, Ci, Cp, NULL, stats) ;
	ok = stats [COLAMD_STATUS] ;
	ok = (ok == COLAMD_OK || ok == COLAMD_OK_BUT_JUMBLED) ;
	/* permutation returned in C->p, if the ordering succeeded */
	for (k = 0 ; k < nrow ; k++)
	{
	    Perm [k] = Cp [k] ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* column etree postordering */
    /* ---------------------------------------------------------------------- */

    /* C->i is larger than A->nzmax + 2*nrow, but only A->nzmax is used for
     * the matrix C below.  Use the rest for Parent and Post */
    Parent = NULL ;
    Post = NULL ;
    if (ok)
    {
	Parent = Ci + A->nzmax ;		/* size nrow workspace */
	Post   = Ci + A->nzmax + nrow ;		/* size nrow workspace */
    }

    /* C = A (p,f)' */
    /* workspace: Iwork (nrow if no fset, MAX (nrow,ncol) if fset) */
    ok = ok && cholmod_transpose_unsym (A, FALSE, Perm, fset, fsize, C, Common);

    /* find the column etree of C and postorder it */
    /* workspace: Iwork (nrow+ncol) */
    ok = ok && cholmod_etree (C, Parent, Common) ;

    /* postorder the column etree */
    /* workspace: Iwork (2*nrow) */
    ok = ok && (cholmod_postorder (Parent, nrow, Post, Common) == nrow) ;

    /* combine the colamd permutation with its postordering */
    if (ok)
    {
	NewPerm = Common->Iwork ;		/* size nrow (i/i/l) */
	for (k = 0 ; k < nrow ; k++)
	{
	    NewPerm [k] = Perm [Post [k]] ;
	}
	for (k = 0 ; k < nrow ; k++)
	{
	    Perm [k] = NewPerm [k] ;
	}
    }

    return (ok) ;
}


/* ========================================================================== */
/* === cholmod_colamd ======================================================= */
/* ========================================================================== */

int cholmod_colamd
(
    /* inputs, not modified */
    cholmod_sparse *A,	    /* user's input matrix, nrow-by-ncol */
    void *fset,
    size_t fsize,

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
	cholmod_error (CHOLMOD_INVALID, "cholmod_colamd"
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

    cholmod_allocate_work (0, nrow + MAX (nrow,ncol), 0, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* allocate workspace for colamd */
    /* ---------------------------------------------------------------------- */

    alen = colamd_recommended (A->nzmax, ncol, nrow) ;
    if (alen < 0)
    {
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_colamd: matrix invalid or too large", Common) ;
	return (FALSE) ;
    }
    C = cholmod_allocate_sparse (ncol, nrow, alen, TRUE, TRUE, 0, FALSE,
	    Common) ;

    /* ---------------------------------------------------------------------- */
    /* order with colamd, followed by coletree postorder */
    /* ---------------------------------------------------------------------- */

    ok = colamd_interface (A, alen, Perm, fset, fsize, C, Common) ;

    /* ---------------------------------------------------------------------- */
    /* free the workspace and return result */
    /* ---------------------------------------------------------------------- */

    cholmod_free_sparse (&C, Common) ;
    return (ok) ;
}
