/* ========================================================================== */
/* === Cholesky/cholmod_factorize =========================================== */
/* ========================================================================== */

/*
 * CHOLMOD/Cholesky version 0.1, May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Authors: Timothy A. Davis
 * and Sivasankaran Rajamanickam
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* Computes the numerical factorization of a symmetric matrix.  The primary
 * inputs to this routine is a sparse matrix A and the symbolic factor L from
 * cholmod_analyze or a prior numerical factor L.  If A is symmetric, this
 * routine factorizes A(p,p)+beta*I (beta can be zero), where p is the
 * fill-reducing permutation (L->Perm).  If A is unsymmetric, either
 * A(p,:)*A(p,:)'+beta*I or A(p,f)*A(p,f)'+beta*I is factorized.  The set f and
 * the nonzero pattern of the matrix A must be the same as the matrix passed to
 * cholmod_analyze for the supernodal case.  For the simplicial case, it can
 * be different, but it should be the same for best performance.
 *
 * A simplicial factorization or supernodal factorization is chosen, based on
 * the type of the factor L.  If CHOLMOD_SYMBOLIC_SUPER or CHOLMOD_LL_SUPER, a
 * supernodal factorization is computed.  If L is of type CHOLMOD_SYMBOLIC,
 * CHOLMOD_LDL_UNPACKED, or CHOLMOD_LDL_DYNAMIC, a simplicial numeric
 * factorization is computed.  A factor of type CHOLMOD_LL_PACKED cannot be
 * factorized.
 *
 * Once the factorization is complete, it can be left as is or optionally
 * converted into any simplicial numeric type, depending on the
 * Common->final_ftype parameter.  If converted from a supernodal to simplicial
 * type, and the Common->final_resymbol parameter is true, then numerically
 * zero entries in L due to relaxed supernodal amalgamation are removed from
 * the simplicial factor (they are always left in the supernodal form of L).
 * Entries that are numerically zero but present in the simplicial symbolic
 * pattern of L are left in place (that is, the graph of L remains chordal).
 * This is required for the update/downdate/rowadd/rowdell routines to work
 * properly.
 *
 * workspace: Flag (nrow), Head (nrow+1),
 *	if symmetric:   Iwork (2*nrow)
 *	if unsymmetric: Iwork (2*nrow+ncol)
 *	if simplicial: W (nrow).
 *	Allocates up to two temporary copies of its input matrix (including
 *	both pattern and numerical values).
 *
 * If the matrix is not positive definite the routine returns TRUE, but
 * sets Common->status to CHOLMOD_NOT_POSDEF and L->minor is set to the
 * column at which the failure occurred.  All columns in the supernode
 * containing column L->minor are set to zero, as are all columns in all
 * subsequent supernodes.
 */

#include "cholmod_cholesky.h"
#include "cholmod_internal.h"

#ifndef NSUPERNODAL
#include "cholmod_supernodal.h"
#endif

/* ========================================================================== */
/* === cholmod_factorize ==================================================== */
/* ========================================================================== */

int cholmod_factorize
(
    cholmod_sparse *A,
    cholmod_factor *L,
    cholmod_common *Common
)
{
    cholmod_scalar zero ;
    zero.x = 0 ;
    zero.z = 0 ;
    return (cholmod_factorize_p (A, zero, NULL, 0, L, Common)) ;
}


/* ========================================================================== */
/* === cholmod_factorize_p ================================================== */
/* ========================================================================== */

int cholmod_factorize_p	    /* returns TRUE if successful, FALSE otherwise */
(
    /* inputs, not modified */
    cholmod_sparse *A,	    /* nrow-by-ncol, stored in column form */
    cholmod_scalar beta,    /* add beta to diagonal of matrix to factorize */
    void *fset,		    /* if A unsymmetric, factoriize A(p,f)*A(p,f)' */
    size_t fsize,	    /* size of the set f */

    /* input/output: */
    cholmod_factor *L,

    cholmod_common *Common
)
{
    cholmod_sparse *S, *F, *A1, *A2 ;
    size_t maxcsize, wsize, maxrank ;
    int nrow, ncol, symmetry, symmetric, supernodal, convert, nf, n, ok ;
    nf = fsize ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    nrow = A->nrow ;
    ncol = A->ncol ;
    n = L->n ;
    symmetry = A->stype ;
    symmetric = (symmetry != 0) ;
    if (L->n != A->nrow)
    {
	cholmod_error (CHOLMOD_INVALID, "A and L dimensions do not match",
		Common) ;
	return (FALSE) ;
    }
    if (symmetric && nrow != ncol)
    {
	cholmod_error (CHOLMOD_INVALID, "matrix invalid", Common) ;
    }

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    supernodal = (L->ftype == CHOLMOD_LL_SUPER
	       || L->ftype == CHOLMOD_SYMBOLIC_SUPER) ;

    /* determine Common->Xwork usage */
    maxrank = cholmod_maxrank (nrow, Common) ;
    if (supernodal)
    {
	/* cholmod_super_numeric uses Xwork if its required workspace C is
	 * smaller in size than maxrank*n */
	maxcsize = L->maxcsize ;
	wsize = (maxcsize <= (maxrank * nrow)) ? maxcsize : 0 ;
    }
    else
    {
	/* cholmod_rowfac uses one column in Xwork */
	wsize = nrow ;
    }

    cholmod_allocate_work (nrow, 2*nrow + (symmetric ? 0 : ncol), wsize,
	    sizeof (double), Common) ;
    if (Common->status < CHOLMOD_OK || maxrank == 0)
    {
	return (FALSE) ;
    }

    S  = NULL ;
    F  = NULL ;
    A1 = NULL ;
    A2 = NULL ;

    /* convert to a simplicial numeric factorization when done, if requested */
    convert = (Common->final_ftype >= CHOLMOD_LDL_PACKED
	    && Common->final_ftype <= CHOLMOD_LL_PACKED) ;

    /* ---------------------------------------------------------------------- */
    /* perform supernodal LL' or simplicial LDL' factorization */
    /* ---------------------------------------------------------------------- */

    if (supernodal)
    {

#ifndef NSUPERNODAL

	/* ------------------------------------------------------------------ */
	/* supernodal factorization */
	/* ------------------------------------------------------------------ */

	if (L->ordering == CHOLMOD_NATURAL)
	{

	    /* -------------------------------------------------------------- */
	    /* natural ordering */
	    /* -------------------------------------------------------------- */

	    if (symmetry > 0)
	    {
		/* S = tril (A'), F not needed */
		/* workspace: Iwork (nrow) */
		A1 = cholmod_transpose (A, TRUE, NULL, NULL, 0, Common) ;
		S = A1 ;
	    }
	    else if (symmetry < 0)
	    {
		/* This is the fastest option for the natural ordering */
		/* S = A; F not needed */
		S = A ;
	    }
	    else
	    {
		/* F = A(:,f)' */
		/* workspace: Iwork (nrow) */
		/* workspace: Iwork (nrow if no fset; MAX (nrow,ncol) if fset)*/
		A1 = cholmod_transpose (A, TRUE, NULL, fset, nf, Common) ;
		F = A1 ;
		/* S = A */
		S = A ;
	    }

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* permute the input matrix before factorization */
	    /* -------------------------------------------------------------- */

	    if (symmetry > 0)
	    {
		/* This is the fastest option for factoring a permuted matrix */
		/* S = tril (PAP'); F not needed */
		/* workspace: Iwork (2*nrow) */
		A1 = cholmod_transpose (A, TRUE, L->Perm, NULL, 0, Common) ;
		S = A1 ;
	    }
	    else if (symmetry < 0)
	    {
		/* A2 = triu (PAP') */
		/* workspace: Iwork (2*nrow) */
		A2 = cholmod_transpose (A, TRUE, L->Perm, NULL, 0, Common) ;
		/* S = tril (A2'); F not needed */
		/* workspace: Iwork (nrow) */
		A1 = cholmod_transpose (A2, TRUE, NULL, NULL, 0, Common) ;
		S = A1 ;
		cholmod_free_sparse (&A2, Common) ;
		ASSERT (A2 == NULL) ;
	    }
	    else
	    {
		/* F = A(p,f)' */
		/* workspace: Iwork (nrow if no fset; MAX (nrow,ncol) if fset)*/
		A1 = cholmod_transpose (A, TRUE, L->Perm, fset, nf, Common) ;
		F = A1 ;
		/* S = F' */
		/* workspace: Iwork (nrow) */
		A2 = cholmod_transpose (F, TRUE, NULL, NULL, 0, Common) ;
		S = A2 ;
	    }
	}

	ok = (Common->status >= CHOLMOD_OK) ;

	/* ------------------------------------------------------------------ */
	/* supernodal factorization */
	/* ------------------------------------------------------------------ */

	/* workspace: Flag (nrow), Head (nrow+1), Iwork (2*nrow) */
	if (ok)
	{
	    ok = cholmod_super_numeric (S, F, beta, L, Common) ;
	}
	ASSERT (IMPLIES (ok, L->ftype == CHOLMOD_LL_SUPER)) ;

	/* ------------------------------------------------------------------ */
	/* convert to simplicial L if requested */
	/* ------------------------------------------------------------------ */

	if (ok && convert)
	{
	    /* workspace: none */
	    ok = cholmod_change_ftype (L, Common->final_ftype, Common) ;
	    if (ok && Common->final_resymbol)
	    {
		/* workspace: Flag (nrow), Head (nrow+1),
		 *	if symmetric:   Iwork (2*nrow)
		 *	if unsymmetric: Iwork (2*nrow+ncol) */
		ok = cholmod_resymbol_noperm (S, fset, nf, L,
		    Common->final_ftype == CHOLMOD_LDL_PACKED ||
		    Common->final_ftype == CHOLMOD_LL_PACKED, Common) ;
	    }
	}

#else

	/* ------------------------------------------------------------------ */
	/* CHOLMOD Supernodal module not installed */
	/* ------------------------------------------------------------------ */

	ok = FALSE ;
	cholmod_error (CHOLMOD_NOT_INSTALLED,
		"cholmod_factorize: Supernodal module not installed", Common) ;

#endif

    }
    else if (L->ftype == CHOLMOD_SYMBOLIC
	  || L->ftype == CHOLMOD_LDL_UNPACKED
	  || L->ftype == CHOLMOD_LDL_DYNAMIC)
    {

	/* ------------------------------------------------------------------ */
	/* simplicial LDL' factorization */
	/* ------------------------------------------------------------------ */

	/* Permute the input matrix A if necessary.  cholmod_rowfac requires
	 * triu(A) in column form for the symmetric case, and A in column form
	 * for the unsymmetric case (the matrix S).  The unsymmetric case
	 * requires A in row form, or equivalently A' in column form (the
	 * matrix F).
	 */

	if (L->ordering == CHOLMOD_NATURAL)
	{

	    /* -------------------------------------------------------------- */
	    /* natural ordering */
	    /* -------------------------------------------------------------- */

	    if (symmetry > 0)
	    {
		/* F is not needed, S = A */
		S = A ;
	    }
	    else if (symmetry < 0)
	    {
		/* F is not needed, S = A' */
		/* workspace: Iwork (nrow) */
		A2 = cholmod_transpose (A, TRUE, NULL, NULL, 0, Common) ;
		S = A2 ;
	    }
	    else
	    {
		/* F = A (:,f)' */
		/* workspace: Iwork (nrow if no fset; MAX (nrow,ncol) if fset)*/
		A1 = cholmod_transpose (A, TRUE, NULL, fset, nf, Common) ;
		F = A1 ;
		S = A ;
	    }

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* permute the input matrix before factorization */
	    /* -------------------------------------------------------------- */

	    if (symmetry > 0)
	    {
		/* F = tril (A (p,p)') */
		/* workspace: Iwork (2*nrow) */
		A1 = cholmod_transpose (A, TRUE, L->Perm, NULL, 0, Common) ;
		/* A2 = triu (F') */
		/* workspace: Iwork (nrow) */
		A2 = cholmod_transpose (A1, TRUE, NULL, NULL, 0, Common) ;
		/* the symmetric case does not need F, free it and set to NULL*/
		cholmod_free_sparse (&A1, Common) ;
	    }
	    else if (symmetry < 0)
	    {
		/* A2 = triu (A (p,p)'), F not needed.  This is the fastest
		 * way to factorize a matrix using the simplicial routine
		 * (cholmod_rowfac). */
		/* workspace: Iwork (2*nrow) */
		A2 = cholmod_transpose (A, TRUE, L->Perm, NULL, 0, Common) ;
	    }
	    else
	    {
		/* F = A (p,f)' */
		/* workspace: Iwork (nrow if no fset; MAX (nrow,ncol) if fset)*/
		A1 = cholmod_transpose (A, TRUE, L->Perm, fset, nf, Common) ;
		F = A1 ;
		/* A2 = F' */
		/* workspace: Iwork (nrow) */
		A2 = cholmod_transpose (F, TRUE, NULL, NULL, 0, Common) ;
	    }
	    S = A2 ;
	}

	ok = (Common->status >= CHOLMOD_OK) ;

	/* ------------------------------------------------------------------ */
	/* simplicial LDL' factorization */
	/* ------------------------------------------------------------------ */

	/* factorize beta*I+S (symmetric) or beta*I+F*F' (unsymmetric) */
	/* workspace: Flag (nrow), W (nrow), Iwork (2*nrow) */
	if (ok)
	{
	    ok = cholmod_rowfac (S, F, beta, 0, nrow, L, Common) ;
	}

	/* convert to final form, if requested */
	if (ok && convert)
	{
	    /* workspace: none */
	    ok = cholmod_change_ftype (L, Common->final_ftype, Common) ;
	}

    }
    else
    {
	/* L invalid.  Cannot factorize a simplicial LL' */
	ok = FALSE ;
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_factorize: Cannot factorize non-supernodal LL'",
		Common) ;
    }

    /* ---------------------------------------------------------------------- */
    /* free A1 and A2 if they exist */
    /* ---------------------------------------------------------------------- */

    cholmod_free_sparse (&A1, Common) ;
    cholmod_free_sparse (&A2, Common) ;
    return (ok) ;
}
