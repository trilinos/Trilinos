/* ========================================================================== */
/* === cholmod_copy ========================================================= */
/* ========================================================================== */

/*
 * CHOLMOD/Core version 0.1. May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* C = A, which allocates C and copies A into C, with possible change of
 * symmetry.  The diagonal can optionally be removed.  The numerical entries
 * can optionally be copied.  This routine differs from cholmod_copy_sparse,
 * which makes an exact copy of a sparse matrix.
 *
 * A can be of any type (packed/unpacked, upper/lower/unsymmetric).  C is
 * packed and can be of any symmetry (upper/lower/unsymmetric), except that if
 * A is rectangular C can only be unsymmetric.  If the symmetry of A and C
 * differ, then the appropriate conversion is made.
 *
 * Symmetry of A (A->stype):
 * <0: lower: assume A is symmetric with just tril(A); the rest of A is ignored
 *  0  unsym: assume A is unsymmetric; consider all entries in A
 * >0  upper: assume A is symmetric with just triu(A); the rest of A is ignored
 *
 * Symmetry of C (symmetric parameter):
 * <0  lower: return just tril(C)
 *  0  unsym: return all of C
 * >0  upper: return just triu(C)
 *
 * In MATLAB:		Using cholmod_copy:
 * ----------		-----------------------------
 * C = A ;		A unsymmetric, C unsymmetric
 *
 * C = tril (A) ;	A unsymmetric, C lower
 *
 * C = triu (A) ;	A unsymmetric, C upper
 *
 * U = triu (A) ;	A upper, C unsymmetric
 * L = tril (U',-1) ;	
 * C = L+U ;
 *
 * C = triu (A)' ;	A upper, C lower
 *
 * C = triu (A) ;	A upper, C upper
 *
 * L = tril (A) ;	A lower, C unsymmetric
 * U = triu (L',1) ;
 * C = L+U ;
 *
 * C = tril (A) ;	A lower, C lower
 *
 * C = tril (A)' ;	A lower, C upper
 *
 * workspace: Iwork (max (nrow,ncol))
 */

#include "cholmod_core.h"
#include "cholmod_internal.h"


/* ========================================================================== */
/* === cholmod_copy_sym_to_unsym ============================================ */
/* ========================================================================== */

/* Construct an unsymmetric copy of a symmetric sparse matrix.  This is the
 * same as C = cholmod_copy (A, 0, mode, Common), except that A must be
 * symmetric (A->stype != 0) and C can be allocated with extra space.
 */

cholmod_sparse *cholmod_copy_sym_to_unsym
(
    /* inputs, not modified on output */
    cholmod_sparse *A,
    size_t extra,	/* extra space to allocate in C */
    int mode,		/* >0: numerical, 0: pattern, <0: pattern (no diag) */
    cholmod_common *Common
)
{
    double aij ;
    double *Ax, *Cx ;
    int *Ap, *Ai, *Anz, *Cp, *Ci, *Wj, *Iwork ;
    cholmod_sparse *C ;
    int nrow, ncol, nz, packed, j, p, pend, i, pc, up, lo, values, diag, asym ;

    DEBUG (int orig) ;
    DEBUG (int nnzdiag) ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    nrow = A->nrow ;
    ncol = A->ncol ;
    if (A->stype == 0 || nrow != ncol)
    {
	/* inputs invalid */
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_copy_sym_to_unsym: matrix invalid", Common);
	return (NULL) ;
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    cholmod_allocate_work (0, MAX (nrow,ncol), 0, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	/* out of memory */
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    DEBUG (orig = Common->malloc_count) ;

    Ap  = A->p ;
    Anz = A->nz ;
    Ai  = A->i ;
    Ax  = A->x ;
    packed = A->packed ;
    values = (mode > 0) && (A->xtype != CHOLMOD_PATTERN) ;
    diag = (mode >= 0) ;

    asym = SIGN (A->stype) ;
    up = (asym > 0) ;
    lo = (asym < 0) ;

    /* ---------------------------------------------------------------------- */
    /* create an unsymmetric copy of a symmetric matrix */
    /* ---------------------------------------------------------------------- */

    Iwork = Common->Iwork ;
    Wj = Iwork ;		    /* size ncol (i/i/l) */

    /* In MATLAB notation, for converting a symmetric/upper matrix:
     *	U = triu (A) ;
     *	L = tril (U',-1) ;
     *	C = L + U ;
     *
     * For converting a symmetric/lower matrix to unsymmetric:
     *	L = tril (A) ;
     *	U = triu (L',1) ;
     *	C = L + U ;
     */
    ASSERT (up || lo) ;
    PRINT1 (("copy: convert symmetric to unsym\n")) ;

    /* count the number of entries in each column of C */
    for (j = 0 ; j < ncol ; j++)
    {
	Wj [j] = 0 ;
    }
    for (j = 0 ; j < ncol ; j++)
    {
	p = Ap [j] ;
	pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	for ( ; p < pend ; p++)
	{
	    i = Ai [p] ;
	    if (i == j)
	    {
		/* the diagonal entry A(i,i) will appear just once
		 * (unless it is excluded with mode < 0) */
		if (diag)
		{
		    Wj [j]++ ;
		}
	    }
	    else if ((up && i < j) || (lo && i > j))
	    {
		/* upper case:  A(i,j) is in the strictly upper part;
		 * A(j,i) will be added to the strictly lower part of C.
		 * lower case is the opposite. */
		Wj [j]++ ;
		Wj [i]++ ;
	    }
	}
    }
    nz = 0 ;
    for (j = 0 ; j < ncol ; j++)
    {
	nz += Wj [j] ;
    }

    /* allocate C.  C is sorted if and only if A is sorted */
    C = cholmod_allocate_sparse (nrow, ncol, nz + extra, A->sorted, TRUE, 0,
	    values, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	ASSERT (Common->malloc_count == orig) ;
	return (NULL) ;
    }

    Cp = C->p ;
    Ci = C->i ;
    Cx = C->x ;

    /* construct the column pointers for C */
    p = 0 ;
    for (j = 0 ; j < ncol ; j++)
    {
	Cp [j] = p ;
	p += Wj [j] ;
    }
    Cp [ncol] = p ;
    for (j = 0 ; j < ncol ; j++)
    {
	Wj [j] = Cp [j] ;
    }

    /* construct C */
    if (values)
    {

	/* pattern and values */
	ASSERT (diag) ;
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		i = Ai [p] ;
		aij = Ax [p] ;
		if (i == j)
		{
		    /* add diagonal entry A(i,i) to column i */
		    pc = Wj [i]++ ;
		    Ci [pc] = i ;
		    Cx [pc] = aij ;
		}
		else if ((up && i < j) || (lo && i > j))
		{
		    /* add A(i,j) to column j */
		    pc = Wj [j]++ ;
		    Ci [pc] = i ;
		    Cx [pc] = aij ;
		    /* add A(j,i) to column i */
		    pc = Wj [i]++ ;
		    Ci [pc] = j ;
		    Cx [pc] = aij ;
		}
	    }
	}

    }
    else
    {

	/* pattern only, possibly excluding the diagonal */
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		i = Ai [p] ;
		if (i == j)
		{
		    /* add diagonal entry A(i,i) to column i
		     * (unless it is excluded with mode < 0) */
		    if (diag)
		    {
			Ci [Wj [i]++] = i ;
		    }
		}
		else if ((up && i < j) || (lo && i > j))
		{
		    /* add A(i,j) to column j */
		    Ci [Wj [j]++] = i ;
		    /* add A(j,i) to column i */
		    Ci [Wj [i]++] = j ;
		}
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* return the result */
    /* ---------------------------------------------------------------------- */

    ASSERT (Common->malloc_count == orig + ((values) ? 4 : 3)) ;
    DEBUG (nnzdiag = cholmod_dump_sparse (C, "copy_sym_to_unsym", Common)) ;
    PRINT1 (("mode %d nnzdiag %d\n", mode, nnzdiag)) ;
    ASSERT (IMPLIES (mode < 0, nnzdiag == 0)) ;
    return (C) ;
}


/* ========================================================================== */
/* === cholmod_copy ========================================================= */
/* ========================================================================== */

cholmod_sparse *cholmod_copy
(
    /* inputs, not modified on output */
    cholmod_sparse *A,
    int symmetric,	/* requested symmetry of C (<0: lo, 0: both, >0 up) */
    int mode,		/* >0: numerical, 0: pattern, <0: pattern (no diag) */

    cholmod_common *Common
)
{
    cholmod_sparse *C ;
    int nrow, ncol, up, lo, values, diag, asym ;

    DEBUG (int orig) ;
    DEBUG (int nnzdiag) ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    nrow = A->nrow ;
    ncol = A->ncol ;
    if ((symmetric || A->stype) && nrow != ncol)
    {
	/* inputs invalid */
	cholmod_error (CHOLMOD_INVALID, "cholmod_copy: matrix invalid", Common);
	return (NULL) ;
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    cholmod_allocate_work (0, MAX (nrow,ncol), 0, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	/* out of memory */
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    DEBUG (orig = Common->malloc_count) ;

    values = (mode > 0) && (A->xtype != CHOLMOD_PATTERN) ;
    diag = (mode >= 0) ;
    asym = SIGN (A->stype) ;
    symmetric = SIGN (symmetric) ;
    up = (asym > 0) ;
    lo = (asym < 0) ;

    /* ---------------------------------------------------------------------- */
    /* copy the matrix */
    /* ---------------------------------------------------------------------- */

    if (asym == symmetric)
    {

	/* ------------------------------------------------------------------ */
	/* symmetry of A and C are the same */
	/* ------------------------------------------------------------------ */

	/* copy A into C, keeping the same symmetry.  If A is symmetric
	 * entries in the ignored part of A are not copied into C */
	C = cholmod_band (A, -nrow, ncol, mode, Common) ;

    }
    else if (!asym)
    {

	/* ------------------------------------------------------------------ */
	/* convert unsymmetric matrix A into a symmetric matrix C */
	/* ------------------------------------------------------------------ */

	if (symmetric > 0)
	{
	    /* C = triu (A) */
	    C = cholmod_band (A, 0, ncol, mode, Common) ;
	}
	else
	{
	    /* C = tril (A) */
	    C = cholmod_band (A, -nrow, 0, mode, Common) ;
	}
	if (Common->status < CHOLMOD_OK)
	{
	    /* out of memory */
	    ASSERT (Common->malloc_count == orig) ;
	    return (NULL) ;
	}
	C->stype = symmetric ;

    }
    else if (asym == -symmetric)
    {

	/* ------------------------------------------------------------------ */
	/* transpose a symmetric matrix */
	/* ------------------------------------------------------------------ */

	/* converting upper to lower or lower to upper */
	/* workspace: Iwork (nrow) */
	C = cholmod_transpose (A, values, NULL, NULL, 0, Common) ;
	if (!diag)
	{
	    /* remove diagonal, if requested */
	    cholmod_band_inplace (C, -nrow, ncol, -1, Common) ;
	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* create an unsymmetric copy of a symmetric matrix */
	/* ------------------------------------------------------------------ */

	C = cholmod_copy_sym_to_unsym (A, 0, mode, Common) ;
    }

    /* ---------------------------------------------------------------------- */
    /* return if error */
    /* ---------------------------------------------------------------------- */

    if (Common->status < CHOLMOD_OK)
    {
	/* out of memory */
	ASSERT (Common->malloc_count == orig) ;
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* return the result */
    /* ---------------------------------------------------------------------- */

    ASSERT (Common->malloc_count == orig + ((values) ? 4 : 3)) ;
    DEBUG (nnzdiag = cholmod_dump_sparse (C, "copy", Common)) ;
    PRINT1 (("mode %d nnzdiag %d\n", mode, nnzdiag)) ;
    ASSERT (IMPLIES (mode < 0, nnzdiag == 0)) ;
    return (C) ;
}
