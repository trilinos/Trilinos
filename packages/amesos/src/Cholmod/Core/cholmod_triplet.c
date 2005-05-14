/* ========================================================================== */
/* === Core/cholmod_triplet ================================================= */
/* ========================================================================== */

/*
 * CHOLMOD/Core version 0.1. May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* Core utility routines for the cholmod_triplet object:
 *
 * A sparse matrix held in triplet form is the simplest one for a user to
 * create.  It consists of a list of nz entries in arbitrary order, held in
 * three arrays: i, j, and x, each of length nk.  The kth entry is in row i[k],
 * column j[k], with value x[k].  There may be duplicate values; if A(i,j)
 * appears more than once, its value is the sum of the entries with those row
 * and column indices.
 *
 * Primary routines:
 * -----------------
 * cholmod_allocate_triplet	allocate a triplet matrix
 * cholmod_free_triplet		free a triplet matrix
 *
 * Secondary routines:
 * -------------------
 * cholmod_reallocate_triplet	reallocate a triplet matrix
 * cholmod_sparse_to_triplet	create a triplet matrix copy of a sparse matrix
 * cholmod_triplet_to_sparse	create a sparse matrix copy of a triplet matrix
 * cholmod_copy_triplet		create a copy of a triplet matrix
 *
 * The relationship between an m-by-n cholmod_sparse matrix A and a
 * cholmod_triplet matrix (i, j, and x) is identical to how they are used in
 * the MATLAB "sparse" and "find" functions:
 *
 *	[i j x] = find (A)
 *	[m n] = size (A)
 *	A = sparse (i,j,x,m,n)
 *
 * with the exception that the cholmod_sparse matrix may be "unpacked", may
 * have either sorted or unsorted columns (depending on the option selected),
 * and may be symmetric with just the upper or lower triangular part stored.
 * Likewise, the cholmod_triplet matrix may contain just the entries in the
 * upper or lower triangular part of a symmetric matrix.
 *
 * MATLAB sparse matrices are always "packed", always have sorted columns,
 * and always store both parts of a symmetric matrix.  In some cases, MATLAB
 * behaves like CHOLMOD by ignoring entries in the upper or lower triangular
 * part of a matrix that is otherwise assumed to be symmetric (such as the
 * input to chol).  In CHOLMOD, that option is a characteristic of the object.
 * In MATLAB, that option is based on how a matrix is used as the input to
 * a function.
 *
 * The triplet matrix is provided to give the user a simple way of constructing
 * a sparse matrix.  There are very few operations supported for triplet
 * matrices.  The assumption is that they will be converted to cholmod_sparse
 * matrix form first.
 *
 * Adding two triplet matrices simply involves concatenating the contents of
 * the three arrays (i, j, and x).   To permute a triplet matrix, just replace
 * the row and column indices with their permuted values.  For example, if
 * P is a permutation vector, then P [k] = j means row/column j is the kth
 * row/column in C=P*A*P'.  In MATLAB notation, C=A(p,p).  If Pinv is an array
 * of size n and T is the triplet form of A, then:
 *
 *	Ti = T->i ;
 *	Tj = T->j ;
 *	for (k = 0 ; k < n  ; k++) Pinv [P [k]] = k ;
 *	for (k = 0 ; k < nz ; k++) Ti [k] = Pinv [Ti [k]] ;
 *	for (k = 0 ; k < nz ; k++) Tj [k] = Pinv [Tj [k]] ;
 *
 * overwrites T with the triplet form of C=P*A*P'.  The conversion
 *
 *	C = cholmod_triplet_to_sparse (T, &Common) ;
 *
 * will then return the matrix C = P*A*P'.
 *
 * Note that T->stype > 0 means that entries in the lower triangular part of
 * T are transposed into the upper triangular part when T is converted to
 * sparse matrix (cholmod_sparse) form with cholmod_triplet_to_sparse.  The
 * opposite is true for T->stype < 0.
 *
 * Since the triplet matrix T is so simple to generate, it's quite easy
 * to remove entries that you do not want, prior to converting T to the
 * cholmod_sparse form.  So if you include these entries in T, CHOLMOD
 * assumes that there must be a reason (such as the one above).  Thus,
 * no entry in a triplet matrix is ever ignored.
 *
 * Other operations, such as extacting a submatrix, horizontal and vertical
 * concatenation, multiply a triplet matrix times a dense matrix, are also
 * simple.  Multiplying two triplet matrices is not trivial; the simplest
 * method is to convert them to cholmod_sparse matrices first.
 */
 
#include "cholmod_core.h"
#include "cholmod_internal.h"

/* ========================================================================== */
/* === cholmod_allocate_triplet ============================================= */
/* ========================================================================== */

/* allocate space for a triplet matrix
 *
 * workspace: none
 */

cholmod_triplet *cholmod_allocate_triplet
(
    size_t nrow,    /* # of rows of T */
    size_t ncol,    /* # of columns of T */
    size_t nzmax,   /* max # of nonzeros of T */
    int symmetry,   /* symmetry of T */
    int values,
    cholmod_common *Common
)
{
    cholmod_triplet *T ;
    DEBUG (int orig) ;

    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = CHOLMOD_OK ;
    DEBUG (orig = Common->malloc_count) ;

    T = cholmod_malloc (1, sizeof (cholmod_triplet), Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	ASSERT (Common->malloc_count == orig) ;
	return (NULL) ;	    /* out of memory */
    }

    PRINT1 (("cholmod_allocate_triplet %d-by-%d nzmax %d values %d\n",
		nrow, ncol, nzmax, values)) ;

    nzmax = MAX (1, nzmax) ;

    T->nrow = nrow ;
    T->ncol = ncol ;
    T->nzmax = nzmax ;
    T->nnz = 0 ;
    T->stype = symmetry ;
    T->itype = Common->itype ;
    T->xtype = values ? (Common->xtype) : CHOLMOD_PATTERN ;
    T->dtype = Common->dtype ;

    T->j = NULL ;
    T->i = NULL ;
    T->x = NULL ;
    T->z = NULL ;

    /* allocate T->i, T->j, and T->x, of size nzmax */
    T->i = cholmod_malloc (nzmax, sizeof (int), Common) ;
    T->j = cholmod_malloc (nzmax, sizeof (int), Common) ;
    if (values)
    {
	T->x = cholmod_malloc (nzmax, sizeof (double), Common) ;
    }

    if (Common->status < CHOLMOD_OK)
    {
	cholmod_free_triplet (&T, Common) ;
	ASSERT (Common->malloc_count == orig) ;
	return (NULL) ;	    /* out of memory */
    }

    ASSERT (Common->malloc_count == orig + 3 + (values ? 1 : 0)) ;
    return (T) ;
}


/* ========================================================================== */
/* === cholmod_free_triplet ================================================= */
/* ========================================================================== */

/* free a triplet matrix
 *
 * workspace: none
 */

int cholmod_free_triplet
(
    cholmod_triplet **THandle,
    cholmod_common *Common
)
{
    int nz ;
    cholmod_triplet *T ;

    RETURN_IF_NULL_COMMON (FALSE) ;

    if (THandle == NULL)
    {
	/* nothing to do */
	return (TRUE) ;
    }
    T = *THandle ;
    if (T == NULL)
    {
	/* nothing to do */
	return (TRUE) ;
    }
    nz = T->nzmax ;
    T->j = cholmod_free (T->j, nz, sizeof (int), Common) ;
    T->i = cholmod_free (T->i, nz, sizeof (int), Common) ;
    T->x = cholmod_free (T->x, nz, sizeof (double), Common) ;
    *THandle = cholmod_free ((*THandle), 1, sizeof (cholmod_triplet), Common) ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_reallocate_triplet =========================================== */
/* ========================================================================== */

/* Change the size of T->i, T->j, and T->x, or allocate them if their current
 * size is zero.  T->x is not modified if T->xtype is CHOLMOD_PATTERN.
 *
 * FUTURE WORK: complex case is not yet supported.
 *
 * workspace: none
 */

int cholmod_reallocate_triplet
(
    cholmod_triplet *T,
    size_t nznew,
    cholmod_common *Common
)
{

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (T, FALSE) ;
    Common->status = CHOLMOD_OK ;
    PRINT1 (("realloc triplet %d to %d, values: %d\n",
		T->nzmax, nznew, T->xtype)) ;

    /* ---------------------------------------------------------------------- */
    /* resize the matrix */
    /* ---------------------------------------------------------------------- */

    if (T->xtype == CHOLMOD_PATTERN)
    {
	cholmod_realloc_multiple (2, 0, &(T->i), &(T->j), NULL, NULL,
		&(T->nzmax), MAX (1,nznew), Common) ;
    }
    else
    {
	ASSERT (T->xtype == CHOLMOD_REAL) ;
	cholmod_realloc_multiple (2, 1, &(T->i), &(T->j), &(T->x), NULL,
		&(T->nzmax), MAX (1,nznew), Common) ;
    }

    return (Common->status == CHOLMOD_OK) ;
}


/* ========================================================================== */
/* === cholmod_triplet_to_sparse ============================================ */
/* ========================================================================== */

/* Convert a set of triplets into a cholmod_sparse matrix.  In MATLAB notation,
 * for unsymmetric matrices:
 *
 *	A = sparse (Ti, Tj, Tx, nrow, ncol) ;
 *
 * For the symmetric upper case:
 *
 *	A = sparse (min(Ti,Tj), max(Ti,Tj), Tx, nrow, ncol) ;
 *
 * For the symmetric lower case:
 *
 *	A = sparse (max(Ti,Tj), min(Ti,Tj), Tx, nrow, ncol) ;
 *
 * If Tx is NULL, then A->x is not allocated, and only the pattern
 * of A is computed.  A is returned in packed form, and can be of any symmetry
 * (upper/lower/unsymmetric).
 *
 * workspace: Iwork (max (nrow,ncol))
 *	allocates a temporary copy of its output matrix.
 */

#include "cholmod_internal.h"

cholmod_sparse *cholmod_triplet_to_sparse
(
    /* inputs, not modified on output */
    cholmod_triplet *T,
    cholmod_common *Common
)
{
    double *Rx, *Tx ;
    int *Wj, *Rp, *Ri, *Rnz, *Ti, *Tj, *Iwork ;
    cholmod_sparse *R, *A ;
    int values, i, j, p, p1, p2, pdest, pj, k, symmetric, nrow, ncol, nz ;

    DEBUG (int orig) ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (T, NULL) ;

    symmetric = SIGN (T->stype) ;
    nrow = T->nrow ;
    ncol = T->ncol ;
    nz = T->nnz ;
    Ti = T->i ;
    Tj = T->j ;
    Tx = T->x ;

    RETURN_IF_NULL (Ti, NULL) ;
    RETURN_IF_NULL (Tj, NULL) ;
    if (symmetric && nrow != ncol)
    {
	/* inputs invalid */
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_triplet_to_sparse: matrix invalid", Common) ;
	return (NULL) ;
    }
    Common->status = CHOLMOD_OK ;
    DEBUG (cholmod_dump_triplet (T, "T", Common)) ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    cholmod_allocate_work (0, MAX (nrow, ncol), 0, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory */
    }

    DEBUG (orig = Common->malloc_count) ;

    /* ---------------------------------------------------------------------- */
    /* allocate temporary matrix R */
    /* ---------------------------------------------------------------------- */

    values = (T->xtype != CHOLMOD_PATTERN && Tx != NULL) ;

    R = cholmod_allocate_sparse (ncol, nrow, nz, FALSE, FALSE, -symmetric,
	    values, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	ASSERT (Common->malloc_count == orig) ;
	return (NULL) ;	    /* out of memory */
    }

    Rp = R->p ;
    Ri = R->i ;
    Rx = R->x ;
    Rnz = R->nz ;

    /* ---------------------------------------------------------------------- */
    /* count the entries in each row of A (also counting duplicates) */
    /* ---------------------------------------------------------------------- */

    for (i = 0 ; i < nrow ; i++)
    {
	Rnz [i] = 0 ;	
    }

    if (symmetric > 0)
    {
	for (k = 0 ; k < nz ; k++)
	{
	    i = Ti [k] ;
	    j = Tj [k] ;
	    if (i < 0 || i > nrow || j < 0 || j > ncol)
	    {
		cholmod_error (CHOLMOD_INVALID,
		    "cholmod_triplet_to_sparse: index out of range", Common) ;
		break ;
	    }
	    /* A will be symmetric with just the upper triangular part stored.
	     * Create a matrix R that is lower triangular.  Entries in the
	     * upper part of R are transposed to the lower part. */
	    Rnz [MIN (i,j)]++ ;
	}
    }
    else if (symmetric < 0)
    {
	for (k = 0 ; k < nz ; k++)
	{
	    i = Ti [k] ;
	    j = Tj [k] ;
	    if (i < 0 || i > nrow || j < 0 || j > ncol)
	    {
		cholmod_error (CHOLMOD_INVALID,
		    "cholmod_triplet_to_sparse: index out of range", Common) ;
		break ;
	    }
	    /* A will be symmetric with just the lower triangular part stored.
	     * Create a matrix R that is upper triangular.  Entries in the
	     * lower part of R are transposed to the upper part. */
	    Rnz [MAX (i,j)]++ ;
	}
    }
    else
    {
	for (k = 0 ; k < nz ; k++)
	{
	    i = Ti [k] ;
	    j = Tj [k] ;
	    if (i < 0 || i > nrow || j < 0 || j > ncol)
	    {
		cholmod_error (CHOLMOD_INVALID,
		    "cholmod_triplet_to_sparse: index out of range", Common) ;
		break ;
	    }
	    /* constructing an unsymmetric matrix */
	    Rnz [i]++ ;
	}
    }

    if (Common->status < CHOLMOD_OK)
    {
	/* triplet matrix is invalid */
	cholmod_free_sparse (&R, Common) ;
	ASSERT (Common->malloc_count == orig) ;
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* construct the row pointers */
    /* ---------------------------------------------------------------------- */

    p = 0 ;
    for (i = 0 ; i < nrow ; i++)
    {
	Rp [i] = p ;
	p += Rnz [i] ;
    }
    Rp [nrow] = p ;

    /* use Wj (i/l/l) as temporary row pointers [ */
    Iwork = Common->Iwork ;
    Wj = Iwork ;		/* size MAX (nrow,ncol) FUTURE WORK: (i/l/l) */
    for (i = 0 ; i < nrow ; i++)
    {
	Wj [i] = Rp [i] ;
    }

    /* ---------------------------------------------------------------------- */
    /* construct the row form */
    /* ---------------------------------------------------------------------- */

    if (symmetric > 0)
    {
	for (k = 0 ; k < nz ; k++)
	{
	    i = Ti [k] ;
	    j = Tj [k] ;
	    if (i < j)
	    {
		/* place triplet (j,i,x) in column i of R */
		p = Wj [i]++ ;
		Ri [p] = j ;
		if (values)
		{
		    Rx [p] = Tx [k] ;
		}
	    }
	    else
	    {
		/* place triplet (i,j,x) in column j of R */
		p = Wj [j]++ ;
		Ri [p] = i ;
		if (values)
		{
		    Rx [p] = Tx [k] ;
		}
	    }
	}
    }
    else if (symmetric < 0)
    {
	for (k = 0 ; k < nz ; k++)
	{
	    i = Ti [k] ;
	    j = Tj [k] ;
	    if (i > j)
	    {
		/* place triplet (j,i,x) in column i of R */
		p = Wj [i]++ ;
		Ri [p] = j ;
		if (values)
		{
		    Rx [p] = Tx [k] ;
		}
	    }
	    else
	    {
		/* place triplet (i,j,x) in column j of R */
		p = Wj [j]++ ;
		Ri [p] = i ;
		if (values)
		{
		    Rx [p] = Tx [k] ;
		}
	    }
	}
    }
    else
    {
	for (k = 0 ; k < nz ; k++)
	{
	    /* place triplet (i,j,x) in column i of R */
	    p = Wj [Ti [k]]++ ;
	    Ri [p] = Tj [k] ;
	    if (values)
	    {
		Rx [p] = Tx [k] ;
	    }
	}
    }

    /* done using Wj (i/l/l) as temporary row pointers ] */

    /* ---------------------------------------------------------------------- */
    /* sum up duplicates */
    /* ---------------------------------------------------------------------- */

    /* use Wj (i/l/l) of size ncol to keep track of duplicates in each row [ */
    for (j = 0 ; j < ncol ; j++)
    {
	Wj [j] = EMPTY ;
    }

    for (i = 0 ; i < nrow ; i++)
    {
	p1 = Rp [i] ;
	p2 = Rp [i+1] ;
	pdest = p1 ;
	/* at this point Wj [j] < p1 holds true for all columns j, because
	 * Ri/Rx is stored in row oriented manner */
	for (p = p1 ; p < p2 ; p++)
	{
	    j = Ri [p] ;
	    pj = Wj [j] ;
	    if (pj >= p1)
	    {
		/* this column index j is already in row i at position pj;
		 * sum up the duplicate entry */
		if (values)
		{
		    Rx [pj] += Rx [p] ;
		}
	    }
	    else
	    {
		/* keep the entry and keep track in Wj [j] for case above */
		Wj [j] = pdest ;
		if (pdest != p)
		{
		    Ri [pdest] = j ;
		    if (values)
		    {
			Rx [pdest] = Rx [p] ;
		    }
		}
		pdest++ ;
	    }
	}
	Rnz [i] = pdest - p1 ;
    }
    /* done using Wj to keep track of duplicate entries in each row ] */

    ASSERT (cholmod_dump_sparse (R, "R", Common) >= 0) ;

    /* ---------------------------------------------------------------------- */
    /* A = R' */
    /* ---------------------------------------------------------------------- */

    /* workspace: Iwork (R->nrow), which is A->ncol */
    A = cholmod_transpose (R, values, NULL, NULL, 0, Common) ;
    cholmod_free_sparse (&R, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	ASSERT (Common->malloc_count == orig) ;
	return (NULL) ;	    /* out of memory */
    }

    /* ---------------------------------------------------------------------- */
    /* return result */
    /* ---------------------------------------------------------------------- */

    ASSERT (Common->malloc_count == orig + ((values) ? 4 : 3)) ;
    return (A) ;
}


/* ========================================================================== */
/* === cholmod_sparse_to_triplet ============================================ */
/* ========================================================================== */

/* Converts a sparse column-oriented matrix to triplet form.
 *
 * workspace: none
 */

cholmod_triplet *cholmod_sparse_to_triplet
(
    /* inputs, not modified on output */
    cholmod_sparse *A,
    cholmod_common *Common
)
{
    double *Ax, *Tx ;
    int *Ap, *Ai, *Ti, *Tj, *Anz ;
    cholmod_triplet *T ;
    int i, values, p, pend, k, j, nrow, ncol, nz, symmetric, packed, up, lo,
	both ;

    DEBUG (int orig) ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    symmetric = SIGN (A->stype) ;
    nrow = A->nrow ;
    ncol = A->ncol ;
    if (symmetric && nrow != ncol)
    {
	/* inputs invalid */
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_sparse_to_triplet: matrix invalid", Common) ;
	return (NULL) ;
    }
    Common->status = CHOLMOD_OK ;

    DEBUG (orig = Common->malloc_count) ;
    ASSERT (cholmod_dump_sparse (A, "A", Common) >= 0) ;

    /* ---------------------------------------------------------------------- */
    /* allocate triplet matrix */
    /* ---------------------------------------------------------------------- */

    nz = cholmod_nnz (A, Common) ;
    values = (A->xtype != CHOLMOD_PATTERN) ;
    T = cholmod_allocate_triplet (nrow, ncol, nz, A->stype, values, Common);

    if (Common->status < CHOLMOD_OK)
    {
	ASSERT (Common->malloc_count == orig) ;
	return (NULL) ;	    /* out of memory */
    }

    /* ---------------------------------------------------------------------- */
    /* convert to a sparse matrix */
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;
    Ai = A->i ;
    Ax = A->x ;
    Anz = A->nz ;
    packed = A->packed ;

    Ti = T->i ;
    Tj = T->j ;
    Tx = T->x ;
    T->stype = A->stype ;

    both = (A->stype == 0) ;
    up = (A->stype > 0) ;
    lo = (A->stype < 0) ;

    k = 0 ;

    for (j = 0 ; j < ncol ; j++)
    {
	p = Ap [j] ;
	pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	for ( ; p < pend ; p++)
	{
	    i = Ai [p] ;
	    if (both || (up && i <= j) || (lo && i >= j))
	    {
		Ti [k] = Ai [p] ;
		Tj [k] = j ;
		if (values)
		{
		    Tx [k] = Ax [p] ;
		}
		k++ ;
		ASSERT (k <= nz) ;
	    }
	}
    }

    T->nnz = k ;

    /* ---------------------------------------------------------------------- */
    /* return result */
    /* ---------------------------------------------------------------------- */

    ASSERT (cholmod_dump_triplet (T, "T", Common)) ;
    ASSERT (Common->malloc_count == orig + ((values) ? 4 : 3)) ;
    return (T) ;
}


/* ========================================================================== */
/* === cholmod_copy_triplet ================================================= */
/* ========================================================================== */

/* Create an exact copy of a triplet matrix, with two exceptions:
 *
 * (1) Entries in unused space are not copied (they might not be initialized,
 *	and copying them would cause program checkers such as purify and
 *	valgrind to complain).
 *
 * (2) C->xtype becomes CHOLMOD_PATTERN and C->x becomes NULL if either T->x
 *	is NULL or T->xtype is CHOLMOD_PATTERN.  The two logical statements
 *	(T->x == NULL) and (T->xtype == CHOLMOD_PATTERN) are kept logically
 *	equivalent by CHOLMOD, but a user might modify one without making a
 *	change to the other, outside of any CHOLMOD routine.
 */

cholmod_triplet *cholmod_copy_triplet
(
    cholmod_triplet *T,
    cholmod_common *Common
)
{
    double *Tx, *Cx ;
    int *Ci, *Cj, *Ti, *Tj ;
    cholmod_triplet *C ;
    int values, k, nz ;

    DEBUG (int orig) ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (T, NULL) ;
    nz = T->nnz ;
    Ti = T->i ;
    Tj = T->j ;
    Tx = T->x ;
    RETURN_IF_NULL (Ti, NULL) ;
    RETURN_IF_NULL (Tj, NULL) ;
    Common->status = CHOLMOD_OK ;
    DEBUG (orig = Common->malloc_count) ;
    DEBUG (cholmod_dump_triplet (T, "T orig", Common)) ;

    /* ---------------------------------------------------------------------- */
    /* allocate copy */
    /* ---------------------------------------------------------------------- */

    values = (Tx != NULL && T->xtype != CHOLMOD_PATTERN) ;

    C = cholmod_allocate_triplet (T->nrow, T->ncol, T->nzmax, T->stype,
	    values, Common) ;

    if (Common->status < CHOLMOD_OK)
    {
	ASSERT (Common->malloc_count == orig) ;
	return (NULL) ;	    /* out of memory */
    }

    /* ---------------------------------------------------------------------- */
    /* copy the triplet matrix */
    /* ---------------------------------------------------------------------- */

    Ci = C->i ;
    Cj = C->j ;
    Cx = C->x ;
    C->nnz = nz ;

    for (k = 0 ; k < nz ; k++)
    {
	Ci [k] = Ti [k] ;
    }
    for (k = 0 ; k < nz ; k++)
    {
	Cj [k] = Tj [k] ;
    }
    if (values)
    {
	for (k = 0 ; k < nz ; k++)
	{
	    Cx [k] = Tx [k] ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* return the result */
    /* ---------------------------------------------------------------------- */

    ASSERT (cholmod_dump_triplet (C, "C triplet copy", Common)) ;
    ASSERT (Common->malloc_count == orig + ((values) ? 4 : 3)) ;
    return (C) ;
}
