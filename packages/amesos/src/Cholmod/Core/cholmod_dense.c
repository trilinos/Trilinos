/* ========================================================================== */
/* === Core/cholmod_dense =================================================== */
/* ========================================================================== */

/*
 * CHOLMOD/Core version 0.1. May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* Core utility routines for the cholmod_dense object:
 *
 * The solve routines and some of the MatrixOps and Modify routines use dense
 * matrices as inputs.  These are held in column-major order.  With a leading
 * dimension of d, the entry in row i and column j is held in x [i+j*d].
 *
 * Primary routines:
 * -----------------
 * cholmod_allocate_dense	allocate a dense matrix
 * cholmod_free_dense		free a dense matrix
 *
 * Secondary routines:
 * -------------------
 * cholmod_sparse_to_dense	create a dense matrix copy of a sparse matrix
 * cholmod_dense_to_sparse	create a sparse matrix copy of a dense matrix
 * cholmod_copy_dense		create a copy of a dense matrix
 * cholmod_copy_dense2		copy a dense matrix (pre-allocated)
 */

#include "cholmod_core.h"
#include "cholmod_internal.h"

/* ========================================================================== */
/* === cholmod_allocate_dense =============================================== */
/* ========================================================================== */

/* Allocate a dense matrix with leading dimension d.  The space is optionally
 * initialized, depending on the values parameter:
 *
 * CHOLMOD_NONE: do not initialize
 * CHOLMOD_ZERO: initialize to zero
 * CHOLMOD_ONE:  initialize to one
 * CHOLMOD_EYE:  initialize to the identity matrix
 */

cholmod_dense *cholmod_allocate_dense
(
    size_t nrow,
    size_t ncol,
    size_t d,
    int values,
    cholmod_common *Common
)
{
    cholmod_dense *X ;
    double *Xx ;
    size_t nzmax ;
    int i, n ;
    DEBUG (int orig) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    if (d < nrow)
    {
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_allocate_dense: leading dimension invalid", Common) ;
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* allocate header */
    /* ---------------------------------------------------------------------- */

    Common->status = CHOLMOD_OK ;
    DEBUG (orig = Common->malloc_count) ;

    X = cholmod_malloc (1, sizeof (cholmod_dense), Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	ASSERT (Common->malloc_count == orig) ;
	return (NULL) ;	    /* out of memory */
    }

    nzmax = MAX (1, d*ncol) ;
    PRINT1 (("cholmod_allocate_dense %d-by-%d nzmax %d values %d\n",
		nrow, ncol, nzmax, values)) ;

    X->nrow = nrow ;
    X->ncol = ncol ;
    X->nzmax = nzmax ;
    X->xtype = Common->xtype ;
    X->dtype = Common->dtype ;
    X->x = NULL ;
    X->z = NULL ;
    X->d = d ;
    X->x = cholmod_malloc (nzmax, sizeof (double), Common) ;

    if (Common->status < CHOLMOD_OK)
    {
	cholmod_free_dense (&X, Common) ;
	ASSERT (Common->malloc_count == orig) ;
	return (NULL) ;	    /* out of memory */
    }

    Xx = X->x ;

    switch (values)
    {

	case CHOLMOD_ZERO:
	    for (i = 0 ; i < ((int) nzmax) ; i++)
	    {
		Xx [i] = 0 ;
	    }
	    break ;

	case CHOLMOD_ONE:
	    for (i = 0 ; i < ((int) nzmax) ; i++)
	    {
		Xx [i] = 1 ;
	    }
	    break ;

	case CHOLMOD_EYE:
	    for (i = 0 ; i < ((int) nzmax) ; i++)
	    {
		Xx [i] = 0 ;
	    }
	    n = MIN (nrow, ncol) ;
	    for (i = 0 ; i < n ; i++)
	    {
		Xx [i + i*d] = 1 ;
	    }
	    break ;

    }

    ASSERT (Common->malloc_count == orig + 2) ;
    return (X) ;
}


/* ========================================================================== */
/* === cholmod_free_dense =================================================== */
/* ========================================================================== */

/* free a dense matrix
 *
 * workspace: none
 */

int cholmod_free_dense
(
    cholmod_dense **XHandle,
    cholmod_common *Common
)
{
    cholmod_dense *X ;

    RETURN_IF_NULL_COMMON (FALSE) ;

    if (XHandle == NULL)
    {
	/* nothing to do */
	return (TRUE) ;
    }
    X = *XHandle ;
    if (X == NULL)
    {
	/* nothing to do */
	return (TRUE) ;
    }
    X->x = cholmod_free (X->x, X->nzmax, sizeof (double), Common) ;
    *XHandle = cholmod_free ((*XHandle), 1, sizeof (cholmod_dense), Common) ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_sparse_to_dense ============================================== */
/* ========================================================================== */

/* Convert a sparse matrix to a dense matrix.  If the values parameter is
 * FALSE, then 1's are placed where entries of A are, regardless of their
 * numerical value (which may in fact be equal to zero).
 */

cholmod_dense *cholmod_sparse_to_dense
(
    /* inputs, not modified on output */
    cholmod_sparse *A,
    int values,
    cholmod_common *Common
)
{
    double aij ;
    double *Ax, *Xx ;
    int *Ap, *Ai, *Anz ;
    cholmod_dense *X ;
    int i, j, p, pend, nrow, ncol, packed ;

    DEBUG (int orig) ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    nrow = A->nrow ;
    ncol = A->ncol ;
    if (A->stype && nrow != ncol)
    {
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_sparse_to_dense: matrix invalid", Common) ;
	return (NULL) ;
    }
    Common->status = CHOLMOD_OK ;
    DEBUG (orig = Common->malloc_count) ;
    ASSERT (cholmod_dump_sparse (A, "A", Common) >= 0) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    nrow = A->nrow ;
    ncol = A->ncol ;
    packed = A->packed ;
    Ap = A->p ;
    Ai = A->i ;
    Ax = A->x ;
    Anz = A->nz ;
    values = values && (A->xtype != CHOLMOD_PATTERN) ;

    /* ---------------------------------------------------------------------- */
    /* allocate result */
    /* ---------------------------------------------------------------------- */

    X = cholmod_allocate_dense (nrow, ncol, nrow, CHOLMOD_ZERO, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory */
    }
    Xx = X->x ;

    /* ---------------------------------------------------------------------- */
    /* copy into dense matrix */
    /* ---------------------------------------------------------------------- */

    if (A->stype < 0)
    {

	/* A is symmetric with lower stored, but both parts of X are present */
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		i = Ai [p] ;
		if (i >= j)
		{
		    aij = (values) ? (Ax [p]) : 1 ;
		    Xx [i + j*nrow] = aij ;
		    Xx [j + i*nrow] = aij ;
		}
	    }
	}

    }
    else if (A->stype > 0)
    {

	/* A is symmetric with upper stored, but both parts of X are present */
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		i = Ai [p] ;
		if (i <= j)
		{
		    aij = (values) ? (Ax [p]) : 1 ;
		    Xx [i + j*nrow] = aij ;
		    Xx [j + i*nrow] = aij ;
		}
	    }
	}

    }
    else
    {

	/* both parts of A and X are present */
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		Xx [Ai [p] + j*nrow] = (values) ? (Ax [p]) : 1 ;
	    }
	}
    }

    ASSERT (Common->malloc_count == orig + 2) ;
    return (X) ;
}


/* ========================================================================== */
/* === cholmod_dense_to_sparse ============================================== */
/* ========================================================================== */

/* Convert a dense matrix to a sparse matrix, similar to the MATLAB statements:
 *
 * C = sparse (X)			symmetric = 0, values = TRUE
 * C = sparse (tril (X))		symmetric > 0, values = TRUE
 * C = sparse (triu (X))		symmetric < 0, values = TRUE
 * C = spones (sparse (X))		symmetric = 0, values = FALSE
 * C = spones (sparse (tril (X)))	symmetric > 0, values = FALSE
 * C = spones (sparse (triu (X)))	symmetric < 0, values = FALSE
 */

cholmod_sparse *cholmod_dense_to_sparse
(
    /* inputs, not modified on output */
    cholmod_dense *X,
    int symmetric,
    int values,
    cholmod_common *Common
)
{
    double xij ;
    double *Xx, *Cx ;
    int *Ci, *Cp ;
    cholmod_sparse *C ;
    int i, j, p, d, nrow, ncol, nz ;

    DEBUG (int orig) ;
    DEBUG (cholmod_dump_dense (X, "X", Common)) ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (X, NULL) ;
    nrow = X->nrow ;
    ncol = X->ncol ;
    d = X->d ;
    Xx = X->x ;
    RETURN_IF_NULL (Xx, NULL) ;
    symmetric = SIGN (symmetric) ;
    if (X->d < X->nrow)
    {
	cholmod_error (CHOLMOD_INVALID,
	    "cholmod_dense_to_sparse: matrix invalid", Common) ;
	return (NULL) ;
    }
    if (symmetric && nrow != ncol)
    {
	cholmod_error (CHOLMOD_INVALID,
	    "cholmod_dense_to_sparse: matrix unsymmetric", Common) ;
	return (NULL) ;
    }
    Common->status = CHOLMOD_OK ;
    DEBUG (orig = Common->malloc_count) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    /* count the number of nonzeros in the result */
    /* ---------------------------------------------------------------------- */

    nz = 0 ;
    if (symmetric > 0)
    {
	/* look at just triu(X) */
	for (j = 0 ; j < ncol ; j++)
	{
	    for (i = 0 ; i <= j ; i++)
	    {
		if (Xx [i+j*d] != 0)
		{
		    nz++ ;
		}
	    }
	}
    }
    else if (symmetric < 0)
    {
	/* look at just tril(X) */
	for (j = 0 ; j < ncol ; j++)
	{
	    for (i = j ; i < nrow ; i++)
	    {
		if (Xx [i+j*d] != 0)
		{
		    nz++ ;
		}
	    }
	}
    }
    else
    {
	/* look at both parts of X */
	for (j = 0 ; j < ncol ; j++)
	{
	    for (i = 0 ; i < nrow ; i++)
	    {
		if (Xx [i+j*d] != 0)
		{
		    nz++ ;
		}
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* allocate the result C */
    /* ---------------------------------------------------------------------- */

    C = cholmod_allocate_sparse (nrow, ncol, nz, TRUE, TRUE, symmetric,
	    values, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	ASSERT (Common->malloc_count == orig) ;
	return (NULL) ;	    /* out of memory */
    }
    Cp = C->p ;
    Ci = C->i ;
    Cx = C->x ;

    /* ---------------------------------------------------------------------- */
    /* copy the dense matrix X into the sparse matrix C */
    /* ---------------------------------------------------------------------- */

    p = 0 ;
    if (symmetric > 0)
    {
	/* look at just triu(X) */
	for (j = 0 ; j < ncol ; j++)
	{
	    Cp [j] = p ;
	    for (i = 0 ; i <= j ; i++)
	    {
		xij = Xx [i+j*d] ;
		if (xij != 0)
		{
		    Ci [p] = i ;
		    if (values)
		    {
			Cx [p] = xij ;
		    }
		    p++ ;
		}
	    }
	}
    }
    else if (symmetric < 0)
    {
	/* look at just tril(X) */
	for (j = 0 ; j < ncol ; j++)
	{
	    Cp [j] = p ;
	    for (i = j ; i < nrow ; i++)
	    {
		xij = Xx [i+j*d] ;
		if (xij != 0)
		{
		    Ci [p] = i ;
		    if (values)
		    {
			Cx [p] = xij ;
		    }
		    p++ ;
		}
	    }
	}
    }
    else
    {
	/* look at all of X */
	for (j = 0 ; j < ncol ; j++)
	{
	    Cp [j] = p ;
	    for (i = 0 ; i < nrow ; i++)
	    {
		xij = Xx [i+j*d] ;
		if (xij != 0)
		{
		    Ci [p] = i ;
		    if (values)
		    {
			Cx [p] = xij ;
		    }
		    p++ ;
		}
	    }
	}
    }
    ASSERT (p == nz) ;
    Cp [ncol] = nz ;

    /* ---------------------------------------------------------------------- */
    /* return result */
    /* ---------------------------------------------------------------------- */

    ASSERT (cholmod_dump_sparse (C, "C", Common) >= 0) ;
    ASSERT (Common->malloc_count == orig + ((values) ? 4 : 3)) ;
    return (C) ;
}


/* ========================================================================== */
/* === cholmod_copy_dense2 ================================================== */
/* ========================================================================== */

/* Y = X, where X and Y are both already allocated.  The leading dimensions of
 * X and Y may differ, but both must be >= the # of rows in X and Y.
 * Entries in rows nrow to d-1 are not copied from X, since the space might not
 * be initialized.  Y->nzmax is unchanged.  X->nzmax is typically
 * (X->d)*(X->ncol), but a user might modify that condition outside of any
 * CHOLMOD routine.
 */

int cholmod_copy_dense2
(
    cholmod_dense *Y,
    cholmod_dense *X,
    cholmod_common *Common
)
{
    double *Xx, *Yx ;
    int i, j, nrow, ncol, dy, dx ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (X, FALSE) ;
    RETURN_IF_NULL (Y, FALSE) ;
    if (X->nrow != Y->nrow || X->ncol != Y->ncol)
    {
	cholmod_error (CHOLMOD_INVALID,
	    "cholmod_copy_dense2: X and Y must have same dimensions", Common) ;
	return (FALSE) ;
    }
    if (X->d < X->nrow || Y->d < Y->nrow || X->x == NULL || Y->x == NULL
	    || (X->d * X->ncol) > X->nzmax || (Y->d * Y->ncol) > Y->nzmax)
    {
	cholmod_error (CHOLMOD_INVALID,
	    "cholmod_copy_dense2: X and/or Y invalid", Common) ;
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* copy */
    /* ---------------------------------------------------------------------- */

    Xx = X->x ;
    Yx = Y->x ;
    dx = X->d ;
    dy = Y->d ;
    nrow = X->nrow ;
    ncol = X->ncol ;

    for (j = 0 ; j < ncol ; j++)
    {
	for (i = 0 ; i < nrow ; i++)
	{
	    Yx [i+j*dy] = Xx [i+j*dx] ;
	}
    }
    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_copy_dense =================================================== */
/* ========================================================================== */

/* Y = X, allocating the result Y.  The leading dimension may be changed.
 * Entries in rows nrow to d-1 are not copied from X, since the space might not
 * be initialized.  Y->nzmax becomes d*ncol, which may differ from X->nzmax even
 * if d is unchanged.  X->nzmax is typically (X->d)*(X->ncol), but a user might
 * modify that condition outside of any CHOLMOD routine.
 */

cholmod_dense *cholmod_copy_dense
(
    /* inputs, not modified on output */
    cholmod_dense *X,
    size_t d,
    cholmod_common *Common
)
{
    cholmod_dense *Y ;
    DEBUG (int orig) ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (X, NULL) ;
    if (d < X->nrow)
    {
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_copy_dense: leading dimension invalid", Common) ;
	return (NULL) ;
    }
    Common->status = CHOLMOD_OK ;
    DEBUG (orig = Common->malloc_count) ;

    /* ---------------------------------------------------------------------- */
    /* allocate result */
    /* ---------------------------------------------------------------------- */

    Y = cholmod_allocate_dense (X->nrow, X->ncol, d, CHOLMOD_NONE, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	ASSERT (Common->malloc_count == orig) ;
	return (NULL) ;	    /* out of memory or d invalid */
    }

    /* ---------------------------------------------------------------------- */
    /* copy */
    /* ---------------------------------------------------------------------- */

    if (!cholmod_copy_dense2 (Y, X, Common))
    {
	/* X and/or Y are invalid */
	cholmod_free_dense (&Y, Common) ;
	ASSERT (Common->malloc_count == orig) ;
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* return result */
    /* ---------------------------------------------------------------------- */

    ASSERT (Common->malloc_count == orig + 2) ;
    return (Y) ;
}
