/* ========================================================================== */
/* === Core/cholmod_sparse ================================================== */
/* ========================================================================== */

/*
 * CHOLMOD/Core version 0.1. May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* Core utility routines for the cholmod_sparse object:
 *
 * A sparse matrix is held in compressed column form.  In the basic type
 * ("packed", which corresponds to a MATLAB sparse matrix), an n-by-n matrix
 * with nz entries is held in three arrays: p of size n+1, i of size nz, and x
 * of size nz.  Row indices of column j are held in i [p [j] ... p [j+1]-1] and
 * in the same locations in x.  There may be no duplicate entries in a column.
 * Row indices in each column may be sorted or unsorted (CHOLMOD keeps track).
 *
 * Primary routines:
 * -----------------
 * cholmod_allocate_sparse	allocate a sparse matrix
 * cholmod_free_sparse		free a sparse matrix
 *
 * Secondary routines:
 * -------------------
 * cholmod_reallocate_sparse	change the size (# entries) of sparse matrix
 * cholmod_nnz			number of nonzeros in a sparse matrix
 * cholmod_speye		sparse identity matrix
 * cholmod_zero			sparse zero matrix
 * cholmod_transpose_unsym	transpose unsymmetric sparse matrix
 * cholmod_transpose_sym	transpose symmetric sparse matrix
 * cholmod_transpose		transpose sparse matrix
 * cholmod_sort			sort row indices in each column of sparse matrix
 * cholmod_copy_sparse		create a copy of a sparse matrix
 */

#include "cholmod_core.h"
#include "cholmod_internal.h"

/* ========================================================================== */
/* === cholmod_allocate_sparse ============================================== */
/* ========================================================================== */

/* Allocate space for a matrix.  The contents are not initialized.
 *
 * workspace: none
 */

cholmod_sparse *cholmod_allocate_sparse
(
    size_t nrow,    /* # of rows of A */
    size_t ncol,    /* # of columns of A */
    size_t nzmax,   /* max # of nonzeros of A */
    int sorted,	    /* TRUE if columns of A will be sorted, FALSE otherwise */
    int packed,	    /* TRUE if A is packed, FALSE otherwise */
    int stype,	    /* stype of A */
    int values,	    /* TRUE: allocate space for numerical values */
    cholmod_common *Common
)
{
    cholmod_sparse *A ;
    DEBUG (int orig) ;

    RETURN_IF_NULL_COMMON (FALSE) ;
    if (stype != 0 && nrow != ncol)
    {
	cholmod_error (CHOLMOD_INVALID,
		"rectangular matrix with stype != 0 invalid", Common) ;
	return (NULL) ;
    }
    Common->status = CHOLMOD_OK ;
    DEBUG (orig = Common->malloc_count) ;

    A = cholmod_malloc (1, sizeof (cholmod_sparse), Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	ASSERT (Common->malloc_count == orig) ;
	return (NULL) ;	    /* out of memory */
    }
    PRINT1 (("cholmod_allocate_sparse %d-by-%d nzmax %d sorted %d packed %d"
		" values %d\n", nrow, ncol, nzmax, sorted, packed, values)) ;

    nzmax = MAX (1, nzmax) ;

    A->nrow = nrow ;
    A->ncol = ncol ;
    A->nzmax = nzmax ;
    A->packed = packed ;    /* default is packed (A->nz not present) */
    A->stype = stype ;
    A->itype = Common->itype ;
    A->xtype = values ? (Common->xtype) : CHOLMOD_PATTERN ;
    A->dtype = Common->dtype ;

    A->nz = NULL ;
    A->p = NULL ;
    A->i = NULL ;
    A->x = NULL ;
    A->z = NULL ;

    /* A 1-by-m matrix always has sorted columns */
    A->sorted = (nrow <= 1) ? TRUE : sorted ;

    /* allocate O(ncol) space */
    A->p = cholmod_malloc (ncol+1, sizeof (int), Common) ;
    if (!packed)
    {
	A->nz = cholmod_malloc (ncol, sizeof (int), Common) ;
    }

    /* allocate A->i and A->x, of size nzmax */
    A->i = cholmod_malloc (nzmax, sizeof (int), Common) ;
    if (values)
    {
	A->x = cholmod_malloc (nzmax, sizeof (double), Common) ;
    }

    if (Common->status < CHOLMOD_OK)
    {
	cholmod_free_sparse (&A, Common) ;
	ASSERT (Common->malloc_count == orig) ;
	return (NULL) ;	    /* out of memory */
    }

    ASSERT (Common->malloc_count == orig +
	    3 + (values ? 1 : 0) + ((packed) ? 0 : 1)) ;
    return (A) ;
}


/* ========================================================================== */
/* === cholmod_free_sparse ================================================== */
/* ========================================================================== */

/* free a sparse matrix
 *
 * workspace: none
 */

int cholmod_free_sparse
(
    cholmod_sparse **AHandle,
    cholmod_common *Common
)
{
    int n, nz ;
    cholmod_sparse *A ;

    RETURN_IF_NULL_COMMON (FALSE) ;

    if (AHandle == NULL)
    {
	/* nothing to do */
	return (TRUE) ;
    }
    A = *AHandle ;
    if (A == NULL)
    {
	/* nothing to do */
	return (TRUE) ;
    }
    n = A->ncol ;
    nz = A->nzmax ;
    A->p     = cholmod_free (A->p,  n+1, sizeof (int), Common) ;
    A->i     = cholmod_free (A->i,  nz,  sizeof (int), Common) ;
    A->x     = cholmod_free (A->x,  nz,  sizeof (double), Common) ;
    A->nz    = cholmod_free (A->nz, n,   sizeof (int), Common) ;
    *AHandle = cholmod_free ((*AHandle), 1, sizeof (cholmod_sparse), Common) ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_reallocate_sparse ============================================ */
/* ========================================================================== */

/* Change the size of A->i and A->x, or allocate them if their current size
 * is zero.  A->x is not modified if A->xtype is CHOLMOD_PATTERN.
 *
 * FUTURE WORK: complex case is not yet supported.
 * 
 * workspace: none
 */

int cholmod_reallocate_sparse
(
    cholmod_sparse *A,
    size_t nznew,
    cholmod_common *Common
)
{

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    Common->status = CHOLMOD_OK ;
    PRINT1 (("realloc matrix %d to %d, values: %d\n",
		A->nzmax, nznew, A->xtype)) ;

    /* ---------------------------------------------------------------------- */
    /* resize the matrix */
    /* ---------------------------------------------------------------------- */

    if (A->xtype == CHOLMOD_PATTERN)
    {
	A->i = cholmod_realloc (A->i, &(A->nzmax), MAX (1,nznew), sizeof (int),
		Common) ;
    }
    else
    {
	ASSERT (A->xtype == CHOLMOD_REAL) ;
	cholmod_realloc_multiple (1, 1, &(A->i), NULL, &(A->x), NULL,
		&(A->nzmax), MAX (1,nznew), Common) ;
    }

    return (Common->status == CHOLMOD_OK) ;
}


/* ========================================================================== */
/* === cholmod_speye ======================================================== */
/* ========================================================================== */

/* Return a sparse identity matrix. */

cholmod_sparse *cholmod_speye
(
    size_t nrow,
    size_t ncol,
    size_t nzmax,
    int values,		    /* with A->x if values is TRUE */
    cholmod_common *Common
)
{
    double *Ax ;
    cholmod_sparse *A ;
    int *Ap, *Ai ;
    int j, n ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate the matrix */
    /* ---------------------------------------------------------------------- */

    n = MIN (nrow, ncol) ;
    nzmax = MAX ((int) nzmax, n) ;
    A = cholmod_allocate_sparse (nrow, ncol, nzmax, TRUE, TRUE, 0, values,
	    Common) ;

    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory */
    }

    /* ---------------------------------------------------------------------- */
    /* create the identity matrix */
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;
    Ai = A->i ;
    Ax = A->x ;

    for (j = 0 ; j < n ; j++)
    {
	Ap [j] = j ;
    }
    for (j = n ; j <= ((int) ncol) ; j++)
    {
	Ap [j] = n ;
    }
    for (j = 0 ; j < n ; j++)
    {
	Ai [j] = j ;
    }
    if (values)
    {
	for (j = 0 ; j < n ; j++)
	{
	    Ax [j] = 1 ;
	}
    }

    return (A) ;
}


/* ========================================================================== */
/* === cholmod_zero ========================================================= */
/* ========================================================================== */

/* Return a sparse zero matrix. */

cholmod_sparse *cholmod_zero
(
    size_t nrow,
    size_t ncol,
    size_t nzmax,
    int values,		    /* with A->x if values is TRUE */
    cholmod_common *Common
)
{
    cholmod_sparse *A ;
    int *Ap ;
    int j ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate the matrix */
    /* ---------------------------------------------------------------------- */

    A = cholmod_allocate_sparse (nrow, ncol, nzmax, TRUE, TRUE, 0, values,
	    Common) ;

    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory */
    }

    /* ---------------------------------------------------------------------- */
    /* create the zero matrix */
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;
    for (j = 0 ; j <= ((int) ncol) ; j++)
    {
	Ap [j] = 0 ;
    }
    return (A) ;
}


/* ========================================================================== */
/* === cholmod_nnz ========================================================== */
/* ========================================================================== */

/* Return the number of entries in a sparse matrix.
 *
 * workspace: none
 * int overflow cannot occur, since the matrix is already allocated.
 */

long cholmod_nnz
(
    cholmod_sparse *A,
    cholmod_common *Common
)
{
    int *Ap, *Anz ;
    size_t nz ;
    int j, ncol ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (EMPTY) ;
    RETURN_IF_NULL (A, EMPTY) ;
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* return nnz (A) */
    /* ---------------------------------------------------------------------- */

    ncol = A->ncol ;
    if (A->packed)
    {
	Ap = A->p ;
	RETURN_IF_NULL (Ap, EMPTY) ;
	nz = Ap [ncol] ;
    }
    else
    {
	Anz = A->nz ;
	RETURN_IF_NULL (Anz, EMPTY) ;
	nz = 0 ;
	for (j = 0 ; j < ncol ; j++)
	{
	    nz += MAX (0, Anz [j]) ;
	}
    }
    return (nz) ;
}


/* ========================================================================== */
/* === cholmod_transpose routines =========================================== */
/* ========================================================================== */

/* Compute the transpose or permuted transpose of a matrix.
 *
 * ---------------------------------------
 * Unsymmetric case: A->stype is zero.
 * ---------------------------------------
 *
 * Computes F = A', F = A (:,f)' or F = A (p,f)', except that the indexing by
 * f does not work the same as the MATLAB notation (see below).  A->stype
 * is zero, which denotes that both the upper and lower triangular parts of
 * A are present (and used).  A may in fact be symmetric in pattern and/or
 * value; A->stype just denotes which part of A are stored.  A may be
 * rectangular.
 *
 * p is a permutation of 0:m-1, and f is a subset of 0:n-1, where A is m-by-n.
 * There can be no duplicate entries in p, but f may include duplicates.
 *
 * The set f is held in fset and fsize.
 *	fset = NULL means ":" in MATLAB. fset is ignored.
 *	fset != NULL means f = fset [0..fset-1].
 *	fset != NULL and fsize = 0 means f is the empty set.
 *	There can be no duplicates in fset (this condition is not checked).
 *
 * Columns not in the set f are considered to be zero.  That is,
 * if A is 5-by-10 then F = A (:,[3 4])' is not 2-by-5, but 10-by-5, and rows
 * 3 and 4 of F are equal to columns 3 and 4 of A (the other rows of F are
 * zero).  More precisely, in MATLAB notation:
 *
 *	[m n] = size (A) ;
 *	F = A ;
 *	notf = ones (1,n) ;
 *	notf (f) = 0 ;
 *	F (:, find (notf)) = 0
 *	F = F'
 *
 * If you want the MATLAB equivalent F=A(p,f) operation, use cholmod_submatrix
 * instead (which does not compute the transpose).
 *
 * F->nzmax must be large enough to hold the matrix F.  It is not modified.
 * If F->nz is present then F->nz [j] = # of entries in column j of F.
 *
 * A can be sorted or unsorted, with packed or unpacked columns.
 *
 * If f is present and not sorted in ascending order, then F is unsorted
 * (that is, it may contain columns whose row indices do not appear in
 * ascending order).  Otherwise, F is sorted (the row indices in each
 * column of F appear in strictly ascending order).
 *
 * F is returned in packed or unpacked form, depending on F->packed on input.
 * If F->packed is false, then F is returned in unpacked form
 * (F->nz must be present).  Each row i of F is large enough to hold all the
 * entries in row i of A, even if f is provided.  That is, F->i and F->x
 * [F->p [i] .. F->p [i] + F->nz [i] - 1] contains all entries in A (i,f),
 * but F->p [i+1] - F->p [i] is equal to the number of nonzeros in A (i,:),
 * not just A (i,f).
 *
 * The cholmod_transpose_unsym routine is the only operation in CHOLMOD
 * that can produce an unpacked matrix.
 *
 * ---------------------------------------
 * Symmetric case: A->stype is nonzero.
 * ---------------------------------------
 *
 * Computes F = A' or F = A(p,p)', the transpose or permuted transpose, where
 * A->stype is nonzero.
 *
 * If A->stype > 0, then A is a symmetric matrix where just the upper part
 * of the matrix is stored.  Entries in the lower triangular part may be
 * present, but are ignored.  A must be square.  If F=A', then F is returned
 * sorted; otherwise F is unsorted for the F=A(p,p)' case.
 *
 * There can be no duplicate entries in p.
 * The fset and fsize parameters are not used.
 *
 * FUTURE WORK: the "values" parameter may be extended to handle both the
 * complex conjugate transpose and the array transpose cases.  Zero will still
 * mean pattern-only.
 */


/* ========================================================================== */
/* === cholmod_transpose_unsym ============================================== */
/* ========================================================================== */

/* F=A' or F=A(:,f)' or F=A(p,f)', where A is unsymmetric and F is already
 * allocated.
 *
 * workspace:
 * Iwork (MAX (nrow,ncol)) if fset is present
 * Iwork (nrow) if fset is NULL
 */

int cholmod_transpose_unsym
(
    /* inputs, not modified */
    cholmod_sparse *A,	/* nrow-by-ncol, stored in column form */
    int values,		/* TRUE if values to be transposed, else FALSE */
    void *Perm_p,	/* size nrow, if present (can be NULL) */
    void *fset_p,	/* size fsize, if present (can be NULL) */
    size_t fsize,

    /* output, allocated but contents undefined on input */
    cholmod_sparse *F,	/* A (Perm,f) stored in row form */

    cholmod_common *Common
)
{
    double *Ax, *Fx ;
    int *Ap, *Anz, *Ai, *Fp, *Fnz, *Fj, *Wi, *Perm, *fset, *Iwork ;
    int i, j, p, k, pend, nrow, ncol, Apacked, use_fset, fp, Fsorted, jlast,
	Fpacked, jj, permute, nf, ineed ;

    Perm = Perm_p ;
    fset = fset_p ;
    nf = fsize ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (F, FALSE) ;
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    use_fset = (fset != NULL) ;
    nrow = A->nrow ;
    ncol = A->ncol ;
    ineed = nrow + (use_fset) ? ncol : 0 ;
    cholmod_allocate_work (0, ineed, 0, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;	/* out of memory */
    }

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;		/* size A->ncol+1, column pointers of A */
    Ai = A->i ;		/* size nz = Ap [A->ncol], row indices of A */
    Ax = A->x ;		/* size nz, numeric values of A */
    Anz = A->nz ;
    Apacked = A->packed ;
    ASSERT (IMPLIES (!Apacked, Anz != NULL)) ;

    permute = (Perm != NULL) ;

    F->ncol = nrow ;
    F->nrow = ncol ;
    Fp = F->p ;		/* size A->nrow+1, row pointers of F */
    Fj = F->i ;		/* size nz, column indices of F */
    Fx = F->x ;		/* size nz, numeric values of F */
    Fnz = F->nz ;
    values = values && (F->xtype != CHOLMOD_PATTERN) ;
    Fpacked = F->packed ;
    ASSERT (IMPLIES (!Fpacked, Fnz != NULL)) ;

    nf = (use_fset) ? nf : ncol ;

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    Iwork = Common->Iwork ;
    Wi = Iwork ;		/* size nrow. FUTURE WORK: Wi is (i/l/l) */

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    if (permute)
    {
	for (i = 0 ; i < nrow ; i++)
	{
	    Wi [i] = 1 ;
	}
	for (k = 0 ; k < nrow ; k++)
	{
	    i = Perm [k] ;
	    if (i < 0 || i > nrow || Wi [i] == 0)
	    {
		cholmod_error (CHOLMOD_INVALID, "invalid permutation", Common) ;
		return (FALSE) ;
	    }
	    Wi [i] = 0 ;
	}
    }

    if (use_fset)
    {
	for (j = 0 ; j < ncol ; j++)
	{
	    Wi [j] = 1 ;
	}
	for (k = 0 ; k < nf ; k++)
	{
	    j = fset [k] ;
	    if (j < 0 || j > ncol || Wi [j] == 0)
	    {
		cholmod_error (CHOLMOD_INVALID, "invalid fset", Common) ;
		return (FALSE) ;
	    }
	    Wi [j] = 0 ;
	}
    }

    /* Perm and fset are now valid */
    ASSERT (cholmod_dump_perm (Perm, nrow, nrow, "Perm", Common)) ;
    ASSERT (cholmod_dump_perm (fset, nf, ncol, "fset", Common)) ;

    /* ---------------------------------------------------------------------- */
    /* count the entries in each row of A or A(:,f) */
    /* ---------------------------------------------------------------------- */

    for (i = 0 ; i < nrow ; i++)
    {
	Wi [i] = 0 ;
    }

    jlast = EMPTY ;
    Fsorted = TRUE ;

    if (use_fset)
    {
	/* count entries in each row of A(:,f) */
	for (jj = 0 ; jj < nf ; jj++)
	{
	    j = fset [jj] ;
	    if (j <= jlast)
	    {
		Fsorted = FALSE ;
	    }
	    p = Ap [j] ;
	    pend = (Apacked) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		Wi [Ai [p]]++ ;
	    }
	    jlast = j ;
	}

	/* save the nz counts if F is unpacked, and recount all of A */
	if (!Fpacked)
	{
	    if (permute)
	    {
		for (i = 0 ; i < nrow ; i++)
		{
		    Fnz [i] = Wi [Perm [i]] ;
		}
	    }
	    else
	    {
		for (i = 0 ; i < nrow ; i++)
		{
		    Fnz [i] = Wi [i] ;
		}
	    }
	    for (i = 0 ; i < nrow ; i++)
	    {
		Wi [i] = 0 ;
	    }

	    /* count entries in each row of A */
	    for (j = 0 ; j < ncol ; j++)
	    {
		p = Ap [j] ;
		pend = (Apacked) ? (Ap [j+1]) : (p + Anz [j]) ;
		for ( ; p < pend ; p++)
		{
		    Wi [Ai [p]]++ ;
		}
	    }
	}

    }
    else
    {

	/* count entries in each row of A */
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (Apacked) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		Wi [Ai [p]]++ ;
	    }
	}

	/* save the nz counts if F is unpacked */
	if (!Fpacked)
	{
	    if (permute)
	    {
		for (i = 0 ; i < nrow ; i++)
		{
		    Fnz [i] = Wi [Perm [i]] ;
		}
	    }
	    else
	    {
		for (i = 0 ; i < nrow ; i++)
		{
		    Fnz [i] = Wi [i] ;
		}
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* compute the row pointers */
    /* ---------------------------------------------------------------------- */

    p = 0 ;
    if (permute)
    {
	for (i = 0 ; i < nrow ; i++)
	{
	    Fp [i] = p ;
	    p += Wi [Perm [i]] ;
	}
	for (i = 0 ; i < nrow ; i++)
	{
	    Wi [Perm [i]] = Fp [i] ;
	}
    }
    else
    {
	for (i = 0 ; i < nrow ; i++)
	{
	    Fp [i] = p ;
	    p += Wi [i] ;
	}
	for (i = 0 ; i < nrow ; i++)
	{
	    Wi [i] = Fp [i] ;
	}
    }
    Fp [nrow] = p ;

    /* ---------------------------------------------------------------------- */
    /* construct the transpose */
    /* ---------------------------------------------------------------------- */

    if (use_fset)
    {
	if (Apacked)
	{
	    if (values)
	    {
		/* use fset, A packed, values */
		for (jj = 0 ; jj < nf ; jj++)
		{
		    j = fset [jj] ;
		    pend = Ap [j+1] ;
		    for (p = Ap [j] ; p < pend ; p++)
		    {
			fp = Wi [Ai [p]]++ ;
			Fj [fp] = j ;
			Fx [fp] = Ax [p] ;
		    }
		    jlast = j ;
		}
	    }
	    else
	    {
		/* use fset, A packed, no values */
		for (jj = 0 ; jj < nf ; jj++)
		{
		    j = fset [jj] ;
		    pend = Ap [j+1] ;
		    for (p = Ap [j] ; p < pend ; p++)
		    {
			Fj [Wi [Ai [p]]++] = j ;
		    }
		}
	    }
	}
	else
	{
	    if (values)
	    {
		/* use fset, A not packed, values */
		for (jj = 0 ; jj < nf ; jj++)
		{
		    j = fset [jj] ;
		    p = Ap [j] ;
		    pend = p + Anz [j] ;
		    for ( ; p < pend ; p++)
		    {
			fp = Wi [Ai [p]]++ ;
			Fj [fp] = j ;
			Fx [fp] = Ax [p] ;
		    }
		    jlast = j ;
		}
	    }
	    else
	    {
		/* use fset, A not packed, no values */
		for (jj = 0 ; jj < nf ; jj++)
		{
		    j = fset [jj] ;
		    p = Ap [j] ;
		    pend = p + Anz [j] ;
		    for ( ; p < pend ; p++)
		    {
			Fj [Wi [Ai [p]]++] = j ;
		    }
		}
	    }
	}
    }
    else
    {
	if (Apacked)
	{
	    if (values)
	    {
		/* A packed, values */
		for (j = 0 ; j < ncol ; j++)
		{
		    pend = Ap [j+1] ;
		    for (p = Ap [j] ; p < pend ; p++)
		    {
			fp = Wi [Ai [p]]++ ;
			Fj [fp] = j ;
			Fx [fp] = Ax [p] ;
		    }
		}
	    }
	    else
	    {
		/* A packed, no values */
		for (j = 0 ; j < ncol ; j++)
		{
		    pend = Ap [j+1] ;
		    for (p = Ap [j] ; p < pend ; p++)
		    {
			Fj [Wi [Ai [p]]++] = j ;
		    }
		}
	    }
	}
	else
	{
	    if (values)
	    {
		/* A unpacked, values */
		for (j = 0 ; j < ncol ; j++)
		{
		    p = Ap [j] ;
		    pend = p + Anz [j] ;
		    for ( ; p < pend ; p++)
		    {
			fp = Wi [Ai [p]]++ ;
			Fj [fp] = j ;
			Fx [fp] = Ax [p] ;
		    }
		}
	    }
	    else
	    {
		/* A unpacked, no values */
		for (j = 0 ; j < ncol ; j++)
		{
		    p = Ap [j] ;
		    pend = p + Anz [j] ;
		    for ( ; p < pend ; p++)
		    {
			Fj [Wi [Ai [p]]++] = j ;
		    }
		}
	    }
	}
    }

    F->sorted = Fsorted ;

    ASSERT (cholmod_dump_sparse (F, "output F unsym", Common) >= 0) ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_transpose_sym ================================================ */
/* ========================================================================== */

/* F = A' or F = A(p,p)' where A is symmetric and F is already allocated.
 *
 * workspace:  Iwork (nrow) if Perm NULL, Iwork (2*nrow) if Perm non-NULL.
 */

int cholmod_transpose_sym
(
    /* inputs, not modified */
    cholmod_sparse *A,
    int values,		/* TRUE if values to be transposed, else FALSE */
    void *Perm_p,	/* size nrow, if present (can be NULL) */

    /* output, allocated but contents undefined on input */
    cholmod_sparse *F,	/* A (Perm,Perm)' */

    cholmod_common *Common  /* Iwork (nrow, or 2*nrow if Perm not NULL) */
)
{
    double *Ax, *Fx ;
    int *Ap, *Anz, *Ai, *Fp, *Fj, *Wi, *Pinv, *Perm, *Iwork ;
    int p, pend, nrow, ncol, packed, fp, upper, permute, jold, n, i, j, k ;

    Perm = Perm_p ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (F, FALSE) ;
    nrow = A->nrow ;
    ncol = A->ncol ;
    if (nrow != ncol || A->stype == 0)
    {
	/* this routine handles square symmetric matrices only */
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_transpose_sym: matrix must be square", Common) ;
	return (FALSE) ;
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    permute = (Perm != NULL) ;
    cholmod_allocate_work (0, (permute) ? 2*nrow : nrow, 0, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;	/* out of memory */
    }

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    n = A->nrow ;

    Ap = A->p ;		/* size A->ncol+1, column pointers of A */
    Ai = A->i ;		/* size nz = Ap [A->ncol], row indices of A */
    Ax = A->x ;		/* size nz, numeric values of A */
    Anz = A->nz ;
    packed = A->packed ;
    ASSERT (IMPLIES (!packed, Anz != NULL)) ;
    upper = (A->stype > 0) ;

    F->ncol = n ;
    F->nrow = n ;
    Fp = F->p ;		/* size A->nrow+1, row pointers of F */
    Fj = F->i ;		/* size nz, column indices of F */
    Fx = F->x ;		/* size nz, numeric values of F */
    values = values && (F->xtype != CHOLMOD_PATTERN) ;
    F->packed = TRUE ;
    F->sorted = TRUE ;
    F->stype = - SIGN (A->stype) ;	/* flip the stype */

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    Iwork = Common->Iwork ;
    Wi = Iwork ;	  /* size nrow. FUTURE WORK: Wi is (i/l/l) */
    Pinv = Iwork + nrow ; /* size nrow (i/i/l) , unused if Perm NULL */

    /* ---------------------------------------------------------------------- */
    /* construct the inverse permutation */
    /* ---------------------------------------------------------------------- */

    if (permute)
    {
	for (i = 0 ; i < n ; i++)
	{
	    Pinv [i] = EMPTY ;
	}
	for (k = 0 ; k < n ; k++)
	{
	    i = Perm [k] ;
	    if (i < 0 || i > n || Pinv [i] != EMPTY)
	    {
		cholmod_error (CHOLMOD_INVALID, "invalid permutation", Common) ;
		return (FALSE) ;
	    }
	    Pinv [i] = k ;
	}
    }

    /* Perm is now valid */
    ASSERT (cholmod_dump_perm (Perm, n, n, "Perm", Common)) ;

    /* ---------------------------------------------------------------------- */
    /* count the entries in each row of F */
    /* ---------------------------------------------------------------------- */

    for (i = 0 ; i < nrow ; i++)
    {
	Wi [i] = 0 ;
    }

    if (packed)
    {
	if (permute)
	{
	    if (upper)
	    {
		/* permuted, upper */
		for (j = 0 ; j < n ; j++)
		{
		    jold = Perm [j] ;
		    pend = Ap [jold+1] ;
		    for (p = Ap [jold] ; p < pend ; p++)
		    {
			if (Ai [p] <= jold)
			{
			    i = Pinv [Ai [p]] ;
			    Wi [MIN (i, j)]++ ;
			}
		    }
		}
	    }
	    else
	    {
		/* permuted, lower */
		for (j = 0 ; j < n ; j++)
		{
		    jold = Perm [j] ;
		    pend = Ap [jold+1] ;
		    for (p = Ap [jold] ; p < pend ; p++)
		    {
			if (Ai [p] >= jold)
			{
			    i = Pinv [Ai [p]] ;
			    Wi [MAX (i, j)]++ ;
			}
		    }
		}
	    }
	}
	else
	{
	    if (upper)
	    {
		/* unpermuted, upper */
		for (j = 0 ; j < n ; j++)
		{
		    pend = Ap [j+1] ;
		    for (p = Ap [j] ; p < pend ; p++)
		    {
			if (Ai [p] <= j)
			{
			    Wi [Ai [p]]++ ;
			}
		    }
		}
	    }
	    else
	    {
		/* unpermuted, lower */
		for (j = 0 ; j < n ; j++)
		{
		    pend = Ap [j+1] ;
		    for (p = Ap [j] ; p < pend ; p++)
		    {
			if (Ai [p] >= j)
			{
			    Wi [Ai [p]]++ ;
			}
		    }
		}
	    }
	}
    }
    else
    {
	if (permute)
	{
	    if (upper)
	    {
		/* permuted, upper */
		for (j = 0 ; j < n ; j++)
		{
		    jold = Perm [j] ;
		    p = Ap [jold] ;
		    pend = p + Anz [jold] ;
		    for ( ; p < pend ; p++)
		    {
			if (Ai [p] <= jold)
			{
			    i = Pinv [Ai [p]] ;
			    Wi [MIN (i, j)]++ ;
			}
		    }
		}
	    }
	    else
	    {
		/* permuted, lower */
		for (j = 0 ; j < n ; j++)
		{
		    jold = Perm [j] ;
		    p = Ap [jold] ;
		    pend = p + Anz [jold] ;
		    for ( ; p < pend ; p++)
		    {
			if (Ai [p] >= jold)
			{
			    i = Pinv [Ai [p]] ;
			    Wi [MAX (i, j)]++ ;
			}
		    }
		}
	    }
	}
	else
	{
	    if (upper)
	    {
		/* unpermuted, upper */
		for (j = 0 ; j < n ; j++)
		{
		    p = Ap [j] ;
		    pend = p + Anz [j] ;
		    for ( ; p < pend ; p++)
		    {
			if (Ai [p] <= j)
			{
			    Wi [Ai [p]]++ ;
			}
		    }
		}
	    }
	    else
	    {
		/* unpermuted, lower */
		for (j = 0 ; j < n ; j++)
		{
		    p = Ap [j] ;
		    pend = p + Anz [j] ;
		    for ( ; p < pend ; p++)
		    {
			if (Ai [p] >= j)
			{
			    Wi [Ai [p]]++ ;
			}
		    }
		}
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* compute the row pointers */
    /* ---------------------------------------------------------------------- */

    p = 0 ;
    for (i = 0 ; i < n ; i++)
    {
	Fp [i] = p ;
	p += Wi [i] ;
    }
    Fp [n] = p ;
    for (i = 0 ; i < n ; i++)
    {
	Wi [i] = Fp [i] ;
    }

    /* ---------------------------------------------------------------------- */
    /* construct the transpose */
    /* ---------------------------------------------------------------------- */

    if (packed)
    {
	if (permute)
	{
	    if (upper)
	    {
		if (values)
		{
		    /* packed, permuted, upper, values */
		    for (j = 0 ; j < n ; j++)
		    {
			jold = Perm [j] ;
			pend = Ap [jold+1] ;
			for (p = Ap [jold] ; p < pend ; p++)
			{
			    if (Ai [p] <= jold)
			    {
				i = Pinv [Ai [p]] ;
				if (i < j)
				{
				    fp = Wi [i]++ ;
				    Fj [fp] = j ;
				}
				else
				{
				    fp = Wi [j]++ ;
				    Fj [fp] = i ;
				}
				Fx [fp] = Ax [p] ;
			    }
			}
		    }
		}
		else
		{
		    /* packed, permuted, upper, no values */
		    for (j = 0 ; j < n ; j++)
		    {
			jold = Perm [j] ;
			pend = Ap [jold+1] ;
			for (p = Ap [jold] ; p < pend ; p++)
			{
			    if (Ai [p] <= jold)
			    {
				i = Pinv [Ai [p]] ;
				if (i < j)
				{
				    Fj [Wi [i]++] = j ;
				}
				else
				{
				    Fj [Wi [j]++] = i ;
				}
			    }
			}
		    }
		}
	    }
	    else
	    {
		if (values)
		{
		    /* packed, permuted, lower, values */
		    for (j = 0 ; j < n ; j++)
		    {
			jold = Perm [j] ;
			pend = Ap [jold+1] ;
			for (p = Ap [jold] ; p < pend ; p++)
			{
			    if (Ai [p] >= jold)
			    {
				i = Pinv [Ai [p]] ;
				if (i > j)
				{
				    fp = Wi [i]++ ;
				    Fj [fp] = j ;
				}
				else
				{
				    fp = Wi [j]++ ;
				    Fj [fp] = i ;
				}
				Fx [fp] = Ax [p] ;
			    }
			}
		    }
		}
		else
		{
		    /* packed, permuted, lower, no values */
		    for (j = 0 ; j < n ; j++)
		    {
			jold = Perm [j] ;
			pend = Ap [jold+1] ;
			for (p = Ap [jold] ; p < pend ; p++)
			{
			    if (Ai [p] >= jold)
			    {
				i = Pinv [Ai [p]] ;
				if (i > j)
				{
				    Fj [Wi [i]++] = j ;
				}
				else
				{
				    Fj [Wi [j]++] = i ;
				}
			    }
			}
		    }
		}
	    }
	}
	else
	{
	    if (upper)
	    {
		if (values)
		{
		    /* packed, unpermuted, upper, values */
		    for (j = 0 ; j < n ; j++)
		    {
			pend = Ap [j+1] ;
			for (p = Ap [j] ; p < pend ; p++)
			{
			    if (Ai [p] <= j)
			    {
				fp = Wi [Ai [p]]++ ;
				Fj [fp] = j ;
				Fx [fp] = Ax [p] ;
			    }
			}
		    }
		}
		else
		{
		    /* packed, unperm., upper, no values */
		    for (j = 0 ; j < n ; j++)
		    {
			pend = Ap [j+1] ;
			for (p = Ap [j] ; p < pend ; p++)
			{
			    if (Ai [p] <= j)
			    {
				Fj [Wi [Ai [p]]++] = j ;
			    }
			}
		    }
		}
	    }
	    else
	    {
		if (values)
		{
		    /* packed, unpermuted, lower, values */
		    for (j = 0 ; j < n ; j++)
		    {
			pend = Ap [j+1] ;
			for (p = Ap [j] ; p < pend ; p++)
			{
			    if (Ai [p] >= j)
			    {
				fp = Wi [Ai [p]]++ ;
				Fj [fp] = j ;
				Fx [fp] = Ax [p] ;
			    }
			}
		    }
		}
		else
		{
		    /* packed, unperm., lower, no values */
		    for (j = 0 ; j < n ; j++)
		    {
			pend = Ap [j+1] ;
			for (p = Ap [j] ; p < pend ; p++)
			{
			    if (Ai [p] >= j)
			    {
				Fj [Wi [Ai [p]]++] = j ;
			    }
			}
		    }
		}
	    }
	}
    }
    else
    {
	if (permute)
	{
	    if (upper)
	    {
		if (values)
		{
		    /* unpacked, permuted, upper, values */
		    for (j = 0 ; j < n ; j++)
		    {
			jold = Perm [j] ;
			p = Ap [jold] ;
			pend = p + Anz [jold] ;
			for ( ; p < pend ; p++)
			{
			    if (Ai [p] <= jold)
			    {
				i = Pinv [Ai [p]] ;
				if (i < j)
				{
				    fp = Wi [i]++ ;
				    Fj [fp] = j ;
				}
				else
				{
				    fp = Wi [j]++ ;
				    Fj [fp] = i ;
				}
				Fx [fp] = Ax [p] ;
			    }
			}
		    }
		}
		else
		{
		    /* unpacked, perm., upper, no values */
		    for (j = 0 ; j < n ; j++)
		    {
			jold = Perm [j] ;
			p = Ap [jold] ;
			pend = p + Anz [jold] ;
			for ( ; p < pend ; p++)
			{
			    if (Ai [p] <= jold)
			    {
				i = Pinv [Ai [p]] ;
				if (i < j)
				{
				    Fj [Wi [i]++] = j ;
				}
				else
				{
				    Fj [Wi [j]++] = i ;
				}
			    }
			}
		    }
		}
	    }
	    else
	    {
		if (values)
		{
		    /* unpacked, permuted, lower, values */
		    for (j = 0 ; j < n ; j++)
		    {
			jold = Perm [j] ;
			p = Ap [jold] ;
			pend = p + Anz [jold] ;
			for ( ; p < pend ; p++)
			{
			    if (Ai [p] >= jold)
			    {
				i = Pinv [Ai [p]] ;
				if (i > j)
				{
				    fp = Wi [i]++ ;
				    Fj [fp] = j ;
				}
				else
				{
				    fp = Wi [j]++ ;
				    Fj [fp] = i ;
				}
				Fx [fp] = Ax [p] ;
			    }
			}
		    }
		}
		else
		{
		    /* unpacked, perm., lower, no values */
		    for (j = 0 ; j < n ; j++)
		    {
			jold = Perm [j] ;
			p = Ap [jold] ;
			pend = p + Anz [jold] ;
			for ( ; p < pend ; p++)
			{
			    if (Ai [p] >= jold)
			    {
				i = Pinv [Ai [p]] ;
				if (i > j)
				{
				    Fj [Wi [i]++] = j ;
				}
				else
				{
				    Fj [Wi [j]++] = i ;
				}
			    }
			}
		    }
		}
	    }
	}
	else
	{
	    if (upper)
	    {
		if (values)
		{
		    /* unpacked, unperm., upper, values */
		    for (j = 0 ; j < n ; j++)
		    {
			p = Ap [j] ;
			pend = p + Anz [j] ;
			for ( ; p < pend ; p++)
			{
			    if (Ai [p] <= j)
			    {
				fp = Wi [Ai [p]]++ ;
				Fj [fp] = j ;
				Fx [fp] = Ax [p] ;
			    }
			}
		    }
		}
		else
		{
		    /* unpacked, unperm, upper, no values */
		    for (j = 0 ; j < n ; j++)
		    {
			p = Ap [j] ;
			pend = p + Anz [j] ;
			for ( ; p < pend ; p++)
			{
			    if (Ai [p] <= j)
			    {
				Fj [Wi [Ai [p]]++] = j ;
			    }
			}
		    }
		}
	    }
	    else
	    {
		if (values)
		{
		    /* unpacked, unperm., lower, values */
		    for (j = 0 ; j < n ; j++)
		    {
			p = Ap [j] ;
			pend = p + Anz [j] ;
			for ( ; p < pend ; p++)
			{
			    if (Ai [p] >= j)
			    {
				fp = Wi [Ai [p]]++ ;
				Fj [fp] = j ;
				Fx [fp] = Ax [p] ;
			    }
			}
		    }
		}
		else
		{
		    /* unpacked, unperm, lower, no values */
		    for (j = 0 ; j < n ; j++)
		    {
			p = Ap [j] ;
			pend = p + Anz [j] ;
			for ( ; p < pend ; p++)
			{
			    if (Ai [p] >= j)
			    {
				Fj [Wi [Ai [p]]++] = j ;
			    }
			}
		    }
		}
	    }
	}
    }

    /* F is sorted if there is no permutation vector */
    F->sorted = !permute ;

    ASSERT (cholmod_dump_sparse (F, "output F sym", Common) >= 0) ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_transpose ==================================================== */
/* ========================================================================== */

/* workspace:
 * Iwork (MAX (nrow,ncol)) if unsymmetric and fset is non-NULL
 * Iwork (nrow) if unsymmetric and fset is NULL
 * Iwork (2*nrow) if symmetric and Perm is non-NULL.
 * Iwork (nrow) if symmetric and Perm is NULL.
 *
 * A simple worst-case upper bound on the workspace is nrow+ncol.
 */

cholmod_sparse *cholmod_transpose
(
    /* inputs, not modified: */
    cholmod_sparse *A,	/* nrow-by-ncol */
    int values,		/* TRUE if values to be transposed, else FALSE */
    void *Perm_p,	/* if non-NULL, F = A(p,f)' or A(p,p)' */
    void *fset_p,	/* if non-NULL, F = A(*,f)' */
    size_t fsize,	/* fset and fsize are used for unsymmetric case only */

    cholmod_common *Common  /* Iwork (nrow, or 2*nrow if Perm not NULL) */
)
{
    int *Ap, *Anz, *Perm, *fset ;
    cholmod_sparse *F ;
    int nrow, ncol, use_fset, j, jj, fnz, packed, symmetric, nf, ineed, ok ;
    DEBUG (int orig) ;

    Perm = Perm_p ;
    fset = fset_p ;
    nf = fsize ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, FALSE) ;
    symmetric = A->stype ;
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    nrow = A->nrow ;
    ncol = A->ncol ;

    if (symmetric)
    {
	use_fset = FALSE ;
	if (Perm != NULL)
	{
	    ineed = 2*nrow ;
	}
	else
	{
	    ineed = nrow ;
	}
    }
    else
    {
	use_fset = (fset != NULL) ;
	if (use_fset)
	{
	    ineed = MAX (nrow, ncol) ;
	}
	else
	{
	    ineed = nrow ;
	}
    }

    cholmod_allocate_work (0, ineed, 0, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory */
    }

    DEBUG (orig = Common->malloc_count) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;
    Anz = A->nz ;
    packed = A->packed ;
    ASSERT (IMPLIES (!packed, Anz != NULL)) ;
    values = values && (A->xtype != CHOLMOD_PATTERN) ;

    /* ---------------------------------------------------------------------- */
    /* allocate F */
    /* ---------------------------------------------------------------------- */

    /* determine # of nonzeros in F */
    if (symmetric)
    {
	/* F=A' or F=A(p,p)', fset is ignored */
	fnz = cholmod_nnz (A, Common) ;
    }
    else
    {
	nf = (use_fset) ? nf : ncol ;
	if (use_fset)
	{
	    fnz = 0 ;
	    /* F=A(:,f)' or F=A(p,f)' */
	    for (jj = 0 ; jj < nf ; jj++)
	    {
		/* The fset is not yet checked; it will be thoroughly checked
		 * in cholmod_transpose_unsym.  For now, just make sure we don't
		 * access Ap and Anz out of bounds. */
		j = fset [jj] ;
		if (j >= 0 && j < ncol)
		{
		    fnz += packed ? (Ap [j+1] - Ap [j]) : MAX (0, Anz [j]) ;
		}
	    }
	}
	else
	{
	    /* F=A' or F=A(p,:)' */
	    fnz = cholmod_nnz (A, Common) ;
	}
    }

    /* F is ncol-by-nrow, fnz nonzeros, sorted unless f is present and unsorted,
     * packed, of opposite stype as A, and with/without numerical values */
    F = cholmod_allocate_sparse (ncol, nrow, fnz, TRUE, TRUE, -SIGN(symmetric),
	    values, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	ASSERT (Common->malloc_count == orig) ;
	return (NULL) ;	    /* out of memory */
    }

    /* ---------------------------------------------------------------------- */
    /* transpose and optionally permute the matrix A */
    /* ---------------------------------------------------------------------- */

    if (symmetric)
    {
	/* F = A (p,p)', using upper or lower triangular part of A only */
	ok = cholmod_transpose_sym (A, values, Perm, F, Common) ;
    }
    else
    {
	/* F = A (p,f)' */
	ok = cholmod_transpose_unsym (A, values, Perm, fset, nf, F, Common) ;
    }

    /* ---------------------------------------------------------------------- */
    /* return the matrix F, or NULL if an error occured */
    /* ---------------------------------------------------------------------- */

    if (!ok)
    {
	cholmod_free_sparse (&F, Common) ;
    }
    ASSERT (Common->malloc_count == orig + (ok ? ((values) ? 4 : 3) : 0)) ;
    return (F) ;
}


/* ========================================================================== */
/* === cholmod_sort ========================================================= */
/* ========================================================================== */

/* Sort the columns of A, in place.  Returns A in packed form, even if it
 * starts as unpacked.  Removes entries in the ignored part of a symmetric
 * matrix.
 *
 * FUTURE WORK: this could be faster.  The 2nd transpose does not need the row
 * counts of its input (the col counts of the final output).  These can be
 * determined directly from the original matrix.
 *
 * workspace: Iwork (max (nrow,ncol)).  Allocates additional workspace for a
 * temporary copy of A'.
 */

int cholmod_sort
(
    cholmod_sparse *A,
    cholmod_common *Common
)
{
    int *Ap ;
    cholmod_sparse *F ;
    int anz, ncol, nrow, symmetric, values ;
    DEBUG (int orig) ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    Common->status = CHOLMOD_OK ;
    nrow = A->nrow ;
    if (nrow <= 1)
    {
	/* a 1-by-n sparse matrix must be sorted */
	A->sorted = TRUE ;
	return (TRUE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    ncol = A->ncol ;
    cholmod_allocate_work (0, MAX (nrow, ncol), 0, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;	/* out of memory */
    }

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    anz = cholmod_nnz (A, Common) ;
    values = (A->xtype != CHOLMOD_PATTERN) ;
    symmetric = A->stype ;
    DEBUG (orig = Common->malloc_count) ;

    /* ---------------------------------------------------------------------- */
    /* sort the columns of the matrix */
    /* ---------------------------------------------------------------------- */

    /* allocate workspace for transpose: ncol-by-nrow, same # of nonzeros as A,
     * sorted, packed, same stype as A, and of the same numeric type as A. */
    F = cholmod_allocate_sparse (ncol, nrow, anz, TRUE, TRUE, symmetric,
	    values, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	ASSERT (Common->malloc_count == orig) ;
	return (FALSE) ;	/* out of memory */
    }

    if (symmetric)
    {
	/* F = A', upper or lower triangular part only */
	cholmod_transpose_sym (A, values, NULL, F, Common) ;
	A->packed = TRUE ;
	/* A = F' */
	cholmod_transpose_sym (F, values, NULL, A, Common) ;
    }
    else
    {
	/* F = A' */
	cholmod_transpose_unsym (A, values, NULL, NULL, 0, F, Common) ;
	A->packed = TRUE ;
	/* A = F' */
	cholmod_transpose_unsym (F, values, NULL, NULL, 0, A, Common) ;
    }

    ASSERT (A->sorted && A->packed) ;
    ASSERT (cholmod_dump_sparse (A, "Asorted", Common) >= 0) ;

    /* ---------------------------------------------------------------------- */
    /* reduce A in size, if needed.  This must succeed. */
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;
    anz = Ap [ncol] ;
    ASSERT ((size_t) anz <= A->nzmax) ;
    (void) cholmod_reallocate_sparse (A, anz, Common) ;
    ASSERT (Common->status >= CHOLMOD_OK) ;

    /* ---------------------------------------------------------------------- */
    /* free workspace */
    /* ---------------------------------------------------------------------- */

    cholmod_free_sparse (&F, Common) ;
    ASSERT (Common->malloc_count == orig) ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_copy_sparse ================================================== */
/* ========================================================================== */

/* Create an exact copy of a sparse matrix, with two exceptions:
 *
 * (1) Entries in unused space are not copied (they might not be initialized,
 *	and copying them would cause program checkers such as purify and
 *	valgrind to complain).
 *
 * (2) C->xtype becomes CHOLMOD_PATTERN and C->x becomes NULL if either A->x
 *	is NULL or A->xtype is CHOLMOD_PATTERN.  The two logical statements
 *	(A->x == NULL) and (A->xtype == CHOLMOD_PATTERN) are kept logically
 *	equivalent by CHOLMOD, but a user might modify one without making a
 *	change to the other, outside of any CHOLMOD routine.
 *
 * See also MatrixOps/cholmod_copy, which copies a matrix with possible changes
 * in stype, presence of diagonal entries, pattern vs. numerical values,
 * real and/or imaginary parts, and so on.
 */

cholmod_sparse *cholmod_copy_sparse
(
    cholmod_sparse *A,
    cholmod_common *Common
)
{
    double *Ax, *Cx ;
    int *Ap, *Ai, *Anz, *Cp, *Ci, *Cnz ;
    cholmod_sparse *C ;
    int p, pend, j, ncol, packed, nzmax, nz, values ;
    DEBUG (int orig) ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    DEBUG (orig = Common->malloc_count) ;
    if (A->stype != 0 && A->nrow != A->ncol)
    {
	cholmod_error (CHOLMOD_INVALID,
		"rectangular matrix with stype != 0 invalid", Common) ;
	return (NULL) ;
    }
    Common->status = CHOLMOD_OK ;
    ASSERT (cholmod_dump_sparse (A, "A original", Common) >= 0) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    ncol = A->ncol ;
    nzmax = A->nzmax ;
    packed = A->packed ;
    Ap = A->p ;
    Ai = A->i ;
    Ax = A->x ;
    Anz = A->nz ;
    values = (Ax != NULL && A->xtype != CHOLMOD_PATTERN) ;

    /* ---------------------------------------------------------------------- */
    /* allocate the copy */
    /* ---------------------------------------------------------------------- */

    C = cholmod_allocate_sparse (A->nrow, A->ncol, A->nzmax, A->sorted,
	    A->packed, A->stype, values, Common) ;

    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory */
    }

    Cp = C->p ;
    Ci = C->i ;
    Cx = C->x ;
    Cnz = C->nz ;

    /* ---------------------------------------------------------------------- */
    /* copy the matrix */
    /* ---------------------------------------------------------------------- */

    for (j = 0 ; j <= ncol ; j++)
    {
	Cp [j] = Ap [j] ;
    }

    if (packed)
    {
	nz = Ap [ncol] ;
	for (p = 0 ; p < nz ; p++)
	{
	    Ci [p] = Ai [p] ;
	}
	if (values)
	{
	    for (p = 0 ; p < nz ; p++)
	    {
		Cx [p] = Ax [p] ;
	    }
	}
    }
    else
    {
	for (j = 0 ; j < ncol ; j++)
	{
	    Cnz [j] = Anz [j] ;
	}
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = p + Anz [j] ;
	    for ( ; p < pend ; p++)
	    {
		Ci [p] = Ai [p] ;
		if (values)
		{
		    Cx [p] = Ax [p] ;
		}
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* return the result */
    /* ---------------------------------------------------------------------- */

    ASSERT (Common->malloc_count == orig +
	    3 + (values ? 1 : 0) + ((packed) ? 0 : 1)) ;
    ASSERT (cholmod_dump_sparse (C, "C copy", Common) >= 0) ;
    return (C) ;
}
