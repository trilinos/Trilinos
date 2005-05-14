/* ========================================================================== */
/* === cholmod_add ========================================================== */
/* ========================================================================== */

/*
 * CHOLMOD/Core version 0.1. May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* C = alpha*A + beta*B, or spones(A+B).  Result is packed, with sorted or
 * unsorted columns.  This routine is much faster and takes less memory if C
 * is allowed to have unsorted columns.
 *
 * If A and B are both symmetric (in upper form) then C is the same.  Likewise,
 * if A and B are both symmetric (in lower form) then C is the same.
 * Otherwise, C is unsymmetric.  A and B must have the same dimension.
 *
 * workspace: Flag (nrow), W (nrow) if values, Iwork (max (nrow,ncol)).
 *	allocates temporary copies for A and B if they are symmetric.
 *	allocates temporary copy of C if it is to be returned sorted.
 */

#include "cholmod_core.h"
#include "cholmod_internal.h"

cholmod_sparse *cholmod_add		/* returns C, or NULL if failure */
(
    /* inputs, not modified on output */
    cholmod_sparse *A,
    cholmod_sparse *B,
    cholmod_scalar alpha,
    cholmod_scalar beta,
    int values,		    /* if TRUE compute the numerical values of C */
    int sorted,		    /* if TRUE then return C with sorted columns */

    cholmod_common *Common
)
{
    double *Ax, *Bx, *Cx, *W ;
    int apacked, up, lo, nrow, ncol, bpacked, nzmax, pa, paend, pb, pbend, i,
	j, p, mark, nz ;
    int *Ap, *Ai, *Anz, *Bp, *Bi, *Bnz, *Flag, *Cp, *Ci ;
    cholmod_sparse *A2, *B2, *C ;

    DEBUG (int orig) ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    RETURN_IF_NULL (B, NULL) ;
    nrow = A->nrow ;
    ncol = A->ncol ;
    if (nrow != (int) (B->nrow) || ncol != (int) (B->ncol))
    {
	/* A and B must have the same dimensions */
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_add: A and B dimesions do not match", Common) ;
	return (NULL) ;
    }
    values = values &&
	(A->xtype != CHOLMOD_PATTERN) && (B->xtype != CHOLMOD_PATTERN) ;

    if (values && A->xtype != B->xtype)
    {
	/* A and B must have the same numerical type */
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_add: A and B must have same numerical type", Common) ;
	return (NULL) ;
    }

    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    cholmod_allocate_work (nrow, MAX (nrow,ncol), values ? nrow : 0,
	    sizeof (double), Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory */
    }

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    if (nrow <= 1)
    {
	/* C will be implicitly sorted, so no need to sort it here */
	sorted = FALSE ;
    }

    DEBUG (orig = Common->malloc_count) ;

    /* convert A or B to unsymmetric, if necessary */
    A2 = NULL ;
    B2 = NULL ;

    if (A->stype != B->stype)
    {
	if (A->stype)
	{
	    /* workspace: Iwork (max (nrow,ncol)) */
	    A2 = cholmod_copy (A, 0, values, Common) ;
	    if (Common->status < CHOLMOD_OK)
	    {
		ASSERT (Common->malloc_count == orig) ;
		return (NULL) ;	    /* out of memory */
	    }
	    A = A2 ;
	}
	if (B->stype)
	{
	    /* workspace: Iwork (max (nrow,ncol)) */
	    B2 = cholmod_copy (B, 0, values, Common) ;
	    if (Common->status < CHOLMOD_OK)
	    {
		cholmod_free_sparse (&A2, Common) ;
		ASSERT (Common->malloc_count == orig) ;
		return (NULL) ;	    /* out of memory */
	    }
	    B = B2 ;
	}
    }

    /* get the A matrix */
    ASSERT (A->stype == B->stype) ;
    up = (A->stype > 0) ;
    lo = (A->stype < 0) ;

    Ap  = A->p ;
    Anz = A->nz ;
    Ai  = A->i ;
    Ax  = A->x ;
    apacked = A->packed ;

    /* get the B matrix */
    Bp  = B->p ;
    Bnz = B->nz ;
    Bi  = B->i ;
    Bx  = B->x ;
    bpacked = B->packed ;

    /* get workspace */
    W = Common->Xwork ;	    /* size nrow, used if values is TRUE */
    Flag = Common->Flag ;   /* size nrow, Flag [0..nrow-1] < mark on input */

    /* ---------------------------------------------------------------------- */
    /* allocate the result C */
    /* ---------------------------------------------------------------------- */

    /* If int overflow occurs, nzmax < 0 and the allocate fails properly
     * (likewise in most other matrix manipulation routines). */
    nzmax = A->nzmax + B->nzmax ;
    C = cholmod_allocate_sparse (nrow, ncol, nzmax, FALSE, TRUE,
	    SIGN (A->stype), values, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	cholmod_free_sparse (&A2, Common) ;
	cholmod_free_sparse (&B2, Common) ;
	ASSERT (Common->malloc_count == orig) ;
	return (NULL) ;	    /* out of memory */
    }

    Cp = C->p ;
    Ci = C->i ;
    Cx = C->x ;

    /* ---------------------------------------------------------------------- */
    /* compute C = alpha*A + beta*B */
    /* ---------------------------------------------------------------------- */

    nz = 0 ;
    for (j = 0 ; j < ncol ; j++)
    {
	Cp [j] = nz ;

	/* clear the Flag array */
	mark = cholmod_clear_flag (Common) ;

	/* scatter B into W */
	pb = Bp [j] ;
	pbend = (bpacked) ? (Bp [j+1]) : (pb + Bnz [j]) ;
	for (p = pb ; p < pbend ; p++)
	{
	    i = Bi [p] ;
	    if ((up && i > j) || (lo && i < j))
	    {
		continue ;
	    }
	    Flag [i] = mark ;
	    if (values)
	    {
		W [i] = beta.x * Bx [p] ;
	    }
	}

	/* add A and gather from W into C(:,j) */
	pa = Ap [j] ;
	paend = (apacked) ? (Ap [j+1]) : (pa + Anz [j]) ;
	for (p = pa ; p < paend ; p++)
	{
	    i = Ai [p] ;
	    if ((up && i > j) || (lo && i < j))
	    {
		continue ;
	    }
	    Flag [i] = EMPTY ;
	    Ci [nz] = i ;
	    if (values)
	    {
		Cx [nz] = W [i] + alpha.x * Ax [p] ;
		W [i] = 0 ;
	    }
	    nz++ ;
	}

	/* gather remaining entries into C(:,j), using pattern of B */
	for (p = pb ; p < pbend ; p++)
	{
	    i = Bi [p] ;
	    if ((up && i > j) || (lo && i < j))
	    {
		continue ;
	    }
	    if (Flag [i] == mark)
	    {
		Ci [nz] = i ;
		if (values)
		{
		    Cx [nz] = W [i] ;
		    W [i] = 0 ;
		}
		nz++ ;
	    }
	}
    }

    Cp [ncol] = nz ;

    /* ---------------------------------------------------------------------- */
    /* reduce C in size and free temporary matrices */
    /* ---------------------------------------------------------------------- */

    ASSERT (MAX (1,nz) <= C->nzmax) ;
    (void) cholmod_reallocate_sparse (C, nz, Common) ;
    ASSERT (Common->status >= CHOLMOD_OK) ;

    /* clear the Flag array */
    mark = cholmod_clear_flag (Common) ;

    cholmod_free_sparse (&A2, Common) ;
    cholmod_free_sparse (&B2, Common) ;

    /* ---------------------------------------------------------------------- */
    /* sort C, if requested */
    /* ---------------------------------------------------------------------- */

    if (sorted)
    {
	/* workspace: Iwork (max (nrow,ncol)) */
	if (!cholmod_sort (C, Common))
	{
	    cholmod_free_sparse (&C, Common) ;
	    if (Common->status < CHOLMOD_OK)
	    {
		ASSERT (Common->malloc_count == orig) ;
		return (NULL) ;		/* out of memory */
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* return result */
    /* ---------------------------------------------------------------------- */

    ASSERT (Common->malloc_count == orig + ((values) ? 4 : 3)) ;
    ASSERT (cholmod_dump_sparse (C, "add", Common) >= 0) ;
    return (C) ;
}
