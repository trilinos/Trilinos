/* ========================================================================== */
/* === Core/cholmod_aat ===================================================== */
/* ========================================================================== */

/*
 * CHOLMOD/Core version 0.1. May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* C = A*A' or C = A(:,f)*A(:,f)'
 *
 * A can be packed or unpacked, sorted or unsorted, but must be stored with
 * both upper and lower parts.  C is returned as packed, "unsymmetric" (both
 * upper and lower parts present), and unsorted.  See cholmod_ssmult in the
 * MatrixOps package for a more general matrix-matrix multiply.
 *
 * workspace:
 *	Flag (A->nrow)
 *	Iwork (max (A->nrow, A->ncol)) if fset
 *	Iwork (A->nrow) if no fset
 *	W (A->nrow) if mode > 0
 *	allocates temporary copy for A'.
 */

#include "cholmod_core.h"
#include "cholmod_internal.h"

cholmod_sparse *cholmod_aat
(
    /* inputs, not modified on output */
    cholmod_sparse *A,
    void *fset,
    size_t fsize,
    int mode,		    /* >0: compute numerical values of C,
			       0:  pattern of C only,
			       <0: pattern of C, excluding the diagonal */
    cholmod_common *Common
)
{
    double fjt ;
    double *Ax, *Fx, *Cx, *W ;
    int *Ap, *Anz, *Ai, *Fp, *Fi, *Cp, *Ci, *Flag ;
    cholmod_sparse *C, *F ;
    int packed, j, i, pa, paend, pf, pfend, n, mark, cnz, t, p, values, diag ;

    DEBUG (int orig) ;
    DEBUG (int nnzdiag) ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    if (A->stype)
    {
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_aat: matrix cannot be symmetric", Common) ;
	return (NULL) ;
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    values = (mode > 0) && (A->xtype != CHOLMOD_PATTERN) ;
    diag = (mode >= 0) ;
    n = A->nrow ;
    cholmod_allocate_work (n, MAX (A->ncol, A->nrow), values ? n : 0,
	    sizeof (double), Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory */
    }
    ASSERT (cholmod_dump_work (TRUE, TRUE, values ? n : 0, Common)) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    DEBUG (orig = Common->malloc_count) ;
    ASSERT (cholmod_dump_sparse (A, "A", Common) >= 0) ;

    /* get the A matrix */
    Ap  = A->p ;
    Anz = A->nz ;
    Ai  = A->i ;
    Ax  = A->x ;
    packed = A->packed ;

    /* get workspace */
    W = Common->Xwork ;		/* size n, unused if values is FALSE */
    Flag = Common->Flag ;	/* size n, Flag [0..n-1] < mark on input*/

    /* ---------------------------------------------------------------------- */
    /* F = A' or A(:,f)' */
    /* ---------------------------------------------------------------------- */

    /* workspace: Iwork (nrow if no fset; MAX (nrow,ncol) if fset)*/
    F = cholmod_transpose (A, values, NULL, fset, fsize, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	ASSERT (Common->malloc_count == orig) ;
	return (NULL) ;	    /* out of memory */
    }

    Fp = F->p ;
    Fi = F->i ;
    Fx = F->x ;

    /* ---------------------------------------------------------------------- */
    /* count the number of entries in the result C */
    /* ---------------------------------------------------------------------- */

    cnz = 0 ;
    for (j = 0 ; j < n ; j++)
    {
	/* clear the Flag array */
	mark = cholmod_clear_flag (Common) ;

	/* exclude the diagonal, if requested */
	if (!diag)
	{
	    Flag [j] = mark ;
	}

	/* for each nonzero F(t,j) in column j, do: */
	pfend = Fp [j+1] ;
	for (pf = Fp [j] ; pf < pfend ; pf++)
	{
	    /* F(t,j) is nonzero */
	    t = Fi [pf] ;

	    /* add the nonzero pattern of A(:,t) to the pattern of C(:,j) */
	    pa = Ap [t] ;
	    paend = (packed) ? (Ap [t+1]) : (pa + Anz [t]) ;
	    for ( ; pa < paend ; pa++)
	    {
		i = Ai [pa] ;
		if (Flag [i] != mark)
		{
		    Flag [i] = mark ;
		    cnz++ ;
		}
	    }
	}
	if (cnz < 0)
	{
	    break ;	    /* int overflow case */
	}
    }

    mark = cholmod_clear_flag (Common) ;

    /* ---------------------------------------------------------------------- */
    /* check for int overflow */
    /* ---------------------------------------------------------------------- */

    if (cnz < 0)
    {
	cholmod_error (CHOLMOD_TOO_LARGE, "problem too large", Common) ;
	(void) cholmod_clear_flag (Common) ;
	cholmod_free_sparse (&F, Common) ;
	return (NULL) ;	    /* problem too large */
    }

    /* ---------------------------------------------------------------------- */
    /* allocate C */
    /* ---------------------------------------------------------------------- */

    C = cholmod_allocate_sparse (n, n, cnz, FALSE, TRUE, 0, values, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	cholmod_free_sparse (&F, Common) ;
	ASSERT (Common->malloc_count == orig) ;
	return (NULL) ;	    /* out of memory */
    }

    Cp = C->p ;
    Ci = C->i ;
    Cx = C->x ;

    /* ---------------------------------------------------------------------- */
    /* C = A*A' */
    /* ---------------------------------------------------------------------- */

    cnz = 0 ;

    if (values)
    {

	/* pattern and values */
	for (j = 0 ; j < n ; j++)
	{
	    /* clear the Flag array */
	    mark = cholmod_clear_flag (Common) ;

	    /* start column j of C */
	    Cp [j] = cnz ;

	    /* for each nonzero F(t,j) in column j, do: */
	    pfend = Fp [j+1] ;
	    for (pf = Fp [j] ; pf < pfend ; pf++)
	    {
		/* F(t,j) is nonzero */
		t = Fi [pf] ;
		fjt = Fx [pf] ;

		/* add the nonzero pattern of A(:,t) to the pattern of C(:,j)
		 * and scatter the values into W */
		pa = Ap [t] ;
		paend = (packed) ? (Ap [t+1]) : (pa + Anz [t]) ;
		for ( ; pa < paend ; pa++)
		{
		    i = Ai [pa] ;
		    if (Flag [i] != mark)
		    {
			Flag [i] = mark ;
			Ci [cnz++] = i ;
		    }
		    W [i] += Ax [pa] * fjt ;
		}
	    }

	    /* gather the values into C(:,j) */
	    for (p = Cp [j] ; p < cnz ; p++)
	    {
		i = Ci [p] ;
		Cx [p] = W [i] ;
		W [i] = 0 ;
	    }
	}

    }
    else
    {

	/* pattern only */
	for (j = 0 ; j < n ; j++)
	{
	    /* clear the Flag array */
	    mark = cholmod_clear_flag (Common) ;

	    /* exclude the diagonal, if requested */
	    if (!diag)
	    {
		Flag [j] = mark ;
	    }

	    /* start column j of C */
	    Cp [j] = cnz ;

	    /* for each nonzero F(t,j) in column j, do: */
	    pfend = Fp [j+1] ;
	    for (pf = Fp [j] ; pf < pfend ; pf++)
	    {
		/* F(t,j) is nonzero */
		t = Fi [pf] ;

		/* add the nonzero pattern of A(:,t) to the pattern of C(:,j) */
		pa = Ap [t] ;
		paend = (packed) ? (Ap [t+1]) : (pa + Anz [t]) ;
		for ( ; pa < paend ; pa++)
		{
		    i = Ai [pa] ;
		    if (Flag [i] != mark)
		    {
			Flag [i] = mark ;
			Ci [cnz++] = i ;
		    }
		}
	    }
	}
    }

    Cp [n] = cnz ;
    ASSERT (MAX (1,cnz) == C->nzmax) ;

    /* ---------------------------------------------------------------------- */
    /* clear workspace and free temporary matrices and return result */
    /* ---------------------------------------------------------------------- */

    cholmod_free_sparse (&F, Common) ;
    cholmod_clear_flag (Common) ;
    ASSERT (cholmod_dump_work (TRUE, TRUE, values ? n : 0, Common)) ;
    ASSERT (Common->malloc_count == orig + ((values) ? 4 : 3)) ;
    DEBUG (nnzdiag = cholmod_dump_sparse (C, "aat", Common)) ;
    ASSERT (IMPLIES (mode < 0, nnzdiag == 0)) ;
    return (C) ;
}
