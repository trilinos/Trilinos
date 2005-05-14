/* ========================================================================== */
/* === Cholesky/cholmod_solve =============================================== */
/* ========================================================================== */

/*
 * CHOLMOD/Cholesky version 0.1, May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Authors: Timothy A. Davis
 * and Sivasankaran Rajamanickam
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* Solve one of the following systems:
 *
 *	Ax=b	    0: CHOLMOD_A	also applies the permutation L->Perm
 *	LDL'x=b	    1: CHOLMOD_LDLt	does not apply L->Perm
 *	LDx=b	    2: CHOLMOD_LD
 *	DL'x=b	    3: CHOLMOD_DLt
 *	Lx=b	    4: CHOLMOD_L
 *	L'x=b	    5: CHOLMOD_Lt
 *	Dx=b	    6: CHOLMOD_D
 *	x=Pb	    7: CHOLMOD_P	apply a permutation (P is L->Perm)
 *	x=P'b	    8: CHOLMOD_Pt	apply an inverse permutation
 *
 * The factorization can be simplicial LDL', simplicial LL', or supernodal LL'.
 * For an LL' factorization, D is the identity matrix.  Thus CHOLMOD_LD and
 * CHOLMOD_L solve the same system, for example.
 *
 * workspace: uses either Common->Xwork, or temporary space which is freed when
 * done.  If the space required is <= Common->maxrank * nrhs, then Common->Xwork
 * is used.  It's faster and less stressful to your memory space if Xwork is
 * used, because multiple solves can share the same workspace.
 *
 *	Solving Dx=b, x=Pb, or x=P'b never requires temporary workspace.
 *
 *	Solving Ax=b, Lx=b, LDx=b, L'x=b, or DL'x=b with a supernodal
 *	    factorization: the space required is nrhs * L->maxesize, where b is
 *	    n-by-nrhs and L->maxesize is bounded above by n.   Thus, if nrhs is
 *	    less than or equal to maxrank, Xwork is always used.
 *
 *	Solving Ax=b, Lx=b, or L'x=b with a simplicial factorization, the space
 *	    required is MIN (nrhs,4) * n.  Thus, if maxrank is 4 or greater
 *	    (or greater than maxrank), Xwork is always used.  If you use a lot
 *	    of simplicial solves, choosing a maxrank of 4 is a good option.
 *
 * The default value of Common->maxrank is 2, so that if you only solve systems
 * where b has at most two columns, no temporary space will ever be allocated.
 * Valid values of Common->maxrank are 2, 4, and 8.  Values outside that
 * range are converted to that range.
 *
 * The supernodal solver uses BLAS routines dtrsv, dgemv, dtrsm, and dgemm.
 */

#include "cholmod_cholesky.h"
#include "cholmod_internal.h"

#ifndef NSUPERNODAL
#include "cholmod_supernodal.h"
#endif

/* ========================================================================== */
/* === ptrans =============================================================== */
/* ========================================================================== */

/* x = P*b, permute and transpose the column-major right-hand-side into the
 * row-major work vector x.  Operates on at most 4 columns at a time.
 */

static void ptrans
(
    int P [ ],
    int n,
    double X [ ],		    /* n-by-nr in row-major order */
    double B [ ],		    /* n-by-nr in col-major order */
    int d,			    /* leading dimension of B */
    int nr			    /* nr = 1, 2, 3, or 4 */
)
{
    int i, k ;

    switch (nr)
    {
	case 1:
	    for (k = 0 ; k < n ; k++)
	    {
		X [k] = B [P [k]] ;
	    }
	    break ;

	case 2:
	    for (k = 0 ; k < n ; k++)
	    {
		i = P [k] ;
		X [2*k    ] = B [i      ] ;
		X [2*k + 1] = B [i + d  ] ;
	    }
	    break ;

	case 3:
	    for (k = 0 ; k < n ; k++)
	    {
		i = P [k] ;
		X [3*k    ] = B [i      ] ;
		X [3*k + 1] = B [i + d  ] ;
		X [3*k + 2] = B [i + d*2] ;
	    }
	    break ;

	case 4:
	    for (k = 0 ; k < n ; k++)
	    {
		i = P [k] ;
		X [4*k    ] = B [i      ] ;
		X [4*k + 1] = B [i + d  ] ;
		X [4*k + 2] = B [i + d*2] ;
		X [4*k + 3] = B [i + d*3] ;
	    }
	    break ;

    }
}


/* ========================================================================== */
/* === trans ================================================================ */
/* ========================================================================== */

/* x = b, transpose column-major right-hand-side into the row-major work vector
 * x.  Operates on 2, 3, or 4 columns at a time.
 */

static void trans
(
    int n,
    double X [ ],		    /* n-by-nr in row-major order */
    double B [ ],		    /* n-by-nr in col-major order */
    int d,			    /* leading dimension of B */
    int nr			    /* nr = 2, 3, or 4 */
)
{
    int k ;

    switch (nr)
    {

	case 2:
	    for (k = 0 ; k < n ; k++)
	    {
		X [2*k    ] = B [k      ] ;
		X [2*k + 1] = B [k + d  ] ;
	    }
	    break ;

	case 3:
	    for (k = 0 ; k < n ; k++)
	    {
		X [3*k    ] = B [k      ] ;
		X [3*k + 1] = B [k + d  ] ;
		X [3*k + 2] = B [k + d*2] ;
	    }
	    break ;

	case 4:
	    for (k = 0 ; k < n ; k++)
	    {
		X [4*k    ] = B [k      ] ;
		X [4*k + 1] = B [k + d  ] ;
		X [4*k + 2] = B [k + d*2] ;
		X [4*k + 3] = B [k + d*3] ;
	    }
	    break ;

    }
}


/* ========================================================================== */
/* === inv_ptrans =========================================================== */
/* ========================================================================== */

/* b = P'*x, permute and transpose the row-major work vector x back into the
 * column-major right-hand-side.  Operates on at most 4 columns at a time.
 */

static void inv_ptrans
(
    int P [ ],
    int n,
    double X [ ],		    /* n-by-nr in row-major order */
    double B [ ],		    /* n-by-nr in col-major order */
    int d,			    /* leading dimension of B */
    int nr			    /* nr = 1, 2, 3, or 4 */
)
{
    int i, k ;

    switch (nr)
    {

	case 1:
	    for (k = 0 ; k < n ; k++)
	    {
		B [P [k]] = X [k] ;
	    }
	    break ;

	case 2:
	    for (k = 0 ; k < n ; k++)
	    {
		i = P [k] ;
		B [i      ] = X [2*k    ] ;
		B [i + d  ] = X [2*k + 1] ;
	    }
	    break ;

	case 3:
	    for (k = 0 ; k < n ; k++)
	    {
		i = P [k] ;
		B [i      ] = X [3*k    ] ;
		B [i + d  ] = X [3*k + 1] ;
		B [i + d*2] = X [3*k + 2] ;
	    }
	    break ;

	case 4:
	    for (k = 0 ; k < n ; k++)
	    {
		i = P [k] ;
		B [i      ] = X [4*k    ] ;
		B [i + d  ] = X [4*k + 1] ;
		B [i + d*2] = X [4*k + 2] ;
		B [i + d*3] = X [4*k + 3] ;
	    }
	    break ;

    }
}


/* ========================================================================== */
/* === inv_trans ============================================================ */
/* ========================================================================== */

/* b = x, transpose the row-major work vector x back into the column-major
 * right-hand-side.  Operates on 2, 3, or 4 columns at a time.
 */

static void inv_trans
(
    int n,
    double X [ ],		    /* n-by-nr in row-major order */
    double B [ ],		    /* n-by-nr in col-major order */
    int d,			    /* leading dimension of B */
    int nr			    /* nr = 2, 3, or 4 */
)
{
    int k ;

    switch (nr)
    {

	case 2:
	    for (k = 0 ; k < n ; k++)
	    {
		B [k      ] = X [2*k    ] ;
		B [k + d  ] = X [2*k + 1] ;
	    }
	    break ;

	case 3:
	    for (k = 0 ; k < n ; k++)
	    {
		B [k      ] = X [3*k    ] ;
		B [k + d  ] = X [3*k + 1] ;
		B [k + d*2] = X [3*k + 2] ;
	    }
	    break ;

	case 4:
	    for (k = 0 ; k < n ; k++)
	    {
		B [k      ] = X [4*k    ] ;
		B [k + d  ] = X [4*k + 1] ;
		B [k + d*2] = X [4*k + 2] ;
		B [k + d*3] = X [4*k + 3] ;
	    }
	    break ;

    }
}


/* ========================================================================== */
/* === perm ================================================================= */
/* ========================================================================== */

/* b = P*b, permute the right-hand-side, but keep in column order.
 */

static void perm
(
    int P [ ],
    int n,
    double X [ ],		    /* workspace of size n */
    double B [ ],		    /* n-by-nrhs in col-major order */
    int d,			    /* leading dimension of B */
    int nrhs,			    /* nrhs >= 1 */
    int ordering
)
{
    int j, k ;
    if (P == NULL || ordering == CHOLMOD_NATURAL)
    {
	/* P = NULL means that P is identity: nothing to do */
	return ;
    }
    for (j = 0 ; j < nrhs ; j++)
    {
	for (k = 0 ; k < n ; k++)
	{
	    X [k] = B [P [k] + d*j] ;
	}
	for (k = 0 ; k < n ; k++)
	{
	    B [k + d*j] = X [k] ;
	}
    }
}


/* ========================================================================== */
/* === inv_perm ============================================================= */
/* ========================================================================== */

/* b = P'*b, inverse permute the right-hand-side
 */

static void inv_perm
(
    int P [ ],
    int n,
    double X [ ],		    /* workspace of size n */
    double B [ ],		    /* n-by-nrhs in col-major order */
    int d,			    /* leading dimension of B */
    int nrhs,			    /* nrhs >= 1 */
    int ordering
)
{
    int j, k ;
    if (P == NULL || ordering == CHOLMOD_NATURAL)
    {
	/* P = NULL means that P is identity: nothing to do */
	return ;
    }
    for (j = 0 ; j < nrhs ; j++)
    {
	for (k = 0 ; k < n ; k++)
	{
	    X [P [k]] = B [k + d*j] ;
	}
	for (k = 0 ; k < n ; k++)
	{
	    B [k + d*j] = X [k] ;
	}
    }
}


/* ========================================================================== */
/* === ldl_solve_d ========================================================== */
/* ========================================================================== */

/* solve Dx=b, x in column form, for a simplicial LDL' factorization.
 */

static void ldl_solve_d
(
    cholmod_factor *L,
    double X [ ],		    /* n-by-nrhs in column form */
    int d,			    /* leading dimension of X */
    int nrhs			    /* nrhs >= 1 */
)
{
    double dk ;
    double *Lx ;
    int *Li, *Lp ;
    int n, j, k ;
    Lp = L->p ;
    Li = L->i ;
    Lx = L->x ;
    n = L->n ;
    ASSERT (L->ftype >= CHOLMOD_LDL_PACKED
	 && L->ftype <= CHOLMOD_LDL_DYNAMIC) ;

    if (nrhs == 1)
    {
	for (k = 0 ; k < n ; k++)
	{
	    X [k] /= Lx [Lp [k]] ;
	}
    }
    else
    {
	for (k = 0 ; k < n ; k++)
	{
	    dk = Lx [Lp [k]] ;
	    for (j = 0 ; j < nrhs ; j++)
	    {
		X [k + j*d] /= dk ;
	    }
	}
    }
}




/* ========================================================================== */
/* === simplicial solvers =================================================== */
/* ========================================================================== */

/* LL' packed, solve Lx=b with non-unit diagonal */
#define LL
#include "cholmod_lsolve.c"

/* LDL' unpacked or dynamic, solve LDx=b */
#define LD
#include "cholmod_lsolve.c"

/* LDL' packed, solve LDx=b */
#define LD
#define PACKED
#include "cholmod_lsolve.c"

/* LDL' packed, solve Lx=b with unit diagonal */
#define PACKED
#include "cholmod_lsolve.c"

/* LDL' unpacked or dynamic, solve Lx=b with unit diagonal */
#include "cholmod_lsolve.c"




/* LL' packed, solve L'x=b with non-unit diagonal */
#define LL
#include "cholmod_ltsolve.c"

/* LDL' unpacked or dynamic, solve DL'x=b */
#define LD
#include "cholmod_ltsolve.c"

/* LDL' packed, solve DL'x=b */
#define LD
#define PACKED
#include "cholmod_ltsolve.c"

/* LDL' packed, solve L'x=b with unit diagonal */
#define PACKED
#include "cholmod_ltsolve.c"

/* LDL' unpacked or dynamic, solve L'x=b with unit diagonal */
#include "cholmod_ltsolve.c"



/* ========================================================================== */
/* === cholmod_solve ======================================================== */
/* ========================================================================== */

/* Solve a linear system.
 *
 * The factorization can be simplicial LDL', simplicial LL', or supernodal LL'.
 * The Dx=b solve returns silently for the LL' factorizations (it is implicitly
 * identity).
 */

int cholmod_solve
(
    int sys,
    cholmod_factor *L,
    cholmod_dense *X,
    cholmod_common *Common
)
{
    double *W, *Xx ;
    int *P ;
    size_t wsize, maxrank ;
    int chunk, nr, n, use_Xwork, i, d, nrhs ;

    DEBUG (int orig) ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    RETURN_IF_NULL (X, FALSE) ;
    n = L->n ;
    d = X->d ;
    nrhs = X->ncol ;

    if (sys < CHOLMOD_A || sys > CHOLMOD_Pt)
    {
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_solve: invalid system", Common) ;
	return (FALSE) ;
    }

    if (d < n || X->nrow != L->n)
    {
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_solve: dimensions of L and X do not match", Common) ;
	return (FALSE) ;
    }

    DEBUG (cholmod_dump_factor (L, "L", Common)) ;
    DEBUG (cholmod_dump_dense  (X, "X", Common)) ;
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* quick return */
    /* ---------------------------------------------------------------------- */

    if (n == 0 || nrhs == 0)
    {
	/* nothing to do */
	return (TRUE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    wsize = 0 ;
    if (sys == CHOLMOD_D)
    {
	wsize = 0 ;
    }
    else if (sys == CHOLMOD_P || sys == CHOLMOD_Pt)
    {
	wsize = n ;
    }
    else if (L->ftype == CHOLMOD_LL_SUPER)
    {
	/* workspace for E in cholmod_super_lsolve and cholmod_super_ltsolve */
	wsize = nrhs * (L->maxesize) ;
	if (sys == CHOLMOD_A && L->ordering != CHOLMOD_NATURAL)
	{
	    /* workspace for W in perm and inv_perm */
	    wsize = MAX (L->n, wsize) ;
	}
    }
    else if (nrhs > 1 || (sys == CHOLMOD_A && L->ordering != CHOLMOD_NATURAL))
    {
	/* workspace for row-form X in simplicial LL' and LDL' solves, but
	 * not needed for an unpermuted single right-hand-side */
	wsize = MIN (4,nrhs) * n ;
    }
    wsize = MAX (1, wsize) ;

    /* make sure maxrank is in the proper range */
    maxrank = cholmod_maxrank (n, Common) ;

    use_Xwork = (wsize <= maxrank * n) ;
    if (use_Xwork)
    {
	/* allocate workspace in Common, and keep it there when done */
	cholmod_allocate_work (0, 0, wsize, sizeof (double), Common) ;
	W = Common->Xwork ;
	DEBUG (orig = Common->malloc_count) ;
    }
    else
    {
	/* allocate temporary space and free it when done */
	DEBUG (orig = Common->malloc_count) ;
	W = cholmod_malloc (wsize, sizeof (double), Common) ;
    }

    if (Common->status < CHOLMOD_OK || maxrank == 0)
    {
	/* out of memory for workspace */
	ASSERT (Common->malloc_count == orig) ;
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    DEBUG (cholmod_dump_factor (L, "L for solve", Common)) ;
    P = L->Perm ;
    Xx = X->x ;

    /* ---------------------------------------------------------------------- */
    /* solve using L, D, L', P, or some combination */
    /* ---------------------------------------------------------------------- */

    if (sys == CHOLMOD_P)
    {

	/* ------------------------------------------------------------------ */
	/* x = P*b */
	/* ------------------------------------------------------------------ */

	perm (P, n, W, Xx, d, nrhs, L->ordering) ;

    }
    else if (sys == CHOLMOD_Pt)
    {

	/* ------------------------------------------------------------------ */
	/* x = P'*b */
	/* ------------------------------------------------------------------ */

	inv_perm (P, n, W, Xx, d, nrhs, L->ordering) ;

    }
    else if (L->ftype == CHOLMOD_LDL_PACKED)
    {

	/* ------------------------------------------------------------------ */
	/* solve using a packed simplicial LDL' factorization */
	/* ------------------------------------------------------------------ */

	if ((sys == CHOLMOD_A && L->ordering == CHOLMOD_NATURAL) ||
	     sys == CHOLMOD_LDLt)
	{
	    for (chunk = 0 ; chunk < nrhs ; chunk += 4)
	    {
		nr = MIN (nrhs - chunk, 4) ;
		if (nr == 1)
		{
		    /* no need to copy Xx into W */
		    ldl_packed_lsolve_k (L, Xx, nr) ;
		    ldl_packed_dltsolve_k (L, Xx, nr) ;
		}
		else
		{
		    /* tranpose Xx into W, solve, then transpose back again */
		    trans (n, W, Xx, d, nr) ;
		    ldl_packed_lsolve_k (L, W, nr) ;
		    ldl_packed_dltsolve_k (L, W, nr) ;
		    inv_trans (n, W, Xx, d, nr) ;
		}
		Xx += 4*d ;
	    }
	}
	else if (sys == CHOLMOD_A)
	{
	    for (chunk = 0 ; chunk < nrhs ; chunk += 4)
	    {
		nr = MIN (nrhs - chunk, 4) ;
		ptrans (P, n, W, Xx, d, nr) ;
		ldl_packed_lsolve_k (L, W, nr) ;
		ldl_packed_dltsolve_k (L, W, nr) ;
		inv_ptrans (P, n, W, Xx, d, nr) ;
		Xx += 4*d ;
	    }
	}
	else if (sys == CHOLMOD_LD)
	{
	    for (chunk = 0 ; chunk < nrhs ; chunk += 4)
	    {
		nr = MIN (nrhs - chunk, 4) ;
		if (nr == 1)
		{
		    ldl_packed_ldsolve_k (L, Xx, nr) ;
		}
		else
		{
		    trans (n, W, Xx, d, nr) ;
		    ldl_packed_ldsolve_k (L, W, nr) ;
		    inv_trans (n, W, Xx, d, nr) ;
		}
		Xx += 4*d ;
	    }
	}
	else if (sys == CHOLMOD_L)
	{
	    for (chunk = 0 ; chunk < nrhs ; chunk += 4)
	    {
		nr = MIN (nrhs - chunk, 4) ;
		if (nr == 1)
		{
		    ldl_packed_lsolve_k (L, Xx, nr) ;
		}
		else
		{
		    trans (n, W, Xx, d, nr) ;
		    ldl_packed_lsolve_k (L, W, nr) ;
		    inv_trans (n, W, Xx, d, nr) ;
		}
		Xx += 4*d ;
	    }
	}
	else if (sys == CHOLMOD_Lt)
	{
	    for (chunk = 0 ; chunk < nrhs ; chunk += 4)
	    {
		nr = MIN (nrhs - chunk, 4) ;
		if (nr == 1)
		{
		    ldl_packed_ltsolve_k (L, Xx, nr) ;
		}
		else
		{
		    trans (n, W, Xx, d, nr) ;
		    ldl_packed_ltsolve_k (L, W, nr) ;
		    inv_trans (n, W, Xx, d, nr) ;
		}
		Xx += 4*d ;
	    }
	}
	else if (sys == CHOLMOD_DLt)
	{
	    for (chunk = 0 ; chunk < nrhs ; chunk += 4)
	    {
		nr = MIN (nrhs - chunk, 4) ;
		if (nr == 1)
		{
		    ldl_packed_dltsolve_k (L, Xx, nr) ;
		}
		else
		{
		    trans (n, W, Xx, d, nr) ;
		    ldl_packed_dltsolve_k (L, W, nr) ;
		    inv_trans (n, W, Xx, d, nr) ;
		}
		Xx += 4*d ;
	    }
	}
	else if (sys == CHOLMOD_D)
	{
	    ldl_solve_d (L, Xx, d, nrhs) ;
	}

    }
    else if (L->ftype == CHOLMOD_LDL_UNPACKED
	  || L->ftype == CHOLMOD_LDL_DYNAMIC)
    {

	/* ------------------------------------------------------------------ */
	/* solve using an unpacked simplicial LDL' factorization */
	/* ------------------------------------------------------------------ */

	if ((sys == CHOLMOD_A && L->ordering == CHOLMOD_NATURAL) ||
	     sys == CHOLMOD_LDLt)
	{
	    for (chunk = 0 ; chunk < nrhs ; chunk += 4)
	    {
		nr = MIN (nrhs - chunk, 4) ;
		if (nr == 1)
		{
		    ldl_unpacked_lsolve_k (L, Xx, nr) ;
		    ldl_unpacked_dltsolve_k (L, Xx, nr) ;
		}
		else
		{
		    trans (n, W, Xx, d, nr) ;
		    ldl_unpacked_lsolve_k (L, W, nr) ;
		    ldl_unpacked_dltsolve_k (L, W, nr) ;
		    inv_trans (n, W, Xx, d, nr) ;
		}
		Xx += 4*d ;
	    }
	}
	else if (sys == CHOLMOD_A)
	{
	    for (chunk = 0 ; chunk < nrhs ; chunk += 4)
	    {
		nr = MIN (nrhs - chunk, 4) ;
		ptrans (P, n, W, Xx, d, nr) ;
		ldl_unpacked_lsolve_k (L, W, nr) ;
		ldl_unpacked_dltsolve_k (L, W, nr) ;
		inv_ptrans (P, n, W, Xx, d, nr) ;
		Xx += 4*d ;
	    }
	}
	else if (sys == CHOLMOD_LD)
	{
	    for (chunk = 0 ; chunk < nrhs ; chunk += 4)
	    {
		nr = MIN (nrhs - chunk, 4) ;
		if (nr == 1)
		{
		    ldl_unpacked_ldsolve_k (L, Xx, nr) ;
		}
		else
		{
		    trans (n, W, Xx, d, nr) ;
		    ldl_unpacked_ldsolve_k (L, W, nr) ;
		    inv_trans (n, W, Xx, d, nr) ;
		}
		Xx += 4*d ;
	    }
	}
	else if (sys == CHOLMOD_L)
	{
	    for (chunk = 0 ; chunk < nrhs ; chunk += 4)
	    {
		nr = MIN (nrhs - chunk, 4) ;
		if (nr == 1)
		{
		    ldl_unpacked_lsolve_k (L, Xx, nr) ;
		}
		else
		{
		    trans (n, W, Xx, d, nr) ;
		    ldl_unpacked_lsolve_k (L, W, nr) ;
		    inv_trans (n, W, Xx, d, nr) ;
		}
		Xx += 4*d ;
	    }
	}
	else if (sys == CHOLMOD_Lt)
	{
	    for (chunk = 0 ; chunk < nrhs ; chunk += 4)
	    {
		nr = MIN (nrhs - chunk, 4) ;
		if (nr == 1)
		{
		    ldl_unpacked_ltsolve_k (L, Xx, nr) ;
		}
		else
		{
		    trans (n, W, Xx, d, nr) ;
		    ldl_unpacked_ltsolve_k (L, W, nr) ;
		    inv_trans (n, W, Xx, d, nr) ;
		}
		Xx += 4*d ;
	    }
	}
	else if (sys == CHOLMOD_DLt)
	{
	    for (chunk = 0 ; chunk < nrhs ; chunk += 4)
	    {
		nr = MIN (nrhs - chunk, 4) ;
		if (nr == 1)
		{
		    ldl_unpacked_dltsolve_k (L, Xx, nr) ;
		}
		else
		{
		    trans (n, W, Xx, d, nr) ;
		    ldl_unpacked_dltsolve_k (L, W, nr) ;
		    inv_trans (n, W, Xx, d, nr) ;
		}
		Xx += 4*d ;
	    }
	}
	else if (sys == CHOLMOD_D)
	{
	    ldl_solve_d (L, Xx, d, nrhs) ;
	}

    }
    else if (L->ftype == CHOLMOD_LL_PACKED)
    {

	/* ------------------------------------------------------------------ */
	/* solve using a simplicial LL' factorization */
	/* ------------------------------------------------------------------ */

	if ((sys == CHOLMOD_A && L->ordering == CHOLMOD_NATURAL) ||
	     sys == CHOLMOD_LDLt)
	{
	    for (chunk = 0 ; chunk < nrhs ; chunk += 4)
	    {
		nr = MIN (nrhs - chunk, 4) ;
		if (nr == 1)
		{
		    ll_lsolve_k  (L, Xx, nr) ;
		    ll_ltsolve_k (L, Xx, nr) ;
		}
		else
		{
		    trans (n, W, Xx, d, nr) ;
		    ll_lsolve_k  (L, W, nr) ;
		    ll_ltsolve_k (L, W, nr) ;
		    inv_trans (n, W, Xx, d, nr) ;
		}
		Xx += 4*d ;
	    }
	}
	else if (sys == CHOLMOD_A)
	{
	    for (chunk = 0 ; chunk < nrhs ; chunk += 4)
	    {
		nr = MIN (nrhs - chunk, 4) ;
		ptrans (P, n, W, Xx, d, nr) ;
		ll_lsolve_k  (L, W, nr) ;
		ll_ltsolve_k (L, W, nr) ;
		inv_ptrans (P, n, W, Xx, d, nr) ;
		Xx += 4*d ;
	    }
	}
	else if (sys == CHOLMOD_L || sys == CHOLMOD_LD)
	{
	    for (chunk = 0 ; chunk < nrhs ; chunk += 4)
	    {
		nr = MIN (nrhs - chunk, 4) ;
		if (nr == 1)
		{
		    ll_lsolve_k  (L, Xx, nr) ;
		}
		else
		{
		    trans (n, W, Xx, d, nr) ;
		    ll_lsolve_k  (L, W, nr) ;
		    inv_trans (n, W, Xx, d, nr) ;
		}
		Xx += 4*d ;
	    }
	}
	else if (sys == CHOLMOD_Lt || sys == CHOLMOD_DLt)
	{
	    for (chunk = 0 ; chunk < nrhs ; chunk += 4)
	    {
		nr = MIN (nrhs - chunk, 4) ;
		if (nr == 1)
		{
		    ll_ltsolve_k (L, Xx, nr) ;
		}
		else
		{
		    trans (n, W, Xx, d, nr) ;
		    ll_ltsolve_k (L, W, nr) ;
		    inv_trans (n, W, Xx, d, nr) ;
		}
		Xx += 4*d ;
	    }
	}

    }
    else if (L->ftype == CHOLMOD_LL_SUPER)
    {

	/* ------------------------------------------------------------------ */
	/* solve using a supernodal LL' factorization */
	/* ------------------------------------------------------------------ */

#ifndef NSUPERNODAL

	cholmod_dense *E, Ematrix ;
	E = &Ematrix ;
	E->nrow = nrhs ;
	E->ncol = L->maxesize ;
	E->nzmax = nrhs * (L->maxesize) ;
	E->d = nrhs ;
	E->x = W ;
	E->z = NULL ;
	E->xtype = CHOLMOD_REAL ;
	E->dtype = CHOLMOD_DOUBLE ;

	if ((sys == CHOLMOD_A && L->ordering == CHOLMOD_NATURAL) ||
	     sys == CHOLMOD_LDLt)
	{
	    cholmod_super_lsolve  (L, X, E, Common) ;
	    cholmod_super_ltsolve (L, X, E, Common) ;
	}
	else if (sys == CHOLMOD_A)
	{
	    perm (P, n, W, Xx, d, nrhs, L->ordering) ;
	    cholmod_super_lsolve  (L, X, E, Common) ;
	    cholmod_super_ltsolve (L, X, E, Common) ;
	    inv_perm (P, n, W, Xx, d, nrhs, L->ordering) ;
	}
	else if (sys == CHOLMOD_L || sys == CHOLMOD_LD)
	{
	    cholmod_super_lsolve (L, X, E, Common) ;
	}
	else if (sys == CHOLMOD_Lt || sys == CHOLMOD_DLt)
	{
	    cholmod_super_ltsolve (L, X, E, Common) ;
	}
#else
	/* CHOLMOD Supernodal module not installed */
	cholmod_error (CHOLMOD_NOT_INSTALLED,
		"cholmod_solve: Supernodal module not installed", Common) ;
#endif

    }

    /* ---------------------------------------------------------------------- */
    /* free or clear workspace */
    /* ---------------------------------------------------------------------- */

    if (use_Xwork)
    {
	for (i = 0 ; i < ((int) wsize) ; i++)
	{
	    W [i] = 0 ;
	}
    }
    else
    {
	cholmod_free (W, wsize, sizeof (double), Common) ;
    }

    DEBUG (cholmod_dump_dense (X, "X", Common)) ;
    ASSERT (Common->malloc_count == orig) ;
    return (Common->status == CHOLMOD_OK) ;
}
