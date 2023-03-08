/* ========================================================================= */
/* === TRILINOS_AMD_order =========================================================== */
/* ========================================================================= */

/* ------------------------------------------------------------------------- */
/* AMD, Copyright (c) Timothy A. Davis,					     */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README.txt for License.     */
/* email: davis at cise.ufl.edu    CISE Department, Univ. of Florida.        */
/* web: http://www.cise.ufl.edu/research/sparse/amd                          */
/* ------------------------------------------------------------------------- */

/* User-callable AMD minimum degree ordering routine.  See amd.h for
 * documentation.
 */

#include "trilinos_amd_internal.h"

/* ========================================================================= */
/* === TRILINOS_AMD_order =========================================================== */
/* ========================================================================= */

GLOBAL Int TRILINOS_AMD_order
(
    Int n,
    const Int Ap [ ],
    const Int Ai [ ],
    Int P [ ],
    double Control [ ],
    double Info [ ]
)
{
    Int *Len, *S, nz, i, *Pinv, info, status, *Rp, *Ri, *Cp, *Ci, ok ;
    size_t nzaat, slen ;
    double mem = 0 ;

#ifndef NDEBUG
    TRILINOS_AMD_debug_init ("amd") ;
#endif

    /* clear the Info array, if it exists */
    info = Info != (double *) NULL ;
    if (info)
    {
	for (i = 0 ; i < TRILINOS_AMD_INFO ; i++)
	{
	    Info [i] = TRILINOS_AMD_EMPTY ;
	}
	Info [TRILINOS_AMD_N] = n ;
	Info [TRILINOS_AMD_STATUS] = TRILINOS_AMD_OK ;
    }

    /* make sure inputs exist and n is >= 0 */
    if (Ai == (Int *) NULL || Ap == (Int *) NULL || P == (Int *) NULL || n < 0)
    {
	if (info) Info [TRILINOS_AMD_STATUS] = TRILINOS_AMD_INVALID ;
	return (TRILINOS_AMD_INVALID) ;	    /* arguments are invalid */
    }

    if (n == 0)
    {
	return (TRILINOS_AMD_OK) ;	    /* n is 0 so there's nothing to do */
    }

    nz = Ap [n] ;
    if (info)
    {
	Info [TRILINOS_AMD_NZ] = nz ;
    }
    if (nz < 0)
    {
	if (info) Info [TRILINOS_AMD_STATUS] = TRILINOS_AMD_INVALID ;
	return (TRILINOS_AMD_INVALID) ;
    }

    /* check if n or nz will cause size_t overflow */
    if (((size_t) n) >= SIZE_T_MAX / sizeof (Int)
     || ((size_t) nz) >= SIZE_T_MAX / sizeof (Int))
    {
	if (info) Info [TRILINOS_AMD_STATUS] = TRILINOS_AMD_OUT_OF_MEMORY ;
	return (TRILINOS_AMD_OUT_OF_MEMORY) ;	    /* problem too large */
    }

    /* check the input matrix:	TRILINOS_AMD_OK, TRILINOS_AMD_INVALID, or TRILINOS_AMD_OK_BUT_JUMBLED */
    status = TRILINOS_AMD_valid (n, n, Ap, Ai) ;

    if (status == TRILINOS_AMD_INVALID)
    {
	if (info) Info [TRILINOS_AMD_STATUS] = TRILINOS_AMD_INVALID ;
	return (TRILINOS_AMD_INVALID) ;	    /* matrix is invalid */
    }

    /* allocate two size-n integer workspaces */
    Len = (Int*) trilinos_amd_malloc (n * sizeof (Int)) ;
    Pinv = (Int*) trilinos_amd_malloc (n * sizeof (Int)) ;
    mem += n ;
    mem += n ;
    if (!Len || !Pinv)
    {
	/* :: out of memory :: */
	trilinos_amd_free (Len) ;
	trilinos_amd_free (Pinv) ;
	if (info) Info [TRILINOS_AMD_STATUS] = TRILINOS_AMD_OUT_OF_MEMORY ;
	return (TRILINOS_AMD_OUT_OF_MEMORY) ;
    }

    if (status == TRILINOS_AMD_OK_BUT_JUMBLED)
    {
	/* sort the input matrix and remove duplicate entries */
	TRILINOS_AMD_DEBUG1 (("Matrix is jumbled\n")) ;
	Rp = (Int*) trilinos_amd_malloc ((n+1) * sizeof (Int)) ;
	Ri = (Int*) trilinos_amd_malloc (MAX (nz,1) * sizeof (Int)) ;
	mem += (n+1) ;
	mem += MAX (nz,1) ;
	if (!Rp || !Ri)
	{
	    /* :: out of memory :: */
	    trilinos_amd_free (Rp) ;
	    trilinos_amd_free (Ri) ;
	    trilinos_amd_free (Len) ;
	    trilinos_amd_free (Pinv) ;
	    if (info) Info [TRILINOS_AMD_STATUS] = TRILINOS_AMD_OUT_OF_MEMORY ;
	    return (TRILINOS_AMD_OUT_OF_MEMORY) ;
	}
	/* use Len and Pinv as workspace to create R = A' */
	TRILINOS_AMD_preprocess (n, Ap, Ai, Rp, Ri, Len, Pinv) ;
	Cp = Rp ;
	Ci = Ri ;
    }
    else
    {
	/* order the input matrix as-is.  No need to compute R = A' first */
	Rp = NULL ;
	Ri = NULL ;
	Cp = (Int *) Ap ;
	Ci = (Int *) Ai ;
    }

    /* --------------------------------------------------------------------- */
    /* determine the symmetry and count off-diagonal nonzeros in A+A' */
    /* --------------------------------------------------------------------- */

    nzaat = TRILINOS_AMD_aat (n, Cp, Ci, Len, P, Info) ;
    TRILINOS_AMD_DEBUG1 (("nzaat: %g\n", (double) nzaat)) ;
    ASSERT ((MAX (nz-n, 0) <= nzaat) && (nzaat <= 2 * (size_t) nz)) ;

    /* --------------------------------------------------------------------- */
    /* allocate workspace for matrix, elbow room, and 6 size-n vectors */
    /* --------------------------------------------------------------------- */

    S = NULL ;
    slen = nzaat ;			/* space for matrix */
    ok = ((slen + nzaat/5) >= slen) ;	/* check for size_t overflow */
    slen += nzaat/5 ;			/* add elbow room */
    for (i = 0 ; ok && i < 7 ; i++)
    {
	ok = ((slen + n) > slen) ;	/* check for size_t overflow */
	slen += n ;			/* size-n elbow room, 6 size-n work */
    }
    mem += slen ;
    ok = ok && (slen < SIZE_T_MAX / sizeof (Int)) ; /* check for overflow */
    ok = ok && (slen < Int_MAX) ;	/* S[i] for Int i must be OK */
    if (ok)
    {
	S = (Int*) trilinos_amd_malloc (slen * sizeof (Int)) ;
    }
    TRILINOS_AMD_DEBUG1 (("slen %g\n", (double) slen)) ;
    if (!S)
    {
	/* :: out of memory :: (or problem too large) */
	trilinos_amd_free (Rp) ;
	trilinos_amd_free (Ri) ;
	trilinos_amd_free (Len) ;
	trilinos_amd_free (Pinv) ;
	if (info) Info [TRILINOS_AMD_STATUS] = TRILINOS_AMD_OUT_OF_MEMORY ;
	return (TRILINOS_AMD_OUT_OF_MEMORY) ;
    }
    if (info)
    {
	/* memory usage, in bytes. */
	Info [TRILINOS_AMD_MEMORY] = mem * sizeof (Int) ;
    }

    /* --------------------------------------------------------------------- */
    /* order the matrix */
    /* --------------------------------------------------------------------- */

    TRILINOS_AMD_1 (n, Cp, Ci, P, Pinv, Len, slen, S, Control, Info) ;

    /* --------------------------------------------------------------------- */
    /* free the workspace */
    /* --------------------------------------------------------------------- */

    trilinos_amd_free (Rp) ;
    trilinos_amd_free (Ri) ;
    trilinos_amd_free (Len) ;
    trilinos_amd_free (Pinv) ;
    trilinos_amd_free (S) ;
    if (info) Info [TRILINOS_AMD_STATUS] = status ;
    return (status) ;	    /* successful ordering */
}
