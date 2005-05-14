/* ========================================================================== */
/* === Core/cholmod_factor ================================================== */
/* ========================================================================== */

/*
 * CHOLMOD/Core version 0.1. May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* Core utility routines for the cholmod_factor object:
 *
 * The data structure for an LL' or LDL' factorization is too complex to
 * describe in one sentence.  This object can hold the symbolic analysis alone,
 * or in combination with a "simplicial" (similar to a sparse matrix) or
 * "supernodal" form of the numerical factorization.  Only the routine to free
 * a factor is primary, since a factor object is created by the factorization
 * routine (cholmod_factorize).  It must be freed with cholmod_free_factor.
 *
 * Primary routine:
 * ----------------
 * cholmod_free_factor		free a factor
 *
 * Secondary routines:
 * -------------------
 * cholmod_allocate_factor	allocate a symbolic factor (LL' or LDL')
 * cholmod_reallocate_factor	change the # entries in a factor 
 * cholmod_change_ftype		change the type of factor (e.g., LDL' to LL')
 * cholmod_pack_factor		pack the columns of a factor
 * cholmod_reallocate_column	resize a single column of a factor
 * cholmod_factor_to_sparse	create a sparse matrix copy of a factor
 * cholmod_copy_factor		create a copy of a factor
 *
 * Note that there is no cholmod_sparse_to_factor routine to create a factor
 * as a copy of a sparse matrix.  It could be done, after a fashion, but a
 * lower triangular sparse matrix would not necessarily have a chordal graph,
 * which would break the many CHOLMOD routines that rely on this property.
 *
 * The cholmod_factor_to_sparse routine is provided so that matrix operations
 * in the MatrixOps module may be applied to L.  Those operations operate on
 * cholmod_sparse objects, and they are not guaranteed to maintain the chordal
 * property of L.  Such a modified L cannot be safely convert back to a
 * cholmod_factor object.
 */

#include "cholmod_core.h"
#include "cholmod_internal.h"

/* ========================================================================== */
/* === cholmod_allocate_factor ============================================== */
/* ========================================================================== */

/* Allocate a CHOLMOD_SYMBOLIC factor, with L->Perm and L->ColCount allocated
 * and initialized to "empty" values (Perm [k] = k, and ColCount[k] = 1).
 * The integer and numerical parts of L are not allocated.
 *
 * This is sufficient (but far from ideal) for input to cholmod_factorize,
 * since the dynamic LDL' factorization can reallocate the columns of L as
 * needed.  The primary purpose of this routine is to allocate space for a
 * symbolic factorization, for the "expert" user to do his or her own
 * symbolic analysis.  The typical user should use cholmod_analyze instead of
 * this routine.
 *
 * workspace: none
 */

cholmod_factor *cholmod_allocate_factor
(
    size_t n,
    cholmod_common *Common
)
{
    int j ;
    int *Perm, *ColCount ;
    cholmod_factor *L ;
    DEBUG (int orig) ;

    RETURN_IF_NULL_COMMON (FALSE) ;
    Common->status = CHOLMOD_OK ;

    DEBUG (orig = Common->malloc_count) ;

    L = cholmod_malloc (1, sizeof (cholmod_factor), Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	ASSERT (Common->malloc_count == orig) ;
	return (NULL) ;	    /* out of memory */
    }
    L->n = n ;
    L->ftype = CHOLMOD_SYMBOLIC ;
    L->itype = Common->itype ;
    L->xtype = CHOLMOD_PATTERN ;
    L->dtype = Common->dtype ;

    /* allocate the purely symbolic part of L */
    L->ordering = CHOLMOD_NATURAL ;
    Perm = cholmod_malloc (n, sizeof (int), Common) ;
    ColCount = cholmod_malloc (n, sizeof (int), Common) ;
    L->Perm = Perm ;
    L->ColCount = ColCount ;

    /* simplicial part of L is empty */
    L->nzmax = 0 ;
    L->p = NULL ;
    L->i = NULL ;
    L->x = NULL ;
    L->z = NULL ;
    L->nz = NULL ;
    L->next = NULL ;
    L->prev = NULL ;

    /* supernodal part of L is also empty */
    L->nsuper = 0 ;
    L->ssize = 0 ;
    L->xsize = 0 ;
    L->maxesize = 0 ;
    L->maxcsize = 0 ;
    L->super = NULL ;
    L->pi = NULL ;
    L->px = NULL ;
    L->s = NULL ;

    /* L has not been factorized */
    L->minor = n ;

    if (Common->status < CHOLMOD_OK)
    {
	cholmod_free_factor (&L, Common) ;
	ASSERT (Common->malloc_count == orig) ;
	return (NULL) ;		/* out of memory */
    }

    /* initialize Perm and ColCount */
    for (j = 0 ; j < ((int) n) ; j++)
    {
	Perm [j] = j ;
    }
    for (j = 0 ; j < ((int) n) ; j++)
    {
	ColCount [j] = 1 ;
    }

    ASSERT (Common->malloc_count == orig + 3) ;
    return (L) ;
}


/* ========================================================================== */
/* === cholmod_free_factor ================================================== */
/* ========================================================================== */

/* Free a factor object.
 *
 * workspace: none
 */

int cholmod_free_factor
(
    cholmod_factor **LHandle,
    cholmod_common *Common
)
{
    int n, lnz, xs, ss, s ;
    cholmod_factor *L ;

    RETURN_IF_NULL_COMMON (FALSE) ;

    if (LHandle == NULL)
    {
	/* nothing to do */
	return (TRUE) ;
    }
    L = *LHandle ;
    if (L == NULL)
    {
	/* nothing to do */
	return (TRUE) ;
    }

    n = L->n ;
    lnz = L->nzmax ;
    s = L->nsuper + 1 ;
    xs = (L->ftype == CHOLMOD_LL_SUPER) ? ((int) (L->xsize)) : (lnz) ;
    ss = L->ssize ;

    /* symbolic part of L */
    L->Perm = cholmod_free (L->Perm,     n, sizeof (int), Common);
    L->ColCount = cholmod_free (L->ColCount, n, sizeof (int), Common);

    /* simplicial form of L */
    L->p    = cholmod_free (L->p,    n+1, sizeof (int), Common) ;
    L->i    = cholmod_free (L->i,    lnz, sizeof (int), Common) ;
    L->x    = cholmod_free (L->x,    xs,  sizeof (double), Common) ;
    L->nz   = cholmod_free (L->nz,   n,   sizeof (int), Common) ;
    L->next = cholmod_free (L->next, n+2, sizeof (int), Common) ;
    L->prev = cholmod_free (L->prev, n+2, sizeof (int), Common) ;

    /* supernodal form of L */
    L->pi    = cholmod_free (L->pi,   s,   sizeof (int),    Common) ;
    L->px    = cholmod_free (L->px,   s,   sizeof (int),    Common) ;
    L->super = cholmod_free (L->super,s,   sizeof (int),    Common) ;
    L->s     = cholmod_free (L->s,    ss,  sizeof (int),    Common) ;

    *LHandle = cholmod_free ((*LHandle), 1, sizeof (cholmod_factor), Common) ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_reallocate_factor ============================================ */
/* ========================================================================== */

/* Change the size of L->i and L->x, or allocate them if their current size
 * is zero.  L must be simplicial.
 *
 * workspace: none
 */

int cholmod_reallocate_factor
(
    cholmod_factor *L,
    size_t nznew,
    cholmod_common *Common
)
{
    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    if (L->ftype < CHOLMOD_SYMBOLIC || L->ftype > CHOLMOD_LL_PACKED)
    {
	/* L must be simplicial, and not symbolic */
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_reallocate_factor: L invalid", Common) ;
	return (FALSE) ;
    }
    Common->status = CHOLMOD_OK ;
    PRINT1 (("realloc factor %d to %d\n", L->nzmax, nznew)) ;

    /* ---------------------------------------------------------------------- */
    /* resize the factor */
    /* ---------------------------------------------------------------------- */

    cholmod_realloc_multiple (1, 1, &(L->i), NULL, &(L->x), NULL, &(L->nzmax),
	    nznew, Common) ;
    return (Common->status == CHOLMOD_OK) ;
}


/* ========================================================================== */
/* === cholmod_reallocate_column =========================================== */
/* ========================================================================== */

/* Column j needs more space, reallocate it at the end of L->i and L->x.
 * The factor is converted to a dynamic LDL' if it is not one already.
 * If the reallocation fails, the factor is converted to a simplicial
 * symbolic ftype (no pattern, just L->Perm and L->ColCount).
 *
 * workspace: none
 */

int cholmod_reallocate_column	    /* returns TRUE if OK, FALSE otherwise */
(
    cholmod_factor *L,
    size_t j,		/* the column to reallocate (in range 0 to L->n-1)*/
    size_t need,	/* required size of column j (must be >= 0) */
    cholmod_common *Common
)
{
    double xneed ;
    double *Lx ;
    int n, pold, pnew, len, k, tail ;
    int *Lp, *Lprev, *Lnext, *Li, *Lnz ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    n = L->n ;
    if (j >= L->n || need == 0)
    {
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_reallocate_column: j invalid", Common) ;
	return (FALSE) ;	    /* j out of range */
    }
    if (L->ftype != CHOLMOD_LDL_DYNAMIC)
    {
	if (!cholmod_change_ftype (L, CHOLMOD_LDL_DYNAMIC, Common))
	{
	    /* out of memory, convert to simplicial symbolic */
	    (void) cholmod_change_ftype (L, CHOLMOD_SYMBOLIC, Common) ;
	    cholmod_error (CHOLMOD_OUT_OF_MEMORY,
		"out of memory, L converted to symbolic", Common) ;
	    return (FALSE) ;	    /* out of memory */
	}
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* increase the size of L if needed */
    /* ---------------------------------------------------------------------- */

    /* head = n+1 ; */
    tail = n ;
    Lp = L->p ;
    Lnz = L->nz ;
    Lprev = L->prev ;
    Lnext = L->next ;

    ASSERT (Lnz != NULL) ;
    ASSERT (Lnext != NULL && Lprev != NULL) ;
    PRINT1 (("col %d need %d\n", j, need)) ;

    /* column j cannot have more than n-j entries if all entries are present */
    need = MIN (need, n-j) ;

    /* compute need in double to avoid int overflow */
    if (Common->grow1 >= 1.0)
    {
	xneed = (double) need ;
	xneed = Common->grow1 * xneed + Common->grow2 ;
	xneed = MIN (xneed, n-j) ;
	need = (int) xneed ;
    }
    PRINT1 (("new need %d\n", need)) ;
    ASSERT (need >= 1 && need <= n-j) ;

    DEBUG (cholmod_dump_factor (L, "start colrealloc", Common)) ;

    if (Lp [tail] + need > L->nzmax)
    {
	/* use double to avoid int overflow */
	xneed = (double) need ;
	if (Common->grow0 < 1.2)	    /* fl. pt. compare, false if NaN */
	{
	    /* if grow0 is less than 1.2 or NaN, don't use it */
	    xneed = 1.2 * (((double) L->nzmax) + xneed + 1) ;
	}
	else
	{
	    xneed = Common->grow0 * (((double) L->nzmax) + xneed + 1) ;
	}
	if (xneed > INT_MAX ||
		!cholmod_reallocate_factor (L, (int) xneed, Common))
	{
	    /* out of memory, convert to simplicial symbolic */
	    cholmod_change_ftype (L, CHOLMOD_SYMBOLIC, Common) ;
	    cholmod_error (CHOLMOD_OUT_OF_MEMORY,
		"out of memory, L converted to symbolic", Common) ;
	    return (FALSE) ;	    /* out of memory */
	}
	PRINT1 (("\n=== GROW L from %d to %d\n", L->nzmax, (int) xneed)) ;
	/* pack all columns so that each column has at most grow2 free space */
	cholmod_pack_factor (L, Common) ;
	ASSERT (Common->status == CHOLMOD_OK) ;
    }

    /* ---------------------------------------------------------------------- */
    /* reallocate the column */
    /* ---------------------------------------------------------------------- */

    Li = L->i ;
    Lx = L->x ;

    /* remove j from its current position in the list */
    Lnext [Lprev [j]] = Lnext [j] ;
    Lprev [Lnext [j]] = Lprev [j] ;

    /* place j at the end of the list */
    Lnext [Lprev [tail]] = j ;
    Lprev [j] = Lprev [tail] ;
    Lnext [j] = n ;
    Lprev [tail] = j ;

    /* allocate space for column j */
    pold = Lp [j] ;
    pnew = Lp [tail] ;
    Lp [j] = pnew  ;
    Lp [tail] += need ;

    /* copy column j to the new space */
    len = Lnz [j] ;
    for (k = 0 ; k < len ; k++)
    {
	Li [pnew + k] = Li [pold + k] ;
	Lx [pnew + k] = Lx [pold + k] ;
    }

    DEBUG (cholmod_dump_factor (L, "colrealloc done", Common)) ;

    /* successful reallocation of column j of L */
    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_change_ftype ================================================= */
/* ========================================================================== */

/* Change the ftype of a cholmod_factor object.
 *
 * There are four basic classes of factor ftypes:
 *
 * (1) simplicial symbolic:  Consists of two size-n arrays: the fill-reducing
 *	permutation (L->Perm) and the nonzero count for each column of L
 *	(L->ColCount).  All other factor ftypes also include this information.
 *	L->ColCount may be exact (obtained from the analysis routines), or
 *	it may be a guess.  During factorization, and certainly after update/
 *	downdate, the columns of L can have a different number of nonzeros.
 *	L->ColCount is used to allocate space.  L->ColCount is exact for the
 *	supernodal factorizations.  The nonzero pattern of L is not kept.
 *
 * (2) simplicial numeric:  These represent L in a compressed column form.  The
 *	variants of this ftype are:
 *
 *	packed LL':	L has a non-unit diagonal.  Row indices in column j are
 *	    located in L->i [L->p [j] ... L->p [j+1]-1], and corresponding
 *	    numeric values are in the same locations in L->x.  The total number
 *	    of entries in L is L->p [n].
 *
 *	packed LDL':	L is unit diagonal.  Row indices in column j are
 *	    located in L->i [L->p [j] ... L->p [j+1]-1], and corresponding
 *	    numeric values are in the same locations in L->x.  The total number
 *	    of entries in L is L->p [n].  The unit diagonal is not stored; D is
 *	    stored on the diagonal of L instead.
 *
 *	unpacked LDL':	L is unit diagonal.  Row indices in column j are
 *	    located in L->i [L->p [j] ... L->p [j] + L->nz [j]], and
 *	    corresponding numeric values are in the same locations in L->x.
 *	    The total number of entries is the sum of L->nz [j].
 *	    The unit diagonal is not stored; D is stored on the diagonal of L
 *	    instead.  L->p is monotonic (that is, L->p [j] <= L->p [j+1],
 *	    except that L->p [n] is not defined).
 *
 *	dynamic LDL':	L is unit diagonal.  Row indices in column j are
 *	    located in L->i [L->p [j] ... L->p [j] + L->nz [j]], and
 *	    corresponding numeric values are in the same locations in L->x.
 *	    The total number of entries is the sum of L->nz [j].
 *	    The unit diagonal is not stored; D is stored on the diagonal of L
 *	    instead.  L->p may or may not be monotonic.  The order of
 *	    storage of the columns in L->i and L->x is given by a doubly-linked
 *	    list (L->prev and L->next).
 *
 * (3) supernodal symbolic:  A representation of the nonzero pattern of the
 *	supernodes for a supernodal factorization.  There are L->nsuper
 *	supernodes.  Columns L->super [k] to L->super [k+1]-1 are in the kth
 *	supernode.  The row indices for the kth supernode are in
 *	L->s [L->pi [k] ... L->pi [k+1]-1].  The numerical values are not
 *	allocated (L->x), but when they are they will be located in
 *	L->x [L->px [k] ... L->px [k+1]-1], and the L->px array is defined
 *	in this factor ftype.
 *
 * (4) supernodal numeric:  L is non-unit diagonal.  L->x contains the numerical
 *	values of the supernodes, as described above.
 *
 * The cholmod_change_ftype routine can do almost all possible conversions:
 *
 *	Overall, the simplicial ftypes are very flexible.  Converting to them
 *	and among them, any conversions can be made and when they are converted
 *	the contents of the factor are always well defined.  The supernodal
 *	ftypes are less flexible.
 *
 *	(1) Convert any factor to a simplicial symbolic ftype.  All information
 *	    other than L->Perm and L->ColCount is deallocated.
 *
 *	(2) Convert between any simplicial numeric ftype, or convert a
 *	    supernodal numeric ftype to any simplicial numeric ftype.  The
 *	    numerical values and nonzero pattern of L and D are preserved.
 *
 *	(3) Convert any symbolic ftype to any simplicial numeric ftype.
 *	    Space for the pattern and values of L (and D for an LDL' ftype)
 *	    is allocated based on L->ColCount, and the contents of L (and D
 *	    if applicable) are set to the identity matrix.  Supernodal
 *	    information is not preserved if converting from a supernodal
 *	    symbolic ftype.
 *
 *	(4) Convert from a symplicial symbolic to a supernodal symbolic.  Space
 *	    is allocated but not initialized (this conversion is used by the
 *	    supernodal analysis routine and unlikely to be needed by the user).
 *
 *	(5) Convert from supernodal symbolic to a supernodal numeric.  Space
 *	    is allocated but not initialized (this conversion is used by the
 *	    supernodal factorization routine and unlikely to be needed by the
 *	    user).
 *
 *	(6) Convert from a supernodal numeric to supernodal symbolic ftype.
 *	    The numerical values are deallocated, but the supernodal pattern
 *	    is kept.  This is useful if you wish to free some space now, but do
 *	    a supernodal factorization of another matrix later on with the same
 *	    nonzero pattern as the one you already factorized, using the same
 *	    fill-reducing ordering.
 *
 *	(7) Converting a factor to its current ftype does nothing, except for
 *	    the packed LL' and LDL' ftypes, in which case the space is reduced
 *	    to exactly what is needed (L->i and L->x become size L->p [n]).
 *
 * The list of what this routine cannot do is shorter:
 *
 *	(1) Simplicial numeric ftypes cannot be converted to a supernodal
 *	    symbolic ftypes.  This would simultaneously deallocate the
 *	    simplicial pattern and numeric values and reallocate uninitialized
 *	    space for the supernodal pattern.  This isn't useful for the user,
 *	    and not needed by CHOLMOD's own routines either.
 *
 *	(2) The only ftype that can be converted to a supernodal numeric factor
 *	    is a supernodal symbolic ftype.
 *
 * workspace: no conversion routine uses workspace in Common.  No temporary
 *	workspace is allocated.
 */

/* ========================================================================== */
/* === is_monotonic ========================================================= */
/* ========================================================================== */

/* Determine if columns of L are monotonic or not.  */

static int is_monotonic
(
    cholmod_factor *L
)
{
    int j, n ;
    int *Lp, *Lnz ;
    ASSERT (L->ftype == CHOLMOD_LDL_UNPACKED
	 || L->ftype == CHOLMOD_LDL_DYNAMIC) ;
    n = L->n ;
    Lp = L->p ;
    Lnz = L->nz ;
    for (j = 0 ; j < n-1 ; j++)
    {
	ASSERT (Lp [j] != Lp [j+1] && Lnz [j] > 0) ;
	if (Lp [j+1] < Lp [j])
	{
	    /* column j+1 comes before column j, thus not monotonic */
	    return (FALSE) ;
	}
    }
    return (TRUE) ;	/* all columns are in their natural order */
}


/* ========================================================================== */
/* === natural_list ========================================================= */
/* ========================================================================== */

/* Create a naturally-ordered doubly-linked list of columns. */

static void natural_list (cholmod_factor *L)
{
    int head, tail, n, j ;
    int *Lnext, *Lprev ;
    PRINT1 (("create natural list, factor ftype is %d\n", L->ftype)) ;
    Lnext = L->next ;
    Lprev = L->prev ;
    ASSERT (Lprev != NULL && Lnext != NULL) ;
    n = L->n ;
    head = n+1 ;
    tail = n ;
    Lnext [head] = 0 ;
    Lprev [head] = EMPTY ;
    Lnext [tail] = EMPTY ;
    Lprev [tail] = n-1 ;
    for (j = 0 ; j < n ; j++)
    {
	Lnext [j] = j+1 ;
	Lprev [j] = j-1 ;
    }
    Lprev [0] = head ;
}


/* ========================================================================== */
/* === allocate_simplicial_numeric ========================================== */
/* ========================================================================== */

/* Allocate O(n) arrays for simplicial numeric factorization.  Initializes
 * the link lists only.  Does not change L->ftype. */

static int allocate_simplicial_numeric
(
    cholmod_factor *L,
    int ftype,
    cholmod_common *Common
)
{
    int packed, dynamic, n ;
    int *Lp, *Lnz, *Lprev, *Lnext ;

    PRINT1 (("Allocate simplicial\n")) ;

    ASSERT (L->ftype == CHOLMOD_SYMBOLIC ||
	    L->ftype == CHOLMOD_SYMBOLIC_SUPER ||
	    L->ftype == CHOLMOD_LL_SUPER) ;
    ASSERT (ftype >= CHOLMOD_LDL_PACKED && ftype <= CHOLMOD_LL_PACKED) ;
    ASSERT (L->p == NULL) ;
    ASSERT (L->nz == NULL) ;
    ASSERT (L->prev == NULL) ;
    ASSERT (L->next == NULL) ;

    n = L->n ;

    /* all L matrices require an L->p array */
    Lp = cholmod_malloc (n+1, sizeof (int), Common) ;

    packed  = (ftype == CHOLMOD_LDL_PACKED || ftype == CHOLMOD_LL_PACKED) ;
    dynamic = (ftype == CHOLMOD_LDL_DYNAMIC) ;

    if (packed)
    {
	/* a packed matrix has no L->nz */
	Lnz = NULL ;
    }
    else
    {
	/* for the unpacked or dynamic cases */
	Lnz = cholmod_malloc (n, sizeof (int), Common) ;
    }

    if (dynamic)
    {
	/* for the dynamic case */
	Lprev = cholmod_malloc (n+2, sizeof (int), Common) ;
	Lnext = cholmod_malloc (n+2, sizeof (int), Common) ;
    }
    else
    {
	/* for the packed or unpacked cases (not dynamic) */
	Lnext = NULL ;
	Lprev = NULL ;
    }

    if (Common->status < CHOLMOD_OK)
    {
	cholmod_free (Lp   , n+1, sizeof (int), Common) ;
	cholmod_free (Lnz  , n,   sizeof (int), Common) ;
	cholmod_free (Lprev, n+2, sizeof (int), Common) ;
	cholmod_free (Lnext, n+2, sizeof (int), Common) ;
	PRINT1 (("Allocate simplicial failed, %d\n", L->ftype)) ;
	return (FALSE) ;	/* out of memory */
    }

    /* ---------------------------------------------- commit the changes to L */

    L->p = Lp ;
    L->nz = Lnz ;
    L->prev = Lprev ;
    L->next = Lnext ;
    if (dynamic)
    {
	/* initialize a doubly linked list for columns in natural order */
	natural_list (L) ;
    }
    PRINT1 (("Allocate simplicial done %d\n", L->ftype)) ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === simplicial_symbolic_to_super_symbolic ================================ */
/* ========================================================================== */

/* Convert a simplicial symbolic factor supernodal symbolic factor.  Does not
 * initialize the new space. */

static int simplicial_symbolic_to_super_symbolic
(
    cholmod_factor *L,
    cholmod_common *Common
)
{
    int nsuper, xsize, ssize ;
    int *Lsuper, *Lpi, *Lpx, *Ls ;

    ASSERT (L->ftype == CHOLMOD_SYMBOLIC) ;

    xsize  = L->xsize ;
    ssize  = L->ssize ;
    nsuper = L->nsuper ;

    PRINT1 (("simple sym to super sym: ssize %d xsize %d nsuper %d status %d\n",
	    ssize, xsize, nsuper, Common->status)) ;

    /* O(nsuper) arrays, where nsuper <= n */
    Lsuper = cholmod_malloc (nsuper+1, sizeof (int), Common) ;
    Lpi    = cholmod_malloc (nsuper+1, sizeof (int), Common) ;
    Lpx    = cholmod_malloc (nsuper+1, sizeof (int), Common) ;

    /* O(ssize) array, where ssize <= nnz(L), and usually much smaller */
    Ls = cholmod_malloc (ssize, sizeof (int), Common) ;

    if (Common->status < CHOLMOD_OK)
    {
	cholmod_free (Lsuper, nsuper+1, sizeof (int), Common) ;
	cholmod_free (Lpi,    nsuper+1, sizeof (int), Common) ;
	cholmod_free (Lpx,    nsuper+1, sizeof (int), Common) ;
	cholmod_free (Ls,     ssize,    sizeof (int), Common) ;
	return (FALSE) ;	/* out of memory */
    }

    /* ---------------------------------------------- commit the changes to L */

    ASSERT (Lsuper != NULL && Lpi != NULL && Lpx != NULL && Ls != NULL) ;

    L->maxcsize = 0 ;
    L->maxesize = 0 ;

    L->super = Lsuper ;
    L->pi = Lpi ;
    L->px = Lpx ;
    L->s  = Ls ;
    Ls [0] = EMPTY ;	    /* supernodal pattern undefined */

    L->ftype = CHOLMOD_SYMBOLIC_SUPER ;
    L->xtype = CHOLMOD_PATTERN ;
    L->dtype = Common->dtype ;
    L->minor = L->n ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === any_to_simplicial_symbolic =========================================== */
/* ========================================================================== */

/* Convert any factor L to a simplicial symbolic factor. */

static void any_to_simplicial_symbolic
(
    cholmod_factor *L,
    cholmod_common *Common
)
{
    int n, lnz, xs, ss, s ;

    n = L->n ;
    lnz = L->nzmax ;
    s = L->nsuper + 1 ;
    xs = (L->ftype == CHOLMOD_LL_SUPER) ? ((int) (L->xsize)) : (lnz) ;
    ss = L->ssize ;

    /* ---------------------------------------------- commit the changes to L */

    /* free all but the symbolic analysis; L->ftype can be anything */
    L->p     = cholmod_free (L->p    , n+1, sizeof (int),    Common) ;
    L->i     = cholmod_free (L->i    , lnz, sizeof (int),    Common) ;
    L->x     = cholmod_free (L->x    , xs,  sizeof (double), Common) ;
    L->nz    = cholmod_free (L->nz   , n,   sizeof (int),    Common) ;
    L->next  = cholmod_free (L->next , n+2, sizeof (int),    Common) ;
    L->prev  = cholmod_free (L->prev , n+2, sizeof (int),    Common) ;
    L->super = cholmod_free (L->super, s,   sizeof (int),    Common) ;
    L->pi    = cholmod_free (L->pi   , s,   sizeof (int),    Common) ;
    L->px    = cholmod_free (L->px   , s,   sizeof (int),    Common) ;
    L->s     = cholmod_free (L->s    , ss,  sizeof (int),    Common) ;
    L->nzmax = 0 ;
    L->ftype = CHOLMOD_SYMBOLIC ;
    L->xtype = CHOLMOD_PATTERN ;
    L->dtype = Common->dtype ;
    L->minor = n ;
    DEBUG (cholmod_dump_factor (L, "done  to simplicial symbolic", Common)) ;
}


/* ========================================================================== */
/* === ll_super_to_super_symbolic =========================================== */
/* ========================================================================== */

/* Convert a numerical supernodal L to symbolic supernodal.  Cannot fail. */

static void ll_super_to_super_symbolic
(
    cholmod_factor *L,
    cholmod_common *Common
)
{

    /* ---------------------------------------------- commit the changes to L */

    /* free all but the supernodal numerical factor */
    ASSERT (L->ftype == CHOLMOD_LL_SUPER) ;
    DEBUG (cholmod_dump_factor (L, "start to super symbolic", Common)) ;
    L->x = cholmod_free (L->x, L->xsize, sizeof (double), Common) ;
    L->ftype = CHOLMOD_SYMBOLIC_SUPER ;
    L->xtype = CHOLMOD_PATTERN ;
    L->dtype = Common->dtype ;
    L->minor = L->n ;
    DEBUG (cholmod_dump_factor (L, "done  to super symbolic", Common)) ;
}


/* ========================================================================== */
/* === simplicial_symbolic_to_simplicial_numeric ============================ */
/* ========================================================================== */

/* Convert a simplicial symbolic L to a simplicial numeric L; allocate space
 * for L using L->ColCount from symbolic analysis, and set L to identity. */

static void simplicial_symbolic_to_simplicial_numeric
(
    cholmod_factor *L,
    int ftype,
    cholmod_common *Common
)
{
    double grow0, grow1, xlen, xlnz ;
    double *Lx ;
    int *Li, *Lp, *Lnz, *ColCount ;
    int grow2, n, grow, p, j, lnz, packed, len, ok ;

    n = L->n ;

    ASSERT (L->ftype == CHOLMOD_SYMBOLIC) ;
    ASSERT (ftype >= CHOLMOD_LDL_PACKED && ftype <= CHOLMOD_LL_PACKED) ;
    packed = (ftype == CHOLMOD_LDL_PACKED || ftype == CHOLMOD_LL_PACKED) ;

    if (!allocate_simplicial_numeric (L, ftype, Common))
    {
	PRINT1 (("out of memory, allocate simplicial numeric %d\n", L->ftype)) ;
	return ;	/* out of memory */
    }

    ASSERT (L->ColCount != NULL) ;
    ASSERT (IMPLIES (!packed, L->nz != NULL)) ;

    ColCount = L->ColCount ;
    Lnz = L->nz ;
    Lp = L->p ;
    ok = TRUE ;

    if (packed)
    {

	/* ------------------------------------------------------------------ */
	/* LDL' or LL' packed */
	/* ------------------------------------------------------------------ */

	PRINT1 (("convert to packed LL' or LDL'\n")) ;
	lnz = 0 ;
	for (j = 0 ; ok && j < n ; j++)
	{
	    /* ensure len is in the range 1 to n-j */
	    len = ColCount [j] ;
	    len = MAX (1, len) ;
	    len = MIN (len, n-j) ;
	    lnz += len ;
	    ok = (lnz >= 0) ;
	}
	for (j = 0 ; j <= n ; j++)
	{
	    Lp [j] = j ;
	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* LDL' unpacked or LDL' dynamic */
	/* ------------------------------------------------------------------ */

	PRINT1 (("convert to unpacked or dynamic\n")) ;
	/* compute new lnzmax */
	/* if any parameter is NaN, grow is false */
	grow0 = Common->grow0 ;
	grow1 = Common->grow1 ;
	grow2 = Common->grow2 ;
	/* fl.pt. comparisons below are false if any parameter is NaN */
	grow = (grow0 >= 1.0) && (grow1 >= 1.0) && (grow2 > 0) ;
	PRINT1 (("init, grow1 %g grow2 %d\n", grow1, grow2)) ;
	/* initialize Lp and Lnz for each column */
	lnz = 0 ;
	for (j = 0 ; ok && j < n ; j++)
	{
	    Lp [j] = lnz ;
	    Lnz [j] = 1 ;

	    /* ensure len is in the range 1 to n-j */
	    len = ColCount [j] ;
	    len = MAX (1, len) ;
	    len = MIN (len, n-j) ;

	    /* compute len in double to avoid int overflow */
	    PRINT1 (("ColCount [%d] = %d\n", j, len)) ;
	    if (grow)
	    {
		xlen = (double) len ;
		xlen = grow1 * xlen + grow2 ;
		xlen = MIN (xlen, n-j) ;
		len = (int) xlen ;
	    }
	    ASSERT (len >= 1 && len <= n-j) ;
	    lnz += len ;
	    ok = (lnz >= 0) ;
	}
	if (ok)
	{
	    Lp [n] = lnz ;
	    if (grow)
	    {
		/* add extra space */
		xlnz = (double) lnz ;
		xlnz *= grow0 ;
		xlnz = MIN (xlnz, INT_MAX) ;
		xlnz = MIN (xlnz, ((double) n * (double) n + (double) n) / 2) ;
		lnz = (int) xlnz ;
	    }
	}
    }

    if (!ok)
    {
	cholmod_error (CHOLMOD_TOO_LARGE, "problem too large", Common) ;
    }

    /* allocate L->i and L->x */
    PRINT1 (("resizing to lnz %d\n", lnz)) ;
    if (!ok || !cholmod_reallocate_factor (L, lnz, Common))
    {
	L->p    = cholmod_free (L->p   , n+1, sizeof (int),    Common) ;
	L->nz   = cholmod_free (L->nz  , n,   sizeof (int),    Common) ;
	L->prev = cholmod_free (L->prev, n+2, sizeof (int),    Common) ;
	L->next = cholmod_free (L->next, n+2, sizeof (int),    Common) ;
	L->i    = cholmod_free (L->i   , lnz, sizeof (int),    Common) ;
	L->x    = cholmod_free (L->x   , lnz, sizeof (double), Common) ;
	PRINT1 (("cannot realloc simplicial numeric %d\n", L->ftype)) ;
	return ;	/* out of memory */
    }

    /* ---------------------------------------------- commit the changes to L */

    /* initialize L to be the identity matrix */
    L->ftype = ftype ;
    L->xtype = Common->xtype ;
    L->dtype = Common->dtype ;
    Li = L->i ;
    Lx = L->x ;

    /* create the unit diagonal for either the LL' or LDL' case */
    for (j = 0 ; j < n ; j++)
    {
	ASSERT (Lp [j] < Lp [j+1]) ;
	p = Lp [j] ;
	Li [p] = j ;
	Lx [p] = 1.0 ;
    }

    PRINT1 (("done convert simplicial symbolic to numeric %d\n", L->ftype)) ;
}


/* ========================================================================== */
/* === pack_ldl_unpacked ==================================================== */
/* ========================================================================== */

/* Pack the columns of an LDL' unpacked factor, in place, so that each column
 * has at most grow2 free space.  Does not reallocate anything.  Cannot fail.
 */

static void pack_ldl_unpacked
(
    cholmod_factor *L,
    int grow2
)
{
    double *Lx ;
    int pnew, j, k, pold, len, n ;
    int *Li, *Lp, *Lnz ;

    grow2 = MAX (0, grow2) ;
    PRINT1 (("\nPACK monotonic: grow2 %d\n", grow2)) ;
    ASSERT (L->ftype == CHOLMOD_LDL_UNPACKED ||
	    L->ftype == CHOLMOD_LDL_DYNAMIC) ;
    ASSERT (is_monotonic (L)) ;

    pnew = 0 ;
    n = L->n ;
    Lp = L->p ;
    Li = L->i ;
    Lnz = L->nz ;
    Lx = L->x ;
    Lp [n] = L->nzmax ;

    for (j = 0 ; j < n ; j++)
    {
	/* pack column j */
	pold = Lp [j] ;
	len = Lnz [j] ;
	ASSERT (len > 0) ;
	PRINT2 (("col %d pnew %d pold %d\n", j, pnew, pold)) ;
	if (pnew < pold)
	{
	    PRINT2 (("    pack this column\n")) ;
	    for (k = 0 ; k < len ; k++)
	    {
		Li [pnew + k] = Li [pold + k] ;
		Lx [pnew + k] = Lx [pold + k] ;
	    }
	    Lp [j] = pnew ;
	}
	len = MIN (len + grow2, n - j) ;
	pnew = MIN (Lp [j] + len, Lp [j+1]) ;
    }
    Lp [n] = pnew ;
    PRINT2 (("Lp [n] = %d\n", pnew)) ;
}


/* ========================================================================== */
/* === ldl_dynamic_to_ldl =================================================== */
/* ========================================================================== */

/* Convert a dynamic LDL' to packed or unpacked LDL'.  Always succeeds if
 * the columns are monotonic. */

static void ldl_dynamic_to_ldl
(
    cholmod_factor *L,
    int ftype,
    cholmod_common *Common
)
{
    double grow0, grow1, xlen, xlnz ;
    double *newLx, *Lx ;
    int monotonic, grow2, n, j, lnz, len, grow, pnew, pold, k, packed, ok ;
    int *newLi, *Lp, *Li, *Lnz ;

    PRINT1 (("\n===Convert LDL dyn to %d\n", ftype)) ;
    DEBUG (cholmod_dump_factor (L, "start LDL dyn to any LDL", Common)) ;
    ASSERT (L->ftype == CHOLMOD_LDL_DYNAMIC) ;
    ASSERT (ftype == CHOLMOD_LDL_UNPACKED || ftype == CHOLMOD_LDL_PACKED) ;

    monotonic = is_monotonic (L) ;
    packed = (ftype == CHOLMOD_LDL_PACKED) ;

    n = L->n ;
    Lp = L->p ;
    Li = L->i ;
    Lx = L->x ;
    Lnz = L->nz ;
    grow = FALSE ;
    grow0 = Common->grow0 ;
    grow1 = Common->grow1 ;
    grow2 = Common->grow2 ;
    ok = TRUE ;

    if (!monotonic)
    {
	/* Columns out of order.  Need to copy L into new space */
	PRINT1 (("L dynamic is non-monotonic\n")) ;

	/* compute new L->nzmax */
	if (!packed)
	{
	    /* if any parameter is NaN, grow is false */
	    /* fl.pt. comparisons below are false if any parameter is NaN */
	    grow = (grow0 >= 1.0) && (grow1 >= 1.0) && (grow2 >= 0) ;
	}
	lnz = 0 ;
	for (j = 0 ; ok && j < n ; j++)
	{
	    len = Lnz [j] ;
	    ASSERT (len >= 1 && len <= n-j) ;

	    /* compute len in double to avoid int overflow */
	    if (grow)
	    {
		xlen = (double) len ;
		xlen = grow1 * xlen + grow2 ;
		xlen = MIN (xlen, n-j) ;
		len = (int) xlen ;
	    }
	    ASSERT (len >= Lnz [j] && len <= n-j) ;

	    PRINT2 (("j: %d Lnz[j] %d len %d p %d\n", j, Lnz [j], len, lnz)) ;

	    lnz += len ;
	    ok = (lnz >= 0) ;
	}

	if (!ok)
	{
	    cholmod_error (CHOLMOD_TOO_LARGE, "problem too large", Common) ;
	    return ;
	}

	if (grow)
	{
	    xlnz = (double) lnz ;
	    xlnz *= grow0 ;
	    xlnz = MIN (xlnz, INT_MAX) ;
	    xlnz = MIN (xlnz, ((double) n * (double) n + (double) n) / 2) ;
	    lnz = (int) xlnz ;
	}

	PRINT1 (("final lnz %d\n", lnz)) ;

	/* allocate new space */
	newLi = cholmod_malloc (lnz, sizeof (int), Common) ;
	newLx = cholmod_malloc (lnz, sizeof (double), Common) ;
	if (Common->status < CHOLMOD_OK)
	{
	    cholmod_free (newLi, lnz, sizeof (int),    Common) ;
	    cholmod_free (newLx, lnz, sizeof (double), Common) ;
	    return ;	    /* out of memory */
	}
    }

    /* ---------------------------------------------- commit the changes to L */

    if (!monotonic)
    {
	pnew = 0 ;
	for (j = 0 ; j < n ; j++)
	{
	    /* copy and pack column j */
	    len = Lnz [j] ;
	    PRINT2 (("j: %d Lnz[j] %d len %d p %d\n", j, Lnz [j], len, pnew)) ;
	    pold = Lp [j] ;
	    for (k = 0 ; k < len ; k++)
	    {
		newLi [pnew + k] = Li [pold + k] ;
		newLx [pnew + k] = Lx [pold + k] ;
	    }
	    Lp [j] = pnew ;

	    /* compute len in double to avoid int overflow */
	    if (grow)
	    {
		xlen = (double) len ;
		xlen = grow1 * xlen + grow2 ;
		xlen = MIN (xlen, n-j) ;
		len = (int) xlen ;
	    }
	    ASSERT (len >= Lnz [j] && len <= n-j) ;
	    pnew += len ;
	    ASSERT (pnew > 0) ;	    /* int overflow case already covered */
	}
	Lp [n] = pnew ;
	PRINT1 (("final pnew = %d, lnz %d lnzmax %d\n", pnew, lnz, L->nzmax)) ;
	ASSERT (pnew <= lnz) ;

	/* free the old L->i and L->x and replace with the new ones */
	L->i = cholmod_free (L->i, L->nzmax, sizeof (int), Common) ;
	L->x = cholmod_free (L->x, L->nzmax, sizeof (double), Common) ;
	L->i = newLi ;
	L->x = newLx ;
	L->nzmax = lnz ;
    }
    else if (packed)
    {
	/* convert to packed LDL' */
	pack_ldl_unpacked (L, 0) ;
    }

    /* free the link list */
    L->next = cholmod_free (L->next, n+2, sizeof (int), Common) ;
    L->prev = cholmod_free (L->prev, n+2, sizeof (int), Common) ;

    if (ftype == CHOLMOD_LDL_PACKED)
    {
	/* free Lnz if L is being converted to LDL' packed */
	L->nz = cholmod_free (L->nz, n, sizeof (int), Common) ;
    }

    L->ftype = ftype ;

    DEBUG (cholmod_dump_factor (L, "done  LDL dyn to any LDL", Common)) ;
}


/* ========================================================================== */
/* === ldl_packed_to_ldl_unpacked =========================================== */
/* ========================================================================== */

/* Convert a packed LDL' to an unpacked LDL' by allocating and
 * computing L->nz */

static void ldl_packed_to_ldl_unpacked
(
    cholmod_factor *L,
    cholmod_common *Common
)
{
    int j, n ;
    int *Lnz, *Lp ;

    DEBUG (cholmod_dump_factor (L, "start LDL pack to LDL unpack", Common)) ;
    ASSERT (L->ftype == CHOLMOD_LDL_PACKED) ;
    ASSERT (L->nz == NULL) ;
    n = L->n ;
    Lnz = cholmod_malloc (n, sizeof (int), Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return ;	/* out of memory */
    }

    /* ---------------------------------------------- commit the changes to L */

    L->nz = Lnz ;
    Lp = L->p ;
    for (j = 0 ; j < n ; j++)
    {
	Lnz [j] = Lp [j+1] - Lp [j] ;
	ASSERT (Lnz [j] > 0) ;
    }
    ASSERT (Lp [n] <= (int) (L->nzmax)) ;
    L->ftype = CHOLMOD_LDL_UNPACKED ;
    DEBUG (cholmod_dump_factor (L, "done  LDL pack to LDL unpack", Common)) ;
}


/* ========================================================================== */
/* === pack_ldl_dynamic ===================================================== */
/* ========================================================================== */

/* Pack the columns of an LDL' dynamic factor.  Does not change the storage
 * order of the columns of L.  Cannot fail. */

static void pack_ldl_dynamic
(
    cholmod_factor *L,
    int grow2
)
{
    double *Lx ;
    int pnew, j, k, pold, len, n, head, tail ;
    int *Lp, *Li, *Lnz, *Lnext ;

    grow2 = MAX (0, grow2) ;
    PRINT1 (("\nPACK grow2 %d\n", grow2)) ;
    ASSERT (L->ftype == CHOLMOD_LDL_DYNAMIC) ;

    pnew = 0 ;
    n = L->n ;
    Lp = L->p ;
    Li = L->i ;
    Lx = L->x ;
    Lnz = L->nz ;
    Lnext = L->next ;

    head = n+1 ;
    tail = n ;

    for (j = Lnext [head] ; j != tail ; j = Lnext [j])
    {
	/* pack column j */
	pold = Lp [j] ;
	len = Lnz [j] ;
	ASSERT (len > 0) ;
	PRINT2 (("col %d pnew %d pold %d\n", j, pnew, pold)) ;
	if (pnew < pold)
	{
	    PRINT2 (("    pack this column\n")) ;
	    for (k = 0 ; k < len ; k++)
	    {
		Li [pnew + k] = Li [pold + k] ;
		Lx [pnew + k] = Lx [pold + k] ;
	    }
	    Lp [j] = pnew ;
	}
	len = MIN (len + grow2, n - j) ;
	pnew = MIN (Lp [j] + len, Lp [Lnext [j]]) ;
    }
    PRINT2 (("final pnew = %d\n", pnew)) ;
}


/* ========================================================================== */
/* === ldl_unpacked_to_ldl_packed =========================================== */
/* ========================================================================== */

/* Convert an unpacked LDL' to a packed LDL'.   Cannot fail.  */

static void ldl_unpacked_to_ldl_packed
(
    cholmod_factor *L,
    cholmod_common *Common
)
{
    int j, n, lnz ;
    int *Lnz ;

    DEBUG (cholmod_dump_factor (L, "start LDL unpack to LDL pack", Common)) ;
    ASSERT (L->ftype == CHOLMOD_LDL_UNPACKED) ;

    Lnz = L->nz ;
    n = L->n ;
    ASSERT (Lnz != NULL) ;


    lnz = 0 ;
    for (j = 0 ; j < n ; j++)
    {
	PRINT2 (("Lnz [%d] = %d\n", j, Lnz [j])) ;
	lnz += Lnz [j] ;
	ASSERT (Lnz [j] > 0) ;
	ASSERT (lnz > 0) ;  /* no int overflow since L is already allocated */
    }

    PRINT1 (("lnz %d L->nzmax %d\n", lnz, L->nzmax)) ;
    ASSERT (is_monotonic (L)) ;

    /* pack the columns of L so that each column has no free space at all. */

    pack_ldl_unpacked (L, 0) ;
    ASSERT (is_monotonic (L)) ;
    PRINT1 (("lnz %d L->nzmax %d n %d\n", lnz, L->nzmax, n)) ;
    ASSERT (lnz <= (int) (L->nzmax)) ;

    L->nz = cholmod_free (L->nz, n, sizeof (int), Common) ;

    /* reduce the size of L to just what is needed */
    cholmod_reallocate_factor (L, lnz, Common) ;
    Common->status = CHOLMOD_OK ;

    L->ftype = CHOLMOD_LDL_PACKED ;
    ASSERT ((int) (L->nzmax) == lnz) ;
    DEBUG (cholmod_dump_factor (L, "done  LDL unpack to LDL pack", Common)) ;
}


/* ========================================================================== */
/* === ldl_packed_to_ldl_dynamic ============================================ */
/* ========================================================================== */

/* Convert a packed LDL' to a dynamic LDL' by allocating and
 * computing L->nz, and creating the link lists */

static void ldl_packed_to_ldl_dynamic
(
    cholmod_factor *L,
    cholmod_common *Common
)
{
    int j, n ;
    int *Lnz, *Lp, *Lnext, *Lprev ;

    DEBUG (cholmod_dump_factor (L, "start LDL pack to LDL dyn", Common)) ;
    ASSERT (L->ftype == CHOLMOD_LDL_PACKED) ;
    ASSERT (L->nz == NULL && L->next == NULL && L->prev == NULL) ;
    n = L->n ;
    Lnz = cholmod_malloc (n, sizeof (int), Common) ;
    Lnext = cholmod_malloc (n+2, sizeof (int), Common) ;
    Lprev = cholmod_malloc (n+2, sizeof (int), Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	cholmod_free (Lnz, n, sizeof (int), Common) ;
	cholmod_free (Lnext, n+2, sizeof (int), Common) ;
	cholmod_free (Lprev, n+2, sizeof (int), Common) ;
	return ;	/* out of memory */
    }

    /* ---------------------------------------------- commit the changes to L */

    L->nz = Lnz ;
    Lp = L->p ;
    for (j = 0 ; j < n ; j++)
    {
	Lnz [j] = Lp [j+1] - Lp [j] ;
	ASSERT (Lnz [j] > 0) ;
    }

    L->next = Lnext ;
    L->prev = Lprev ;

    /* initialize a doubly linked list for columns in natural order */
    natural_list (L) ;

    L->ftype = CHOLMOD_LDL_DYNAMIC ;
    DEBUG (cholmod_dump_factor (L, "done  LDL pack to LDL dyn", Common)) ;
}


/* ========================================================================== */
/* === ldl_unpacked_to_ldl_dynamic ========================================== */
/* ========================================================================== */

/* Convert an unpacked LDL' to a dynamic LDL', by constructing the link list
 * L->prev and L->next.  */

static void ldl_unpacked_to_ldl_dynamic
(
    cholmod_factor *L,
    cholmod_common *Common
)
{
    int n ;
    int *Lnext, *Lprev ;

    DEBUG (cholmod_dump_factor (L, "start LDL unpack to LDL dyn", Common)) ;
    ASSERT (L->ftype == CHOLMOD_LDL_UNPACKED) ;
    ASSERT (is_monotonic (L)) ;
    n = L->n ;

    ASSERT (L->next == NULL && L->prev == NULL) ;
    Lnext = cholmod_malloc (n+2, sizeof (int), Common) ;
    Lprev = cholmod_malloc (n+2, sizeof (int), Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	cholmod_free (Lnext, n+2, sizeof (int), Common) ;
	cholmod_free (Lprev, n+2, sizeof (int), Common) ;
	return ;	/* out of memory */
    }

    /* ---------------------------------------------- commit the changes to L */

    L->next = Lnext ;
    L->prev = Lprev ;
    L->ftype = CHOLMOD_LDL_DYNAMIC ;

    /* initialize a doubly linked list for columns in natural order */
    natural_list (L) ;

    DEBUG (cholmod_dump_factor (L, "done  LDL unpack to LDL dyn", Common)) ;
}


/* ========================================================================== */
/* === ldl_packed_to_ll ===================================================== */
/* ========================================================================== */

/* Convert a packed LDL' to an LL' factorization.   Any column with a negative
 * or zero diagonal entry is not modified so that conversion back to LDL' will
 * succeed.  Cannot fail, but can result in a matrix L with a negative entry
 * on the diagonal.  If the kth entry on the diagonal of D is negative, the
 * it and the kth column of L are left unchanged.  A subsequent conversion back
 * to an LDL' form will also leave the column unchanged, so the correct LDL'
 * factorization will be restored.  L->minor is set to the smallest k for which
 * D(k,k) is negative.
 */

static void ldl_packed_to_ll
(
    cholmod_factor *L,
    cholmod_common *Common
)
{
    double dj, ljj ;
    double *Lx ;
    int j, n, p, pend ;
    int *Lp ;

    /* ---------------------------------------------- commit the changes to L */

    DEBUG (cholmod_dump_factor (L, "start LDL unpack to LL", Common)) ;
    ASSERT (L->ftype == CHOLMOD_LDL_PACKED) ;
    n = L->n ;
    Lp = L->p ;
    Lx = L->x ;

    L->minor = n ;
    for (j = 0 ; j < n ; j++)
    {
	p = Lp [j] ;
	pend = Lp [j+1] ;
	dj = Lx [p] ;
	if (dj <= 0)
	{
	    /* Conversion has failed; matrix is not positive definite.  Do not
	     * modify the column so that the LDL' factorization can be restored
	     * if desired. */
	    cholmod_error (CHOLMOD_NOT_POSDEF, "L not positive definite",
		    Common) ;
	    L->minor = MIN (L->minor, (size_t) j) ;
	}
	else
	{
	    ljj = sqrt (dj) ;
	    Lx [p++] = ljj ;
	    for ( ; p < pend ; p++)
	    {
		Lx [p] *= ljj ;
	    }
	}
    }

    L->ftype = CHOLMOD_LL_PACKED ;
    DEBUG (if ((int) (L->minor) != n) PRINT1 (("L not posdef %d\n", L->minor)));
    DEBUG (cholmod_dump_factor (L, "done LDL unpack to LL", Common)) ;
}


/* ========================================================================== */
/* === ll_to_ldl ============================================================ */
/* ========================================================================== */

/* Convert an LL' factorization to any LDL'.  Conversion to LDL' packed
 * cannot fail. */

static void ll_to_ldl
(
    cholmod_factor *L,
    int ftype,
    cholmod_common *Common
)
{
    int n, j, p, pend, packed, dynamic ;
    int *Lp, *Lnz, *Lnext, *Lprev ;
    double ljj ;
    double *Lx ;

    DEBUG (cholmod_dump_factor (L, "start LL to LDL", Common)) ;
    ASSERT (L->ftype == CHOLMOD_LL_PACKED) ;
    ASSERT (ftype == CHOLMOD_LDL_PACKED || ftype == CHOLMOD_LDL_UNPACKED ||
	    ftype == CHOLMOD_LDL_DYNAMIC) ;
    ASSERT (L->next == NULL && L->prev == NULL && L->nz == NULL) ;

    packed = (ftype == CHOLMOD_LDL_PACKED) ;
    dynamic = (ftype == CHOLMOD_LDL_DYNAMIC) ;

    n = L->n ;
    Lp = L->p ;
    Lx = L->x ;
    Lnz = NULL ;
    Lprev = NULL ;
    Lnext = NULL ;

    if (!packed)
    {
	Lnz = cholmod_malloc (n, sizeof (int), Common) ;
    }
    if (dynamic)
    {
	Lnext = cholmod_malloc (n+2, sizeof (int), Common) ;
	Lprev = cholmod_malloc (n+2, sizeof (int), Common) ;
    }

    if (Common->status < CHOLMOD_OK)
    {
	cholmod_free (Lnz  , n,   sizeof (int), Common) ;
	cholmod_free (Lnext, n+2, sizeof (int), Common) ;
	cholmod_free (Lprev, n+2, sizeof (int), Common) ;
	return ;	/* out of memory */
    }

    /* ---------------------------------------------- commit the changes to L */

    L->nz = Lnz ;

    for (j = 0 ; j < n ; j++)
    {
	p = Lp [j] ;
	pend = Lp [j+1] ;
	if (!packed)
	{
	    Lnz [j] = pend - p ;
	}
	ljj = Lx [p] ;
	if (ljj > 0)
	{
	    Lx [p++] = ljj*ljj ;
	    for ( ; p < pend ; p++)
	    {
		Lx [p] /= ljj ;
	    }
	}
    }

    L->next = Lnext ;
    L->prev = Lprev ;
    L->ftype = ftype ;

    if (dynamic)
    {
	/* initialize a doubly linked list for columns in natural order */
	natural_list (L) ;
    }

    DEBUG (cholmod_dump_factor (L, "done  LL to LDL", Common)) ;
}


/* ========================================================================== */
/* === ll_super_to_simplicial_numeric ======================================= */
/* ========================================================================== */

/* Convert a supernodal numeric factorization to any simplicial numeric one */

static void ll_super_to_simplicial_numeric
(
    cholmod_factor *L,
    int ftype,
    cholmod_common *Common
)
{
    double ljj ;
    double *Lx ;
    int n, lnz, packed, ldl, s, nsuper, p, psi, psx, psend, nsrow, nscol, ii,
	jj, j, k1, k2, erows, ll ;
    int *Ls, *Lpi, *Lpx, *Super, *Lp, *Li, *Lnz ;

    DEBUG (cholmod_dump_factor (L, "start LL super to simplicial", Common)) ;
    PRINT1 (("super -> simplicial (%d)\n", ftype)) ;
    ASSERT (L->ftype == CHOLMOD_LL_SUPER) ;
    ASSERT (L->x != NULL && L->i == NULL) ;

    n = L->n ;
    ldl = (ftype <= CHOLMOD_LDL_DYNAMIC) ;
    packed = (ftype == CHOLMOD_LDL_PACKED || ftype == CHOLMOD_LL_PACKED) ;
    ll = (ftype == CHOLMOD_LL_PACKED) ;

    nsuper = L->nsuper ;
    Lpi = L->pi ;
    Lpx = L->px ;
    Ls = L->s ;
    Super = L->super ;

    /* int overflow cannot occur since supernodal L already exists */

    if (packed)
    {
	/* count the number of nonzeros in L.  Each supernode is of the form
	 *
	 *    d		    For this example, nscol = 4 (# columns).
	 *    l d	    nsrow = 9.  The "l" entries are in an LDL' factor,
	 *    l l d	    with the "d" appearing in D.  Both "l" and "d" are
	 *    l l l d	    in the LL' factor.  For both LDL' and LL', the "e"
	 *    e e e e	    entries in the erows-by-nscol are placed in the
	 *    e e e e	    simplicial L.  Note that some "l" and "e" entries
	 *    e e e e	    may be numerically zero and even symbolically zero
	 *    e e e e	    if a tight simplicial factorization or resymbol
	 *    e e e e	    were done.
	 */
	lnz = 0 ;
	for (s = 0 ; s < nsuper ; s++)
	{
	    k1 = Super [s] ;
	    k2 = Super [s+1] ;
	    psi = Lpi [s] ;
	    psend = Lpi [s+1] ;
	    nsrow = psend - psi ;
	    nscol = k2 - k1 ;
	    ASSERT (nsrow >= nscol) ;
	    erows = nsrow - nscol ;

	    /* lower triangular part, including the diagonal,
	     * counting the "l" and "d" terms in the figure above. */
	    lnz += nscol * (nscol+1) / 2 ;

	    /* rectangular part, below the diagonal block (the "e" terms) */
	    lnz += nscol * erows ;
	}
	ASSERT (lnz <= (int) (L->xsize)) ;
    }
    else
    {
	/* Li will be the same size as Lx */
	lnz = L->xsize ;
    }
    ASSERT (lnz >= 0) ;
    PRINT1 (("simplicial lnz = %d  packed: %d  ldl: %d L->xsize %d\n",
		lnz, ldl, packed, L->xsize)) ;

    Li = cholmod_malloc (lnz, sizeof (int), Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return ;	/* out of memory */
    }

    if (!allocate_simplicial_numeric (L, ftype, Common))
    {
	cholmod_free (Li, lnz, sizeof (int), Common) ;
	return ;	/* out of memory */
    }

    /* ---------------------------------------------- commit the changes to L */

    L->i = Li ;
    L->ftype = ftype ;
    L->nzmax = lnz ;

    Lp = L->p ;
    Li = L->i ;
    Lx = L->x ;
    Lnz = L->nz ;

    p = 0 ;

    for (s = 0 ; s < nsuper ; s++)
    {
	k1 = Super [s] ;
	k2 = Super [s+1] ;
	psi = Lpi [s] ;
	psend = Lpi [s+1] ;
	psx = Lpx [s] ;
	nsrow = psend - psi ;
	nscol = k2 - k1 ;

	for (jj = 0 ; jj < nscol ; jj++)
	{
	    /* column j of L starts here */
	    j = jj + k1 ;

	    if (ll)
	    {

		/* convert to LL' packed */
		Lp [j] = p ;
		PRINT2 (("Col j %d p %d\n", Lp [j], p)) ;
		for (ii = jj ; ii < nsrow ; ii++)
		{
		    /* get L(i,j) from supernode and store in column j */
		    ASSERT (p < (int) (L->xsize) && p <= psx + ii + jj*nsrow) ;
		    Li [p] = Ls [psi + ii] ;
		    Lx [p] = Lx [psx + ii + jj*nsrow] ;
		    PRINT2 (("  i %d %g\n", Li [p], Lx [p])) ;
		    p++ ;
		}

	    }
	    else if (packed)
	    {

		/* convert to LDL' packed */
		Lp [j] = p ;
		PRINT2 (("Col j %d p %d\n", Lp [j], p)) ;
		ljj = Lx [psx + jj + jj*nsrow] ;

		if (ljj <= 0)
		{
		    /* the matrix is not positive definite; do not divide */
		    Lx [p] = ljj ;
		    ljj = 1 ;
		}
		else
		{
		    Lx [p] = ljj*ljj ;
		}

		Li [p] = j ;
		p++ ;
		for (ii = jj + 1 ; ii < nsrow ; ii++)
		{
		    /* get L(i,j) from supernode and store in column j */
		    ASSERT (p < (int) (L->xsize) && p <= psx + ii + jj*nsrow) ;
		    Li [p] = Ls [psi + ii] ;
		    Lx [p] = Lx [psx + ii + jj*nsrow] / ljj ;
		    PRINT2 (("  i %d %g\n", Li [p], Lx [p])) ;
		    p++ ;
		}

	    }
	    else
	    {

		/* convert to LDL' unpacked or dynamic */
		ASSERT (ldl) ;
		p = psx + jj + jj*nsrow ;
		Lp [j] = p ;
		ljj = Lx [p] ;

		if (ljj <= 0)
		{
		    /* the matrix is not positive definite; do not divide */
		    Lx [p] = ljj ;
		    ljj = 1 ;
		}
		else
		{
		    Lx [p] = ljj*ljj ;
		}

		Li [p] = j ;

		Lnz [j] = nsrow - jj ;
		p++ ;
		for (ii = jj + 1 ; ii < nsrow ; ii++)
		{
		    /* get L(i,j) from supernode and store in column j */
		    Li [psx + ii + jj*nsrow] = Ls [psi + ii] ;
		    Lx [psx + ii + jj*nsrow] /= ljj ;
		}
	    }
	}
    }

    if (packed)
    {
	Lp [n] = p ;
	PRINT1 (("Final Lp %d n %d lnz %d\n", p, n, lnz)) ;
	ASSERT (Lp [n] == lnz) ;
	ASSERT (lnz <= (int) (L->xsize)) ;
	/* reduce size of L->x to match L->i.  This cannot fail. */
	L->x = cholmod_realloc (L->x, &(L->xsize), lnz, sizeof (double),
		Common) ;
	ASSERT (lnz == (int) (L->xsize)) ;
	Common->status = CHOLMOD_OK ;
    }
    else
    {
	Lp [n] = Lpx [nsuper] ;
	ASSERT (Lp [n] == (int) (L->xsize)) ;
	ASSERT (Lp [n] == (int) (L->nzmax)) ;
    }

    /* free unused parts of L */
    L->super = cholmod_free (L->super, nsuper+1, sizeof (int), Common) ;
    L->pi    = cholmod_free (L->pi   , nsuper+1, sizeof (int), Common) ;
    L->px    = cholmod_free (L->px   , nsuper+1, sizeof (int), Common) ;
    L->s     = cholmod_free (L->s    , L->ssize, sizeof (int), Common) ;

    L->ssize = 0 ;
    L->xsize = 0 ;
    L->nsuper = 0 ;
    L->maxesize = 0 ;
    L->maxcsize = 0 ;

    DEBUG (cholmod_dump_factor (L, "done  LL super to simplicial", Common)) ;
}


/* ========================================================================== */
/* === super_symbolic_to_ll_super =========================================== */
/* ========================================================================== */

/* Convert a supernodal symbolic factorization to a supernodal numeric
 * factorization by allocating L->x.  Contents of L->x are undefined.
 */

static int super_symbolic_to_ll_super
(
    cholmod_factor *L,
    cholmod_common *Common
)
{
    double *Lx ;
    PRINT1 (("convert super sym to num\n")) ;
    ASSERT (L->ftype == CHOLMOD_SYMBOLIC_SUPER) ;
    /* don't use cholmod_calloc, even though L->x will be set to zero.
     * L->x is large, and it's probably better to incrementally set it to
     * zero, a single supernode at a time. */
    Lx = cholmod_malloc (L->xsize, sizeof (double), Common) ;
    PRINT1 (("xsize %d\n", L->xsize)) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;	/* out of memory */
    }

    /* ---------------------------------------------- commit the changes to L */

    L->x = Lx ;
    L->ftype = CHOLMOD_LL_SUPER ;
    L->xtype = Common->xtype ;
    L->dtype = Common->dtype ;
    L->minor = L->n ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_change_ftype ================================================= */
/* ========================================================================== */

/* Convert a factor L from any ftype to any other ftype (with a few exceptions).
 * Some conversions simply allocate uninitialized space that meant to be filled
 * later.
 *
 * If the conversion fails, the factor is left in its original form, with one
 * exception.  Converting a supernodal symbolic factor to a simplicial numeric
 * one (with L=D=I) may leave the factor in simplicial symbolic form.
 *
 * Memory allocated for each conversion is listed below.
 */

int cholmod_change_ftype
(
    /* inputs, modified on output: */
    cholmod_factor *L,

    /* inputs */
    int ftype,		/* what ftype of factor to convert L to */

    cholmod_common *Common
)
{
    int *Lp ;
    DEBUG (int orig) ;
    DEBUG (int forig) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    Common->status = CHOLMOD_OK ;

    DEBUG (forig = L->ftype) ;
    DEBUG (cholmod_dump_factor (L, "start change ftype", Common)) ;
    PRINT1 (("-----------------------------------convert from %d to ftype %d\n",
	    L->ftype, ftype)) ;

    if (ftype < CHOLMOD_SYMBOLIC_SUPER || ftype > CHOLMOD_LL_SUPER)
    {
	/* invalid ftype */
	cholmod_error (CHOLMOD_INVALID, "ftype invalid", Common) ;
	return (FALSE) ;
    }
    DEBUG (orig = Common->malloc_count) ;
    PRINT1 ((">>convert %d to %d\n", L->ftype, ftype)) ;

    /* ---------------------------------------------------------------------- */
    /* convert to the same ftype, or to symbolic or supernodal ftypes */
    /* ---------------------------------------------------------------------- */

    if (ftype == L->ftype)
    {

	/* ------------------------------------------------------------------ */
	/* converting to same ftype */
	/* ------------------------------------------------------------------ */

	if (ftype == CHOLMOD_LDL_PACKED || ftype == CHOLMOD_LL_PACKED)
	{
	    /* reduce L in size so that it is just big enough */
	    Lp = L->p ;
	    cholmod_reallocate_factor (L, Lp [L->n], Common) ;
	}
	ASSERT (L->ftype == ftype) ;			/* always succeeds */

    }
    else if (ftype == CHOLMOD_SYMBOLIC)
    {

	/* ------------------------------------------------------------------ */
	/* convert any factor into a simplicial symbolic factor */
	/* ------------------------------------------------------------------ */

	any_to_simplicial_symbolic (L, Common) ;
	ASSERT (L->ftype == ftype) ;			/* always succeeds */

    }
    else if (ftype == CHOLMOD_SYMBOLIC_SUPER)
    {

	/* ------------------------------------------------------------------ */
	/* convert to a supernodal symbolic factor */
	/* ------------------------------------------------------------------ */

	if (L->ftype == CHOLMOD_LL_SUPER)
	{
	    /* this preserves the symbolic pattern of L, discards numeric val */
	    ll_super_to_super_symbolic (L, Common) ;
	    ASSERT (L->ftype == ftype) ;		/* always succeeds */
	}
	else if (L->ftype == CHOLMOD_SYMBOLIC)
	{
	    /* contents of supernodal pattern are uninitialized */
	    simplicial_symbolic_to_super_symbolic  (L, Common) ;
	    ASSERT (L->ftype == forig			/* fully fails */
		 || L->ftype == ftype) ;		/* or succeeds */
	}
	else
	{
	    /* cannot convert from simplicial numeric to supernodal symbolic */
	    cholmod_error (CHOLMOD_INVALID,
		"cannot convert L to supernodal symbolic", Common) ;
	}

    }
    else if (ftype == CHOLMOD_LL_SUPER)
    {

	/* ------------------------------------------------------------------ */
	/* convert supernodal symbolic factor to a supernodal numeric factor */
	/* ------------------------------------------------------------------ */

	if (L->ftype == CHOLMOD_SYMBOLIC_SUPER)
	{
	    /* Contents of supernodal numeric values are uninitialized.  This
	     * is used by cholmod_super_numeric.  Not meant for the end user. */
	    super_symbolic_to_ll_super (L, Common) ;
	    ASSERT (L->ftype == forig			/* fully fails */
		 || L->ftype == ftype) ;		/* or succeeds */
	}
	else
	{
	    /* can only convert to numeric supernodal LL' from a symbolic
	     * supernodal LL' */
	    cholmod_error (CHOLMOD_INVALID,
		"cannot convert L to supernodal", Common) ;
	}

    }

    /* ---------------------------------------------------------------------- */
    /* convert any ftype to simplicial numeric ftypes */
    /* ---------------------------------------------------------------------- */

    else if (L->ftype == CHOLMOD_SYMBOLIC)
    {

	/* ------------------------------------------------------------------ */
	/* convert simplicial symbolic factorization to numeric (L=I,D=I) */
	/* ------------------------------------------------------------------ */

	/* contents of numerical values reflect L and D identity */
	simplicial_symbolic_to_simplicial_numeric (L, ftype, Common) ;
	ASSERT (L->ftype == forig			/* fully fails */
	     || L->ftype == ftype) ;			/* or succeeds */

    }
    else if (L->ftype == CHOLMOD_LL_SUPER)
    {

	/* ------------------------------------------------------------------ */
	/* convert a supernodal LL' to simplicial numeric */
	/* ------------------------------------------------------------------ */

	ll_super_to_simplicial_numeric (L, ftype, Common) ;
	ASSERT (L->ftype == forig			/* fully fails */
	     || L->ftype == ftype) ;			/* or succeeds */

    }
    else if (L->ftype == CHOLMOD_SYMBOLIC_SUPER)
    {

	/* ------------------------------------------------------------------ */
	/* convert a supernodal symbolic to simplicial numeric */
	/* ------------------------------------------------------------------ */

	/* contents of numerical values reflect identity */
	any_to_simplicial_symbolic (L, Common) ;	/* always succeeds */
	ASSERT (L->ftype == CHOLMOD_SYMBOLIC) ;
	simplicial_symbolic_to_simplicial_numeric (L, ftype, Common) ;
	ASSERT (L->ftype == forig			/* fully fails */
	     || L->ftype == CHOLMOD_SYMBOLIC		/* partial failure */
	     || L->ftype == ftype) ;			/* or succeeds */

    }
    else if (L->ftype == CHOLMOD_LDL_PACKED)
    {

	/* ------------------------------------------------------------------ */
	/* convert a packed LDL' to simplicial numeric */
	/* ------------------------------------------------------------------ */

	if (ftype == CHOLMOD_LDL_UNPACKED)
	{
	    ldl_packed_to_ldl_unpacked (L, Common) ;
	}
	else if (ftype == CHOLMOD_LDL_DYNAMIC)
	{
	    ldl_packed_to_ldl_dynamic (L, Common) ;
	}
	else if (ftype == CHOLMOD_LL_PACKED)
	{
	    ldl_packed_to_ll (L, Common) ;
	}
	ASSERT (L->ftype == forig			/* fully fails */
	     || L->ftype == ftype) ;			/* or succeeds */

    }
    else if (L->ftype == CHOLMOD_LDL_UNPACKED)
    {

	/* ------------------------------------------------------------------ */
	/* convert an unpacked LDL' to simplicial numeric */
	/* ------------------------------------------------------------------ */

	if (ftype == CHOLMOD_LDL_PACKED)
	{
	    ldl_unpacked_to_ldl_packed (L, Common) ;
	}
	else if (ftype == CHOLMOD_LDL_DYNAMIC)
	{
	    ldl_unpacked_to_ldl_dynamic (L, Common) ;
	}
	else if (ftype == CHOLMOD_LL_PACKED)
	{
	    ldl_unpacked_to_ldl_packed (L, Common) ;
	    ldl_packed_to_ll (L, Common) ;
	}
	ASSERT (L->ftype == forig			/* fully fails */
	     || L->ftype == ftype) ;			/* or succeeds */

    }
    else if (L->ftype == CHOLMOD_LDL_DYNAMIC)
    {

	/* ------------------------------------------------------------------ */
	/* convert a dynamic LDL' to simplicial numeric */
	/* ------------------------------------------------------------------ */

	if (ftype == CHOLMOD_LDL_UNPACKED)
	{
	    ldl_dynamic_to_ldl (L, ftype, Common) ;
	}
	else if (ftype == CHOLMOD_LDL_PACKED)
	{
	    ldl_dynamic_to_ldl (L, ftype, Common) ;
	}
	else if (ftype == CHOLMOD_LL_PACKED)
	{
	    /* this conversion cannot leave L in LDL' packed form */
	    ldl_dynamic_to_ldl (L, CHOLMOD_LDL_PACKED, Common) ;
	    if (Common->status < CHOLMOD_OK)
	    {
		return (FALSE) ;	/* out of memory */
	    }
	    ldl_packed_to_ll (L, Common) ;
	}
	ASSERT (L->ftype == forig			/* fully fails */
	     || L->ftype == ftype) ;			/* or succeeds */

    }
    else if (L->ftype == CHOLMOD_LL_PACKED)
    {

	/* ------------------------------------------------------------------ */
	/* convert a packed LL' to simplicial LDL' numeric */
	/* ------------------------------------------------------------------ */

	ll_to_ldl (L, ftype, Common) ;
	ASSERT (L->ftype == forig			/* fully fails */
	     || L->ftype == ftype) ;			/* or succeeds */
    }

    /* ---------------------------------------------------------------------- */
    /* return result */
    /* ---------------------------------------------------------------------- */

    PRINT1 (("Common->status %d L->ftype %d forig %d ftype %d\n"
		"orig: %d current %d mdelta %d\n",
		Common->status, L->ftype, forig, ftype, orig,
		Common->malloc_count, cholmod_dump_mdelta (forig, L->ftype))) ;
    ASSERT (Common->malloc_count == orig +
	    cholmod_dump_mdelta (forig, L->ftype)) ;
    return (Common->status >= CHOLMOD_OK) ;
}


/* ========================================================================== */
/* === cholmod_pack_factor ================================================== */
/* ========================================================================== */

/* Pack the columns of an LDL' unpacked or dynamic factor.  This can be followed
 * by a call to cholmod_reallocate_factor to reduce the size of L to the exact
 * size required by the factor, if desired.  Alternatively, you can leave the
 * size of L->i and L->x the same, to allow space for future updates/rowadds.
 *
 * Each column is reduced in size so that it has at most Common->grow2 free
 * space at the end of the column.
 *
 * Does nothing and returns silently if given any other ftype of factor
 * (LL' packed, LDL' packed, and LL' supernodal are always kept packed).
 */

int cholmod_pack_factor
(
    /* inputs, modified on output: */
    cholmod_factor *L,

    cholmod_common *Common
)
{
    int grow2 ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    Common->status = CHOLMOD_OK ;
    DEBUG (cholmod_dump_factor (L, "start pack", Common)) ;
    PRINT1 (("PACK factor %d\n", L->ftype)) ;

    grow2 = Common->grow2 ;

    /* ---------------------------------------------------------------------- */
    /* pack */
    /* ---------------------------------------------------------------------- */

    if (L->ftype == CHOLMOD_LDL_UNPACKED)
    {
	pack_ldl_unpacked (L, grow2) ;
    }
    else if (L->ftype == CHOLMOD_LDL_DYNAMIC)
    {
	pack_ldl_dynamic (L, grow2) ;
    }

    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_factor_to_sparse ============================================= */
/* ========================================================================== */

/* Constructs a column-oriented sparse matrix containing the pattern and values
 * of a simplicial numerical factor, and then converts the factor into a
 * simplicial symbolic factor.
 *
 * Can only convert factors of ftype LDL' or LL' packed.  If you want to convert
 * other ftypes, use cholmod_change_ftype first.
 *
 * Only a single malloc of size equal to the struct required for Lsparse is
 * allocated (about 52 bytes independent of the size of Lsparse), which is very
 * unlikely to fail.  This routine takes O(1) time, since all it does is
 * pointer manipulation on the L and Lsparse objects.  If this routine does
 * fail, L is left unmodified.
 */

cholmod_sparse *cholmod_factor_to_sparse
(
    /* input, converted to CHOLMOD_SYMBOLIC on output */
    cholmod_factor *L,
    cholmod_common *Common
)
{
    cholmod_sparse *Lsparse ;
    DEBUG (int orig) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (L, NULL) ;
    if (!(L->ftype == CHOLMOD_LL_PACKED || L->ftype == CHOLMOD_LDL_PACKED))
    {
	cholmod_error (CHOLMOD_INVALID,
	    "cholmod_factor_to_sparse: cannot convert L", Common) ;
	return (NULL) ;
    }
    Common->status = CHOLMOD_OK ;
    DEBUG (cholmod_dump_factor (L, "start convert to matrix", Common)) ;
    DEBUG (orig = Common->malloc_count) ;

    /* ---------------------------------------------------------------------- */
    /* create A */
    /* ---------------------------------------------------------------------- */

    /* allocate the header for Lsparse, the sparse matrix version of L */
    Lsparse = cholmod_malloc (1, sizeof (cholmod_sparse), Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	ASSERT (Common->malloc_count == orig) ;
	return (NULL) ;		/* out of memory */
    }

    /* transfer the contents from L to Lsparse */
    Lsparse->nrow = L->n ;
    Lsparse->ncol = L->n ;
    Lsparse->p = L->p ;
    Lsparse->i = L->i ;
    Lsparse->x = L->x ;
    Lsparse->z = L->z ;
    Lsparse->nz = NULL ;
    Lsparse->stype = 0 ;
    Lsparse->itype = L->itype ;
    Lsparse->xtype = L->xtype ;
    Lsparse->dtype = L->dtype ;
    Lsparse->sorted = TRUE ;
    Lsparse->packed = TRUE ;
    Lsparse->nzmax = L->nzmax ;
    ASSERT (L->nz == NULL) ;
    ASSERT (cholmod_dump_sparse (Lsparse, "Lsparse", Common) >= 0) ;

    /* convert L to symbolic, but do not free contents transfered to Lsparse */
    L->p = NULL ;
    L->i = NULL ;
    L->x = NULL ;
    L->z = NULL ;
    any_to_simplicial_symbolic (L, Common) ;

    ASSERT (Common->malloc_count == orig + 1) ;
    return (Lsparse) ;
}


/* ========================================================================== */
/* === cholmod_copy_factor ================================================== */
/* ========================================================================== */

/* Create an exact copy of a factor, with one exception:
 *
 * Entries in unused space are not copied (they might not be initialized,
 *	and copying them would cause program checkers such as purify and
 *	valgrind to complain).
 */

cholmod_factor *cholmod_copy_factor
(
    cholmod_factor *L,
    cholmod_common *Common
)
{
    cholmod_factor *L2 ;
    double *Lx, *L2x ;
    int *Perm, *ColCount, *Lp, *Li, *Lnz, *Lnext, *Lprev, *Lsuper, *Lpi, *Lpx,
	*Ls, *Perm2, *ColCount2, *L2p, *L2i, *L2nz, *L2next, *L2prev, *L2super,
	*L2pi, *L2px, *L2s ;
    int n, j, p, pend, s, nz, nzmax, ftype, xsize, ssize, nsuper ;
    DEBUG (int orig) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (L, NULL) ;

    Common->status = CHOLMOD_OK ;
    DEBUG (cholmod_dump_factor (L, "start copy", Common)) ;
    DEBUG (orig = Common->malloc_count) ;

    n = L->n ;
    ftype = L->ftype ;

    /* ---------------------------------------------------------------------- */
    /* allocate a simplicial symbolic factor  */
    /* ---------------------------------------------------------------------- */

    /* allocates L2->Perm and L2->ColCount and sets ftype to CHOLMOD_SYMBOLIC */
    L2 = cholmod_allocate_factor (n, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	ASSERT (Common->malloc_count == orig) ;
	return (NULL) ;	    /* out of memory */
    }

    Perm = L->Perm ;
    ColCount = L->ColCount ;
    Perm2 = L2->Perm ;
    ColCount2 = L2->ColCount ;
    L2->ordering = L->ordering ;

    for (j = 0 ; j < n ; j++)
    {
	Perm2 [j] = Perm [j] ;
    }
    for (j = 0 ; j < n ; j++)
    {
	ColCount2 [j] = ColCount [j] ;
    }

    /* ---------------------------------------------------------------------- */
    /* convert and copy the rest of the factor */
    /* ---------------------------------------------------------------------- */

    if (ftype >= CHOLMOD_LDL_PACKED && ftype <= CHOLMOD_LL_PACKED)
    {

	/* ------------------------------------------------------------------ */
	/* convert to a simplicial numeric factor */
	/* ------------------------------------------------------------------ */

	/* allocate L2->p, L2->nz, L2->prev, and L2->next */
	if (!allocate_simplicial_numeric (L2, ftype, Common))
	{
	    cholmod_free_factor (&L2, Common) ;
	    ASSERT (Common->malloc_count == orig) ;
	    return (NULL) ;	/* out of memory */
	}

	/* allocate L2->i and L2->x, of size nzmax */
	nzmax = L->nzmax ;
	if (!cholmod_reallocate_factor (L2, nzmax, Common))
	{
	    cholmod_free_factor (&L2, Common) ;
	    ASSERT (Common->malloc_count == orig) ;
	    return (NULL) ;	/* out of memory */
	}

	/* ------------------------------------------------------------------ */
	/* copy the contents of a simplicial numeric factor */
	/* ------------------------------------------------------------------ */

	Lp = L->p ;
	Li = L->i ;
	Lx = L->x ;
	Lnz = L->nz ;
	Lnext = L->next ;
	Lprev = L->prev ;

	L2p = L2->p ;
	L2i = L2->i ;
	L2x = L2->x ;
	L2nz = L2->nz ;
	L2next = L2->next ;
	L2prev = L2->prev ;
	L2->ftype = ftype ;
	L2->xtype = L->xtype ;
	L2->dtype = L->dtype ;

	for (j = 0 ; j <= n ; j++)
	{
	    L2p [j] = Lp [j] ;
	}

	if (ftype == CHOLMOD_LDL_DYNAMIC)
	{
	    for (j = 0 ; j < n+2 ; j++)
	    {
		L2prev [j] = Lprev [j] ;
	    }
	    for (j = 0 ; j < n+2 ; j++)
	    {
		L2next [j] = Lnext [j] ;
	    }
	}

	if (ftype == CHOLMOD_LDL_UNPACKED || ftype == CHOLMOD_LDL_DYNAMIC)
	{
	    for (j = 0 ; j < n ; j++)
	    {
		L2nz [j] = Lnz [j] ;
	    }
	    for (j = 0 ; j < n ; j++)
	    {
		p = Lp [j] ;
		pend = p + Lnz [j] ;
		for ( ; p < pend ; p++)
		{
		    L2i [p] = Li [p] ;
		    L2x [p] = Lx [p] ;
		}
	    }
	}
	else
	{
	    nz = Lp [n] ;
	    for (p = 0 ; p < nz ; p++)
	    {
		L2i [p] = Li [p] ;
	    }
	    for (p = 0 ; p < nz ; p++)
	    {
		L2x [p] = Lx [p] ;
	    }
	}

    }
    else if (ftype == CHOLMOD_SYMBOLIC_SUPER || ftype == CHOLMOD_LL_SUPER)
    {

	/* ------------------------------------------------------------------ */
	/* convert to a supernodal factor */
	/* ------------------------------------------------------------------ */

	xsize = L->xsize ;
	ssize = L->ssize ;
	nsuper = L->nsuper ;

	L2->xsize = xsize ;
	L2->ssize = ssize ;
	L2->nsuper = nsuper ;

	/* allocate L2->super, L2->pi, L2->px, and L2->s */
	if (!simplicial_symbolic_to_super_symbolic (L2, Common))
	{
	    cholmod_free_factor (&L2, Common) ;
	    ASSERT (Common->malloc_count == orig) ;
	    return (NULL) ;	/* out of memory */
	}

	/* allocate L2->x */
	if (ftype == CHOLMOD_LL_SUPER)
	{
	    if (!super_symbolic_to_ll_super (L2, Common))
	    {
		cholmod_free_factor (&L2, Common) ;
		ASSERT (Common->malloc_count == orig) ;
		return (NULL) ;	    /* out of memory */
	    }
	}

	/* ------------------------------------------------------------------ */
	/* copy the contents of a supernodal factor */
	/* ------------------------------------------------------------------ */

	Lsuper = L->super ;
	Lpi = L->pi ;
	Lpx = L->px ;
	Ls = L->s ;
	Lx = L->x ;

	L2super = L2->super ;
	L2pi = L2->pi ;
	L2px = L2->px ;
	L2s = L2->s ;
	L2x = L2->x ;

	L2->maxcsize = L->maxcsize ;
	L2->maxesize = L->maxesize ;

	for (s = 0 ; s <= nsuper ; s++)
	{
	    L2super [s] = Lsuper [s] ;
	}
	for (s = 0 ; s <= nsuper ; s++)
	{
	    L2pi [s] = Lpi [s] ;
	}
	for (s = 0 ; s <= nsuper ; s++)
	{
	    L2px [s] = Lpx [s] ;
	}

	L2s [0] = 0 ;
	for (p = 0 ; p < ssize ; p++)
	{
	    L2s [p] = Ls [p] ;
	}

	if (ftype == CHOLMOD_LL_SUPER)
	{
	    for (p = 0 ; p < xsize ; p++)
	    {
		L2x [p] = Lx [p] ;
	    }
	}
    }

    L2->minor = L->minor ;

    DEBUG (cholmod_dump_factor (L2, "L2", Common)) ;
    return (L2) ;
}
