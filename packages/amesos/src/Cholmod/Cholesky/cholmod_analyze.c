/* ========================================================================== */
/* === Cholesky/cholmod_analyze ============================================= */
/* ========================================================================== */

/*
 * CHOLMOD/Cholesky version 0.1, May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Authors: Timothy A. Davis
 * and Sivasankaran Rajamanickam
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* Order and analyze a matrix (either simplicial or supernodal), in prepartion
 * for numerical factorization via cholmod_factorize or via the "expert"
 * routines cholmod_rowfac and cholmod_super_numeric.
 *
 * symmetric case:    A or A(p,p)
 * unsymmetric case:  AA', A(p,:)*A(p,:)', A(:,f)*A(:,f)', or A(p,f)*A(p,f)'
 *
 * For the symmetric case, only the upper or lower triangular part of A is
 * accessed (depending on the type of A).  LL'=A (or permuted A) is analzed.
 * For the unsymmetric case (LL'=AA' or permuted A).
 *
 * There can be no duplicate entries in p or f.  p is of length m if A is
 * m-by-n.  f can be length 0 to n.
 *
 * In both cases, the columns of A need not be sorted.  A can be in packed
 * or unpacked form.
 *
 * Ordering options include:
 *
 *	natural:    A is not permuted to reduce fill-in
 *	given:	    a permutation can be provided to this routine (UserPerm)
 *	AMD:	    approximate minumum degree (AMD for the symmetric case,
 *		    COLAMD for the AA' case).
 *	METIS:	    nested dissection with METIS_NodeND
 *	ND:	    nested dissection using METIS_NodeComputeSeparator,
 *		    typically followed by a constrained minimum degree
 *		    (CSYMAMD for the symmetric case, CCOLAMD for the AA' case).
 *
 * Multiple ordering options can be tried (up to 8 of them), and the best one
 * is selected (the one that gives the smallest number of nonzeros in the
 * simplicial factor L).  If one method fails, cholmod_analyze keeps going, and
 * picks the best among the methods that succeeded.  This routine fails (and
 * returns NULL) if either initial memory allocation fails, all ordering methods
 * fail, or the supernodal analysis (if requested) fails.  By default, the 8
 * methods available are:
 *
 *	0) given permutation (skipped if UserPerm is NULL)
 *	1) AMD (symmetric case) or COLAMD (unsymmetric case)
 *	2) METIS with default parameters
 *	3) ND with default parameters (stopping the partitioning when
 *	    the graph is of size nd_small = 200 or less, remove nodes with
 *	    more than max (16, nd_prune * sqrt (n)) nodes where nd_prune = 10,
 *	    and follow partitioning with a constrained min. degree ordering).
 *	4) natural
 *	5) ND, nd_small = 20000, nd_prune = 10
 *	6) ND, nd_small =     4, nd_prune = 10, no min degree
 *	7) ND, nd_small =   200, nd_prune = 0
 *
 * By default, the first three are tried.  If you do not have METIS, only
 * the first two will be tried.  You can modify these 8 methods and the number
 * of methods tried by changing parameters in the Common argument.  If you know
 * the best ordering for your matrix, set Common->nmethods to 1 and set
 * Common->method[0].ordering to the requested ordering method.  Parameters
 * for each method can also be modified (refer to cholmod.h for details).
 *
 * Note that it is possible for METIS to terminate your program if it runs out
 * of memory.  This is not the case for any CHOLMOD or minimum degree ordering
 * routine (AMD, COLAMD, CCOLAMD, or CSYMAMD).
 *
 * The factor L is returned as type CHOLMOD_SYMBOLIC (a simplicial analysis) if
 * Common->supernodal is FALSE or as type CHOLMOD_SYMBOLIC_SUPER otherwise.
 * The default is to perform a supernodal analysis if the Supernodal module is
 * installed, or a simplicial analysis otherwise.  A subsequent call to
 * cholmod_factorize will perform a simplicial or supernodal factorization,
 * depending on the type of L.
 *
 * For the simplicial case, L contains the fill-reducing permutation (L->Perm)
 * and the counts of nonzeros in each column of L (L->ColCount).  For the
 * supernodal case, L also contains the nonzero pattern of each supernode.
 *
 * workspace: Flag (nrow), Head (nrow+1), Xwork (4*n*sizeof(int))
 *	if symmetric:   Iwork (2*nrow)
 *	if unsymmetric: Iwork (2*nrow+ncol).
 *	calls various ordering routines, which typically allocate O(nnz(A))
 *	temporary workspace ((2 to 3)*nnz(A) * sizeof (int) is typical, but it
 *	can be much higher if A*A' must be explicitly formed for METIS).  Also
 *	allocates up to 4 temporary (permuted/transpose) copies of the nonzero
 *	pattern of A, and up to 3*n*sizeof(int) additional workspace.
 */

#include "cholmod_cholesky.h"
#include "cholmod_internal.h"

#ifndef NSUPERNODAL
#include "cholmod_supernodal.h"
#endif

#ifndef NPARTITION
#include "cholmod_partition.h"
#endif


/* free the current temporary matrices, restore status, and try next method */
#define ERROR_CONTINUE \
{ \
    status = MIN (status, Common->status) ; \
    Common->status = CHOLMOD_OK ; \
    cholmod_free_sparse (&A1, Common) ; \
    cholmod_free_sparse (&A2, Common) ; \
    continue ; \
}


/* ========================================================================== */
/* === cholmod_analyze ====================================================== */
/* ========================================================================== */

/* Order and analyze A or AA'. */

cholmod_factor *cholmod_analyze
(
    cholmod_sparse *A,
    cholmod_common *Common
)
{
    return (cholmod_analyze_p (A, NULL, NULL, 0, Common)) ;
}


/* ========================================================================== */
/* === cholmod_analyze_p ==================================================== */
/* ========================================================================== */

/* Analyze A, AA', PAP', PAA'P, or A(p,f)*A(p,f)' */

cholmod_factor *cholmod_analyze_p
(
    /* inputs, not modified: */
    cholmod_sparse *A,
    void *UserPerm_p,
    void *fset,			/* the set f.  NULL means use all cols of A */
    size_t fsize,		/* the size of fset */

    /* input parameters, workspace, and output statistics */
    cholmod_common *Common
)
{
    double lnz_best ;
    double *W ;
    int *First, *Level, *Work4n, *Cmember, *CParent, *ColCount, *Lperm,
	*Parent, *Post, *Perm, *Parent_best, *UserPerm, *Lcolcount ;
    cholmod_sparse *S, *F, *A2, *A1, *A1_best, *A2_best, *F_best, *S_best ;
    cholmod_factor *L ;
    int k, symmetric, n, ordering, ncol, method, wsize, nmethods, status, nf,
	ok ;

    DEBUG (int orig) ;
    nf = fsize ;
    UserPerm = UserPerm_p ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    Common->status = CHOLMOD_OK ;
    status = CHOLMOD_OK ;
    Common->selected = EMPTY ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    n = A->nrow ;
    ncol = A->ncol ;
    symmetric = A->stype ;
    lnz_best = (double) EMPTY ;
    nmethods = MIN (Common->nmethods, CHOLMOD_MAXMETHODS) ;
    nmethods = MAX (0, nmethods) ;
    PRINT1 (("nmethods %d\n", nmethods)) ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    /* Note: enough space needs to be allocated here so that routines called by
     * cholmod_analyze do not reallocate the space.
     */

    wsize = cholmod_allocate_work (n, 2*n + (symmetric ? 0: ncol), 4*n,
	    sizeof (int), Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory */
    }
    ASSERT (cholmod_dump_work (TRUE, TRUE, 0, Common)) ;

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    /* Do not use Iwork because etree, postorder, rowcolcounts, and many
     * of the ordering methods need it.  Use Xwork instead.  These arguments
     * are passed to other CHOLMOD routines, which themselves use space from
     * Common.
     */

    W = Common->Xwork ;
    Work4n = Common->Xwork ;
    Parent = Work4n ;
    First  = Work4n + n ;
    Level  = Work4n + 2*n ;
    Post   = Work4n + 3*n ;
    Cmember = Post ;
    CParent = Level ;

    DEBUG (orig = Common->malloc_count) ;
    PRINT1 (("orig: %d\n", orig)) ;

    /* ---------------------------------------------------------------------- */
    /* permuted and/or transposed copies of the input matrix (pattern only) */
    /* ---------------------------------------------------------------------- */

    S = NULL ;
    F = NULL ;
    A1 = NULL ;
    A2 = NULL ;
    A1_best = NULL ;
    A2_best = NULL ;
    S_best = NULL ;
    F_best = NULL ;

    /* ---------------------------------------------------------------------- */
    /* allocate an empty factor object of type CHOLMOD_SYMBOLIC */
    /* ---------------------------------------------------------------------- */

    L = cholmod_allocate_factor (n, Common) ;
    Perm = cholmod_malloc (n, sizeof (int), Common) ;
    ColCount = cholmod_malloc (n, sizeof (int), Common) ;
    Parent_best = NULL ;
#ifndef NSUPERNODAL
    if (Common->supernodal)
    {
	Parent_best = cholmod_malloc (n, sizeof (int), Common) ;
    }
#else
    /* CHOLMOD Supernodal module not installed, just do simplicial analysis */
    Common->supernodal = FALSE ;
#endif
    if (Common->status < CHOLMOD_OK)
    {
	/* out of memory: set nmethods to EMPTY so no methods are tried */
	nmethods = EMPTY ;
	status = Common->status ;	/* out of memory */
    }

    /* ---------------------------------------------------------------------- */
    /* try all the requested ordering options and backup to AMD if needed */
    /* ---------------------------------------------------------------------- */

    /* turn off error printing */
    Common->try_catch = TRUE ;

    for (method = 0 ; method <= nmethods ; method++)
    {

	/* ------------------------------------------------------------------ */
	/* determine the method to try */
	/* ------------------------------------------------------------------ */

	if (method == nmethods)
	{
	    /* All methods failed: backup to AMD */
	    if (Common->selected == EMPTY)
	    {
		PRINT1 (("All methods requested failed: backup to AMD\n")) ;
		ordering = CHOLMOD_AMD ;
	    }
	    else
	    {
		break ;
	    }
	}
	else
	{
	    ordering = Common->method [method].ordering ;
	}
	Common->current = method ;
	PRINT1 (("method %d: Try ordering method: %d\n", method, ordering)) ;
	ASSERT (A1 == NULL && A2 == NULL) ;

	/* ------------------------------------------------------------------ */
	/* find the fill-reducing permutation */
	/* ------------------------------------------------------------------ */

	if (ordering == CHOLMOD_NATURAL)
	{

	    /* -------------------------------------------------------------- */
	    /* natural ordering */
	    /* -------------------------------------------------------------- */

	    for (k = 0 ; k < n ; k++)
	    {
		Perm [k] = k ;
	    }

	}
	else if (ordering == CHOLMOD_GIVEN)
	{

	    /* -------------------------------------------------------------- */
	    /* use given ordering of A, if provided */
	    /* -------------------------------------------------------------- */

	    if (UserPerm == NULL)
	    {
		/* this is not an error condition */
		PRINT1 (("skip, no user perm given\n")) ;
		continue ;
	    }
	    for (k = 0 ; k < n ; k++)
	    {
		/* Perm is checked in cholmod_transpose */
		Perm [k] = UserPerm [k] ;
	    }

	}
	else if (ordering == CHOLMOD_AMD)
	{

	    /* -------------------------------------------------------------- */
	    /* AMD or COLAMD ordering of A */
	    /* -------------------------------------------------------------- */

	    if (symmetric)
	    {
		/* use AMD with no constraints */
		/* workspace: Iwork (2*nrow), Xwork (4*nrow*sizeof(int)),
		 * Head (nrow+1).  Note that this routine uses the same space
		 * (Common->Xwork) used here for the Parent, First, Level, and
		 * Post arrays. */
		(void) cholmod_amd (A, Perm, Common) ;
	    }
	    else
	    {
		/* use COLAMD with no constraints */
		/* workspace: Iwork (nrow+MAX(nrow,ncol)) */
		(void) cholmod_colamd (A, fset, nf, Perm, Common) ;
	    }

	}
	else if (ordering == CHOLMOD_METIS)
	{

	    /* -------------------------------------------------------------- */
	    /* use METIS_NodeND directly (via a CHOLMOD wrapper) */
	    /* -------------------------------------------------------------- */

#ifndef NPARTITION
	    /* workspace: Iwork (nrow).  If A unsymmetric: Flag (nrow) */
	    (void) cholmod_metis (A, fset, nf, Perm, Common) ;
#else
	    Common->status = CHOLMOD_NOT_INSTALLED ;
#endif

	}
	else if (ordering == CHOLMOD_ND)
	{

	    /* -------------------------------------------------------------- */
	    /* use CHOLMOD's nested dissection */
	    /* -------------------------------------------------------------- */

	    /* this method is based on METIS' node bissection routine
	     * (METIS_NodeComputeSeparator).  In contrast to METIS_NodeND,
	     * it calls CSYMAMD or CCOLAMD on the whole graph, instead of MMD
	     * on just the leaves. */

#ifndef NPARTITION
	    /* workspace: Flag (nrow), Head (nrow+1), Iwork (2*nrow) */
	    cholmod_nested_dissection (A, fset, nf, Perm, CParent, Cmember,
		    Common) ;
#else
	    Common->status = CHOLMOD_NOT_INSTALLED ;
#endif

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* invalid ordering method */
	    /* -------------------------------------------------------------- */

	    Common->status = CHOLMOD_INVALID ;
	    PRINT1 (("No such ordering: %d\n", ordering)) ;
	}

	ASSERT (cholmod_dump_work (TRUE, TRUE, 0, Common)) ;

	if (Common->status < CHOLMOD_OK)
	{
	    ERROR_CONTINUE ;	    /* out of memory, or method failed */
	}

	/* ------------------------------------------------------------------ */
	/* permute the matrix and get its permuted transpose */
	/* ------------------------------------------------------------------ */

	if (ordering == CHOLMOD_NATURAL)
	{

	    /* -------------------------------------------------------------- */
	    /* natural ordering of A */
	    /* -------------------------------------------------------------- */

	    if (symmetric < 0)
	    {
		/* symmetric lower case: A already in lower form, so let S=A' */
		A1 = NULL ;
		/* workspace: Iwork (nrow) */
		A2 = cholmod_transpose (A, FALSE, NULL, NULL, 0, Common) ;
		F = A ;
		S = A2 ;
	    }
	    else
	    {
		/* symmetric upper case: F = pattern of triu (A)' and S = A */
		/* unsymmetric case:     F = pattern of A (:,f)'  and S = A */
		/* workspace:
		 * unsym: Iwork (nrow if no fset, MAX(nrow,ncol) if fset)
		 * sym:   Iwork (nrow) */
		A1 = cholmod_transpose (A, FALSE, NULL, fset, nf, Common) ;
		A2 = NULL ;
		F = A1 ;
		S = A ;
	    }

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* A is permuted */
	    /* -------------------------------------------------------------- */

	    if (symmetric < 0)
	    {
		/* symmetric lower case: S = tril (A (p,p))' and F = S' */
		/* workspace: Iwork (2*nrow) */
		A2 = cholmod_transpose (A, FALSE, Perm, NULL, 0, Common) ;
		S = A2 ;
		/* workspace: Iwork (nrow) */
		A1 = cholmod_transpose (A2, FALSE, NULL, NULL, 0, Common) ;
		F = A1 ;
	    }
	    else
	    {
		/* symmetric upper case: F = triu (A (p,p))'  and S = F' */
		/* unsymmetric case:     F = A (p,f)'         and S = F' */
		/* workspace: symmetric: Iwork (2*nrow),
		 * unsym: Iwork (nrow if no fset, MAX(nrow,ncol) if fset)
		 */
		A1 = cholmod_transpose (A, FALSE, Perm, fset, nf, Common) ;
		F = A1 ;
		/* workspace: Iwork (nrow) */
		A2 = cholmod_transpose (A1, FALSE, NULL, NULL, 0, Common) ;
		S = A2 ;
	    }
	}

	/* transpose may have failed */
	ok = (Common->status == CHOLMOD_OK) ;

	/* ------------------------------------------------------------------ */
	/* find etree of S (symmetric upper/lower case) or F (unsym. case) */
	/* ------------------------------------------------------------------ */

	/* workspace: symmmetric: Iwork (nrow), unsymmetric: Iwork (nrow+ncol)*/
	ok = ok && cholmod_etree (symmetric ? S:F, Parent, Common) ;

	/* ------------------------------------------------------------------ */
	/* postorder the etree */
	/* ------------------------------------------------------------------ */

	/* workspace: Iwork (2*nrow) */
	ok = ok && (cholmod_postorder (Parent, n, Post, Common) == n) ;

	/* ------------------------------------------------------------------ */
	/* analyze LL'=S or SS' or S(:,f)*S(:,f)' */
	/* ------------------------------------------------------------------ */

	/* workspace:
	 *	if symmetric:   Flag (nrow), Iwork (2*nrow)
	 *	if unsymmetric: Flag (nrow), Iwork (2*nrow+ncol), Head (nrow+1)
	 */
	ok = ok && cholmod_rowcolcounts (symmetric ? F:S, fset, nf, Parent,
		Post, NULL, ColCount, First, Level, Common) ;

	if (!ok)
	{
	    ERROR_CONTINUE ;  /* out of memory, Perm invalid, or fset invalid */
	}

	Common->method [method].fl  = Common->fl ;
	Common->method [method].lnz = Common->lnz ;
	PRINT1 (("lnz %g fl %g\n", Common->lnz, Common->fl)) ;
	ASSERT (Common->lnz >= 0.0) ;

	/* ------------------------------------------------------------------ */
	/* pick the best method */
	/* ------------------------------------------------------------------ */

	/* fl.pt. compare, but lnz can never be NaN */
	if (Common->selected == EMPTY || Common->lnz < lnz_best)
	{
	    Common->selected = method ;
	    PRINT1 (("this is best so far, method %d\n", method)) ;
	    L->ordering = ordering ;
	    Lperm = L->Perm ;
	    Lcolcount = L->ColCount ;
	    for (k = 0 ; k < n ; k++)
	    {
		Lperm [k] = Perm [k] ;
	    }
	    for (k = 0 ; k < n ; k++)
	    {
		Lcolcount [k] = ColCount [k] ;
	    }
	    lnz_best = Common->lnz ;

	    if (Common->supernodal)
	    {
		/* keep a copy of Parent, S, and F for supernodal analysis */
		for (k = 0 ; k < n ; k++)
		{
		    Parent_best [k] = Parent [k] ;
		}

		/* free the prior permuted copies of the input matrix, if any */
		cholmod_free_sparse (&A1_best, Common) ;
		cholmod_free_sparse (&A2_best, Common) ;

		/* A1 and A2 are not freed, but copied to A1_best and A2_best */
		A1_best = A1 ;
		A2_best = A2 ;
		F_best = F ;
		S_best = S ;
		A1 = NULL ;
		A2 = NULL ;
	    }
	}

	/* A1, A2 matrices no longer needed (does nothing if A1, A2 are NULL) */
	cholmod_free_sparse (&A1, Common) ;
	cholmod_free_sparse (&A2, Common) ;
    }

    /* turn error printing back on */
    Common->try_catch = FALSE ;

    /* ---------------------------------------------------------------------- */
    /* supernodal analysis, if requested and installed */
    /* ---------------------------------------------------------------------- */

    if (Common->selected == EMPTY)
    {
	/* All methods failed.  
	 * If two or more methods failed, they may have failed for different
	 * reasons.  Both would clear Common->status and skip to the next method
	 * via ERROR_CONTINUE.  Common->status needs to be restored here to the
	 * worst error obtained in any of the methods.  CHOLMOD_INVALID is worse
	 * than CHOLMOD_OUT_OF_MEMORY, since the former implies something may
	 * be wrong with the user's input.  CHOLMOD_OUT_OF_MEMORY is simply an
	 * indication of lack of resources. */
	ASSERT (status < CHOLMOD_OK) ;
	cholmod_error (status, "cholmod_analyze: all methods failed", Common) ;
    }
    else
    {
	/* At least one method succeeded. */
	Common->fl  = Common->method [Common->selected].fl  ;
	Common->lnz = Common->method [Common->selected].lnz ;
	ASSERT (Common->lnz >= 0) ;
#ifndef NSUPERNODAL
	if (Common->supernodal)
	{
	    /* workspace: Flag (nrow), Head (nrow), Iwork (2*nrow) */
	    cholmod_super_symbolic (S_best, F_best, Parent_best, L, Common) ;
	    PRINT1 (("status %d L->ftype %d\n", Common->status, L->ftype)) ;
	}
#endif
    }

    /* ---------------------------------------------------------------------- */
    /* free temporary matrices and workspace, and return result */
    /* ---------------------------------------------------------------------- */

    cholmod_free_sparse (&A1, Common) ;
    cholmod_free_sparse (&A2, Common) ;
    cholmod_free_sparse (&A1_best, Common) ;
    cholmod_free_sparse (&A2_best, Common) ;
    cholmod_free (Parent_best, n, sizeof (int), Common) ;
    cholmod_free (Perm,        n, sizeof (int), Common) ;
    cholmod_free (ColCount,    n, sizeof (int), Common) ;
    ASSERT (W == Common->Xwork) ;
    for (k = 0 ; k < wsize ; k++)
    {
	W [k] = 0. ;
    }
    if (Common->status < CHOLMOD_OK)
    {
	cholmod_free_factor (&L, Common) ;
    }
    PRINT1 (("orig: %d malloc_count %d\n", orig, Common->malloc_count)) ;
    ASSERT (Common->malloc_count == orig +
	    (L == NULL) ? 0 : (3 + ((Common->supernodal ? 4:0)))) ;
    ASSERT (cholmod_dump_work (TRUE, TRUE, -1, Common)) ;
    return (L) ;
}
