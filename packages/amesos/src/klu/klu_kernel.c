/* ========================================================================== */
/* === klu kernel =========================================================== */
/* ========================================================================== */

/* Sparse left-looking LU factorization, with partial pivoting.  Based on
 * Gilbert & Peierl's method, with or without Eisenstat & Liu's symmetric
 * pruning.  No user-callable routines are in this file.
 *
 * TODO: NaN handling.
 * TODO: Check for tiny pivots (not just zero pivots).
 * TODO: make GROWTH a run-time control parameter
 * TODO: keep track of diagonal when pivoting (= DiagonalMap in UMFPACK)
 */

#include "klu.h"
#include "klu_kernel.h"

#define GROWTH 1.2
#define PIVOTAL(i) (Pinv [i] >= 0)

/* force pruning on, and use non-recursive depth-first-search */
#ifndef PRUNE
#define PRUNE
#endif
#ifdef RECURSIVE
#undef RECURSIVE
#endif

/* ========================================================================== */
/* === GET_COLUMN =========================================================== */
/* ========================================================================== */

/* Given a row index j, find the corresponding pivotal column jnew, and return
 * the start and end indices into Li that contains the pattern of the column
 * of L.  With pruning, some of the column of L is not needed. */

#ifdef PRUNE
#define GET_COLUMN(j,p,pend) \
{ \
    int jnew = Pinv [j] ; \
    pend = Lpend [jnew] ; \
    p = Lp [jnew] ; \
}
#else
#define GET_COLUMN(j,p,pend) \
{ \
    int jnew = Pinv [j] ; \
    pend = Lp [jnew+1] ; \
    p = Lp [jnew] ; \
}
#endif

/* ========================================================================== */
/* === dfs ================================================================== */
/* ========================================================================== */

/* Does a depth-first-search, starting at node j.  The dfs routine is private
 * to this file.  L is unit diagonal, but this routine also works if the
 * diagonal of L is not present.  The columns of L need not be sorted. */

static int dfs		/* returns top of stack on output */
(
    /* input, not modified on output: */
    int j,		/* node at which to start the DFS */
    int k,		/* mark value, for the Flag array */
    int n,		/* L is n-by-n */
    int Lp [ ],		/* size n+1, column pointers of L */
    int Li [ ],		/* size lnz = Lp [n], row indices of L.  Diagonal (if
			 * present) is ignored.  L must be lower triangular. */
    int Pinv [ ],	/* Pinv [i] = k if row i is kth pivot row, or EMPTY if
			 * row i is not yet pivotal. */
    int top,		/* top of stack on input */

    /* input/output: */
    int Stack [ ],	/* stack containing the pattern of x */
    int Flag [ ],	/* Flag [i] == k means i is marked */

    /* other */
    WorkStackType WorkStack [ ],    /* workspace for non-recursive version */
    int Lpend [ ]		    /* for symmetryc pruning */
)
{
    int i, p, pend ;

#ifdef RECURSIVE

    /* ====================================================================== */
    /* === recursive DFS ==================================================== */
    /* ====================================================================== */

    /* mark j as having been visited */
    PRINTF (("[ start dfs at %d : new %d\n", j, Pinv [j])) ;
    Flag [j] = k ;

    /* start a DFS at this node j.  L (:,j) doesn't exist if j not pivotal*/
    if (PIVOTAL (j))
    {
	GET_COLUMN (j, p, pend) ;
	for ( ; p < pend ; p++)
	{
	    i = Li [p] ;
	    if (Flag [i] != k)
	    {
		top = dfs (i, k, n, Lp, Li, Pinv, top, Stack, Flag,
			WorkStack, Lpend) ;
	    }
	}
    }

    /* push j on the stack and return the top of the stack */
    PRINTF (("  end   dfs at %d ] top: %d\n", j, top-1)) ;
    Stack [--top] = j ;

#else

    /* ====================================================================== */
    /* === non-recursive DFS ================================================ */
    /* ====================================================================== */

    /* WorkStack is empty */
    int wtop = n ;

    PRINTF (("[ start dfs at %d : new %d\n", j, Pinv [j])) ;
    Flag [j] = k ;
    if (PIVOTAL (j))
    {
	/* push [j p pend] onto work stack to start the recursion */
	GET_COLUMN (j, p, pend) ;
	wtop-- ;
	WorkStack [wtop].j = j ;
	WorkStack [wtop].p = p ;
	WorkStack [wtop].pend = pend ;

	/* perform the recursion until the workstack is empty */
	while (wtop < n)
	{
	    /* continue iteration for column j */
	    p    = WorkStack [wtop].p ;
	    pend = WorkStack [wtop].pend ;

	    for ( ; p < pend ; p++)
	    {
		j = Li [p] ;
		if (Flag [j] != k)
		{
		    PRINTF (("[ start dfs at %d : new %d\n", j, Pinv [j])) ;
		    Flag [j] = k ;
		    if (PIVOTAL (j))
		    {
			/* copy p back into the top of workstack */
			WorkStack [wtop].p = p ;
			ASSERT (pend == WorkStack [wtop].pend) ;

			/* begin recursion for new column j of L */
			GET_COLUMN (j, p, pend) ;
			wtop-- ;
			WorkStack [wtop].j = j ;
			WorkStack [wtop].p = p ;
			WorkStack [wtop].pend = pend ;
		    }
		    else
		    {
			/* j is not pivotal, just push j onto output stack */
			PRINTF (("  end   dfs at %d ] top: %d\n", j, top-1)) ;
			Stack [--top] = j ;
		    }
		}
	    }

	    /* pivotal column j, on top of stack, is done. pop from WorkStack */
	    ASSERT (pend == WorkStack [wtop].pend) ;
	    j = WorkStack [wtop++].j ;

	    /* push j onto output stack */
	    PRINTF (("  end   dfs at %d ] top: %d\n", j, top-1)) ;
	    Stack [--top] = j ;

	    /* continue prior iteration, for predecessor of j */
	}
    }
    else
    {
	/* j is not pivotal, push j on the stack and finish */
	PRINTF (("  end   dfs at %d ] top: %d\n", j, top-1)) ;
	Stack [--top] = j ;
    }

#endif

    /* return top of output stack */
    return (top) ;
}

/* ========================================================================== */
/* === lsolve_symbolic ====================================================== */
/* ========================================================================== */

/* Finds the pattern of x, for the solution of Lx=b */

static int lsolve_symbolic	/* returns top of stack */
(
    /* input, not modified on output: */
    int k,		/* mark value, for the Flag array */
    int n,		/* L is n-by-n, where n >= 0 */
    int Lp [ ],		/* of size n+1, column pointers of L */
    int Li [ ],		/* size lnz=Lp[n], row indices of L */
    int Pinv [ ],	/* Pinv [i] = k if i is kth pivot row, or EMPTY if row i
			 * is not yet pivotal */
    int nzb,		/* nz in b */
    int Bi [ ],		/* pattern of b, of size nzb */

    /* output, not defined on input */
    int Stack [ ],	/* size n.  Pattern of x on output is in
			 * Stack [top .. n-1]. */

    /* workspace, defined on input and output */
    int Flag [ ],	/* size n.  Initially, all of Flag [0..n-1] < k.  After
			 * lsolve_symbolic is done, Flag [i] == k if i is in
			 * the pattern of the output, Stack [top .. n-1],
			 * and Flag [0..n-1] <= k. */

    /* other */
    WorkStackType WorkStack [ ],    /* workspace for non-recursive DFS */
    int Lpend [ ]		    /* for symmetric pruning */
)
{
    int i, p, top ;

#ifndef NDEBUG
    for (i = 0 ; i < n ; i++) ASSERT (Flag [i] < k) ;
#endif

    /* compute the pattern of X */
    top = n ;
    for (p = 0 ; p < nzb ; p++)
    {
	i = Bi [p] ;
	PRINTF (("\n ===== DFS at node %d in b, Pinv[i]: %d\n", i, Pinv [i])) ;
	if (Flag [i] != k)
	{
	    top = dfs (i, k, n, Lp, Li, Pinv, top, Stack, Flag,
		WorkStack, Lpend) ;
	}
    }

#ifndef NDEBUG
    for (p = top ; p < n ; p++)
    {
	PRINTF (("pattern of X: %d : %d pivot row: %d\n",
	    p, Stack [p], Pinv [Stack [p]])) ;
	ASSERT (Flag [Stack [p]] == k) ;
    }
    p = 0 ;
    for (i = 0 ; i < n ; i++)
    {
	ASSERT (Flag [i] <= k) ;
	if (Flag [i] == k) p++ ;
    }
    ASSERT (p == n - top) ;
#endif

    return (top) ;
}

/* ========================================================================== */
/* === lsolve_numeric ======================================================= */
/* ========================================================================== */

/* Computes the numerical values of x, for the solution of Lx=b.  Note that x
 * may include explicit zeros if numerical cancelation occurs.  L is assumed
 * to be unit-diagonal, with possibly unsorted columns (but the first entry in
 * the column must always be the diagonal entry). */

static void lsolve_numeric
(

    /* input, not modified on output: */
    int n,		/* L is n-by-n, where n >= 0 */
    int Lp [ ],		/* of size n+1, column pointers of L */
    int Li [ ],		/* size lnz=Lp[n], row indices of L */
    double Lx [ ],	/* size lnz=Lp[n], values of L  */
    int Pinv [ ],	/* Pinv [i] = k if i is kth pivot row, or EMPTY if row i
			 * is not yet pivotal */
    int nzb,		/* nz in b */
    int Bi [ ],		/* pattern of b, of size nzb */
    double Bx [ ],	/* values of b, of size nzb */
    int top,		/* top of stack containing the pattern of X */
    int Stack [ ],	/* pattern of x is in Stack [top .. n-1] */

    /* output, must be zero on input: */
    double X [ ]	/* size n, initially zero.  On output,
			 * X [Stack [top..n-1]] contains the solution. */
)
{
    int p, s, j, jnew, pend ;
    double xj ;

#ifndef NDEBUG
    for (j = 0 ; j < n ; j++) ASSERT (X [j] == 0) ;
#endif

    /* scatter b into x */
    for (p = 0 ; p < nzb ; p++)
    {
	X [Bi [p]] = Bx [p] ;
    }
    /* solve Lx=b */
    for (s = top ; s < n ; s++)
    {
	/* forward solve with column j of L (skip the unit diagonal of L) */
	j = Stack [s] ;
	if (!PIVOTAL (j)) continue ;
	xj = X [j] ;
	jnew = Pinv [j] ;
	ASSERT (Lp [jnew] < Lp [jnew+1] && Li [Lp [jnew]] == j) ;
	pend = Lp [jnew+1] ;
	for (p = Lp [jnew] + 1 ; p < pend ; p++)
	{
	    X [Li [p]] -= Lx [p] * xj ;
	}
    }

#ifndef NDEBUG
    for (p = top ; p < n ; p++) PRINTF (("%d : %g\n", Stack [p], X [Stack[p]]));
#endif
}

/* ========================================================================== */
/* === klu_kernel =========================================================== */
/* ========================================================================== */

/* returns KLU_OK (0), KLU_SINGULAR (-1), or KLU_OUT_OF_MEMORY (-2) */

#if 0
#ifdef RECURSIVE
#ifdef PRUNE
int klu_prune			/* with symmetric pruning, recursive DFS */
#else
int klu_noprune			/* no symmetric pruning, recursive DFS */
#endif
#else
#ifdef PRUNE
int klu_prune_nonrecursive	/* with symmetric pruning, non-recursive DFS */
				/* this is the default method */
#else
int klu_noprune_nonrecursive	/* no symmetric pruning, non-recursive DFS */
#endif
#endif
#else

int klu_kernel	    /* with pruning, and non-recursive depth-first-search */

#endif
(
    /* input, not modified */
    int n,	    /* A is n-by-n */
    int Ap [ ],	    /* size n+1, column pointers for A */
    int Ai [ ],	    /* size nz = Ap [n], row indices for A */
    double Ax [ ],  /* size nz, values of A */
    int Q [ ],	    /* size n, optional input permutation */
    double tol,	    /* partial pivoting tolerance parameter */

    /* input and output: */
    int *p_lsize,   /* input: initial size of Li and Lx, output: final size */
    int *p_usize,   /* input: initial size of Ui and Ux, output: final size */

    /* output, not defined on input */
    int Lp [ ],	    /* size n+1, column pointers for L */
    int **p_Li,	    /* size lsize, row indices of L */
    double **p_Lx,  /* size lsize, values of L */
    int Up [ ],	    /* size n+1, column pointers for U */
    int **p_Ui,	    /* size usize, row indices of U */
    double **p_Ux,  /* size usize, values of U */
    int Pinv [ ],   /* size n, inverse row permutation, where Pinv [i] = k if
		     * row i is the kth pivot row */
    int P [ ],	    /* size n, row permutation, where P [k] = i if row i is the
		     * kth pivot row. */
    int *p_noffdiag,	/* # of off-diagonal pivots chosen */

    /* workspace, not defined on input or output */
    double X [ ],   /* size n */
    int Stack [ ],  /* size n */
    int Flag [ ],   /* size n */

    /* other workspace: */
    WorkStackType WorkStack [ ],    /* size n, for non-recursive DFS only */
    int Lpend [ ],		    /* size n workspace, for pruning only */
    int Lpruned [ ]		    /* size n workspace, for pruning only */
)
{ 
    int lp, up, k, nzb, *Bi, top, p, lnz, unz, i, pivrow, found, kbar,
	*Li, *Ui, lsize, usize, ok, oki, okx, kcol, diagrow, noffdiag ;
    double pivot, *Bx, *Lx, *Ux, x, diag, abs_pivot, abs_diag, abs_x, xsize ;

    /* ---------------------------------------------------------------------- */
    /* get initial Li, Lx, Ui, and Ux */
    /* ---------------------------------------------------------------------- */

    lsize = *p_lsize ;
    usize = *p_usize ;
    PRINTF (("input: lsize %d usize %d\n", lsize, usize)) ;

    if (lsize < 0 || usize < 0)
    {
	return (KLU_OUT_OF_MEMORY) ;
    }

    Li = *p_Li ;
    Lx = *p_Lx ;
    Ui = *p_Ui ;
    Ux = *p_Ux ;

    /* ---------------------------------------------------------------------- */
    /* initializations */
    /* ---------------------------------------------------------------------- */

    noffdiag = 0 ;
    lp = 0 ;
    up = 0 ;
    for (k = 0 ; k < n ; k++)
    {
	X [k] = 0 ;
	Flag [k] = EMPTY ;
#ifdef PRUNE
	Lpruned [k] = FALSE ;	/* flag k as not pruned */
#endif
    }

    /* mark all rows as non-pivotal, and find the row indices of the initial
     * diagonal entries for each column */
    if (Q == (int *) NULL)
    {
	for (k = 0 ; k < n ; k++)
	{
	    P [k] = k ;
	    Pinv [k] = FLIP (k) ;	/* mark all rows as non-pivotal */
	}
    }
    else
    {
	for (k = 0 ; k < n ; k++)
	{
	    /* Assume Q is applied symmetrically to the system, so the initial
	     * P is equal to Q.  This can change via partial pivoting. */
	    P [k] = Q [k] ;
	    Pinv [Q [k]] = FLIP (k) ;	/* mark all rows as non-pivotal */
	}
    }

    /* P [k] = row means that UNFLIP (Pinv [row]) = k, and visa versa.
     * If row is pivotal, then Pinv [row] >= 0.  A row is initially "flipped"
     * (Pinv [k] < EMPTY), and then marked "unflipped" when it becomes
     * pivotal. */

#ifndef NDEBUG
    for (k = 0 ; k < n ; k++)
    {
	PRINTF (("Initial P [%d] = %d\n", k, P [k])) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* factorize */
    /* ---------------------------------------------------------------------- */

    for (k = 0 ; k < n ; k++)
    {
	PRINTF (("\n\n==================================== k: %d\n", k)) ;
	Lp [k] = lp ;
	Up [k] = up ;

	/* ------------------------------------------------------------------ */
	/* get the kth column of A and find the "diagonal" */
	/* ------------------------------------------------------------------ */

	if (Q == (int *) NULL)
	{
	    kcol = k ;
	}
	else
	{
	    kcol = Q [k] ;
	}

	nzb = Ap [kcol+1] - Ap [kcol] ;
	Bi = Ai + Ap [kcol] ;
	Bx = Ax + Ap [kcol] ;

	diagrow = P [k] ;   /* might already be pivotal */
	PRINTF (("k %d kcol %d, diagrow = %d, UNFLIP(diagrow) = %d\n",
	    k, kcol, diagrow, UNFLIP (diagrow))) ;

	/* ------------------------------------------------------------------ */
	/* compute the kth column of L and U */
	/* ------------------------------------------------------------------ */

	top = lsolve_symbolic (k, n, Lp, Li, Pinv, nzb, Bi, Stack, Flag,
		WorkStack, Lpend) ;
	lsolve_numeric (n, Lp, Li, Lx, Pinv, nzb, Bi, Bx, top, Stack, X) ;

	/* ------------------------------------------------------------------ */
	/* partial pivoting */
	/* ------------------------------------------------------------------ */

	/* TODO: NaN case will fail */

	/* determine # nonzeros in L and U, and find the pivot */
	lnz = 0 ;
	pivot = 0 ;
	abs_pivot = 0 ;
	pivrow = EMPTY ;

	/* give the diagonal entry some preference */
	diag = 0 ;
	abs_diag = 0 ;
	found = FALSE ;

	for (p = top ; p < n ; p++)
	{
	    i = Stack [p] ;
	    if (!PIVOTAL (i))
	    {
		lnz++ ;
		x = X [i] ;
		abs_x = ABS (x) ;
		/* find the diagonal, if it is nonzero */
		if (!found && i == diagrow && abs_x > 0.)
		{
		    found = TRUE ;
		    diag = x ;
		    abs_diag = abs_x ;
		    PRINTF (("diagonal : %g\n", diag)) ;
		}
		/* find the partial-pivoting choice */
		if (abs_x > abs_pivot)
		{
		    pivrow = i ;
		    pivot = x ;
		    abs_pivot = abs_x ;
		}
	    }
	}

	PRINTF (("Partial pivot:  pivrow %d %g\n", pivrow, pivot)) ;
	if (found)
	{
	    /* compare the diagonal with the largest entry */
	    if (abs_diag >= tol * abs_pivot)
	    {
		/* the diagonal is large enough */
		pivrow = diagrow ;
		pivot = diag ;
		PRINTF (("diagonal pivot: pivrow %d %g\n", pivrow, pivot)) ;
	    }
	}
	ASSERT (!PIVOTAL (pivrow)) ;

	/* ------------------------------------------------------------------ */
	/* determine if LU factors have grown too big */
	/* ------------------------------------------------------------------ */

	unz = (n - top) - lnz + 1 ;

	PRINTF (("lp %d lnz %d lsize %d\n", lp, lnz, lsize)) ;
	if (lp + lnz > lsize)
	{
	    xsize = (GROWTH * ((double) lsize) + 2*n + 1) * sizeof (double) ;
	    if (INT_OVERFLOW (xsize))
	    {
		PRINTF (("Matrix is too large (L, int overflow)\n")) ;
		return (KLU_OUT_OF_MEMORY) ;
	    }
	    lsize = GROWTH * lsize + n + 1 ;
	    PRINTF (("inc L to %d\n", lsize)) ;
	    REALLOCATE (Li, int,    lsize, oki) ;
	    PRINTF (("did Li: %d\n", oki)) ;
	    REALLOCATE (Lx, double, lsize, okx) ;
	    PRINTF (("did Lx: %d\n", okx)) ;
	    *p_Li = Li ;
	    *p_Lx = Lx ;
	    if (!oki || !okx)
	    {
		PRINTF (("Matrix is too large (L)\n")) ;
		return (KLU_OUT_OF_MEMORY) ;
	    }
	    *p_lsize = lsize ;
	    PRINTF (("inc L to %d done\n", lsize)) ;
	}

	PRINTF (("up %d unz %d usize %d\n", up, unz, usize)) ;
	if (up + unz > usize)
	{
	    xsize = (GROWTH * ((double) usize) + 2*n + 1) * sizeof (double) ;
	    if (INT_OVERFLOW (xsize))
	    {
		PRINTF (("Matrix is too large (U, int overflow)\n")) ;
		return (KLU_OUT_OF_MEMORY) ;
	    }
	    usize = GROWTH * usize + n + 1 ;
	    PRINTF (("inc U to %d\n", usize)) ;
	    REALLOCATE (Ui, int,    usize, oki) ;
	    PRINTF (("did Ui: %d\n", oki)) ;
	    REALLOCATE (Ux, double, usize, okx) ;
	    PRINTF (("did Ux: %d\n", okx)) ;
	    *p_Ui = Ui ;
	    *p_Ux = Ux ;
	    if (!oki || !okx)
	    {
		PRINTF (("Matrix is too large (U)\n")) ;
		return (KLU_OUT_OF_MEMORY) ;
	    }
	    *p_usize = usize ;
	    PRINTF (("inc U to %d done\n", usize)) ;
	}

	/* ------------------------------------------------------------------ */
	/* singular case */
	/* ------------------------------------------------------------------ */

	PRINTF (("\nk %d : Pivot row %d : %g\n", k, pivrow, pivot)) ;
	if (pivrow == EMPTY || pivot == 0.0)
	{
	    PRINTF (("Matrix is singular\n")) ;
	    return (KLU_SINGULAR) ;
	}

	/* ------------------------------------------------------------------ */
	/* create the unit diagonal of L as first entry in the column of L */
	/* ------------------------------------------------------------------ */

	ASSERT (lp < lsize) ;
	Li [lp] = pivrow ;
	Lx [lp] = 1 ;
	lp++ ;

	/* ------------------------------------------------------------------ */
	/* copy into L and U, clear X, and scale the pivot column */
	/* ------------------------------------------------------------------ */

	for (p = top ; p < n ; p++)
	{
	    i = Stack [p] ;
	    /* NOTE: symmetric pruning and dropping zeros are not compatible */
	    if (PIVOTAL (i))
	    {
		ASSERT (up < usize) ;
		Ui [up] = Pinv [i] ;
		Ux [up] = X [i] ;
		up++ ;
	    }
	    else if (i != pivrow)
	    {
		ASSERT (lp < lsize) ;
		Li [lp] = i ;
		Lx [lp] = X [i] / pivot ;
		lp++ ;
	    }
	    X [i] = 0 ;
	}

	/* ------------------------------------------------------------------ */
	/* copy the pivot value as the last entry in the column of U */
	/* ------------------------------------------------------------------ */

	ASSERT (up < usize) ;
	Ui [up] = k ;
	Ux [up] = pivot ;
	up++ ;
	ASSERT (lp <= lsize) ;
	ASSERT (up <= usize) ;

	/* ------------------------------------------------------------------ */
	/* log the pivot permutation */
	/* ------------------------------------------------------------------ */

#ifndef NDEBUG
	kbar = UNFLIP (Pinv [diagrow]) ;
	ASSERT (kbar < n) ;
	ASSERT (P [kbar] == diagrow) ;
#endif
	if (pivrow != diagrow)
	{
	    /* an off-diagonal pivot has been chosen */
	    noffdiag++ ;
	    PRINTF ((">>>>>>>>>>>>>>>>> pivrow %d k %d off-diagonal\n",
			pivrow, k)) ;
	    if (!PIVOTAL (diagrow))
	    {
		/* the former diagonal row index, diagrow, has not yet been
		 * chosen as a pivot row.  Log this diagrow as the "diagonal"
		 * entry in the column kbar for which the chosen pivot row,
		 * pivrow, was originally logged as the "diagonal" */
		kbar = FLIP (Pinv [pivrow]) ;
		P [kbar] = diagrow ;
		Pinv [diagrow] = FLIP (kbar) ;
	    }
	}
	P [k] = pivrow ;
	Pinv [pivrow] = k ;

#ifndef NDEBUG
	for (i = 0 ; i < n ; i++) ASSERT (X [i] == 0) ;
	for (p = Up [k] ; p < up ; p++)
	{
	    PRINTF (("Column %d of U: %d : %g\n", k, Ui [p], Ux [p])) ;
	}
	for (p = Lp [k] ; p < lp ; p++)
	{
	    PRINTF (("Column %d of L: %d : %g\n", k, Li [p], Lx [p])) ;
	}
#endif

#ifdef PRUNE
	/* ------------------------------------------------------------------ */
	/* symmetric pruning */
	/* ------------------------------------------------------------------ */

	/* mark end of column k of L */
	/* this will change if and when column k of L is pruned. */
	Lpend [k] = lp ;

	/* check to see if any column of L can be pruned */
	for (p = Up [k] ; p < up-1 ; p++)   /* skip the diagonal entry */
	{
	    int j = Ui [p] ;
	    ASSERT (j < k) ;
	    if (!Lpruned [j])
	    {
		/* scan column j of L for the pivot row */
		int p2 ;
		int pend = Lpend [j] ;
		double x ;
		ASSERT (Lpend [j] == Lp [j+1]) ;
		for (p2 = Lp [j] ; p2 < pend ; p2++)
		{
		    if (pivrow == Li [p2])
		    {
			/* found it!  This column can be pruned */
			int phead, ptail ;

#ifndef NDEBUG
			PRINTF (("==== PRUNE: col j %d of L\n", j)) ;
			{
			    int p3 ;
			    for (p3 = Lp [j] ; p3 < Lp [j+1] ; p3++)
			    {
				PRINTF (("before: %i  pivotal: %d\n", Li [p3],
					    PIVOTAL (Li [p3]))) ;
			    }
			}
#endif

			/* flag j as pruned */
			Lpruned [j] = TRUE ;

			/* partition column j of L.  The unit diagonal of L
			 * remains as the first entry in the column of L. */
			phead = Lp [j] ;
			ptail = Lp [j+1] ;
			while (phead < ptail)
			{
			    i = Li [phead] ;
			    if (PIVOTAL (i))
			    {
				/* leave at the head */
				phead++ ;
			    }
			    else
			    {
				/* swap with the tail */
				ptail-- ;
				Li [phead] = Li [ptail] ;
				Li [ptail] = i ;
				x = Lx [phead] ;
				Lx [phead] = Lx [ptail] ;
				Lx [ptail] = x ;

			    }
			}
			Lpend [j] = ptail ;

#ifndef NDEBUG
			{
			    int p3 ;
			    for (p3 = Lp [j] ; p3 < Lp [j+1] ; p3++)
			    {
				if (p3 == Lpend [j]) PRINTF (("----\n")) ;
				PRINTF (("after: %i  pivotal: %d\n", Li [p3],
					    PIVOTAL (Li [p3]))) ;
			    }
			}
#endif

			break ;
		    }
		}
	    }
	}
#endif

    }

    /* ---------------------------------------------------------------------- */
    /* finalize column pointers for L and U, and put L in the pivotal order */
    /* ---------------------------------------------------------------------- */

    Lp [n] = lp ;
    Up [n] = up ;
    for (p = 0 ; p < lp ; p++)
    {
	ASSERT (PIVOTAL (Li [p])) ;
	Li [p] = Pinv [Li [p]] ;
    }

#ifndef NDEBUG
    for (i = 0 ; i < n ; i++)
    {
	PRINTF (("P[%d] = %d   Pinv[%d] = %d\n", i, P [i], i, Pinv [i])) ;
    }
    for (i = 0 ; i < n ; i++)
    {
	ASSERT (Pinv [i] >= 0 && Pinv [i] < n) ;
	ASSERT (P [i] >= 0 && P [i] < n) ;
	ASSERT (P [Pinv [i]] == i) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* shrink the LU factors to just the required size */
    /* ---------------------------------------------------------------------- */

    lsize = lp ;
    usize = up ;
    *p_lsize = lsize ;
    *p_usize = usize ;
    *p_noffdiag = noffdiag ;
    PRINTF (("noffdiag %d\n", noffdiag)) ;

    REALLOCATE (Li, int,    lsize, ok) ;
    REALLOCATE (Lx, double, lsize, ok) ;
    REALLOCATE (Ui, int,    usize, ok) ;
    REALLOCATE (Ux, double, usize, ok) ;
    *p_Li = Li ;
    *p_Lx = Lx ;
    *p_Ui = Ui ;
    *p_Ux = Ux ;

    return (KLU_OK) ;
}
