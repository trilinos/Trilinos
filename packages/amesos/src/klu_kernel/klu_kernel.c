/* ========================================================================== */
/* === klu kernel =========================================================== */
/* ========================================================================== */

/* Sparse left-looking LU factorization, with partial pivoting.  Based on
 * Gilbert & Peierl's method, with a non-recursive DFS and with Eisenstat &
 * Liu's symmetric pruning.  No user-callable routines are in this file.
 */

#include "klu.h"
#include "klu_kernel.h"

/* ========================================================================== */
/* === dfs ================================================================== */
/* ========================================================================== */

/* Does a depth-first-search, starting at node j. */

static void dfs
(
    /* input, not modified on output: */
    int j,		/* node at which to start the DFS */
    int k,		/* mark value, for the Flag array */
    int n,		/* L is n-by-n */
    int Lp [ ],		/* size n+1, column pointers of L */
    int Li [ ],		/* size lnz = Lp [n], row indices of L.  Diagonal (if
			 * present) is ignored.  L must be lower triangular. */
    int Pinv [ ],	/* Pinv [i] = k if row i is kth pivot row, or EMPTY if
			 * row i is not yet pivotal.  */

    /* workspace, not defined on input or output */
    int Stack [ ],	/* size n */

    /* input/output: */
    int Flag [ ],	/* Flag [i] == k means i is marked */
    int Lpend [ ],	/* for symmetric pruning */

    /* other, not defined on input or output */
    int adj_pos [ ]	/* keeps track of position in adj list during DFS */

    , int *plp
    , int Ui [ ]
    , int *pup
)
{
    int i, pos, jnew, head, pstart, lp, up ;

    lp = *plp ;
    up = *pup ;

    head = 0;
    Stack [0] = j ;
    ASSERT (Flag [j] != k) ;

    while (head >= 0)
    {
	j = Stack [head] ;
	jnew = Pinv [j] ;
	ASSERT (jnew >= 0 && jnew < k) ;	/* j is pivotal */

	if (Flag [j] != k)	    /* a node is not yet visited */
	{
	    /* first time that j has been visited */
	    Flag [j] = k;
	    PRINTF (("[ start dfs at %d : new %d\n", j, jnew)) ;
	    /* set adj_pos [head] to one past the last entry in col j to scan */
	    adj_pos [head] = (Lpend[jnew] == EMPTY) ? Lp[jnew+1] : Lpend[jnew] ;
	}

	/* add the adjacent nodes to the recursive stack by iterating through
	 * until finding another non-visited pivotal node */
	pstart = Lp [jnew] ;
	for (pos = --adj_pos [head] ; pos > pstart ; --pos)
	{
	    i = Li [pos] ;
	    if (Flag [i] != k)
	    {
		/* node i is not yet visited */
		if (Pinv [i] >= 0)
		{
		    /* keep track of where we left off in the scan of the
		     * adjacency list of node j so we can restart j where we
		     * left off. */
		    adj_pos [head] = pos ;

		    /* node i is pivotal; push it onto the recursive stack
		     * and immediately break so we can recurse on node i. */
		    Stack [++head] = i ;
		    break ;
		}
		else
		{
		    /* node i is not pivotal (no outgoing edges). */
		    /* Flag as visited and store directly into L,
		     * and continue with current node j. */
		    Flag [i] = k ;
		    Li [lp++] = i ;
		}
	    }
	}

	if (pos == pstart)
	{
	    /* if all adjacent nodes of j are already visited, pop j from
	     * recursive stack and push j onto output stack */
	    head-- ;
	    Ui [up++] = j ;
	    PRINTF (("  end   dfs at %d ] head : %d\n", j, head)) ;
	}
    }

    *plp = lp ;
    *pup = up ;
}

/* ========================================================================== */
/* === lsolve_symbolic ====================================================== */
/* ========================================================================== */

/* Finds the pattern of x, for the solution of Lx=b */

static void lsolve_symbolic
(
    /* input, not modified on output: */
    int k,		/* mark value, for the Flag array */
    int n,		/* L is n-by-n, where n >= 0 */
    int Lp [ ],		/* of size n+1, column pointers of L */
    int Li [ ],		/* size lnz=Lp [n], row indices of L */
    int Pinv [ ],	/* Pinv [i] = k if i is kth pivot row, or EMPTY if row i
			 * is not yet pivotal.  */
    int nzb,		/* nz in b */
    int Bi [ ],		/* pattern of b, of size nzb */

    /* workspace, not defined on input or output */
    int Stack [ ],	/* size n */

    /* workspace, defined on input and output */
    int Flag [ ],	/* size n.  Initially, all of Flag [0..n-1] < k.  After
			 * lsolve_symbolic is done, Flag [i] == k if i is in
			 * the pattern of the output, and Flag [0..n-1] <= k. */

    /* other */
    int Lpend [ ],	/* for symmetric pruning */
    int adj_pos [ ],	/* workspace used in dfs */

    int *plp,
    int Ui [ ],
    int *pup
)
{
    int i, p ;

    /* compute the pattern of X */
    for (p = 0 ; p < nzb ; p++)
    {
	i = Bi [p] ;
	PRINTF (("\n ===== DFS at node %d in b, inew: %d\n", i, Pinv [i])) ;
	if (Flag [i] != k)
	{
	    if (Pinv [i] >= 0)
	    {
		dfs (i, k, n, Lp, Li, Pinv, Stack, Flag,
		       Lpend, adj_pos, plp, Ui, pup) ;
	    }
	    else
	    {
		/* i is not pivotal, and not flagged. Flag it and put in L */
		Flag [i] = k ;
		Li [(*plp)++] = i ;
	    }
	}
    }
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
    int Li [ ],		/* size lnz=Lp [n], row indices of L */
    double Lx [ ],	/* size lnz=Lp [n], values of L  */
    int Pinv [ ],	/* Pinv [i] = k if i is kth pivot row, or EMPTY if row i
			 * is not yet pivotal.  */
    int nzb,		/* nz in b */
    int Bi [ ],		/* pattern of b, of size nzb */
    double Bx [ ],	/* values of b, of size nzb */

    /* output, must be zero on input: */
    double X [ ]	/* size n, initially zero.  On output,
			 * X [Ui [up1..up-1]] and X [Li [lp1..lp-1]]
			 * contains the solution. */

    , int Ui [ ]
    , int up1
    , int up
)
{
    double xj ;
    int p, s, j, jnew, pend ;

    /* scatter b into x */
    for (p = 0 ; p < nzb ; p++)
    {
	X [Bi [p]] = Bx [p] ;
    }
    /* solve Lx=b */
    for (s = up-1 ; s >= up1 ; s--)
    {
	/* forward solve with column j of L (skip the unit diagonal of L) */
	j = Ui [s] ;
	jnew = Pinv [j] ;
	ASSERT (jnew >= 0) ;
	xj = X [j] ;
	ASSERT (Lp [jnew] < Lp [jnew+1] && Li [Lp [jnew]] == j) ;
	pend = Lp [jnew+1] ;
	for (p = Lp [jnew] + 1 ; p < pend ; p++)
	{
	    X [Li [p]] -= Lx [p] * xj ;
	}
    }
}

/* ========================================================================== */
/* === lpivot =============================================================== */
/* ========================================================================== */

/* Find a pivot via partial pivoting, and scale the column of L. */

static void lpivot
(
    int lp,
    int lp1,
    int diagrow,
    int Li [ ],
    double Lx [ ],
    int *p_pivrow,
    double *p_pivot,
    double *p_abs_pivot,
    double tol,
    double X [ ]
)
{
    double x, pivot, abs_pivot ;
    int p, i, ppivrow, pdiag, pivrow ;

    /* look in Li [lp1 .. lp-1] for a pivot row */

    /* gather the first entry from X and store in L */
    p = lp1 ;
    i = Li [p] ;
    x = X [i] ;
    X [i] = 0 ;
    Lx [p] = x ;
    x = fabs (x) ;

    /* find the diagonal */
    pdiag = (i == diagrow) ? p : EMPTY ;


    /* find the partial-pivoting choice */
    abs_pivot = x ;
    ppivrow = p ;

    for (p = lp1 + 1 ; p < lp ; p++)
    {
	/* gather the entry from X and store in L */
	i = Li [p] ;
	x = X [i] ;
	X [i] = 0 ;
	Lx [p] = x ;
	x = fabs (x) ;

	/* find the diagonal */
	if (i == diagrow)
	{
	    pdiag = p ;
	}

	/* find the partial-pivoting choice */
	if (x > abs_pivot)
	{
	    abs_pivot = x ;
	    ppivrow = p ;
	}
    }

    /* compare the diagonal with the largest entry */
    if (pdiag != EMPTY)
    {
	x = fabs (Lx [pdiag]) ;
	if (x >= tol * abs_pivot)
	{
	    /* the diagonal is large enough */
	    abs_pivot = x ;
	    ppivrow = pdiag ;
	}
    }

    /* swap pivot row to first entry in column k of L */
    pivrow = Li [ppivrow] ;
    pivot  = Lx [ppivrow] ;
    Li [ppivrow] = Li [lp1] ;
    Lx [ppivrow] = Lx [lp1] ;
    Li [lp1] = pivrow ;
    Lx [lp1] = 1.0 ;

    /* divide L by the pivot value */
    for (p = lp1 + 1 ; p < lp ; p++)
    {
	Lx [p] /= pivot ;
    }

    *p_pivrow = pivrow ;
    *p_pivot = pivot ;
    *p_abs_pivot = abs_pivot ;
}


/* ========================================================================== */
/* === prune ================================================================ */
/* ========================================================================== */

/* Prune the columns of L to reduce work in subsequent depth-first searches */

static void prune
(
    int Up [ ],
    int Ui [ ],
    int Lp [ ],
    int Li [ ],
    double Lx [ ],
    int Lpend [ ],
    int Pinv [ ],
    int up,
    int k,
    int pivrow
)
{
    double x ;
    int p, i, j, p2, pend, phead, ptail ;

    /* check to see if any column of L can be pruned */
    for (p = Up [k] ; p < up-1 ; p++)   /* skip the diagonal entry */
    {
	int j = Ui [p] ;
	ASSERT (j < k) ;
	PRINTF (("%d is pruned: %d. Lpend[j] %d Lp[j+1] %d\n",
	    j, Lpend [j] != EMPTY, Lpend [j], Lp [j+1])) ;
	if (Lpend [j] == EMPTY)
	{
	    /* scan column j of L for the pivot row */
	    pend = Lp [j+1] ;
	    for (p2 = Lp [j] ; p2 < pend ; p2++)
	    {
		if (pivrow == Li [p2])
		{
		    /* found it!  This column can be pruned */

#ifndef NDEBUG
		    PRINTF (("==== PRUNE: col j %d of L\n", j)) ;
		    {
			int p3 ;
			for (p3 = Lp [j] ; p3 < Lp [j+1] ; p3++)
			{
			    PRINTF (("before: %i  pivotal: %d\n", Li [p3],
					Pinv [Li [p3]] >= 0)) ;
			}
		    }
#endif

		    /* partition column j of L.  The unit diagonal of L
		     * remains as the first entry in the column of L. */
		    phead = Lp [j] ;
		    ptail = Lp [j+1] ;
		    while (phead < ptail)
		    {
			i = Li [phead] ;
			if (Pinv [i] >= 0)
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

		    /* set Lpend to one past the last entry in the
		     * first part of the column of L.  Entries in
		     * Li [Lp [j] ... Lpend [j]-1] are the only part of
		     * column j of L that needs to be scanned in the DFS.
		     * Lpend [j] was EMPTY; setting it >= 0 also flags
		     * column j as pruned. */
		    Lpend [j] = ptail ;

#ifndef NDEBUG
		    {
			int p3 ;
			for (p3 = Lp [j] ; p3 < Lp [j+1] ; p3++)
			{
			    if (p3 == Lpend [j]) PRINTF (("----\n")) ;
			    PRINTF (("after: %i  pivotal: %d\n", Li [p3],
					Pinv [Li [p3]] >= 0)) ;
			}
		    }
#endif

		    break ;
		}
	    }
	}
    }
}


/* ========================================================================== */
/* === klu_kernel =========================================================== */
/* ========================================================================== */

int klu_kernel	    /* returns KLU_OK (0) or KLU_OUT_OF_MEMORY (-2) */
(
    /* input, not modified */
    int n,	    /* A is n-by-n */
    int Ap [ ],	    /* size n+1, column pointers for A */
    int Ai [ ],	    /* size nz = Ap [n], row indices for A */
    double Ax [ ],  /* size nz, values of A */
    int Q [ ],	    /* size n, optional input permutation */
    double tol,	    /* partial pivoting tolerance parameter */
    double growth,  /* memory growth factor */
    int lsize,	    /* initial size of Li and Lx, final size is Lp[n] */
    int usize,	    /* initial size of Ui and Ux, final size is Up[n] */

    /* output, not defined on input */
    int Lp [ ],	    /* size n+1, column pointers for L */
    int **p_Li,	    /* size lsize on input, Lp[n] on output. row indices of L */
    double **p_Lx,  /* size lsize on input, Lp[n] on output. values of L */
    int Up [ ],	    /* size n+1, column pointers for U */
    int **p_Ui,	    /* size usize on input, Up[n] on input. row indices of U */
    double **p_Ux,  /* size usize on input, Up[n] on input. values of U */
    int Pinv [ ],   /* size n, inverse row permutation, where Pinv [i] = k if
		     * row i is the kth pivot row */
    int P [ ],	    /* size n, row permutation, where P [k] = i if row i is the
		     * kth pivot row. */
    int *p_noffdiag,	/* # of off-diagonal pivots chosen */
    double *p_umin, /* smallest entry on the diagonal of U */
    double *p_umax, /* largest entry on the diagonal of U */

    /* workspace, not defined on input */
    double X [ ],   /* size n, undefined on input, zero on output */

    /* workspace, not defined on input or output */
    int Stack [ ],  /* size n */
    int Flag [ ],   /* size n */
    int adj_pos [ ],	/* size n */

    /* other workspace: */
    int Lpend [ ]		    /* size n workspace, for pruning only */
)
{ 
    double pivot, *Bx, *Lx, *Ux, abs_pivot, xsize, umin, umax ;
    int lp, up, k, nzb, *Bi, p, i, pivrow, kbar,
	*Li, *Ui, ok, oki, okx, kcol, diagrow, noffdiag, firstrow,
	jnew, lp1, up1 ;

    /* ---------------------------------------------------------------------- */
    /* get initial Li, Lx, Ui, and Ux */
    /* ---------------------------------------------------------------------- */

    PRINTF (("input: lsize %d usize %d\n", lsize, usize)) ;

    if (lsize <= 0 || usize <= 0)
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

    umin = 0 ;
    umax = 0 ;
    firstrow = 0 ;
    noffdiag = 0 ;
    lp = 0 ;
    up = 0 ;
    for (k = 0 ; k < n ; k++)
    {
	X [k] = 0 ;
	Flag [k] = EMPTY ;
	Lpend [k] = EMPTY ;	/* flag k as not pruned */
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

	/* ------------------------------------------------------------------ */
	/* determine if LU factors have grown too big */
	/* ------------------------------------------------------------------ */

	PRINTF (("lp %d lsize %d\n", lp, lsize)) ;
	if (lp + (n-k) > lsize)
	{
	    xsize = (growth * ((double) lsize) + 2*n + 1) * sizeof (double) ;
	    if (INT_OVERFLOW (xsize))
	    {
		PRINTF (("Matrix is too large (L, int overflow)\n")) ;
		return (KLU_OUT_OF_MEMORY) ;
	    }
	    lsize = growth * lsize + n + 1 ;
	    REALLOCATE (Li, int,    lsize, oki) ;
	    REALLOCATE (Lx, double, lsize, okx) ;
	    *p_Li = Li ;
	    *p_Lx = Lx ;
	    if (!oki || !okx)
	    {
		PRINTF (("Matrix is too large (L)\n")) ;
		return (KLU_OUT_OF_MEMORY) ;
	    }
	    PRINTF (("inc L to %d done\n", lsize)) ;
	}

	PRINTF (("up %d usize %d\n", up, usize)) ;
	if (up + k > usize)
	{
	    xsize = (growth * ((double) usize) + 2*n + 1) * sizeof (double) ;
	    if (INT_OVERFLOW (xsize))
	    {
		PRINTF (("Matrix is too large (U, int overflow)\n")) ;
		return (KLU_OUT_OF_MEMORY) ;
	    }
	    usize = growth * usize + n + 1 ;
	    REALLOCATE (Ui, int,    usize, oki) ;
	    REALLOCATE (Ux, double, usize, okx) ;
	    *p_Ui = Ui ;
	    *p_Ux = Ux ;
	    if (!oki || !okx)
	    {
		PRINTF (("Matrix is too large (U)\n")) ;
		return (KLU_OUT_OF_MEMORY) ;
	    }
	    PRINTF (("inc U to %d done\n", usize)) ;
	}

	/* ------------------------------------------------------------------ */
	/* start the kth column of L and U */
	/* ------------------------------------------------------------------ */

	PRINTF (("\n\n==================================== k: %d\n", k)) ;
	Lp [k] = lp ;
	Up [k] = up ;
	lp1 = lp ;
	up1 = up ;

	/* ------------------------------------------------------------------ */
	/* get the kth column of A and find the "diagonal" */
	/* ------------------------------------------------------------------ */

	kcol = (Q == (int *) NULL) ? (k) : (Q [k]) ;
	nzb = Ap [kcol+1] - Ap [kcol] ;
	Bi = Ai + Ap [kcol] ;
	Bx = Ax + Ap [kcol] ;

	diagrow = P [k] ;   /* might already be pivotal */
	PRINTF (("k %d kcol %d, diagrow = %d, UNFLIP(diagrow) = %d\n",
	    k, kcol, diagrow, UNFLIP (diagrow))) ;

	/* ------------------------------------------------------------------ */
	/* compute the nonzero pattern of the kth column of L and U */
	/* ------------------------------------------------------------------ */

#ifndef NDEBUG
	for (i = 0 ; i < n ; i++)
	{
	    ASSERT (Flag [i] < k) ;
	    ASSERT (X [i] == 0) ;
	}
#endif

	lsolve_symbolic (k, n, Lp, Li, Pinv, nzb, Bi, Stack, Flag,
		               Lpend, adj_pos, &lp, Ui, &up) ;

#ifndef NDEBUG
	PRINTF (("--- in U:\n")) ;
	for (p = up-1 ; p >= up1 ; p--)
	{
	    PRINTF (("pattern of X for U: %d : %d pivot row: %d\n",
		p, Ui [p], Pinv [Ui [p]])) ;
	    ASSERT (Flag [Ui [p]] == k) ;
	}
	PRINTF (("--- in L:\n")) ;
	for (p = lp1 ; p < lp ; p++)
	{
	    PRINTF (("pattern of X in L: %d : %d pivot row: %d\n",
		p, Li [p], Pinv [Li [p]])) ;
	    ASSERT (Flag [Li [p]] == k) ;
	}
	p = 0 ;
	for (i = 0 ; i < n ; i++)
	{
	    ASSERT (Flag [i] <= k) ;
	    if (Flag [i] == k) p++ ;
	}
	ASSERT (p == (up-up1) + (lp-lp1)) ;
#endif

	/* ------------------------------------------------------------------ */
	/* compute the numerical values of the kth column (s = L \ A (:,k)) */
	/* ------------------------------------------------------------------ */

	lsolve_numeric (n, Lp, Li, Lx, Pinv, nzb, Bi, Bx, X, Ui, up1, up) ; 

#ifndef NDEBUG
	for (p = up-1 ; p >= up1 ; p--)
	{
	    PRINTF (("X for U %d : %g\n", Ui [p], X [Ui [p]])) ;
	}
	for (p = lp1 ; p < lp ; p++)
	{
	    PRINTF (("X for L %d : %g\n", Li [p], X [Li [p]])) ;
	}
#endif

	/* ------------------------------------------------------------------ */
	/* partial pivoting with diagonal preference */
	/* ------------------------------------------------------------------ */

	if (lp == lp1)
	{
	    /* matrix is structurally singular */
	    PRINTF (("Matrix is singular\n")) ;
	    for ( ; firstrow < n ; firstrow++)
	    {
		PRINTF (("check %d\n", firstrow)) ;
		if (Pinv [firstrow] < 0)
		{
		    /* found the lowest-numbered non-pivotal row.  Pick it. */
		    pivrow = firstrow ;
		    PRINTF (("Got pivotal row: %d\n", pivrow)) ;
		    break ;
		}
	    }
	    Li [lp1] = pivrow ;
	    Lx [lp1] = 1.0 ;
	    lp++ ;
	    pivot = 0 ;
	    abs_pivot = 0 ;
	}
	else
	{
	    /* find a pivot and scale the pivot column */
	    lpivot (lp, lp1, diagrow, Li, Lx, &pivrow, &pivot, &abs_pivot,
		   tol, X) ;
	}

	/* we now have a valid pivot row, even if the column has NaN's or
	 * has no entries on or below the diagonal at all. */
	PRINTF (("\nk %d : Pivot row %d : %g\n", k, pivrow, pivot)) ;
	ASSERT (pivrow >= 0 && pivrow < n) ;
	ASSERT (Pinv [pivrow] < 0) ;

	/* ------------------------------------------------------------------ */
	/* keep track of the smallest and largest entry on the diagonal of U */
	/* ------------------------------------------------------------------ */

	if (k == 0)
	{
	    umin = abs_pivot ;
	    umax = abs_pivot ;
	}
	else if (abs_pivot < umin)
	{
	    umin = abs_pivot ;
	}
	else if (abs_pivot > umax)
	{
	    umax = abs_pivot ;
	}

	/* ------------------------------------------------------------------ */
	/* finalize the kth column of U and clear X */
	/* ------------------------------------------------------------------ */

	for (p = up-1 ; p >= up1 ; p--)
	{
	    i = Ui [p] ;
	    Ui [p] = Pinv [i] ;
	    Ux [p] = X [i] ;
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
	    if (Pinv [diagrow] < 0)
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
	for (p = lp1 ; p < lp ; p++)
	{
	    PRINTF (("Column %d of L: %d : %g\n", k, Li [p], Lx [p])) ;
	}
#endif

	/* ------------------------------------------------------------------ */
	/* symmetric pruning */
	/* ------------------------------------------------------------------ */

	prune (Up, Ui, Lp, Li, Lx, Lpend, Pinv, up, k, pivrow) ;

    }

    /* ---------------------------------------------------------------------- */
    /* finalize column pointers for L and U, and put L in the pivotal order */
    /* ---------------------------------------------------------------------- */

    Lp [n] = lp ;
    Up [n] = up ;
    for (p = 0 ; p < lp ; p++)
    {
	Li [p] = Pinv [Li [p]] ;
    }

#ifndef NDEBUG
    for (i = 0 ; i < n ; i++)
    {
	PRINTF (("P [%d] = %d   Pinv [%d] = %d\n", i, P [i], i, Pinv [i])) ;
    }
    for (i = 0 ; i < n ; i++)
    {
	ASSERT (Pinv [i] >= 0 && Pinv [i] < n) ;
	ASSERT (P [i] >= 0 && P [i] < n) ;
	ASSERT (P [Pinv [i]] == i) ;
	ASSERT (X [i] == 0) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* shrink the LU factors to just the required size */
    /* ---------------------------------------------------------------------- */

    lsize = lp ;
    usize = up ;
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
