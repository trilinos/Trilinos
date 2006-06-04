/* ========================================================================== */
/* === paraklete_kernel ===================================================== */
/* ========================================================================== */

#include "paraklete.h"

/* Sparse left-looking LU factorization, with partial pivoting.  Based on
 * Gilbert & Peierl's method, with a non-recursive DFS and with Eisenstat &
 * Liu's symmetric pruning.  No user-callable routines are in this file.
 *
 * PARAKLETE version 0.1: parallel sparse LU factorization.  May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

#ifndef NDEBUG
/* global variables for debugging only */
static int debug_k, debug_nfound, debug_n, debug_npiv ;
#endif

/* ========================================================================== */
/* === dfs ================================================================== */
/* ========================================================================== */

/* Does a depth-first-search, starting at node j. */

static int dfs		/* returns new top of output stack */
(
    /* input, not modified on output: */
    int npiv,
    int j,		/* node at which to start the DFS */
    int mark,		/* mark value for Flag */
    int Pinv [ ],	/* Pinv [i] = k if row i is kth pivot row, or EMPTY if
			 * row i is not yet pivotal.  Only used in phase 1. */
    int Llen [ ],
    int Lip [ ],

    int phase1,		/* TRUE if in phase 1, FALSE if in phase 2 */

    /* workspace, not defined on input or output */
    int Stack [ ],	/* size n */

    /* input/output: */
    int Flag [ ],	/* Flag [i] >= mark means i is marked */
    int Lprune [ ],	/* for symmetric pruning */
    int top,		/* top of stack on input */
    double LU [ ],
    int Li [ ],		/* resulting column of L */
    int *plength,

    /* other, not defined on input or output */
    int Pstack [ ]	/* keeps track of position in adj list during DFS */
)
{
    int i, pos, jnew, head, llen ;
    int *Lj ;

    llen = *plength ;
    head = 0 ;
    Stack [0] = j ;
    ASSERT (Flag [j] < mark) ;

    while (head >= 0)
    {
	/* In phase 1, j is in the original row index space of A.  In
	 * phase 2, j is in the final pivotal order.  In both phases, jnew is
	 * in the pivotal row order. */
	j = Stack [head] ;
	jnew = (phase1 && j < npiv) ? (Pinv [j]) : (j) ;

	ASSERT (IMPLIES ( phase1, jnew >= 0 && jnew < debug_k)) ;
	ASSERT (IMPLIES (!phase1, jnew >= 0 && jnew < debug_nfound)) ;

	if (Flag [j] < mark)	    /* a node is not yet visited */
	{
	    /* first time that j has been visited */
	    Flag [j] = mark ;
	    /* set Pstack [head] to one past the last entry in col j to scan */
	    Pstack [head] = (Lprune [jnew] == EMPTY) ?
		(Llen [jnew]) : (Lprune [jnew]) ;
	}

	/* add the adjacent nodes to the recursive stack by iterating through
	 * until finding another non-visited pivotal node */
	Lj = (int *) (LU + Lip [jnew]) ;
	for (pos = --Pstack [head] ; pos > 0 ; --pos)
	{
	    /* In phase 1, i is in the original row index space of A.  In
	     * phase 2, i is in the final pivotal order. */
	    i = Lj [pos] ;

	    ASSERT (IMPLIES ( phase1, i >= 0 && i < debug_n)) ;
	    ASSERT (IMPLIES (!phase1, i >= 0 && i < debug_nfound)) ;

	    if (Flag [i] < mark)
	    {
		/* node i is not yet visited */
		if (i < npiv && (!phase1 || Pinv [i] >= 0))
		{
		    ASSERT (i < npiv) ;
		    /* keep track of where we left off in the scan of the
		     * adjacency list of node j so we can restart j where we
		     * left off. */
		    Pstack [head] = pos ;

		    /* node i is pivotal; push it onto the recursive stack
		     * and immediately break so we can recurse on node i. */
		    Stack [++head] = i ;
		    break ;
		}
		else
		{
		    /* node i is not pivotal (no outgoing edges).
		     * Flag as visited and store directly into L,
		     * and continue with current node j.
		     * This cannot occur during phase 2. */
		    Flag [i] = mark ;
		    Li [llen++] = i ;
		}
	    }
	}

	if (pos == 0)
	{
	    /* if all adjacent nodes of j are already visited, pop j from
	     * recursive stack and push j onto output stack */
	    head-- ;
	    Stack [--top] = j ;
	}
    }

    *plength = llen ;
    return (top) ;	    /* return top of output stack */
}


/* ========================================================================== */
/* === lsolve_symbolic ====================================================== */
/* ========================================================================== */

/* Finds the pattern of x, for the solution of Lx=b */

static void lsolve_symbolic
(
    /* input, not modified on output: */
    int n,              /* L is n-by-n, where n >= 0 */
    int k,		/* kth step of factorization */
    int mark,		/* mark for Flag */
    int kcol,		/* b = A (:,kcol) */
    int Ap [ ],
    int Ai [ ],
    int Pinv [ ],	/* Pinv [i] = k if i is kth pivot row, or EMPTY if row i
			 * is not yet pivotal.  In phase 2, all rows are in
			 * their final order, and Pinv is a complete inverse
			 * permutation. */

    int phase1,		/* TRUE if in phase 1, FALSE if in phase 2 */
    int nfound,		/* used in phase 2 only */
    int npiv,

    /* workspace, not defined on input or output */
    int Stack [ ],	/* size n */

    /* workspace, defined on input and output */
    int Flag [ ],	/* size n.  Initially, all of Flag [0..n-1] < mark.
			 * After lsolve_symbolic is done, Flag [i] == mark if i
			 * is in the pattern of the output, and
			 * Flag [0..n-1] <= mark. */

    /* other */
    int Lprune [ ],	/* for symmetric pruning */
    int Pstack [ ],	/* workspace used in dfs */

    /* TODO: comment this */
    double LU [ ],
    int *plu,
    int Llen [ ],
    int Ulen [ ],
    int Lip [ ],
    int Uip [ ]
)
{
    int i, p, p2, pend, top, llen, ulen, lup ;
    int *Li, *Ui ;

    top = n ;
    llen = 0 ;
    lup = *plu ;

    if (phase1)
    {
	Lip [k] = lup ;
    }
    Li = (int *) (LU + lup) ;
    pend = Ap [kcol+1] ;

    for (p = Ap [kcol] ; p < pend ; p++)
    {
	/* Ai is always in original row order, since it is never modified.
	 * In phase 1, we need to leave i in its original row index space.
	 * In phase 2, i needs to refer to the ith pivot row instead. */
	if (phase1)
	{
	    i = Ai [p] ;
	}
	else
	{
	    i = Ai [p] ;
	    i = (i < npiv) ? Pinv [i] : i ;
	    if (i >= nfound) continue ;
	}

	/* (i,k) is an entry in the block.  start a DFS at node i */
	if (Flag [i] < mark)
	{
	    if (i < npiv && (!phase1 || Pinv [i] >= 0))
	    {
		top = dfs (npiv, i, mark, Pinv, Llen, Lip, phase1, Stack, Flag,
			   Lprune, top, LU, Li, &llen, Pstack) ;
	    }
	    else
	    {
		/* i is not pivotal, and not flagged. Flag and put in L.
		 * This cannot occur during phase 2. */
		Flag [i] = mark ;
		Li [llen++] = i ;
	    }
	}
    }

    /* for phase 2, llen will be zero */
    if (phase1)
    {
	Llen [k] = llen ;
    }
    ASSERT (IMPLIES (!phase1, llen == 0)) ;

    /* advance LU pointer past the pattern and values of the new column of L */ 
    lup += ((llen + 1) / 2) + llen ;

    /* find the length of the column of U */
    ulen = n - top ;
    if (phase1)
    {
	/* add one for the pivot itself */
	ulen++ ;
    }
    Ulen [k] = ulen ;
    Ui = (int *) (LU + lup) ;
    Uip [k] = lup ;

    /* advance LU pointer past the pattern and values of the new column of U */ 
    lup += ((ulen + 1) / 2) + ulen ;

    /* extract Stack [top..n-1] to Ui */
    for (p = top, p2 = 0 ; p < n ; p++, p2++)
    {
	Ui [p2] = Stack [p] ;
	ASSERT (IMPLIES (!phase1, Ui [p2] < debug_nfound)) ;
    }

    /* position the lu index at the starting point for next column */
    *plu = lup ;
}


/* ========================================================================== */
/* === construct_column ===================================================== */
/* ========================================================================== */

/* Scatter the column of A into the workspace X */

static void construct_column
(
    /* inputs, not modified on output */
    int kcol,	    /* the column of A to construct */
    int Ap [ ],
    int Ai [ ],
    double Ax [ ],
    int phase1,	    /* phase 1: computing L1, L2 and U1.  phase 2: just U2 */
    int nfound,	    /* only used in phase 2 */
    int npiv,
    int Pinv [ ],   /* only used in phase 2 */

    /* zero on input, modified on output */
    double X [ ]
)
{
    int p, pend, i ;

    pend = Ap [kcol+1] ;
    if (phase1)
    {
	/* scatter the entire column of A */
	for (p = Ap [kcol] ; p < pend ; p++)
	{
	    X [Ai [p]] = Ax [p] ;
	}
    }
    else
    {
	/* scatter only the pivotal rows of A (for the U2 part only),
	 * and permute to their final pivot order as well. */
	for (p = Ap [kcol] ; p < pend ; p++)
	{
	    i = Ai [p] ;
	    i = (i < npiv) ? Pinv [i] : i ;
	    if (i < nfound)
	    {
		X [i] = Ax [p] ;
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
    int npiv,
    int Pinv [ ],	/* Pinv [i] = k if i is kth pivot row, or EMPTY if row i
			 * is not yet pivotal.  */

    /* TODO comment this */
    double *LU,
    int Uip [ ],
    int Lip [ ],
    int Ulen [ ],
    int Llen [ ],
    int k,

    int phase1,
    int Lphase_len [ ],	/* phase 1: Llen, phase 2: L1_len */

    /* output, must be zero on input: */
    double X [ ]	/* size n, initially zero.  On output,
			 * X [Ui [up1..up-1]] and X [Li [lp1..lp-1]]
			 * contains the solution. */
)
{
    double xj ;
    double *Lx ;
    int p, s, j, jnew, llen, ulen ;
    int *Ui, *Li ;  

    Ui = (int *) (LU + Uip [k]) ;
    ulen = Ulen [k] ;

    if (phase1) ulen-- ;    /* skip the diagonal */

    /* solve Lx=b */
    for (s = 0 ; s < ulen ; s++)
    {
	/* forward solve with column j of L (skip the unit diagonal of L).
	 * In phase 1, all of L1 and L2 is used.  Ui not yet in pivotal order.
	 * In phase 2, only L1 is used.  Ui already in pivotal order. */
 	j = Ui [s] ;
	jnew = (phase1 && j < npiv) ? (Pinv [j]) : j ;

	ASSERT (IMPLIES ( phase1, jnew >= 0 && jnew < debug_n)) ;
	ASSERT (IMPLIES (!phase1, jnew >= 0 && jnew < debug_nfound)) ;

	xj = X [j] ;
	GET_COLUMN (Lip, Llen, LU, jnew, Li, Lx, llen) ;

	/* phase 1: solve with entire column of L (L1 and L2).
	 * phase 2: solve with L1 only */
	llen = Lphase_len [jnew] ;

	for (p = 1 ; p < llen ; p++)
	{
	    ASSERT (IMPLIES (!phase1, Li [p] > j && Li [p] < debug_nfound)) ;
	    X [Li [p]] -= Lx [p] * xj ;
	}
    }
}


/* ========================================================================== */
/* === lsolve =============================================================== */
/* ========================================================================== */

/* Sparse symbolic and numeric solve of Lx=b to compute the kth column of L and
 * U.  In phase 1: L1, L2, and U1 are computed.  In phase 2: only U2 is
 * computed. */

static void lsolve
(
   int phase1,
   int nfound,
   int npiv,
   int n,
   int k,
   int kcol,
   int Ap [ ],
   int Ai [ ],
   double Ax [ ],

   double *LU,
   int *lup,
   int Lip [ ],
   int Uip [ ],
   int Llen [ ],
   int Ulen [ ],

   int Lphase_len [ ],

   int Pinv [ ],
   int Stack [ ],
   int Lprune [ ],
   int Pstack [ ],
   int Flag [ ],
   int mark,
   double X [ ]
)
{
    double *Ux ;
    int i, p, ulen ;
    int *Ui ;

    /* ---------------------------------------------------------------------- */
    /* compute the nonzero pattern of the kth column of L and U */
    /* ---------------------------------------------------------------------- */

    lsolve_symbolic (n, k, mark, kcol, Ap, Ai, Pinv, phase1, nfound, npiv,
	    Stack, Flag, Lprune, Pstack, LU, lup, Llen, Ulen, Lip, Uip) ;

    /* ---------------------------------------------------------------------- */
    /* get the column of the matrix to factorize and scatter into X */
    /* ---------------------------------------------------------------------- */

    construct_column (kcol, Ap, Ai, Ax, phase1, nfound, npiv, Pinv, X) ;

    /* ---------------------------------------------------------------------- */
    /* compute the numerical values of the kth column (s = L \ A(:,kcol)) */
    /* ---------------------------------------------------------------------- */

    lsolve_numeric (npiv,
	    Pinv, LU, Uip, Lip, Ulen, Llen, k, phase1, Lphase_len, X) ;

    /* --------------------------------------------------------------------- */
    /* finalize the kth column of U (pivot excluded) and clear X */
    /* --------------------------------------------------------------------- */

    GET_COLUMN (Uip, Ulen, LU, k, Ui, Ux, ulen) ;

    if (phase1)
    {
	/* phase 1: modify row indices in Ui to be in their final pivot order */
	ulen-- ;	    /* skip the diagonal */
	for (p = 0 ; p < ulen ; p++)
	{
	    i = Ui [p] ;
	    ASSERT (i >= 0 && i < npiv) ;
	    ASSERT (Pinv [i] >= 0 && Pinv [i] < npiv) ;
	    Ui [p] = Pinv [i] ;
	    Ux [p] = X [i] ;
	    X [i] = 0 ;
	}
    }
    else
    {
	/* phase 2: Ui is already in pivotal order */
	for (p = 0 ; p < ulen ; p++)
	{
	    i = Ui [p] ;
	    ASSERT (i >= 0 && i < nfound) ;
	    Ux [p] = X [i] ;
	    X [i] = 0 ;
	}
    }
}


/* ========================================================================== */
/* === lpivot =============================================================== */
/* ========================================================================== */

/* Find a pivot via threshold partial pivoting, and scale the column of L. 
 * This routine is active during phase 1 only.  Returns TRUE if pivot found,
 * FALSE otherwise. */

static int lpivot
(
    int diagrow,
    int *p_pivrow,
    double *p_pivot,
    double *p_abs_pivot,
    double tol_diag,
    double tol_offdiag,
    double X [ ],
    double *LU,
    int Lip [ ],
    int Llen [ ],
    int k,
    int npiv
)
{
    double x, pivot, abs_pivot, max_entry ;
    double *Lx ;
    int p, i, ppivrow, pdiag, pivrow, pivot_found, llen ;
    int *Li ;

    /* ---------------------------------------------------------------------- */
    /* look in Li [0 ... Llen [k] - 1 ] for a pivot row */
    /* ---------------------------------------------------------------------- */

    GET_COLUMN (Lip, Llen, LU, k, Li, Lx, llen) ;

    pdiag = EMPTY ;
    abs_pivot = -1 ;
    max_entry = -1 ;
    ppivrow = EMPTY ;

    for (p = 0 ; p < llen ; p++)
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

	/* find the partial-pivoting choice (constrained to rows 0..npiv-1) */
	if (x > abs_pivot && i < npiv)
	{
	    abs_pivot = x ;
	    ppivrow = p ;
	}

	/* find the max entry in this column (any row) */
	max_entry = MAX (max_entry, x) ;
    }

    /* compare the diagonal with the largest entry */
    if (pdiag != EMPTY)
    {
	x = fabs (Lx [pdiag]) ;
	if (x >= tol_diag * max_entry)
	{
	    /* the diagonal is large enough */
	    abs_pivot = x ;
	    ppivrow = pdiag ;
	}
	else
	{
	    /* diagonal entry is too small, discard it */
	    pdiag = EMPTY ;
	}
    }

    /* if diagonal not found or not chosen, try for an off-diagonal pivot */
    if (pdiag == EMPTY && ppivrow != EMPTY)
    {
	/* no diagonal pivot.  Try to pick the off-diagonal pivot */
	if (abs_pivot < tol_offdiag * max_entry)
	{
	    /* off-diagonal pivot is too small, discard it */
	    ppivrow = EMPTY ;
	}
    }

    pivot_found = (ppivrow != EMPTY && abs_pivot > 0) ;

    if (pivot_found)
    {
	/* swap pivot row to first entry in column k of L */
	pivrow = Li [ppivrow] ;
	pivot  = Lx [ppivrow] ;
	Li [ppivrow] = Li [0] ;
	Lx [ppivrow] = Lx [0] ;
	Li [0] = pivrow ;
	Lx [0] = 1.0 ;

	ASSERT (pivrow >= 0 && pivrow < npiv) ;
	ASSERT (pivot != 0) ;

	/* divide L by the pivot value */
	for (p = 1 ; p < Llen [k] ; p++)
	{
	    Lx [p] /= pivot ;
	}

	*p_pivrow = pivrow ;
	*p_pivot = pivot ;
	*p_abs_pivot = abs_pivot ;
    }

    return (pivot_found) ;
}


/* ========================================================================== */
/* === prune ================================================================ */
/* ========================================================================== */

/* Prune the columns of L to reduce work in subsequent depth-first searches.
 * This routine is not active during phase 2. */

static void prune
(
    /* TODO comment this */
    int npiv,
    int Lprune [ ],
    int Pinv [ ],
    int k,
    int pivrow,
    double *LU,
    int Uip [ ],
    int Lip [ ],
    int Ulen [ ],
    int Llen [ ]
)
{
    double x ;
    double *Lx, *Ux ;
    int p, i, j, p2, phead, ptail, ulen, llen ;
    int *Li, *Ui ;

    /* ---------------------------------------------------------------------- */
    /* check to see if any column of L can be pruned */
    /* ---------------------------------------------------------------------- */

    GET_COLUMN (Uip, Ulen, LU, k, Ui, Ux, ulen) ;
    ulen-- ;	/* skip the diagonal entry */

    for (p = 0 ; p < ulen ; p++)
    {
	j = Ui [p] ;
	ASSERT (j >= 0 && j < k) ;
	if (Lprune [j] == EMPTY)
	{
	    /* scan column j of L for the pivot row */
	    GET_COLUMN (Lip, Llen, LU, j, Li, Lx, llen) ;

	    for (p2 = 0 ; p2 < Llen [j] ; p2++)
	    {
		if (pivrow == Li [p2])
		{
		    /* found it!  This column can be pruned.
		     * partition column j of L.  The unit diagonal of L
		     * remains as the first entry in the column of L. */
		    phead = 0 ;
		    ptail = Llen [j] ;
		    while (phead < ptail)
		    {
			i = Li [phead] ;
			if (i < npiv && Pinv [i] >= 0)
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

		    /* set Lprune to one past the last entry in the
		     * first part of the column of L.  Entries in
		     * Li [0 ... Lprune [j]-1] are the only part of
		     * column j of L that needs to be scanned in the DFS.
		     * Lprune [j] was EMPTY; setting it >= 0 also flags
		     * column j as pruned. */
		    Lprune [j] = ptail ;
		    break ;
		}
	    }
	}
    }
}


/* ========================================================================== */
/* === paraklete_kernel ===================================================== */
/* ========================================================================== */

/* Factorize P*A*Q into L*U+S.  Returns TRUE if successful, FALSE if out of
 * memory.  If memory runs out, the factors of this node (LUnode) are not
 * freed; that is done by the caller.   If the matrix is singular, the
 * Schur complement (S) has a larger dimension than expected. */

int paraklete_kernel
(
    cholmod_sparse *A,
    paraklete_node *LUnode,
    paraklete_common *Common
)
{ 
    cholmod_common *cm ;
    double pivot, abs_pivot, umin, umax, ujk, x, tol_diag, tol_offdiag, growth ;
    double *Lx, *Ux, *LU, *S, *Sx, *Ax, *X ;
    int *Li, *Ui, *Si, *Qinv, *L1_len, *Ap, *Ai, *Llen, *Ulen, *Lip, *Uip, *P,
	*Q, *Pinv, *Sip, *Slen, *Flag, *Queue, *Iwork, *Stack, *Pstack,
	*Lprune ;
    int k, p, i, j, pivrow, kbar, diagrow, noffdiag, lup, mark, unpruned,
	phead, ptail, sp, slen, llen, ulen, p2, pend, ngrow, qhead, qtail, sn,
	found, kcol, second_try, nfound, n, npiv, lnz, unz, snz, lup_orig ;
    size_t lusize, ssize ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    cm = &(Common->cm) ;	/* CHOLMOD Common */

    /* Matrix to factorize */
    n = A->nrow ;		/* A is n-by-n */
    Ap = A->p ;			/* size n+1, column pointers for A */
    Ai = A->i ;			/* size nz = Ap [n], row indices for A */
    Ax = A->x ;			/* size nz, values of A */

    /* LU factors of this node */
    Llen = LUnode->llen ;	/* size n, column length of L */ 
    Ulen = LUnode->ulen ;	/* size n, column length of U */
    Lip = LUnode->lp ;		/* size n, column pointers for L */
    Uip = LUnode->up ;		/* size n, column pointers for U */
    LU = LUnode->ix ;		/* size LU->size on input and output */
    lusize = LUnode->lusize ;	/* initial and final size of LU */
    npiv = LUnode->PK_NPIV	;	/* number of pivots to find */

    P = LUnode->Plocal ;	/* size npiv, row permutation, where P [k] = i
				 * if row i is the kth pivot row. */
    Q = LUnode->Qlocal ; 	/* size npiv, column permutation */
    Pinv = LUnode->Pinv ;	/* size npiv, inverse row permutation, where
				 * Pinv [i] = k if row i is the kth pivot row */
    Qinv = LUnode->Qinv ;	/* size npiv, inverse of Q */

    tol_diag = Common->tol_diag ;	    /* partial pivoting tolerance */
    tol_offdiag = Common->tol_offdiag ;    /* partial pivoting tolerance */
    growth = Common->growth ;		    /* memory growth factor */

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    cholmod_allocate_work (n, 3*n, n, cm) ;
    if (cm->status != CHOLMOD_OK)
    {
	/* out of memory */
	return (FALSE) ;
    }

    /* workspace, not defined on input or output: */
    Flag = cm->Flag ;		/* size n */
    Queue = cm->Head ;		/* size n+1 */
    Iwork = cm->Iwork ;
    Stack  = Iwork ;		/* size n */
    Pstack = Iwork + n ;	/* size n */
    Lprune = Iwork + 2*n ;	/* size n */

    /* workspace, not defined on input */
    X = cm->Xwork ;		/* size n, undefined on input, zero on output */

    /* ====================================================================== */
    /* PHASE 1: compute L1, L2, and U1 ====================================== */
    /* ====================================================================== */

    /* ---------------------------------------------------------------------- */
    /* initializations */
    /* ---------------------------------------------------------------------- */

    npiv = MAX (0, MIN (npiv, n)) ;
    nfound = 0 ;

    DEBUG (debug_n = n) ;
    DEBUG (debug_npiv = npiv) ;

    lnz = 0 ;
    unz = 0 ;
    snz = 0 ;

    umin = 0 ;
    umax = 0 ;
    noffdiag = 0 ;
    lup = 0 ;
    LUnode->nrealloc = 0 ;
    for (k = 0 ; k < n ; k++)
    {
	X [k] = 0 ;
	Flag [k] = EMPTY ;	/* flag k as unmarked */
	Lprune [k] = EMPTY ;	/* flag k as not pruned */
    }
    mark = 0 ;

    ngrow = 
	(n+1)		/* reals, in L and U (both diag L and U are stored) */
	+ (n+1)/2	/* ints */
	+ 2 ;		/* int to real may be rounded up */

	/*
	(3*n/2) + 2 ;
	*/

    /* ---------------------------------------------------------------------- */
    /* mark all rows as non-pivotal and determine initial diagonal mapping */
    /* ---------------------------------------------------------------------- */

    /* initial Q is identity, so the diagonal entries are the natural ones */
    for (k = 0 ; k < npiv ; k++)
    {
	P [k] = k ;
	Pinv [k] = FLIP (k) ;	/* mark all candidiate rows as non-pivotal */
    }

    /* (P [k] = row) means that (UNFLIP (Pinv [row]) = k), and visa versa.
     * If row is pivotal, then Pinv [row] >= 0.  A row is initially "flipped"
     * (Pinv [k] < EMPTY), and then marked "unflipped" when it becomes
     * pivotal. */

    /* ---------------------------------------------------------------------- */
    /* place all potential pivot columns in the candidate column queue */
    /* ---------------------------------------------------------------------- */

    for (k = 0 ; k < npiv ; k++)
    {
	Queue [k] = k ;
    }
    qhead = 0 ;
    qtail = npiv ;

    /* ---------------------------------------------------------------------- */
    /* factorize up to npiv columns */
    /* ---------------------------------------------------------------------- */

    for (k = 0 ; k < npiv ; k++)
    {

	PR1 ((Common->file, "\n\n==========================L,U1:  k: %d\n", k)) ;
	DEBUG (debug_k = k) ;

	/* ------------------------------------------------------------------ */
	/* determine if LU factors have grown too big */
	/* ------------------------------------------------------------------ */

        /* LU can grow by at most 3n/2 + 2 entries if the column is dense */
	if (lup + ngrow > (int) lusize)
	{
	    size_t new = (growth * lusize + ngrow + 2) ;
	    PR1 ((Common->file, "growing LU from %d to %d\n", lusize, new)) ;
	    LU = cholmod_realloc (new, sizeof (double), LU, &lusize, cm) ;
	    LUnode->ix = LU ;
	    LUnode->lusize = lusize ;
	    if (cm->status != CHOLMOD_OK)
	    {
		/* out of memory */
		return (FALSE) ;
	    }
	}

	/* ------------------------------------------------------------------ */
	/* find a good pivot column, if possible */
	/* ------------------------------------------------------------------ */

	found = FALSE ;
	diagrow = EMPTY ;
	kcol = EMPTY ;
	while (!found && qhead != qtail)
	{

	    /* -------------------------------------------------------------- */
	    /* grab a pivot column candidate from the head of the queue */
	    /* -------------------------------------------------------------- */

	    kcol = Queue [qhead] ;
	    qhead = (qhead + 1) % (npiv + 1) ;
	    second_try = (kcol < 0) ;
	    kcol = UNFLIP (kcol) ;

	    /* -------------------------------------------------------------- */
	    /* s = L \ A (:, kcol) ; store result in L(:,k) and U(:,k)  */
	    /* -------------------------------------------------------------- */

	    lup_orig = lup ;
	    lsolve (TRUE, nfound, npiv,
		n, k, kcol, Ap, Ai, Ax, LU, &lup, Lip, Uip,
		Llen, Ulen, Llen, Pinv, Stack, Lprune, Pstack, Flag, mark, X) ;

	    /* -------------------------------------------------------------- */
	    /* partial pivoting with diagonal preference */
	    /* -------------------------------------------------------------- */

	    /* determine what the "diagonal" is */
	    diagrow = P [k] ;   /* might already be pivotal */
	    PR1 ((Common->file, "k %d, diagrow = %d, UNFLIP(diagrow) = %d\n",
		k, diagrow, UNFLIP (diagrow))) ;

	    /* find a pivot and scale the pivot column */
	    found = lpivot (diagrow, &pivrow, &pivot, &abs_pivot, tol_diag,
		    tol_offdiag, X, LU, Lip, Llen, k, npiv) ;

	    if (!found)
	    {
		/* pivot column candidate failed.  If this was already its
		 * second chance, discard it.  Otherwise, flag it and put it
		 * at the end of the queue.  The kth column of L and U 
		 * computed by lsolve, above, is discarded. */
		if (!second_try)
		{
		    PR1 ((Common->file, "requeue kcol %d for 2nd try\n", kcol)) ;
		    Queue [qtail] = FLIP (kcol) ;
		    qtail = (qtail + 1) % (npiv + 1) ;
		}
		else
		{
		    PR1 ((Common->file, "discarding kcol %d\n", kcol)) ;
		}
		PR1 ((Common->file, "restore lup, %d to %d\n", lup_orig, lup)) ;
		lup = lup_orig ;
	    }
	    else
	    {
		/* we found a pivot column */
		Q [nfound++] = kcol ;
		PR1 ((Common->file, "found kcol: %d nfound %d\n", kcol, nfound));
	    }

	    /* clear the flag array */
	    mark++ ;
	}

	DEBUG (debug_nfound = nfound) ;

	if (!found)
	{
	    /* No more pivots to be found. Go to phase 2 to compute U2 */
	    break ;
	}

	/* ------------------------------------------------------------------ */
	/* we now have a valid pivot row and column */
	/* ------------------------------------------------------------------ */

	PR1 ((Common->file, "\nk %d: pivcol %d pivrow %d: %g\n",
	    k, kcol, pivrow, pivot)) ;
	ASSERT (pivrow >= 0 && pivrow < npiv) ;
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
	/* copy the pivot value as the last entry in the column of U */
	/* ------------------------------------------------------------------ */

	GET_COLUMN (Uip, Ulen, LU, k, Ui, Ux, ulen) ;
	PR1 ((Common->file, "k %d Ulen[k] %d\n", k, Ulen [k])) ;
	ASSERT (ulen > 0) ;
	Ui [ulen-1] = k ;
	Ux [ulen-1] = pivot ;
	PR1 ((Common->file, "pivot: Ui %d Ux %g\n", Ui [ulen-1], Ux [ulen-1])) ;

	/* ------------------------------------------------------------------ */
	/* log the pivot permutation */
	/* ------------------------------------------------------------------ */

	DEBUG (kbar = UNFLIP (Pinv [diagrow])) ;
	ASSERT (kbar < n) ;
	ASSERT (P [kbar] == diagrow) ;

	if (pivrow != diagrow)
	{
	    /* an off-diagonal pivot has been chosen */
	    noffdiag++ ;
	    PR1 ((Common->file, ">>>>>>>>>>>>>>>> pivrow %d k %d off-diagonal\n",
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
	for (p = 0 ; p < ulen ; p++)
	{
	    PR2 ((Common->file, "Column %d of U: %d: %g\n", k, Ui [p], Ux [p])) ;
	}
	GET_COLUMN (Lip, Llen, LU, k, Li, Lx, llen) ;
	for (p = 0 ; p < llen ; p++)
	{
	    PR2 ((Common->file, "Column %d of L: %d: %g\n", k, Li [p], Lx [p])) ;
	}
#endif

	/* ------------------------------------------------------------------ */
	/* symmetric pruning */
	/* ------------------------------------------------------------------ */

	prune (npiv, Lprune, Pinv, k, pivrow, LU, Uip, Lip, Ulen, Llen) ;

	lnz += Llen [k] ;
	unz += Ulen [k] ;
    }

    LUnode->noffdiag = noffdiag ;	/* # of off-diagonal pivots chosen */
    LUnode->umin = umin ;		/* smallest entry on diagonal of U */
    LUnode->umax = umax ;		/* largest entry on the diagonal of U */
    LUnode->PK_NFOUND = nfound ;		/* number of pivots found */
    LUnode->lnz = lnz ;			/* nnz in L */

    PR1 ((Common->file, "noffdiag %d\n", noffdiag)) ;

    /* ====================================================================== */
    /* PHASE 2: compute U2 and S ============================================ */
    /* ====================================================================== */

    /* ---------------------------------------------------------------------- */
    /* finalize the row and column permutations */
    /* ---------------------------------------------------------------------- */

    DEBUG (for (i = 0 ; i < n ; i++) ASSERT (Flag [i] < mark)) ;

    k = nfound ;
    for (i = 0 ; i < npiv ; i++)
    {
	if (Pinv [i] < 0)
	{
	    /* a non-pivotal row */
	    P [k] = i ;
	    Pinv [i] = k++ ;
	}
    }

    for (j = 0 ; j < npiv ; j++)
    {
	Qinv [j] = EMPTY ;
    }
    for (k = 0 ; k < nfound ; k++)
    {
	PR2 ((Common->file, "k %d Q [k] %d npiv %d\n", k, Q [k], npiv)) ;
	ASSERT (Q [k] >= 0 && Q [k] < npiv) ;
	Qinv [Q [k]] = k ;
    }

    k = nfound ;
    for (j = 0 ; j < npiv ; j++)
    {
	if (Qinv [j] == EMPTY)
	{
	    /* a non-pivotal column */
	    Q [k] = j ;
	    Qinv [j] = k++ ;
	}
    }
    ASSERT (k == npiv) ;

#ifndef NDEBUG
    for (i = 0 ; i < npiv ; i++)
    {
	if (i == nfound) PR1 ((Common->file, "-------- nfound = %d\n", nfound)) ;
	if (i == npiv  ) PR1 ((Common->file, "-------- npiv   = %d\n", npiv  )) ;
	PR2 ((Common->file, "P[%d] = %d Pinv[%d] = %d\n", i, P[i], i, Pinv[i])) ;
	PR2 ((Common->file, "Q[%d] = %d Qinv[%d] = %d\n", i, Q[i], i, Qinv[i])) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* partition the columns of L into L1 (pivotal) and L2 (non-pivotal) */
    /* ---------------------------------------------------------------------- */

    /* use Queue as workspace for L1_len */
    L1_len = Queue ;

    for (j = 0 ; j < nfound ; j++)
    {
	GET_COLUMN (Lip, Llen, LU, j, Li, Lx, llen) ;

	/* change all row indices in L1 and L2 to final row indices */
	for (p = 0 ; p < llen ; p++) 
	{
	    i = Li [p] ;
	    if (i < npiv)
	    {
		Li [p] = Pinv [i] ;
	    }
	}

	unpruned = (Lprune [j] == EMPTY) ;

	phead = (unpruned) ? 1 : (Lprune [j]) ;
	ptail = llen ;

#ifndef NDEBUG
	/* any pruned part of the column is already pivotal */
	for (p = 0 ; p < phead ; p++)
	{
	    ASSERT (Li [p] < nfound) ;
	}
#endif

	/* partition row indices in L1 (< nfound) and L2 (>= nfound) */
	while (phead < ptail)
	{
	    i = Li [phead] ;
	    if (i < nfound)
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

	L1_len [j] = ptail ;

	/* prune the column for the next lsolve */
	if (unpruned)
	{
	    Lprune [j] = L1_len [j] ;
	}

#ifndef NDEBUG
	/* the column is now partitioned */
	for (p = 0 ; p < Lprune [j] ; p++)
	{
	    ASSERT (Li [p] < nfound) ;
	}
	ASSERT (Lprune [j] <= L1_len [j]) ;
	for ( ; p < L1_len [j] ; p++)
	{
	    ASSERT (Li [p] < nfound) ;
	}
	for ( ; p < Llen [j] ; p++)
	{
	    ASSERT (Li [p] >= nfound) ;
	}
#endif

    }

    /* ---------------------------------------------------------------------- */
    /* compute the U2 block, U2 = L1 \ A12 */
    /* ---------------------------------------------------------------------- */

    ngrow = 
	(nfound+1)	/* reals, in L and U (both diag L and U are stored) */
	+ (nfound+1)/2	/* ints */
	+ 2 ;		/* int to real may be rounded up */

    for (k = nfound ; k < n ; k++)
    {
	PR1 ((Common->file, "\n\n=============================U2: k: %d\n", k)) ;
	DEBUG (debug_k = k) ;

	/* ------------------------------------------------------------------ */
	/* determine if LU factors have grown too big */
	/* ------------------------------------------------------------------ */

        /* LU can grow by at most 3*nfound/2+2 entries if the column is dense */
	if (lup + ngrow > (int) lusize)
	{
	    size_t new = (growth * lusize + ngrow + 2) ;
	    LU = cholmod_realloc (new, sizeof (double), LU, &lusize, cm) ;
	    LUnode->ix = LU ;
	    LUnode->lusize = lusize ;
	    if (cm->status != CHOLMOD_OK)
	    {
		/* out of memory */
		return (FALSE) ;
	    }
	}

	/* ------------------------------------------------------------------ */
	/* determine which column of A this is */
	/* ------------------------------------------------------------------ */

	kcol = (k < npiv) ? Q [k] : k ;

	/* ------------------------------------------------------------------ */
	/* s = L \ A (:, kcol) and store result in U(:,k)  */
	/* ------------------------------------------------------------------ */

	lsolve (FALSE, nfound, npiv,
	    n, k, kcol, Ap, Ai, Ax, LU, &lup, Lip, Uip,
	    Llen, Ulen, L1_len, Pinv, Stack, Lprune, Pstack, Flag, mark, X) ;

	PR2 ((Common->file, "lup is now %d\n", lup)) ;

#ifndef NDEBUG
	GET_COLUMN (Uip, Ulen, LU, k, Ui, Ux, ulen) ;
	for (p = 0 ; p < ulen ; p++)
	{
	    PR2 ((Common->file, "Column %d of U2: %d: %g\n", k, Ui[p], Ux[p])) ;
	}
#endif

	unz += Ulen [k] ;

	/* clear the flag array */
	mark++ ;
    }

    LUnode->unz = unz ;			/* nnz in U */

    /* ---------------------------------------------------------------------- */
    /* shrink the LU factors to just the required size */
    /* ---------------------------------------------------------------------- */

    LU = cholmod_realloc (lup, sizeof (double), LU, &lusize, cm) ;
    LUnode->ix = LU ; 
    LUnode->lusize = lusize ;
    ASSERT (cm->status == CHOLMOD_OK)

    /* ---------------------------------------------------------------------- */
    /* compute the Schur complement, S = A22 - L2*U2 */
    /* ---------------------------------------------------------------------- */

    sn = n - nfound ;
    ssize = lusize ;
    LUnode->PK_SSIZE = ssize ;
    PR1 ((Common->file, "allocate S, : ssize %d sn %d \n", ssize, sn)) ;

    S = cholmod_malloc (ssize, sizeof (double), cm) ;
    Sip = cholmod_malloc (sn, sizeof (int), cm) ;
    Slen = cholmod_malloc (sn, sizeof (int), cm) ;
    LUnode->sx = S ;
    LUnode->sp = Sip ;
    LUnode->slen = Slen ;
    LUnode->PK_SSIZE = ssize ; 
    LUnode->PK_SN = sn ;
    if (cm->status != CHOLMOD_OK)
    {
	/* out of memory */
	return (FALSE) ;
    }

    sp = 0 ;
    ngrow = (3 * (n - nfound) / 2) + 2 ;

    /* S is square, of order n-nfound, with rows/cols in range 0 : n-nfound-1 */

    for (k = nfound ; k < n ; k++)
    {
	PR2 ((Common->file, "\n\n========================== S:  k: %d sk: %d\n",
		    k, k-nfound)) ;
	DEBUG (for (i = 0 ; i < n ; i++) ASSERT (Flag [i] < mark)) ;

	/* ------------------------------------------------------------------ */
	/* determine if the Schur complement needs to grow */
	/* ------------------------------------------------------------------ */

	if (sp + ngrow > (int) ssize)
	{
	    size_t new = (growth * ssize + ngrow + 2) ;
	    S = cholmod_realloc (new, sizeof (double), S, &ssize, cm) ;
	    LUnode->sx = S ;
	    LUnode->PK_SSIZE = ssize ;
	    if (cm->status != CHOLMOD_OK)
	    {
		return (FALSE) ;
	    }
	}

	/* ------------------------------------------------------------------ */
	/* compute the kth column of S */
	/* ------------------------------------------------------------------ */

	Sip [k - nfound] = sp ;
	Si = (int *) (S + sp) ;

	/* ------------------------------------------------------------------ */
	/* scatter the kth column of A*Q (namely, A(:,kcol)) into X */
	/* ------------------------------------------------------------------ */

	kcol = (k < npiv) ? Q [k] : k ;
	slen = 0 ;
	DEBUG (for (i = 0 ; i < n ; i++) ASSERT (X [i] == 0)) ;

	/* scatter only the non-pivotal rows of A (for the S part only) */
	pend = Ap [kcol+1] ;
	for (p = Ap [kcol] ; p < pend ; p++)
	{
	    i = Ai [p] ;
	    i = (i < npiv) ? Pinv [i] : i ;
	    if (i >= nfound)
	    {
		PR3 ((Common->file, "S: Entry in A: %d %g\n", i, Ax [p])) ;
		X [i] = Ax [p] ;
		Flag [i] = mark ;
		Si [slen++] = i - nfound ;
	    }
	}

	/* ------------------------------------------------------------------ */
	/* X = X - L2 * U2 (:,k) */
	/* ------------------------------------------------------------------ */

	/* get the pointers for the kth column of U */
	GET_COLUMN (Uip, Ulen, LU, k, Ui, Ux, ulen) ;
	PR2 ((Common->file, "Multiply L2 by U2(:,%d) :\n", k)) ;

	/* for each j for which U (j,k) is nonzero, do: */
	for (p = 0 ; p < ulen ; p++)
	{
	    /* X = X - L (nfound:n, j) * U (j,k) */
	    j = Ui [p] ;
	    ASSERT (j >= 0 && j < nfound) ;
	    ujk = Ux [p] ;
	    PR2 ((Common->file, "U(%d,%d) = %g\n", j, k, ujk)) ;
	    GET_COLUMN (Lip, Llen, LU, j, Li, Lx, llen) ;
	    for (p2 = L1_len [j] ; p2 < llen ; p2++)
	    {
		i = Li [p2] ;
		ASSERT (i >= nfound && i < n) ;
		PR3 ((Common->file, "    L(%d,%d) = %g\n", i, j, Lx [p2])) ;
		if (Flag [i] < mark)
		{
		    PR3 ((Common->file, "    new entry in Si, slen %d\n", slen));
		    Si [slen++] = i - nfound ;
		    Flag [i] = mark ;
		}
		X [i] -= Lx [p2] * ujk ;
	    }
	}

	/* ------------------------------------------------------------------ */
	/* allocate space for S(:,k) */
	/* ------------------------------------------------------------------ */

	Slen [k - nfound] = slen ;
	PR2 ((Common->file, "Slen [%d] = %d\n", k - nfound, Slen [k - nfound])) ;
	snz += slen ;
	sp += ((slen + 1) / 2) ;
	Sx = (double *) (S + sp) ;
	sp += slen ;

	/* ------------------------------------------------------------------ */
	/* gather X into S and clear X */
	/* ------------------------------------------------------------------ */

	PR2 ((Common->file, "snz so far: %d. Final col of S(:,%d)=\n", snz, k)) ;
	for (p = 0 ; p < slen ; p++)
	{
	    i = Si [p] + nfound ;
	    Sx [p] = X [i] ;
	    X [i] = 0 ;
	    PR3 ((Common->file, " S(%d,%d) = %g\n", i-nfound, k-nfound, Sx [p])) ;
	}
	DEBUG (for (i = 0 ; i < n ; i++) ASSERT (X [i] == 0)) ;

	/* clear the flag array */
	mark++ ;
    }

    /* ---------------------------------------------------------------------- */
    /* shrink the Schur complement to just the required size */
    /* ---------------------------------------------------------------------- */

    PR0 ((Common->file, "Final ssize = %d sn = %d snz: %d\n", sp, sn, snz)) ;
    S = cholmod_realloc (sp, sizeof (double), S, &ssize, cm) ;
    LUnode->sx = S ; 
    LUnode->PK_SNZ = snz ; 
    LUnode->PK_SSIZE = ssize ;

    return (TRUE) ;
}
