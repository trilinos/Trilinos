/* ========================================================================== */
/* === Core/cholmod_postorder =============================================== */
/* ========================================================================== */

/*
 * CHOLMOD/Core version 0.1. May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* Compute the postorder of a tree.  This is part of the Core module because
 * several unrelated CHOLMOD modules rely on it.  It is not a primary routine.
 */

#include "cholmod_core.h"
#include "cholmod_internal.h"

/* ========================================================================== */
/* === dfs ================================================================== */
/* ========================================================================== */

/* The code below includes both a recursive and non-recursive depth-first-search
 * of a tree.  The recursive code is simpler, but can lead to stack overflow.
 * It is left here for reference, to understand what the non-recursive code
 * is computing.  To try the recursive version, uncomment the following
 * #define, or compile the code with -DRECURSIVE.  Be aware that stack
 * overflow may occur.
#define RECURSIVE
 */

#ifdef RECURSIVE

/* recursive version: a working code for reference only, not actual use */

static int dfs			/* return the new value of k */
(
    int p,		/* start a DFS at node p */
    int k,		/* start the node numbering at k */
    int Post [ ],	/* Post ordering, modified on output */
    int Head [ ],	/* Head [p] = youngest child of p; EMPTY on output */
    int Next [ ],	/* Next [j] = sibling of j; unmodified */
    int Pstack [ ]	/* unused */
)
{
    int j ;
    /* start a DFS at each child of node p */
    for (j = Head [p] ; j != EMPTY ; j = Next [j])
    {
	/* start a DFS at child node j */
	k = dfs (j, k, Post, Head, Next, Pstack) ;
    }
    Post [k++] = p ;	/* order node p as the kth node */
    Head [p] = EMPTY ;	/* link list p no longer needed */
    return (k) ;	/* the next node will be numbered k */
}

#else

/* non-recursive version */

static int dfs			/* return the new value of k */
(
    int p,		/* start the DFS at a root node p */
    int k,		/* start the node numbering at k */
    int Post [ ],	/* Post ordering, modified on output */
    int Head [ ],	/* Head [p] = youngest child of p; EMPTY on output */
    int Next [ ],	/* Next [j] = sibling of j; unmodified */
    int Pstack [ ]	/* workspace of size n, undefined on input or output */
)
{
    int j, phead ;

    /* put the root node on the stack */
    Pstack [0] = p ;
    phead = 0 ;

    /* while the stack is not empty, do: */
    while (phead >= 0)
    {
	/* grab the node p from top of the stack and get its youngest child j */
	p = Pstack [phead] ;
	j = Head [p] ;
	if (j == EMPTY)
	{
	    /* all children of p ordered.  remove p from stack and order it */
	    phead-- ;
	    Post [k++] = p ;	/* order node p as the kth node */
	}
	else
	{
	    /* leave p on the stack.  Start a DFS at child node j by putting
	     * j on the stack and removing j from the list of children of p. */
	    Head [p] = Next [j] ;
	    Pstack [++phead] = j ;
	}
    }
    return (k) ;	/* the next node will be numbered k */
}

#endif

/* ========================================================================== */
/* === cholmod_postorder ==================================================== */
/* ========================================================================== */

/* Postorder a tree.  The tree is either an elimination tree (the output from
 * from cholmod_etree) or a component tree (from cholmod_nd).
 *
 * An elimination tree is a complete tree of n nodes with Parent [j] > j or
 * Parent [j] = EMPTY if j is a root.  On output Post [0..n-1] is a complete
 * permutation vector.
 *
 * A component tree is a subset of 0..n-1.  Parent [j] = -2 if node j is not
 * in the component tree.  Parent [j] = EMPTY if j is a root of the component
 * tree, and Parent [j] is in the range 0 to n-1 if j is in the component
 * tree but not a root.  On output, Post [j] is defined only for nodes in
 * the component tree.  Post [j] = k if node j is the kth node in the
 * postordered component tree, where k is in the range 0 to the number of
 * components minus 1.
 *
 * As a result, check_parent (Parent, n,...) may fail on input, since
 * cholmod_check_parent assumes Parent is an elimination tree.  Similarly,
 * cholmod_check_perm (Post, ...) may fail on output, since Post is a partial
 * permutation if Parent is a component tree.
 *
 * workspace: Head (n), Iwork (2*n)
 */

long cholmod_postorder	/* returns # of nodes postordered, or -1 if error */
(
    /* inputs, not modified on output: */
    void *Parent_p,	/* size n.  Parent [j] = p if p is the parent of j.
			 * node j is a root if Parent [j] is EMPTY.
			 * node j is ignored and not included in the postorder
			 * if Parent [j] < EMPTY. */
    size_t n,		/* nodes are in the range 0 to n-1 */

    /* outputs, not defined on input */
    void *Post_p,	/* size n.  Post [k] = j if j is kth node in
			 * postordered tree */

    cholmod_common *Common
)
{
    int *Head, *Next, *Pstack, *Parent, *Post, *Iwork ;
    int j, p, k ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (EMPTY) ;
    Parent = Parent_p ;
    Post = Post_p ;
    RETURN_IF_NULL (Parent, EMPTY) ;
    RETURN_IF_NULL (Post, EMPTY) ;
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    cholmod_allocate_work (n, 2*n, 0, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (EMPTY) ;
    }
    ASSERT (cholmod_dump_work (TRUE, TRUE, 0, Common)) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    Head  = Common->Head ;	/* size n+1, initially all EMPTY */
    Iwork = Common->Iwork ;
    Next  = Iwork ;		/* size n (i/i/l) */
    Pstack = Iwork + n ;	/* size n (i/i/l) */

    /* ---------------------------------------------------------------------- */
    /* construct a link list of children for each node */
    /* ---------------------------------------------------------------------- */

    /* in reverse order so the children are in ascending order in each list */
    for (j = n-1 ; j >= 0 ; j--)
    {
	p = Parent [j] ;
	if (p >= 0 && p < ((int) n))
	{
	    /* add j to the list of children for node p */
	    Next [j] = Head [p] ;
	    Head [p] = j ;
	}
    }

    /* Head [p] = j if j is the youngest (least-numbered) child of p */
    /* Next [j1] = j2 if j2 is the next-oldest sibling of j1 */

    /* ---------------------------------------------------------------------- */
    /* start a DFS at each root node of the etree */
    /* ---------------------------------------------------------------------- */

    k = 0 ;
    for (j = 0 ; j < ((int) n) ; j++)
    {
	if (Parent [j] == EMPTY)
	{
	    /* j is the root of a tree; start a DFS here */
	    k = dfs (j, k, Post, Head, Next, Pstack) ;
	}
    }

    /* this would normally be EMPTY already, unless Parent is invalid */
    for (j = 0 ; j < ((int) n) ; j++)
    {
	Head [j] = EMPTY ;
    }

    PRINT1 (("postordered %d nodes\n", k)) ;
    ASSERT (cholmod_dump_work (TRUE, TRUE, 0, Common)) ;
    return (k) ;
}
