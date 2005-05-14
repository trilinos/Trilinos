/* ========================================================================== */
/* === Partition/cholmod_nd ================================================= */
/* ========================================================================== */

/*
 * CHOLMOD/Partition version 0.1. May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* CHOLMOD nested dissection and graph partitioning.
 *
 * cholmod_bisect:
 *
 *	Finds a set of nodes that partitions the graph into two parts.
 *	Compresses the graph first.  Requires METIS.
 *
 * cholmod_nested_dissection:
 *
 *	Nested dissection, using its own compression and connected-commponents
 *	algorithms, an external graph partitioner (METIS), and a constrained
 *	minimum degree ordering algorithm (CCOLAMD or CSYMAMD).  Typically
 *	gives better orderings than METIS_NodeND (about 5% to 10% fewer
 *	nonzeros in L).
 *
 *	FUTURE WORK: CSYMAMD is not a "native" symmetric ordering, but recasts
 *	the problem for CCOLAMD.  It gives good orderings, but the related
 *	SYMAMD can be up to 4 times slower than COLAMD.  A constrained AMD
 *	would be faster.

 * cholmod_metis:
 *
 *	A wrapper for METIS_NodeND.

 * This file contains several routines private to this file:
 *
 *	partition	compress and partition a graph
 *	clear_flag	clear Common->Flag, but do not modify negative entries
 *	find_components	find the connected components of a graph
 */

#include "cholmod_partition.h"
#include "cholmod_internal.h"

/* ========================================================================== */
/* === partition ============================================================ */
/* ========================================================================== */

/* Find a set of nodes that partition a graph.  The graph must be symmetric
 * with no diagonal entries.  To compress the graph first, compress is TRUE
 * and on input Hash [j] holds the hash key for node j, which must be in the
 * range 0 to csize-1. The input graph (Cp, Ci) is destroyed.  Cew is all 1's
 * on input and output.  Cnw [j] > 0 is the initial weight of node j.  On
 * output, Cnw [i] = 0 if node i is absorbed into j and the original weight
 * Cnw [i] is added to Cnw [j].  If compress is FALSE, the graph is not
 * compressed and Cnw and Hash are unmodified.  The partition itself is held in
 * the output array Part of size n.  Part [j] is 0, 1, or 2, depending on
 * whether node j is in the left part of the graph, the right part, or the
 * separator, respectively.  Note that the input graph need not be connected,
 * and the output subgraphs (the three parts) may also be unconnected.
 *
 * The size of the separator is guaranteed to be between 1 and n nodes.
 * If it is of size less than n, then both the left and right parts are
 * guaranteed to be non-empty.
 */

static int partition	/* size of separator or -1 if failure */
(
    /* inputs, not modified on output */
#ifndef NDEBUG
    int csize,		/* upper bound on # of edges in the graph;
			 * csize >= MAX (n, nnz(C)) must hold. */
#endif
    int compress,	/* if TRUE the compress the graph first */

    /* input/output */
    int Hash [ ],	/* Hash [i] = hash >= 0 is the hash function for node
			 * i on input.  On output, Hash [i] = FLIP (j) if node
			 * i is absorbed into j.  Hash [i] >= 0 if i has not
			 * been absorbed. */

    /* input graph, compressed graph of cn nodes on output */
    cholmod_sparse *C,

    /* input/output */
    int Cnw [ ],	/* size n.  Cnw [j] > 0 is the weight of node j on
			 * input.  On output, if node i is absorbed into
			 * node j, then Cnw [i] = 0 and the original weight of
			 * node i is added to Cnw [j].  The sum of Cnw [0..n-1]
			 * is not modified. */

    /* workspace */
    int Cew [ ],	/* size csize, all 1's on input and output */

    /* more workspace, undefined on input and output */
    int Cmap [ ],	/* size n (i/i/l) */

    /* output */
    int Part [ ],	/* size n, Part [j] = 0, 1, or 2. */

    cholmod_common *Common
)
{
    int n, hash, head, i, j, k, p, pend, ilen, ilast, pi, piend,
	jlen, ok, cn, csep, pdest, nodes_pruned, nz, total_weight,
	jscattered ;
    int *Cp, *Ci, *Next, *Hhead ;

#ifndef NDEBUG
    int cnt, pruned ;
    double work = 0, goodwork = 0 ;
#endif

    /* ---------------------------------------------------------------------- */
    /* quick return for empty graphs */
    /* ---------------------------------------------------------------------- */

    n = C->nrow ;
    if (n <= 0)
    {
	/* nothing to do */
	return (0) ;
    }

    Cp = C->p ;
    Ci = C->i ;
    nz = Cp [n] ;

    PRINT1 (("Partition start, n %d nz %d\n", n, nz)) ;

    if (n == 1 || nz <= 0)
    {
	/* no edges, this is easy */
	PRINT1 (("diagonal matrix\n")) ;
	k = n/2 ;
	for (j = 0 ; j < k ; j++)
	{
	    Part [j] = 0 ;
	}
	for ( ; j < n ; j++)
	{
	    Part [j] = 1 ;
	}
	/* ensure the separator is not empty (required by nested dissection) */
	Part [n-1] = 2 ;
	return (1) ;
    }

    total_weight = 0 ;
    for (j = 0 ; j < n ; j++)
    {
	ASSERT (Cnw [j] > 0) ;
	total_weight += Cnw [j] ;
    }

#ifndef NDEBUG
    ASSERT (n > 1 && nz > 0) ;
    PRINT1 (("original graph:\n")) ;
    for (j = 0 ; j < n ; j++)
    {
	PRINT2 (("%d: ", j)) ;
	for (p = Cp [j] ; p < Cp [j+1] ; p++)
	{
	    i = Ci [p] ;
	    PRINT3 (("%d ", i)) ;
	    ASSERT (i >= 0 && i < n && i != j) ;
	}
	PRINT2 (("hash: %d\n", Hash [j])) ;
    }
    DEBUG (for (p = 0 ; p < csize ; p++) ASSERT (Cew [p] == 1)) ;
#endif

    nodes_pruned = 0 ;

    if (compress)
    {

	/* ------------------------------------------------------------------ */
	/* get workspace */
	/* ------------------------------------------------------------------ */

	Next = Part ;	/* use Part as workspace for Next [ */
	Hhead = Cew ;	/* use Cew as workspace for Hhead [ */

	/* ------------------------------------------------------------------ */
	/* create the hash buckets */
	/* ------------------------------------------------------------------ */

	for (j = 0 ; j < n ; j++)
	{
	    /* get the hash key for node j */
	    hash = Hash [j] ;
	    ASSERT (hash >= 0 && hash < csize) ;
	    head = Hhead [hash] ;
	    if (head > EMPTY)
	    {
		/* hash bucket for this hash key is empty. */
		head = EMPTY ;
	    }
	    else
	    {
		/* hash bucket for this hash key is not empty.  get old head */
		head = FLIP (head) ;
		ASSERT (head >= 0 && head < n) ;
	    }
	    /* node j becomes the new head of the hash bucket.  FLIP it so that
	     * we can tell the difference between an empty or non-empty hash
	     * bucket. */
	    Hhead [hash] = FLIP (j) ;
	    Next [j] = head ;
	    ASSERT (head >= EMPTY && head < n) ;
	}

#ifndef NDEBUG
	for (cnt = 0, k = 0 ; k < n ; k++)
	{
	    ASSERT (Hash [k] >= 0 && Hash [k] < csize) ;    /* k is alive */
	    hash = Hash [k] ;
	    ASSERT (hash >= 0 && hash < csize) ;
	    head = Hhead [hash] ;
	    ASSERT (head < EMPTY) ;	/* hash bucket not empty */
	    j = FLIP (head) ;
	    ASSERT (j >= 0 && j < n) ;
	    if (j == k)
	    {
		PRINT2 (("hash %d: ", hash)) ;
		for ( ; j != EMPTY ; j = Next [j])
		{
		    PRINT3 ((" %d", j)) ;
		    ASSERT (j >= 0 && j < n) ;
		    ASSERT (Hash [j] == hash) ;
		    cnt++ ;
		    ASSERT (cnt <= n) ;
		}
		PRINT2 (("\n")) ;
	    }
	}
	ASSERT (cnt == n) ;
#endif

	/* ------------------------------------------------------------------ */
	/* scan the non-empty hash buckets for indistinguishable nodes */
	/* ------------------------------------------------------------------ */

	/* If there are no hash collisions and no compression occurs, this takes
	 * O(n) time.  If no hash collisions, but some nodes are removed, this
	 * takes time O(n+e) where e is the sum of the degress of the nodes
	 * that are removed.  Even with many hash collisions (a rare case),
	 * this algorithm has never been observed to perform more than nnz(A)
	 * useless work.
	 *
	 * Cmap is used as workspace to mark nodes of the graph, [
	 * for comparing the nonzero patterns of two nodes i and j.
	 */

#define MARK(i)   Cmap [i] = j
#define MARKED(i) (Cmap [i] == j)

	for (i = 0 ; i < n ; i++)
	{
	    Cmap [i] = EMPTY ;
	}

	for (k = 0 ; k < n ; k++)
	{
	    hash = Hash [k] ;
	    ASSERT (hash >= FLIP (n-1) && hash < csize) ;
	    if (hash < 0)
	    {
		/* node k has already been absorbed into some other node */
		ASSERT (FLIP (Hash [k]) >= 0 && FLIP (Hash [k] < n)) ;
		continue ;
	    }
	    head = Hhead [hash] ;
	    ASSERT (head < EMPTY || head == 1) ;
	    if (head == 1)
	    {
		/* hash bucket is already empty */
		continue ;
	    }
	    PRINT2 (("\n--------------------hash %d:\n", hash)) ;
	    for (j = FLIP (head) ; j != EMPTY && Next[j] > EMPTY ; j = Next [j])
	    {
		/* compare j with all nodes i following it in hash bucket */
		ASSERT (j >= 0 && j < n && Hash [j] == hash) ;
		p = Cp [j] ;
		pend = Cp [j+1] ;
		jlen = pend - p ;
		jscattered = FALSE ;
		DEBUG (for (i = 0 ; i < n ; i++) ASSERT (!MARKED (i))) ;
		DEBUG (pruned = FALSE) ;
		ilast = j ;
		for (i = Next [j] ; i != EMPTY ; i = Next [i])
		{
		    ASSERT (i >= 0 && i < n && Hash [i] == hash && i != j) ;
		    pi = Cp [i] ;
		    piend = Cp [i+1] ;
		    ilen = piend - pi ;
		    DEBUG (work++) ;
		    if (ilen != jlen)
		    {
			/* i and j have different degrees */
			ilast = i ;
			continue ;
		    }
		    /* scatter the pattern of node j, if not already */
		    if (!jscattered)
		    {
			MARK (j) ;
			for ( ; p < pend ; p++)
			{
			    MARK (Ci [p]) ;
			}
			jscattered = TRUE ;
			DEBUG (work += jlen) ;
		    }
		    for (ok = MARKED (i) ; ok && pi < piend ; pi++)
		    {
			ok = MARKED (Ci [pi]) ;
			DEBUG (work++) ;
		    }
		    if (ok)
		    {
			/* found it.  kill node i and merge it into j */
			PRINT2 (("found %d absorbed into %d\n", i, j)) ;
			Hash [i] = FLIP (j) ;
			Cnw [j] += Cnw [i] ;
			Cnw [i] = 0 ;
			ASSERT (ilast != i && ilast >= 0 && ilast < n) ;
			Next [ilast] = Next [i] ; /* delete i from bucket */
			nodes_pruned++ ;
			DEBUG (goodwork += (ilen+1)) ;
			DEBUG (pruned = TRUE) ;
		    }
		    else
		    {
			/* i and j are different */
			ilast = i ;
		    }
		}
		DEBUG (if (pruned) goodwork += jlen) ;
	    }
	    /* empty the hash bucket, restoring Cew */
	    Hhead [hash] = 1 ;
	}

	DEBUG (if (((work - goodwork) / (double) nz) > 0.20) PRINT0 ((
	    "work %12g good %12g nz %12d (wasted work/nz: %6.2f )\n",
	    work, goodwork, nz, (work - goodwork) / ((double) nz)))) ;

	/* All hash buckets now empty.  Cmap no longer needed as workspace. ]
	 * Cew no longer needed as Hhead; Cew is now restored to all ones. ]
	 * Part no longer needed as workspace for Next. ] */
    }

    /* Edge weights are all one, node weights reflect node absorption */
    DEBUG (for (p = 0 ; p < csize ; p++) ASSERT (Cew [p] == 1)) ;
    DEBUG (for (cnt = 0, j = 0 ; j < n ; j++) cnt += Cnw [j]) ;
    ASSERT (cnt == total_weight) ;

    /* ---------------------------------------------------------------------- */
    /* compress and partition the graph */
    /* ---------------------------------------------------------------------- */

    if (nodes_pruned == 0)
    {

	/* ------------------------------------------------------------------ */
	/* no pruning done at all.  Do not create the compressed graph */
	/* ------------------------------------------------------------------ */

	/* FUTURE WORK: could call CHACO here too */
	csep = cholmod_metis_bisector (C, Cnw, Cew, Part, Common) ;

    }
    else if (nodes_pruned == n-1)
    {

	/* ------------------------------------------------------------------ */
	/* only one node left.  This is a dense graph */
	/* ------------------------------------------------------------------ */

	PRINT1 (("completely dense graph\n")) ;
	csep = total_weight ;
	for (j = 0 ; j < n ; j++)
	{
	    Part [j] = 2 ;
	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* compress the graph and partition the compressed graph */
	/* ------------------------------------------------------------------ */

	/* ------------------------------------------------------------------ */
	/* create the map from the uncompressed graph to the compressed graph */
	/* ------------------------------------------------------------------ */

	/* Cmap [j] = k if node j is alive and the kth node of compressed graph.
	 * The mapping is done monotonically (that is, k <= j) to simplify the
	 * uncompression later on.  Cmap [j] = EMPTY if node j is dead. */

	for (j = 0 ; j < n ; j++)
	{
	    Cmap [j] = EMPTY ;
	}
	k = 0 ;
	for (j = 0 ; j < n ; j++)
	{
	    if (Cnw [j] > 0)
	    {
		ASSERT (k <= j) ;
		Cmap [j] = k++ ;
	    }
	}
	cn = k ;	    /* # of nodes in compressed graph */
	PRINT1 (("compressed graph from %d to %d nodes\n", n, cn)) ;
	ASSERT (cn > 1 && cn == n - nodes_pruned) ;

	/* ------------------------------------------------------------------ */
	/* create the compressed graph */
	/* ------------------------------------------------------------------ */

	k = 0 ;
	pdest = 0 ;
	for (j = 0 ; j < n ; j++)
	{
	    if (Cnw [j] > 0)
	    {
		/* node j in the full graph is node k in the compressed graph */
		ASSERT (k <= j && Cmap [j] == k) ;
		p = Cp [j] ;
		pend = Cp [j+1] ;
		Cp [k] = pdest ;
		Cnw [k] = Cnw [j] ;
		for ( ; p < pend ; p++)
		{
		    /* prune dead nodes, and remap to new node numbering */
		    i = Ci [p] ;
		    ASSERT (i >= 0 && i < n && i != j) ;
		    i = Cmap [i] ;
		    ASSERT (i >= EMPTY && i < cn && i != k) ;
		    if (i > EMPTY)
		    {
			ASSERT (pdest <= p) ;
			Ci [pdest++] = i ;
		    }
		}
		k++ ;
	    }
	}
	Cp [cn] = pdest ;
	C->nrow = cn ;
	C->ncol = cn ;	/* affects mem stats unless restored when C free'd */

#ifndef NDEBUG
	PRINT1 (("pruned graph (%d/%d) nodes, (%d/%d) edges\n",
		    cn, n, pdest, nz)) ;
	PRINT1 (("compressed graph:\n")) ;
	for (cnt = 0, j = 0 ; j < cn ; j++)
	{
	    PRINT2 (("%d: ", j)) ;
	    for (p = Cp [j] ; p < Cp [j+1] ; p++)
	    {
		i = Ci [p] ;
		PRINT3 (("%d ", i)) ;
		ASSERT (i >= 0 && i < cn && i != j) ;
	    }
	    PRINT2 (("weight: %d\n", Cnw [j])) ;
	    ASSERT (Cnw [j] > 0) ;
	    cnt += Cnw [j] ;
	}
	ASSERT (cnt == total_weight) ;
	for (j = 0 ; j < n ; j++) PRINT2 (("Cmap [%d] = %d\n", j, Cmap [j])) ;
	ASSERT (k == cn) ;
#endif

	/* ------------------------------------------------------------------ */
	/* find the separator of the compressed graph */
	/* ------------------------------------------------------------------ */

	/* FUTURE WORK: could call CHACO here too */
	csep = cholmod_metis_bisector (C, Cnw, Cew, Part, Common) ;

	if (csep < 0)
	{
	    /* failed */
	    return (-1) ;
	}

	PRINT1 (("Part: ")) ;
	DEBUG (for (j = 0 ; j < cn ; j++) PRINT2 (("%d ", Part [j]))) ;
	PRINT1 (("\n")) ;

	/* Cp and Ci no longer needed */

	/* ------------------------------------------------------------------ */
	/* find the separator of the uncompressed graph */
	/* ------------------------------------------------------------------ */

	/* expand the separator to live nodes in the uncompressed graph */
	for (j = n-1 ; j >= 0 ; j--)
	{
	    /* do this in reverse order so that Cnw can be expanded in place */
	    k = Cmap [j] ;
	    ASSERT (k >= EMPTY && k < n) ;
	    if (k > EMPTY)
	    {
		/* node k in compressed graph and is node j in full graph */
		ASSERT (k <= j) ;
		ASSERT (Hash [j] >= EMPTY) ;
		Part [j] = Part [k] ;
		Cnw [j] = Cnw [k] ;
	    }
	    else
	    {
		/* node j is a dead node */
		Cnw [j] = 0 ;
		DEBUG (Part [j] = EMPTY) ;
		ASSERT (Hash [j] < EMPTY) ;
	    }
	}

	/* find the components for the dead nodes */
	for (i = 0 ; i < n ; i++)
	{
	    if (Hash [i] < EMPTY)
	    {
		/* node i has been absorbed into node j */
		j = FLIP (Hash [i]) ;
		ASSERT (Part [i] == EMPTY && j >= 0 && j < n && Cnw [i] == 0) ;
		Part [i] = Part [j] ;
	    }
	    ASSERT (Part [i] >= 0 && Part [i] <= 2) ;
	}

#ifndef NDEBUG
	PRINT1 (("Part: ")) ;
	for (cnt = 0, j = 0 ; j < n ; j++)
	{
	    ASSERT (Part [j] != EMPTY) ;
	    PRINT2 (("%d ", Part [j])) ;
	    if (Part [j] == 2) cnt += Cnw [j] ;
	}
	PRINT1 (("\n")) ;
	PRINT1 (("csep %d %d\n", cnt, csep)) ;
	ASSERT (cnt == csep) ;
	for (cnt = 0, j = 0 ; j < n ; j++) cnt += Cnw [j] ;
	ASSERT (cnt == total_weight) ;
#endif

    }

    /* ---------------------------------------------------------------------- */
    /* return the separator (or -1 if error) */
    /* ---------------------------------------------------------------------- */

    PRINT1 (("Partition done, n %d csep %d\n", n, csep)) ;
    return (csep) ;
}


/* ========================================================================== */
/* === clear_flag =========================================================== */
/* ========================================================================== */

/* A node j has been removed from the graph if Flag [j] < EMPTY.
 * If Flag [j] >= EMPTY && Flag [j] < mark, then node j is alive but unmarked.
 * Flag [j] == mark means that node j is alive and marked.  Incrementing mark
 * means that all nodes are either (still) dead, or live but unmarked.
 *
 * On output, Common->mark < Common->Flag [i] for all i from 0 to Common->nrow.
 * This is the same output condition as cholmod_clear_flag, except that this
 * routine maintains the Flag [i] < EMPTY condition as well, if that condition
 * was true on input.
 *
 * workspace: Flag (nrow)
 */

static long clear_flag (cholmod_common *Common)
{
    int nrow, i ;
    int *Flag ;
    PRINT2 (("old mark %ld\n", Common->mark)) ;
    Common->mark++ ;
    PRINT2 (("new mark %ld\n", Common->mark)) ;
    if (Common->mark <= 0)
    {
	nrow = Common->nrow ;
	Flag = Common->Flag ;
	for (i = 0 ; i < nrow ; i++)
	{
	    /* if Flag [i] < EMPTY, leave it alone */
	    if (Flag [i] >= EMPTY)
	    {
		Flag [i] = EMPTY ;
	    }
	}
	/* now Flag [i] <= EMPTY for all i */
	Common->mark = 0 ;
    }
    return (Common->mark) ;
}


/* ========================================================================== */
/* === find_components ====================================================== */
/* ========================================================================== */

/* Find all connected components of the current subgraph C.  The subgraph C
 * consists of the nodes of B that appear in the set Map [0..cn-1].  If Map
 * is NULL, then it is assumed to be the identity mapping
 * (Map [0..cn-1] = 0..cn-1).
 *
 * A node j does not appear in B if it has been ordered (Flag [j] < EMPTY,
 * which means that j has been ordered and is "deleted" from B).
 *
 * If the size of a component is large, it is placed on the component stack,
 * Cstack.  Otherwise, its nodes are ordered and it is not placed on the Cstack.
 *
 * A component S is defined by a "representative node" (repnode for short)
 * called the snode, which is one of the nodes in the subgraph.  Likewise, the
 * subgraph C is defined by its repnode, called cnode.
 *
 * workspace: Flag (nrow)
 */

static void find_components
(
    /* inputs, not modified on output */
    cholmod_sparse *B,
    int Map [ ],	    /* size n, only Map [0..cn-1] used */
    int cn,		    /* # of nodes in C */
    int cnode,		    /* root node of component C, or EMPTY if C is the
			     * entire graph B */

    /* input/output */
    int Bnz [ ],	    /* size n.  Bnz [j] = # nonzeros in column j of B.
			     * Reduce since B is pruned of dead nodes. */

    int CParent [ ],	    /* CParent [i] = j if component with repnode j is
			     * the parent of the component with repnode i.
			     * CParent [i] = EMPTY if the component with
			     * repnode i is a root of the separator tree.
			     * CParent [i] is -2 if i is not a repnode. */
    int Cstack [ ],	    /* component stack for nested dissection */
    int *top,		    /* Cstack [0..top] contains root nodes of the
			     * the components currently in the stack */

    /* workspace, undefined on input and output: */
    int Queue [ ],	    /* size n, for breadth-first search */

    cholmod_common *Common
)
{

    int n, mark, cj, j, sj, sn, p, i, snode, pstart, pdest, pend, nd_small ;
    int *Bp, *Bi, *Flag ;

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    Flag = Common->Flag ;	    /* size n */
    Common->mark = EMPTY ;	    /* force initialization of Flag array */
    mark = clear_flag (Common) ;    /* clear Flag but preserve Flag [i]<EMPTY */
    Bp = B->p ;
    Bi = B->i ;
    n = B->nrow ;
    ASSERT (cnode >= EMPTY && cnode < n) ;
    ASSERT (IMPLIES (cnode >= 0, Flag [cnode] < EMPTY)) ;

    /* get ordering parameters */
    nd_small = Common->method [Common->current].nd_small ;
    nd_small = MAX (4, nd_small) ;

    /* ---------------------------------------------------------------------- */
    /* find the connected components of C via a breadth-first search */
    /* ---------------------------------------------------------------------- */

    for (cj = 0 ; cj < cn ; cj++)
    {
	/* get node snode, which is node cj of C.  It might already be in the
	 * separator of C (and thus ordered, with Flag [snode] < EMPTY) */
	snode = (Map == NULL) ? (cj) : (Map [cj]) ;
	ASSERT (snode >= 0 && snode < n) ;

	if (Flag [snode] >= EMPTY && Flag [snode] < mark)
	{

	    /* -------------------------------------------------------------- */
	    /* find new connected component S */
	    /* -------------------------------------------------------------- */

	    /* node snode is the repnode of a connected component S, the
	     * parent of which is cnode, the repnode of C.  If cnode is
	     * EMPTY then C is the original graph B. */
	    PRINT1 (("--------------:::snode %d cnode %d\n", snode, cnode)) ;
	    ASSERT (CParent [snode] == -2) ;
	    CParent [snode] = cnode ;

	    /* place j in the queue and mark it */
	    sj = 0 ;
	    Queue [0] = snode ;
	    Flag [snode] = mark ;
	    sn = 1 ;

	    /* breadth-first traversal, starting at node j */
	    for (sj = 0 ; sj < sn ; sj++)
	    {
		/* get node j from head of Queue and traverse its edges */
		j = Queue [sj] ;
		PRINT2 (("    j: %d\n", j)) ;
		ASSERT (j >= 0 && j < n) ;
		ASSERT (Flag [j] == mark) ;
		pstart = Bp [j] ;
		pdest = pstart ;
		pend = pstart + Bnz [j] ;
		for (p = pstart ; p < pend ; p++)
		{
		    i = Bi [p] ;
		    if (i != j && Flag [i] >= EMPTY)
		    {
			/* node is still in the graph */
			Bi [pdest++] = i ;
			if (Flag [i] < mark)
			{
			    /* node i is in this component S, and is unflagged
			     * (first time node i has been seen in this BFS).
			     * place node i in the queue and mark it */
			    Queue [sn++] = i ;
			    Flag [i] = mark ;
			}
		    }
		}
		/* edges to dead nodes have been removed */
		Bnz [j] = pdest - pstart ;
	    }

	    /* -------------------------------------------------------------- */
	    /* order S if it is small; place it on Cstack otherwise */
	    /* -------------------------------------------------------------- */

	    PRINT2 (("sn %d\n", sn)) ;
	    if (sn <= nd_small)
	    {
		/* the S component is tiny, order it now */
		for (sj = 0 ; sj < sn ; sj++)
		{
		    j = Queue [sj] ;
		    ASSERT (j >= 0 && j < n) ;
		    ASSERT (Flag [j] == mark) ;
		    Flag [j] = FLIP (snode) ;
		}
		ASSERT (Flag [snode] == FLIP (snode)) ;
	    }
	    else
	    {
		/* place the new component on the Cstack */
		Cstack [++(*top)] = snode ;
	    }
	}
    }

    /* clear Flag array, but preserve Flag [i] < EMPTY */
    (void) clear_flag (Common) ;
}


/* ========================================================================== */
/* === cholmod_bisect ======================================================= */
/* ========================================================================== */

/* Finds a node bisector of A, A*A', A(:,f)*A(:,f)'.
 *
 * workspace: Flag (nrow), Iwork (nrow + max (nrow,ncol)).
 *	Allocates a temporary matrix B=A*A' or B=A,
 *	and O(nnz(A)) temporary memory space.
 */

long cholmod_bisect	/* # of nodes in separator or -1 if error */
(
    /* input only, not modified on output: */
    cholmod_sparse *A,		/* must be square, need not be sorted */
    void *fset,
    size_t fsize,
    int compress,		/* if TRUE, compress the graph first */

    /* output, allocated but contents undefined on input */
    void *Partition,		/* size A->nrow.  Node i is in the left graph
				 * if Partition [i] = 0, the right graph if 1,
				 * and in the separator if 2. */

    cholmod_common *Common
)
{
    int *Bp, *Bi, *Hash, *Cmap, *Bnw, *Bew, *Iwork ;
    cholmod_sparse *B ;
    unsigned int hash ;
    int j, n, bnz, csize, sepsize, p, pend ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (EMPTY) ;
    RETURN_IF_NULL (A, EMPTY) ;
    RETURN_IF_NULL (Partition, EMPTY) ;
    n = A->nrow ;
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    cholmod_allocate_work (n, n + MAX (n, ((int) (A->ncol))), 0, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (-1) ;
    }
    ASSERT (cholmod_dump_work (TRUE, TRUE, 0, Common)) ;

    Iwork = Common->Iwork ;
    Hash = Iwork ;		/* size n, FUTURE WORK: (i/l/l) */
    Cmap = Iwork + n ;		/* size n, (i/i/l) */

    /* ---------------------------------------------------------------------- */
    /* convert the matrix to adjacency list form */
    /* ---------------------------------------------------------------------- */

    /* The input graph to must be symmetric, with no diagonal entries
     * present.  The columns need not be sorted. */

    /* B = A+A', A*A', or A(:,f)*A(:,f)', upper and lower parts present */

    if (A->stype)
    {
	/* Add the upper/lower part to a symmetric lower/upper matrix by
	 * converting to unsymmetric mode */
	/* workspace: Iwork (max (nrow,ncol)) */
	B = cholmod_copy (A, 0, -1, Common) ;
    }
    else
    {
	/* B = A*A' or A(:,f)*A(:,f)', no diagonal */
	/* workspace: Flag (nrow), Iwork (max (nrow,ncol)) */
	B = cholmod_aat (A, fset, fsize, -1, Common) ;
    }

    if (Common->status < CHOLMOD_OK)
    {
	return (EMPTY) ;
    }
    Bp = B->p ;
    Bi = B->i ;
    bnz = Bp [n] ;
    ASSERT ((int) (B->nrow) == n && (int) (B->ncol) == n) ;

    /* Bew should be at least size n for the hash function to work well */
    csize = MAX (n+1, bnz) ;

    /* create the graph using Flag as workspace for node weights [ */
    Bnw = Common->Flag ;    /* size n workspace */

    /* compute hash for each node if compression requested */
    if (compress)
    {
	for (j = 0 ; j < n ; j++)
	{
	    hash = j ;
	    pend = Bp [j+1] ;
	    for (p = Bp [j] ; p < pend ; p++)
	    {
		hash += Bi [p] ;
		ASSERT (Bi [p] != j) ;
	    }
	    /* finalize the hash key for node j */
	    hash %= csize ;
	    Hash [j] = (int) hash ;
	    ASSERT (Hash [j] >= 0 && Hash [j] < csize) ;
	}
    }

    /* allocate edge weights */
    Bew = cholmod_malloc (csize, sizeof (int), Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	/* out of memory */
	cholmod_free_sparse (&B, Common) ;
	cholmod_free (Bew, csize, sizeof (int), Common) ;
	return (EMPTY) ;
    }

    /* graph has unit node and edge weights */
    for (j = 0 ; j < n ; j++)
    {
	Bnw [j] = 1 ;
    }
    for (p = 0 ; p < csize ; p++)
    {
	Bew [p] = 1 ;
    }

    /* ---------------------------------------------------------------------- */
    /* compress and partition the graph */
    /* ---------------------------------------------------------------------- */

    sepsize = partition (
#ifndef NDEBUG
	    csize,
#endif
	    compress, Hash, B, Bnw, Bew, Cmap, Partition, Common) ;

    /* contents of Bp, Bi, Bnw, and Bew no longer needed ] */

    /* If partition fails, free the workspace below and return sepsize < 0 */

    /* ---------------------------------------------------------------------- */
    /* free workspace */
    /* ---------------------------------------------------------------------- */

    B->ncol = n ;   /* restore size for memory usage statistics */
    cholmod_free_sparse (&B, Common) ;
    Common->mark = EMPTY ;
    cholmod_clear_flag (Common) ;
    cholmod_free (Bew, csize, sizeof (int), Common) ;
    return (sepsize) ;
}


/* ========================================================================== */
/* === cholmod_nested_dissection -=========================================== */
/* ========================================================================== */

/* This method uses a node bisector, applied recursively (but using a
 * non-recursive algorithm).  Once the graph is partitioned, it calls a
 * constrained min degree code (CSYMAMD for A+A', and CCOLAMD for A*A') to
 * order all the nodes in the graph - but obeying the constraints determined
 * by the separators.  This routine is similar to METIS_NodeND, except for how
 * it treats the leaf nodes.  METIS_NodeND orders the leaves of the separator
 * tree with MMD, ignoring the rest of the matrix when ordering a single leaf.
 * This routine orders the whole matrix with CSYMAMD or CCOLAMD, all at once,
 * when the graph partitioning is done.
 *
 * workspace: Flag (nrow), Head (nrow+1), Iwork (nrow + max (nrow,ncol)).
 *	Allocates a temporary matrix B=A*A' or B=A,
 *	and O(nnz(A)) temporary memory space.
 */

long cholmod_nested_dissection	/* returns # of components, or -1 if error */
(
    cholmod_sparse *A,
    void *fset,
    size_t fsize,

    /* outputs, contents undefined on input */
    void *Perm_p,	/* size n.  Perm [k] = j if node j is kth node in the
			 * permuted matrix */
    void *CParent_p,	/* size n.  On output, CParent [c] is the parent
			 * of component c, or EMPTY if c is a root.  c is in
			 * the range 0 to the # of components minus 1. */
    void *Cmember_p,	/* size n.  Cmember [j] = c if node j of A is
			 * in component c */

    cholmod_common *Common
)
{
    double nd_prune ;
    int *Bp, *Bi, *Bnz, *Cstack, *Imap, *Map, *Flag, *Head, *Next, *Bnw, *Iwork,
	*Ipost, *NewParent, *Hash, *Cmap, *Cp, *Ci, *Cew, *Cnw, *Part, *Post,
	*Perm, *CParent, *Cmember ;
    unsigned int hash ;
    int n, bnz, top, i, j, k, cnode, p, cj, cn, ci, cnz, mark, c,
	sepsize, parent, ncomponents, threshold, ndense, pstart, pdest, pend,
	symmetric, ok, nd_compress, nd_camd, csize, jnext ;
    cholmod_sparse *B, *C ;

    DEBUG (int cnt) ;
    DEBUG (int orig) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (EMPTY) ;
    RETURN_IF_NULL (A, EMPTY) ;
    Perm = Perm_p ;
    CParent = CParent_p ;
    Cmember = Cmember_p ;
    RETURN_IF_NULL (Perm, EMPTY) ;
    RETURN_IF_NULL (CParent, EMPTY) ;
    RETURN_IF_NULL (Cmember, EMPTY) ;
    n = A->nrow ;
    Common->status = CHOLMOD_OK ;
    if (A->nrow == 0)
    {
	return (1) ;	    /* nothing to do */
    }
    ASSERT (cholmod_dump_sparse (A, "A_ND:", Common) >= 0) ;
    symmetric = A->stype ;

    /* get ordering parameters */
    nd_prune = Common->method [Common->current].nd_prune ;
    nd_compress = Common->method [Common->current].nd_compress ;
    nd_camd = Common->method [Common->current].nd_camd ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    cholmod_allocate_work (n, n + MAX (n, ((int) (A->ncol))), 0, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (EMPTY) ;
    }
    DEBUG (orig = Common->malloc_count) ;
    PRINT1 (("orig in nd: %d\n", orig)) ;
    ASSERT (cholmod_dump_work (TRUE, TRUE, 0, Common)) ;

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    /* Note:  cholmod_analyze uses Common->Xwork for Cmember and CParent,
     * so do not use it here. */

    Cstack = Perm ;		/* use Perm as workspace for Cstack */
    Flag = Common->Flag ;	/* size n */
    Head = Common->Head ;	/* size n+1, all equal to -1 */
    Iwork = Common->Iwork ;
    Imap = Iwork ;		/* size n, same as Queue in find_components */
    Map  = Iwork + n ;		/* size n */

    Bnz  = cholmod_malloc (n, sizeof (int), Common) ;
    Hash = cholmod_malloc (n, sizeof (int), Common) ;
    Cmap = cholmod_malloc (n, sizeof (int), Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	cholmod_free (Bnz,  n, sizeof (int), Common) ;
	cholmod_free (Hash, n, sizeof (int), Common) ;
	cholmod_free (Cmap, n, sizeof (int), Common) ;
	ASSERT (Common->malloc_count == orig) ;
	return (EMPTY) ;
    }

    /* ---------------------------------------------------------------------- */
    /* convert B to symmetric form with both upper/lower parts present */
    /* ---------------------------------------------------------------------- */

    /* B = A+A', A*A', or A(:,f)*A(:,f)', upper and lower parts present */

    if (A->stype)
    {
	/* Add the upper/lower part to a symmetric lower/upper matrix by
	 * converting to unsymmetric mode */
	/* workspace: Iwork (max (nrow,ncol)) */
	B = cholmod_copy (A, 0, -1, Common) ;
    }
    else
    {
	/* B = A*A' or A(:,f)*A(:,f)', no diagonal */
	/* workspace: Flag (nrow), Iwork (max (nrow,ncol)) */
	B = cholmod_aat (A, fset, fsize, -1, Common) ;
    }

    if (Common->status < CHOLMOD_OK)
    {
	cholmod_free (Bnz,  n, sizeof (int), Common) ;
	cholmod_free (Hash, n, sizeof (int), Common) ;
	cholmod_free (Cmap, n, sizeof (int), Common) ;
	ASSERT (Common->malloc_count == orig) ;
	return (EMPTY) ;
    }
    Bp = B->p ;
    Bi = B->i ;
    bnz = cholmod_nnz (B, Common) ;
    ASSERT ((int) (B->nrow) == n && (int) (B->ncol) == n) ;
    csize = MAX (n, bnz) ;
    ASSERT (cholmod_dump_sparse (B, "B for nd:", Common) >= 0) ;

    /* ---------------------------------------------------------------------- */
    /* initializations */
    /* ---------------------------------------------------------------------- */

    /* all nodes start out unmarked and unordered (Type 4, see below) */
    Common->mark = EMPTY ;
    cholmod_clear_flag (Common) ;
    ASSERT (Flag == Common->Flag) ;
    ASSERT (cholmod_dump_work (TRUE, TRUE, 0, Common)) ;

    for (j = 0 ; j < n ; j++)
    {
	CParent [j] = -2 ;
    }

    /* prune dense nodes from B */
    if (ISNAN (nd_prune) || nd_prune < 0)
    {
	threshold = n-1 ;
    }
    else
    {
	threshold = (int) (MAX (16, nd_prune * sqrt ((double) (n)))) ;
	threshold = MIN (n-1, threshold) ;
    }
    ndense = 0 ;
    cnode = EMPTY ;

    for (j = 0 ; j < n ; j++)
    {
	Bnz [j] = Bp [j+1] - Bp [j] ;
	if (Bnz [j] >= threshold)
	{
	    /* node j is dense, prune it from B */
	    ndense++ ;
	    if (cnode == EMPTY)
	    {
		/* first dense node found becomes root of this component,
		 * which contains all of the dense nodes found here */
		cnode = j ;
		CParent [cnode] = EMPTY ;
	    }
	    Flag [j] = FLIP (cnode) ;
	}
    }
    B->packed = FALSE ;

    if (ndense == n)
    {
	/* all nodes removed: Perm is identity, all nodes in component zero,
	 * and the separator tree has just one node. */
	PRINT1 (("all nodes are dense\n")) ;
	for (k = 0 ; k < n ; k++)
	{
	    Perm [k] = k ;
	    Cmember [k] = 0 ;
	}
	CParent [0] = EMPTY ;
	cholmod_free_sparse (&B, Common) ;
	cholmod_free (Bnz,  n, sizeof (int), Common) ;
	cholmod_free (Hash, n, sizeof (int), Common) ;
	cholmod_free (Cmap, n, sizeof (int), Common) ;
	ASSERT (Common->malloc_count == orig) ;
	Common->mark = EMPTY ;
	cholmod_clear_flag (Common) ;
	return (1) ;
    }

    /* Cp and Ci are workspace to construct the subgraphs to partition */
    C = cholmod_allocate_sparse (n, n, csize, FALSE, TRUE, 0, FALSE, Common) ;
    Cew  = cholmod_malloc (csize, sizeof (int), Common) ;
    Cnw  = cholmod_malloc (n, sizeof (int), Common) ;
    Part = cholmod_malloc (n, sizeof (int), Common) ;
    Bnw  = cholmod_malloc (n, sizeof (int), Common) ;

    if (Common->status < CHOLMOD_OK)
    {
	/* out of memory */
	cholmod_free_sparse (&C, Common) ;
	cholmod_free_sparse (&B, Common) ;
	cholmod_free (Cew, csize, sizeof (int), Common) ;
	cholmod_free (Bnz,  n, sizeof (int), Common) ;
	cholmod_free (Hash, n, sizeof (int), Common) ;
	cholmod_free (Cmap, n, sizeof (int), Common) ;
	cholmod_free (Cnw,  n, sizeof (int), Common) ;
	cholmod_free (Part, n, sizeof (int), Common) ;
	cholmod_free (Bnw,  n, sizeof (int), Common) ;
	ASSERT (Common->malloc_count == orig) ;
	Common->mark = EMPTY ;
	cholmod_clear_flag (Common) ;
	PRINT1 (("out of memory for C, etc\n")) ;
	return (EMPTY) ;
    }

    Cp = C->p ;
    Ci = C->i ;

    /* create initial unit node and edge weights */
    for (j = 0 ; j < n ; j++)
    {
	Bnw [j] = 1 ;
    }
    for (p = 0 ; p < csize ; p++)
    {
	Cew [p] = 1 ;
    }

    /* push the initial connnected components of B onto the Cstack */
    top = EMPTY ;	/* Cstack is empty */
    /* workspace: Flag (nrow), Iwork (nrow); use Imap as workspace for Queue [*/
    find_components (B, NULL, n, cnode, Bnz, CParent, Cstack, &top, Imap,
	    Common) ;
    /* done using Imap as workspace for Queue ] */

    /* Nodes can now be of Type 0, 1, 2, or 4 (see definition below) */

    /* ---------------------------------------------------------------------- */
    /* while Cstack is not empty, do: */
    /* ---------------------------------------------------------------------- */

    while (top >= 0)
    {

	/* ------------------------------------------------------------------ */
	/* get a node from the top of the Cstack */
	/* ------------------------------------------------------------------ */

	/* cnode is the repnode of its (unordered) connected component. */
	cnode = Cstack [top--] ;
	ASSERT (cnode >= 0 && cnode < n && Flag [cnode] >= EMPTY) ;
	ASSERT (CParent [cnode] >= EMPTY && CParent [cnode] < n) ;

	/* During ordering, there are five kinds of nodes in the graph of B,
	 * based on Flag [j] and CParent [j] for nodes j = 0 to n-1:
	 *
	 * Type 0: If cnode is a repnode of an unordered component, then
	 * CParent [cnode] is in the range EMPTY to n-1 and
	 * Flag [cnode] >= EMPTY.  This is a "live" node.
	 *
	 * Type 1: If cnode is a repnode of an ordered separator component,
	 * then Flag [cnode] < EMPTY and FLAG [cnode] = FLIP (cnode).
	 * CParent [cnode] is in the range EMPTY to n-1.  cnode is a root of
	 * the separator tree if CParent [cnode] == EMPTY.  This node is dead.
	 *
	 * Type 2: If node j isn't a repnode, has not been absorbed via
	 * graph compression into another node, but is in an ordered separator
	 * component, then cnode = FLIP (Flag [j]) gives the repnode of the
	 * component that contains j and CParent [j]  is -2.  This node is dead.
	 * Note that Flag [j] < EMPTY.
	 *
	 * Type 3: If node i has been absorbed via graph compression into some
	 * other node j = FLIP (Flag [i]) where j is not a repnode.
	 * CParent [j] is -2.  Node i may or may not be in an ordered
	 * component.  This node is dead.  Note that Flag [j] < EMPTY.
	 *
	 * Type 4: If node j is "live" (not in an ordered component, and not
	 * absorbed into any other node), then Flag [j] >= EMPTY.
	 *
	 * Only "live" nodes (of type 0 or 4) are placed in a subgraph to be
	 * partitioned.  Node j is alive if Flag [j] >= EMPTY, and dead if
	 * Flag [j] < EMPTY.
	 */

	/* ------------------------------------------------------------------ */
	/* create the subgraph for this connected component C */
	/* ------------------------------------------------------------------ */

	/* Do a breadth-first search of the graph starting at cnode.
	 * use Map [0..cn-1] for nodes in the component C [
	 * use Cnw and Cew for node and edge weights of the resulting subgraph [
	 * use Cp and Ci for the resulting subgraph [
	 * use Imap [i] for all nodes i in B that are in the component C [
	 */

	/* clear the Flag array, but do not modify negative entries in Flag  */
	mark = clear_flag (Common) ;
	DEBUG (for (i = 0 ; i < n ; i++) Imap [i] = EMPTY) ;

	/* place cnode in the queue and mark it */
	Map [0] = cnode ;
	Flag [cnode] = mark ;
	Imap [cnode] = 0 ;
	cn = 1 ;

	cnz = 0 ;
	for (cj = 0 ; cj < cn ; cj++)
	{
	    /* get node j from the head of the queue; it is node cj of C */
	    j = Map [cj] ;
	    ASSERT (Flag [j] == mark) ;
	    Cp [cj] = cnz ;
	    Cnw [cj] = Bnw [j] ;
	    ASSERT (Cnw [cj] >= 0) ;
	    pstart = Bp [j] ;
	    pdest = pstart ;
	    pend = pstart + Bnz [j] ;
	    hash = cj ;
	    for (p = pstart ; p < pend ; p++)
	    {
		i = Bi [p] ;
		/* prune diagonal entries and dead edges from B */
		if (i != j && Flag [i] >= EMPTY)
		{
		    /* live node i is in the current component */
		    Bi [pdest++] = i ;
		    if (Flag [i] != mark)
		    {
			/* First time node i has been seen, it is a new node
			 * of C.  place node i in the queue and mark it */
			Map [cn] = i ;
			Flag [i] = mark ;
			Imap [i] = cn ;
			cn++ ;
		    }
		    /* place the edge (cj,ci) in the adjacency list of cj */
		    ci = Imap [i] ;
		    ASSERT (ci >= 0 && ci < cn && ci != cj && cnz < csize) ;
		    Ci [cnz++] = ci ;
		    hash += ci ;
		}
	    }
	    /* edges to dead nodes have been removed */
	    Bnz [j] = pdest - pstart ;
	    /* finalize the hash key for column j */
	    hash %= csize ;
	    Hash [cj] = (int) hash ;
	    ASSERT (Hash [cj] >= 0 && Hash [cj] < csize) ;
	}
	Cp [cn] = cnz ;
	C->nrow = cn ;
	C->ncol = cn ;	/* affects mem stats unless restored when C free'd */

	/* contents of Imap no longer needed ] */

#ifndef NDEBUG
	for (cj = 0 ; cj < cn ; cj++)
	{
	    j = Map [cj] ;
	    PRINT2 (("---------------------------------C column cj: %d j: %d\n",
		cj, j)) ;
	    ASSERT (j >= 0 && j < n) ;
	    ASSERT (Flag [j] >= EMPTY) ;
	    for (p = Cp [cj] ; p < Cp [cj+1] ; p++)
	    {
		ci = Ci [p] ;
		i = Map [ci] ;
		PRINT3 (("ci: %d i: %d\n", ci, i)) ;
		ASSERT (ci != cj && ci >= 0 && ci < cn) ;
		ASSERT (i != j && i >= 0 && i < n) ;
		ASSERT (Flag [i] >= EMPTY) ;
	    }
	}
#endif

	/* small components are never placed on the stack */
	ASSERT (cn > MAX (4, Common->method [Common->current].nd_small)) ;

	/* Cp and Ci now contain the component, with cn nodes and cnz nonzeros.
	 * The mapping of a node cj into node j the main graph B is given by
	 * Map [cj] = j */

	/* ------------------------------------------------------------------ */
	/* compress and partition the graph C */
	/* ------------------------------------------------------------------ */

	/* The edge weights Cew [0..csize-1] are all 1's on input to and output
	 * from the partition routine. */

	sepsize = partition (
#ifndef NDEBUG
		csize,
#endif
		nd_compress, Hash, C, Cnw, Cew,
		Cmap, Part, Common) ;

	/* contents of Cp and Ci no longer needed ] */

	if (sepsize < 0)
	{
	    /* failed */
	    C->ncol = n ;   /* restore size for memory usage statistics */
	    cholmod_free_sparse (&C, Common) ;
	    cholmod_free_sparse (&B, Common) ;
	    cholmod_free (Cew, csize, sizeof (int), Common) ;
	    cholmod_free (Bnz,  n, sizeof (int), Common) ;
	    cholmod_free (Hash, n, sizeof (int), Common) ;
	    cholmod_free (Cmap, n, sizeof (int), Common) ;
	    cholmod_free (Cnw,  n, sizeof (int), Common) ;
	    cholmod_free (Part, n, sizeof (int), Common) ;
	    cholmod_free (Bnw,  n, sizeof (int), Common) ;
	    Common->mark = EMPTY ;
	    cholmod_clear_flag (Common) ;
	    ASSERT (Common->malloc_count == orig) ;
	    return (EMPTY) ;
	}

	/* ------------------------------------------------------------------ */
	/* compress B based on how C was compressed */
	/* ------------------------------------------------------------------ */

	for (ci = 0 ; ci < cn ; ci++)
	{
	    if (Hash [ci] < EMPTY)
	    {
		/* ci is dead in C, having been absorbed into cj */
		cj = FLIP (Hash [ci]) ;
		PRINT1 (("In C, %d absorbed into %d (wgt now %d)\n", ci, cj,
			    Cnw [cj])) ;
		/* i is dead in B, having been absorbed into j */
		i = Map [ci] ;
		j = Map [cj] ;
		PRINT1 (("In B, %d (wgt %d) absorbed into %d (wgt %d -> %d)\n",
			    i, Bnw [i], j, Bnw [j], Cnw [cj])) ;
		/* more than one node may be absorbed into j.  This is accounted
		 * for in Cnw [cj].  Assign it here rather than += Bnw [i] */
		Bnw [i] = 0 ;
		Bnw [j] = Cnw [cj] ;
		Flag [i] = FLIP (j) ;
	    }
	}

	DEBUG (for (cnt = 0, j = 0 ; j < n ; j++) cnt += Bnw [j]) ;
	ASSERT (cnt == n) ;

	/* contents of Cnw [0..cn-1] no longer needed ] */

	/* ------------------------------------------------------------------ */
	/* order the separator, and stack the components when C is split */
	/* ------------------------------------------------------------------ */

	/* one more component has been found: either the separator of C,
	 * or all of C */

	if (sepsize == cn || sepsize == 0)
	{
	    /* Order the nodes in the component.  The separator is too large,
	     * or empty.  Note that the partition routine cannot return a
	     * sepsize of zero, but it can return a separator consisting of the
	     * whole graph.  The "sepsize == 0" test is kept, above, in case the
	     * partition routine changes.  In either case, this is component
	     * remains unsplit, and becomes a leaf of the separator tree. */
	    PRINT1 (("sepsize zero or all of graph: %d\n", sepsize)) ;
	    for (cj = 0 ; cj < cn ; cj++)
	    {
		j = Map [cj] ;
		Flag [j] = FLIP (cnode) ;
		PRINT2 (("      node cj: %d j: %d ordered\n", cj, j)) ;
	    }
	    ASSERT (cnode == Map [0]) ;
	    ASSERT (cnode != EMPTY && Flag [cnode] < EMPTY) ;
	}
	else
	{
	    /* Order the nodes in the separator of C and find a new repnode
	     * cnode that is in the separator of C.  This requires the separator
	     * to be non-empty. */
	    PRINT1 (("sepsize not tiny: %d\n", sepsize)) ;
	    parent = CParent [cnode] ;
	    ASSERT (parent >= EMPTY && parent < n) ;
	    CParent [cnode] = -2 ;
	    cnode = EMPTY ;
	    for (cj = 0 ; cj < cn ; cj++)
	    {
		j = Map [cj] ;
		if (Part [cj] == 2)
		{
		    /* All nodes in the separator become part of a component
		     * whose repnode is cnode */
		    PRINT2 (("node cj: %d j: %d ordered\n", cj, j)) ;
		    if (cnode == EMPTY)
		    {
			PRINT2(("------------new cnode: cj %d j %d\n", cj, j)) ;
			cnode = j ;
		    }
		    Flag [j] = FLIP (cnode) ;
		}
		else
		{
		    PRINT2 (("      node cj: %d j: %d not ordered\n", cj, j)) ;
		}
	    }
	    ASSERT (cnode != EMPTY && Flag [cnode] < EMPTY) ;
	    ASSERT (CParent [cnode] == -2) ;
	    CParent [cnode] = parent ;

	    /* find the connected components when C is split, and push
	     * then on the Cstack.  Use Imap as workspace for Queue. [ */
	    /* workspace: Flag (nrow) */
	    find_components (B, Map, cn, cnode, Bnz,
		    CParent, Cstack, &top, Imap, Common) ;
	    /* done using Imap as workspace for Queue ] */
	}
	/* contents of Map [0..cn-1] no longer needed ] */
    }

    /* ---------------------------------------------------------------------- */
    /* place nodes removed via compression into their proper component */
    /* ---------------------------------------------------------------------- */

    /* At this point, all nodes are of Type 1, 2, or 3, as defined above. */

    for (i = 0 ; i < n ; i++)
    {
	/* find the repnode cnode that contains node i */
	j = FLIP (Flag [i]) ;
	PRINT2 (("\nfind component for %d, in: %d\n", i, j)) ;
	ASSERT (j >= 0 && j < n) ;
	DEBUG (cnt = 0) ;
	while (CParent [j] == -2)
	{
	    j = FLIP (Flag [j]) ;
	    PRINT2 (("    walk up to %d ", j)) ;
	    ASSERT (j >= 0 && j < n) ;
	    PRINT2 ((" CParent %d\n", CParent [j])) ;
	    ASSERT (cnt < n) ;
	    DEBUG (cnt++) ;
	}
	cnode = j ;
	ASSERT (cnode >= 0 && cnode < n) ;
	ASSERT (CParent [cnode] >= EMPTY && CParent [cnode] < n) ;
	PRINT2 (("i %d is in component with cnode %d\n", i, cnode)) ;
	ASSERT (Flag [cnode] == FLIP (cnode)) ;

	/* Mark all nodes along the path from i to cnode as being in the
	 * component whos repnode is cnode.  Perform path compression.  */
	j = FLIP (Flag [i]) ;
	Flag [i] = FLIP (cnode) ;
	DEBUG (cnt = 0) ;
	while (CParent [j] == -2)
	{
	    ASSERT (j >= 0 && j < n) ;
	    jnext = FLIP (Flag [j]) ;
	    PRINT1 (("    %d walk %d set cnode to %d\n", i, j, cnode)) ;
	    ASSERT (cnt < n) ;
	    DEBUG (cnt++) ;
	    Flag [j] = FLIP (cnode) ;
	    j = jnext ;
	}
    }

    /* At this point, all nodes fall into Types 1 or 2, as defined above. */

#ifndef NDEBUG
    for (j = 0 ; j < n ; j++)
    {
	if (CParent [j] >= EMPTY && CParent [j] < n)
	{
	    /* case 1: j is a repnode of a component */
	    cnode = j ;
	}
	else
	{
	    /* case 2: j is not a repnode of a component */
	    cnode = FLIP (Flag [j]) ;
	    ASSERT (cnode >= 0 && cnode < n) ;
	    ASSERT (CParent [cnode] >= EMPTY && CParent [cnode] < n) ;
	}
	ASSERT (Flag [cnode] == FLIP (cnode)) ;
	/* case 3 no longer holds */
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* free workspace */
    /* ---------------------------------------------------------------------- */

    C->ncol = n ;   /* restore size for memory usage statistics */
    cholmod_free_sparse (&C, Common) ;
    cholmod_free_sparse (&B, Common) ;
    cholmod_free (Cew, csize, sizeof (int), Common) ;
    cholmod_free (Bnz,  n, sizeof (int), Common) ;
    cholmod_free (Hash, n, sizeof (int), Common) ;
    cholmod_free (Cmap, n, sizeof (int), Common) ;
    cholmod_free (Cnw,  n, sizeof (int), Common) ;
    cholmod_free (Part, n, sizeof (int), Common) ;
    cholmod_free (Bnw,  n, sizeof (int), Common) ;
    ASSERT (Common->malloc_count == orig) ;

    /* ---------------------------------------------------------------------- */
    /* postorder the components */
    /* ---------------------------------------------------------------------- */

    DEBUG (for (cnt = 0, j = 0 ; j < n ; j++) if (CParent [j] != -2) cnt++) ;

    /* use Cmember as workspace for Post [ */
    Post = Cmember ;

    /* cholmod_postorder uses Head and Iwork [0..2n].  It does not use Flag,
     * which here holds the mapping of nodes to repnodes.  It ignores all nodes
     * for which CParent [j] < -1, so it operates just on the repnodes. */
    /* workspace: Head (n), Iwork (2*n) */
    ncomponents = cholmod_postorder (CParent, n, Post, Common) ;
    ASSERT (cnt == ncomponents) ;

    /* use Iwork [0..n-1] as workspace for Ipost ( */
    Ipost = Iwork ;
    DEBUG (for (j = 0 ; j < n ; j++) Ipost [j] = EMPTY) ;

    /* compute inverse postorder */
    for (c = 0 ; c < ncomponents ; c++)
    {
	cnode = Post [c] ;
	ASSERT (cnode >= 0 && cnode < n) ;
	Ipost [cnode] = c ;
	ASSERT (Head [c] == EMPTY) ;
    }

    /* adjust the parent array */
    /* Iwork [n..2n-1] used for NewParent [ */
    NewParent = Iwork + n ;
    for (c = 0 ; c < ncomponents ; c++)
    {
	parent = CParent [Post [c]] ;
	NewParent [c] = (parent == EMPTY) ? EMPTY : (Ipost [parent]) ;
    }
    for (c = 0 ; c < ncomponents ; c++)
    {
	CParent [c] = NewParent [c] ;
    }
    ASSERT (cholmod_dump_parent (CParent, ncomponents, "CParent", Common)) ;

    /* Iwork [n..2n-1] no longer needed for NewParent ] */
    /* Cmember no longer needed for Post ] */

    /* ---------------------------------------------------------------------- */
    /* place each node in its component */
    /* ---------------------------------------------------------------------- */

    for (j = 0 ; j < n ; j++)
    {
	/* node j is in the cth component, whose root node is cnode */
	cnode = FLIP (Flag [j]) ;
	PRINT2 (("j %d  flag %d cnode %d\n", j, Flag [j], FLIP (Flag [j]))) ;
	ASSERT (cnode >= 0 && cnode < n) ;
	c = Ipost [cnode] ;
	ASSERT (c >= 0 && c < ncomponents) ;
	Cmember [j] = c ;
    }

    /* Flag no longer needed for the node-to-component mapping */

    /* done using Iwork [0..n-1] as workspace for Ipost ) */

    /* ---------------------------------------------------------------------- */
    /* clear the Flag array */
    /* ---------------------------------------------------------------------- */

    Common->mark = EMPTY ;
    cholmod_clear_flag (Common) ;
    ASSERT (cholmod_dump_work (TRUE, TRUE, 0, Common)) ;

    /* ---------------------------------------------------------------------- */
    /* find the permutation */
    /* ---------------------------------------------------------------------- */

    if (nd_camd)
    {

	/* ------------------------------------------------------------------ */
	/* apply csymamd or ccolamd using the Cmember constraints */
	/* ------------------------------------------------------------------ */

	if (symmetric)
	{
	    /* ordering A+A', so fset and fsize are ignored.
	     * Add the upper/lower part to a symmetric lower/upper matrix by
	     * converting to unsymmetric mode
	     * workspace: Iwork (max (nrow,ncol)) */
	    B = cholmod_copy (A, 0, -1, Common) ;
	    if (Common->status < CHOLMOD_OK)
	    {
		ASSERT (Common->malloc_count == orig) ;
		PRINT0 (("make symmetric failed\n")) ;
		return (EMPTY) ;
	    }
	    ASSERT ((int) (B->nrow) == n && (int) (B->ncol) == n) ;
	    PRINT1 (("nested dissection (2)\n")) ;
	    B->stype = -1 ;
	    /* workspace:  Head (nrow+1), Iwork (nrow) if symmetric-upper */
	    ok = cholmod_csymamd (B, Cmember, Perm, Common) ;
	    cholmod_free_sparse (&B, Common) ;
	    if (!ok)
	    {
		/* csymamd failed */
		ASSERT (Common->malloc_count == orig) ;
		PRINT0 (("csymamd failed\n")) ;
		return (EMPTY) ;
	    }
	}
	else
	{
	    /* ordering A*A' or A(:,f)*A(:,f)' */
	    /* workspace: Iwork (nrow if no fset; MAX(nrow,ncol) if fset) */
	    if (!cholmod_ccolamd (A, fset, fsize, Cmember, Perm, Common))
	    {
		/* ccolamd failed */
		ASSERT (Common->malloc_count == orig) ;
		PRINT1 (("ccolamd failed\n")) ;
		return (EMPTY) ;
	    }
	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* natural ordering of each component */
	/* ------------------------------------------------------------------ */

	/* use Iwork [0..n-1] for Next [ */
	Next = Iwork  ;

	/* ------------------------------------------------------------------ */
	/* place the nodes in link lists, one list per component */
	/* ------------------------------------------------------------------ */

	/* do so in reverse order, to preserve original ordering */
	for (j = n-1 ; j >= 0 ; j--)
	{
	    /* node j is in the cth component */
	    c = Cmember [j] ;
	    ASSERT (c >= 0 && c < ncomponents) ;
	    /* place node j in link list for component c */
	    Next [j] = Head [c] ;
	    Head [c] = j ;
	}

	/* ------------------------------------------------------------------ */
	/* order each node in each component */
	/* ------------------------------------------------------------------ */

	k = 0 ;
	for (c = 0 ; c < ncomponents ; c++)
	{
	    for (j = Head [c] ; j != EMPTY ; j = Next [j])
	    {
		Perm [k++] = j ;
	    }
	    Head [c] = EMPTY ;
	}
	ASSERT (k == n) ;

	/* done using Iwork [0..n-1] for Next ] */
    }

    /* ---------------------------------------------------------------------- */
    /* clear workspace and return number of components */
    /* ---------------------------------------------------------------------- */

    ASSERT (Common->malloc_count == orig) ;
    ASSERT (cholmod_dump_work (TRUE, TRUE, 0, Common)) ;
    return (ncomponents) ;
}
