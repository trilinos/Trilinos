/* ========================================================================== */
/* === paraklete_analyze ==================================================== */
/* ========================================================================== */

#include "paraklete.h"

/* LUsymbolic = paraklete_analyze (A, Common) finds a fill-reducing permutation
 * of A and its separator tree, and returns a paraklete_symbolic object
 * containing the symbolic analysis.
 *
 * TODO: check return values of MPI
 *
 * PARAKLETE version 0.1: parallel sparse LU factorization.  May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* ========================================================================== */
/* === paraklete_bcast_symbolic ============================================= */
/* ========================================================================== */

/* Process 0 has computed the symbolic object; broadcast it to all processes.
 * Returns TRUE and a non-NULL LUsymbolic object if successful.  Otherwise,
 * returns FALSE and a NULL LUsymbolic object.  This routine is not used for
 * the sequential case.
 *
 * If the symbolic analysis fails, all processes receive an object with
 * n=-1, which denotes a failure.  All processes then return NULL.
 */

#ifndef NMPI

static int paraklete_bcast_symbolic
(
    paraklete_symbolic **LUsymbolicHandle,
    paraklete_common *Common
)

{
    paraklete_symbolic *LUsymbolic ;
    int n, ncomponents, ok, all_ok, header [2] ;
    cholmod_common *cm ;

    cm = &(Common->cm) ;
    n = EMPTY ;
    ncomponents = EMPTY ;

    if (Common->myid == 0)
    {
	LUsymbolic = *LUsymbolicHandle ;
	if (LUsymbolic != NULL)
	{
	    n = LUsymbolic->n ;
	    ncomponents = LUsymbolic->ncomponents ;
	}
    }
    else
    {
	/* other processes do not yet have the symbolic object */
	*LUsymbolicHandle = NULL ;
    }

    /* broadcast the size of the object, or -1 if a failure occured */
    header [0] = n ;
    header [1] = ncomponents ;
    MPI_Bcast (&header, 2, MPI_INT, TAG0, Common->mpicomm) ;
    n = header [0] ;
    ncomponents = header [1] ;
    if (n == EMPTY)
    {
	/* the analysis in the root process failed */
	PR0 ((Common->file, "proc %d root analyze fails\n", Common->myid)) ;
	return (FALSE) ;
    }

    PR1 ((Common->file, "proc %d in bcast symbolic: status %d header %d %d\n",
	    Common->myid, cm->status, header [0], header [1])) ;

    if (Common->myid != 0)
    {
	LUsymbolic = cholmod_malloc (1, sizeof (paraklete_symbolic), cm) ;
	if (LUsymbolic != NULL)
	{
	    LUsymbolic->Mem_n  = cholmod_malloc (3*n, sizeof (int), cm) ;
	    LUsymbolic->Cperm  = LUsymbolic->Mem_n ;
	    LUsymbolic->Cinv   = LUsymbolic->Mem_n + n ;
	    LUsymbolic->Cparent= LUsymbolic->Mem_n + 2*n ;

	    LUsymbolic->Mem_c  = cholmod_malloc (7*ncomponents+2,
				    sizeof(int), cm) ;

	    /* each of size ncomponents: */
	    LUsymbolic->Child  = LUsymbolic->Mem_c ;
	    LUsymbolic->Clnz   = LUsymbolic->Mem_c + ncomponents ;
	    LUsymbolic->Cn     = LUsymbolic->Mem_c + 2*ncomponents ;
	    LUsymbolic->Cnz    = LUsymbolic->Mem_c + 3*ncomponents ;
	    LUsymbolic->Sched  = LUsymbolic->Mem_c + 4*ncomponents ;

	    /* each of size ncomponents+1: */
	    LUsymbolic->Cstart = LUsymbolic->Mem_c + 5*ncomponents ;
	    LUsymbolic->Childp = LUsymbolic->Mem_c + 6*ncomponents + 1 ;

	    LUsymbolic->n = n ;
	    LUsymbolic->ncomponents = ncomponents ;
	}
	*LUsymbolicHandle = LUsymbolic ;
    }

    ok = (cm->status == CHOLMOD_OK) ;

    /* all processes find out if any one process fails to allocate memory */
    MPI_Allreduce (&ok, &all_ok, 1, MPI_INT, MPI_LAND, Common->mpicomm) ;
    if (!all_ok)
    {
	/* out of memory; inform all processes */
	PR0 ((Common->file, "proc %d all fail in analyze\n", Common->myid)) ;
	paraklete_free_symbolic (&LUsymbolic, Common) ;
	*LUsymbolicHandle = NULL ;
	return (FALSE) ;
    }

    /* broadcast the contents of the symbolic object */
    MPI_Bcast (LUsymbolic->Mem_n, 3*n, MPI_INT, TAG0, Common->mpicomm) ;
    MPI_Bcast (LUsymbolic->Mem_c, 7*ncomponents+2, MPI_INT, TAG0, Common->mpicomm) ;

    return (TRUE) ;
}
#endif


/* ========================================================================== */
/* === paraklete_analyze ==================================================== */
/* ========================================================================== */

paraklete_symbolic *paraklete_analyze
(
    /* matrix to analyze */ 
    cholmod_sparse *A,
    paraklete_common *Common
)
{
    double work, cnt ;
    double *Cwork ;
    cholmod_common *cm ;
    paraklete_symbolic *LUsymbolic ;
    cholmod_sparse *C, *AT, *Elo, *Eup ;
    int *Cperm, *Cinv, *Cparent, *Cmember, *ColCount, *Cstart, *Childp, *Clnz,
	*Cn, *Cnz, *Parent, *Post, *First, *Level, *Child, *Ap, *Ai, *Sched,
	*Leaves, *Merged, *Nchildren, *NewNode, *Cparent2, *W ;
    double one [2] = {1,1} ;
    int p, k, n, ncomponents, ci, cj, i, j, clast, c, parent, nparent, nroots,
	nproc, nleaves, cp, cmerge, cc, c2, ncomp2, parent2 ;

#ifndef HACK
    int proc, nchild ;
#endif

    /* ---------------------------------------------------------------------- */
    /* all processes except process 0 get the symbolic object from process 0 */
    /* ---------------------------------------------------------------------- */

    LUsymbolic = NULL ;
    if (Common->myid != 0)
    {
	MPI (paraklete_bcast_symbolic (&LUsymbolic, Common)) ;
	return (LUsymbolic) ;
    }

    /* ---------------------------------------------------------------------- */
    /* process 0 does the analysis */
    /* ---------------------------------------------------------------------- */

    nproc = Common->nproc ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

#if 0
    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    if (A->nrow != A->ncol || A->stype)
    {
	cholmod_error (CHOLMOD_INVALID, "paraklete: invalid matrix", cm) ;
	return (NULL) ;
    }
#endif

    n = A->nrow ;
    C = NULL ;
    W = NULL ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    cm = &(Common->cm) ;
    cholmod_allocate_work (n, 2*n, n, cm) ;
    if (cm->status != CHOLMOD_OK)
    {
	/* out of memory; inform all processes */
	MPI (paraklete_bcast_symbolic (&LUsymbolic, Common)) ;
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* allocate first part of symbolic factor */
    /* ---------------------------------------------------------------------- */

    LUsymbolic = cholmod_malloc (1, sizeof (paraklete_symbolic), cm) ;

    if (cm->status != CHOLMOD_OK)
    {
	/* out of memory; inform all processes */
	MPI (paraklete_bcast_symbolic (&LUsymbolic, Common)) ;
	return (NULL) ;
    }

    LUsymbolic->Mem_n = cholmod_malloc (3*n, sizeof (int), cm) ;
    Cperm      = LUsymbolic->Mem_n ;		/* size n */
    Cinv       = LUsymbolic->Mem_n + n ;	/* size n */
    Cparent    = LUsymbolic->Mem_n + 2*n ;	/* size n */

    LUsymbolic->n = n ;
    LUsymbolic->Cperm = Cperm ;
    LUsymbolic->Cinv = Cinv ;
    LUsymbolic->Cparent = Cparent ;

    LUsymbolic->Mem_c = NULL ;
    LUsymbolic->Cstart = NULL ;
    LUsymbolic->Child = NULL ;
    LUsymbolic->Childp = NULL ;
    LUsymbolic->Clnz = NULL ;
    LUsymbolic->Cn = NULL ;
    LUsymbolic->Cnz = NULL ;
    LUsymbolic->Sched = NULL ;

    LUsymbolic->ncomponents = 0 ;

    if (cm->status != CHOLMOD_OK)
    {
	/* out of memory; inform all processes */
	PR0 ((Common->file, "oops, proc 0 ran out\n")) ;
	paraklete_free_symbolic (&LUsymbolic, Common) ;
	MPI (paraklete_bcast_symbolic (&LUsymbolic, Common)) ;
	return (NULL) ;
    }

    Cmember = Cinv ;	    /* use Cinv as workspace for Cmember */

    /* ---------------------------------------------------------------------- */
    /* C = pattern of triu (A+A'), in symmetric/upper form */
    /* ---------------------------------------------------------------------- */

    AT = cholmod_transpose (A, FALSE, cm) ;
    C = cholmod_add (A, AT, one, one, FALSE, FALSE, cm) ;
    cholmod_free_sparse (&AT, cm) ;
    cholmod_band_inplace (0, n, 0, C, cm) ;

    if (cm->status != CHOLMOD_OK)
    {
	/* out of memory; inform all processes */
	PR0 ((Common->file, "oops, proc 0 ran out here2\n")) ;
	cholmod_free_sparse (&C, cm) ;
	paraklete_free_symbolic (&LUsymbolic, Common) ;
	MPI (paraklete_bcast_symbolic (&LUsymbolic, Common)) ;
	return (NULL) ;
    }

    C->stype = 1 ;

    /* ---------------------------------------------------------------------- */
    /* fill-reducing nested dissection ordering of C */
    /* ---------------------------------------------------------------------- */

    /* Cperm [k] = i if row/col i of A is the kth row/col of A(p,p)
     * ncomponents = # of components in separator tree
     * Cparent [c] is parent of c in separator tree, or EMPTY if c is a root
     * Cmember [i] = c if row/col i of A is in component c
     */

    cm->current = 0 ;
    cm->method [0].nd_small = 4 ;
    ncomponents = cholmod_nested_dissection (C, NULL, 0, Cperm, Cparent,
	    Cmember, cm) ;

    if (cm->status != CHOLMOD_OK)
    {
	/* out of memory; inform all processes */
	PR0 ((Common->file, "oops, proc 0 ran out here3\n")) ;
	paraklete_free_symbolic (&LUsymbolic, Common) ;
	cholmod_free_sparse (&C, cm) ;
	MPI (paraklete_bcast_symbolic (&LUsymbolic, Common)) ;
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* Elo = C (p,p)', Eup = Elo' */
    /* ---------------------------------------------------------------------- */

    Elo = cholmod_ptranspose (C, FALSE, Cperm, NULL, 0, cm) ;
    cholmod_free_sparse (&C, cm) ;
    Eup = cholmod_transpose (Elo, FALSE, cm) ;

    /* ---------------------------------------------------------------------- */
    /* allocate more workspace */
    /* ---------------------------------------------------------------------- */

    W = cholmod_malloc (5*(n+1), sizeof (int), cm) ;
    ColCount = W ;	    /* size n [ */
    Parent   = W + n ;	    /* size n [ */
    Post     = W + 2*n ;    /* size n [ */
    First    = W + 3*n ;    /* size n [ */
    Level    = W + 4*n ;    /* size n [ */

    if (cm->status != CHOLMOD_OK)
    {
	/* out of memory; inform all processes */
	PR0 ((Common->file, "oops, proc 0 ran out here4\n")) ;
	paraklete_free_symbolic (&LUsymbolic, Common) ;
	cholmod_free_sparse (&Eup, cm) ;
	cholmod_free_sparse (&Elo, cm) ;
	cholmod_free (5*(n+1), sizeof (int), W, cm) ;
	MPI (paraklete_bcast_symbolic (&LUsymbolic, Common)) ;
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get an estimate of the # entries in L and U */
    /* ---------------------------------------------------------------------- */

    /* This assumes an LU factorization of C, with no partial pivoting */
    cholmod_etree (Eup, Parent, cm) ;
    cholmod_postorder (Parent, n, NULL, Post, cm) ;
    cholmod_rowcolcounts (Elo, NULL, 0, Parent, Post, NULL, ColCount,
	    First, Level, cm) ;

    cholmod_free_sparse (&Eup, cm) ;
    cholmod_free_sparse (&Elo, cm) ;

    /* Parent, Post, First, Level no longer needed ]]]] */

    if (cm->status != CHOLMOD_OK)
    {
	/* out of memory or other error; inform all processes */
	PR0 ((Common->file, "oops, proc 0 ran out here5\n")) ;
	paraklete_free_symbolic (&LUsymbolic, Common) ;
	cholmod_free (5*(n+1), sizeof (int), W, cm) ;
	MPI (paraklete_bcast_symbolic (&LUsymbolic, Common)) ;
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* compute Cwork [c] = flops to be done at component c */
    /* ---------------------------------------------------------------------- */

    Cwork = cm->Xwork ;
    for (c = 0 ; c < ncomponents ; c++)
    {
	Cwork [c] = 0 ;
    }
    for (k = 0 ; k < n ; k++)
    {
	c = Cmember [Cperm [k]] ;
	cnt = ColCount [k] ;
	Cwork [c] += cnt * cnt ;
    }

#ifndef NDEBUG
    for (c = 0 ; c < ncomponents ; c++)
    {
	PR1 ((Common->file, "Node %d work %g parent %d\n",
		    c, Cwork [c], Cparent [c])) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* compress the tree until it has <= nproc leaves */
    /* ---------------------------------------------------------------------- */

    /* Nchildren [c] = the number of children of node c */
    /* Merged [c] = node that c is merged into, or EMPTY if c is not merged */
    /* Leaves [0..nleaves-1] is a list of the leaves of the current tree */

    /* Note that W [0..n-1] is still in use for ColCount [0..n-1] */
    Merged = W + n ;				/* size ncomponents+1 [ */
    Nchildren = W + n + ncomponents + 1 ;	/* size ncomponents+1 [ */
    Leaves = W + n + 2 * (ncomponents+1) ;	/* size ncomponents [ */

    for (c = 0 ; c <= ncomponents ; c++)
    {
	Nchildren [c] = 0 ;
	Merged [c] = EMPTY ;
    }
    for (c = 0 ; c < ncomponents ; c++)
    {
	parent = Cparent [c] ;
	if (parent == EMPTY)
	{
	    parent = ncomponents ;
	}
	ASSERT (parent > c && parent <= ncomponents) ;
	Nchildren [parent]++ ;
    }

    /* make a list of all leaves */
    nleaves = 0 ;
    for (c = 0 ; c < ncomponents ; c++)
    {
	PR1 ((Common->file, "Node %d has %d children\n", c, Nchildren [c])) ;
	if (Nchildren [c] == 0)
	{
	    PR1 ((Common->file, "Leaf: %d\n", c)) ;
	    Leaves [nleaves++] = c ;
	}
    }

    /* TODO find out why cholmod_nested_dissection returns a graph with
     * 2 nodes (one parent and one child) for vanHeukelum/cage3
     *
     * TODO: use a heap for the leaves
     */

    while (nleaves > nproc)
    {
	PR1 ((Common->file, "\n------------ nleaves: %d\n", nleaves)) ;

	/* find the lightest leaf (skip node ncomponents-1) */
	work = EMPTY ;
	c = EMPTY ;
	cp = EMPTY ;
	for (p = 0 ; p < nleaves ; p++)
	{
	    PR2 ((Common->file, "Leaf %d work %g\n",
			Leaves [p], Cwork [Leaves [p]])) ;
	    ASSERT (Merged [Leaves [p]] == EMPTY) ;
	    if (Leaves [p] == ncomponents-1)
	    {
		/* node ncomponents-1 has no live node to its right (that is,
		 * a higher-numbered node), so skip it.  The alternative is to
		 * merge a node cmerge to its left into node ncomponents-1.
		 * This may lead to a better balance of work, but is more
		 * complicated.  The parent of the children of cmerge would
		 * have to be updated.  FUTURE WORK: consider handling this
		 * case. */
		continue ;
	    }
	    if (work == EMPTY || Cwork [Leaves [p]] < work)
	    {
		c = Leaves [p] ;
		work = Cwork [c] ;
		cp = p ;
	    }
	}
	ASSERT (c != EMPTY) ;
	PR2 ((Common->file, "Lightest leaf is %d with work %g\n", c, Cwork [c]));
	ASSERT (c < ncomponents-1) ;
	ASSERT (Nchildren [c] == 0) ;

	/* find the live node to the right of this node */
	for (cmerge = c+1 ; cmerge < ncomponents ; cmerge++)
	{
	    if (Merged [cmerge] == EMPTY)
	    {
		break ;
	    }
	}

	/* find the parent of c */
	parent = Cparent [c] ;
	if (parent == EMPTY)
	{
	    parent = ncomponents ;
	}

	/* merge c into cmerge node, where c is a leaf */
	PR1 ((Common->file, "merge %d into %d, parent %d\n", c, cmerge, parent));
	ASSERT (cmerge < ncomponents) ;
	Cwork [cmerge] += Cwork [c] ;
	Cwork [c] = 0 ;
	Merged [c] = cmerge ;
	Leaves [cp] = Leaves [--nleaves] ;
	Nchildren [parent]-- ;
	ASSERT (Merged [parent] == EMPTY) ;

	if (Nchildren [parent] == 0 && parent != ncomponents)
	{
	    /* parent is a new leaf, add it to the list of Leaves */
	    PR1 ((Common->file, "parent is new leaf: %d\n", parent)) ;
	    Leaves [nleaves++] = parent ;
	}
    }

    /* Leaves no longer needed ] */

    PR1 ((Common->file, "\n--------------------------- done merging leaves\n")) ;

    /* merge nodes that have just one child, with the one child */
    for (c = 0 ; c < ncomponents-1 ; c++)
    {
	if (Merged [c] == EMPTY)
	{
	    parent = Cparent [c] ;
	    if (parent == EMPTY) continue ;
	    if (Nchildren [parent] == 1)
	    {
		PR1 ((Common->file, "\nparent %d of c %d has one child\n",
			    parent, c));
		Cwork [parent] += Cwork [c] ;
		Cwork [c] = 0 ;
		Merged [c] = parent ;
		for (cc = c+1 ; cc < parent ; cc++)
		{
		    PR1 ((Common->file, "merge %d into %d\n", cc, parent)) ;
		    ASSERT (Merged [cc] != EMPTY) ;
		    Merged [cc] = parent ;
		}
	    }
	}
    }

    /* Nchildren no longer needed ] */

    /* compress the paths in Merged */
    for (c = 0 ; c < ncomponents ; c++)
    {
	/* find the ultimate node that node c was merged into */
	PR1 ((Common->file, "\nFind ultimate node for %d\n", c)) ;
	for (cc = c ; Merged [cc] != EMPTY ; cc = Merged [cc]) ;
	for (c2 = c ; Merged [c2] != EMPTY ; c2 = Merged [c2])
	{
	    PR1 ((Common->file, "   merge %d into %d\n", c2, cc)) ;
	    Merged [c2] = cc ;
	}
    }

    /* find the new node numbering, using Leaves as workspace */
    NewNode = Leaves ;
    ncomp2 = 0 ;
    for (c = 0 ; c < ncomponents ; c++)
    {
	if (Merged [c] == EMPTY)
	{
	    PR1 ((Common->file, "Live node %d becomes node %d\n", c, ncomp2)) ;
	    NewNode [c] = ncomp2++ ;
	}
    }
    for (c = 0 ; c < ncomponents ; c++)
    {
	if (Merged [c] != EMPTY)
	{
	    NewNode [c] = NewNode [Merged [c]] ;
	    PR1 ((Common->file, "Dead node %d becomes part of node %d\n",
			c, NewNode [c])) ;
	}
    }

    /* fix Cmember */
    for (k = 0 ; k < n ; k++)
    {
	c = Cmember [k] ;
	c = NewNode [c] ;
	Cmember [k] = c ;
    }

    /* fix Cparent, using Nchildren as workspace */
    Cparent2 = Nchildren ;
    for (c = 0 ; c < ncomponents ; c++)
    {
	if (Merged [c] == EMPTY)
	{
	    c2 = NewNode [c] ;
	    parent = Cparent [c] ;
	    parent2 = (parent == EMPTY) ? EMPTY : (NewNode [parent]) ;
	    Cparent2 [c2] = parent2 ;
	}
    }

    /* Merged no longer needed ] */

    for (c = 0 ; c < ncomponents ; c++)
    {
	Cwork [c] = 0 ;
    }

    ncomponents = ncomp2 ;
    for (c = 0 ; c < ncomponents ; c++)
    {
	Cparent [c] = Cparent2 [c] ;
    }

    for (c = 0 ; c < ncomponents ; c++)
    {
	PR1 ((Common->file, "New node: %d new parent %d\n", c, Cparent [c])) ;
    }

    /* ---------------------------------------------------------------------- */
    /* allocate remainder of symbolic factor */
    /* ---------------------------------------------------------------------- */

    LUsymbolic->Mem_c = cholmod_malloc (7*ncomponents+2, sizeof (int), cm) ;

    /* each of size ncomponents: */
    Child = LUsymbolic->Child  = LUsymbolic->Mem_c ;
    Clnz  = LUsymbolic->Clnz   = LUsymbolic->Mem_c + ncomponents ;
    Cn    = LUsymbolic->Cn     = LUsymbolic->Mem_c + 2*ncomponents ;
    Cnz   = LUsymbolic->Cnz    = LUsymbolic->Mem_c + 3*ncomponents ;
    Sched = LUsymbolic->Sched  = LUsymbolic->Mem_c + 4*ncomponents ;

    /* each of size ncomponents+1: */
    Cstart = LUsymbolic->Cstart = LUsymbolic->Mem_c + 5*ncomponents ;
    Childp = LUsymbolic->Childp = LUsymbolic->Mem_c + 6*ncomponents + 1 ;

    LUsymbolic->ncomponents = ncomponents ;

    if (cm->status != CHOLMOD_OK)
    {
	/* out of memory or other error; inform all processes */
	paraklete_free_symbolic (&LUsymbolic, Common) ;
	cholmod_free (5*(n+1), sizeof (int), W, cm) ;
	MPI (paraklete_bcast_symbolic (&LUsymbolic, Common)) ;
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* Cstart = start and end nodes of each component */
    /* ---------------------------------------------------------------------- */

    clast = EMPTY ;
    for (k = 0 ; k < n ; k++)
    {
	c = Cmember [Cperm [k]] ;
	if (c != clast)
	{
	    ASSERT (c == clast + 1) ;
	    Cstart [c] = k ;
	}
	clast = c ;
    }
    Cstart [ncomponents] = n ;

    /* ---------------------------------------------------------------------- */
    /* Clnz = estimate of # of entries in L for each component */
    /* ---------------------------------------------------------------------- */

    for (c = 0 ; c < ncomponents ; c++)
    {
	Clnz [c] = 0 ;
	for (k = Cstart [c] ; k < Cstart [c+1] ; k++)
	{
	    Clnz [c] += ColCount [k] ;
	}
    }

    /* ColCount no longer needed ] */

    /* ---------------------------------------------------------------------- */
    /* Child, Childp: list of children of each component */
    /* ---------------------------------------------------------------------- */

    /* count the number of children of each node, using Cn as workspace */
    for (c = 0 ; c < ncomponents ; c++)
    {
	Cn [c] = 0 ;
    }
    nroots = 0 ;
    for (c = 0 ; c < ncomponents ; c++)
    {
	parent = Cparent [c] ;
	PR1 ((Common->file, "node %d: parent %d\n", c, parent)) ;
	if (parent == EMPTY)
	{
	    nroots++ ;
	}
	else
	{
	    ASSERT (parent > 0 && parent < ncomponents) ;
	    Cn [parent]++ ;
	}
    }

#ifndef NDEBUG
    for (c = 0 ; c < ncomponents ; c++)
    {
	PR1 ((Common->file, "node %d: %d children\n", c, Cn [c])) ;
    }
#endif

    /* find the cumulative sum of the children */
    k = 0 ;
    for (c = 0 ; c < ncomponents ; c++)
    {
	Childp [c] = k ;
	k += Cn [c] ;
    }
    PR1 ((Common->file, "k %d ncomponents %d\n", k, ncomponents)) ;
    ASSERT (k == ncomponents - nroots) ;
    Childp [ncomponents] = k ;

    /* create a list of children for each node */
    for (c = 0 ; c < ncomponents ; c++)
    {
	Cn [c] = Childp [c] ;
    }
    for (c = 0 ; c < ncomponents ; c++)
    {
	parent = Cparent [c] ;
	if (parent != EMPTY)
	{
	    Child [Cn [parent]++] = c ;
	}
    }

#ifndef NDEBUG
    for (c = 0 ; c < ncomponents ; c++)
    {
	PR1 ((Common->file, "Node %d children: ", c)) ;
	for (cp = Childp [c] ; cp < Childp [c+1] ; cp++)
	{
	    PR1 ((Common->file, "%d ", Child [cp])) ;
	}
	PR1 ((Common->file, "\n")) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* find the nominal dimensions of each node (assuming no pivot delays) */
    /* ---------------------------------------------------------------------- */

    for (c = ncomponents - 1 ; c >= 0 ; c--)
    {
	parent = Cparent [c] ;
	nparent = (parent == EMPTY) ? 0 : Cn [parent] ;
	Cn [c] = (Cstart [c+1] - Cstart [c]) + nparent ;
	PR1 ((Common->file, "node %d Cn: %d\n", c, Cn [c])) ;
    }

    /* ---------------------------------------------------------------------- */
    /* count the nonzeros in A that map to node c */
    /* ---------------------------------------------------------------------- */

    for (c = 0 ; c < ncomponents ; c++)
    {
	Cnz [c] = 0 ;
    }
    Ap = A->p ;
    Ai = A->i ;
    for (j = 0 ; j < n ; j++)
    {
	cj = Cmember [j] ;
	for (p = Ap [j] ; p < Ap [j+1] ; p++)
	{
	    i = Ai [p] ;
	    ci = Cmember [i] ;
	    c = MIN (ci, cj) ;
	    Cnz [c]++ ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* compute Cinv = inverse of Cperm */
    /* ---------------------------------------------------------------------- */

    for (k = 0 ; k < n ; k++)
    {
	Cinv [Cperm [k]] = k ;
    }

    /* ---------------------------------------------------------------------- */
    /* schedule the processes to the nodes */
    /* ---------------------------------------------------------------------- */

#ifdef HACK
    for (c = 0 ; c < ncomponents ; c++)
    {
	Sched [c] = 0 ;
    }
#else
    for (c = 0 ; c < ncomponents ; c++)
    {
	Sched [c] = EMPTY ;
    }
    proc = 0 ;
    for (c = 0 ; c < ncomponents ; c++)
    {
	nchild = Childp [c+1] - Childp [c] ;
	if (nchild == 0)
	{
	    PR1 ((Common->file, "\nSchedule child %d to process %d\n", c, proc));
	    ASSERT (proc < nproc) ;
	    for (cc = c ; cc != EMPTY && Sched [cc] == EMPTY ; cc = Cparent[cc])
	    {
		PR1 ((Common->file, "  node %d to process %d\n", cc, proc)) ;
		Sched [cc] = proc ;
	    }
	    proc++ ;
	}
    }
#endif

#ifndef NDEBUG
    PR0 ((Common->file, "\nncomponents:: %d\n", ncomponents)) ;
    for (c = 0 ; c < ncomponents ; c++)
    {
	PR0 ((Common->file, "    node %d Sched %d : Cparent %d proc %d\n",
		    c, Sched [c], Cparent [c],
		    (Cparent [c] == EMPTY) ? EMPTY : Sched [Cparent [c]])) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* return the symbolic factorization, and broadcast it to all processes */
    /* ---------------------------------------------------------------------- */

    cholmod_free (5*(n+1), sizeof (int), W, cm) ;
    PR1 ((Common->file, "analysis done\n")) ;
    MPI (paraklete_bcast_symbolic (&LUsymbolic, Common)) ;
    return (LUsymbolic) ;
}


/* ========================================================================== */
/* === paraklete_free_symbolic ============================================== */
/* ========================================================================== */

/* Free the symbolic object.  All processes own a copy after it's broadcast. */

void paraklete_free_symbolic
(
    paraklete_symbolic **LUsymbolicHandle,
    paraklete_common *Common
)
{
    paraklete_symbolic *LUsymbolic ;
    cholmod_common *cm ;
    int n, ncomponents ;

    if (LUsymbolicHandle == NULL)
    {
	/* nothing to do */
	return ;
    }
    LUsymbolic = *LUsymbolicHandle ;
    if (LUsymbolic == NULL)
    {
	/* nothing to do */
	return ;
    }

    cm = &(Common->cm) ;
    ncomponents = LUsymbolic->ncomponents ;
    n = LUsymbolic->n ;

    cholmod_free (3*n, sizeof (int), LUsymbolic->Mem_n, cm) ;
    LUsymbolic->Cperm = NULL ; 
    LUsymbolic->Cinv = NULL ; 
    LUsymbolic->Cparent = NULL ; 

    cholmod_free (7*ncomponents + 2, sizeof (int), LUsymbolic->Mem_c, cm) ;
    LUsymbolic->Cstart = NULL ;
    LUsymbolic->Child = NULL ;
    LUsymbolic->Childp = NULL ;
    LUsymbolic->Clnz = NULL ;
    LUsymbolic->Cn = NULL ;
    LUsymbolic->Cnz = NULL ;
    LUsymbolic->Sched = NULL ;

    *LUsymbolicHandle = cholmod_free (
	    1, sizeof (paraklete_symbolic), (*LUsymbolicHandle), cm) ;
}
