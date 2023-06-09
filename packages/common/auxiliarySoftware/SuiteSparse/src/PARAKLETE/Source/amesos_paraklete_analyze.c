/* ========================================================================== */
/* === paraklete_analyze ==================================================== */
/* ========================================================================== */

#include "amesos_paraklete_decl.h"

/* LUsymbolic = paraklete_analyze (A, Common) finds a fill-reducing permutation
 * of A and its separator tree, and returns a paraklete_symbolic object
 * containing the symbolic analysis.
 *
 * TODO: check return values of MPI
 *
 * PARAKLETE version 0.3: parallel sparse LU factorization.  Nov 13, 2007
 * Copyright (C) 2007, Univ. of Florida.  Author: Timothy A. Davis
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

static Int paraklete_bcast_symbolic
(
    paraklete_symbolic **LUsymbolicHandle,
    paraklete_common *Common
)

{
    paraklete_symbolic *LUsymbolic = NULL ;
    Int n, ncomponents, header [2] ;
    int ok, all_ok ;
    cholmod_common *cm ;

    cm = &(Common->cm) ;
    n = TRILINOS_CHOLMOD_EMPTY ;
    ncomponents = TRILINOS_CHOLMOD_EMPTY ;

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
    MPI_Bcast (&header, 2, MPI_Int, TAG0, MPI_COMM_WORLD) ;
    n = header [0] ;
    ncomponents = header [1] ;
    if (n == TRILINOS_CHOLMOD_EMPTY)
    {
	/* the analysis in the root process failed */
	PR0 ((Common->file, "proc "ID" root analyze fails\n", Common->myid)) ;
	return (FALSE) ;
    }

    PR1 ((Common->file, "proc "ID" in bcast symbolic: status "ID" header "ID" "ID"\n",
	    Common->myid, cm->status, header [0], header [1])) ;

    if (Common->myid != 0)
    {
	LUsymbolic = amesos_paraklete_alloc_symbolic (n, ncomponents, FALSE, Common) ;
	*LUsymbolicHandle = LUsymbolic ;
    }

    ok = (cm->status == CHOLMOD_OK) && (LUsymbolic != NULL) ;
    all_ok = ok ;

    /* all processes find out if any one process fails to allocate memory */
    MPI_Allreduce (&ok, &all_ok, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD) ;
    if (!all_ok)
    {
	/* out of memory; inform all processes */
	PR0 ((Common->file, "proc "ID" all fail in analyze\n", Common->myid)) ;
	amesos_paraklete_free_symbolic (&LUsymbolic, Common) ;
	*LUsymbolicHandle = NULL ;
	return (FALSE) ;
    }

    /* broadcast the contents of the symbolic object */
    MPI_Bcast (LUsymbolic->Mem_n, 3*n, MPI_Int, TAG0, MPI_COMM_WORLD) ;
    MPI_Bcast (LUsymbolic->Mem_c, 7*ncomponents+2, MPI_Int, TAG0,
	MPI_COMM_WORLD) ;

#if 0
    {
	/* each of size ncomponents: */
	Int *Child = LUsymbolic->Child ;
	Int *Clnz  = LUsymbolic->Clnz ;
	Int *Cn    = LUsymbolic->Cn ;
	Int *Cnz   = LUsymbolic->Cnz ;
	Int *Sched = LUsymbolic->Sched ;

	/* each of size ncomponents+1: */
	Int *Cstart = LUsymbolic->Cstart ;
	Int *Childp = LUsymbolic->Childp ;
	Int cc ;

	for (cc = 0 ; cc < ncomponents ; cc++)
	{
	    printf ("component "ID"\n", cc) ;
	    printf ("Child "ID"\n", Child [cc]) ;
	    printf ("Clnz "ID"\n", Clnz [cc]) ;
	    printf ("Cn "ID"\n", Cn [cc]) ;
	    printf ("Cnz "ID"\n", Cnz [cc]) ;
	    printf ("Sched "ID"\n", Sched [cc]) ;
	    printf ("Cstart "ID"\n", Cstart [cc]) ;
	    printf ("Childp "ID"\n", Childp [cc]) ;
	}
	    printf ("Cstart "ID"\n", Cstart [ncomponents]) ;
	    printf ("Childp "ID"\n", Childp [ncomponents]) ;

    }
#endif

    return (TRUE) ;
}
#endif


/* ========================================================================== */
/* === paraklete_analyze ==================================================== */
/* ========================================================================== */

paraklete_symbolic *amesos_paraklete_analyze
(
    /* matrix to analyze */ 
    cholmod_sparse *A,
    paraklete_common *Common
)
{
    double cnt ;
    double *Cwork ;
    cholmod_common *cm ;
    paraklete_symbolic *LUsymbolic ;
    cholmod_sparse *C, *AT, *Elo, *Eup ;
    Int *Cperm, *RpermInv, *Cparent, *Cmember, *ColCount, *Cstart, *Childp,
	*Clnz, *Cn, *Cnz, *Parent, *Post, *First, *Level, *Child, *Ap, *Ai,
	*Sched, *W, *Rperm,
        *Lo_id, *Hi_id ;
    double one [2] = {1,1} ;
    Int p, k, n, ncomponents, ci, cj, i, j, clast, c, parent, nparent, nroots,
	nproc, cp, nc0, nc1 ;
    int ok = TRUE ;
    size_t n5, nc7 ;

#if 0
    double work ;
    Int proc = 0, parent2, ncomp2, c2, cc, cmerge, nleaves,
    Int *NewNode, *Cparent2, *Merged, *Leaves, *Nchildren ;
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
	CHOLMOD (error) (CHOLMOD_INVALID, "paraklete: invalid matrix", cm) ;
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
    CHOLMOD (allocate_work) (n, 2*n, n, cm) ;
    if (cm->status != CHOLMOD_OK)
    {
	/* out of memory; inform all processes */
	MPI (paraklete_bcast_symbolic (&LUsymbolic, Common)) ;
        /*
        printf ("   analyze failed 1!\n") ;
        */
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* allocate first part of symbolic factor */
    /* ---------------------------------------------------------------------- */

    LUsymbolic = amesos_paraklete_alloc_symbolic (n, 0, TRUE, Common) ;

    if (LUsymbolic == NULL)
    {
	/* out of memory; inform all processes */
	PR0 ((Common->file, "oops, proc 0 ran out\n")) ;
	MPI (paraklete_bcast_symbolic (&LUsymbolic, Common)) ;
        /*
        printf ("   analyze failed 2!\n") ;
        */
	return (NULL) ;
    }

    Cperm    = LUsymbolic->Cperm ;
    RpermInv = LUsymbolic->RpermInv ;
    Cparent  = LUsymbolic->Cparent ;

    Cmember = RpermInv ;	    /* use RpermInv as workspace for Cmember */

    /* ---------------------------------------------------------------------- */
    /* C = pattern of triu (A+A'), in symmetric/upper form */
    /* ---------------------------------------------------------------------- */

    /*
    printf ("pattern of A+A', n = "ID" ("ID")\n", A->nrow, cm->status) ;
    printf ("A "ID" by "ID", nzmax "ID"\n", A->nrow, A->ncol, A->nzmax) ;
    */
    AT = CHOLMOD (transpose) (A, FALSE, cm) ;
    /*
    printf ("AT is %p ("ID")\n", (void *) AT, cm->status) ;
    */
    C = CHOLMOD (add) (A, AT, one, one, FALSE, FALSE, cm) ;
    /*
    printf ("C is %p ("ID")\n", (void *) C, cm->status) ;
    */
    CHOLMOD (free_sparse) (&AT, cm) ;
    CHOLMOD (band_inplace) (0, n, 0, C, cm) ;
    /*
    printf ("status is ("ID")\n", cm->status) ;
    */

    if (cm->status != CHOLMOD_OK)
    {
	/* out of memory; inform all processes */
	PR0 ((Common->file, "oops, proc 0 ran out here2\n")) ;
	CHOLMOD (free_sparse) (&C, cm) ;
	amesos_paraklete_free_symbolic (&LUsymbolic, Common) ;
	MPI (paraklete_bcast_symbolic (&LUsymbolic, Common)) ;
        /*
        printf ("   analyze failed 3!\n") ;
        */
	return (NULL) ;
    }

    C->stype = 1 ;

    /* ---------------------------------------------------------------------- */
    /* fill-reducing nested dissection ordering of C */
    /* ---------------------------------------------------------------------- */

    /* Cperm [k] = i if row/col i of A is the kth row/col of A(p,p)
     * ncomponents = # of components in separator tree
     * Cparent [c] is parent of c in separator tree, or TRILINOS_CHOLMOD_EMPTY if c is a root
     * Cmember [i] = c if row/col i of A is in component c
     */

    /* TODO rename Common->nleaves to be something else */
    cm->method [0].nd_oksep = 0.1 ;

    if (Common->nleaves <= 0)
    {
        cm->method [0].nd_small = MAX (1000, -(Common->nleaves)) ;
    }
    else
    {
        cm->method [0].nd_small = n / Common->nleaves ;
    }

    cm->current = 0 ;
    cm->method [0].nd_components = 0 ;  /* default value */
    /*
    printf ("nd_components "ID"\n", cm->method [0].nd_components) ;
    */

    ncomponents = CHOLMOD (nested_dissection) (C, NULL, 0, Cperm, Cparent,
        Cmember, cm) ;

    nc0 = ncomponents ; /* from CHOLMOD (nested_dissection) */
    nc1 = ncomponents ; /* after collapsing */

#ifndef NDEBUG
    /* check results: */
    clast = TRILINOS_CHOLMOD_EMPTY ;
    for (k = 0 ; k < n ; k++)
    {
	c = Cmember [Cperm [k]] ;
	if (c != clast)
	{
	    /*
	    printf ("Cmember ["ID"] = "ID"\n", k, Cmember [Cperm [k]]) ;
	    printf ("start of component\n") ;
	    */
	    /* if (c != clast+1) { printf ("ERROR!\n") ; exit (1) ; } */
	    ASSERT (c == clast + 1) ;
	}
	clast = c ;
    }
#endif

    if (cm->status != CHOLMOD_OK)
    {
	/* out of memory; inform all processes */
	PR0 ((Common->file, "oops, proc 0 ran out here3\n")) ;
	amesos_paraklete_free_symbolic (&LUsymbolic, Common) ;
	CHOLMOD (free_sparse) (&C, cm) ;
	MPI (paraklete_bcast_symbolic (&LUsymbolic, Common)) ;
        /*
        printf ("   analyze failed 4!\n") ;
        */
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* Elo = C (p,p)', Eup = Elo' */
    /* ---------------------------------------------------------------------- */

    Elo = CHOLMOD (ptranspose) (C, FALSE, Cperm, NULL, 0, cm) ;
    CHOLMOD (free_sparse) (&C, cm) ;
    Eup = CHOLMOD (transpose) (Elo, FALSE, cm) ;

    /* ---------------------------------------------------------------------- */
    /* allocate more workspace */
    /* ---------------------------------------------------------------------- */

    /* n5 = 5*(n+1) */
    n5 = CHOLMOD (mult_size_t) (n+1, 5, &ok) ;
    if (!ok) PARAKLETE_ERROR (PK_TOO_LARGE, "problem too large") ;

    W = CHOLMOD (malloc) (n5, sizeof (Int), cm) ;
    ColCount = W ;	    /* size n [ */
    Parent   = W + n ;	    /* size n [ */
    Post     = W + 2*n ;    /* size n [ */
    First    = W + 3*n ;    /* size n [ */
    Level    = W + 4*n ;    /* size n [ */

    if (cm->status != CHOLMOD_OK)
    {
	/* out of memory; inform all processes */
	PR0 ((Common->file, "oops, proc 0 ran out here4\n")) ;
	amesos_paraklete_free_symbolic (&LUsymbolic, Common) ;
	CHOLMOD (free_sparse) (&Eup, cm) ;
	CHOLMOD (free_sparse) (&Elo, cm) ;
	CHOLMOD (free) (n5, sizeof (Int), W, cm) ;
	MPI (paraklete_bcast_symbolic (&LUsymbolic, Common)) ;
        PARAKLETE_ERROR (PK_OUT_OF_MEMORY, "out of memory") ;
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get an estimate of the # entries in L and U */
    /* ---------------------------------------------------------------------- */

    /* This assumes an LU factorization of C, with no partial pivoting */
    CHOLMOD (etree) (Eup, Parent, cm) ;
    CHOLMOD (postorder) (Parent, n, NULL, Post, cm) ;
    CHOLMOD (rowcolcounts) (Elo, NULL, 0, Parent, Post, NULL, ColCount,
	    First, Level, cm) ;

    CHOLMOD (free_sparse) (&Eup, cm) ;
    CHOLMOD (free_sparse) (&Elo, cm) ;

    /* Parent, Post, First, Level no longer needed ]]]] */

    if (cm->status != CHOLMOD_OK)
    {
	/* out of memory or other error; inform all processes */
        PARAKLETE_ERROR (PK_UNKNOWN, "out of memory or other error") ;
	amesos_paraklete_free_symbolic (&LUsymbolic, Common) ;
	CHOLMOD (free) (n5, sizeof (Int), W, cm) ;
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
	PR1 ((Common->file, "Node "ID" work %g parent "ID"\n",
		    c, Cwork [c], Cparent [c])) ;
    }
#endif

#if 0

    /* ---------------------------------------------------------------------- */
    /* compress the tree until it has <= nproc leaves */
    /* ---------------------------------------------------------------------- */

    /* Nchildren [c] = the number of children of node c */
    /* Merged [c] = node that c is merged into, or TRILINOS_CHOLMOD_EMPTY if c is not merged */
    /* Leaves [0..nleaves-1] is a list of the leaves of the current tree */

    /* Note that W [0..n-1] is still in use for ColCount [0..n-1] */
    Merged = W + n ;				/* size ncomponents+1 [ */
    Nchildren = W + n + ncomponents + 1 ;	/* size ncomponents+1 [ */
    Leaves = W + n + 2 * (ncomponents+1) ;	/* size ncomponents [ */

    for (c = 0 ; c <= ncomponents ; c++)
    {
	Nchildren [c] = 0 ;
	Merged [c] = TRILINOS_CHOLMOD_EMPTY ;
    }
    for (c = 0 ; c < ncomponents ; c++)
    {
	parent = Cparent [c] ;
	if (parent == TRILINOS_CHOLMOD_EMPTY)
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
	PR1 ((Common->file, "Node "ID" has "ID" children\n", c, Nchildren [c])) ;
	if (Nchildren [c] == 0)
	{
	    PR1 ((Common->file, "Leaf: "ID"\n", c)) ;
	    Leaves [nleaves++] = c ;
	}
    }

    /* CHOLMOD (nested_dissection) returns a graph with
     * 2 nodes (one parent and one child) for vanHeukelum/cage3
     *
     * could use a heap for the leaves
     */

    while (nleaves > target_nleaves)
    {
	PR1 ((Common->file, "\n------------ nleaves: "ID"\n", nleaves)) ;

	/* find the lightest leaf (skip node ncomponents-1) */
	work = TRILINOS_CHOLMOD_EMPTY ;
	c = TRILINOS_CHOLMOD_EMPTY ;
	cp = TRILINOS_CHOLMOD_EMPTY ;
	for (p = 0 ; p < nleaves ; p++)
	{
	    PR2 ((Common->file, "Leaf "ID" work %g\n",
			Leaves [p], Cwork [Leaves [p]])) ;
	    ASSERT (Merged [Leaves [p]] == TRILINOS_CHOLMOD_EMPTY) ;
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
	    if (work == TRILINOS_CHOLMOD_EMPTY || Cwork [Leaves [p]] < work)
	    {
		c = Leaves [p] ;
		work = Cwork [c] ;
		cp = p ;
	    }
	}
	ASSERT (c != TRILINOS_CHOLMOD_EMPTY) ;
	PR2 ((Common->file,"Lightest leaf is "ID" with work %g\n", c, Cwork [c]));
	ASSERT (c < ncomponents-1) ;
	ASSERT (Nchildren [c] == 0) ;

	/* find the live node to the right of this node */
	for (cmerge = c+1 ; cmerge < ncomponents ; cmerge++)
	{
	    if (Merged [cmerge] == TRILINOS_CHOLMOD_EMPTY)
	    {
		break ;
	    }
	}

	/* find the parent of c */
	parent = Cparent [c] ;
	if (parent == TRILINOS_CHOLMOD_EMPTY)
	{
	    parent = ncomponents ;
	}

	/* merge c into cmerge node, where c is a leaf */
	PR1 ((Common->file,"merge "ID" into "ID", parent "ID"\n", c, cmerge, parent));
	ASSERT (cmerge < ncomponents) ;
	Cwork [cmerge] += Cwork [c] ;
	Cwork [c] = 0 ;
	Merged [c] = cmerge ;
	Leaves [cp] = Leaves [--nleaves] ;
	Nchildren [parent]-- ;
	ASSERT (Merged [parent] == TRILINOS_CHOLMOD_EMPTY) ;

	if (Nchildren [parent] == 0 && parent != ncomponents)
	{
	    /* parent is a new leaf, add it to the list of Leaves */
	    PR1 ((Common->file, "parent is new leaf: "ID"\n", parent)) ;
	    Leaves [nleaves++] = parent ;
	}
    }

    /* Leaves no longer needed ] */

    PR1 ((Common->file, "\n-------------------------- done merging leaves\n")) ;

    /* merge nodes that have just one child, with the one child */
    for (c = 0 ; c < ncomponents-1 ; c++)
    {
	if (Merged [c] == TRILINOS_CHOLMOD_EMPTY)
	{
	    parent = Cparent [c] ;
	    if (parent == TRILINOS_CHOLMOD_EMPTY) continue ;
	    if (Nchildren [parent] == 1)
	    {
		PR1 ((Common->file, "\nparent "ID" of c "ID" has one child\n",
			    parent, c));
		Cwork [parent] += Cwork [c] ;
		Cwork [c] = 0 ;
		Merged [c] = parent ;
		for (cc = c+1 ; cc < parent ; cc++)
		{
		    PR1 ((Common->file, "merge "ID" into "ID"\n", cc, parent)) ;
		    ASSERT (Merged [cc] != TRILINOS_CHOLMOD_EMPTY) ;
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
	PR1 ((Common->file, "\nFind ultimate node for "ID"\n", c)) ;
	for (cc = c ; Merged [cc] != TRILINOS_CHOLMOD_EMPTY ; cc = Merged [cc]) ;
	for (c2 = c ; Merged [c2] != TRILINOS_CHOLMOD_EMPTY ; c2 = Merged [c2])
	{
	    PR1 ((Common->file, "   merge "ID" into "ID"\n", c2, cc)) ;
	    Merged [c2] = cc ;
	}
    }

    /* find the new node numbering, using Leaves as workspace */
    NewNode = Leaves ;
    ncomp2 = 0 ;
    for (c = 0 ; c < ncomponents ; c++)
    {
	if (Merged [c] == TRILINOS_CHOLMOD_EMPTY)
	{
	    PR1 ((Common->file, "Live node "ID" becomes node "ID"\n", c, ncomp2)) ;
	    NewNode [c] = ncomp2++ ;
	}
    }
    for (c = 0 ; c < ncomponents ; c++)
    {
	if (Merged [c] != TRILINOS_CHOLMOD_EMPTY)
	{
	    NewNode [c] = NewNode [Merged [c]] ;
	    PR1 ((Common->file, "Dead node "ID" becomes part of node "ID"\n",
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
	if (Merged [c] == TRILINOS_CHOLMOD_EMPTY)
	{
	    c2 = NewNode [c] ;
	    parent = Cparent [c] ;
	    parent2 = (parent == TRILINOS_CHOLMOD_EMPTY) ? TRILINOS_CHOLMOD_EMPTY : (NewNode [parent]) ;
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

#ifndef NDEBUG
    printf ("Final components: "ID" leaves: "ID"\n", ncomponents, nleaves) ;
    for (c = 0 ; c < ncomponents ; c++)
    {
	PR1 ((Common->file, "New node: "ID" new parent "ID"\n", c, Cparent [c])) ;
    }
#endif

#endif

    /* ---------------------------------------------------------------------- */
    /* allocate remainder of symbolic factor */
    /* ---------------------------------------------------------------------- */

    /* nc7 = 7*ncomponents + 2 */
    nc7 = CHOLMOD (mult_size_t) (ncomponents, 7, &ok) ;
    nc7 = CHOLMOD (add_size_t) (nc7, 2, &ok) ;
    if (!ok) PARAKLETE_ERROR (PK_TOO_LARGE, "problem too large") ;

    LUsymbolic->Mem_c = CHOLMOD (malloc) (nc7, sizeof (Int), cm) ;

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
    LUsymbolic->ncomp0 = nc0 ;
    LUsymbolic->ncomp1 = nc1 ;

    if (cm->status != CHOLMOD_OK)
    {
	/* out of memory or other error; inform all processes */
	amesos_paraklete_free_symbolic (&LUsymbolic, Common) ;
	CHOLMOD (free) (n5, sizeof (Int), W, cm) ;
	MPI (paraklete_bcast_symbolic (&LUsymbolic, Common)) ;
        PARAKLETE_ERROR (PK_OUT_OF_MEMORY, "out of memory") ;
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* Cstart = start and end nodes of each component */
    /* ---------------------------------------------------------------------- */

    clast = TRILINOS_CHOLMOD_EMPTY ;
    for (k = 0 ; k < n ; k++)
    {
	c = Cmember [Cperm [k]] ;
	if (c != clast)
	{
	    /* if (c != clast+1) { printf ("Error!\n") ; exit (1) ; } */
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
        size_t s = 0 ;
	for (k = Cstart [c] ; k < Cstart [c+1] ; k++)
	{
	    s = CHOLMOD (add_size_t) (s, ColCount [k], &ok) ;
	}
        if (!ok)
        {
            /* TODO return NULL, and broadcast error to all processes */
            PARAKLETE_ERROR (PK_TOO_LARGE, "problem too large") ;
        }
	/*
	printf ("Clnz ["ID"] = "ID", cols "ID" to "ID"\n", c, Clnz [c],
	    Cstart [c], Cstart [c+1]-1) ;
	*/
	Clnz [c] = s ;
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
	PR1 ((Common->file, "node "ID": parent "ID"\n", c, parent)) ;
	if (parent == TRILINOS_CHOLMOD_EMPTY)
	{
	    nroots++ ;
	}
	else
	{
	    ASSERT (parent > 0 && parent < ncomponents) ;
	    Cn [parent]++ ;
	}
    }

    if (nroots > 1)
    {
        /* TODO - this is an assertion */
        PARAKLETE_ERROR (PK_UNKNOWN, "separator tree cannot be forest") ;
        abort ( ) ;
    }

#ifndef NDEBUG
    for (c = 0 ; c < ncomponents ; c++)
    {
	PR1 ((Common->file, "node "ID": "ID" children\n", c, Cn [c])) ;
    }
#endif

    /* find the cumulative sum of the children */
    k = 0 ;
    for (c = 0 ; c < ncomponents ; c++)
    {
	Childp [c] = k ;
	k += Cn [c] ;
    }
    PR1 ((Common->file, "k "ID" ncomponents "ID"\n", k, ncomponents)) ;
    ASSERT (k == ncomponents - nroots) ;
    Childp [ncomponents] = k ;

    /* create a list of children for each node */
    for (c = 0 ; c < ncomponents ; c++)
    {
	Cn [c] = Childp [c] ;
    }
    for (k = 0 ; k < ncomponents ; k++)
    {
	Child [k] = -1 ;
    }
    for (c = 0 ; c < ncomponents ; c++)
    {
	parent = Cparent [c] ;
	if (parent != TRILINOS_CHOLMOD_EMPTY)
	{
	    Child [Cn [parent]++] = c ;
	}
    }

#ifndef NDEBUG
    for (c = 0 ; c < ncomponents ; c++)
    {
	PR1 ((Common->file, "Node "ID" children: ", c)) ;
	for (cp = Childp [c] ; cp < Childp [c+1] ; cp++)
	{
	    PR1 ((Common->file, ""ID" ", Child [cp])) ;
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
	nparent = (parent == TRILINOS_CHOLMOD_EMPTY) ? 0 : Cn [parent] ;
	Cn [c] = (Cstart [c+1] - Cstart [c]) + nparent ;
	PR1 ((Common->file, "node "ID" Cn: "ID"\n", c, Cn [c])) ;
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
	    PR3 ((Common->file, "A(1+"ID",1+"ID") = %g ; c(1+"ID",1+"ID") = "ID" ;\n",
		i,j, 
                A->x ? 0 : (((double *) (A->x)) [p]),
                i,j, c)) ;
	}
    }

#ifndef NDEBUG
    for (c = 0 ; c < ncomponents ; c++)
    {
	PR0 ((Common->file, "Cnz ["ID"] = "ID"\n", c, Cnz [c])) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* compute RpermInv = inverse of Cperm */
    /* ---------------------------------------------------------------------- */

    Rperm = LUsymbolic->Rperm ;
    ASSERT (Rperm != NULL) ;

    /* Rperm starts out equal to Cperm */
    for (k = 0 ; k < n ; k++)
    {
	Rperm [k] = Cperm [k] ;
    }
    for (k = 0 ; k < n ; k++)
    {
	RpermInv [Rperm [k]] = k ;
    }

    /* ---------------------------------------------------------------------- */
    /* schedule the processes to the nodes */
    /* ---------------------------------------------------------------------- */

#ifndef NMPI
    /* processes Lo_id [c] to Hi_id [c] are assigned to node c or descendants */
    Lo_id = W ;
    Hi_id = W + ncomponents ; 

    for (c = 0 ; c < ncomponents ; c++)
    {
	Sched [c] = TRILINOS_CHOLMOD_EMPTY ;
        Lo_id [c] = -1 ;
        Hi_id [c] = -1 ;
    }

    Lo_id [ncomponents-1] = 0 ;
    Hi_id [ncomponents-1] = nproc-1 ;

    for (c = ncomponents - 1 ; c >= 0 ; c--)
    {
        Int nchild, child, child_left, child_right, c_nproc ;

        /* first available process does this node */
        Sched [c] = Lo_id [c] ;

        /* split the processes amongst the children */
	nchild = Childp [c+1] - Childp [c] ;
        cp = Childp [c] ;

        if (nchild == 0)
        {
            /* nothing else to do */
        }
        else if (nchild == 1)
        {
            /* all processes go to the one child */
            child = Child [cp] ;
            Lo_id [child] = Lo_id [c] ;
            Hi_id [child] = Hi_id [c] ;
        }
        else if (nchild == 2)
        {
            /* two children; split the processors between them */
            child_left  = Child [cp] ;
            child_right = Child [cp+1] ;

            c_nproc = Hi_id [c] - Lo_id [c] + 1 ;
            if (c_nproc > 1)
            {
                Lo_id [child_left ] = Lo_id [c] ;
                Hi_id [child_left ] = Lo_id [c] + c_nproc / 2 - 1 ;
                Lo_id [child_right] = Lo_id [c] + c_nproc / 2 ;
                Hi_id [child_right] = Hi_id [c] ;
            }
            else
            {
                Lo_id [child_left ] = Lo_id [c] ;
                Hi_id [child_left ] = Lo_id [c] ;
                Lo_id [child_right] = Lo_id [c] ;
                Hi_id [child_right] = Lo_id [c] ;
            }
        }
        else
        {
            /* TODO this is an assertion - it cannot occur */
            PARAKLETE_ERROR (PK_UNKNOWN, "invalid separator tree") ;
            abort ( ) ;
        }
    }

#if 0
    proc = 0 ;
    for (c = 0 ; c < ncomponents ; c++)
    {
	Int nchild = Childp [c+1] - Childp [c] ;
	if (nchild == 0)
	{
	    PR1 ((Common->file,"\nSchedule child "ID" to process "ID"\n", c, proc));
	    for (cc = c ; cc != TRILINOS_CHOLMOD_EMPTY && Sched [cc] == TRILINOS_CHOLMOD_EMPTY ; cc = Cparent[cc])
	    {
		PR1 ((Common->file, "  node "ID" to process "ID"\n", cc, proc)) ;
		Sched [cc] = proc ;
	    }
            /* advance to the next process */
	    proc = (proc + 1) % nproc ;
	}
    }
#endif

#else
    /* all components are done by process zero */
    proc = 0 ;
    for (c = 0 ; c < ncomponents ; c++)
    {
	Sched [c] = proc ;
    }
#endif

#ifndef NDEBUG
    PR0 ((Common->file, "\nncomponents:: "ID"\n", ncomponents)) ;
    for (c = 0 ; c < ncomponents ; c++)
    {
	PR0 ((Common->file, "    node "ID" Sched "ID" : Cparent "ID" proc "ID"\n",
		    c, Sched [c], Cparent [c],
		    (Cparent [c] == TRILINOS_CHOLMOD_EMPTY) ? TRILINOS_CHOLMOD_EMPTY : Sched [Cparent [c]])) ;
    }
#endif

#if 0
    for (c = 0 ; c < ncomponents ; c++)
    {
        printf ("   node "ID" on "ID" : Cparent "ID" on "ID"\n",
		    c, Sched [c], Cparent [c],
		    (Cparent [c] == TRILINOS_CHOLMOD_EMPTY) ? TRILINOS_CHOLMOD_EMPTY : Sched [Cparent [c]]) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* return the symbolic factorization, and broadcast it to all processes */
    /* ---------------------------------------------------------------------- */

    CHOLMOD (free) (n5, sizeof (Int), W, cm) ;
    PR1 ((Common->file, "analysis done\n")) ;
    MPI (paraklete_bcast_symbolic (&LUsymbolic, Common)) ;
    return (LUsymbolic) ;
}


/* ========================================================================== */
/* === paraklete_alloc_symbolic ============================================= */
/* ========================================================================== */

/* allocate a symbolic object */

paraklete_symbolic *amesos_paraklete_alloc_symbolic
(
    Int n,
    Int ncomponents,
    Int do_Rperm,
    paraklete_common *Common
)
{
    paraklete_symbolic *LUsymbolic ;
    cholmod_common *cm ;
    size_t n3, nc7 ;
    int ok = TRUE ;

    cm = &(Common->cm) ;

    /* n3 = 3*n */
    n3 = CHOLMOD (mult_size_t) (n, 3, &ok) ;
    if (!ok) PARAKLETE_ERROR (PK_TOO_LARGE, "problem too large") ;

    /* nc7 = 7*ncomponents + 2 */
    nc7 = CHOLMOD (mult_size_t) (ncomponents, 7, &ok) ;
    nc7 = CHOLMOD (add_size_t) (nc7, 2, &ok) ;
    if (!ok) PARAKLETE_ERROR (PK_TOO_LARGE, "problem too large") ;

    LUsymbolic = CHOLMOD (malloc) (1, sizeof (paraklete_symbolic), cm) ;

    if (cm->status != CHOLMOD_OK)
    {
        PARAKLETE_ERROR (PK_OUT_OF_MEMORY, "out of memory") ;
	return (NULL) ;
    }

    if (n > 0)
    {
	LUsymbolic->Mem_n    = CHOLMOD (malloc) (n3, sizeof (Int), cm) ;
	LUsymbolic->Cperm    = LUsymbolic->Mem_n ;
	LUsymbolic->RpermInv = LUsymbolic->Mem_n + n ;
	LUsymbolic->Cparent  = LUsymbolic->Mem_n + 2*n ;

	if (do_Rperm)
	{
	    LUsymbolic->Rperm = CHOLMOD (malloc) (n, sizeof (Int), cm) ;
	}
	else
	{
	    /* fill-reducing ordering is symmetric, Rperm is implicitly
	     * equal to Cperm */
	    LUsymbolic->Rperm = NULL ;
	}
    }
    else
    {
	LUsymbolic->Mem_n    = NULL ;
	LUsymbolic->Cperm    = NULL ;
	LUsymbolic->RpermInv = NULL ;
	LUsymbolic->Cparent  = NULL ;
	LUsymbolic->Rperm    = NULL ;
    }

    if (ncomponents > 0)
    {

	LUsymbolic->Mem_c  = CHOLMOD (malloc) (nc7, sizeof(Int), cm) ;

	/* each of size ncomponents: */
	LUsymbolic->Child  = LUsymbolic->Mem_c ;
	LUsymbolic->Clnz   = LUsymbolic->Mem_c + ncomponents ;
	LUsymbolic->Cn     = LUsymbolic->Mem_c + 2*ncomponents ;
	LUsymbolic->Cnz    = LUsymbolic->Mem_c + 3*ncomponents ;
	LUsymbolic->Sched  = LUsymbolic->Mem_c + 4*ncomponents ;

	/* each of size ncomponents+1: */
	LUsymbolic->Cstart = LUsymbolic->Mem_c + 5*ncomponents ;
	LUsymbolic->Childp = LUsymbolic->Mem_c + 6*ncomponents + 1 ;

    }
    else
    {
	LUsymbolic->Mem_c  = NULL ;
	LUsymbolic->Child  = NULL ;
	LUsymbolic->Clnz   = NULL ;
	LUsymbolic->Cn     = NULL ;
	LUsymbolic->Cnz    = NULL ;
	LUsymbolic->Sched  = NULL ;
	LUsymbolic->Cstart = NULL ;
	LUsymbolic->Childp = NULL ;
    }

    LUsymbolic->n = n ;
    LUsymbolic->ncomponents = ncomponents ;

    if (cm->status != CHOLMOD_OK)
    {
        /* out of memory */
        PARAKLETE_ERROR (PK_OUT_OF_MEMORY, "out of memory") ;
	amesos_paraklete_free_symbolic (&LUsymbolic, Common) ;
    }

    return (LUsymbolic) ;
}


/* ========================================================================== */
/* === paraklete_free_symbolic ============================================== */
/* ========================================================================== */

/* Free the symbolic object.  All processes own a copy after it is broadcast. */

void amesos_paraklete_free_symbolic
(
    paraklete_symbolic **LUsymbolicHandle,
    paraklete_common *Common
)
{
    paraklete_symbolic *LUsymbolic ;
    cholmod_common *cm ;
    Int n, ncomponents ;

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

    /* size-3n space, in Mem_n */
    CHOLMOD (free) (3*n, sizeof (Int), LUsymbolic->Mem_n, cm) ;
    LUsymbolic->Cperm = NULL ; 
    LUsymbolic->RpermInv = NULL ; 
    LUsymbolic->Cparent = NULL ; 

    /* size-n space, only used for reanalyze/refactorize, or if the fill-
     * reducing ordering is unsymmetric.  Otherwise, Rperm is implicitly
     * equal to Cperm. */
    CHOLMOD (free) (n, sizeof (Int), LUsymbolic->Rperm, cm) ;
    LUsymbolic->Rperm = NULL ; 

    /* size-(7*components+2) space, in Mem_c */
    CHOLMOD (free) (7*ncomponents + 2, sizeof (Int), LUsymbolic->Mem_c, cm) ;
    LUsymbolic->Cstart = NULL ;
    LUsymbolic->Child = NULL ;
    LUsymbolic->Childp = NULL ;
    LUsymbolic->Clnz = NULL ;
    LUsymbolic->Cn = NULL ;
    LUsymbolic->Cnz = NULL ;
    LUsymbolic->Sched = NULL ;

    *LUsymbolicHandle = CHOLMOD (free) (
	    1, sizeof (paraklete_symbolic), (*LUsymbolicHandle), cm) ;
}
