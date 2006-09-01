/* ========================================================================== */
/* === paraklete_lsolve_node ================================================ */
/* ========================================================================== */

#include "paraklete.h"

/* Solve Lx=b with node c of the separator tree.
 *
 * PARAKLETE version 0.1: parallel sparse LU factorization.  May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

int paraklete_lsolve_node
(
    int c,
    paraklete_numeric *LU,
    paraklete_symbolic *LUsymbolic,
    paraklete_common *Common
)

{
    double xj ;
    paraklete_node *LUnode ;
    cholmod_common *cm ;
    double *X, *B, *LUix, *Xchild, *Lx ;
    int *Child, *Childp, *Lost, *Lostp, *Cstart, *Lip, *Llen,
	*Cn, *Pinv, *Li, *Cparent, *Sched ;
    int cp, nchild, child, cn, k, j, p, i, nfound, k1, k2, llen,
	ci, nlost_in, cn_child, cn_nfound, parent, myid, pass ;
    MPI (MPI_Status ms) ;
    MPI (MPI_Request req) ;

    /* ---------------------------------------------------------------------- */
    /* get local workspace */
    /* ---------------------------------------------------------------------- */

    cm = &(Common->cm) ;
    PR0 ((Common->file, "\n\n########################## Lsolve NODE %d\n", c)) ;

    /* ---------------------------------------------------------------------- */
    /* get the symbolic analysis of this node */
    /* ---------------------------------------------------------------------- */

    myid = Common->myid ;
    Childp = LUsymbolic->Childp ;	/* Child [Childp [c]..Childp[c+1]-1] */
    Child  = LUsymbolic->Child ;	/* is list of children of node c */
    nchild = Childp [c+1] - Childp [c] ;    /* # of children of node c */
    Cn     = LUsymbolic->Cn ;		/* dimension of each node */
    Cstart = LUsymbolic->Cstart ;
    Sched = LUsymbolic->Sched ;
    Cparent = LUsymbolic->Cparent ;

    /* ---------------------------------------------------------------------- */
    /* get factors at this node */
    /* ---------------------------------------------------------------------- */

    LUnode = LU->LUnode [c] ;

    Lip = LUnode->lp ;
    Llen = LUnode->llen ;
    LUix = LUnode->ix ;
    nfound = LUnode->PK_NFOUND ;

    nlost_in = 0 ;
    Lost  = LUnode->Lost ;
    Lostp = LUnode->Lostp ;
    nlost_in = LUnode->nlost_in ;

    /* The matrix factorized at this node is cn-by-cn. */
    cn = nlost_in + Cn [c] ;
    k1 = Cstart [c] ;
    k2 = Cstart [c+1] ;

    Pinv = LUnode->Pinv ;

#ifndef NDEBUG
    {
	int npiv = LUnode->PK_NPIV ;
	int *Plocal = LUnode->Plocal ;
	ASSERT (cn == LUnode->PK_NN) ;
	ASSERT (npiv == nlost_in + (k2 - k1)) ;
	for (k = 0 ; k < npiv ; k++)
	{
	    i = Plocal [k] ;
	    ASSERT (i >= 0 && i < npiv) ;
	    ASSERT (Pinv [i] == k) ;
	}
	DEBUG (cholmod_print_perm (Plocal, npiv, npiv, "Plocal at node c", cm));
	DEBUG (cholmod_print_perm (Pinv,   npiv, npiv, "Pinv at node c", cm)) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* add up all the contributions to the right-hand-side */
    /* ---------------------------------------------------------------------- */

    X = LU->LUnode [c]->X ;
    PR1 ((Common->file, "Lsolve at Node %d, cn %d\n", c, cn)) ;
    for (i = 0 ; i < cn ; i++)
    {
	X [i] = 0 ;
    }

    /* B [0..(k2-k1)-1] contains the right-hand-side corresponding to rows
     * Cstart [c] to Cstart [c+1]-1 of the global system (after permuted by
     * Cperm).
     */
    B = LUnode->B ;

    /* wait for B to arrive at this node */
    MPI (MPI_Wait (&(LU->LUnode [c]->req), &ms)) ;

    for (i = 0 ; i < k2-k1 ; i++)
    {
	ci = i + nlost_in ;
	k = Pinv [ci] ;
	PR2 ((Common->file, "orig B [%d] = %g goes to X [%d]\n", i, B [i], k)) ;
	ASSERT (k >= 0 && k < cn) ;
	X [k] = B [i] ;
    }

    /* scatter and add the contributions from each child */
    for (pass = 1 ; pass <= 2 ; pass++)
    {
	for (cp = 0 ; cp < nchild ; cp++)
	{
	    child = Child [Childp [c] + cp] ;
	    if (pass == 1)
	    {
		/* do local contributions on first pass */
		if (Sched [child] != myid) continue ;
	    }
	    else
	    {
		/* do non-local contributions on second pass */
		if (Sched [child] == myid) continue ;
		MPI (MPI_Wait (&(LUnode->Req [cp]), &ms)) ;
	    }

	    cn_child = LU->LUnode [child]->PK_NN ;
	    cn_nfound = LU->LUnode [child]->PK_NFOUND ;
	    Xchild = LU->LUnode [child]->X + cn_nfound ;

	    PR1 ((Common->file, "child %d cn %d nfound %d\n",
		    child, cn_child, cn_nfound)) ;

	    /* get the contributions of child to its lost pivot rows */
	    for (i = 0 ; i < Lost [cp] ; i++)
	    {
		ci = i + Lostp [cp] ;	/* ci is now "original" local index */
		k = Pinv [ci] ;		/* kth local pivot row */
		PR2 ((Common->file, "Xchild [%d] = %g goes to X [%d] (lost)\n",
			i, Xchild [i], k)) ;
		ASSERT (k >= 0 && k < cn) ;
		X [k] += Xchild [i] ;
	    }

	    /* get the contributions to candidate pivot rows of node c */
	    for ( ; i < Lost [cp] + (k2-k1) ; i++)
	    {
		ci = i + (nlost_in - Lost [cp]) ;
		k = Pinv [ci] ;		/* kth local pivot row */
		PR2 ((Common->file, "Xchild [%d] = %g goes to X [%d] (cand)\n",
			i, Xchild [i], k)) ;
		ASSERT (k >= 0 && k < cn) ;
		X [k] += Xchild [i] ;
	    }

	    /* get contributions to candidate pivot rows of ancestors of c */
	    for ( ; i < cn_child - cn_nfound ; i++)
	    {
		k = i + (nlost_in - Lost [cp]) ;
		PR2 ((Common->file, "Xchild [%d] = %g goes to X [%d] (anc)\n",
			i, Xchild [i], k)) ;
		ASSERT (k >= 0 && k < cn) ;
		X [k] += Xchild [i] ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* solve Lx=b */
    /* ---------------------------------------------------------------------- */

    for (j = 0 ; j < nfound ; j++)
    {
	xj = X [j] ;
	GET_COLUMN (Lip, Llen, LUix, j, Li, Lx, llen) ;
	for (p = 1 ; p < llen ; p++)
	{
	    X [Li [p]] -= Lx [p] * xj ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* send results to parent */
    /* ---------------------------------------------------------------------- */

    /* X [0..nfound-1] is part of the global solution, X [nfound..nn-1] is
     * the contribution to the parent and ancestors */

    parent = Cparent [c] ;
    if (parent != EMPTY && Sched [parent] != myid)
    {
	MPI (MPI_Isend (X, cn, MPI_DOUBLE, Sched [parent], TAG0, Common->mpicomm,
	    &req)) ;
    }

    return (TRUE) ;
}
