/* ========================================================================== */
/* === paraklete_usolve_node ================================================ */
/* ========================================================================== */

#include "paraklete.h"

/* Solve Ux=y with node c of the separator tree
 *
 * PARAKLETE version 0.1: parallel sparse LU factorization.  May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

int paraklete_usolve_node
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
    double *LUix, *Xchild, *Ux, *X ;
    int *Child, *Childp, *Lost, *Lostp, *Cstart, *Uip, *Ulen,
	*Cn, *Qinv, *Ui, *Sched, *Cparent ;
    int cp, nchild, child, cn, k, j, p, i, nfound, k1, k2, ulen,
	ci, nlost_in, cn_child, cn_nfound, myid, parent ;
    MPI (MPI_Status ms) ;
    MPI (MPI_Request req) ;

    /* ---------------------------------------------------------------------- */
    /* get local workspace */
    /* ---------------------------------------------------------------------- */

    cm = &(Common->cm) ;
    PR0 ((Common->file, "\n\n########################### Usolve NODE %d\n", c)) ;

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
    /* get solution to Ly=Pb, and solution to Ux=y from parent */
    /* ---------------------------------------------------------------------- */

    LUnode = LU->LUnode [c] ;
    cn = LUnode->PK_NN ;
    nfound = LUnode->PK_NFOUND ;
    X = LUnode->X ;
    parent = Cparent [c] ;

    PR1 ((Common->file, "Usolve node %d cn %d nfound %d\n", c, cn, nfound)) ;

    if (parent != EMPTY && Sched [parent] != myid)
    {
	PR1 ((Common->file, "Recving usolve from %d, size %d\n", Sched [parent],
		cn - nfound)) ;
	MPI (MPI_Irecv (X + nfound, cn-nfound, MPI_DOUBLE, Sched [parent],
		TAG0, MPI_COMM_WORLD, &req)) ;
	MPI (MPI_Wait (&req, &ms)) ;
    }

    /* ---------------------------------------------------------------------- */
    /* solve Ux=y at this node */
    /* ---------------------------------------------------------------------- */

    Uip = LUnode->up ;
    Ulen = LUnode->ulen ;
    LUix = LUnode->ix ;

    /* X1 = U2*X2 - X1 */
    for (j = cn-1 ; j >= nfound ; j--)
    {
	xj = X [j] ;
	GET_COLUMN (Uip, Ulen, LUix, j, Ui, Ux, ulen) ;
	for (p = 0 ; p < ulen ; p++)
	{
	    X [Ui [p]] -= Ux [p] * xj ;
	}
    }

    /* solve U*X1 */
    for ( ; j >= 0 ; j--)
    {
	GET_COLUMN (Uip, Ulen, LUix, j, Ui, Ux, ulen) ;
	xj = X [j] / Ux [ulen-1] ;
	X [j] = xj ;
	for (p = 0 ; p < ulen-1 ; p++)
	{
	    X [Ui [p]] -= Ux [p] * xj ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* get factors at this node */
    /* ---------------------------------------------------------------------- */

    nlost_in = 0 ;
    Lost  = LUnode->Lost ;
    Lostp = LUnode->Lostp ;
    nlost_in = LUnode->nlost_in ;

    /* The matrix factorized at this node is cn-by-cn. */
    cn = nlost_in + Cn [c] ;
    k1 = Cstart [c] ;
    k2 = Cstart [c+1] ;

    PR1 ((Common->file, "NODE %d nlost_in: %d\n", c, nlost_in)) ;

    Qinv = LUnode->Qinv ;

#ifndef NDEBUG
    {
	int npiv = LUnode->PK_NPIV ;
	int *Qlocal = LUnode->Qlocal ;
	ASSERT (cn == LUnode->PK_NN) ;
	ASSERT (npiv == nlost_in + (k2 - k1)) ;
	for (k = 0 ; k < npiv ; k++)
	{
	    i = Qlocal [k] ;
	    ASSERT (i >= 0 && i < npiv) ;
	    ASSERT (Qinv [i] == k) ;
	}
	DEBUG (cholmod_print_perm (Qlocal, npiv, npiv, "Qlocal at node c", cm));
	DEBUG (cholmod_print_perm (Qinv,   npiv, npiv, "Qinv at node c", cm)) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* send solutions to each child */
    /* ---------------------------------------------------------------------- */

    for (cp = 0 ; cp < nchild ; cp++)
    {
	child = Child [Childp [c] + cp] ;
	cn_child = LU->LUnode [child]->PK_NN ;
	cn_nfound = LU->LUnode [child]->PK_NFOUND ;
	Xchild = LU->LUnode [child]->X + cn_nfound ;

	PR1 ((Common->file, "child %d cn %d nfound %d\n",
		child, cn_child, cn_nfound)) ;

	/* send the child its lost pivot cols */
	for (i = 0 ; i < Lost [cp] ; i++)
	{
	    ci = i + Lostp [cp] ;	/* ci is now "original" local index */
	    k = Qinv [ci] ;		/* kth local pivot col */
	    PR2 ((Common->file, "Xchild [%d] = %g from X [%d] (lost)\n",
		    i, Xchild [i], k)) ;
	    Xchild [i] = X [k] ;
	}

	/* send the child the rest of the solution from c */
	for ( ; i < Lost [cp] + (k2-k1) ; i++)
	{
	    ci = i + (nlost_in - Lost [cp]) ;
	    k = Qinv [ci] ;		/* kth local pivot col */
	    PR2 ((Common->file, "Xchild [%d] = %g from X [%d] (cand)\n",
		    i, Xchild [i], k)) ;
	    Xchild [i] = X [k] ;
	}

	/* get the solutions of ancestors of c */
	for ( ; i < cn_child - cn_nfound ; i++)
	{
	    k = i + (nlost_in - Lost [cp]) ;
	    PR2 ((Common->file, "Xchild [%d] = %g from X [%d] (anc)\n",
		    i, Xchild [i], k)) ;
	    Xchild [i] = X [k] ;
	}

	PR1 ((Common->file, "Usolve: Sending from %d to child %d\n", c, child)) ;
	if (Sched [child] != myid)
	{
	    PR1 ((Common->file, "Sending to %d, size %d\n", Sched [child], 
		    cn_child - cn_nfound)) ;
	    MPI (MPI_Isend (Xchild, cn_child - cn_nfound, MPI_DOUBLE,
			Sched [child], TAG0, MPI_COMM_WORLD, &req)) ;
	}
    }

    return (TRUE) ;
}
