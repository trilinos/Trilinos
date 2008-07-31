/* ========================================================================== */
/* === paraklete_solve ====================================================== */
/* ========================================================================== */

#include "amesos_paraklete_decl.h"

/* paraklete_solve (LU, LUsymbolic, B, Common) solves Ly=Pb, where P is
 * the initial fill-reducing ordering (P=Cperm).
 *
 * The solution is left in LU->LUnode [c]->X, for each node c.  Next, it solves
 * UQ'x=y, where the solution is written into the input array B,
 * and where Q is the final column ordering.
 *
 * This solve is performed using the original LUsymbolic object, not the one
 * returned from paraklete_reanalyze.
 *
 * PARAKLETE version 0.3: parallel sparse LU factorization.  Nov 13, 2007
 * Copyright (C) 2007, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

Int amesos_paraklete_solve
(
    /* inputs, not modified */
    paraklete_numeric *LU,
    paraklete_symbolic *LUsymbolic,
    double *B,
    paraklete_common *Common
)
{
    cholmod_common *cm ;
    double *Blocal, *X, *W, *X2 ;
    Int *Cperm, *Rperm, *Cstart, *Q, *Sched, *Child, *Childp ;
    Int i, n, ncomponents, k, c, k1, k2, nfound, myid, cp, nchild, child ;
    MPI (MPI_Status ms) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    cm = &(Common->cm) ;
    n = LU->n ;
    W = LU->W ;	    /* only non-NULL on the root process 0 */
    myid = Common->myid ;

    ncomponents = LUsymbolic->ncomponents ;
    Cperm  = LUsymbolic->Cperm ;
    Rperm  = LUsymbolic->Rperm ;
    Rperm  = (Rperm == NULL) ? Cperm : Rperm ;
    Cstart = LUsymbolic->Cstart ;
    Sched  = LUsymbolic->Sched ;

    /* ---------------------------------------------------------------------- */
    /* W = Rperm*B, on process 0 only */
    /* ---------------------------------------------------------------------- */

    if (myid == 0)
    {
	for (k = 0 ; k < n ; k++)
	{
	    W [k] = B [Rperm [k]] ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* distribute the permuted B to the nodes */
    /* ---------------------------------------------------------------------- */

    for (c = 0 ; c < ncomponents ; c++)
    {
	k1 = Cstart [c] ;
	k2 = Cstart [c+1] ;

	Blocal = LU->LUnode [c]->B ;

	/* send Blocal to node c */
	MPI (LU->LUnode [c]->req = MPI_REQUEST_NULL) ;
	if (myid == 0)
	{
	    if (Sched [c] == myid)
	    {
		for (i = 0 ; i < k2-k1 ; i++)
		{
		    Blocal [i] = W [k1 + i] ;
		}
	    }
	    else
	    {
		MPI (MPI_Isend (W + k1, k2-k1, MPI_DOUBLE, Sched [c],
			/* TAG: */ c, MPI_COMM_WORLD, &(LU->LUnode [c]->req))) ;
	    }
	}
	else
	{
	    if (Sched [c] == myid)
	    {
		MPI (MPI_Irecv (Blocal, k2-k1, MPI_DOUBLE, 0,
			/* TAG: */ c, MPI_COMM_WORLD, &(LU->LUnode [c]->req))) ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* forward solve: post a receive for X from each non-local child */
    /* ---------------------------------------------------------------------- */

    Childp = LUsymbolic->Childp ;	/* Child [Childp [c]..Childp[c+1]-1] */
    Child  = LUsymbolic->Child ;	/* is list of children of node c */

    for (c = 0 ; c < ncomponents ; c++)
    {
	if (Sched [c] == myid)
	{
	    nchild = Childp [c+1] - Childp [c] ;
	    for (cp = 0 ; cp < nchild ; cp++)
	    {
		child = Child [Childp [c] + cp] ;
		if (Sched [child] != myid)
		{
		    MPI (MPI_Irecv (LU->LUnode [child]->X,
			LU->LUnode [child]->PK_NN, MPI_DOUBLE, Sched [child],
			TAG0, MPI_COMM_WORLD, &(LU->LUnode [c]->Req [cp]))) ;
		}
		else
		{
		    MPI (LU->LUnode [c]->Req [cp] = MPI_REQUEST_NULL) ;
		}
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* forward solve: Ly=b */
    /* ---------------------------------------------------------------------- */

    for (c = 0 ; c < ncomponents ; c++)
    {
	if (Sched [c] == myid)
	{
	    amesos_paraklete_lsolve_node (c, LU, LUsymbolic, Common) ;
	}
    }

    MPI (MPI_Barrier (MPI_COMM_WORLD)) ;

    /* ---------------------------------------------------------------------- */
    /* backsolve: Ux=y */
    /* ---------------------------------------------------------------------- */

    for (c = ncomponents-1 ; c >= 0 ; c--)
    {
	if (Sched [c] == myid)
	{
	    amesos_paraklete_usolve_node (c, LU, LUsymbolic, Common) ;
	}
    }

    MPI (MPI_Barrier (MPI_COMM_WORLD)) ;

    /* ---------------------------------------------------------------------- */
    /* get the permuted solution from each node */
    /* ---------------------------------------------------------------------- */

    Q = LU->Q ;
    k = 0 ;
    for (c = 0 ; c < ncomponents ; c++)
    {
	X = LU->LUnode [c]->X ;
	nfound = LU->LUnode [c]->PK_NFOUND ;
	if (myid == 0)
	{
	    PR1 ((Common->file, "get soln, node c="ID", nfound "ID"\n", c, nfound));
	    /* get X from node c */
	    if (Sched [c] != myid)
	    {
		PR1 ((Common->file, "recv node "ID" from "ID"\n", c, Sched [c])) ;
		MPI (MPI_Recv (W, nfound, MPI_DOUBLE, Sched [c],
			    TAG0, MPI_COMM_WORLD, &ms)) ;
		X2 = W ;
	    }
	    else
	    {
		PR1 ((Common->file, "I own it already\n")) ;
		X2 = X ;
	    }
	    PR1 ((Common->file, "got it from Sched [c] = "ID"\n", Sched [c])) ;
	    for (i = 0 ; i < nfound ; i++)
	    {
		B [Q [k]] = X2 [i] ;
		PR2 ((Common->file, "X ["ID"] is global B ["ID"] %g\n",
		    i, k, X2 [i])) ;
		k++ ;
	    }
	}
	else
	{
	    if (Sched [c] == myid)
	    {
		PR1 ((Common->file,
		    "send soln, node c = "ID", myid "ID" nfound "ID"\n",
		    c, myid, nfound)) ;
		MPI (MPI_Send (X, nfound, MPI_DOUBLE, 0, TAG0, MPI_COMM_WORLD));
	    }
	}
    }

    return (TRUE) ;
}
