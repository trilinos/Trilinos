/* ========================================================================== */
/* === paraklete_reanalyze ================================================== */
/* ========================================================================== */

/* paraklete_analyze and paraklete_factorize have already been called.
 * This function constructs a new symbolic object that includes the original
 * fill-reducing nested dissection ordering (Cperm) and separator tree
 * from paraklete_analyze, combined with the partial-pivoting permutation
 * and lost-pivots from paraklete_factorize.
 *
 * All processes should do this, independently, since all processes have their
 * own complete copy of the LUsymbolic object.  Each process also has its own
 * copy of LU->P and LU->Q, the final row and column permutations from
 * paraklete_factorize.
 *
 * Everyone has the global P, Q, and LU->LUnode [*]->header.  The assignment
 * of rows/columns to the nodes can be finalized, to reflect the final
 * permutation.  All processes compute the new schedule.
 *
 * PARAKLETE version 0.3: parallel sparse LU factorization.  Nov 13, 2007
 * Copyright (C) 2007, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

#include "amesos_paraklete_decl.h"

paraklete_symbolic *amesos_paraklete_reanalyze
(
    cholmod_sparse *A,	    /* only root processor owns this */
    paraklete_numeric *LU,
    paraklete_symbolic *LUsymbolic,
    paraklete_common *Common
)
{
    paraklete_symbolic *LUsymbolic_new ;
    Int *Cparent, *Cstart, *Cn, *Cmember, *Pmember, *P, *Q, *Cnz, *Ap, *Ai ;
    Int i, j, k, c, n, ncomponents, jold, iold, parent, nparent, ci, cj, p ;

    /* TODO check inputs to make sure they are not NULL */

    /* ---------------------------------------------------------------------- */
    /* allocate the new symbolic object */
    /* ---------------------------------------------------------------------- */

    n = LUsymbolic->n ;
    ncomponents = LUsymbolic->ncomponents ;
    LUsymbolic_new = amesos_paraklete_alloc_symbolic (n, ncomponents, TRUE, Common) ;

    /* TODO check for out-of-memory here */

    /* ---------------------------------------------------------------------- */
    /* copy the parts of the separator tree that do not change with pivoting */
    /* ---------------------------------------------------------------------- */

    for (c = 0 ; c < ncomponents ; c++)
    {
	LUsymbolic_new->Child [c]  = LUsymbolic->Child [c] ;
    }

    for (c = 0 ; c <= ncomponents ; c++)
    {
	LUsymbolic_new->Childp [c] = LUsymbolic->Childp [c] ;
    }

    for (c = 0 ; c < ncomponents ; c++)
    {
	LUsymbolic_new->Sched [c] = LUsymbolic->Sched [c] ;
    }

    Cparent = LUsymbolic_new->Cparent ;
    for (c = 0 ; c < ncomponents ; c++)
    {
	Cparent [c] = LUsymbolic->Cparent [c] ;
    }

    /* ---------------------------------------------------------------------- */
    /* determine actual number of entries in LU factors at each node */
    /* ---------------------------------------------------------------------- */

    for (c = 0 ; c < ncomponents ; c++)
    {
	LUsymbolic_new->Clnz [c] = LU->LUnode [c]->lnz + LU->LUnode [c]->unz ;
    }

    /* ---------------------------------------------------------------------- */
    /* find the range of nodes in P*A*Q assigned to each node */
    /* ---------------------------------------------------------------------- */

    Cstart = LUsymbolic_new->Cstart ;
    Cstart [0] = 0 ;
    for (c = 0 ; c < ncomponents ; c++)
    {
	Cstart [c+1] = Cstart [c] + LU->LUnode [c]-> PK_NFOUND ;
    }

    /*
    if (Common->myid == 0)
	for (c = 0 ; c <= ncomponents ; c++)
	    printf ("new Cstart "ID"\n", Cstart [c]) ;
    */

    /* ---------------------------------------------------------------------- */
    /* find size of C matrix at each node (all pivoting accounted for) */
    /* ---------------------------------------------------------------------- */

    Cn = LUsymbolic_new->Cn ;
    for (c = ncomponents - 1 ; c >= 0 ; c--)
    {
	parent = Cparent [c] ;
	nparent = (parent == TRILINOS_CHOLMOD_EMPTY) ? 0 : Cn [parent] ;
	Cn [c] = (Cstart [c+1] - Cstart [c]) + nparent ;
	PR1 ((Common->file, "node "ID" new Cn: "ID"\n", c, Cn [c])) ;
    }

    /* ---------------------------------------------------------------------- */
    /* determine Cnz for each node */
    /* ---------------------------------------------------------------------- */

    Cmember = LUsymbolic_new->Cperm ;	/* use Cperm as workspace */
    Pmember = LUsymbolic_new->Rperm ;	/* use Rperm as workspace */

    Cnz = LUsymbolic_new->Cnz ;		/* new nonzero count for each node */
    Q = LU->Q ;				/* final column ordering, P*A*Q=L*U */
    P = LU->P ;				/* final row ordering */

    if (Common->myid == 0)
    {
	for (c = 0 ; c < ncomponents ; c++)
	{
	    Cnz [c] = 0 ;
	}

	/* Cmember [k] = c if the kth diagonal entry of P*A*Q is in node c */
	for (c = 0 ; c < ncomponents ; c++)
	{
	    for (k = Cstart [c] ; k < Cstart [c+1] ; k++)
	    {
		Cmember [k] = c ;
	    }
	}

	/*
	for (k = 0 ; k < n ; k++)
	{
	    printf ("Cmember ["ID"] = "ID"\n", k, Cmember [k]) ;
	}
	for (k = 0 ; k < n ; k++)
	{
	    printf ("P ["ID"] = "ID"\n", k, P [k]) ;
	}
	*/

	/* Pmember [i] = c if A(i,i) becomes a pivot in node c */
	for (k = 0 ; k < n ; k++)
	{
	    i = P [k] ;		/* kth row of P*A*Q is ith row of A */
	    c = Cmember [k] ;	/* kth row of P*A*Q is in node c */
	    Pmember [i] = c ;	/* A(i,i) becomes (P*A*Q)_kk, in node c */
	}

	/*
	for (i = 0 ; i < n ; i++)
	{
	    printf ("Pmember ["ID"] = "ID"\n", i, Pmember [i]) ;
	}
	*/

	/* count the entries in each node */
	Ap = A->p ;
	Ai = A->i ;
	for (j = 0 ; j < n ; j++)
	{
	    jold = Q [j] ;	/* jth column of P*A*Q is column jold of A */
	    cj = Cmember [j] ;	/* jth diagonal entry of P*A*Q in node cj */
	    for (p = Ap [jold] ; p < Ap [jold+1] ; p++)
	    {
		iold = Ai [p] ;		/* A(iold,jold) in original matrix */
		ci = Pmember [iold] ;	/* A(iold,iold) in node ci */
		c = MIN (ci, cj) ;	/* determine node of (P*A*Q)_ij */

		/*
		printf ("A("ID","ID") PAQ("ID","ID") ci "ID" cj "ID" c "ID"\n",
			iold, jold, Ai [p], j, ci, cj, c) ;
		*/
		Cnz [c]++ ;
	    }
	}

	/*
	for (c = 0 ; c < ncomponents ; c++)
	{
	    printf ("component "ID"\n", c) ;
	    printf ("   new Cn "ID" Cnz "ID" Cstart "ID"\n",
			Cn [c], Cnz [c], Cstart [c]) ;
	}
	*/

    }

    /* Cmember and Pmember workspace no longer needed */

    /* broadcast results to all processors */
    MPI (MPI_Bcast (Cnz, ncomponents, MPI_Int, TAG0, MPI_COMM_WORLD)) ;

    /* ---------------------------------------------------------------------- */
    /* copy the final P and Q */
    /* ---------------------------------------------------------------------- */

    /* column permutation Q */
    for (k = 0 ; k < n ; k++)
    {
	LUsymbolic_new->Cperm [k] = Q [k] ;
    }

    /* row permutation P and its inverse */
    for (k = 0 ; k < n ; k++)
    {
	i = P [k] ;
	LUsymbolic_new->RpermInv [i] = k ;
	LUsymbolic_new->Rperm [k] = i ;
    }

    return (LUsymbolic_new) ;
}
