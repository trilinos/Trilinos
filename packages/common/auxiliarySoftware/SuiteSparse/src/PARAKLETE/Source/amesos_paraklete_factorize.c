/* ========================================================================== */
/* === paraklete_factorize ================================================== */
/* ========================================================================== */

#include "amesos_paraklete_decl.h"

/* LU = paraklete_factorize (A, LUsymbolic, Common) factorizes P*A*Q into L*U.
 * Returns NULL if A is singular or if memory is exhausted.
 *
 * PARAKLETE version 0.3: parallel sparse LU factorization.  Nov 13, 2007
 * Copyright (C) 2007, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* ========================================================================== */
/* === paraklete_allocate_numeric =========================================== */
/* ========================================================================== */

/* Allocate initial part of LU factors */

static paraklete_numeric *amesos_paraklete_allocate_numeric
(
    paraklete_symbolic *LUsymbolic,
    paraklete_common *Common
)
{
    paraklete_numeric *LU ;
    paraklete_node *LUnode ;
    cholmod_common *cm ;
    Int *Cstart, *Sched, *Childp ;
    Int c, n, ncomponents, myid ;

    cm = &(Common->cm) ;
    myid = Common->myid ;
    LU = CHOLMOD (malloc) (1, sizeof (paraklete_numeric), cm) ;

    if (cm->status != CHOLMOD_OK)
    {
        /* TODO return NULL, and broadcast error to all processes */
        PARAKLETE_ERROR (PK_OUT_OF_MEMORY, "out of memory") ;
    }

    n = LUsymbolic->n ;
    ncomponents = LUsymbolic->ncomponents ;

    LU->magic = PARAKLETE_MAGIC ;
    LU->n = n ;
    LU->ncomponents = ncomponents ;
    LU->LUnode = CHOLMOD (malloc) (ncomponents, sizeof (paraklete_node *), cm) ;
    LU->P = CHOLMOD (malloc) (n, sizeof (Int), cm) ;
    LU->Q = CHOLMOD (malloc) (n, sizeof (Int), cm) ;
    LU->Pinv = CHOLMOD (malloc) (n, sizeof (Int), cm) ;
    LU->Qinv = CHOLMOD (malloc) (n, sizeof (Int), cm) ;
    LU->W = NULL ;
    LU->Ep2 = NULL ;
    LU->E = NULL ;

    if (myid == 0)
    {
	/* allocate workspace for subsequent solve */
	LU->W = CHOLMOD (malloc) (n, sizeof (double), cm) ;
	/* allocate workspace for distributing the input matrix to the nodes */
	LU->Ep2 = CHOLMOD (malloc) (n+1, sizeof (Int), cm) ;
    }

    if (LU->LUnode != NULL)
    {
	for (c = 0 ; c < ncomponents ; c++)
	{
	    LU->LUnode [c] = CHOLMOD (malloc) (1, sizeof (paraklete_node), cm) ;
	}
    }

    Cstart = LUsymbolic->Cstart ;
    Sched = LUsymbolic->Sched ;
    Childp = LUsymbolic->Childp ;

    if (cm->status != CHOLMOD_OK)
    {
        /* TODO return NULL, and broadcast error to all processes, and
         * free the LU object */
        PARAKLETE_ERROR (PK_OUT_OF_MEMORY, "out of memory") ;
    }

    for (c = 0 ; c < ncomponents ; c++)
    {
	/* Each process has an LUnode [c] for each node in the tree, but it
	 * will be populated only with the parts this process needs */
	LUnode = (LU->LUnode == NULL) ? NULL : (LU->LUnode [c]) ;

	if (LUnode != NULL)
	{

	    LUnode->nk = Cstart [c+1] - Cstart [c] ;
	    LUnode->nchild = Childp [c+1] - Childp [c] ;

	    /* no LU factorization of this node yet */
	    LUnode->PK_STATUS = PK_UNKNOWN ;
	    LUnode->PK_NN = 0 ;
	    LUnode->PK_NPIV = 0 ;
	    LUnode->PK_NFOUND = 0 ;
	    LUnode->PK_NLOST = 0 ;
	    LUnode->lusize = 0 ;
	    LUnode->llen = NULL ;
	    LUnode->lp = NULL ;
	    LUnode->ulen = NULL ;
	    LUnode->up = NULL ;
	    LUnode->ix = NULL ;

	    /* information and messages from each child of node c */
	    if (Sched [c] == myid)
	    {
		LUnode->Lost = CHOLMOD (malloc) (LUnode->nchild,
			sizeof (Int), cm) ;
		LUnode->Lostp = CHOLMOD (malloc) (LUnode->nchild+1,
			sizeof (Int), cm);
		MPI (LUnode->Req = CHOLMOD (malloc) (LUnode->nchild,
			sizeof (MPI_Request), cm)) ;
	    }
	    else
	    {
		LUnode->Lost = NULL ;
		LUnode->Lostp = NULL ;
		MPI (LUnode->Req = NULL) ;
	    }

	    /* no permutation vectors yet */
	    LUnode->Pglobal = NULL ;
	    LUnode->Qglobal = NULL ;
	    LUnode->Plocal = NULL ;
	    LUnode->Qlocal = NULL ;
	    LUnode->Pinv = NULL ;
	    LUnode->Qinv = NULL ;

	    /* no Schur complement yet */
	    LUnode->PK_SSIZE = 0 ;
	    LUnode->PK_SNZ = 0 ;
	    LUnode->PK_SN = 0 ;
	    LUnode->slen = NULL ;
	    LUnode->sp = NULL ;
	    LUnode->sx = NULL ;

	    /* no solution and right-hand-side yet */
	    LUnode->B = NULL ;
	    LUnode->X = NULL ;

	    /* no input matrix yet */
	    LUnode->A = NULL ;
	    LUnode->C = NULL ;

	    /* no workspace yet */
	    LUnode->W2 = NULL ;
	}
    }

    if (cm->status != CHOLMOD_OK)
    {
        /* TODO return NULL, and broadcast error to all processes, and */
        /* free the LU object */
        PARAKLETE_ERROR (PK_OUT_OF_MEMORY, "out of memory") ;
    }

    PR1 ((Common->file, "proc "ID" LU done\n", myid)) ;
    return (LU) ;
}


/* ========================================================================== */
/* === paraklete_free_numeric =============================================== */
/* ========================================================================== */

/* Free the numeric object on all processors */

void amesos_paraklete_free_numeric
(
    paraklete_numeric **LUHandle,
    paraklete_common *Common
)
{
    paraklete_numeric *LU ;
    paraklete_node *LUnode ;
    cholmod_common *cm ;
    Int c ;

    if (LUHandle == NULL)
    {
	/* nothing to do */
	return ;
    }
    LU = *LUHandle ;
    if (LU == NULL)
    {
	/* nothing to do */
	return ;
    }

    cm = &(Common->cm) ;

    /* global P and Q, broadcast to all processors */
    CHOLMOD (free) (LU->n, sizeof (Int), LU->P,    cm) ;
    CHOLMOD (free) (LU->n, sizeof (Int), LU->Q,    cm) ;
    CHOLMOD (free) (LU->n, sizeof (Int), LU->Pinv, cm) ;
    CHOLMOD (free) (LU->n, sizeof (Int), LU->Qinv, cm) ;
    CHOLMOD (free) (LU->n, sizeof (double), LU->W,    cm) ;
    CHOLMOD (free) (LU->n+1, sizeof (Int), LU->Ep2,  cm) ;
    CHOLMOD (free_sparse) (&(LU->E), cm) ;

    if (LU->LUnode != NULL)
    {
	for (c = 0 ; c < LU->ncomponents ; c++)
	{
	    LUnode = LU->LUnode [c] ;
	    if (LUnode != NULL)
	    {
		/* solution and right-hand-side at this node */
		PR2 ((Common->file,"proc "ID" node "ID" free numeric, nk "ID" B %p\n",
			Common->myid, c, LUnode->nk, (void *) (LUnode->B))) ;
		CHOLMOD (free) (LUnode->nk, sizeof (double), LUnode->B, cm) ;
		CHOLMOD (free) (LUnode->PK_NN, sizeof (double), LUnode->X, cm) ;

		/* LU factors at this node */
		CHOLMOD (free) (LUnode->PK_NPIV, sizeof (Int), LUnode->llen, cm) ;
		CHOLMOD (free) (LUnode->PK_NPIV, sizeof (Int), LUnode->lp, cm) ;
		CHOLMOD (free) (LUnode->PK_NN, sizeof (Int), LUnode->ulen, cm) ;
		CHOLMOD (free) (LUnode->PK_NN, sizeof (Int), LUnode->up, cm) ;
		CHOLMOD (free) (LUnode->lusize, sizeof (double), LUnode->ix, cm) ;
		CHOLMOD (free) (LUnode->nchild, sizeof (Int), LUnode->Lost, cm) ;
		CHOLMOD (free) (LUnode->nchild+1, sizeof (Int),
			LUnode->Lostp, cm) ;
		MPI (CHOLMOD (free) (LUnode->nchild,
			    sizeof (MPI_Request), LUnode->Req, cm)) ;

		/* P and Q at this node */
		CHOLMOD (free) (LUnode->PK_NPIV, sizeof (Int),
			LUnode->Pglobal, cm) ;
		CHOLMOD (free) (LUnode->PK_NPIV, sizeof (Int),
			LUnode->Qglobal, cm) ;
		CHOLMOD (free) (LUnode->PK_NPIV, sizeof (Int),
			LUnode->Plocal, cm) ;
		CHOLMOD (free) (LUnode->PK_NPIV, sizeof (Int),
			LUnode->Qlocal, cm) ;
		CHOLMOD (free) (LUnode->PK_NPIV, sizeof (Int), LUnode->Pinv, cm) ;
		CHOLMOD (free) (LUnode->PK_NPIV, sizeof (Int), LUnode->Qinv, cm) ;

		/* Schur complement of this node */
		CHOLMOD (free) (LUnode->PK_SN, sizeof (Int), LUnode->sp, cm) ;
		CHOLMOD (free) (LUnode->PK_SN, sizeof (Int), LUnode->slen, cm) ;
		CHOLMOD (free) (LUnode->PK_SSIZE, sizeof (double),
			LUnode->sx, cm) ;

		/* input matrix and sum of Schur complements at this node */
		CHOLMOD (free_sparse) (&(LUnode->A), cm) ;
		CHOLMOD (free_sparse) (&(LUnode->C), cm) ;

		CHOLMOD (free) (2*(LUnode->nlost_in), sizeof (Int),
			LUnode->W2, cm) ;

		/* free the LUnode itself */
		CHOLMOD (free) (1, sizeof (paraklete_node), LUnode, cm) ;
	    }
	}
    }

    CHOLMOD (free) (LU->ncomponents, sizeof (paraklete_node *), LU->LUnode, cm) ;
    *LUHandle = CHOLMOD (free) (1, sizeof (paraklete_numeric), (*LUHandle), cm) ;
}


/* ========================================================================== */
/* === paraklete_permute ==================================================== */
/* ========================================================================== */

/* E = A(p,q) */

static cholmod_sparse *paraklete_permute
(
    cholmod_sparse *A,
    paraklete_numeric *LU,
    paraklete_symbolic *LUsymbolic,
    paraklete_common *Common
)
{
    double *Ax, *Ex ;
    Int *Ap, *Ai, *Ep2, *Ep, *Ei, *Cperm, *RpermInv, *Rperm ;
    cholmod_common *cm ;
    cholmod_sparse *E ;
    Int n, enz, anz, k, j, p ;

    cm = &(Common->cm) ;
    ASSERT (Common->myid == 0) ;

    Cperm = LUsymbolic->Cperm ;
    Rperm = LUsymbolic->Rperm ;
    Rperm = (Rperm == NULL) ? Cperm : Rperm ;
    RpermInv = LUsymbolic->RpermInv ;

    n = A->nrow ;
    ASSERT (n == LUsymbolic->n) ;
    Ap = A->p ;
    Ai = A->i ;
    Ax = A->x ;
    anz = Ap [n] ;
    Ep2 = LU->Ep2 ;

#ifndef NDEBUG
    for (k = 0 ; k < n ; k++)
    {
        PR3 ((Common->file, "Cperm ("ID") = "ID" ;\n", k+1, Cperm[k]+1));
    }
    for (k = 0 ; k < n ; k++)
    {
        PR3 ((Common->file, "Rperm ("ID") = "ID" ;\n", k+1, Rperm[k]+1));
    }
    for (k = 0 ; k < n ; k++)
    {
        PR3 ((Common->file, "RpermInv ("ID") = "ID" ;\n", k+1, RpermInv [k]+1)) ;
    }
#endif

    E = CHOLMOD (allocate_sparse) (n, n, anz, FALSE, TRUE, 0, TRUE, cm) ;

    if (cm->status != CHOLMOD_OK)
    {
	/* out of memory */
	PR0 ((Common->file, "failed to allocate E\n")) ;
	return (NULL) ;
    }

    Ep = E->p ;
    Ei = E->i ;
    Ex = E->x ;
    enz = 0 ;
    for (k = 0 ; k < n ; k++)
    {
	/* column j of A becomes column k of E */
	j = Cperm [k] ;
	PR1 ((Common->file, "Col "ID" of A becomes col "ID" of E\n", j, k)) ;
	Ep [k] = enz ;
	for (p = Ap [j] ; p < Ap [j+1] ; p++)
	{
	    /* row i = Ai [p] becomes row RpermInv[i] of E */
	    PR2 ((Common->file, "    "ID" %g -> "ID"\n",
                Ai [p], Ax [p], RpermInv [Ai [p]])) ;
	    Ei [enz] = RpermInv [Ai [p]] ;
	    Ex [enz] = Ax [p] ;
	    enz++ ;
	    ASSERT (enz <= anz) ;
	}
    }
    Ep [n] = enz ;
    DEBUG (CHOLMOD (print_sparse) (E, "E = A(p,p)", cm)) ;
    CHOLMOD (sort) (E, cm) ;

    if (cm->status != CHOLMOD_OK)
    {
	/* out of memory */
	PR0 ((Common->file, "out of memory to sort E\n")) ;
	CHOLMOD (free_sparse) (&E, cm) ;
	return (NULL) ;
    }

    DEBUG (CHOLMOD (print_sparse) (E, "E = A(p,p), sorted", cm)) ;
    for (k = 0 ; k <= n ; k++)
    {
	Ep2 [k] = Ep [k] ;
    }

    return (E) ;
}


/* ========================================================================== */
/* === paraklete_factorize ================================================== */
/* ========================================================================== */

paraklete_numeric *amesos_paraklete_factorize
(
    /* inputs, not modified */
    cholmod_sparse *A,
    paraklete_symbolic *LUsymbolic,
    paraklete_common *Common
)
{
    paraklete_numeric *LU ;
    cholmod_common *cm ;
    cholmod_sparse *E, *C ;
    double *Ex, *Cx ;
    Int *Cperm, *Cn, *Cnz, *Ep, *Ei, *Ep2, *Map, *Cparent, *Rperm,
	*Cstart, *P, *Q, *Cp, *Ci, *Pc, *Qc, *Pinv, *Qinv, *Sched ;
    Int cj, i, n, ncomponents, k, p, c, a, cn, k1, k2, kk, cnz,
	nfound, myid, npiv ;
    int ok, all_ok ;

    MPI (MPI_Status ms) ;
    MPI (MPI_Request req) ;

    /* TODO: return NULL if any input argument is NULL */

    Common->status = PK_OK ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    cm = &(Common->cm) ;
    cm->status = CHOLMOD_OK ;

    myid = Common->myid ;

    PR0 ((Common->file, "proc "ID" start factorize\n", myid)) ;

    n = LUsymbolic->n ;
    /*
    printf ("   factorize n "ID" LUsymbolic %p\n", n, (void *) LUsymbolic) ;
    */
    ncomponents = LUsymbolic->ncomponents ;
    Cperm   = LUsymbolic->Cperm ;
    Rperm   = LUsymbolic->Rperm ;
    Rperm   = (Rperm == NULL) ? Cperm : Rperm ;
    Cn      = LUsymbolic->Cn ;
    Cnz     = LUsymbolic->Cnz ;
    Cparent = LUsymbolic->Cparent ;
    Cstart  = LUsymbolic->Cstart  ;
    Sched   = LUsymbolic->Sched ;

    PR0 ((Common->file, "in factor: proc "ID" my_tries "ID"\n", myid, my_tries)) ;

#ifndef NDEBUG
    for (c = 0 ; c < ncomponents ; c++)
    {
	PR1 ((Common->file, "proc: "ID" node "ID": "ID" to "ID", parent "ID" Sched "ID"\n",
		Common->myid, c, Cstart [c], Cstart [c+1]-1, Cparent [c],
		Sched [c])) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* allocate and initialize the global LU factors */
    /* ---------------------------------------------------------------------- */

    LU = amesos_paraklete_allocate_numeric (LUsymbolic, Common) ;
    ok = (LU != NULL) ;
    PR1 ((Common->file, "factor proc "ID" alloc ok: "ID"\n", myid, ok)) ;

    /* ---------------------------------------------------------------------- */
    /* E = A(p,p), with sorted column indices */
    /* ---------------------------------------------------------------------- */

    E = NULL ;
    Ep2 = NULL ;
    if (ok && myid == 0)
    {
	/* only the root constructs the matrix E = A (p,p) */
	Ep2 = LU->Ep2 ;
	E = LU->E = paraklete_permute (A, LU, LUsymbolic, Common) ;
	ok = (E != NULL) ;
	DEBUG (CHOLMOD (print_sparse) (E, "E = A(p,p) on master node", cm)) ;
    }

    /* ---------------------------------------------------------------------- */
    /* allocate space for submatrices at each node */
    /* ---------------------------------------------------------------------- */

    for (c = 0 ; ok && c < ncomponents ; c++)
    {
	if (myid == 0 || Sched [c] == myid)
	{
	    /* Two copies are made, one on the root process and one on the
	     * process that factorizes that submatrix.  If the root process
	     * is factorizing the node, then only one copy is made. */
	    C = CHOLMOD (allocate_sparse) (Cn [c], Cn [c], Cnz [c], TRUE, TRUE,
		    0, TRUE, cm) ;
	    LU->LUnode [c]->A = C ;
	    ok = (C != NULL) ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* all processes find out if initial allocation fails for any process */
    /* ---------------------------------------------------------------------- */

    all_ok = ok ;
    MPI (MPI_Allreduce (&ok, &all_ok, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD)) ;
    if (!all_ok)
    {
	/* out of memory; inform all processes */
	PR0 ((Common->file, "proc "ID" all fail factorize\n", Common->myid)) ;
	Common->status = PK_OUT_OF_MEMORY ;
	amesos_paraklete_free_numeric (&LU, Common) ;
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* distribute the permuted matrix to the nodes */
    /* ---------------------------------------------------------------------- */

    if (myid == 0)
    {

	/* ------------------------------------------------------------------ */
	/* process 0 partitions the input matrix and sends to all processes */
	/* ------------------------------------------------------------------ */

	Ep = E->p ;
	Ei = E->i ;
	Ex = E->x ;

#ifndef NDEBUG
	for (k = 0 ; k < n ; k++)
	{
	    for (p = Ep [k] ; p < Ep [k+1] ; p++)
	    {
		PR3 ((Common->file, "E ("ID","ID") = %g ;\n", 1+Ei[p], 1+k, Ex[p]));
	    }
	}
	for (c = 0 ; c <= ncomponents ; c++)
        {
	    PR3 ((Common->file, "Cstart ("ID") = "ID" ;\n", 1+c, 1+Cstart [c])) ;
        }
#endif

	Map = cm->Iwork ;

	for (c = 0 ; c < ncomponents ; c++)
	{
	    PR1 ((Common->file, "distribute node c = "ID"\n", c)) ;
	    cn = Cn [c] ;

	    /* find mapping of global rows/columns to node c's local indices */
	    cj = 0 ;
	    for (a = c ; a != TRILINOS_CHOLMOD_EMPTY ; a = Cparent [a])
	    {
		PR2 ((Common->file, "ancestor "ID", for c "ID", ncomp "ID"\n",
			a, c, ncomponents)) ;
		ASSERT (a >= c && a < ncomponents) ;
		k1 = Cstart [a] ;
		k2 = Cstart [a+1] ;
		PR2 ((Common->file, "k1 "ID" k2 "ID"\n", k1, k2)) ;
		for (k = k1 ; k < k2 ; k++)
		{
		    /* global index k becomes local index j */
		    PR3 ((Common->file, "  global: "ID" local "ID"\n", k, cj)) ;
		    Map [k] = cj++ ;
		}
	    }
	    ASSERT (cj == cn) ;

	    /* get the local matrix for node c */
	    C = LU->LUnode [c]->A ;
	    Cp = C->p ;
	    Ci = C->i ;
	    Cx = C->x ;

	    cj = 0 ;
	    cnz = 0 ;

	    /* create columns 0:(k2-k1-1) of C, containing candidate pivot
	     * columns for node c */
	    k1 = Cstart [c] ;
	    k2 = Cstart [c+1] ;
	    PR2 ((Common->file, "c: "ID" k1 to k2-1: "ID" "ID"\n", c, k1, k2-1)) ;
	    for (k = k1 ; k < k2 ; k++)
	    {
		Cp [cj++] = cnz ;
		for (p = Ep2 [k] ; p < Ep [k+1] ; p++)
		{
		    i = Ei [p] ;
		    ASSERT (Ei [p] >= k1) ;
		    PR3 ((Common->file,
                        "E(1+"ID",1+"ID") = %g ; ce(1+"ID",1+"ID") = "ID"\n",
		        i,k, Ex[p], i,k,c)) ;
		    Ci [cnz] = Map [i] ;
		    Cx [cnz] = Ex [p] ;
		    cnz++ ;
		    ASSERT (cnz <= Cnz [c]) ;
		}
	    }

	    /* create columns for the ancestors of C.  These will become columns
	     * of U2 and S for node c. */
	    for (a = Cparent [c] ; a != TRILINOS_CHOLMOD_EMPTY ; a = Cparent [a])
	    {
		PR2 ((Common->file, "ancestor "ID"\n", a)) ;
		PR2 ((Common->file, "k1 "ID" k2 "ID"\n", Cstart [a], Cstart [a+1]));
		for (k = Cstart [a] ; k < Cstart [a+1] ; k++)
		{
		    Cp [cj++] = cnz ;
		    ASSERT (cj <= cn) ;
		    for (p = Ep2 [k] ; p < Ep [k+1] ; p++)
		    {
			i = Ei [p] ;
			ASSERT (i >= k1) ;
			if (i >= k2)
			{
			    /* only get rows Cstart [c] to Cstart [c+1]-1 */
			    break ;
			}
			PR3 ((Common->file,
                            "E(1+"ID",1+"ID") = %g ; ce(1+"ID",1+"ID") = "ID"\n",
			    i,k, Ex[p], i,k,c)) ;
			Ci [cnz] = Map [i] ;
			Cx [cnz] = Ex [p] ;
			cnz++ ;
			ASSERT (cnz <= Cnz [c]) ;
		    }
		    /* next component will start here */
		    Ep2 [k] = p ;
		}
	    }
	    ASSERT (cj == cn) ;
	    ASSERT (cnz == Cnz [c]) ;
	    Cp [cn] = cnz ;

	    /* place matrix C in node c of the tree */
	    PR2 ((Common->file, "node: "ID" sched: "ID"\n", c, Sched [c])) ;

	    if (Sched [c] != myid)
	    {
		/* send this matrix to the process that owns node c.
		   Note that the Isend requests are immediately free'd,
		   because the sends are synchronized with an barrier,
		   later on. */

		/*
		if (cnz == 4040)
		{
		    Int k3 = 1084 ;
		    fprintf (Common->file,
			"sending "ID": ["ID" %g]\n", k3, Ci [k3], Cx [k3]) ;
		}
		*/

		PR2 ((Common->file, "n "ID" nz "ID"\n", cn, cnz)) ;
		DEBUG (CHOLMOD (print_sparse) (C, "send C to a node", cm)) ;

		MPI (MPI_Isend (Cp, cn+1, MPI_Int, Sched [c],
			    TAG0, MPI_COMM_WORLD, &req)) ;
		MPI (MPI_Request_free (&req)) ;
		MPI (MPI_Isend (Ci, cnz,  MPI_Int, Sched [c],
			    TAG0, MPI_COMM_WORLD, &req)) ;
		MPI (MPI_Request_free (&req)) ;
		MPI (MPI_Isend (Cx, cnz,  MPI_DOUBLE, Sched [c],
			    TAG0, MPI_COMM_WORLD, &req)) ;
		MPI (MPI_Request_free (&req)) ;
	    }
	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* all other processes acquire their input matrix from process 0 */
	/* ------------------------------------------------------------------ */

	for (c = 0 ; c < ncomponents ; c++)
	{
	    if (Sched [c] == myid)
	    {
		/*
		{
		    Int *Mi = LU->LUnode [c]->A->i ;
		    double *Mx = LU->LUnode [c]->A->x ;
		    Mi [1084] = -9999999 ;
		    Mx [1084] = -9999999 ;
		    fprintf (Common->file,
			"pre ["ID" %g]\n", Mi [1084], Mx [1084]) ;
		}
		*/

		MPI (MPI_Recv (LU->LUnode [c]->A->p, Cn [c] +1, MPI_Int, 0,
			    TAG0, MPI_COMM_WORLD, &ms)) ;
		MPI (MPI_Recv (LU->LUnode [c]->A->i, Cnz [c], MPI_Int, 0,
			    TAG0, MPI_COMM_WORLD, &ms)) ;
		MPI (MPI_Recv (LU->LUnode [c]->A->x, Cnz [c], MPI_DOUBLE, 0,
			    TAG0, MPI_COMM_WORLD, &ms)) ;

		/*
		{
		    Int *Mi = LU->LUnode [c]->A->i ;
		    double *Mx = LU->LUnode [c]->A->x ;
		    char *CC = (char *) (LU->LUnode [c]->A->x) ;
		    Int k3 = 1084 ;
		    fprintf (Common->file,
			"got "ID" ["ID" %g]\n", k3, Mi [k3], Mx [k3]) ;
		    fprintf (Common->file, "byte "ID"\n", (Int) (CC [8*k3])) ;
		    fprintf (Common->file, "byte "ID"\n", (Int) (CC [8*k3+1])) ;
		    fprintf (Common->file, "byte "ID"\n", (Int) (CC [8*k3+2])) ;
		    fprintf (Common->file, "byte "ID"\n", (Int) (CC [8*k3+3])) ;
		    fprintf (Common->file, "byte "ID"\n", (Int) (CC [8*k3+4])) ;
		    fprintf (Common->file, "byte "ID"\n", (Int) (CC [8*k3+5])) ;
		    fprintf (Common->file, "byte "ID"\n", (Int) (CC [8*k3+6])) ;
		    fprintf (Common->file, "byte "ID"\n", (Int) (CC [8*k3+7])) ;
		}
		*/

		PR1 ((Common->file, "proc "ID" Node "ID" got orig:\n", myid, c));
		PR2 ((Common->file, "n "ID" nz "ID"\n", Cn [c], Cnz [c])) ;
		DEBUG (CHOLMOD (print_sparse) (LU->LUnode [c]->A, "got A", cm)) ;
#ifndef NDEBUG
		/* {
		    Int *Mp, *Mi, j ;
		    double *Mx ;
		    Mp = LU->LUnode [c]->A->p ;
		    Mi = LU->LUnode [c]->A->i ;
		    Mx = LU->LUnode [c]->A->x ;
		    for (j = 0 ; j < Cn [c] ; j++)
		    {
			fprintf (Common->file, "my own col "ID"\n", j) ;
			for (k = Mp [j] ; k < Mp [j+1] ; k++)
			{
			    Int ii = Mi [k] ;
			    double x = Mx [k] ;
			    fprintf (Common->file, " is: "ID" %g\n", ii, x) ;
			}
		    }
		    k = 0 ;
		} */
#endif
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* free temporary copy of A(p,p) */
    /* ---------------------------------------------------------------------- */

    /* only the root needs to do this */
    if (myid == 0)
    {
	LU->Ep2 = CHOLMOD (free) (n+1, sizeof (Int), LU->Ep2, cm) ;
	CHOLMOD (free_sparse) (&(LU->E), cm) ;
    }

    /* process 0 frees LUnode [c]->A, but only when recvd by c */
    MPI (MPI_Barrier (MPI_COMM_WORLD)) ;
    PR1 ((Common->file, "proc "ID" everybody OK so far2\n", Common->myid)) ;
    if (myid == 0)
    {
	for (c = 0 ; c < ncomponents ; c++)
	{
	    if (Sched [c] != 0)
	    {
		/* root process no longer needs submatrices sent to others */
		CHOLMOD (free_sparse) (&(LU->LUnode [c]->A), cm) ;
	    }
	}
    }

#ifndef NDEBUG
    for (c = 0 ; c < ncomponents ; c++)
    {
	if (Sched [c] == myid)
	{
	    PR1 ((Common->file, "proc "ID" Node "ID" original matrix:\n", myid, c));
	    DEBUG (CHOLMOD (print_sparse) (LU->LUnode [c]->A, "A", cm)) ;
	}
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* factorize each node */
    /* ---------------------------------------------------------------------- */

    ok = TRUE ;
    for (c = 0 ; c < ncomponents ; c++)
    {
	if (Sched [c] == myid)
	{
	    PR1 ((Common->file, "proc "ID" doing node "ID"\n", myid, c)) ;
	    if (!amesos_paraklete_factorize_node (c, LU, LUsymbolic, Common))
	    {
		ok = FALSE ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* determine global ordering, including Cperm and numerical pivoting */
    /* ---------------------------------------------------------------------- */

    P = LU->P ;
    Q = LU->Q ;
    kk = 0 ;
    Common->status = TRUE ;
    for (c = 0 ; c < ncomponents ; c++)
    {
	/* Give all processors nfound and npiv for each node of the tree.
	 * This also allows us to determine if the matrix is singular, or if
	 * any process ran out of memory. */
	MPI (MPI_Bcast (LU->LUnode [c]->header, PK_HEADER, MPI_Int,
		    Sched [c], MPI_COMM_WORLD)) ;
	kk += LU->LUnode [c]->PK_NFOUND ;
	Common->status = MIN (Common->status, LU->LUnode [c]->PK_STATUS) ;
    }

    if (kk < n)
    {
	PR0 ((Common->file, "proc "ID" Singular: "ID" of "ID"\n", myid, kk, n)) ;
	/*
        printf ("proc "ID" Singular: "ID" of "ID"\n", myid, kk, n) ;
        */
	Common->status = PK_SINGULAR ;
    }

    if (Common->status != PK_OK)
    {
	/* out of memory, singular matrix, or unknown error */
	amesos_paraklete_free_numeric (&LU, Common) ;
	return (NULL) ;
    }

    /* all processes allocate space for the global permutation */
    for (c = 0 ; c < ncomponents ; c++)
    {
	npiv = LU->LUnode [c]->PK_NPIV ;
	if (LU->LUnode [c]->Pglobal == NULL)
	{
	    LU->LUnode [c]->Pglobal = CHOLMOD (malloc) (npiv, sizeof (Int), cm) ;
	    LU->LUnode [c]->Qglobal = CHOLMOD (malloc) (npiv, sizeof (Int), cm) ;
	}
    }

    ok = ok && (kk == n) && (cm->status == CHOLMOD_OK) ;

    /* check to see if all processes had space for P and Q */
    PR1 ((Common->file, "proc "ID" here in PQ: ok "ID"\n", Common->myid, ok)) ;
    all_ok = ok ;
    MPI (MPI_Allreduce (&ok, &all_ok, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD)) ;
    if (!all_ok)
    {
        /* TODO return NULL, broadcast error to all processes, and free LU */
        PARAKLETE_ERROR (PK_OUT_OF_MEMORY, "out of memory") ;
    }

    kk = 0 ;
    for (c = 0 ; c < ncomponents ; c++)
    {
	npiv = LU->LUnode [c]->PK_NPIV ;

	/* give Pglobal and Qglobal to all processors */
	MPI (MPI_Bcast (LU->LUnode [c]->Pglobal, npiv, MPI_Int, Sched [c],
		    MPI_COMM_WORLD)) ;
	MPI (MPI_Bcast (LU->LUnode [c]->Qglobal, npiv, MPI_Int, Sched [c],
		    MPI_COMM_WORLD)) ;

	Pc = LU->LUnode [c]->Pglobal ;
	Qc = LU->LUnode [c]->Qglobal ;

	nfound = LU->LUnode [c]->PK_NFOUND ;
	for (k = 0 ; k < nfound ; k++)
	{
	    /* row Pc[k] and column Qc[k] of E are the kth pivot row/column of
	     * node c.  They are the kk-th pivot row/column of the global LU
	     * factors. */
	    P [kk] = Rperm [Pc [k]] ;
	    Q [kk] = Cperm [Qc [k]] ;
	    kk++ ;
	}
    }
    ASSERT (kk == n) ;

    /* compute Pinv and Qinv.  TODO: this is not needed */
    Pinv = LU->Pinv ;
    Qinv = LU->Qinv ;
    for (k = 0 ; k < n ; k++)
    {
	Pinv [P [k]] = k ;
	Qinv [Q [k]] = k ;
    }

#ifndef NDEBUG

    /* ---------------------------------------------------------------------- */
    /* dump the matrix */
    /* ---------------------------------------------------------------------- */

    if (Common->dump > 1)
    {

	if (myid == 0)
	{
	    Int *Ap, *Ai ;
	    double *Ax ;
	    Int j ;

	    Ap = A->p ;
	    Ai = A->i ;
	    Ax = A->x ;

	    for (k = 0 ; k < n ; k++) printf ("P (1+"ID") = "ID" ;\n", k, P [k]) ;
	    for (k = 0 ; k < n ; k++) printf ("Q (1+"ID") = "ID" ;\n", k, Q [k]) ;

	    printf ("P = 1+P ;\n") ;
	    printf ("Q = 1+Q ;\n") ;

	    for (j = 0 ; j < n ; j++)
	    {
		for (p = Ap [j] ; p < Ap [j+1] ; p++)
		{
		    printf ("A (1+"ID",1+"ID") = %.16g ;\n", Ai [p], j, Ax [p]) ;
		}
	    }
	}

	for (c = 0 ; c < ncomponents ; c++)
	{
	    paraklete_node *LUnode ;
	    Int *Lip, *Llen, *Li, *Uip, *Ulen, *Ui ;
	    double *LUix, *Lx, *Ux ;
	    Int llen, ulen, j ;

	    MPI (MPI_Barrier (MPI_COMM_WORLD)) ;
	    if (Sched [c] != myid) continue ;
	    printf ("\n%% ---- node "ID"\n", c) ;

	    /* dump L */
	    LUnode = LU->LUnode [c] ;
	    Lip = LUnode->lp ;
	    Llen = LUnode->llen ;
	    LUix = LUnode->ix ;
	    nfound = LUnode->PK_NFOUND ;
	    for (j = 0 ; j < nfound ; j++)
	    {
		GET_COLUMN (Lip, Llen, LUix, j, Li, Lx, llen) ;
		printf ("\nL (1+"ID",1+"ID") = 1 ;\n", j, j) ;
		for (p = 1 ; p < llen ; p++)
		{
		    printf ("L (1+"ID",1+"ID") = %.16g ;\n", Li [p], j, Lx [p]) ;
		}
	    }

	    /* dump U */
	    cn = LUnode->PK_NN ;
	    Uip = LUnode->up ;
	    Ulen = LUnode->ulen ;
	    for (j = cn-1 ; j >= nfound ; j--)
	    {
		printf ("\n") ;
		GET_COLUMN (Uip, Ulen, LUix, j, Ui, Ux, ulen) ;
		for (p = 0 ; p < ulen ; p++)
		{
		    printf ("U (1+"ID",1+"ID") = %.16g ;\n", Ui [p], j, Ux [p]) ;
		}
	    }
	    for ( ; j >= 0 ; j--)
	    {
		GET_COLUMN (Uip, Ulen, LUix, j, Ui, Ux, ulen) ;
		printf ("\nU (1+"ID",1+"ID") = %.16g ; %% pivot\n",
		    j, j, Ux [ulen-1]) ;
		for (p = 0 ; p < ulen-1 ; p++)
		{
		    printf ("U (1+"ID",1+"ID") = %.16g ;\n", Ui [p], j, Ux [p]) ;
		}
	    }
	}

	if (Common->nproc == 1 && Common->dump > 1)
	{
	    printf ("norm (L*U-A(P,Q))\n") ;
	}
    }
#endif

    return (LU) ;
}
