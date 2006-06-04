/* ========================================================================== */
/* === paraklete_factorize ================================================== */
/* ========================================================================== */

#include "paraklete.h"

/* LU = paraklete_factorize (A, LUsymbolic, Common) factorizes P*A*Q into L*U.
 * Returns NULL if A is singular or if memory is exhausted.
 *
 * PARAKLETE version 0.1: parallel sparse LU factorization.  May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* ========================================================================== */
/* === paraklete_allocate_numeric =========================================== */
/* ========================================================================== */

/* Allocate initial part of LU factors */

static paraklete_numeric *paraklete_allocate_numeric
(
    paraklete_symbolic *LUsymbolic,
    paraklete_common *Common
)
{
    paraklete_numeric *LU ;
    paraklete_node *LUnode ;
    cholmod_common *cm ;
    int *Cstart, *Sched, *Childp ;
    int c, n, ncomponents, myid ;

    cm = &(Common->cm) ;
    myid = Common->myid ;
    LU = cholmod_malloc (1, sizeof (paraklete_numeric), cm) ;

    if (cm->status != CHOLMOD_OK)
    {
	/* out of memory */
	PR0 ((Common->file, "proc %d LU header failed\n", myid)) ;
	return (NULL) ;
    }

    n = LUsymbolic->n ;
    ncomponents = LUsymbolic->ncomponents ;

    LU->n = n ;
    LU->ncomponents = ncomponents ;
    LU->LUnode = cholmod_malloc (ncomponents, sizeof (paraklete_node *), cm) ;
    LU->P = cholmod_malloc (n, sizeof (int), cm) ;
    LU->Q = cholmod_malloc (n, sizeof (int), cm) ;
    LU->Pinv = cholmod_malloc (n, sizeof (int), cm) ;
    LU->Qinv = cholmod_malloc (n, sizeof (int), cm) ;
    LU->W = NULL ;
    LU->Ep2 = NULL ;
    LU->E = NULL ;
    if (myid == 0)
    {
	/* allocate workspace for subsequent solve */
	LU->W = cholmod_malloc (n, sizeof (double), cm) ;
	/* allocate workspace for distributing the input matrix to the nodes */
	LU->Ep2 = cholmod_malloc (n+1, sizeof (int), cm) ;
    }

    if (LU->LUnode != NULL)
    {
	for (c = 0 ; c < ncomponents ; c++)
	{
	    LU->LUnode [c] = cholmod_malloc (1, sizeof (paraklete_node), cm) ;
	}
    }

    Cstart = LUsymbolic->Cstart ;
    Sched = LUsymbolic->Sched ;
    Childp = LUsymbolic->Childp ;

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
		LUnode->Lost = cholmod_malloc (LUnode->nchild,
			sizeof (int), cm) ;
		LUnode->Lostp = cholmod_malloc (LUnode->nchild+1,
			sizeof (int), cm);
		MPI (LUnode->Req = cholmod_malloc (LUnode->nchild,
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
	/* out of memory */
	PR0 ((Common->file, "proc %d LU contents failed\n", myid)) ;
	paraklete_free_numeric (&LU, Common) ;
    }

    PR1 ((Common->file, "proc %d LU ok\n", myid)) ;
    return (LU) ;
}


/* ========================================================================== */
/* === paraklete_free_numeric =============================================== */
/* ========================================================================== */

/* Free the numeric object on all processors */

void paraklete_free_numeric
(
    paraklete_numeric **LUHandle,
    paraklete_common *Common
)
{
    paraklete_numeric *LU ;
    paraklete_node *LUnode ;
    cholmod_common *cm ;
    int c ;

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
    cholmod_free (LU->n, sizeof (int), LU->P,    cm) ;
    cholmod_free (LU->n, sizeof (int), LU->Q,    cm) ;
    cholmod_free (LU->n, sizeof (int), LU->Pinv, cm) ;
    cholmod_free (LU->n, sizeof (int), LU->Qinv, cm) ;
    cholmod_free (LU->n, sizeof (double), LU->W,    cm) ;
    cholmod_free (LU->n+1, sizeof (int), LU->Ep2,  cm) ;
    cholmod_free_sparse (&(LU->E), cm) ;

    if (LU->LUnode != NULL)
    {
	for (c = 0 ; c < LU->ncomponents ; c++)
	{
	    LUnode = LU->LUnode [c] ;
	    if (LUnode != NULL)
	    {
		/* solution and right-hand-side at this node */
		PR2 ((Common->file, "proc %d node %d free numeric, nk %d B %p\n",
			Common->myid, c, LUnode->nk, (void *) (LUnode->B))) ;
		cholmod_free (LUnode->nk, sizeof (double), LUnode->B, cm) ;
		cholmod_free (LUnode->PK_NN, sizeof (double), LUnode->X, cm) ;

		/* LU factors at this node */
		cholmod_free (LUnode->PK_NPIV, sizeof (int), LUnode->llen, cm) ;
		cholmod_free (LUnode->PK_NPIV, sizeof (int), LUnode->lp, cm) ;
		cholmod_free (LUnode->PK_NN, sizeof (int), LUnode->ulen, cm) ;
		cholmod_free (LUnode->PK_NN, sizeof (int), LUnode->up, cm) ;
		cholmod_free (LUnode->lusize, sizeof (double), LUnode->ix, cm) ;
		cholmod_free (LUnode->nchild, sizeof (int), LUnode->Lost, cm) ;
		cholmod_free (LUnode->nchild+1, sizeof (int),
			LUnode->Lostp, cm) ;
		MPI (cholmod_free (LUnode->nchild,
			    sizeof (MPI_Request), LUnode->Req, cm)) ;

		/* P and Q at this node */
		cholmod_free (LUnode->PK_NPIV, sizeof (int),
			LUnode->Pglobal, cm) ;
		cholmod_free (LUnode->PK_NPIV, sizeof (int),
			LUnode->Qglobal, cm) ;
		cholmod_free (LUnode->PK_NPIV, sizeof (int),
			LUnode->Plocal, cm) ;
		cholmod_free (LUnode->PK_NPIV, sizeof (int),
			LUnode->Qlocal, cm) ;
		cholmod_free (LUnode->PK_NPIV, sizeof (int), LUnode->Pinv, cm) ;
		cholmod_free (LUnode->PK_NPIV, sizeof (int), LUnode->Qinv, cm) ;

		/* Schur complement of this node */
		cholmod_free (LUnode->PK_SN, sizeof (int), LUnode->sp, cm) ;
		cholmod_free (LUnode->PK_SN, sizeof (int), LUnode->slen, cm) ;
		cholmod_free (LUnode->PK_SSIZE, sizeof (double),
			LUnode->sx, cm) ;

		/* input matrix and sum of Schur complements at this node */
		cholmod_free_sparse (&(LUnode->A), cm) ;
		cholmod_free_sparse (&(LUnode->C), cm) ;

		cholmod_free (2*(LUnode->nlost_in), sizeof (int),
			LUnode->W2, cm) ;

		/* free the LUnode itself */
		cholmod_free (1, sizeof (paraklete_node), LUnode, cm) ;
	    }
	}
    }

    cholmod_free (LU->ncomponents, sizeof (paraklete_node *), LU->LUnode, cm) ;
    *LUHandle = cholmod_free (1, sizeof (paraklete_numeric), (*LUHandle), cm) ;
}


/* ========================================================================== */
/* === paraklete_permute ==================================================== */
/* ========================================================================== */

/* E = A(p,p) */

static cholmod_sparse *paraklete_permute
(
    cholmod_sparse *A,
    paraklete_numeric *LU,
    paraklete_symbolic *LUsymbolic,
    paraklete_common *Common
)
{
    double *Ax, *Ex ;
    int *Ap, *Ai, *Ep2, *Ep, *Ei, *Cperm, *Cinv ;
    cholmod_common *cm ;
    cholmod_sparse *E ;
    int n, enz, anz, k, j, p ;

    cm = &(Common->cm) ;
    ASSERT (Common->myid == 0) ;

    Cperm = LUsymbolic->Cperm ;
    Cinv = LUsymbolic->Cinv ;

    n = A->nrow ;
    ASSERT (n == LUsymbolic->n) ;
    Ap = A->p ;
    Ai = A->i ;
    Ax = A->x ;
    anz = Ap [n] ;
    Ep2 = LU->Ep2 ;

    E = cholmod_allocate_sparse (n, n, anz, FALSE, TRUE, 0, TRUE, cm) ;

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
	PRINT1 (("Col %d of A becomes col %d of E\n", j, k)) ;
	Ep [k] = enz ;
	for (p = Ap [j] ; p < Ap [j+1] ; p++)
	{
	    /* row i = Ai [p] becomes row Cinv[i] of E */
	    PRINT2 (("    %d %g -> %d\n", Ai [p], Ax [p], Cinv [Ai [p]])) ;
	    Ei [enz] = Cinv [Ai [p]] ;
	    Ex [enz] = Ax [p] ;
	    enz++ ;
	    ASSERT (enz <= anz) ;
	}
    }
    Ep [n] = enz ;
    DEBUG (cholmod_print_sparse (E, "E = A(p,p)", cm)) ;
    cholmod_sort (E, cm) ;

    if (cm->status != CHOLMOD_OK)
    {
	/* out of memory */
	PR0 ((Common->file, "out of memory to sort E\n")) ;
	cholmod_free_sparse (&E, cm) ;
	return (NULL) ;
    }

    DEBUG (cholmod_print_sparse (E, "E = A(p,p), sorted", cm)) ;
    for (k = 0 ; k <= n ; k++)
    {
	Ep2 [k] = Ep [k] ;
    }

    return (E) ;
}


/* ========================================================================== */
/* === paraklete_factorize ================================================== */
/* ========================================================================== */

paraklete_numeric *paraklete_factorize
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
    int *Cperm, *Cn, *Cnz, *Ep, *Ei, *Ep2, *Map, *Cparent,
	*Cstart, *P, *Q, *Cp, *Ci, *Pc, *Qc, *Pinv, *Qinv, *Sched ;
    int cj, i, n, ncomponents, k, p, c, a, cn, k1, k2, kk, cnz,
	nfound, myid, npiv, ok, all_ok ;

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
    n = LUsymbolic->n ;
    ncomponents = LUsymbolic->ncomponents ;
    Cperm   = LUsymbolic->Cperm ;
    Cn      = LUsymbolic->Cn ;
    Cnz     = LUsymbolic->Cnz ;
    Cparent = LUsymbolic->Cparent ;
    Cstart  = LUsymbolic->Cstart  ;
    Sched   = LUsymbolic->Sched ;

    PR0 ((Common->file, "in factor: proc %d my_tries %d\n", myid, my_tries)) ;

#ifndef NDEBUG
    for (c = 0 ; c < ncomponents ; c++)
    {
	PR1 ((Common->file, "proc: %d node %d: %d to %d, parent %d Sched %d\n",
		Common->myid, c, Cstart [c], Cstart [c+1]-1, Cparent [c],
		Sched [c])) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* allocate and initialize the global LU factors */
    /* ---------------------------------------------------------------------- */

    LU = paraklete_allocate_numeric (LUsymbolic, Common) ;
    ok = (LU != NULL) ;
    PR1 ((Common->file, "factor proc %d alloc ok: %d\n", myid, ok)) ;

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
	    C = cholmod_allocate_sparse (Cn [c], Cn [c], Cnz [c], TRUE, TRUE,
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
	PR0 ((Common->file, "proc %d all fail factorize\n", Common->myid)) ;
	Common->status = PK_OUT_OF_MEMORY ;
	paraklete_free_numeric (&LU, Common) ;
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

	Map = cm->Iwork ;

	for (c = 0 ; c < ncomponents ; c++)
	{
	    PR1 ((Common->file, "distribute node c = %d\n", c)) ;
	    cn = Cn [c] ;

	    /* find mapping of global rows/columns to node c's local indices */
	    cj = 0 ;
	    for (a = c ; a != EMPTY ; a = Cparent [a])
	    {
		PR2 ((Common->file, "ancestor %d, for c %d, ncomp %d\n",
			a, c, ncomponents)) ;
		ASSERT (a >= c && a < ncomponents) ;
		k1 = Cstart [a] ;
		k2 = Cstart [a+1] ;
		PR2 ((Common->file, "k1 %d k2 %d\n", k1, k2)) ;
		for (k = k1 ; k < k2 ; k++)
		{
		    /* global index k becomes local index j */
		    PR3 ((Common->file, "  global: %d local %d\n", k, cj)) ;
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
	    for (k = k1 ; k < k2 ; k++)
	    {
		Cp [cj++] = cnz ;
		for (p = Ep2 [k] ; p < Ep [k+1] ; p++)
		{
		    ASSERT (Ei [p] >= k1) ;
		    Ci [cnz] = Map [Ei [p]] ;
		    Cx [cnz] = Ex [p] ;
		    cnz++ ;
		    ASSERT (cnz <= Cnz [c]) ;
		}
	    }

	    /* create columns for the ancestors of C.  These will become columns
	     * of U2 and S for node c. */
	    for (a = Cparent [c] ; a != EMPTY ; a = Cparent [a])
	    {
		PR2 ((Common->file, "ancestor %d\n", a)) ;
		PR2 ((Common->file, "k1 %d k2 %d\n", Cstart [a], Cstart [a+1])) ;
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
	    PR2 ((Common->file, "node: %d\n", c)) ;
	    DEBUG (cholmod_print_sparse (C, "send C to a node", cm)) ;

	    if (Sched [c] != myid)
	    {
		/* send this matrix to the process that owns node c */
		MPI (MPI_Isend (Cp, cn+1, MPI_INT, Sched [c],
			    TAG0, MPI_COMM_WORLD, &req)) ;
		MPI (MPI_Isend (Ci, cnz,  MPI_INT, Sched [c],
			    TAG0, MPI_COMM_WORLD, &req)) ;
		MPI (MPI_Isend (Cx, cnz,  MPI_DOUBLE, Sched [c],
			    TAG0, MPI_COMM_WORLD, &req)) ;
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
		MPI (MPI_Recv (LU->LUnode [c]->A->p, Cn [c] +1, MPI_INT, 0,
			    TAG0, MPI_COMM_WORLD, &ms)) ;
		MPI (MPI_Recv (LU->LUnode [c]->A->i, Cnz [c], MPI_INT, 0,
			    TAG0, MPI_COMM_WORLD, &ms)) ;
		MPI (MPI_Recv (LU->LUnode [c]->A->x, Cnz [c], MPI_DOUBLE, 0,
			    TAG0, MPI_COMM_WORLD, &ms)) ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* free temporary copy of A(p,p) */
    /* ---------------------------------------------------------------------- */

    /* only the root needs to do this */
    if (myid == 0)
    {
	LU->Ep2 = cholmod_free (n+1, sizeof (int), LU->Ep2, cm) ;
	cholmod_free_sparse (&(LU->E), cm) ;
    }


#if 0
    /* TODO: process 0 should free LUnode [c]->A, but only when recvd by c */
    MPI (MPI_Barrier (MPI_COMM_WORLD)) ;
    PR1 ((Common->file, "proc %d everybody OK so far2\n", Common->myid)) ;
    cm->print = 5 ;
    if (myid == 0)
    {
	for (c = 0 ; c < ncomponents ; c++)
	{
	    if (Sched [c] != 0)
	    {
		/* root process no longer needs submatrices sent to others */
		cholmod_free_sparse (&(LU->LUnode [c]->A), cm) ;
	    }
	}
    }
#endif

#ifndef NDEBUG
    for (c = 0 ; c < ncomponents ; c++)
    {
	if (Sched [c] == myid)
	{
	    PR1 ((Common->file, "proc %d Node %d original matrix:\n", myid, c)) ;
	    DEBUG (cholmod_print_sparse (LU->LUnode [c]->A, "A", cm)) ;
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
	    PR1 ((Common->file, "proc %d doing node %d\n", myid, c)) ;
	    if (!paraklete_factorize_node (c, LU, LUsymbolic, Common))
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
	MPI (MPI_Bcast (LU->LUnode [c]->header, PK_HEADER, MPI_INT,
		    Sched [c], MPI_COMM_WORLD)) ;
	kk += LU->LUnode [c]->PK_NFOUND ;
	Common->status = MIN (Common->status, LU->LUnode [c]->PK_STATUS) ;
    }

    if (kk < n)
    {
	PR0 ((Common->file, "proc %d Singular: %d of %d\n", myid, kk, n)) ;
	Common->status = PK_SINGULAR ;
    }

    if (Common->status != PK_OK)
    {
	/* out of memory, singular matrix, or unknown error */
	paraklete_free_numeric (&LU, Common) ;
	return (NULL) ;
    }

    /* all processes allocate space for the global permutation */
    for (c = 0 ; c < ncomponents ; c++)
    {
	npiv = LU->LUnode [c]->PK_NPIV ;
	if (LU->LUnode [c]->Pglobal == NULL)
	{
	    LU->LUnode [c]->Pglobal = cholmod_malloc (npiv, sizeof (int), cm) ;
	    LU->LUnode [c]->Qglobal = cholmod_malloc (npiv, sizeof (int), cm) ;
	}
    }

    ok = ok && (kk == n) && (cm->status == CHOLMOD_OK) ;

    /* check to see if all processes had space for P and Q */
    PR1 ((Common->file, "proc %d here in PQ: ok %d\n", Common->myid, ok)) ;
    all_ok = ok ;
    MPI (MPI_Allreduce (&ok, &all_ok, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD)) ;
    if (!all_ok)
    {
	/* out of memory */
	PR0 ((Common->file, "proc %d everybody fails in PQ\n", Common->myid)) ;
	Common->status = PK_OUT_OF_MEMORY ;
	paraklete_free_numeric (&LU, Common) ;
	return (NULL) ;
    }

    kk = 0 ;
    for (c = 0 ; c < ncomponents ; c++)
    {
	npiv = LU->LUnode [c]->PK_NPIV ;

	/* give Pglobal and Qglobal to all processors */
	MPI (MPI_Bcast (LU->LUnode [c]->Pglobal, npiv, MPI_INT, Sched [c],
		    MPI_COMM_WORLD)) ;
	MPI (MPI_Bcast (LU->LUnode [c]->Qglobal, npiv, MPI_INT, Sched [c],
		    MPI_COMM_WORLD)) ;

	Pc = LU->LUnode [c]->Pglobal ;
	Qc = LU->LUnode [c]->Qglobal ;

	nfound = LU->LUnode [c]->PK_NFOUND ;
	for (k = 0 ; k < nfound ; k++)
	{
	    /* row Pc[k] and column Qc[k] of E are the kth pivot row/column of
	     * node c.  They are the kk-th pivot row/column of the global LU
	     * factors. */
	    P [kk] = Cperm [Pc [k]] ;
	    Q [kk] = Cperm [Qc [k]] ;
	    kk++ ;
	}
    }
    ASSERT (kk == n) ;

    Pinv = LU->Pinv ;
    Qinv = LU->Qinv ;
    for (k = 0 ; k < n ; k++)
    {
	Pinv [P [k]] = k ;
	Qinv [Q [k]] = k ;
    }

    return (LU) ;
}
