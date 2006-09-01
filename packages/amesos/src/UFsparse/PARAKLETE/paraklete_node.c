/* ========================================================================== */
/* === paraklete_factorize_node ============================================= */
/* ========================================================================== */

#include "paraklete.h"

/* Factorize node c of the separator tree:
 *  (1) obtain the Schur complements from the children of node c
 *  (2) sum up the Schur complements, and the input matrix at this node c
 *  (3) factorize node c
 *  (4) send the Schur complement of node c to its parent
 *
 * PARAKLETE version 0.1: parallel sparse LU factorization.  May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* ========================================================================== */
/* === paraklete_free_children ============================================== */
/* ========================================================================== */

/* Free the Schur complement of each child of node c.  Free the input matrix
 * A for this node c. */

static void paraklete_free_children
(
    int c,
    paraklete_numeric *LU,
    paraklete_symbolic *LUsymbolic,
    paraklete_common *Common
)
{
    int *Child, *Childp ;
    int nchild, child, sn, cp ;
    cholmod_common *cm ;

    cm = &(Common->cm) ;

    /* free A for this node */
    cholmod_free_sparse (&(LU->LUnode [c]->A), cm) ;

    Child = LUsymbolic->Child ;
    Childp = LUsymbolic->Childp ;
    nchild = LU->LUnode [c]->nchild ;
    ASSERT (nchild == Childp [c+1] - Childp [c]) ;

    for (cp = 0 ; cp < nchild ; cp++)
    {
	/* free the Schur complement of the child */
	child = Child [Childp [c] + cp] ;
	sn = LU->LUnode [child]->PK_SN ;

	LU->LUnode [child]->sx = cholmod_free (
		LU->LUnode [child]->PK_SSIZE, sizeof (double),
		LU->LUnode [child]->sx, cm) ;
	LU->LUnode [child]->sp = cholmod_free (
		sn, sizeof (int), LU->LUnode [child]->sp, cm) ;
	LU->LUnode [child]->slen = cholmod_free (
		sn, sizeof (int), LU->LUnode [child]->slen, cm) ;
    }
}


/* ========================================================================== */
/* === paraklete_send_to_parent ============================================= */
/* ========================================================================== */

/* Send the Schur complement or an error message to the parent */

static void paraklete_send_to_parent
(
    int c,
    int status,
    int parent_id,
    paraklete_numeric *LU,
    paraklete_symbolic *LUsymbolic,
    paraklete_common *Common
)
{
    int ok ;
    cholmod_common *cm ;
    paraklete_node *LUnode ;
    MPI (MPI_Status ms) ;

    cm = &(Common->cm) ;
    LUnode = LU->LUnode [c] ;

    /* workspace W2 and C no longer needed */
    LUnode->W2 = cholmod_free (2*(LUnode->nlost_in), sizeof (int),
	    LUnode->W2, cm) ;
    cholmod_free_sparse (&(LUnode->C), cm) ;

    LUnode->PK_STATUS = status ;

    ok = TRUE ;
    if (parent_id != EMPTY && parent_id != Common->myid)
    {
	/* send header to process that owns the parent node (parent_id) */
	PR1 ((Common->file, "proc %d sends header parent proc %d\n",
		Common->myid, parent_id)) ;
	MPI (MPI_Send (LUnode->header, PK_HEADER, MPI_INT,
		    parent_id, TAG0, Common->mpicomm)) ;

	/* wait to see if parent_id is OK */
	PR1 ((Common->file, "proc %d wait for parent %d\n",
		Common->myid, parent_id)) ;
	MPI (MPI_Recv (&ok, 1, MPI_INT, parent_id, TAG0, Common->mpicomm, &ms)) ;

	/* if status is not PK_OK, then parent_id will send back ok = FALSE */
	ASSERT (IMPLIES (status != PK_OK, !ok)) ;

	if (ok)
	{
	    /* both parent_id and myid agree to send the Schur complement */
	    /* TODO: send this as one or two messages */
	    MPI (MPI_Rsend (LUnode->sp, LUnode->PK_SN, MPI_INT,
			parent_id, TAG0, Common->mpicomm)) ;
	    MPI (MPI_Rsend (LUnode->slen, LUnode->PK_SN, MPI_INT,
			parent_id, TAG0, Common->mpicomm)) ;
	    MPI (MPI_Rsend (LUnode->sx, LUnode->PK_SSIZE, MPI_DOUBLE,
			parent_id, TAG0, Common->mpicomm)) ;
	    MPI (MPI_Rsend (LUnode->Pglobal, LUnode->PK_NPIV, MPI_INT,
			parent_id, TAG0, Common->mpicomm)) ;
	    MPI (MPI_Rsend (LUnode->Qglobal, LUnode->PK_NPIV, MPI_INT,
			parent_id, TAG0, Common->mpicomm)) ;
	}
    }

    if (!ok || parent_id == EMPTY || parent_id != Common->myid)
    {
	/* free the Schur complement of node c if a failure occured, if node
	 * c has no parent, or if the Schur complement was sent to a
	 * different process (parent_id). */
	LUnode->sx = cholmod_free (LUnode->PK_SSIZE,
		sizeof (double), LUnode->sx, cm) ;
	LUnode->sp = cholmod_free (LUnode->PK_SN,
		sizeof (int), LUnode->sp, cm) ;
	LUnode->slen = cholmod_free (LUnode->PK_SN,
		sizeof (int), LUnode->slen, cm) ;
    }

    if (!ok)
    {
	/* free the Schur complements of the children if a failure occured */
	paraklete_free_children (c, LU, LUsymbolic, Common) ;
    }
}


/* ========================================================================== */
/* === paraklete_factorize_node ============================================= */
/* ========================================================================== */

int paraklete_factorize_node
(
    int c,
    paraklete_numeric *LU,
    paraklete_symbolic *LUsymbolic,
    paraklete_common *Common
)
{
    paraklete_node *LUnode ;
    cholmod_common *cm ;
    cholmod_sparse *C ;
    double *Cx, *W, *Sx, *Ax, *S ;
    int *Child, *Childp, *Lost, *Lostp, *Plost_in, *Qlost_in, *Cp, *Ci, *Flag,
	*Sp, *Si, *Slen, *Ap, *Ai, *Cn, *Plocal, *Qlocal, *Pglobal, *Qglobal,
	*Sched, *Pinv ;
    int cp, cnz, clnz, nchild, k1, k2, child, cn, npiv, k, j, p, i, nfound,
	mark, cj, ci, nlost_in, len, sn, myid, ssize, parent, ok, nmessages,
	status, parent_id ;

    MPI (MPI_Status ms) ;
    MPI (MPI_Request req [5]) ;

    /* ---------------------------------------------------------------------- */
    /* get local workspace */
    /* ---------------------------------------------------------------------- */

    cm = &(Common->cm) ;
    PR0 ((Common->file, "\n\n######################## FACTORIZE NODE %d\n", c)) ;

    /* ---------------------------------------------------------------------- */
    /* get the symbolic analysis of this node */
    /* ---------------------------------------------------------------------- */

    cnz    = LUsymbolic->Cnz [c] ;	/* # entries in A for this node c */
    Childp = LUsymbolic->Childp ;	/* Child [Childp [c]..Childp[c+1]-1] */
    Child  = LUsymbolic->Child ;	/* is list of children of node c */
    clnz   = LUsymbolic->Clnz [c] ;	/* est. # entries in L for node c */
    nchild = Childp [c+1] - Childp [c] ;    /* # of children of node c */
    k1     = LUsymbolic->Cstart [c] ;	/* global index of first pivot cand. */
    k2     = LUsymbolic->Cstart [c+1] ;	/* global index of last piv plus one */
    Cn     = LUsymbolic->Cn ;		/* dimension of each node */
    Sched  = LUsymbolic->Sched ;
    parent = LUsymbolic->Cparent [c] ;
    myid = Common->myid ;
    parent_id = (parent == EMPTY) ? EMPTY : Sched [parent] ;

    /* ---------------------------------------------------------------------- */
    /* get the arrowhead of the input matrix to factorize at this node */
    /* ---------------------------------------------------------------------- */

    LUnode = LU->LUnode [c] ;
    ASSERT (nchild == LUnode->nchild) ;
    Ap = LUnode->A->p ;	    /* A matrix of dimension Cn [c] */
    Ai = LUnode->A->i ;
    Ax = LUnode->A->x ;
    DEBUG (cholmod_print_sparse (LUnode->A, "Arrowhead", cm)) ;

#ifndef NDEBUG
    PR1 ((Common->file, "proc %d node %d nchild %d\n", myid, c, nchild)) ;
    {
	int cc ;
	for (cc = 0 ; cc < LUsymbolic->ncomponents ; cc++)
	{
	    PR2 ((Common->file,
		"proc %d: node %d sched %d cparent %d childp [%d:%d]\n",
		myid, cc, Sched [cc], LUsymbolic->Cparent [cc], Childp [cc],
		Childp [cc+1]-1)) ;
	    for (cp = 0 ; cp < Childp [cc+1] - Childp [cc] ; cp++)
	    {
		PR2 ((Common->file, "node %d: child node %d\n",
			    cc, Child [Childp[c] + cp])) ;
	    }
	}
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* post a non-blocking receive for the header information from each child */
    /* ---------------------------------------------------------------------- */

    PR1 ((Common->file, "proc %d posting recv at node %d, nchild %d status %d\n",
	    myid, c, nchild, cm->status)) ;
    nmessages = 0 ;
    for (cp = 0 ; cp < nchild ; cp++)
    {
	child = Child [Childp [c] + cp] ;
	PR2 ((Common->file, "proc %d child %d owned by %d\n",
		    myid, child, Sched [child])) ;
	if (Sched [child] != myid)
	{
	    PR2 ((Common->file, "parent proc %d awaits header from %d\n",
		    myid, Sched [child])) ;
	    MPI (MPI_Irecv (LU->LUnode [child]->header, PK_HEADER, MPI_INT,
		    Sched [child], TAG0, Common->mpicomm, &(LUnode->Req [cp]))) ;
	    nmessages++ ;
	}
	else
	{
	    MPI (LUnode->Req [cp] = MPI_REQUEST_NULL) ;
	}
    }

#ifdef NMPI
    ASSERT (nmessages == 0) ;
#endif

    /* ---------------------------------------------------------------------- */
    /* get header and Schur complement from each child */
    /* ---------------------------------------------------------------------- */

    ok = TRUE ;
    for (k = 0 ; k < nmessages ; k++)
    {

	/* ------------------------------------------------------------------ */
	/* get header from a child who is ready to send its Schur complement*/
	/* ------------------------------------------------------------------ */

	MPI (MPI_Waitany (nchild, LUnode->Req, &cp, &ms)) ;

	child = Child [Childp [c] + cp] ;
	ASSERT (Sched [child] != myid) ;

	status = LU->LUnode [child]->PK_STATUS ;
	sn = LU->LUnode [child]->PK_SN ;
	nfound = LU->LUnode [child]->PK_NFOUND ;
	npiv = LU->LUnode [child]->PK_NPIV ;
	ssize = LU->LUnode [child]->PK_SSIZE ;
	cn = LU->LUnode [child]->PK_NN ;

	ok = ok && (status == PK_OK) ;

	/* ------------------------------------------------------------------ */
	/* allocate space for Schur complement of the child */
	/* ------------------------------------------------------------------ */

	if (ok)
	{
	    /* allocate space for next message: the Schur complement */
	    LU->LUnode [child]->sp = cholmod_malloc (sn,
		    sizeof (int), cm) ;
	    LU->LUnode [child]->slen = cholmod_malloc (sn,
		    sizeof (int), cm) ;
	    LU->LUnode [child]->Pglobal = cholmod_malloc (npiv,
		    sizeof (int), cm) ;
	    LU->LUnode [child]->Qglobal = cholmod_malloc (npiv,
		    sizeof (int), cm) ;
	    LU->LUnode [child]->sx = cholmod_malloc (ssize,
		    sizeof (double), cm) ;

	    /* allocate space for foward/backsolve */
	    LU->LUnode [child]->X = cholmod_malloc (cn,
		    sizeof (double), cm) ;
	}

	/* check if we ran out of memory */
	ok = ok && (cm->status == CHOLMOD_OK) ;

	/* ------------------------------------------------------------------ */
	/* post receives for the Schur complement from the child */
	/* ------------------------------------------------------------------ */

	if (ok)
	{
	    /* both parent and child agree to send the Schur complement */
	    /* TODO: fold this information into fewer messages */
	    MPI (MPI_Irecv (LU->LUnode [child]->sp, sn, MPI_INT,
			Sched [child], TAG0, Common->mpicomm, &(req [0]))) ;
	    MPI (MPI_Irecv (LU->LUnode [child]->slen, sn, MPI_INT,
			Sched [child], TAG0, Common->mpicomm, &(req [1]))) ;
	    MPI (MPI_Irecv (LU->LUnode [child]->sx, ssize, MPI_DOUBLE,
			Sched [child], TAG0, Common->mpicomm, &(req [2]))) ;
	    MPI (MPI_Irecv (LU->LUnode [child]->Pglobal, npiv, MPI_INT,
			Sched [child], TAG0, Common->mpicomm, &(req [3]))) ;
	    MPI (MPI_Irecv (LU->LUnode [child]->Qglobal, npiv, MPI_INT,
			Sched [child], TAG0, Common->mpicomm, &(req [4]))) ;
	}

	/* ------------------------------------------------------------------ */
	/* tell child that we are ready to receive its Schur complement */
	/* ------------------------------------------------------------------ */

	PR1 ((Common->file, "parent proc %d replies to proc %d\n",
		    myid, Sched [child])) ;
	MPI (MPI_Send (&ok, 1, MPI_INT, Sched [child], TAG0, Common->mpicomm)) ;

	/* ------------------------------------------------------------------ */
	/* wait for the Schur complement to be received */
	/* ------------------------------------------------------------------ */

	if (ok)
	{
	    MPI (MPI_Waitall (5, req, MPI_STATUSES_IGNORE)) ;
	}
    }

    PR1 ((Common->file, "proc %d node %d arrowhead again\n", myid, c)) ;
    DEBUG (cholmod_print_sparse (LUnode->A, "Arrowhead again", cm)) ;

    /* ---------------------------------------------------------------------- */
    /* report failure to parent of c, if a failure occured */
    /* ---------------------------------------------------------------------- */

    if (!ok)
    {
	PR0 ((Common->file, "proc %d, node %d, report failure to parent: %d\n",
		myid, c, parent_id)) ;
	paraklete_send_to_parent (c, PK_OUT_OF_MEMORY, parent_id,
		LU, LUsymbolic, Common) ;
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get lost pivots and Schur complement from each child */
    /* ---------------------------------------------------------------------- */

    nlost_in = 0 ;

    Lost  = LUnode->Lost ;
    Lostp = LUnode->Lostp ;
    LUnode->nchild = nchild ;

    for (cp = 0 ; cp < nchild ; cp++)
    {
	/* find the failed pivot rows/cols of this child and add them
	 * to the pivot candidate sets of this node s */
	child = Child [Childp [c] + cp] ;
	Lostp [cp] = nlost_in ;
	Lost [cp] = LU->LUnode [child]->PK_NLOST ;
	PR1 ((Common->file, "child %d lost %d \n", child, Lost [cp])) ;
	nlost_in += Lost [cp] ;
	cnz += LU->LUnode [child]->PK_SNZ ;
    }
    Lostp [nchild] = nlost_in ;
    LUnode->nlost_in = nlost_in ;

    /* The matrix to factorize at this node is cn-by-cn.  Up to npiv pivots
     * can be found in this node */
    cn = nlost_in + Cn [c] ;
    npiv = nlost_in + (k2 - k1) ;
    LUnode->PK_NN = cn ;

    LUnode->W2 = cholmod_malloc (2*nlost_in, sizeof (int), cm) ;
    Plost_in = LUnode->W2 ;		/* size nlost_in */
    Qlost_in = LUnode->W2 + nlost_in ;	/* size nlost_in */

    /* get workspace */
    cholmod_allocate_work (cn, 3*cn, cn, cm) ;

    /* ---------------------------------------------------------------------- */
    /* C = original entries in this node, plus Schur complements of children */
    /* ---------------------------------------------------------------------- */

    /* TODO: assemble the Schur complements but do not compute C, just to
     * determine cnz.  Then do the cholmod_allocate_sparse, below: */

    C = cholmod_allocate_sparse (cn, cn, cnz, FALSE, TRUE, 0, TRUE, cm) ;
    LUnode->C = C ;

    if (cm->status != CHOLMOD_OK)
    {
	/* out of memory; tell the parent that we failed */
	PR0 ((Common->file, "node %d proc %d fails here\n", c, myid)) ;
	paraklete_send_to_parent (c, PK_OUT_OF_MEMORY, parent_id,
		LU, LUsymbolic, Common) ;
	return (FALSE) ;
    }

    Cp = C->p ;
    Ci = C->i ;
    Cx = C->x ;

    /* ---------------------------------------------------------------------- */
    /* concatenate and remap the lost pivot columns of each child */
    /* ---------------------------------------------------------------------- */

    k = 0 ;
    cnz = 0 ;
    for (cp = 0 ; cp < nchild ; cp++)
    {
	/* get the Schur complement of the child */
	child = Child [Childp [c] + cp] ;
	sn = LU->LUnode [child]->PK_SN ;

	PR1 ((Common->file, "\nconcatenate lost child %d, Lost [cp] %d\n",
		    child, Lost [cp])) ;
	S    = LU->LUnode [child]->sx ;
	Sp   = LU->LUnode [child]->sp ;
	Slen = LU->LUnode [child]->slen ;
	nfound = LU->LUnode [child]->PK_NFOUND ;

	PR1 ((Common->file, "Cn [c] is %d\n", Cn [c])) ;
	PR1 ((Common->file, "child is %d\n", child)) ;
	PR1 ((Common->file, "Lost[cp] is %d\n", Lost [cp])) ;
	PR1 ((Common->file, "sn %d  Lost[cp]+Cn[cn] %d\n",
	    sn, Lost [cp] + Cn [c]));

	ASSERT (sn == Lost [cp] + Cn [c]) ;

	for (j = 0 ; j < Lost [cp] ; j++)
	{
	    /* column j of the Schur complement of the child becomes column
	     * k of C */
	    PR2 ((Common->file, "lost child col %d becomes col %d of C\n",
		j, k)) ;
	    Cp [k] = cnz ;
	    GET_COLUMN (Sp, Slen, S, j, Si, Sx, len) ;
	    for (p = 0 ; p < len ; p++)
	    {
		i = Si [p] ;
		ci = i + ((i < Lost [cp]) ? Lostp [cp] : (nlost_in - Lost[cp]));
		Ci [cnz] = ci ;
		Cx [cnz] = Sx [p] ;
		PR3 ((Common->file,
			"  Lost child col: row %d newrow %d value %g\n",
			i, ci, Sx [p])) ;
		cnz++ ;
	    }
	    /* get the lost pivot row/column from the child */
	    Plost_in [k] = LU->LUnode [child]->Pglobal [nfound+j] ;
	    Qlost_in [k] = LU->LUnode [child]->Qglobal [nfound+j] ;
	    k++ ;
	}
    }
    ASSERT (k == nlost_in) ;

    /* ---------------------------------------------------------------------- */
    /* assemble original entries and Schur complements of each child */
    /* ---------------------------------------------------------------------- */

    Flag = cm->Flag ;
    W = cm->Xwork ;
    mark = cholmod_clear_flag (cm) ;
    ASSERT (cm->xworksize >= cn * sizeof (double)) ;
    DEBUG (for (cj = 0 ; cj < cn ; cj++) ASSERT (W [cj] == 0)) ;

    for (j = 0 ; j < Cn [c] ; j++)
    {

	/* Column j is nominal column index, if there were no failed pivots.
	 * The column is shifted over to accomodate incoming lost pivot columns
	 * from the children, to obtain the new column cj. */

	cj = j + nlost_in ;
	Cp [cj] = cnz ;
	PR2 ((Common->file, "\nAdd column %d of Schur\n", cj)) ;

	/* ------------------------------------------------------------------ */
	/* add up all the contributions to column cj of C */
	/* ------------------------------------------------------------------ */

	/* scatter the original entries of column j.  A is the input arrowhead
	 * for this node, with row/col indices already mapped to the local
	 * row/col index space (in range 0 to Cn [c]). */
	for (p = Ap [j] ; p < Ap [j+1] ; p++)
	{
	    i = Ai [p] ;
	    ci = i + nlost_in ;
	    PR3 ((Common->file, "scatter original (%d,%d) %g to (%d,%d)\n",
			i, j, Ax [p], ci, cj)) ;
	    Flag [ci] = mark ;
	    Ci [cnz++] = ci ;
	    W [ci] = Ax [p] ;
	}

	/* scatter and add the contributions from each child */
	for (cp = 0 ; cp < nchild ; cp++)
	{
	    /* get column j+Lost[cp] of the Schur complement of the child */
	    child = Child [Childp [c] + cp] ;
	    S    = LU->LUnode [child]->sx ;
	    Sp   = LU->LUnode [child]->sp ;
	    Slen = LU->LUnode [child]->slen ;
	    GET_COLUMN (Sp, Slen, S, j + Lost [cp], Si, Sx, len) ;
	    for (p = 0 ; p < len ; p++)
	    {
		i = Si [p] ;
		ci = i + ((i < Lost [cp]) ? Lostp [cp] : (nlost_in - Lost[cp]));
		if (Flag [ci] < mark)
		{
		    Flag [ci] = mark ;
		    Ci [cnz++] = ci ;
		}
		W [ci] += Sx [p] ;
	    }
	}

	/* gather the results */
	PR2 ((Common->file, "gather %d to %d-1\n", Cp [cj], cnz)) ;
	for (p = Cp [cj] ; p < cnz ; p++)
	{
	    ci = Ci [p] ;
	    Cx [p] = W [ci] ;
	    W [ci] = 0 ;
	    PR3 ((Common->file, "%d: %g\n", ci, Cx [p])) ;
	}
	DEBUG (for (cj = 0 ; cj < cn ; cj++) ASSERT (W [cj] == 0)) ;

	/* clear Flag array */
	mark = cholmod_clear_flag (cm) ;
    }

    Cp [cn] = cnz ;

    /* shrink C to be just large enough */
    cholmod_reallocate_sparse (cnz, C, cm) ;
    ASSERT (cm->status == CHOLMOD_OK) ;

    /* ---------------------------------------------------------------------- */
    /* A for this node, and Schur complements of children, no longer needed */
    /* ---------------------------------------------------------------------- */

    paraklete_free_children (c, LU, LUsymbolic, Common) ;

    /* ---------------------------------------------------------------------- */
    /* allocate the LU factors for this node */
    /* ---------------------------------------------------------------------- */

    LUnode->lusize = 1.2 * (2 * ((3*clnz)/2 + 2)) ;
    LUnode->PK_NPIV = npiv ;

    LUnode->llen = cholmod_malloc (npiv, sizeof (int), cm) ;
    LUnode->lp   = cholmod_malloc (npiv, sizeof (int), cm) ;

    LUnode->ulen = cholmod_malloc (cn, sizeof (int), cm) ;
    LUnode->up   = cholmod_malloc (cn, sizeof (int), cm) ;

    LUnode->Plocal  = cholmod_malloc (npiv, sizeof (int), cm) ;
    LUnode->Pglobal = cholmod_malloc (npiv, sizeof (int), cm) ;
    LUnode->Qlocal  = cholmod_malloc (npiv, sizeof (int), cm) ;
    LUnode->Qglobal = cholmod_malloc (npiv, sizeof (int), cm) ;

    LUnode->Pinv = cholmod_malloc (npiv, sizeof (int), cm) ;
    LUnode->Qinv = cholmod_malloc (npiv, sizeof (int), cm) ;

    /* this block of memory may grow in size */
    LUnode->ix = cholmod_malloc (LUnode->lusize, sizeof (double), cm) ;
    PR1 ((Common->file, "allocated LUnode->ix: size %d\n", LUnode->lusize)) ;

    if (cm->status != CHOLMOD_OK)
    {
	/* out of memory */
	paraklete_send_to_parent (c, PK_OUT_OF_MEMORY, parent_id,
		LU, LUsymbolic, Common) ;
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* factorize P*C*Q into L*U+S */
    /* ---------------------------------------------------------------------- */

    ASSERT (cholmod_print_sparse (C, "C = A + sum (Schur)", cm)) ;
    ok = paraklete_kernel (C, LUnode, Common) ;

    if (!ok)
    {
	/* out of memory */
	paraklete_send_to_parent (c, PK_OUT_OF_MEMORY, parent_id,
		LU, LUsymbolic, Common) ;
	return (FALSE) ;
    }

    PR0 ((Common->file, "============= NODE %d NPIV %d FOUND %d LOST %d\n", c,
	    LUnode->PK_NPIV, LUnode->PK_NFOUND,
	    LUnode->PK_NPIV - LUnode->PK_NFOUND)) ;

    /* TODO ensure the kernel does this */
    cm->mark = EMPTY ;
    cholmod_clear_flag (cm) ;

    Plocal = LUnode->Plocal ;
    Pinv = LUnode->Pinv ;
    for (k = 0 ; k < npiv ; k++)
    {
	i = Plocal [k] ;
	ASSERT (i >= 0 && i < npiv) ;
	ASSERT (Pinv [i] == k) ;
    }

    /* ---------------------------------------------------------------------- */
    /* determine the global pivot ordering (not including LUsymbolic->Cperm) */
    /* ---------------------------------------------------------------------- */

    Plocal = LUnode->Plocal ;
    Qlocal = LUnode->Qlocal ;

    Pglobal = LUnode->Pglobal ;
    Qglobal = LUnode->Qglobal ;

    DEBUG (cholmod_print_perm (Plocal, npiv, npiv, "Node P, local", cm)) ;
    DEBUG (cholmod_print_perm (Qlocal, npiv, npiv, "Node Q, local", cm)) ;

    nfound = LUnode->PK_NFOUND ;

    for (k = 0 ; k < npiv ; k++)
    {
	i = Plocal [k] ;
	if (i < nlost_in)
	{
	    /* row i was the ith incoming lost row, from the children */
	    Pglobal [k] = Plost_in [i] ;
	}
	else
	{
	    /* row i was an original candidate for this node.  It was shifted
	     * by nlost_in to obtain the local index.  Then k1 needs to be added
	     * to obtain the global index. */
	    Pglobal [k] = i - nlost_in + k1 ;
	}
    }
    for (k = 0 ; k < npiv ; k++)
    {
	j = Qlocal [k] ;
	if (j < nlost_in)
	{
	    /* col j was the jth incoming lost col, from the children */
	    Qglobal [k] = Qlost_in [j] ;
	}
	else
	{
	    /* col j was an original candidate for this node.  It was shifted
	     * by nlost_in to obtain the local index.  Then k1 needs to be added
	     * to obtain the global index. */
	    Qglobal [k] = j - nlost_in + k1 ;
	}
    }

    LUnode->PK_NLOST = npiv - nfound ;

    DEBUG (cholmod_print_perm (Pglobal, npiv, LU->n, "Node P, global", cm)) ;
    DEBUG (cholmod_print_perm (Qglobal, npiv, LU->n, "Node Q, global", cm)) ;

    /* ---------------------------------------------------------------------- */
    /* allocate space for subsequent forward/backsolve */
    /* ---------------------------------------------------------------------- */

    LU->LUnode [c]->X = cholmod_malloc (cn, sizeof (double), cm) ;
    ASSERT (k2-k1 == LUnode->nk) ;
    LU->LUnode [c]->B = cholmod_malloc (k2-k1, sizeof (double), cm) ;

    if (cm->status != CHOLMOD_OK)
    {
	/* out of memory */
	paraklete_send_to_parent (c, PK_OUT_OF_MEMORY, parent_id,
		LU, LUsymbolic, Common) ;
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* send Schur complement to the parent */
    /* ---------------------------------------------------------------------- */

    paraklete_send_to_parent (c, PK_OK, parent_id, LU, LUsymbolic, Common) ;

    return (TRUE) ;
}
