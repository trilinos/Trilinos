/* ========================================================================== */
/* === Paraklete/paraklete_btf.c ============================================ */
/* ========================================================================== */

/* User-callable functions that combine the functions of Paraklete,
 * KLU, and BTF.  See pb.c for examples on how to use these functions.
 *
 * PARAKLETE version 0.3: parallel sparse LU factorization.  Nov 13, 2007
 * Copyright (C) 2007, Univ. of Florida.  Author: Timothy A. Davis
 */

#include "amesos_paraklete_decl.h"

/* ========================================================================== */
/* === paraklete_btf_bcast_symbolic ========================================= */
/* ========================================================================== */

#ifndef NMPI
static Int paraklete_btf_bcast_symbolic
(
    paraklete_btf_symbolic **LU_btf_symbolic_handle,
    paraklete_common *Common
)
{
    paraklete_btf_symbolic *LU_btf_symbolic = NULL ;
    Int n, nblocks, header [4] ;
    int ok, all_ok ;
    cholmod_common *cm ;
    KLU_common *km ;

    cm = &(Common->cm) ;
    km = &(Common->km) ;

    n = EMPTY ;
    nblocks = EMPTY ;

    /* broadcast number of diagonal blocks in the BTF form, or -1 if failure */
    if (Common->myid == 0)
    {
	LU_btf_symbolic = *LU_btf_symbolic_handle ;
	if (LU_btf_symbolic != NULL)
	{
	    n = LU_btf_symbolic->n ;
	    nblocks = LU_btf_symbolic->nblocks ;
	}
    }
    else
    {
	/* other processes do not yet have the BTF symbolic object */
	*LU_btf_symbolic_handle = NULL ;
    }

    /* broadcast the size of the object, or -1 if a failure occurred */
    header [0] = n ;
    header [1] = nblocks ;
    MPI (MPI_Bcast (&header, 2, MPI_Int, TAG0, MPI_COMM_WORLD)) ;
    n = header [0] ;
    nblocks = header [1] ;
    if (n == EMPTY)
    {
	return (FALSE) ;
    }

    if (Common->myid != 0)
    {
	/* allocate the copy of this analysis on the slave */
	LU_btf_symbolic = amesos_paraklete_btf_alloc_symbolic (n, nblocks, Common) ;
	*LU_btf_symbolic_handle = LU_btf_symbolic ;
        /* TODO return if failed, and broadcast error to all processes */
        /* PARAKLETE_ERROR already called */
    }

    /* all processes find out if any one process fails to allocate memory */
    ok = (cm->status == CHOLMOD_OK) && (LU_btf_symbolic != NULL) ;
    all_ok = ok ;
    MPI (MPI_Allreduce (&ok, &all_ok, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD)) ;
    if (!all_ok)
    {
	/* out of memory; inform all processes */
	PR0 ((Common->file, "proc "ID" all fail in analyze\n", Common->myid)) ;
	amesos_paraklete_btf_free_symbolic (&LU_btf_symbolic, Common) ;
	*LU_btf_symbolic_handle = NULL ;
	return (FALSE) ;
    }

    /* broadcast the contents of the LU_btf_symbolic object (excluding
     * the LUsymbolic analysis of each diagonal block) */
    MPI (MPI_Bcast (LU_btf_symbolic->Mem_n, 3*n+1, MPI_Int, TAG0,
        MPI_COMM_WORLD)) ;
    return (TRUE) ;
}
#endif


/* ========================================================================== */
/* === paraklete_btf_alloc_symbolic ========================================= */
/* ========================================================================== */

paraklete_btf_symbolic *amesos_paraklete_btf_alloc_symbolic
(
    Int n,
    Int nblocks,
    paraklete_common *Common
)
{
    cholmod_common *cm ;
    KLU_common *km ;
    paraklete_btf_symbolic *LU_btf_symbolic ;
    Int *p ;
    int ok = TRUE ;
    size_t n3 ;

    cm = &(Common->cm) ;
    km = &(Common->km) ;

    /* n3 = 3*n+1 */
    n3 = CHOLMOD (mult_size_t) (n, 3, &ok) ;
    n3 = CHOLMOD (add_size_t) (n3, 1, &ok) ;
    if (!ok) PARAKLETE_ERROR (PK_TOO_LARGE, "problem too large") ;

    LU_btf_symbolic = CHOLMOD (malloc) (1, sizeof (paraklete_btf_symbolic), cm) ;
    if (cm->status != CHOLMOD_OK)
    {
        PARAKLETE_ERROR (PK_OUT_OF_MEMORY, "out of memory") ;
	return (NULL) ;
    }

    LU_btf_symbolic->n = n ;
    LU_btf_symbolic->nblocks = nblocks ;
    LU_btf_symbolic->cnz = 0 ;
    LU_btf_symbolic->fnz = 0 ;
    p = LU_btf_symbolic->Mem_n = CHOLMOD (calloc) (n3, sizeof (Int), cm) ;
    LU_btf_symbolic->Qbtf = p ;			    /* size n */
    LU_btf_symbolic->Pbinv = p + n ;		    /* size n */
    LU_btf_symbolic->Rbtf = p + 2*n ;		    /* size n+1 */
    /* only nblocks of the LUsymbolic array is used: */
    LU_btf_symbolic->LUsymbolic = CHOLMOD (calloc) (n, sizeof (void *), cm) ;

    if (cm->status != CHOLMOD_OK)
    {
        /* TODO free object and return NULL if malloc fails */
        PARAKLETE_ERROR (PK_OUT_OF_MEMORY, "out of memory") ;
    }
    return (LU_btf_symbolic) ;
}


/* ========================================================================== */
/* === paraklete_btf_analyze ================================================ */
/* ========================================================================== */

paraklete_btf_symbolic *amesos_paraklete_btf_analyze
(
    cholmod_sparse *A,		    /* matrix to analyze */ 
    paraklete_common *Common
)
{
    double work = 0 ;
    Int *Pbtf = NULL, *Qbtf = NULL, *Rbtf = NULL, *Work = NULL, *Ap = NULL,
	*Ai = NULL, *Cp = NULL, *Ci = NULL, *Pbinv = NULL ;
    cholmod_sparse *C = NULL ;
    cholmod_common *cm ;
    KLU_common *km ;
    paraklete_btf_symbolic *LU_btf_symbolic = NULL ;
    void **LUsymbolic ;
    Int nblocks, n = 0, nmatch, block, k, k1, k2, inew, jold, anz = 0, fnz,
	cnz, p, nk ;

    /* ---------------------------------------------------------------------- */
    /* root processor finds the BTF ordering */
    /* ---------------------------------------------------------------------- */

    cm = &(Common->cm) ;
    km = &(Common->km) ;

    if (Common->myid == 0)
    {
	Ap = A->p ;
	Ai = A->i ;
	n = A->nrow ;
	anz = Ap [n] ;

	LU_btf_symbolic = amesos_paraklete_btf_alloc_symbolic (n, n, Common) ;
        if (LU_btf_symbolic == NULL)
        {
            /* TODO return NULL if malloc fails */
            /* PARAKLETE_ERROR already called */ ;
        }

	Qbtf = LU_btf_symbolic->Qbtf ;		    /* size n */
	Rbtf = LU_btf_symbolic->Rbtf ;		    /* size n+1 */
	Pbinv = LU_btf_symbolic->Pbinv ;	    /* size n */
	LUsymbolic = LU_btf_symbolic->LUsymbolic ;  /* size n */

	/* ------------------------------------------------------------------ */
	/* find the BTF ordering */
	/* ------------------------------------------------------------------ */

	Work = CHOLMOD (malloc) (n, 6*sizeof (Int), cm) ;
	Pbtf = Work + 5*n ;
        if (cm->status != CHOLMOD_OK)
        {
            /* TODO return if malloc fails */
            PARAKLETE_ERROR (PK_OUT_OF_MEMORY, "out of memory") ;
        }

	nblocks = BTF_order (n, Ap, Ai, 0, &work, Pbtf, Qbtf, Rbtf, &nmatch,
	    Work) ;

	/* Pbinv = inverse of Pbtf */
	for (k = 0 ; k < n ; k++)
	{
	    Pbinv [Pbtf [k]] = k ;
	}

	CHOLMOD (free) (n, 6*sizeof (Int), Work, cm) ;
	LU_btf_symbolic->nblocks = nblocks ;
    }

    /* ---------------------------------------------------------------------- */
    /* broadcast the BTF information */
    /* ---------------------------------------------------------------------- */

    MPI (paraklete_btf_bcast_symbolic (&LU_btf_symbolic, Common)) ;
    /* TODO return if broadcast fails */
    Rbtf = LU_btf_symbolic->Rbtf ;
    LUsymbolic = LU_btf_symbolic->LUsymbolic ;
    nblocks = LU_btf_symbolic->nblocks ;

    /* ---------------------------------------------------------------------- */
    /* symbolic analysis of diagonal blocks of A(p,q) */
    /* ---------------------------------------------------------------------- */

    if (Common->myid == 0)
    {
	/* C = empty n-by-n sparse matrix */
	fnz = 0 ;
	C = CHOLMOD (allocate_sparse) (n, n, anz, FALSE, TRUE, 0, CHOLMOD_PATTERN,
	    cm) ;
        if (cm->status != CHOLMOD_OK)
        {
            /* TODO return if malloc fails */
            PARAKLETE_ERROR (PK_OUT_OF_MEMORY, "out of memory") ;
        }
	Cp = C->p ;
	Ci = C->i ;
	LU_btf_symbolic->cnz = 0 ;
    }

    for (block = 0 ; block < nblocks ; block++)
    {
	k1 = Rbtf [block] ;
	k2 = Rbtf [block+1] ;
	nk = k2 - k1 ;

	if (Common->myid == 0)
	{
	    cnz = 0 ;
	    for (k = k1 ; k < k2 ; k++)
	    {
		Cp [k - k1] = cnz ;		    /* note change of index */
		jold = Qbtf [k] ;
		for (p = Ap [jold] ; p < Ap [jold+1] ; p++)
		{
		    inew = Pbinv [Ai [p]] ;
		    if (inew < k1)
		    {
			/* this entry goes in the off-diagonal part */
			fnz++ ;
		    }
		    else
		    {
			/* this entry goes in the diagonal block */
			Ci [cnz++] = inew - k1 ;  /* note the change of index */
		    }
		}
	    }
	    LU_btf_symbolic->cnz = MAX (LU_btf_symbolic->cnz, cnz) ;
	    Cp [nk] = cnz ;
	    C->nrow = nk ;
	    C->ncol = nk ;
	}

	if (nk == 1)
	{
	    /* singleton; nothing to do */
	    LUsymbolic [block] = NULL ;
	}
        else if (nk < 1000)
        {
            if (Common->myid == 0)
            {
                /* use KLU on process 0 */
                LUsymbolic [block] = (void *) KLU_analyze (nk, Cp, Ci, km) ;
                if (km->status != KLU_OK)
                {
                    /* TODO return if KLU_analyze failed; broadcast the error */
                    PARAKLETE_ERROR (PK_UNKNOWN, "KLU analyze failed") ;
                }
            }
            else
            {
                /* small block; nothing to do for other processors*/
                LUsymbolic [block] = NULL ;
            }
        }
	else
	{
            /* parallel factorization of a large block */
	    /* process 0 does the work and broadcasts it to all the others */
	    LUsymbolic [block] = (void *) amesos_paraklete_analyze (C, Common) ;
            if (LUsymbolic [block] == NULL)
            {
	        /* TODO return if analyze fails and broadcast the error */
                /* note that PARAKLETE_ERROR has already been called */
            }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* free workspace sparse matrix C */
    /* ---------------------------------------------------------------------- */

    if (Common->myid == 0)
    {
        /* process 0 finalizes the object and frees its workspace */
	LU_btf_symbolic->fnz = fnz ;
	C->nrow = n ;
	C->ncol = n ;
	CHOLMOD (free_sparse) (&C, cm) ;
    }

    /* ---------------------------------------------------------------------- */
    /* return the result */
    /* ---------------------------------------------------------------------- */

    return (LU_btf_symbolic) ;
}


/* ========================================================================== */
/* === paraklete_btf_factorize ============================================== */
/* ========================================================================== */

paraklete_btf_numeric *amesos_paraklete_btf_factorize
(
    cholmod_sparse *A,				/* matrix to analyze */ 
    paraklete_btf_symbolic *LU_btf_symbolic,	/* symbolic analysis */
    paraklete_common *Common
)
{
    cholmod_sparse *F = NULL, *C = NULL ;
    cholmod_common *cm ;
    KLU_common *km ;
    double *Ax = NULL, *Cx = NULL, *Fx = NULL, *Singleton = NULL ;
    Int *Qbtf = NULL, *Rbtf = NULL, *Ap = NULL, *Ai = NULL, *Cp = NULL,
	*Ci = NULL, *Pbinv = NULL, *Fp = NULL, *Fi = NULL ;
    Int fnz = 0, cnz = 0, nblocks, n, k, k1, k2, block, jold, inew, p, nk ;
    int ok = TRUE, all_ok ;
    void **LUsymbolic = NULL;
    paraklete_btf_numeric *LU_btf_numeric = NULL ;
    void **LUnumeric = NULL ;

    /* ---------------------------------------------------------------------- */
    /* get inputs and allocate result */
    /* ---------------------------------------------------------------------- */

    cm = &(Common->cm) ;
    km = &(Common->km) ;

    LU_btf_numeric = CHOLMOD (malloc) (1, sizeof (paraklete_btf_numeric), cm) ;
    if (cm->status != CHOLMOD_OK)
    {
        PARAKLETE_ERROR (PK_OUT_OF_MEMORY, "out of memory") ;
	return (NULL) ;
    }

    if (Common->myid == 0)
    {
	fnz = LU_btf_symbolic->fnz ;
	cnz = LU_btf_symbolic->cnz ;
	n = A->nrow ;
	F = CHOLMOD (allocate_sparse) (n, n, fnz, FALSE, TRUE, 0,
	    CHOLMOD_REAL, cm) ;
	C = CHOLMOD (allocate_sparse) (n, n, cnz, FALSE, TRUE, 0,
	    CHOLMOD_REAL, cm) ;
        if (cm->status != CHOLMOD_OK)
        {
	    /* TODO free LU_btf_numeric and return NULL;
             * broadcast error to all processes */
            PARAKLETE_ERROR (PK_OUT_OF_MEMORY, "out of memory") ;
        }
    }

    n = LU_btf_symbolic->n ;
    nblocks = LU_btf_symbolic->nblocks ;
    Rbtf = LU_btf_symbolic->Rbtf ;
    Qbtf = LU_btf_symbolic->Qbtf ;
    Pbinv = LU_btf_symbolic->Pbinv ;
    LUsymbolic = LU_btf_symbolic->LUsymbolic ;

    if (Common->myid == 0)
    {
	fnz = 0 ;
	Fp = F->p ;
	Fi = F->i ;
	Fx = F->x ;
	Ap = A->p ;
	Ai = A->i ;
	Ax = A->x ;
	Cp = C->p ;
	Ci = C->i ;
	Cx = C->x ;
	Singleton = CHOLMOD (calloc) (nblocks, sizeof (double), cm) ;
        if (cm->status != CHOLMOD_OK)
        {
	    /* TODO free LU_btf_numeric and return NULL;
             * broadcast error to all processes */
            PARAKLETE_ERROR (PK_OUT_OF_MEMORY, "out of memory") ;
        }
    }

    /* all processes do this */
    LU_btf_numeric->Singleton = Singleton ;
    LU_btf_numeric->LUnumeric = LUnumeric =
	CHOLMOD (calloc) (nblocks, sizeof (void *), cm) ;
    LU_btf_numeric->F = F ;
    LU_btf_numeric->nblocks = nblocks ;

    if (cm->status != CHOLMOD_OK)
    {
        /* TODO free LU_btf_numeric and return NULL;
         * broadcast error to all processes */
        PARAKLETE_ERROR (PK_OUT_OF_MEMORY, "out of memory") ;
    }

    /* ---------------------------------------------------------------------- */
    /* factorize each diagonal block, and construct F */
    /* ---------------------------------------------------------------------- */

    for (block = 0 ; block < nblocks ; block++)
    {
	k1 = Rbtf [block] ;
	k2 = Rbtf [block+1] ;
	nk = k2 - k1 ;

	if (Common->myid == 0)
	{
	    cnz = 0 ;
	    for (k = k1 ; k < k2 ; k++)
	    {
		Fp [k] = fnz ;
		Cp [k - k1] = cnz ;		    /* note change of index */
		jold = Qbtf [k] ;
		for (p = Ap [jold] ; p < Ap [jold+1] ; p++)
		{
		    inew = Pbinv [Ai [p]] ;
		    if (inew < k1)
		    {
			/* this entry goes in the off-diagonal part */
			Fi [fnz] = inew ;
			Fx [fnz] = Ax [p] ;
			fnz++ ;
		    }
		    else
		    {
			/* this entry goes in the diagonal block */
			Ci [cnz] = inew - k1 ;  /* note the change of index */
			Cx [cnz] = Ax [p] ;
			cnz++ ;
		    }
		}
	    }
            Cp [nk] = cnz ;
            C->nrow = nk ;
            C->ncol = nk ;
	}

	if (nk == 1)
	{
	    /* singleton */
	    if (Common->myid == 0)
	    {
		Singleton [block] = Cx [0] ;
	    }
            LUnumeric [block] = NULL ;
	}
        else if (nk < 1000)
        {
            if (Common->myid == 0)
            {
                /* use KLU on process 0 */
                LUnumeric [block] = (void *) KLU_factor (Cp, Ci, Cx,
                    (KLU_symbolic *) LUsymbolic [block], km) ;
                if (km->status != KLU_OK)
                {
                    /* TODO return if KLU failed; broadcast the error */
                    PARAKLETE_ERROR (PK_UNKNOWN, "KLU factorize failed") ;
                    ok = FALSE ;
                }
            }
            else
            {
                /* small block; nothing to do for other processors */
                LUnumeric [block] = NULL ;
            }
        }
	else
	{
            /* parallel factorization of a large block */
	    LUnumeric [block] = (void *) amesos_paraklete_factorize (C,
                    (paraklete_symbolic *) LUsymbolic [block], Common) ;
            if (LUnumeric [block] == NULL)
            {
                /* PARAKLETE failed */
                ok = FALSE ;
            }
	}
    }

    if (Common->myid == 0)
    {
        /* process 0 frees its workspace */
	Fp [n] = fnz ;
	C->nrow = n ;
	C->ncol = n ;
	CHOLMOD (free_sparse) (&C, cm) ;
    }

    /* all processes find out if any one process fails */
    all_ok = ok ;
    MPI (MPI_Allreduce (&ok, &all_ok, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD)) ;
    if (!all_ok)
    {
	/* one of the blocks failed; all processes free result */
	amesos_paraklete_btf_free_numeric (&LU_btf_numeric, Common) ;
    }
    return (LU_btf_numeric) ;
}


/* ========================================================================== */
/* === paraklete_btf_solve ================================================== */
/* ========================================================================== */

Int amesos_paraklete_btf_solve                 /* TRUE if OK, FALSE otherwise */
(
    paraklete_btf_numeric *LU_btf_numeric,	/* numeric factorization */
    paraklete_btf_symbolic *LU_btf_symbolic,	/* symbolic analysis */
    double *B,				/* right-hand-side; soln on output */
    paraklete_common *Common
)
{
    double wk ;
    cholmod_common *cm ;
    KLU_common *km ;
    cholmod_sparse *F ;
    void **LUsymbolic ;
    void **LUnumeric ;
    double *Fx = NULL, *W = NULL, *Singleton ;
    Int *Qbtf, *Pbinv, *Rbtf, *Fp = NULL, *Fi = NULL ;
    Int iold, n, jnew, k, k1, k2, nk, p, nblocks, block ;

    /* ---------------------------------------------------------------------- */
    /* get inputs and allocate workspace */
    /* ---------------------------------------------------------------------- */

    cm = &(Common->cm) ;
    km = &(Common->km) ;

    nblocks = LU_btf_symbolic->nblocks ;
    n = LU_btf_symbolic->n ;
    Rbtf = LU_btf_symbolic->Rbtf ;
    Qbtf = LU_btf_symbolic->Qbtf ;
    Pbinv = LU_btf_symbolic->Pbinv ;
    LUsymbolic = LU_btf_symbolic->LUsymbolic ;
    LUnumeric = LU_btf_numeric->LUnumeric ;
    Singleton = LU_btf_numeric->Singleton ;
    F = LU_btf_numeric->F ;
    if (Common->myid == 0)
    {
        /* only the master process has the off-diagonal blocks F */
        Fp = F->p ;
        Fi = F->i ;
        Fx = F->x ;
    }

    /* ---------------------------------------------------------------------- */
    /* permute right-hand side; W = Pbtf*B */
    /* ---------------------------------------------------------------------- */

    if (Common->myid == 0)
    {
        /* only the master process does this */
        W = CHOLMOD (malloc) (n, sizeof (double), cm) ;
        if (cm->status != CHOLMOD_OK)
        {
            /* TODO return NULL, and broadcast error to all processes */
            PARAKLETE_ERROR (PK_OUT_OF_MEMORY, "out of memory") ;
        }
        for (iold = 0 ; iold < n ; iold++)
        {
            W [Pbinv [iold]] = B [iold] ;
        }
    }

    /* ---------------------------------------------------------------------- */
    /* block backsolve */
    /* ---------------------------------------------------------------------- */

    for (block = nblocks-1 ; block >= 0 ; block--)
    {
	k1 = Rbtf [block] ;
	k2 = Rbtf [block+1] ;
	nk = k2 - k1 ;
	if (nk == 1)
	{
	    /* singleton */
            if (Common->myid == 0)
            {
                W [k1] /= Singleton [block] ;
            }
	}
        else if (nk < 1000)
        {
            if (Common->myid == 0)
            {
                /* use KLU on process 0 */
                KLU_solve (
                    (KLU_symbolic *) LUsymbolic [block],
                    (KLU_numeric  *) LUnumeric  [block], nk, 1, W+k1, km) ;

                if (km->status != KLU_OK)
                {
                    /* TODO return NULL, and broadcast error to all processes */
                    PARAKLETE_ERROR (PK_UNKNOWN, "KLU solve failed") ;
                }

            }
            else
            {
                /* small block; nothing to do for other processors*/
                LUnumeric [block] = NULL ;
            }
        }
	else
	{
	    /* solve the diagonal block */
	    amesos_paraklete_solve (
                (paraklete_numeric *)  LUnumeric  [block], 
                (paraklete_symbolic *) LUsymbolic [block], W+k1, Common) ;
            /* TODO check for error here */
	}
	/* backsolve for off-diagonal entries */
        if (Common->myid == 0)
        {
            for (k = k1 ; k < k2 ; k++)
            {
                wk = W [k] ;
                for (p = Fp [k] ; p < Fp [k+1] ; p++)
                {
                    W [Fi [p]] -= Fx [p] * wk ;
                }
            }
        }
    }

    /* ---------------------------------------------------------------------- */
    /* permute solution; B = Qbtf'*W */
    /* ---------------------------------------------------------------------- */

    if (Common->myid == 0)
    {
        for (jnew = 0 ; jnew < n ; jnew++)
        {
            B [Qbtf [jnew]] = W [jnew] ;
        }
        CHOLMOD (free) (n, sizeof (double), W, cm) ;
    }

    return (TRUE) ;
}


/* ========================================================================== */
/* === paraklete_btf_free_symbolic ========================================== */
/* ========================================================================== */

void amesos_paraklete_btf_free_symbolic
(
    paraklete_btf_symbolic **LU_btf_symbolic_handle,
    paraklete_common *Common
)
{
    paraklete_btf_symbolic *LU_btf_symbolic ;
    void **LUsymbolic ;
    Int n, block, nblocks, k1, k2, nk, *Rbtf ;
    cholmod_common *cm ;
    KLU_common *km ;

    if (LU_btf_symbolic_handle == NULL)
    {
	return ;
    }
    LU_btf_symbolic = *LU_btf_symbolic_handle ; 
    if (LU_btf_symbolic == NULL)
    {
	*LU_btf_symbolic_handle = NULL ;
	return ;
    }
    cm = &(Common->cm) ;
    km = &(Common->km) ;

    LUsymbolic = LU_btf_symbolic->LUsymbolic ;
    Rbtf = LU_btf_symbolic->Rbtf ;
    n = LU_btf_symbolic->n ;
    nblocks = LU_btf_symbolic->nblocks ;
    for (block = 0 ; block < nblocks ; block++)
    {
	k2 = Rbtf [block+1] ;
	k1 = Rbtf [block] ;
	nk = k2 - k1 ;
        if (nk < 1000)
        {
            KLU_symbolic *S ;
            S = (KLU_symbolic *) LUsymbolic [block] ;
            KLU_free_symbolic (&S, &(Common->km)) ;
        }
        else
        {
            paraklete_symbolic *S ;
            S = (paraklete_symbolic *) LUsymbolic [block] ;
            amesos_paraklete_free_symbolic (&S, Common) ;
        }
        LUsymbolic [block] = NULL ;
    }
    CHOLMOD (free) (n, sizeof (void *), LU_btf_symbolic->LUsymbolic, cm) ;
    CHOLMOD (free) (3*n+1, sizeof (Int), LU_btf_symbolic->Mem_n, cm) ;
    *LU_btf_symbolic_handle = CHOLMOD (free) (1, sizeof (paraklete_btf_symbolic),
	LU_btf_symbolic, cm) ;
}


/* ========================================================================== */
/* === paraklete_btf_free_numeric =========================================== */
/* ========================================================================== */

void amesos_paraklete_btf_free_numeric
(
    paraklete_btf_numeric **LU_btf_numeric_handle,
    paraklete_common *Common
)
{
    paraklete_btf_numeric *LU_btf_numeric ;
    void **LUnumeric ;
    Int block, nblocks ;
    cholmod_common *cm ;

    if (LU_btf_numeric_handle == NULL)
    {
	return ;
    }
    LU_btf_numeric = *LU_btf_numeric_handle ; 
    if (LU_btf_numeric == NULL)
    {
	*LU_btf_numeric_handle = NULL ;
	return ;
    }
    cm = &(Common->cm) ;

    LUnumeric = LU_btf_numeric->LUnumeric ;
    nblocks = LU_btf_numeric->nblocks ;

    for (block = 0 ; block < nblocks ; block++)
    {
        paraklete_numeric *N ;
        N = (paraklete_numeric *) LUnumeric [block] ;
        if (N != NULL)
        {
            if (N->magic == PARAKLETE_MAGIC)
            {
                /* this block was factorized with Paraklete */
                amesos_paraklete_free_numeric (&N, Common) ;
            }
            else
            {
                /* this block was factorized with KLU */
                KLU_numeric *K ;
                K = (KLU_numeric *) LUnumeric [block] ;
                KLU_free_numeric (&K, &(Common->km)) ;
            }
        }
        LUnumeric [block] = NULL ;
    }

    CHOLMOD (free) (nblocks, sizeof (double), LU_btf_numeric->Singleton, cm) ;
    CHOLMOD (free) (nblocks, sizeof (void *), LUnumeric, cm) ;
    CHOLMOD (free_sparse) (&(LU_btf_numeric->F), cm) ;

    *LU_btf_numeric_handle = CHOLMOD (free) (1, sizeof (paraklete_btf_numeric),
	LU_btf_numeric, cm) ;
}
