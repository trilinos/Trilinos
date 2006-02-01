/* ========================================================================== */
/* === Paraklete/pk ========================================================= */
/* ========================================================================== */

/* Demo program for Paraklete.
 *
 * Usage:
 *
 *	pk matrix rhs id tries
 *
 * id and tries are used for malloc-error testing.  Not used if not present.
 * id=-1 means all processes use tries.  If id or tries are missing, then
 *
 * matrix and rhs are the files used for A and b.  If b is "-", a default
 * b = 0:n-1 is used.  If matrix is not present, stdin is used.
 *
 * PARAKLETE version 0.1: parallel sparse LU factorization.  May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 *
 * To turn on debugging, edit the Cholmod/Include/cholmod_internal.h file and
 * uncomment the "#undef NDEBUG".
 */

#include "pk.h"

/* ========================================================================== */
/* === read_triplet ========================================================= */
/* ========================================================================== */

cholmod_triplet *read_triplet
(
    FILE *f,
    cholmod_common *cm
)
{
    cholmod_triplet *T ;
    double *Tx ;
    int *Ti, *Tj ;
    int n, k, nrow, ncol, nz, stype ;

    /* ---------------------------------------------------------------------- */
    /* read in a triplet matrix from a file */
    /* ---------------------------------------------------------------------- */

    if (f == NULL)
    {
	return (NULL) ;
    }

    if (fscanf (f, "%d %d %d %d\n", &nrow, &ncol, &nz, &stype) == EOF)
    {
	return (NULL) ;
    }

    printf ("nrow %d ncol %d nz %d stype %d\n", nrow, ncol, nz, stype) ;

    n = MAX (nrow, ncol) ;
    if (stype != 0)
    {
	cholmod_error (CHOLMOD_INVALID, "pk: stype must be 0", cm);
	return (NULL) ;
    }

    T = cholmod_allocate_triplet (nrow, ncol, nz, stype, TRUE, cm) ;
    if (T == NULL)
    {
	cholmod_error (CHOLMOD_INVALID, "pk: cannot create triplet matrix", cm);
	return (NULL) ;
    }
    Ti = T->i ;
    Tj = T->j ;
    Tx = T->x ;

    for (k = 0 ; k < nz ; k++)
    {
	if (fscanf (f, "%d %d %lg\n", Ti+k, Tj+k, Tx+k) == EOF)
	{
	    cholmod_error (CHOLMOD_INVALID, "pk: cannot read triplet", cm);
	    return (NULL) ;
	}
    }
    T->nnz = nz ;
    return (T) ;
}


void my_handler (int status, char *msg)
{
    printf ("Error handler: %d: %s\n", status, msg) ;
    if (status != CHOLMOD_OK)
    {
	fprintf (stderr, "\n\n*********************************************\n");
	fprintf (stderr, "**** Test failure: %d %s\n", status, msg) ;
	fprintf (stderr, "*********************************************\n\n");
	fflush (stderr) ;
	fflush (stdout) ;
    }
}

/* ========================================================================== */
/* === FINISH =============================================================== */
/* ========================================================================== */

/* Free workspace, matrices, and factorization.  Finalize CHOLMOD and MPI.
 * Finally, print "1" to the file "done", so that MATLAB can read this file in
 * to see if the code finished correctly. */

#define FINISH \
{ \
    cholmod_free_dense (&R, cm) ; \
    cholmod_free_dense (&B, cm) ; \
    cholmod_free_triplet (&T, cm) ; \
    cholmod_free_sparse (&A, cm) ; \
    paraklete_free_symbolic (&LUsymbolic, Common) ; \
    paraklete_free_numeric (&LU, Common) ; \
    if (cm->file != stdout) \
    { \
	fclose (cm->file) ; \
    } \
    cholmod_finish (cm) ; \
    if (cm->malloc_count != 0 || cm->memory_inuse != 0) \
    { \
	printf ("finish: %d %d\n", cm->malloc_count, cm->memory_inuse) ; \
	MPI (MPI_Finalize ( )) ; \
	exit (1) ; \
    } \
    ASSERT (cm->malloc_count == 0 && cm->memory_inuse == 0) ; \
    if (myid == 0) \
    { \
	f = fopen ("done", "w") ; \
	fprintf (f, "1\n") ; \
	fclose (f) ; \
    } \
    MPI (MPI_Finalize ( )) ; \
    return (0) ; \
}


/* ========================================================================== */
/* === pk =================================================================== */
/* ========================================================================== */

int PK_main(int argc, char **argv)
{
    cholmod_sparse *A ;
    paraklete_common *Common, pcommon ;
    paraklete_symbolic *LUsymbolic ;
    paraklete_numeric *LU ;
    cholmod_common *cm ;
    cholmod_dense *R, *B ;
    cholmod_scalar one, minusone ;
    cholmod_triplet *T ;
    double *Bx ;
    int *P, *Q ;
    int i, n, nproc, myid, ok ;
    FILE *f ;
    /*
    char filename [100] ;
    */

    MPI (MPI_Init (NULL, NULL)) ;

    /* ---------------------------------------------------------------------- */

    R = NULL ;
    B = NULL ;
    T = NULL ;
    A = NULL ;
    LUsymbolic = NULL ;
    LU = NULL ;
    n = 0 ;

    /* ---------------------------------------------------------------------- */
    /* start MPI */
    /* ---------------------------------------------------------------------- */

    myid = 0 ;
    nproc = 1 ;
#ifdef HACK
    printf ("hacked case\n") ;
    nproc = HACK ;
#endif

    MPI (MPI_Comm_size (MPI_COMM_WORLD, &nproc)) ;
    MPI (MPI_Comm_rank (MPI_COMM_WORLD, &myid)) ;

    if (myid == 0)
    {
	printf ("\n-------------------------------------------------------\n") ;
	printf ("mpirun -np %d ", nproc) ;
	for (i = 0 ; i < argc ; i++)
	{
	    printf ("%s ", argv [i]) ;
	}
	printf ("\n")  ;
    }

    MPI (MPI_Barrier (MPI_COMM_WORLD)) ;

    printf ("proc %d pid %d\n", myid, getpid ( )) ;

    DEBUG (cholmod_dump_malloc = FALSE) ;

    /* ---------------------------------------------------------------------- */
    /* initialize workspace and parameters */
    /* ---------------------------------------------------------------------- */

    one.x = 1 ;
    one.z = 0 ;
    minusone.x = -1 ;
    minusone.z = 0 ;

    Common = &pcommon ;
    cm = &(Common->cm) ;
    cholmod_start (CHOLMOD_INT, CHOLMOD_REAL, CHOLMOD_DOUBLE, cm) ;
    DEBUG_INIT ("pk") ;
    Common->nproc = nproc ;
    Common->myid = myid ;
    cm->print = 1 ;
    cm->precise = TRUE ;
    cm->ieee = FALSE ;
    cm->blas_conform = FALSE ;
    cm->error_handler = my_handler ;
    my_tries = -1 ;

    Common->tol_diag = 0.001 ;
    Common->tol_offdiag = 0.1 ;
    Common->growth = 2. ;
#ifndef NDEBUG
    Common->dump = cholmod_dump ;
#else
    Common->dump = 0 ;
#endif

    /* ---------------------------------------------------------------------- */
    /* get input matrix */
    /* ---------------------------------------------------------------------- */

    ok = TRUE ;
    if (myid == 0)
    {
	if (argc > 1)
	{
	    f = fopen (argv [1], "r") ;
	}
	else
	{
	    f = stdin ;
	}
	ok = (f != NULL) ;
	if (ok)
	{
	    T = read_triplet (f, cm) ;
	    A = cholmod_triplet_to_sparse (T, cm) ;
	    cholmod_free_triplet (&T, cm) ;
	    if (cm->status != CHOLMOD_OK)
	    {
		printf ("Failed to read input matrix\n") ;
		ok = FALSE ;
	    }
	    if (argc > 1)
	    {
		fclose (f) ;
	    }
	}
	else
	{
	    printf ("Unable to open matrix file: %s\n", argv [1]) ;
	}
    }

#if 0
    sprintf (filename, "proc%d", myid) ;
#ifdef NMPI
    cm->file = stdout ;
#else
    cm->file = fopen (filename, "w") ;
    ok = ok && (cm->file != NULL) ;
#endif
#endif

    cm->file = stdout ;

    /* ---------------------------------------------------------------------- */
    /* load or construct B */
    /* ---------------------------------------------------------------------- */

    Bx = NULL ;
    if (ok && myid == 0)
    {
	n = A->nrow ;
	B = cholmod_allocate_dense (n, 1, n, CHOLMOD_NONE, cm) ;
	if (cm->status != CHOLMOD_OK)
	{
	    /* out of memroy */
	    printf ("Failed to create right-hand-side\n") ;
	    ok = FALSE ;
	}
	if (ok)
	{
	    Bx = B->x ;
	    if (argc > 2 && argv [2][0] != '-')
	    {
		f = fopen (argv [2], "r") ;
		if (f == NULL)
		{
		    printf ("Failed to read right-hand-side: %s\n", argv [2]) ;
		    ok = FALSE ;
		}
		if (ok)
		{
		    for (i = 0 ; i < n ; i++)
		    {
			fscanf (f, "%lg\n", &Bx [i]) ;
		    }
		    fclose (f) ;
		}
	    }
	    else
	    {
		for (i = 0 ; i < n ; i++)
		{
		    Bx [i] = i ;
		}
	    }
	    R = cholmod_copy_dense (B, n, cm) ;
	    if (cm->status != CHOLMOD_OK)
	    {
		/* out of memroy */
		ok = FALSE ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* get memory handler options */
    /* ---------------------------------------------------------------------- */

    /* Usage: pk matrix rhs id tries, id=-1 means all processes use tries */

    i = 0 ;
    if (argc > 3)
    {
	i = atoi (argv [3]) ;
    }
    if (argc > 4)
    {
	if (i == -1 || i == myid)
	{
	    test_memory_handler (cm) ;
	    my_tries = atoi (argv [4]) ;
	    printf ("proc %d will try malloc %d times\n", myid, my_tries) ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* determine if initializations succeeded */
    /* ---------------------------------------------------------------------- */

    MPI (MPI_Bcast (&ok, 1, MPI_INT, 0, MPI_COMM_WORLD)) ;
    if (!ok)
    {
	printf ("proc %d: initializations failed\n", myid) ;
	FINISH ;
    }

    /* ---------------------------------------------------------------------- */
    /* analyze, factorize, and solve */
    /* ---------------------------------------------------------------------- */

    /* analyze */
    LUsymbolic = paraklete_analyze (A, Common) ;
    if (LUsymbolic == NULL)
    {
	/* out of memory; all processors get a NULL return value */
	printf ("proc %d: paraklete_analyze failed\n", myid) ;
	FINISH ;
    }

    /* factorize L*U = P*A*Q */
    LU = paraklete_factorize (A, LUsymbolic, Common) ;
    if (LU == NULL)
    {
	/* out of memory; all processors get a NULL return value */
	printf ("proc %d: paraklete_factorize failed\n", myid) ;
	FINISH ;
    }

    /* solve Ly=Pb then UQ'x=y, overwritting B with the final solution X */
    paraklete_solve (LU, LUsymbolic, Bx, Common) ;

    /* ---------------------------------------------------------------------- */
    /* compute residual */
    /* ---------------------------------------------------------------------- */

    normal_memory_handler (cm) ;

    if (myid == 0)
    {

#ifndef NMATRIXOPS
	/* compute the residual */
	double r ;
	cholmod_sdmult (A, FALSE, one, minusone, B, R, cm) ;
	r = cholmod_norm_dense (R, 1, cm) ;
	printf ("resid %g\n", r) ;
#endif

	/* print the solution to the file "x" */
	f = fopen ("x", "w") ;
	if (f != NULL)
	{
	    for (i = 0 ; i < n ; i++)
	    {
		fprintf (f, "%30.20e\n", Bx [i]) ;
	    }
	    fclose (f) ;
	}

	/* print P */
	P = LU->P ;
	f = fopen ("P", "w") ;
	if (f != NULL)
	{
	    for (i = 0 ; i < n ; i++)
	    {
		fprintf (f, "%d\n", P [i]) ;
	    }
	    fclose (f) ;
	}

	/* print Q */
	Q = LU->Q ;
	f = fopen ("Q", "w") ;
	if (f != NULL)
	{
	    for (i = 0 ; i < n ; i++)
	    {
		fprintf (f, "%d\n", Q [i]) ;
	    }
	    fclose (f) ;
	}

    }

    /* ---------------------------------------------------------------------- */
    /* done: free everything, finalize CHOLMOD and MPI */
    /* ---------------------------------------------------------------------- */

    FINISH ;
}
