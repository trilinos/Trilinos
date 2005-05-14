/* ========================================================================== */
/* === Partition/cholmod_metis_interface ==================================== */
/* ========================================================================== */

/*
 * CHOLMOD/Partition version 0.1.  May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* CHOLMOD interface to the METIS package (Version 4.0.1):
 *
 * cholmod_metis_bisector:
 *
 *	Interface to METIS_NodeComputeSeparator.  Finds a set of nodes that
 *	partitions the graph into two parts.
 *
 * cholmod_metis_nodend:
 *
 *	Interface to METIS_NodeND, METIS's own nested dissection algorithm.
 *	Typically faster than cholmod_nested_dissection, mostly because it
 *	uses minimum degree on just the leaves of the separator tree, rather
 *	than the whole matrix.
 *
 * Note that METIS does not return an error if it runs out of memory.  Instead,
 * it terminates the program.  This interface attempts to avoid that problem
 * by preallocating space that should be large enough for any memory allocations
 * within METIS, and then freeing that space, just before the call to METIS.
 * While this is not guaranteed to work, it is very unlikely to fail.  If you
 * still encounter a problem, increase Common->metis_memory.  If you don't mind
 * having your program terminated, set Common->metis_memory to zero.  Several
 * other METIS workarounds are made in the routines in this file.
 *
 * FUTURE WORK: an interface to other graph partitioners (CHACO, SCOTCH, ...).
 *
 * workspace: several size-nz and size-n temporary arrays.  Uses no workspace
 * in Common.
 *
 * If CHOLMOD is compiled with -DNMEIS, it makes no calls to routines in this
 * file and this file can also be excluded from your compiled CHOLMOD library.
 */

#include "metis.h"
/* METIS has its own ASSERT that it reveals to the user, so remove it here: */
#undef ASSERT

#include "cholmod_partition.h"
#include "cholmod_internal.h"

/* ========================================================================== */
/* === metis_memory_ok ====================================================== */
/* ========================================================================== */

/* METIS_NodeND and METIS_NodeComputeSeparator will terminate your program it
 * they run out of memory.  In an attempt to workaround METIS' behavior, this
 * routine allocates a single block of memory of size equal to an observed
 * upper bound on METIS' memory usage.  It then immediately deallocates the
 * block.  If the allocation fails, METIS is not called.
 *
 * Median memory usage for a graph with n nodes and nz edges (counting each
 * edge twice, or both upper and lower triangular parts of a matrix) is
 * 4*nz + 40*n + 4096 int's.  A "typical" upper bound is 10*nz + 50*n + 4096
 * int's.  Nearly all matrices tested fit within that upper bound, with the
 * exception two: Schenk_IBMNA/c-64 and Gupta/gupta2.  The latter exceeds the
 * "upper bound" by a factor of just less than 2.
 *
 * If you do not mind having your program terminated if it runs out of memory,
 * set Common->metis_memory to zero.  Its default value is 2, which allows for
 * some memory fragmentation.
 */

#define GUESS(nz,n) (10 * (nz) + 50 * (n) + 4096)

static int metis_memory_ok
(
    int n,
    int nz,
    cholmod_common *Common
)
{
    double s ;
    int t ;
    void *p ;

    if (Common->metis_memory <= 0)
    {
	/* do not prevent METIS from running out of memory */
	return (TRUE) ;
    }

    n  = MAX (1, n) ;
    nz = MAX (0, nz) ;

    /* compute in double, to avoid int overflow */
    s = GUESS ((double) nz, (double) n) ;
    s *= Common->metis_memory ;

    if (s * sizeof (int) >= ((double) INT_MAX))
    {
	/* don't even attempt to malloc such a large block */
	cholmod_error (CHOLMOD_TOO_LARGE, "problem too large", Common) ;
	return (FALSE) ;
    }

    /* recompute in size_t */
    t = GUESS ((size_t) nz, (size_t) n) ;
    t *= Common->metis_memory ;

    /* attempt to malloc the block */
    p = cholmod_malloc (t, sizeof (int), Common) ;
    if (p == NULL)
    {
	/* failure - return out-of-memory condition */
	return (FALSE) ;
    }

    /* success - free the block */
    cholmod_free (p, t, sizeof (int), Common) ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_metis_bisector =============================================== */
/* ========================================================================== */

/* A wrapper for METIS_NodeComputeSeparator. 
 *
 * The input matrix A must be square, symmetric (with both upper and lower
 * parts present) and with no diagonal entries.  These conditions are NOT
 * checked.
 */

long cholmod_metis_bisector   /* returns separator size, or -1 if failure. */
(
    /* inputs, not modified on output */
    cholmod_sparse *A,
    void *Anw_p,	    /* size n, node weights */
    void *Aew_p,	    /* size nz, edge weights */

    /* output, undefined on input */
    void *Partition_p,	    /* size n */

    cholmod_common *Common
)
{
    int *Ap, *Ai, *Anw, *Aew, *Partition ;
    idxtype *Mp, *Mi, *Mnw, *Mew, *Mpart ;
    int n, nleft, nright, j, Opt [8], csep, total_weight, lightest, nz ;
    DEBUG (int nsep) ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (EMPTY) ;
    RETURN_IF_NULL (A, EMPTY) ;
    n = A->nrow ;
    Anw = Anw_p ;
    Aew = Aew_p ;
    Partition = Partition_p ;
    RETURN_IF_NULL (Anw, EMPTY) ;
    RETURN_IF_NULL (Aew, EMPTY) ;
    RETURN_IF_NULL (Partition, EMPTY) ;
    if (A->stype || A->nrow != A->ncol)
    {
	/* A must be square, with both upper and lower parts present */
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_metis_bisector: matrix must be square, symmetric,"
		" and with both upper/lower parts present", Common) ;
	return (EMPTY) ;
    }
    Common->status = CHOLMOD_OK ;
    if (n == 0)
    {
	/* nothing to do */
	return (0) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;
    Ai = A->i ;

    /* ---------------------------------------------------------------------- */
    /* set default options */
    /* ---------------------------------------------------------------------- */

    Opt [0] = 0 ;	/* use defaults */
    Opt [1] = 3 ;	/* matching type */
    Opt [2] = 1 ;	/* init. partitioning algo*/
    Opt [3] = 2 ;	/* refinement algorithm */
    Opt [4] = 0 ;	/* no debug */
    Opt [5] = 0 ;	/* unused */
    Opt [6] = 0 ;	/* unused */
    Opt [7] = -1 ;	/* random seed */

    DEBUG (for (j = 0 ; j < n ; j++) ASSERT (Anw [j] > 0)) ;

    /* ---------------------------------------------------------------------- */
    /* copy int to METIS idxtype, if necessary */
    /* ---------------------------------------------------------------------- */

    nz = Ap [n] ;
    DEBUG (for (j = 0 ; j < nz ; j++) ASSERT (Aew [j] > 0)) ;
    if (sizeof (int) == sizeof (idxtype))
    {
	/* this is the typical case */
	Mi    = Ai ;
	Mew   = Aew ;
	Mp    = Ap ;
	Mnw   = Anw ;
	Mpart = Partition ;
    }
    else
    {
	/* idxtype and int differ; copy the graph into the METIS idxtype */
	Mi    = cholmod_malloc (nz,  sizeof (idxtype), Common) ;
	Mew   = cholmod_malloc (nz,  sizeof (idxtype), Common) ;
	Mp    = cholmod_malloc (n+1, sizeof (idxtype), Common) ;
	Mnw   = cholmod_malloc (n,   sizeof (idxtype), Common) ;
	Mpart = cholmod_malloc (n,   sizeof (idxtype), Common) ;
	if (Common->status < CHOLMOD_OK)
	{
	    cholmod_free (Mi,    nz,  sizeof (idxtype), Common) ;
	    cholmod_free (Mew,   nz,  sizeof (idxtype), Common) ;
	    cholmod_free (Mp,    n+1, sizeof (idxtype), Common) ;
	    cholmod_free (Mnw,   n,   sizeof (idxtype), Common) ;
	    cholmod_free (Mpart, n,   sizeof (idxtype), Common) ;
	    return (EMPTY) ;
	}
	for (j = 0 ; j < nz ; j++) Mi  [j] = Ai  [j] ;
	for (j = 0 ; j < nz ; j++) Mew [j] = Aew [j] ;
	for (j = 0 ; j <= n ; j++) Mp  [j] = Ap  [j] ;
	for (j = 0 ; j <  n ; j++) Mnw [j] = Anw [j] ;
    }

    /* ---------------------------------------------------------------------- */
    /* METIS workaround: try to ensure METIS doesn't run out of memory */
    /* ---------------------------------------------------------------------- */

    if (!metis_memory_ok (n, nz, Common))
    {
	/* METIS might ask for too much memory and thus terminate the program */
	return (EMPTY) ;
    }

    /* ---------------------------------------------------------------------- */
    /* partition the graph */
    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    PRINT1 (("Metis graph, n = %d\n", n)) ;
    for (j = 0 ; j < n ; j++)
    {
	int p ;
	PRINT2 (("M(:,%d) node weight %d\n", j, Mnw [j])) ;
	ASSERT (Mnw [j] > 0) ;
	for (p = Mp [j] ; p < Mp [j+1] ; p++)
	{
	    PRINT3 ((" %d : %d\n", Mi [p], Mew [p])) ;
	    ASSERT (Mi [p] != j) ;
	    ASSERT (Mew [p] > 0) ;
	}
    }
#endif

    METIS_NodeComputeSeparator (&n, Mp, Mi, Mnw, Mew, Opt, &csep, Mpart) ;

    PRINT1 (("METIS csep %d\n", csep)) ;

    /* ---------------------------------------------------------------------- */
    /* copy the results back from idxtype, if required */
    /* ---------------------------------------------------------------------- */

    if (sizeof (int) != sizeof (idxtype))
    {
	for (j = 0 ; j < n ; j++) Partition [j] = Mpart [j] ;
	cholmod_free (Mi,    nz,  sizeof (idxtype), Common) ;
	cholmod_free (Mew,   nz,  sizeof (idxtype), Common) ;
	cholmod_free (Mp,    n+1, sizeof (idxtype), Common) ;
	cholmod_free (Mnw,   n,   sizeof (idxtype), Common) ;
	cholmod_free (Mpart, n,   sizeof (idxtype), Common) ;
    }

    /* ---------------------------------------------------------------------- */
    /* ensure a reasonable separator */
    /* ---------------------------------------------------------------------- */

    /* METIS can return a valid separator with no nodes in (for example) the
     * left part.  In this case, there really is no separator.  CHOLMOD
     * prefers, in this case, for all nodes to be in the separator (and both
     * left and right parts to be empty).  Also, if the graph is unconnected,
     * METIS can return a valid empty separator.  CHOLMOD prefers at least one
     * node in the separator.  Note that cholmod_nested_dissection only calls
     * this routine on connected components, but cholmod_bisect can call this
     * routine for any graph. */

    if (csep == 0)
    {
	/* The separator is empty, select lightest node as separator.  If
	 * ties, select the highest numbered node. */
	lightest = 0 ;
	for (j = 0 ; j < n ; j++)
	{
	    if (Anw [j] <= Anw [lightest])
	    {
		lightest = j ;
	    }
	}
	PRINT1 (("Force %d as sep\n", lightest)) ;
	Partition [lightest] = 2 ;
	csep = Anw [lightest] ;
    }

    /* determine the node weights in the left and right part of the graph */
    nleft = 0 ;
    nright = 0 ;
    DEBUG (nsep = 0) ;
    for (j = 0 ; j < n ; j++)
    {
	PRINT1 (("Partition [%d] = %d\n", j, Partition [j])) ;
	if (Partition [j] == 0)
	{
	    nleft += Anw [j] ;
	}
	else if (Partition [j] == 1)
	{
	    nright += Anw [j] ;
	}
#ifndef NDEBUG
	else
	{
	    ASSERT (Partition [j] == 2) ;
	    nsep += Anw [j] ;
	}
#endif
    }
    ASSERT (csep == nsep) ;

    total_weight = nleft + nright + csep ;

    if (csep < total_weight)
    {
	/* The separator is less than the whole graph.  Make sure the left and
	 * right parts are either both empty or both non-empty. */
	PRINT1 (("nleft %d nright %d csep %d tot %d\n",
		nleft, nright, csep, total_weight)) ;
	ASSERT (nleft + nright + csep == total_weight) ;
	ASSERT (nleft > 0 || nright > 0) ;
	if ((nleft == 0 && nright > 0) || (nleft > 0 && nright == 0))
	{
	    /* left or right is empty; put all nodes in the separator */
	    PRINT1 (("Force all in sep\n")) ;
	    csep = total_weight ;
	    for (j = 0 ; j < n ; j++)
	    {
		Partition [j] = 2 ;
	    }
	}
    }

    ASSERT (cholmod_dump_partition (n, Ap, Ai, Anw, Partition, csep)) ;

    /* ---------------------------------------------------------------------- */
    /* return the sum of the weights of nodes in the separator */
    /* ---------------------------------------------------------------------- */

    return (csep) ;
}


/* ========================================================================== */
/* === cholmod_metis ======================================================== */
/* ========================================================================== */

/* CHOLMOD wrapper for the METIS_NodeND ordering routine.  Creates A+A',
 * A*A' or A(:,f)*A(:,f)' and then calls METIS_NodeND on the resulting graph.
 * This routine is comparable to cholmod_nested_dissection, except that it
 * calls METIS_NodeND directly, and it does not return the separator tree.
 *
 * workspace:  Flag (nrow), Iwork (max (nrow,ncol))
 *	Allocates a temporary matrix B=A*A' or B=A.
 */

int cholmod_metis	/* returns TRUE if successful, FALSE otherwise */
(
    cholmod_sparse *A,
    void *fset,
    size_t fsize,

    /* outputs, contents undefined on input */
    void *Perm_p,	/* size n.  Perm [k] = j if node j is kth node in the
			 * permuted matrix */
    cholmod_common *Common
)
{

    double d ;
    int *Perm, *Iperm, *Iwork, *Bp, *Bi ;
    idxtype *Mp, *Mi, *Mperm, *Miperm ;
    int i, j, n, nz, p, Opt [8], zero = 0, identity ;
    cholmod_sparse *B ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    Perm = Perm_p ;
    RETURN_IF_NULL (Perm, FALSE) ;
    n = A->nrow ;
    Common->status = CHOLMOD_OK ;
    if (n == 0)
    {
	return (TRUE) ;	    /* nothing to do */
    }

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    cholmod_allocate_work (n, MAX (n, ((int) A->ncol)), 0, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;
    }
    ASSERT (cholmod_dump_work (TRUE, TRUE, 0, Common)) ;

    /* ---------------------------------------------------------------------- */
    /* convert the matrix to adjacency list form */
    /* ---------------------------------------------------------------------- */

    /* The input graph for METIS must be symmetric, with both upper and lower
     * parts present, and with no diagonal entries.  The columns need not be
     * sorted.
     * B = A+A', A*A', or A(:,f)*A(:,f)', upper and lower parts present */
    if (A->stype)
    {
	/* Add the upper/lower part to a symmetric lower/upper matrix by
	 * converting to unsymmetric mode */
	/* workspace: Iwork (max (nrow,ncol)) */
	B = cholmod_copy (A, 0, -1, Common) ;
    }
    else
    {
	/* B = A*A' or A(:,f)*A(:,f)', no diagonal */
	/* workspace: Flag (nrow), Iwork (max (nrow,ncol)) */
	B = cholmod_aat (A, fset, fsize, -1, Common) ;
    }
    ASSERT (cholmod_dump_sparse (B, "B for NodeND", Common) >= 0) ;
    PRINT1 (("NodeND\n")) ;
    if (Common->status < CHOLMOD_OK)
    {
	PRINT1 (("create B failed\n")) ;
	return (FALSE) ;
    }
    ASSERT (B->nrow == A->nrow) ;

    /* ---------------------------------------------------------------------- */
    /* find the permutation */
    /* ---------------------------------------------------------------------- */

    Iwork = Common->Iwork ;
    Iperm = Iwork ;		/* size n (i/i/l) */

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    Bp = B->p ;
    Bi = B->i ;
    nz = Bp [n] ;

    /* ---------------------------------------------------------------------- */
    /* set control parameters for METIS_NodeND */
    /* ---------------------------------------------------------------------- */

    Opt [0] = 0 ;	/* use defaults */
    Opt [1] = 3 ;	/* matching type */
    Opt [2] = 1 ;	/* init. partitioning algo*/
    Opt [3] = 2 ;	/* refinement algorithm */
    Opt [4] = 0 ;	/* no debug */
    Opt [5] = 1 ;	/* initial compression */
    Opt [6] = 0 ;	/* no dense node removal */
    Opt [7] = 1 ;	/* number of separators @ each step */

    /* ---------------------------------------------------------------------- */
    /* allocate the METIS input arrays, if needed */
    /* ---------------------------------------------------------------------- */

    if (sizeof (int) == sizeof (idxtype))
    {
	/* This is the typical case. */
	Miperm = Iperm ;
	Mperm  = Perm ;
	Mp     = Bp ;
	Mi     = Bi ;
    }
    else
    {
	/* allocate graph for METIS only if int and idxtype differ */
	Miperm = cholmod_malloc (n,   sizeof (idxtype), Common) ;
	Mperm  = cholmod_malloc (n,   sizeof (idxtype), Common) ;
	Mp     = cholmod_malloc (n+1, sizeof (idxtype), Common) ;
	Mi     = cholmod_malloc (nz,  sizeof (idxtype), Common) ;
	if (Common->status < CHOLMOD_OK)
	{
	    /* out of memory */
	    cholmod_free (Miperm, n,   sizeof (idxtype), Common) ;
	    cholmod_free (Mperm , n,   sizeof (idxtype), Common) ;
	    cholmod_free (Mp    , n+1, sizeof (idxtype), Common) ;
	    cholmod_free (Mi    , nz,  sizeof (idxtype), Common) ;
	    return (FALSE) ;
	}
	for (j = 0 ; j <= n ; j++)
	{
	    Mp [j] = Bp [j] ;
	}
	for (p = 0 ; p < nz ; j++)
	{
	    Mi [p] = Bi [p] ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* METIS workarounds */
    /* ---------------------------------------------------------------------- */

    identity = FALSE ;
    if (nz == 0)
    {
	/* The matrix has no off-diagonal entries.  METIS_NodeND fails in this
	 * case, so avoid using it.  The best permutation is identity anyway,
	 * so this is an easy fix. */
	identity = TRUE ;
	PRINT1 (("no nz\n")) ;
    }
    else if (Common->metis_nswitch > 0)
    {
	/* METIS_NodeND in METIS 4.0.1 gives a seg fault with one matrix of
	 * order n = 3005 and nz = 6,036,025, including the diagonal entries.
	 * The workaround is to return the identity permutation instead of using
	 * METIS for matrices of dimension 3000 or more and with density of 66%
	 * or more - admittedly an uncertain fix, but such matrices are so dense
	 * that any reasonable ordering will do, even identity (n^2 is only 50%
	 * higher than nz in this case).  CHOLMOD's nested dissection method
	 * (cholmod_nested_dissection) has no problems with the same matrix,
	 * even though it too uses METIS_NodeComputeSeparator.  The matrix is
	 * derived from LPnetlib/lpi_cplex1.  If C is the lpi_cplex matrix (of
	 * order 3005-by-5224), A = (C*C')^2 results in the seg fault.  The seg
	 * fault also occurs in the stand-alone onmetis program that comes with
	 * METIS. If a future version of METIS fixes this problem, then set 
	 * Common->metis_nswitch to zero.
	 */
	d = ((double) nz) / (((double) n) * ((double) n)) ;
	if (n > (int) (Common->metis_nswitch) && d > Common->metis_dswitch)
	{
	    identity = TRUE ;
	    PRINT1 (("nswitch/dswitch activated\n")) ;
	}
    }

    if (!identity && !metis_memory_ok (n, nz, Common))
    {
	/* METIS might ask for too much memory and thus terminate the program */
	identity = TRUE ;
	PRINT1 (("no metis called\n")) ;
    }

    /* ---------------------------------------------------------------------- */
    /* find the permutation */
    /* ---------------------------------------------------------------------- */

    if (identity)
    {
	for (i = 0 ; i < n ; i++)
	{
	    Mperm [i] = i ;
	}
	PRINT1 (("identity\n")) ;
    }
    else
    {
	PRINT1 (("calling metis\n")) ;
	METIS_NodeND (&n, Mp, Mi, &zero, Opt, Mperm, Miperm) ;
	PRINT1 (("metis done\n")) ;
    }

    /* ---------------------------------------------------------------------- */
    /* free the METIS input arrays and return the result */
    /* ---------------------------------------------------------------------- */

    if (sizeof (int) != sizeof (idxtype))
    {
	for (i = 0 ; i < n ; i++)
	{
	    Perm [i] = (int) (Mperm [i]) ;
	}
	cholmod_free (Miperm, n,   sizeof (idxtype), Common) ;
	cholmod_free (Mperm , n,   sizeof (idxtype), Common) ;
	cholmod_free (Mp    , n+1, sizeof (idxtype), Common) ;
	cholmod_free (Mi    , nz,  sizeof (idxtype), Common) ;
    }

    cholmod_free_sparse (&B, Common) ;
    ASSERT (cholmod_dump_work (TRUE, TRUE, 0, Common)) ;
    PRINT1 (("cholmod_metis done\n")) ;
    return (TRUE) ;
}
