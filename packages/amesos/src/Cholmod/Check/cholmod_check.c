/* ========================================================================== */
/* === Check/cholmod_check ================================================== */
/* ========================================================================== */

/*
 * CHOLMOD/Check version 0.1, May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* Routines to check and print the contents of the 5 CHOLMOD objects:
 *
 * No CHOLMOD routine calls the check or print routines.  If a user wants to
 * check CHOLMOD's input parameters, a separate call to the appropriate check
 * routine should be used before calling other CHOLMOD routines.
 *
 * cholmod_check_common		check statistics and workspace in Common
 * cholmod_check_sparse		check sparse matrix in compressed column form
 * cholmod_check_dense		check dense matrix
 * cholmod_check_factor		check factorization
 * cholmod_check_triplet	check sparse matrix in triplet form
 *
 * cholmod_print_common		print statistics in Common
 * cholmod_print_sparse		print sparse matrix in compressed column form
 * cholmod_print_dense		print dense matrix
 * cholmod_print_factor		print factorization
 * cholmod_print_triplet	print sparse matrix in triplet form
 *
 * In addition, this file contains routines to check and print three types of
 * integer vectors:
 * 
 * cholmod_check_perm		check a permutation of 0:n-1 (no duplicates)
 * cholmod_check_subset		check a subset of 0:n-1 (duplicates OK)
 * cholmod_check_parent		check an elimination tree
 *
 * cholmod_print_perm		print a permutation
 * cholmod_print_subset		print a subset
 * cholmod_print_parent		print an elimination tree
 *
 * Each Common->print level prints the items at or below the given level:
 *
 *	0: nothing
 *	1: error messages
 *	2: warning messages
 *	3: one-line summary of each object printed
 *	4: short summary of each object (first and last few entries)
 *	5: entire contents of the object
 *
 * No CHOLMOD routine calls these routines, so no printing occurs unless
 * the user specifically calls a cholmod_print_* routine.  Thus, the default
 * print level is 3.
 *
 * Common->precise controls the # of digits printed for numerical entries
 * (5 if FALSE, 15 if TRUE).
 *
 * If NPRINT is defined at compile time, then no printing occurs and the
 * standard I/O library is not accessed at all (stdio.h is not included).  The
 * cholmod_check_* and cholmod_print_* routines still check their inputs and
 * return TRUE/FALSE if the object is valid or not.  (Note, however, that the
 * Partition module still includes stdio.h, because METIS includes stdio.h).
 *
 * This file also includes debugging routines that are enabled only when
 * NDEBUG is defined in cholmod_internal.h (cholmod_dump_*).
 */

#include "cholmod_check.h"
#include "cholmod_internal.h"

/* ========================================================================== */
/* === printing definitions ================================================= */
/* ========================================================================== */

/* If NPRINT is defined at compile time, then no CHOLMOD routine prints
 * anything, even the cholmod_print_* routines.
 */

#ifndef NPRINT
#ifdef MATLAB_MEX_FILE
/* MATLAB mexFunction: use mexPrintf */
#define PR(format,arg) { mexPrintf (format, arg) ; }
#define P1(format,arg) { if (print >= 1) mexPrintf (format, arg) ; }
#define P2(format,arg) { if (print >= 2) mexPrintf (format, arg) ; }
#define P3(format,arg) { if (print >= 3) mexPrintf (format, arg) ; }
#define P4(format,arg) { if (print >= 4) mexPrintf (format, arg) ; }
#else
/* ANSI C fprintf */
#include <stdio.h>
#define PR(format,arg) { printf (format, arg) ; }
#define FP ((Common->file == NULL) ? stdout : Common->file)
#define P1(format,arg) { if (print >= 1) fprintf (FP, format, arg) ; }
#define P2(format,arg) { if (print >= 2) fprintf (FP, format, arg) ; }
#define P3(format,arg) { if (print >= 3) fprintf (FP, format, arg) ; }
#define P4(format,arg) { if (print >= 4) fprintf (FP, format, arg) ; }
#endif
#else
/* no printing at all */
#define P1(format,arg)
#define P2(format,arg)
#define P3(format,arg)
#define P4(format,arg)
#endif

#define ERROR(msg) \
{ \
    P1 ("\nCHOLMOD ERROR: %s: ", type) ; \
    if (name != NULL) \
    { \
	P1 ("%s", name) ; \
    } \
    P1 (": %s\n", msg) ; \
    cholmod_error (CHOLMOD_INVALID, "invalid", Common) ; \
    return (FALSE) ; \
}

#ifndef NPRINT

/* print a numerical value */
#define PRINTVALUE(value) \
{ \
    if (precise) \
    { \
	P4 (" % 23.15e", value) ; \
    } \
    else \
    { \
	P4 (" % .5g", value) ; \
    } \
}

/* start printing */
#define ETC_START(count,limit) \
{ \
    count = (init_print == 4) ? (limit) : (-1) ; \
}

/* re-enable printing if condition is met */
#define ETC_ENABLE(condition,count,limit) \
{ \
    if ((condition) && init_print == 4) \
    { \
	count = limit ; \
	print = 4 ; \
    } \
}

/* turn off printing if limit is reached */
#define ETC_DISABLE(count) \
{ \
    if ((count >= 0) && (count-- == 0) && print == 4) \
    { \
	P4 ("%s", "    ...\n")  ; \
	print = 3 ; \
    } \
}

/* re-enable printing, or turn if off after limit is reached */
#define ETC(condition,count,limit) \
{ \
    ETC_ENABLE (condition, count, limit) ; \
    ETC_DISABLE (count) ; \
}

#else
/* no printing */
#define PRINTVALUE(value) ;
#define ETC_START(count,limit) ;
#define ETC_ENABLE(condition,count,limit) ;
#define ETC_DISABLE(count) ;
#endif

#define BOOLEAN(x) ((x) ? "true " : "false")

/* ========================================================================== */
/* === cholmod_check_common ================================================= */
/* ========================================================================== */

/* Print and verify the contents of Common */

static int check_common
(
    int print,
    char *name,
    cholmod_common *Common
)
{
    double fl, lnz ;
    double *Xwork ;
    int *Flag, *Head ;
    long mark ;
    int i, nrow, nmethods, ordering, wsize, amd_printed ;
    char *type = "common" ;

    if (Common == NULL)
    {
	/* this macro does not access Common to determine the print level
	 * or the file to print to. */
	PR ("%s", "\nCHOLMOD ERROR common: ") ;
	if (name != NULL)
	{
	    PR ("%s:", name) ;
	}
	PR ("%s", " null\n") ;
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* print control parameters and statistics */
    /* ---------------------------------------------------------------------- */

    P2 ("%s", "\n") ;
    P1 ("%s", "CHOLMOD common:  ") ;
    if (name != NULL)
    {
	P1 ("%s: ", name) ;
    }
    switch (Common->status)
    {

	case CHOLMOD_OK:
	    P1 ("%s", "status: OK\n") ;
	    break ;

	case CHOLMOD_OUT_OF_MEMORY:
	    P1 ("%s", "status: ERROR, out of memory\n") ;
	    break ;

	case CHOLMOD_INVALID:
	    P1 ("%s", "status: ERROR, invalid parameter\n") ;
	    break ;

	case CHOLMOD_TOO_LARGE:
	    P1 ("%s", "status: ERROR, problem too large\n") ;
	    break ;

	case CHOLMOD_NOT_INSTALLED:
	    P1 ("%s", "status: ERROR, method not installed\n") ;
	    break ;

	case CHOLMOD_NOT_POSDEF:
	    P1 ("%s", "status: warning, matrix not positive definite\n") ;
	    break ;

	default:
	    ERROR ("unknown status") ;
    }

    if (Common->fl != EMPTY)
    {
	P2 ("%s", "  Results from most recent analysis:\n") ;
	P2 ("    Cholesky flop count:        %.5g\n", Common->fl) ;
	P2 ("    Nonzeros in L:              %.5g\n", Common->lnz) ;
    }
    if (Common->modfl != EMPTY)
    {
	P2 ("    Update/downdate flop count: %.5g\n", Common->modfl) ;
    }

    P2 ("  memory blocks in use:    %8d\n", Common->malloc_count) ;
    P2 ("  memory in use (MB):      %8.1f\n", 
	(double) (Common->memory_inuse) / 1048576.) ;
    P2 ("  peak memory usage (MB):  %8.1f\n", 
	(double) (Common->memory_usage) / 1048576.) ;

    /* ---------------------------------------------------------------------- */
    /* primary control parameters and related ordering statistics */
    /* ---------------------------------------------------------------------- */

    P3 ("  maxrank:    update/downdate rank:   %d\n", Common->maxrank) ;
    P3 ("  supernodal: do supernodal Cholesky: %s\n",
	    BOOLEAN (Common->supernodal)) ;
    P3 ("%s", "  nmethods:   number of ordering methods to try: ") ;

    nmethods = MIN (Common->nmethods, CHOLMOD_MAXMETHODS) ;
    nmethods = MAX (nmethods, 0) ;
    P3 ("%d\n", nmethods) ;
    amd_printed = FALSE ;
    for (i = 0 ; i < nmethods ; i++)
    {
	P3 ("    method %d: ", i) ;
	ordering = Common->method [i].ordering ;
	fl = Common->method [i].fl ;
	lnz = Common->method [i].lnz ;
	switch (ordering)
	{

	    case CHOLMOD_NATURAL:
		P3 ("%s", "natural\n") ;
		break ;

	    case CHOLMOD_GIVEN:
		P3 ("%s", "user permutation (if given)\n") ;
		break ;

	    case CHOLMOD_AMD:
		P3 ("%s", "AMD (or COLAMD if factorizing AA')\n") ;
		amd_printed = TRUE ;
		break ;

	    case CHOLMOD_METIS:
		P3 ("%s", "METIS_NodeND nested dissection\n") ;
		break ;

	    case CHOLMOD_ND:
		P3 ("%s", "CHOLMOD nested dissection\n") ;
		if (Common->method [i].nd_prune < 0)
		{
		    P3 ("        nd_prune: for pruning dense nodes:   %s\n",
			    " none pruned") ;
		}
		else
		{
		    P3 ("        nd_prune: for pruning dense nodes:   "
			"%.5g\n",
			Common->method [i].nd_prune) ;
		    P3 ("        a dense node has degree "
			    ">= max(16,(%.5g)*sqrt(n))\n",
			Common->method [i].nd_prune) ;
		}
		P3 ("        nd_small: # nodes in uncut subgraph: %d\n",
			Common->method [i].nd_small) ;
		P3 ("        nd_compress: compress the graph:     %s\n",
			BOOLEAN (Common->method [i].nd_compress)) ;
		P3 ("        nd_camd: use constrained min degree: %s\n",
			BOOLEAN (Common->method [i].nd_camd)) ;
		break ;

	    default:
		P3 ("%d", ordering) ;
		ERROR ("unknown ordering method") ;

	}
	if (fl  != EMPTY) P3 ("        flop count: %.5g\n", fl) ;
	if (lnz != EMPTY) P3 ("        nnz(L):     %.5g\n", lnz) ;
    }

    /* backup AMD results, if any */
    if (!amd_printed)
    {
	P3 ("%s", "    backup method: ") ;
	P3 ("%s", "AMD (or COLAMD if factorizing AA')\n") ;
	fl = Common->method [nmethods].fl ;
	lnz = Common->method [nmethods].lnz ;
	if (fl  != EMPTY) P3 ("        AMD flop count: %.5g\n", fl) ;
	if (lnz != EMPTY) P3 ("        AMD nnz(L):     %.5g\n", lnz) ;
    }

    /* ---------------------------------------------------------------------- */
    /* arcane control parameters */
    /* ---------------------------------------------------------------------- */

    P4 ("%s", "  final_ftype: ") ;
    switch (Common->final_ftype)
    {
	case CHOLMOD_LDL_PACKED:   P4 ("%s", "convert to LDL' packed") ; break ;
	case CHOLMOD_LDL_UNPACKED: P4 ("%s", "convert to LDL' unpacked") ;break;
	case CHOLMOD_LDL_DYNAMIC:  P4 ("%s", "convert to LDL' dynamic") ; break;
	case CHOLMOD_LL_PACKED:	   P4 ("%s", "convert to LL' packed") ; break ;
	default:		   P4 ("%s", "leave as-is") ; break ;
    }
    P4 ("%s", " when done\n") ;
    P4 ("  final_resymbol: remove zeros from amalgamation: %s\n",
	    BOOLEAN (Common->final_resymbol)) ;

    P4 ("  dmin:  LDL' diagonal threshold: % .5g\n    Entries with abs. value"
	    " less than dmin are replaced with +/- dmin.\n", Common->dmin) ;

    P4 ("  grow0: memory reallocation: % .5g\n", Common->grow0) ;
    P4 ("  grow1: memory reallocation: % .5g\n", Common->grow1) ;
    P4 ("  grow2: memory reallocation: %d\n", Common->grow2) ;

    P4 ("%s", "  nrelax, zrelax:  supernodal amalgamation rule:\n") ;
    P4 ("%s", "    s = # columns in two adjacent supernodes\n") ;
    P4 ("%s", "    z = % of zeros in new supernode if they are merged.\n") ;
    P4 ("%s", "    Two supernodes are merged if") ;
    P4 (" (s <= %d) or (no new zero entries) or\n", Common->nrelax [0]) ;
    P4 ("    (s <= %d and ",   Common->nrelax [1]) ;
    P4 ("z < %.5g%%) or",      Common->zrelax [0] * 100) ;
    P4 (" (s <= %d and ",      Common->nrelax [2]) ;
    P4 ("z < %.5g%%) or",      Common->zrelax [1] * 100) ;
    P4 (" (z < %.5g%%)\n",     Common->zrelax [2] * 100) ;

    /* ---------------------------------------------------------------------- */
    /* check workspace */
    /* ---------------------------------------------------------------------- */

    mark = Common->mark ;
    nrow = Common->nrow ;
    Flag = Common->Flag ;
    Head = Common->Head ;
    if (nrow > 0)
    {
	if (mark < 0 || Flag == NULL || Head == NULL)
	{
	    ERROR ("workspace corrupted (Flag and/or Head missing)") ;
	}
	for (i = 0 ; i < nrow ; i++)
	{
	    if (Flag [i] >= mark)
	    {
		PRINT0 (("Flag [%d] = %d, mark = %ld\n", i, Flag [i], mark)) ;
		ERROR ("workspace corrupted (Flag)") ;
	    }
	}
	for (i = 0 ; i <= nrow ; i++)
	{
	    if (Head [i] != EMPTY)
	    {
		PRINT0 (("Head [%d] = %d,\n", i, Head [i])) ;
		ERROR ("workspace corrupted (Head)") ;
	    }
	}
    }
    wsize = Common->xworksize / sizeof (double) ;
    Xwork = Common->Xwork ;
    if (wsize > 0)
    {
	if (Xwork == NULL)
	{
	    ERROR ("workspace corrupted (Xwork missing)") ;
	}
	for (i = 0 ; i < wsize ; i++)
	{
	    if (Xwork [i] != 0.)
	    {
		PRINT0 (("Xwork [%d] = %g\n", i, Xwork [i])) ;
		ERROR ("workspace corrupted (Xwork)") ;
	    }
	}
    }

    /* workspace and parameters are valid */
    P3 ("%s", "  OK\n") ;
    return (TRUE) ;
}


int cholmod_check_common
(
    cholmod_common *Common
)
{
    return (check_common (0, NULL, Common)) ;
}


int cholmod_print_common
(
    char *name,
    cholmod_common *Common
)
{
    int print = (Common == NULL) ? 3 : (Common->print) ;
    return (check_common (print, name, Common)) ;
}


/* ========================================================================== */
/* === cholmod_check_sparse ================================================= */
/* ========================================================================== */

/* Ensure that a sparse matrix in column-oriented form is valid, and optionally
 * print it.  Returns the number of entries on the diagonal or -1 if error.
 *
 * workspace: Iwork (nrow)
 */


static long check_sparse
(
    int *Wi,
    int print,
    char *name,
    cholmod_sparse *A,
    long *nnzdiag,
    cholmod_common *Common
)
{
    double *Ax ;
    int *Ap, *Ai, *Anz ;
    int nrow, ncol, nzmax, sorted, packed, j, p, pend, i, nz, ilast, values,
	space, count, precise, init_print, dnz ;
    char *type = "sparse" ;

    /* ---------------------------------------------------------------------- */
    /* print header information */
    /* ---------------------------------------------------------------------- */

    P4 ("%s", "\n") ;
    P3 ("%s", "CHOLMOD sparse:  ") ;
    if (name != NULL)
    {
	P3 ("%s: ", name) ;
    }

    if (A == NULL)
    {
	ERROR ("null") ;
    }
    Common->status = CHOLMOD_OK ;

    nrow = A->nrow ;
    ncol = A->ncol ;
    nzmax = A->nzmax ;
    sorted = A->sorted ;
    packed = A->packed ;
    values = (A->xtype != CHOLMOD_PATTERN) ;
    Ap = A->p ;
    Ai = A->i ;
    Ax = A->x ;
    Anz = A->nz ;
    nz = cholmod_nnz (A, Common) ;

    precise = Common->precise ;
    P3 (" %d", nrow) ;
    P3 ("-by-%d, ", ncol) ;
    P3 ("nz %d", nz) ;
    P4 ("\n  nzmax %d, ", nzmax) ;
    if (nz > nzmax)
    {
	ERROR ("nzmax too small") ;
    }
    if (!sorted)
    {
	P4 ("%s", "un") ;
    }
    P4 ("%s", "sorted, ") ;
    if (!packed)
    {
	P4 ("%s", "un") ;
    }
    P4 ("%s", "packed, ") ;
    if (A->stype > 0)
    {
	P4 ("%s", "upper, ") ;
    }
    else if (A->stype < 0)
    {
	P4 ("%s", "lower, ") ;
    }
    else
    {
	P4 ("%s", "up/lo, ") ;
    }

    switch (A->itype)
    {
	case CHOLMOD_INT:     P4 ("%s", "int, ") ; break ;
	case CHOLMOD_INTLONG: ERROR ("int/long type unsupported") ;
	case CHOLMOD_LONG:    ERROR ("long type unsupported") ;
	default:	      ERROR ("unknown itype") ;
    }

    /* FUTURE WORK:
    if (A->itype != Common->itype) ERROR ("itype must match Common->itype") ;
    */

    switch (A->xtype)
    {
	case CHOLMOD_PATTERN: P4 ("%s", "pattern") ;	       break ;
	case CHOLMOD_REAL:    P4 ("%s", "real") ;	       break ;
	case CHOLMOD_COMPLEX: ERROR ("complex unsupported") ;
	case CHOLMOD_ZOMPLEX: ERROR ("zomplex unsupported") ;
	default:	      ERROR ("unknown xtype") ;
    }

    switch (A->dtype)
    {
	case CHOLMOD_DOUBLE:  P4 ("%s", ", double\n") ;	       break ;
	case CHOLMOD_FLOAT:   ERROR ("float unsupported") ;
	default:	      ERROR ("unknown dtype") ;
    }

    if (A->stype && nrow != ncol) ERROR ("symmetric but not square") ;

    /* check for existence of Ap, Ai, Anz, and Ax arrays */
    if (Ap == NULL)
    {
	ERROR ("p array not present") ;
    }
    if (Ai == NULL)
    {
	ERROR ("i array not present") ;
    }
    if (!packed && Anz == NULL)
    {
	ERROR ("nz array not present") ;
    }
    if (values && Ax == NULL)
    {
	ERROR ("x array not present") ;
    }

    /* packed matrices must start at Ap [0] = 0 */
    if (packed && Ap [0] != 0)
    {
	ERROR ("p [0] must be zero") ;
    }
    if (Ap [ncol] < Ap [0] || Ap [ncol] > nzmax)
    {
	ERROR ("p [ncol] invalid") ;
    }

    /* ---------------------------------------------------------------------- */
    /* allocate workspace if needed */
    /* ---------------------------------------------------------------------- */

    if (!sorted)
    {
	if (Wi == NULL)
	{
	    cholmod_allocate_work (0, nrow, 0, 0, Common) ;
	    Wi = Common->Iwork ;	/* size nrow, (i/i/l) */
	}
	if (Common->status < CHOLMOD_OK)
	{
	    return (FALSE) ;	    /* out of memory */
	}
	for (i = 0 ; i < nrow ; i++)
	{
	    Wi [i] = EMPTY ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* check and print each column */
    /* ---------------------------------------------------------------------- */

    init_print = print ;
    dnz = 0 ;
    ETC_START (count, 8) ;

    for (j = 0 ; j < ncol ; j++)
    {
	ETC (j == ncol-1, count, 4) ;
	p = Ap [j] ;
	if (packed)
	{
	    pend = Ap [j+1] ;
	    nz = pend - p ;
	}
	else
	{
	    /* Note that Anz [j] < 0 is treated as zero */
	    nz = MAX (0, Anz [j]) ;
	    pend = p + nz ;
	}
	space = Ap [j+1] - p ;
	P4 ("  col %d:", j) ;
	P4 (" nz %d", nz) ;
	P4 (" start %d", p) ;
	P4 (" end %d", pend) ;
	if (!packed)
	{
	    P4 (" space %d", space) ;
	}
	P4 ("%s", ":\n") ;
	if (p < 0 || pend > nzmax || space < 0)
	{
	    ERROR ("pointer invalid") ;
	}
	if (nz < 0 || nz > nrow || nz > space)
	{
	    ERROR ("nz invalid") ;
	}
	ilast = EMPTY ;

	for ( ; p < pend ; p++)
	{
	    ETC (j == ncol-1 && p >= pend-4, count, -1) ;
	    i = Ai [p] ;
	    P4 ("  %8d", i) ;
	    if (values)
	    {
		P4 ("%s", ":") ;
		PRINTVALUE (Ax [p]) ;
	    }
	    if (i == j)
	    {
		dnz++ ;
	    }
	    if (i < 0 || i >= nrow)
	    {
		ERROR ("row index out of range") ;
	    }
	    if (sorted && i <= ilast)
	    {
		ERROR ("row indices out of order") ;
	    }
	    if (!sorted && Wi [i] == j)
	    {
		ERROR ("duplicate row index") ;
	    }
	    P4 ("%s", "\n") ;
	    ilast = i ;
	    if (!sorted)
	    {
		Wi [i] = j ;
	    }
	}
    }

    /* matrix is valid */
    P4 ("  nnz on diagonal: %d\n", dnz) ;
    P3 ("%s", "  OK\n") ;
    *nnzdiag = dnz ;
    return (TRUE) ;
}


int cholmod_check_sparse
(
    cholmod_sparse *A,
    cholmod_common *Common
)
{
    long nnzdiag ;
    RETURN_IF_NULL_COMMON (FALSE) ;
    return (check_sparse (NULL, 0, NULL, A, &nnzdiag, Common)) ;
}


int cholmod_print_sparse
(
    cholmod_sparse *A,
    char *name,
    cholmod_common *Common
)
{
    long nnzdiag ;
    RETURN_IF_NULL_COMMON (FALSE) ;
    return (check_sparse (NULL, Common->print, name, A, &nnzdiag, Common)) ;
}


/* ========================================================================== */
/* === cholmod_check_dense ================================================== */
/* ========================================================================== */

/* Ensure a dense matrix is valid, and optionally print it. */

static int check_dense
(
    int print,
    char *name,
    cholmod_dense *X,
    cholmod_common *Common
)
{
    double *Xx ;
    int i, j, d, nrow, ncol, nzmax, nz, precise, init_print, count ;
    char *type = "dense" ;

    /* ---------------------------------------------------------------------- */
    /* print header information */
    /* ---------------------------------------------------------------------- */

    P4 ("%s", "\n") ;
    P3 ("%s", "CHOLMOD dense:   ") ;
    if (name != NULL)
    {
	P3 ("%s: ", name) ;
    }

    if (X == NULL)
    {
	ERROR ("null") ;
    }
    Common->status = CHOLMOD_OK ;

    nrow = X->nrow ;
    ncol = X->ncol ;
    nzmax = X->nzmax ;
    d = X->d ;
    Xx = X->x ;

    precise = Common->precise ;

    P3 (" %d", nrow) ;
    P3 ("-by-%d, ", ncol) ;
    P4 ("\n  leading dimension %d, ", d) ;
    P4 ("nzmax %d, ", nzmax) ;
    if (d * ncol > nzmax)
    {
	ERROR ("nzmax too small") ;
    }
    if (d < nrow)
    {
	ERROR ("leading dimension must be >= # of rows") ;
    }
    if (Xx == NULL)
    {
	ERROR ("null") ;
    }

    switch (X->xtype)
    {
	case CHOLMOD_PATTERN: ERROR ("pattern unsupported") ;  break ;
	case CHOLMOD_REAL:    P4 ("%s", "real") ;	       break ;
	case CHOLMOD_COMPLEX: ERROR ("complex unsupported") ;
	case CHOLMOD_ZOMPLEX: ERROR ("zomplex unsupported") ;
	default:	      ERROR ("unknown xtype") ;
    }

    switch (X->dtype)
    {
	case CHOLMOD_DOUBLE:  P4 ("%s", ", double\n") ;	       break ;
	case CHOLMOD_FLOAT:   ERROR ("float unsupported") ;
	default:	      ERROR ("unknown dtype") ;
    }

    /* ---------------------------------------------------------------------- */
    /* check and print each entry */
    /* ---------------------------------------------------------------------- */

    if (print >= 4)
    {
	init_print = print ;
	ETC_START (count, 9) ;
	nz = nrow * ncol ;
	for (j = 0 ; j < ncol ; j++)
	{
	    ETC (j == ncol-1, count, 5) ;
	    P4 ("  col %d:\n", j) ;
	    for (i = 0 ; i < nrow ; i++)
	    {
		ETC (j == ncol-1 && i >= nrow-4, count, -1) ;
		P4 ("  %8d:", i) ;
		PRINTVALUE (Xx [i+j*d]) ;
		P4 ("%s", "\n") ;
	    }
	}
    }

    /* dense  is valid */
    P3 ("%s", "  OK\n") ;
    return (TRUE) ;
}


int cholmod_check_dense
(
    cholmod_dense *X,
    cholmod_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    return (check_dense (0, NULL, X, Common)) ;
}


int cholmod_print_dense
(
    cholmod_dense *X,
    char *name,
    cholmod_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    return (check_dense (Common->print, name, X, Common)) ;
}


/* ========================================================================== */
/* === cholmod_check_subset ================================================= */
/* ========================================================================== */

/* Ensure S (0:len-1) is a subset of 0:n-1.  Duplicates are allowed.  S may be
 * NULL.
 *
 * workspace: none
 */

static int check_subset
(
    void *S_p,
    size_t len,
    size_t n,
    int print,
    char *name,
    cholmod_common *Common
)
{
    int *S ;
    int i, k, count, init_print = print ;
    char *type = "subset" ;

    P4 ("%s", "\n") ;
    P3 ("%s", "CHOLMOD subset:  ") ;
    if (name != NULL)
    {
	P3 ("%s: ", name) ;
    }

    P3 (" len: %d ", len) ;
    P3 ("n: %d", n) ;
    P4 ("%s", "\n") ;

    S = S_p ;
    Common->status = CHOLMOD_OK ;
    if (len == 0 || S == NULL)
    {
	P3 ("%s", "  OK\n") ;
	return (TRUE) ;
    }

    if (print >= 4)
    {
	ETC_START (count, 8) ;
	for (k = 0 ; k < ((int) len) ; k++)
	{
	    ETC (k == ((int) len) - 4, count, -1) ;
	    i = S [k] ;
	    P4 ("  %8d:", k) ;
	    P4 (" %d\n", i) ;
	    if (i < 0 || i >= ((int) n))
	    {
		ERROR ("entry out range") ;
	    }
	}
    }
    else
    {
	for (k = 0 ; k < ((int) len) ; k++)
	{
	    i = S [k] ;
	    if (i < 0 || i >= ((int) n))
	    {
		ERROR ("entry out range") ;
	    }
	}
    }
    P3 ("%s", "  OK\n") ;
    return (TRUE) ;
}


int cholmod_check_subset
(
    void *S,
    size_t len,
    size_t n,
    cholmod_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    return (check_subset (S, len, n, 0, NULL, Common)) ;
}


int cholmod_print_subset
(
    void *S,
    size_t len,
    size_t n,
    char *name,
    cholmod_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    return (check_subset (S, len, n, Common->print, name, Common)) ;
}


/* ========================================================================== */
/* === cholmod_check_perm =================================================== */
/* ========================================================================== */

/* Ensure that Perm [0..len-1] is a permutation of a subset of 0:n-1.  Perm
 * may be NULL, which is interpreted as the identity permutation.  There can
 * be no duplicate entries (len must be <= n).
 *
 * If n <= Common->nrow, then this routine takes O(len) time and does not
 * allocate any memory, by using Common->Flag.  Otherwise, it takes O(n) time
 * and ensures that Common->Iwork is at least n*sizeof(int) in size.
 *
 * To check the fset:	    cholmod_check_perm (fset, fsize, ncol, Common) ;
 * To check a permutation:  cholmod_check_perm (Perm, n, n, Common) ;
 *
 * workspace:  Flag (n) if n <= Common->nrow, Iwork (n) otherwise.
 */

static int check_perm
(
    int *Wi,
    int print,
    char *name,
    void *Perm_p,
    size_t len,
    size_t n,
    cholmod_common *Common
)
{
    int *Flag, *Perm ;
    int i, k, mark, init_print, count ;
    char *type = "perm" ;

    /* ---------------------------------------------------------------------- */
    /* checks that take O(1) time */
    /* ---------------------------------------------------------------------- */

    Common->status = CHOLMOD_OK ;
    Perm = Perm_p ;
    if (Perm == NULL || n == 0)
    {
	/* Perm is valid implicit identity, or empty */
	return (TRUE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* checks that take O(n) time or require memory allocation */
    /* ---------------------------------------------------------------------- */

    init_print = print ;
    ETC_START (count, 8) ;

    if (Wi == NULL && n <= Common->nrow)
    {
	/* use the Common->Flag array if it's big enough */
	mark = cholmod_clear_flag (Common) ;
	Flag = Common->Flag ;
	ASSERT (cholmod_dump_work (TRUE, FALSE, 0, Common)) ;
	if (print >= 4)
	{
	    for (k = 0 ; k < ((int) len) ; k++)
	    {
		ETC (k >= ((int) len) - 4, count, -1) ;
		i = Perm [k] ;
		P4 ("    %d:", k) ;
		P4 ("%d\n", i) ;
		if (i < 0 || i >= ((int) n) || Flag [i] == mark)
		{
		    (void) cholmod_clear_flag (Common) ;
		    ERROR ("invalid permutation") ;
		}
		Flag [i] = mark ;
	    }
	}
	else
	{
	    for (k = 0 ; k < ((int) len) ; k++)
	    {
		i = Perm [k] ;
		if (i < 0 || i >= ((int) n) || Flag [i] == mark)
		{
		    (void) cholmod_clear_flag (Common) ;
		    ERROR ("invalid permutation") ;
		}
		Flag [i] = mark ;
	    }
	}
	cholmod_clear_flag (Common) ;
	ASSERT (cholmod_dump_work (TRUE, FALSE, 0, Common)) ;
    }
    else
    {
	if (Wi == NULL)
	{
	    /* use Common->Iwork instead, but initialize it first */
	    cholmod_allocate_work (0, n, 0, 0, Common) ;
	    Wi = Common->Iwork ;		    /* size n, (i/i/i) is OK */
	}
	if (Common->status < CHOLMOD_OK)
	{
	    return (FALSE) ;	    /* out of memory */
	}
	for (i = 0 ; i < ((int) n) ; i++)
	{
	    Wi [i] = FALSE ;
	}
	if (print >= 4)
	{
	    for (k = 0 ; k < ((int) len) ; k++)
	    {
		ETC (k >= ((int) len) - 4, count, -1) ;
		i = Perm [k] ;
		P4 ("    %d:", k) ;
		P4 ("%d\n", i) ;
		if (i < 0 || i >= ((int) n) || Wi [i])
		{
		    ERROR ("invalid permutation") ;
		}
		Wi [i] = TRUE ;
	    }
	}
	else
	{
	    for (k = 0 ; k < ((int) len) ; k++)
	    {
		i = Perm [k] ;
		if (i < 0 || i >= ((int) n) || Wi [i])
		{
		    ERROR ("invalid permutation") ;
		}
		Wi [i] = TRUE ;
	    }
	}
    }

    /* perm is valid */
    return (TRUE) ;
}


int cholmod_check_perm
(
    void *Perm,
    size_t len,
    size_t n,
    cholmod_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    return (check_perm (NULL, 0, NULL, Perm, len, n, Common)) ;
}


int cholmod_print_perm
(
    void *Perm,
    size_t len,
    size_t n,
    char *name,
    cholmod_common *Common
)
{
    int ok, print ;
    RETURN_IF_NULL_COMMON (FALSE) ;
    print = Common->print ;
    P4 ("%s", "\n") ;
    P3 ("%s", "CHOLMOD perm:    ") ;
    if (name != NULL)
    {
	P3 ("%s: ", name) ;
    }
    P3 (" len: %d", len) ;
    P3 (" n: %d", n) ;
    P4 ("%s", "\n") ;
    ok = check_perm (NULL, print, name, Perm, len, n, Common) ;
    if (ok)
    {
	P3 ("%s", "  OK\n") ;
    }
    return (ok) ;
}


/* ========================================================================== */
/* === cholmod_check_parent ================================================= */
/* ========================================================================== */

/* Ensure that Parent is a valid elimination tree of nodes 0 to n-1.
 * If j is a root of the tree then Parent [j] is EMPTY (-1).
 *
 * NOTE: this check will fail if applied to the component tree (CParent) in
 * cholmod_nested_dissection, unless it has been postordered and renumbered.
 *
 * workspace: none
 */

static int check_parent
(
    void *Parent_p,
    size_t n,
    int print,
    char *name,
    cholmod_common *Common
)
{
    int *Parent ;
    int j, p, count, init_print ;
    char *type = "parent" ;

    init_print = print ;
    Parent = Parent_p ;

    P4 ("%s", "\n") ;
    P3 ("%s", "CHOLMOD parent:  ") ;
    if (name != NULL)
    {
	P3 ("%s: ", name) ;
    }

    P3 (" n: %d", n) ;
    P4 ("%s", "\n") ;

    Common->status = CHOLMOD_OK ;
    if (Parent == NULL)
    {
	ERROR ("null") ;
    }

    /* ---------------------------------------------------------------------- */
    /* checks that take O(n) time */
    /* ---------------------------------------------------------------------- */

    ETC_START (count, 8) ;
    for (j = 0 ; j < ((int) n) ; j++)
    {
	ETC (j == ((int) n) - 4, count, -1) ;
	p = Parent [j] ;
	P4 ("  %8d:", j) ;
	P4 (" %d\n", p) ;
	if (!(p == EMPTY || p > j))
	{
	    ERROR ("invalid") ;
	}
    }
    P3 ("%s", "  OK\n") ;
    return (TRUE) ;
}

int cholmod_check_parent
(
    void *Parent,
    size_t n,
    cholmod_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    return (check_parent (Parent, n, 0, NULL, Common)) ;
}


int cholmod_print_parent
(
    void *Parent,
    size_t n,
    char *name,
    cholmod_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    return (check_parent (Parent, n, Common->print, name, Common)) ;
}



/* ========================================================================== */
/* === cholmod_check_factor ================================================= */
/* ========================================================================== */

static int check_factor
(
    int *Wi,
    int print,
    char *name,
    cholmod_factor *L,
    cholmod_common *Common
)
{
    double *Lx ;
    int *Lp, *Li, *Lnz, *Lnext, *Lprev, *Perm, *ColCount, *Lpi, *Lpx, *Super,
	*Ls ;
    int n, nzmax, packed, j, p, pend, i, nz, minor, ordering, ftype, space,
	count, precise, init_print, dynamic, ilast, lnz, head, tail, jprev,
	jnext, examine_super, nsuper, s, k1, k2, psi, psend, psx, nsrow, nscol,
	ps2, psxend, ssize, xsize, maxcsize, maxesize, nsrow2, jj, ii ;
    char *type = "factor" ;

    /* ---------------------------------------------------------------------- */
    /* print header information */
    /* ---------------------------------------------------------------------- */

    P4 ("%s", "\n") ;
    P3 ("%s", "CHOLMOD factor:  ") ;
    if (name != NULL)
    {
	P3 ("%s: ", name) ;
    }

    if (L == NULL)
    {
	ERROR ("null") ;
    }
    Common->status = CHOLMOD_OK ;

    n = L->n ;
    minor = L->minor ;
    ordering = L->ordering ;
    ftype = L->ftype ;

    Perm = L->Perm ;
    ColCount = L->ColCount ;
    lnz = 0 ;

    precise = Common->precise ;

    P3 (" %d", n) ;
    P3 ("-by-%d", n) ;

    switch (L->itype)
    {
	case CHOLMOD_INT:     P4 ("%s", ", int") ; break ;
	case CHOLMOD_INTLONG: ERROR ("int/long type unsupported") ;
	case CHOLMOD_LONG:    ERROR ("long type unsupported") ;
	default:	      ERROR ("unknown itype") ;
    }

    /* FUTURE WORK:
    if (L->itype != Common->itype) ERROR ("itype must match Common->itype") ;
    */

    switch (L->xtype)
    {
	case CHOLMOD_PATTERN: P4 ("%s", ", pattern") ;	       break ;
	case CHOLMOD_REAL:    P4 ("%s", ", real") ;	       break ;
	case CHOLMOD_COMPLEX: ERROR ("complex unsupported") ;
	case CHOLMOD_ZOMPLEX: ERROR ("zomplex unsupported") ;
	default:	      ERROR ("unknown xtype") ;
    }

    switch (L->dtype)
    {
	case CHOLMOD_DOUBLE:  P4 ("%s", ", double") ;	       break ;
	case CHOLMOD_FLOAT:   ERROR ("float unsupported") ;
	default:	      ERROR ("unknown dtype") ;
    }

    switch (ftype)
    {
	case CHOLMOD_SYMBOLIC_SUPER: P3 ("%s", ", supernodal-symbolic.") ;
				     break ;
	case CHOLMOD_SYMBOLIC:       P3 ("%s", ", symbolic.") ;		 break ;
	case CHOLMOD_LDL_PACKED:     P3 ("%s", ", LDL' packed.") ;	 break ;
	case CHOLMOD_LDL_UNPACKED:   P3 ("%s", ", LDL' unpacked.") ;	 break ;
	case CHOLMOD_LDL_DYNAMIC:    P3 ("%s", ", LDL' dynamic.") ;	 break ;
	case CHOLMOD_LL_PACKED:      P3 ("%s", ", LL' packed.") ;	 break ;
	case CHOLMOD_LL_SUPER:       P3 ("%s", ", LL' supernodal.") ;	 break ;
	default:	             ERROR ("unknown factor type") ;
    }

    P4 ("%s", "\n  ordering method used: ") ;
    switch (L->ordering)
    {
	case CHOLMOD_NATURAL:	P4 ("%s", "natural") ;			 break ;
	case CHOLMOD_GIVEN:	P4 ("%s", "user-provided") ;		 break ;
	case CHOLMOD_AMD:	P4 ("%s", "AMD/COLAMD") ;		 break ;
	case CHOLMOD_METIS:	P4 ("%s", "METIS NodeND") ;		 break ;
	case CHOLMOD_ND:	P4 ("%s", "CHOLMOD nested dissection") ; break ;
	default:		ERROR ("unknown ordering") ;
    }

    P4 ("%s", "\n") ;

    if (ftype <= CHOLMOD_SYMBOLIC)
    {
	if (L->xtype != CHOLMOD_PATTERN)
	{
	    ERROR ("symbolic L cannot have numerical values") ;
	}
    }
    else
    {
	if (L->xtype == CHOLMOD_PATTERN)
	{
	    ERROR ("numeric L must have numerical values") ;
	}
    }

    init_print = print ;

    /* ---------------------------------------------------------------------- */
    /* check L->Perm */
    /* ---------------------------------------------------------------------- */

    if (!check_perm (Wi, print, name, Perm, n, n, Common))
    {
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* check L->ColCount */
    /* ---------------------------------------------------------------------- */

    if (ColCount == NULL)
    {
	ERROR ("ColCount vector invalid") ;
    }

    ETC_START (count, 8) ;
    for (j = 0 ; j < n ; j++)
    {
	ETC (j >= n-4, count, -1) ;
	P4 ("  col: %d ", j) ;
	nz = ColCount [j] ;
	P4 ("colcount: %d\n", nz) ;
	if (nz < 0 || nz > n-j)
	{
	    ERROR ("ColCount out of range") ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* check factor */
    /* ---------------------------------------------------------------------- */

    if (ftype == CHOLMOD_SYMBOLIC)
    {

	/* ------------------------------------------------------------------ */
	/* check simplicial symbolic factor */
	/* ------------------------------------------------------------------ */

	/* nothing else to do */ ;

    }
    else if (ftype >= CHOLMOD_LDL_PACKED && ftype <= CHOLMOD_LL_PACKED)
    {

	/* ------------------------------------------------------------------ */
	/* check simplicial numerical factor */
	/* ------------------------------------------------------------------ */

	nzmax = L->nzmax ;
	Lp = L->p ;
	Li = L->i ;
	Lx = L->x ;
	Lnz = L->nz ;
	Lnext = L->next ;
	Lprev = L->prev ;

	packed = (ftype == CHOLMOD_LDL_PACKED || ftype == CHOLMOD_LL_PACKED) ;
	dynamic = (ftype == CHOLMOD_LDL_DYNAMIC) ;

	/* check for existence of Lp, Li, Lnz, and Lx arrays */
	if (Lp == NULL)
	{
	    ERROR ("p array not present") ;
	}
	if (Li == NULL)
	{
	    ERROR ("i array not present") ;
	}
	if (!packed && Lnz == NULL)
	{
	    ERROR ("nz array not present") ;
	}
	if (Lx == NULL)
	{
	    ERROR ("x array not present") ;
	}

	ETC_START (count, 8) ;

	/* check each column of L */
	for (j = 0 ; j < n ; j++)
	{
	    ETC (j >= n-3, count, -1) ;
	    p = Lp [j] ;
	    if (packed)
	    {
		pend = Lp [j+1] ;
		nz = pend - p ;
	    }
	    else
	    {
		nz = Lnz [j] ;
		pend = p + nz ;
	    }
	    lnz += nz ;

	    P4 ("  col %d:", j) ;
	    P4 (" nz %d", nz) ;
	    P4 (" start %d", p) ;
	    P4 (" end %d", pend) ;

	    if (dynamic)
	    {
		if (Lnext [j] < 0 || Lnext [j] > n)
		{
		    ERROR ("invalid link list")  ;
		}
		space = Lp [Lnext [j]] - p ;
	    }
	    else if (j < n-1)
	    {
		space = Lp [j+1] - p ;
	    }
	    else
	    {
		space = nzmax - p ;
	    }

	    if (!packed)
	    {
		P4 (" space %d", space) ;
		P4 (" free %d", space - nz) ;
	    }
	    P4 ("%s", ":\n") ;
	    if (p < 0 || pend > nzmax || space < 1)
	    {
		ERROR ("pointer invalid") ;
	    }
	    if (nz < 1 || nz > (n-j) || nz > space)
	    {
		ERROR ("nz invalid") ;
	    }
	    ilast = j-1 ;

	    i = Li [p] ;
	    P4 ("  %8d:", i) ;
	    if (i != j)
	    {
		ERROR ("diagonal missing") ;
	    }
	    PRINTVALUE (Lx [p]) ;
	    P4 ("%s", "\n") ;
	    ilast = j ;
	    for (p++ ; p < pend ; p++)
	    {
		ETC_DISABLE (count) ;
		i = Li [p] ;
		P4 ("  %8d:", i) ;
		if (i < j || i >= n)
		{
		    ERROR ("row index out of range") ;
		}
		if (i <= ilast)
		{
		    ERROR ("row indices out of order") ;
		}
		PRINTVALUE (Lx [p]) ;
		P4 ("%s", "\n") ;
		ilast = i ;
	    }
	}

	/* check the link list */
	if (dynamic)
	{
	    head = n+1 ;
	    tail = n ;
	    j = head ;
	    jprev = EMPTY ;
	    count = 0 ;
	    for ( ; ; )
	    {
		if (j < 0 || j > n+1 || count > n+2)
		{
		    ERROR ("invalid link list") ;
		}
		jnext = Lnext [j] ;
		if (j >= 0 && j < n)
		{
		    if (jprev != Lprev [j])
		    {
			ERROR ("invalid link list") ;
		    }
		}
		count++ ;
		if (j == tail)
		{
		    break ;
		}
		jprev = j ;
		j = jnext ;
	    }
	    if (Lnext [tail] != EMPTY || count != n+2)
	    {
		ERROR ("invalid link list") ;
	    }
	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* check supernodal numeric or symbolic factor */
	/* ------------------------------------------------------------------ */

	nsuper = L->nsuper ;
	ssize = L->ssize ;
	xsize = L->xsize ;
	maxcsize = L->maxcsize ;
	maxesize = L->maxesize ;
	Ls = L->s ;
	Lpi = L->pi ;
	Lpx = L->px ;
	Super = L->super ;
	Lx = L->x ;
	ETC_START (count, 8) ;

	P4 ("  ssize %d ", ssize) ;
	P4 ("xsize %d ", xsize) ;
	P4 ("maxcsize %d ", maxcsize) ;
	P4 ("maxesize %d\n", maxesize) ;

	if (Ls == NULL)
	{
	    ERROR ("invalid: L->s missing") ;
	}
	if (Lpi == NULL)
	{
	    ERROR ("invalid: L->pi missing") ;
	}
	if (Lpx == NULL)
	{
	    ERROR ("invalid: L->px missing") ;
	}
	if (Super == NULL)
	{
	    ERROR ("invalid: L->super missing") ;
	}

	if (ftype == CHOLMOD_LL_SUPER)
	{
	    /* numerical supernodal factor */
	    if (Lx == NULL)
	    {
		ERROR ("invalid: L->x missing") ;
	    }
	    if (Ls [0] == EMPTY)
	    {
		ERROR ("invalid: L->s not defined") ;
	    }
	    examine_super = TRUE ;
	}
	else
	{
	    /* symbolic supernodal factor, but only if it has been computed */
	    examine_super = (Ls [0] != EMPTY) ;
	}

	if (examine_super)
	{
	    if (Lpi [0] != 0 || Lpi [nsuper] != ssize)
	    {
		printf ("Lpi [0] %d, Lpi [nsuper = %d] = %d\n",
			Lpi [0], nsuper, Lpi [nsuper]) ;
		ERROR ("invalid: L->pi invalid") ;
	    }
	    if (Lpx [0] != 0 || Lpx [nsuper] != xsize)
	    {
		ERROR ("invalid: L->px invalid") ;
	    }

	    /* check and print each supernode */
	    for (s = 0 ; s < nsuper ; s++)
	    {
		k1 = Super [s] ;
		k2 = Super [s+1] ;
		psi = Lpi [s] ;
		psend = Lpi [s+1] ;
		psx = Lpx [s] ;
		nsrow = psend - psi ;
		nscol = k2 - k1 ;
		nsrow2 = nsrow - nscol ;
		ps2 = psi + nscol ;
		psxend = Lpx [s+1] ;

		ETC (s == nsuper-1, count, 4) ;

		P4 ("  supernode %d, ", s) ;
		P4 ("col %d ", k1) ;
		P4 ("to %d. ", k2-1) ;
		P4 ("nz in first col: %d.\n", nsrow) ;
		P4 ("  values start %d, ", psx) ;
		P4 ("end %d\n", psxend) ;

		if (k1 > k2 || k1 < 0 || k2 > n || nsrow < nscol || nsrow2 < 0
		    || psxend - psx != nsrow * nscol)
		{
		    ERROR ("invalid supernode") ;
		}

		lnz += nscol * nsrow - (nscol*nscol - nscol)/2 ;

		if (ftype == CHOLMOD_LL_SUPER)
		{
		    /* print each column of the supernode */
		    for (jj = 0 ; jj < nscol ; jj++)
		    {
			ETC_ENABLE (s == nsuper-1 && jj >= nscol-3, count, -1) ;
			j = k1 + jj ;
			P4 ("  col %d\n", j) ;
			ilast = j ;
			i = Ls [psi + jj] ;
			P4 ("  %8d:", i) ;
			if (i != j)
			{
			    ERROR ("row index invalid") ;
			}
			PRINTVALUE (Lx [psx + jj + jj*nsrow]) ;
			P4 ("%s", "\n") ;
			for (ii = jj + 1 ; ii < nsrow ; ii++)
			{
			    ETC_DISABLE (count) ;
			    i = Ls [psi + ii] ;
			    P4 ("  %8d:", i) ;
			    if (i <= ilast || i > n)
			    {
				ERROR ("row index out of range") ;
			    }
			    PRINTVALUE (Lx [psx + ii + jj*nsrow]) ;
			    P4 ("%s", "\n") ;
			    ilast = i ;
			}
		    }
		}
		else
		{
		    /* just print the leading column of the supernode */
		    P4 ("  col %d\n", k1) ;
		    for (jj = 0 ; jj < nscol ; jj++)
		    {
			ETC (s == nsuper-1 && jj >= nscol-3, count, -1) ;
			j = k1 + jj ;
			i = Ls [psi + jj] ;
			P4 ("  %8d", i) ;
			if (i != j)
			{
			    ERROR ("row index invalid") ;
			}
			P4 ("%s", "\n") ;
		    }
		    ilast = j ;
		    for (ii = nscol ; ii < nsrow ; ii++)
		    {
			ETC_DISABLE (count) ;
			i = Ls [psi + ii] ;
			P4 ("  %8d", i) ;
			if (i <= ilast || i > n)
			{
			    ERROR ("row index out of range") ;
			}
			P4 ("%s", "\n") ;
			ilast = i ;
		    }
		}
	    }
	}
    }

    /* factor is valid */
    P3 ("  nz %d", lnz) ;
    P3 ("%s", "  OK\n") ;
    return (TRUE) ;
}


int cholmod_check_factor
(
    cholmod_factor *L,
    cholmod_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    return (check_factor (NULL, 0, NULL, L, Common)) ;
}


int cholmod_print_factor
(
    cholmod_factor *L,
    char *name,
    cholmod_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    return (check_factor (NULL, Common->print, name, L, Common)) ;
}


/* ========================================================================== */
/* === cholmod_check_triplet ================================================ */
/* ========================================================================== */

/* Ensure a triplet matrix is valid, and optionally print it. */

static int check_triplet
(
    int print,
    char *name,
    cholmod_triplet *T,
    cholmod_common *Common
)
{
    double *Tx ;
    int *Ti, *Tj ;
    int i, j, p, nrow, ncol, nzmax, nz, precise, values, init_print, count ;
    char *type = "triplet" ;

    /* ---------------------------------------------------------------------- */
    /* print header information */
    /* ---------------------------------------------------------------------- */

    P4 ("%s", "\n") ;
    P3 ("%s", "CHOLMOD triplet: ") ;
    if (name != NULL)
    {
	P3 ("%s: ", name) ;
    }

    if (T == NULL)
    {
	ERROR ("null") ;
    }
    Common->status = CHOLMOD_OK ;

    nrow = T->nrow ;
    ncol = T->ncol ;
    nzmax = T->nzmax ;
    nz = T->nnz ;
    Ti = T->i ;
    Tj = T->j ;
    Tx = T->x ;
    values = (T->xtype != CHOLMOD_PATTERN) ;

    precise = Common->precise ;

    P3 (" %d", nrow) ;
    P3 ("-by-%d, ", ncol) ;
    P3 ("nz %d", nz) ;
    P4 ("\n  nzmax %d, ", nzmax) ;
    if (nz > nzmax)
    {
	ERROR ("nzmax too small") ;
    }
    if (T->stype > 0)
    {
	P4 ("%s", "upper, ") ;
    }
    else if (T->stype < 0)
    {
	P4 ("%s", "lower, ") ;
    }
    else
    {
	P4 ("%s", "up/lo, ") ;
    }

    switch (T->itype)
    {
	case CHOLMOD_INT:     P4 ("%s", "int, ") ; break ;
	case CHOLMOD_INTLONG: ERROR ("int/long type unsupported") ;
	case CHOLMOD_LONG:    ERROR ("long type unsupported") ;
	default:	      ERROR ("unknown itype") ;
    }

    /* FUTURE WORK:
    if (T->itype != Common->itype) ERROR ("itype must match Common->itype") ;
    */

    switch (T->xtype)
    {
	case CHOLMOD_PATTERN: P4 ("%s", "pattern") ;	       break ;
	case CHOLMOD_REAL:    P4 ("%s", "real") ;	       break ;
	case CHOLMOD_COMPLEX: ERROR ("complex unsupported") ;
	case CHOLMOD_ZOMPLEX: ERROR ("zomplex unsupported") ;
	default:	      ERROR ("unknown xtype") ;
    }

    switch (T->dtype)
    {
	case CHOLMOD_DOUBLE:  P4 ("%s", ", double\n") ;	       break ;
	case CHOLMOD_FLOAT:   ERROR ("float unsupported") ;
	default:	      ERROR ("unknown dtype") ;
    }

    if (T->stype && nrow != ncol) ERROR ("symmetric but not square") ;

    /* check for existence of Ti, Tj, Tx arrays */
    if (Tj == NULL)
    {
	ERROR ("j array not present") ;
    }
    if (Ti == NULL)
    {
	ERROR ("i array not present") ;
    }
    if (values && Tx == NULL)
    {
	ERROR ("x array not present") ;
    }

    /* ---------------------------------------------------------------------- */
    /* check and print each entry */
    /* ---------------------------------------------------------------------- */

    init_print = print ;
    ETC_START (count, 8) ;

    for (p = 0 ; p < nz ; p++)
    {
	ETC (p >= nz-4, count, -1) ;
	i = Ti [p] ;
	P4 ("  %8d:", p) ;
	P4 (" %-8d", i) ;
	if (i < 0 || i > nrow)
	{
	    ERROR ("row index out of range") ;
	}
	j = Tj [p] ;
	P4 (" %-8d", j) ;
	if (j < 0 || j > ncol)
	{
	    ERROR ("row index out of range") ;
	}
	if (values)
	{
	    PRINTVALUE (Tx [p]) ;
	}
	P4 ("%s", "\n") ;
    }

    /* triplet matrix is valid */
    P3 ("%s", "  OK\n") ;
    return (TRUE) ;
}


int cholmod_check_triplet
(
    cholmod_triplet *T,
    cholmod_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    return (check_triplet (0, NULL, T, Common)) ;
}


int cholmod_print_triplet
(
    cholmod_triplet *T,
    char *name,
    cholmod_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    return (check_triplet (Common->print, name, T, Common)) ;
}



/* ========================================================================== */
/* === CHOLMOD debugging routines =========================================== */
/* ========================================================================== */

#ifndef NDEBUG

/* The global variable cholmod_dump is present only when debugging enabled. */
int cholmod_dump = 0 ;
int cholmod_dump_malloc = 0 ;
int *cholmod_work = NULL ;
size_t cholmod_work_size = 0 ;

/* workspace: no debug routines use workspace in Common */

/* ========================================================================== */
/* === cholmod_dump_init ==================================================== */
/* ========================================================================== */

void cholmod_dump_init (char *s)
{
    FILE *f ;
    f = fopen ("debug", "r") ;
    if (f == NULL)
    {
	cholmod_dump = 0 ;
    }
    else
    {
	fscanf (f, "%d", &cholmod_dump) ;
	fclose (f) ;
    }
    PRINT1 (("%s: cholmod_dump_init, D = %d\n", s, cholmod_dump)) ;
}


/* ========================================================================== */
/* === cholmod_dump_getwork ================================================= */
/* ========================================================================== */

/* get workspace for dump routines */

void cholmod_dump_getwork
(
    size_t n
)
{
    n = MAX (n,1) ;
    if (n > cholmod_work_size)
    {
	cholmod_dump_freework ( ) ;
#ifdef MATLAB_MEX_FILE
	cholmod_work = mxMalloc (n * sizeof (int)) ;
#else
	cholmod_work = malloc (n * sizeof (int)) ;
#endif
	cholmod_work_size = n ;
	ASSERT (cholmod_work) ;
    }
}


/* ========================================================================== */
/* === cholmod_dump_freework ================================================ */
/* ========================================================================== */

/* free workspace for dump routines */

void cholmod_dump_freework
(
    void
)
{
    if (cholmod_work != NULL)
    {
#ifdef MATLAB_MEX_FILE
	mxFree (cholmod_work) ;
#else
	free (cholmod_work) ;
#endif
	cholmod_work_size = 0 ;
	cholmod_work = NULL ;
    }
}


/* ========================================================================== */
/* === cholmod_dump_sparse ================================================== */
/* ========================================================================== */

long cholmod_dump_sparse	/* returns nnz (diag (A)) or EMPTY if error */
(
    cholmod_sparse *A,
    char *name,
    cholmod_common *Common
)
{
    long nnzdiag ;
    int ok ;

    if (cholmod_dump < -1)
    {
	/* no checks if debug level is -2 or less */
	return (0) ;
    }

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    cholmod_dump_getwork (A->nrow) ;
    ok = check_sparse (cholmod_work, cholmod_dump, name, A, &nnzdiag, Common) ;
    return (ok ? nnzdiag : EMPTY) ;
}


/* ========================================================================== */
/* === cholmod_dump_factor ================================================== */
/* ========================================================================== */

int cholmod_dump_factor
(
    cholmod_factor *L,
    char *name,
    cholmod_common *Common
)
{

    if (cholmod_dump < -1)
    {
	/* no checks if debug level is -2 or less */
	return (TRUE) ;
    }
    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    cholmod_dump_getwork (L->n) ;
    return (check_factor (cholmod_work, cholmod_dump, name, L, Common)) ;
}


/* ========================================================================== */
/* === cholmod_dump_perm ==================================================== */
/* ========================================================================== */

int cholmod_dump_perm
(
    void *Perm,
    size_t len,
    size_t n,
    char *name,
    cholmod_common *Common
)
{

    if (cholmod_dump < -1)
    {
	/* no checks if debug level is -2 or less */
	return (TRUE) ;
    }
    RETURN_IF_NULL_COMMON (FALSE) ;
    cholmod_dump_getwork (n) ;
    return (check_perm (cholmod_work, cholmod_dump, name, Perm, len, n,Common));
}


/* ========================================================================== */
/* === cholmod_dump_dense =================================================== */
/* ========================================================================== */

int cholmod_dump_dense
(
    cholmod_dense *X,
    char *name,
    cholmod_common *Common
)
{
    if (cholmod_dump < -1)
    {
	/* no checks if debug level is -2 or less */
	return (TRUE) ;
    }
    RETURN_IF_NULL_COMMON (FALSE) ;
    return (check_dense (cholmod_dump, name, X, Common)) ;
}


/* ========================================================================== */
/* === cholmod_dump_triplet ================================================= */
/* ========================================================================== */

int cholmod_dump_triplet
(
    cholmod_triplet *T,
    char *name,
    cholmod_common *Common
)
{
    if (cholmod_dump < -1)
    {
	/* no checks if debug level is -2 or less */
	return (TRUE) ;
    }
    RETURN_IF_NULL_COMMON (FALSE) ;
    return (check_triplet (cholmod_dump, name, T, Common)) ;
}


/* ========================================================================== */
/* === cholmod_dump_subset ================================================== */
/* ========================================================================== */

int cholmod_dump_subset
(
    void *S,
    size_t len,
    size_t n,
    char *name,
    cholmod_common *Common
)
{
    if (cholmod_dump < -1)
    {
	/* no checks if debug level is -2 or less */
	return (TRUE) ;
    }
    RETURN_IF_NULL_COMMON (FALSE) ;
    return (check_subset (S, len, n, cholmod_dump, name, Common)) ;
}


/* ========================================================================== */
/* === cholmod_dump_parent ================================================== */
/* ========================================================================== */

int cholmod_dump_parent
(
    void *Parent,
    size_t n,
    char *name,
    cholmod_common *Common
)
{
    if (cholmod_dump < -1)
    {
	/* no checks if debug level is -2 or less */
	return (TRUE) ;
    }
    RETURN_IF_NULL_COMMON (FALSE) ;
    return (check_parent (Parent, n, cholmod_dump, name, Common)) ;
}


/* ========================================================================== */
/* === cholmod_dump_real ==================================================== */
/* ========================================================================== */

void cholmod_dump_real
(
    char *name, double *X, int nrow, int ncol
)
{
    /* dump an nrow-by-ncol real dense matrix */
    int i, j ;
    double value ;
    if (cholmod_dump < -1)
    {
	/* no checks if debug level is -2 or less */
	return ;
    }
    PRINT1 (("%s: dump_real, nrow: %d ncol: %d\n", name, nrow, ncol)) ;
    for (j = 0 ; j < ncol ; j++)
    {
	PRINT2 (("    col %d\n", j)) ;
	for (i = 0 ; i < nrow ; i++)
	{
	    /* X is stored in column-major form */
	    value = *X++ ;
	    PRINT2 (("        %5d: %e\n", i, value)) ;
	}
    }
}


/* ========================================================================== */
/* === cholmod_dump_super =================================================== */
/* ========================================================================== */

void cholmod_dump_super
(
    int s, int Super [ ], int Lpi [ ], int Ls [ ], int Lpx [ ], double Lx [ ]
)
{
    int k1, k2, do_values, psi, psx, nsrow, nscol, psend, ilast, p, i ;
    if (cholmod_dump < -1)
    {
	/* no checks if debug level is -2 or less */
	return ;
    }
    k1 = Super [s] ;
    k2 = Super [s+1] ;
    nscol = k2 - k1 ;
    do_values = (Lpx != NULL) && (Lx != NULL) ;
    psi = Lpi [s] ;
    psend = Lpi [s+1] ;
    nsrow = psend - psi ;
    PRINT1 (("\nSuper %d, columns %d to %d, %d rows %d cols\n",
		s, k1, k2-1, nsrow, nscol)) ;
    ilast = -1 ;
    for (p = psi ; p < psend ; p++)
    {
	i = Ls [p] ;
	PRINT2 (("  %d : p-psi %d\n", i, p-psi)) ;
	if (p-psi < nscol) ASSERT (i == k1 + (p-psi)) ;
	if (p-psi == nscol-1) PRINT2 (("------\n")) ;
	ASSERT (i > ilast) ;
	ilast = i ;
    }
    if (do_values)
    {
	psx = Lpx [s] ;
	cholmod_dump_real ("Supernode", Lx + psx, nsrow, nscol) ;
    }
}


/* ========================================================================== */
/* === cholmod_dump_mdelta ================================================== */
/* ========================================================================== */

/* Defines the number of blocks of memory allocated for each type of factor
 * object. */

static int mem_count [7] = {		/* Each has Perm and ColCount, plus: */
/* CHOLMOD_SYMBOLIC_SUPER */    6,	/* super, pi, px, s		*/
/* CHOLMOD_SYMBOLIC	  */    2,	/*				*/
/* CHOLMOD_LDL_PACKED	  */    5,	/* p, i, x			*/
/* CHOLMOD_LDL_UNPACKED	  */    6,	/* p, i, x, nz			*/
/* CHOLMOD_LDL_DYNAMIC	  */    8,	/* p, i, x, nz, next, prev	*/
/* CHOLMOD_LL_PACKED	  */    5,	/* p, i, x			*/
/* CHOLMOD_LL_SUPER	  */    7 } ;	/* super, pi, px, s, x		*/

/* return the change in malloc_count if a factor of type fin is converted to
 * a factor of type fout */
int cholmod_dump_mdelta (int fin, int fout)
{
    ASSERT (fin  >= CHOLMOD_SYMBOLIC_SUPER && fin  <= CHOLMOD_LL_SUPER) ;
    ASSERT (fout >= CHOLMOD_SYMBOLIC_SUPER && fout <= CHOLMOD_LL_SUPER) ;
    return (mem_count [fout+2] - mem_count [fin+2]) ;
}


/* ========================================================================== */
/* === cholmod_dump_mem ===================================================== */
/* ========================================================================== */

int cholmod_dump_mem (char *where, int should, cholmod_common *Common)
{
    int diff = should - Common->memory_inuse ;
    if (diff != 0)
    {
	printf ("mem: %-15s peak %10d inuse %10d should %10d\n",
	    where, Common->memory_usage, Common->memory_inuse, should) ;
	printf ("mem: %s diff %d !\n", where, diff) ;
    }
    return (diff == 0) ;
}


/* ========================================================================== */
/* === cholmod_dump_partition =============================================== */
/* ========================================================================== */

/* make sure we have a proper separator (for debugging only)
 *
 * workspace: none
 */

int cholmod_dump_partition
(
    int n,
    int Cp [ ],
    int Ci [ ],
    int Cnw [ ],
    int Part [ ],
    int sepsize
)
{
    int chek [3], which, ok, i, j, p ;
    PRINT1 (("bisect sepsize %d\n", sepsize)) ;
    ok = TRUE ;
    chek [0] = 0 ;
    chek [1] = 0 ;
    chek [2] = 0 ;
    for (j = 0 ; j < n ; j++)
    {
	PRINT2 (("--------------j %d in part %d nw %d\n", j, Part [j], Cnw[j]));
	which = Part [j] ;
	for (p = Cp [j] ; p < Cp [j+1] ; p++)
	{
	    i = Ci [p] ;
	    PRINT3 (("i %d, part %d\n", i, Part [i])) ;
	    if (which == 0)
	    {
		if (Part [i] == 1)
		{
		    PRINT0 (("Error! %d %d\n", i, j)) ;
		    ok = FALSE ;
		}
	    }
	    else if (which == 1)
	    {
		if (Part [i] == 0)
		{
		    PRINT0 (("Error! %d %d\n", i, j)) ;
		    ok = FALSE ;
		}
	    }
	}
	if (which < 0 || which > 2)
	{
	    PRINT0 (("Part out of range\n")) ;
	    ok = FALSE ;
	}
	chek [which] += Cnw [j] ;
    }
    PRINT1 (("sepsize %d check %d %d %d\n", sepsize, chek[0], chek[1],chek[2]));
    if (sepsize != chek[2])
    {
	PRINT0 (("mismatch!\n")) ;
	ok = FALSE ;
    }
    return (ok) ;
}


/* ========================================================================== */
/* === cholmod_dump_work ==================================================== */
/* ========================================================================== */

int cholmod_dump_work (int flag, int head, int wsize, cholmod_common *Common)
{
    double *W ;
    int *Flag, *Head ;
    int k, nrow, mark ;

    if (cholmod_dump < -1)
    {
	/* no checks if debug level is -2 or less */
	return (TRUE) ;
    }

    RETURN_IF_NULL_COMMON (FALSE) ;
    nrow = Common->nrow ;
    Flag = Common->Flag ;
    Head = Common->Head ;
    W = Common->Xwork ;
    mark = Common->mark ;
    
    if (wsize < 0)
    {
	/* check all of Xwork */
	wsize = Common->xworksize / sizeof (double) ;
    }
    else
    {
	/* check on the first wsize doubles in Xwork */
	wsize = MIN (wsize, (int) (Common->xworksize / sizeof (double))) ;
    }

    if (flag)
    {
	for (k = 0 ; k < nrow ; k++)
	{
	    if (Flag [k] >= mark)
	    {
		PRINT0 (("Flag invalid, Flag [%d] = %d, mark = %d\n",
			    k, Flag [k], mark)) ;
		ASSERT (0) ;
		return (FALSE) ;
	    }
	}
    }

    if (head)
    {
	for (k = 0 ; k < nrow ; k++)
	{
	    if (Head [k] != EMPTY)
	    {
		PRINT0 (("Head invalid, Head [%d] = %d\n", k, Head [k])) ;
		ASSERT (0) ;
		return (FALSE) ;
	    }
	}
    }

    for (k = 0 ; k < wsize ; k++)
    {
	if (W [k] != 0.)
	{
	    PRINT0 (("W invalid, W [%d] = %g\n", k, W [k])) ;
	    ASSERT (0) ;
	    return (FALSE) ;
	}
    }

    return (TRUE) ;
}
#endif
