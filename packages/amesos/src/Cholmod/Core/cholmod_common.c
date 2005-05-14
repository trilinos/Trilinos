/* ========================================================================== */
/* === Core/cholmod_common ================================================== */
/* ========================================================================== */

/*
 * CHOLMOD/Core version 0.1. May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* Core utility routines for the cholmod_common object:
 *
 * Primary routines:
 * -----------------
 * cholmod_start		the first call to CHOLMOD
 * cholmod_finish		the last call to CHOLMOD
 *
 * Secondary routines:
 * -------------------
 * cholmod_defaults		restore (most) default control parameters
 * cholmod_allocate_work	allocate (or reallocate) workspace in Common
 * cholmod_free_work		free workspace in Common
 * cholmod_clear_flag		clear Common->Flag in workspace
 * cholmod_maxrank		column dimension of Common->Xwork workspace
 *
 * The Common object is unique.  It cannot be allocated or deallocated by
 * CHOLMOD, since it contains the definition of the memory management routines
 * used (pointers to malloc, free, realloc, and calloc, or their equivalent).
 * The Common object contains workspace that is used between calls to
 * CHOLMOD routines.  This workspace allocated by CHOLMOD as needed, by
 * cholmod_allocate_work and cholmod_free_work.
 */

#include "cholmod_core.h"
#include "cholmod_internal.h"

/* ========================================================================== */
/* === cholmod_start ======================================================== */
/* ========================================================================== */

/* Initialize Common default parameters and statistics.  Sets workspace
 * pointers to NULL.
 *
 * This routine must be called just once, prior to calling any other CHOLMOD
 * routine.  Do not call this routine after any other CHOLMOD routine (except
 * cholmod_finish, to start a new CHOLMOD session), or a memory leak will
 * occur.
 *
 * workspace: none
 */

int cholmod_start
(
    int itype,
    int xtype,
    int dtype,
    cholmod_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    DEBUG_INIT ("cholmod start") ;

    /* ---------------------------------------------------------------------- */
    /* integer and numerical types */
    /* ---------------------------------------------------------------------- */

    if (itype != CHOLMOD_INT)
    {
	/* only CHOLMOD_INT is currently supported */
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_start: only 'int' type supported", Common) ;
	return (FALSE) ;
    }
    if (xtype != CHOLMOD_REAL)
    {
	/* only CHOLMOD_REAL is currently supported */
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_start: complex type not supported", Common) ;
	return (FALSE) ;
    }
    if (dtype != CHOLMOD_DOUBLE)
    {
	/* only CHOLMOD_DOUBLE is currently supported */
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_start: only 'double' type supported", Common) ;
	return (FALSE) ;
    }

    Common->itype = itype ;
    Common->xtype = xtype ;
    Common->dtype = dtype ;

    /* ---------------------------------------------------------------------- */
    /* default control parameters */
    /* ---------------------------------------------------------------------- */

    cholmod_defaults (Common) ;
    Common->file = NULL ;
    Common->try_catch = FALSE ;

    /* ---------------------------------------------------------------------- */
    /* user error handling routine */
    /* ---------------------------------------------------------------------- */

    Common->error_handler = NULL ;

    /* ---------------------------------------------------------------------- */
    /* memory management routines */
    /* ---------------------------------------------------------------------- */

    /* The user can replace cholmod's memory management routines by redefining
     * these function pointers. */

#ifdef MATLAB_MEX_FILE
    /* MATLAB mexFunction */
    Common->malloc_memory  = mxMalloc ;
    Common->free_memory    = mxFree ;
    Common->realloc_memory = mxRealloc ;
    Common->calloc_memory  = mxCalloc ;
#else
    /* stand-alone ANSI C program */
    Common->malloc_memory  = malloc ;
    Common->free_memory    = free ;
    Common->realloc_memory = realloc ;
    Common->calloc_memory  = calloc ;
#endif

    /* ---------------------------------------------------------------------- */
    /* workspace */
    /* ---------------------------------------------------------------------- */

    /* This code assumes the workspace held in Common is not initialized.  If
     * it is, then a memory leak will occur because the pointers are
     * overwritten with NULL. */

    Common->nrow = 0 ;
    Common->mark = EMPTY ;
    Common->xworksize = 0 ;
    Common->iworksize = 0 ;
    Common->Flag = NULL ;
    Common->Head = NULL ;
    Common->Iwork = NULL ;
    Common->Xwork = NULL ;

    /* ---------------------------------------------------------------------- */
    /* statistics */
    /* ---------------------------------------------------------------------- */

    /* fl and lnz are computed in cholmod_analyze and cholmod_rowcolcounts */
    Common->fl = EMPTY ;
    Common->lnz = EMPTY ;

    /* modfl is computed in cholmod_updown, cholmod_rowadd, and cholmod_rowdel*/
    Common->modfl = EMPTY ;

    /* all routines use status as their error-report code */
    Common->status = CHOLMOD_OK ;

    Common->malloc_count = 0 ;	/* # calls to malloc minus # calls to free */
    Common->memory_usage = 0 ;	/* peak memory usage (in bytes) */
    Common->memory_inuse = 0 ;	/* current memory in use (in bytes) */

    DEBUG (cholmod_dump_getwork (0)) ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_defaults ===================================================== */
/* ========================================================================== */

/* Set Common default parameters, except for the memory management function
 * pointers and the Common->file pointer.
 *
 * workspace: none
 */

int cholmod_defaults
(
    cholmod_common *Common
)
{
    int i ;

    RETURN_IF_NULL_COMMON (FALSE) ;

    /* ---------------------------------------------------------------------- */
    /* default control parameters */
    /* ---------------------------------------------------------------------- */

    Common->dmin = 0.0 ;
    Common->grow0 = 1.2 ;
    Common->grow1 = 1.2 ;
    Common->grow2 = 5 ;
    Common->maxrank = 2 ;
    Common->final_ftype = CHOLMOD_AS_IS ;
    Common->final_resymbol = FALSE ;
    Common->supernodal = TRUE ;
    Common->nrelax [0] = 4 ;
    Common->nrelax [1] = 16 ;
    Common->nrelax [2] = 48 ;
    Common->zrelax [0] = 0.8 ;
    Common->zrelax [1] = 0.1 ;
    Common->zrelax [2] = 0.05 ;
    Common->metis_memory = 2.0 ;
    Common->metis_nswitch = 3000 ;
    Common->metis_dswitch = 0.66 ;
    Common->print = 3 ;
    Common->precise = FALSE ;
    Common->ieee = TRUE ;
    Common->blas_conform = FALSE ;

    /* ---------------------------------------------------------------------- */
    /* default ordering methods */
    /* ---------------------------------------------------------------------- */

    /* Note that if the Partition module is not installed, the CHOLMOD_METIS
     * and CHOLMOD_ND methods will not be available.  cholmod_analyze will
     * report the CHOLMOD_NOT_INSTALLED error, and safely skip over them.
     */

#if (CHOLMOD_MAXMETHODS < 8)
#error "CHOLMOD_MAXMETHODS must be 8 or more (defined in cholmod_core.h)."
#endif

    Common->nmethods = 3 ;	/* try methods 0, 1, and 2 and pick the best */
    Common->current = 0 ;
    Common->selected = 0 ;	/* the best method selected */

    /* first, fill each method with default parameters */
    for (i = 0 ; i <= CHOLMOD_MAXMETHODS ; i++)
    {
	/* CHOLMOD's default method is AMD for A and COLAMD for AA' */
	Common->method [i].ordering = CHOLMOD_AMD ;

	/* CHOLMOD's nested dissection (METIS + constrained AMD) */
	Common->method [i].nd_prune = 10.0 ;	/* dense node control */
	Common->method [i].nd_small = 200 ;	/* small graphs aren't cut */
	Common->method [i].nd_compress = TRUE ;	/* compress graph & subgraphs */
	Common->method [i].nd_camd = TRUE ;	/* use constrained AMD */

	/* statistics for each method */
	Common->method [i].fl = EMPTY ;
	Common->method [i].lnz = EMPTY ;
    }

    /* Next, define some methods.  The first five use default parameters. */
    Common->method [0].ordering = CHOLMOD_GIVEN ;   /* skip if UserPerm NULL */
    Common->method [1].ordering = CHOLMOD_AMD ;
    Common->method [2].ordering = CHOLMOD_METIS ;
    Common->method [3].ordering = CHOLMOD_ND ;
    Common->method [4].ordering = CHOLMOD_NATURAL ;

    /* CHOLMOD's nested dissection with large leaves of separator tree */
    Common->method [5].ordering = CHOLMOD_ND ;
    Common->method [5].nd_small = 20000 ;

    /* CHOLMOD's nested dissection with tiny leaves, and no AMD ordering */
    Common->method [6].ordering = CHOLMOD_ND ;
    Common->method [6].nd_small = 4 ;
    Common->method [6].nd_camd = FALSE ;

    /* CHOLMOD's nested dissection with no dense node removal */
    Common->method [7].ordering = CHOLMOD_ND ;
    Common->method [7].nd_prune = -1. ;

    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_finish ======================================================= */
/* ========================================================================== */

/* The last call to CHOLMOD must be cholmod_finish.  You may call this routine
 * more than once, and can safely call any other CHOLMOD routine after calling
 * it (including cholmod_start).
 *
 * The statistics and parameter settings in Common are preserved.  The
 * workspace in Common is freed.  This routine is just another name for
 * cholmod_free_work.
 */

int cholmod_finish
(
    cholmod_common *Common
)
{
    return (cholmod_free_work (Common)) ;
}


/* ========================================================================== */
/* === cholmod_allocate_work ================================================ */
/* ========================================================================== */

/* Allocate and initialize workspace for CHOLMOD routines, or increase the size
 * of already-allocated workspace.  If enough workspace is already allocated,
 * then nothing happens.
 *
 * workspace: Flag (nrow), Head (nrow+1), Iwork (iworksize), Xwork (xworksize)
 */

long cholmod_allocate_work  /* returns size of Xwork requested in units of
			     * sizeof (double), or -1 if failure */
(
    size_t nrow,		/* the A matrix will be nrow-by-ncol */
    size_t iworksize,		/* size of Iwork in sizeof (int)'s */
    size_t xwork,		/* size of Xwork in xunits */
    size_t xunits,		/* size of each item in Xwork */
    cholmod_common *Common
)
{
    double *W ;
    int *Head ;
    size_t wsize, xworksize ;
    int i ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (EMPTY) ;
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* Allocate Flag (nrow) and Head (nrow+1) */
    /* ---------------------------------------------------------------------- */

    nrow = MAX (1, nrow) ;
    if (nrow > Common->nrow)
    {
	/* free the old workspace (if any) and allocate new space */
	Common->Flag = cholmod_free (Common->Flag, Common->nrow, sizeof (int),
		Common) ;
	Common->Head = cholmod_free (Common->Head, Common->nrow+1, sizeof (int),
		Common) ;
	Common->Flag = cholmod_malloc (nrow,   sizeof (int), Common) ;
	Common->Head = cholmod_malloc (nrow+1, sizeof (int), Common) ;

	/* record the new size of Flag and Head */
	Common->nrow = nrow ;

	if (Common->status < CHOLMOD_OK)
	{
	    cholmod_free_work (Common) ;
	    return (EMPTY) ;
	}

	/* initialize Flag and Head */
	Common->mark = EMPTY ;
	(void) cholmod_clear_flag (Common) ;
	Head = Common->Head ;
	for (i = 0 ; i <= (int) (nrow) ; i++)
	{
	    Head [i] = EMPTY ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* Allocate Iwork (iworksize) */
    /* ---------------------------------------------------------------------- */

    iworksize = MAX (1, iworksize) ;
    if (iworksize > Common->iworksize)
    {
	/* free the old workspace (if any) and allocate new space.
	 * int overflow safely detected in cholmod_malloc */
	Common->Iwork = cholmod_free (Common->Iwork,
		Common->iworksize, sizeof (int), Common) ;
	Common->Iwork = cholmod_malloc (iworksize, sizeof (int), Common) ;

	/* record the new size of Iwork */
	Common->iworksize = iworksize ;

	if (Common->status < CHOLMOD_OK)
	{
	    cholmod_free_work (Common) ;
	    return (EMPTY) ;
	}

	/* note that Iwork does not need to be initialized */
    }

    /* ---------------------------------------------------------------------- */
    /* Allocate Xwork (xworksize) and set it to ((double) 0.) */
    /* ---------------------------------------------------------------------- */

    /* make sure xworksize is >= sizeof (double) */
    xwork = MAX (1, xwork) ;
    xunits = MAX (1, xunits) ;

    if (xwork > INT_MAX / xunits)
    {
	/* int overflow will occur */
	PRINT0 (("int overflow in Xwork\n")) ;
	cholmod_error (CHOLMOD_TOO_LARGE, "problem too large", Common) ;
	cholmod_free_work (Common) ;
	return (EMPTY) ;
    }

    xworksize = xwork * xunits ;

    /* round up xworksize to a multiple of sizeof (double) */
    xworksize = ROUNDUP (xworksize, sizeof (double)) ;
    wsize = xworksize / sizeof (double) ;

    if (((unsigned int) xworksize) > Common->xworksize)
    {
	/* free the old workspace (if any) and allocate new space */
	Common->Xwork = cholmod_free (Common->Xwork, Common->xworksize,
		sizeof (char), Common) ;

	if (Common->ieee)
	{
	    /* (case 1): set Xwork to zero by using calloc.  This is faster,
	     * but assumes ((double) 0.) is an all-zero bit field, which is not
	     * guaranteed.  This is the default, since CHOLMOD assumes IEEE
	     * floating-point. */
	    Common->Xwork = cholmod_calloc (xworksize, sizeof (char), Common) ;
	}
	else
	{
	    /* (case 2): if ((double) 0) is not a bit field of all zeros, use
	     * malloc and then set Xwork to zero below.  Not the default. */
	    Common->Xwork = cholmod_malloc (xworksize, sizeof (char), Common) ;
	}

	/* record the new size of Xwork */
	Common->xworksize = xworksize ;

	if (Common->status < CHOLMOD_OK)
	{
	    cholmod_free_work (Common) ;
	    return (EMPTY) ;	    /* out of memory */
	}

	if (!(Common->ieee))
	{
	    /* (case 2): malloc was used, so initialize Xwork to zero */
	    W = Common->Xwork ;
	    for (i = 0 ; i < (int) wsize ; i++)
	    {
		W [i] = 0. ;
	    }
	}
    }

    return (wsize) ;
}


/* ========================================================================== */
/* === cholmod_free_work ==================================================== */
/* ========================================================================== */

/* Deallocate the CHOLMOD workspace.
 *
 * workspace: deallocates all workspace in Common
 */

int cholmod_free_work	    /* always succeeds, returning TRUE */
(
    cholmod_common *Common
)
{
    RETURN_IF_NULL_COMMON (FALSE) ;
    DEBUG (cholmod_dump_freework ( )) ;
    Common->Flag  = cholmod_free (Common->Flag , Common->nrow, sizeof (int),
	    Common) ;
    Common->Head  = cholmod_free (Common->Head , Common->nrow+1, sizeof (int),
	    Common) ;
    Common->Iwork = cholmod_free (Common->Iwork, Common->iworksize,
	    sizeof (int), Common) ;
    Common->Xwork = cholmod_free (Common->Xwork, Common->xworksize,
	    sizeof (char), Common) ;
    Common->nrow = 0 ;
    Common->iworksize = 0 ;
    Common->xworksize = 0 ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_clear_flag =================================================== */
/* ========================================================================== */

/* Increment mark to ensure Flag [0..nrow-1] < mark.  If int overflow
 * occurs, or mark was initially negative, reset the entire array.  This is
 * not an error condition, but an intended function of the Flag workspace.
 *
 * workspace: Flag (nrow).  Does not modify Flag if nrow is zero.
 */

long cholmod_clear_flag
(
    cholmod_common *Common
)
{
    int i, nrow, *Flag ;

    RETURN_IF_NULL_COMMON (-1) ;

    Common->mark++ ;
    if (Common->mark <= 0)
    {
	nrow = Common->nrow ;
	Flag = Common->Flag ;
	PRINT2 (("reset Flag: nrow %d\n", nrow)) ;
	PRINT2 (("reset Flag: mark %ld\n", Common->mark)) ;
	for (i = 0 ; i < nrow ; i++)
	{
	    Flag [i] = EMPTY ;
	}
	Common->mark = 0 ;
    }
    return (Common->mark) ;
}


/* ========================================================================== */
/* ==== cholmod_maxrank ===================================================== */
/* ========================================================================== */

/* Find a valid value of Common->maxrank.  Returns 0 if error, or 2, 4, or 8
 * if successful. */

size_t cholmod_maxrank
(
    size_t n,
    cholmod_common *Common
)
{
    size_t maxrank ;
    RETURN_IF_NULL_COMMON (0) ;

    maxrank = Common->maxrank ;

    if (n > 0)
    {
	/* Ensure maxrank*n*sizeof(double) does not result in int overflow.
	 * If n is so large that 2*n*sizeof(double) results in int overflow
	 * (n = 268,435,455 if an int is 32 bits), then maxrank will be 0 or 1,
	 * but maxrank will be set to 2 below.  2*n will not result in int
	 * overflow, and CHOLMOD will run out of memory or safely detect int
	 * overflow elsewhere.
	 */
	maxrank = MIN (maxrank, INT_MAX / (n * sizeof (double))) ;
    }
    if (maxrank <= 2)
    {
	maxrank = 2 ;
    }
    else if (maxrank <= 4)
    {
	maxrank = 4 ;
    }
    else
    {
	maxrank = 8 ;
    }
    Common->maxrank = maxrank ;
    return (maxrank) ;
}


/* ========================================================================== */
/* === cholmod_dmin ========================================================= */
/* ========================================================================== */

/* Ensure the absolute value of a diagonal entry, D(j,j), is greater than
 * Common->dmin.  This routine is not meant for the user to call.  It is used
 * by the various LDL' factorization and update/downdate routines.  The
 * default value of Common->dmin is zero, and in that case this routine is not
 * called at all.  No change is made if D(j,j) or dmin are NaN.
 */

/* TODO keep track of the column at which the correction occurs */

double cholmod_dmin
(
    double dj,
    cholmod_common *Common
)
{
    double dmin ;
    RETURN_IF_NULL_COMMON (0) ;
    dmin = Common->dmin ;
    if (dj < 0)
    {
	if (dj > -dmin)
	{
	    dj = -dmin ;
	}
    }
    else
    {
	if (dj < dmin)
	{
	    dj = dmin ;
	}
    }
    return (dj) ;
}
