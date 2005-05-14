/* ========================================================================== */
/* === Core/cholmod_memory ================================================== */
/* ========================================================================== */

/*
 * CHOLMOD/Core version 0.1. May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* Core memory management routines:
 *
 * Primary routines:
 * -----------------
 * cholmod_malloc		malloc wrapper
 * cholmod_free			free wrapper
 *
 * Secondary routines:
 * -------------------
 * cholmod_calloc		calloc wrapper
 * cholmod_realloc		realloc wrapper
 * cholmod_realloc_multiple	realloc wrapper for multiple objects
 *
 * The void * pointers here are truly void *.  For all other routines in
 * CHOLMOD, void * parameters are used to point to integer arrays of type
 * int or long, depending on Common->itype.
 *
 * The user may make use of these, just like malloc and free.  You can even
 * malloc an object and safely free it with cholmod_free, and visa versa
 * (except that the memory usage statistics will be corrupted).  These routines
 * do differ from malloc and free.  If cholmod_free is given a NULL pointer,
 * for example, it does nothing (unlike the ANSI free).  cholmod_realloc does
 * not return NULL if given a non-NULL pointer and a nonzero size, even if it
 * fails (it sets an error code in Common->status instead).
 *
 * CHOLMOD keeps track of the amount of memory it has allocated, and so the
 * cholmod_free routine includes as a parameter the size of the object being
 * freed.  This is only used for memory usage statistics, which are very useful
 * in finding memory leaks in your program.  If you, the user of CHOLMOD, pass
 * the wrong size, the only consequence is that the memory usage statistics
 * will be corrupted.  This will causes assertions to fail if CHOLMOD is
 * compiled with debugging enabled.
 *
 * The cholmod_free_* routines for each CHOLMOD object keep track of the size
 * of the blocks they free, so they do not require you to pass their sizes
 * as a parameter.
 */

#include "cholmod_core.h"
#include "cholmod_internal.h"

/* ========================================================================== */
/* === cholmod_malloc ======================================================= */
/* ========================================================================== */

/* Wrapper around malloc routine (mxMalloc for a mexFunction).  Allocates
 * space of size MAX(1,n)*size, where size is normally a sizeof (...).
 *
 * This routine, cholmod_calloc, and cholmod_realloc do not set Common->status
 * to CHOLMOD_OK on success, so that a sequence of cholmod_malloc's, _calloc's,
 * or _realloc's can be used.  If any of them fails, the Common->status will
 * hold the most recent error status.
 *
 * Usage, for a pointer to int:
 *
 *	p = cholmod_malloc (n, sizeof (int), Common) ;
 *
 * Uses a pointer to the malloc routine (or its equivalent) defined in Common.
 */

void *cholmod_malloc
(
    size_t n,			/* number of items to allocate */
    size_t size,		/* size of each item */
    cholmod_common *Common
)
{
    void *p ;

    RETURN_IF_NULL_COMMON (NULL) ;
    if (size == 0)
    {
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_malloc: sizeof(item) must be > 0", Common) ;
	p = NULL ;
    }
    else if (n >= (INT_MAX / size))
    {
	/* object is too big to allocate without causing int overflow */
	cholmod_error (CHOLMOD_TOO_LARGE, "problem too large", Common) ;
	p = NULL ;
    }
    else
    {
	/* call malloc, or its equivalent */
	p = (Common->malloc_memory) (MAX (1,n) * size) ;
	if (p == NULL)
	{
	    /* failure: out of memory */
	    cholmod_error (CHOLMOD_OUT_OF_MEMORY, "out of memory", Common) ;
	}
	else
	{
	    /* success: increment the count of objects allocated */
	    Common->malloc_count++ ;
	    Common->memory_inuse += (n * size) ;
	    Common->memory_usage =
		MAX (Common->memory_usage, Common->memory_inuse) ;

	    PRINTM (("cholmod_malloc %p %d cnt: %d inuse %d\n",
		    p, n*size, Common->malloc_count, Common->memory_inuse)) ;
	}
    }
    return (p) ;
}


/* ========================================================================== */
/* === cholmod_free ========================================================= */
/* ========================================================================== */

/* Wrapper around free routine (mxFree for a mexFunction).  Returns NULL,
 * which should be assigned to the pointer being freed, as in:
 *
 *	p = cholmod_free (p, n, sizeof (int), Common) ;
 *
 * In CHOLMOD, the syntax:
 *
 *	cholmod_free (p, n, sizeof (int), Common) ;
 *
 * is used if p is a local pointer and the routine is returning shortly.
 * Uses a pointer to the free routine (or its equivalent) defined in Common.
 */

void *cholmod_free
(
    void *p,
    size_t n,
    size_t size,
    cholmod_common *Common
)
{
    RETURN_IF_NULL_COMMON (NULL) ;
    if (p != NULL)
    {
	/* only free the object if the pointer is not NULL */
	/* call free, or its equivalent */
	(Common->free_memory) (p) ;
	Common->malloc_count-- ;
	Common->memory_inuse -= (n * size) ;
	PRINTM (("cholmod_free   %p %d cnt: %d inuse %d\n",
		p, n*size, Common->malloc_count, Common->memory_inuse)) ;
	/* This assertion will fail if the user calls cholmod_malloc and
	 * cholmod_free with mismatched memory sizes.  It shouldn't fail
	 * otherwise. */
	DEBUG (if (Common->malloc_count == 0 && Common->memory_inuse != 0)
	    PRINT0 (("inuse: %d\n", Common->memory_inuse))) ;
	ASSERT (IMPLIES (Common->malloc_count == 0, Common->memory_inuse == 0));
    }
    /* return NULL, and the caller should assign this to p.  This avoids
     * freeing the same pointer twice. */
    return (NULL) ;
}


/* ========================================================================== */
/* === cholmod_calloc ======================================================= */
/* ========================================================================== */

/* Wrapper around calloc routine (mxCalloc for a mexFunction).
 *
 * Uses a pointer to the calloc routine (or its equivalent) defined in Common.
 * This routine is identical to malloc, except that it zeros the newly allocated
 * block to zero.
 */

void *cholmod_calloc
(
    size_t n,			/* number of items to allocate */
    size_t size,		/* size of each item */
    cholmod_common *Common
)
{
    void *p ;

    RETURN_IF_NULL_COMMON (NULL) ;
    if (size == 0)
    {
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_calloc: sizeof(item) must be > 0", Common) ;
	p = NULL ;
    }
    else if (n >= (INT_MAX / size))
    {
	/* object is too big to allocate without causing int overflow */
	cholmod_error (CHOLMOD_TOO_LARGE, "problem too large", Common) ;
	p = NULL ;
    }
    else
    {
	/* call calloc, or its equivalent */
	p = (Common->calloc_memory) (MAX (1,n), size) ;
	if (p == NULL)
	{
	    /* failure: out of memory */
	    cholmod_error (CHOLMOD_OUT_OF_MEMORY, "out of memory", Common) ;
	}
	else
	{
	    /* success: increment the count of objects allocated */
	    Common->malloc_count++ ;
	    Common->memory_inuse += (n * size) ;
	    Common->memory_usage =
		MAX (Common->memory_usage, Common->memory_inuse) ;
	    PRINTM (("cholmod_calloc %p %d cnt: %d inuse %d\n",
		    p, n*size, Common->malloc_count, Common->memory_inuse)) ;
	}
    }
    return (p) ;
}


/* ========================================================================== */
/* === cholmod_realloc ====================================================== */
/* ========================================================================== */

/* Wrapper around realloc routine (mxRealloc for a mexFunction).  Given a
 * pointer p to a block of size (*n)*size memory, it changes the size of the
 * block pointed to by p to be MAX(1,new)*size in size.  It may return a
 * pointer different than p.  This should be used as (for a pointer to int):
 *
 *	p = cholmod_realloc (p, *n, new, sizeof (int), Common) ;
 *
 * If p is NULL, this is the same as p = cholmod_malloc (...).
 * A size of new=0 is treated as new=1.
 *
 * If the realloc fails, p is returned unchanged and Common->status is set
 * to CHOLMOD_OUT_OF_MEMORY.  If successful, Common->status is not modified,
 * and p is returned (possibly changed) and pointing to a large block of memory.
 *
 * Uses a pointer to the realloc routine (or its equivalent) defined in Common.
 */

void *cholmod_realloc
(
    void *p,
    size_t *n,	/* on input, n = current size. on output, n = new size */
    size_t new,	/* requested size of the block of memory */
    size_t size,		/* size of each item */
    cholmod_common *Common
)
{
    size_t old = (*n) ;
    void *pnew ;

    RETURN_IF_NULL_COMMON (NULL) ;
    if (size == 0)
    {
	cholmod_error (CHOLMOD_INVALID,
		"cholmod_realloc: sizeof(item) must be > 0", Common) ;
	p = NULL ;
    }
    else if (p == NULL)
    {
	/* A fresh object is being allocated. */
	PRINT1 (("realloc fresh: %d %d\n", new, size)) ;
	p = cholmod_malloc (new, size, Common) ;
	*n = (p == NULL) ? 0 : new ;
    }
    else if (old == new)
    {
	/* Nothing to do.  Do not change p or n. */
	PRINT1 (("realloc nothing: %d %d\n", new, size)) ;
    }
    else if (new >= (INT_MAX / size))
    {
	/* failure: The new size too big.  Do not change p or n. */
	cholmod_error (CHOLMOD_TOO_LARGE, "problem too large", Common) ;
    }
    else
    {
	/* The object exists, and is changing to some other nonzero size. */
	/* call realloc, or its equivalent */
	PRINT1 (("realloc : %d to %d, %d\n", old, MAX(1,new), size)) ;
	pnew = NULL ;
	/*
	pnew = realloc (p, MAX (1,new) * size) ;
	 */
	pnew = (Common->realloc_memory) (p, MAX (1,new) * size) ;
	if (pnew == NULL)
	{
	    /* Do not change p, since it still points to allocated memory */
	    if (new <= old)
	    {
		/* The attempt to reduce the size of the block from n to
		 * new has failed.  The current block is not modified, so
		 * pretend to succeed, but do not change p.  Do change
		 * CHOLMOD's notion of the size of the block, however. */
		*n = new ;
		PRINTM (("cholmod_realloc_old: %p %d cnt: %d inuse %d\n"
			"cholmod_realloc_new: %p %d cnt: %d inuse %d\n",
		    p, old*size,   Common->malloc_count-1,
				    Common->memory_inuse - old*size,
		    pnew, new*size, Common->malloc_count,
				    Common->memory_inuse + (new-old)*size)) ;
		Common->memory_inuse += ((new-old) * size) ;
	    }
	    else
	    {
		/* Increasing the size of the block has failed.
		 * Do not change n. */
		cholmod_error (CHOLMOD_OUT_OF_MEMORY, "out of memory", Common) ;
	    }
	}
	else
	{
	    /* success: return the new p and change the size of the block */

	    PRINTM (("cholmod_realloc_old: %p %d cnt: %d inuse %d\n"
		    "cholmod_realloc_new: %p %d cnt: %d inuse %d\n",
		p, old*size,   Common->malloc_count-1,
				Common->memory_inuse - old*size,
		pnew, new*size, Common->malloc_count,
				Common->memory_inuse + (new-old)*size)) ;
	    p = pnew ;
	    *n = new ;
	    Common->memory_inuse += ((new-old) * size) ;
	}
	Common->memory_usage = MAX (Common->memory_usage, Common->memory_inuse);
    }

    return (p) ;
}


/* ========================================================================== */
/* === cholmod_realloc_multiple ============================================= */
/* ========================================================================== */

/* Reallocate multiple blocks of memory, all of the same size (up to two integer
 * and two real blocks).  Either reallocations all succeed, or all are returned
 * in the original size (they are freed if the original size is zero).
 *
 * FUTURE WORK: only int and double are currently supported.
 */

int cholmod_realloc_multiple
(
    int nint,	    /* number int/long blocks to allocate */
    int nreal,	    /* number of double/complex/float blocks to realloc */

    void **I,	    /* int or long block */
    void **J,	    /* int or long block */

    void **X,	    /* complex, double, or float block */
    void **Z,	    /* double or float block */

    size_t *old_p,  /* current size of both blocks */
    size_t new,	    /* requested size */
    cholmod_common *Common
)
{
    size_t i, j, x, z, old ;

    RETURN_IF_NULL_COMMON (FALSE) ;
    old = *old_p ;

    if (nint < 1 && nreal < 1)
    {
	/* nothing to do */
	return (TRUE) ;
    }

    if (nreal > 1 || Z != NULL)
    {
	cholmod_error (CHOLMOD_INVALID, "complex case not supported", Common) ;
	return (FALSE) ;
    }

    i = old ;
    j = old ;
    x = old ;
    z = old ;

    if (nint > 0)
    {
	*I = cholmod_realloc (*I, &i, new, sizeof (int), Common) ;
    }
    if (nint > 1)
    {
	*J = cholmod_realloc (*J, &j, new, sizeof (int), Common) ;
    }
    if (nreal > 0)
    {
	*X = cholmod_realloc (*X, &x, new, sizeof (double), Common) ;
    }

    /* FUTURE WORK:
    if (nreal > 1)
    {
	*Z = cholmod_realloc (*Z, &z, new, sizeof (double), Common) ;
    }
    */

    if (Common->status < CHOLMOD_OK)
    {
	/* one or more realloc's failed.  Resize all back down to old. */
	ASSERT (new > old) ;

	if (old == 0)
	{

	    if (nint > 0)
	    {
		*I = cholmod_free (*I, i, sizeof (int), Common) ;
	    }
	    if (nint > 1)
	    {
		*J = cholmod_free (*J, j, sizeof (int), Common) ;
	    }
	    if (nreal > 0)
	    {
		*X = cholmod_free (*X, x, sizeof (double), Common) ;
	    }
	    /*
	    if (nreal > 1)
	    {
		*Z = cholmod_free (*Z, z, sizeof (double), Common) ;
	    }
	    */

	}
	else
	{

	    if (nint > 0)
	    {
		*I = cholmod_realloc (*I, &i, old, sizeof (int), Common) ;
	    }
	    if (nint > 1)
	    {
		*J = cholmod_realloc (*J, &j, old, sizeof (int), Common) ;
	    }
	    if (nreal > 0)
	    {
		*X = cholmod_realloc (*X, &x, old, sizeof (double), Common) ;
	    }
	    /*
	    if (nreal > 1)
	    {
		*Z = cholmod_realloc (*Z, &z, old, sizeof (double), Common) ;
	    }
	    */

	}

	return (FALSE) ;
    }
    else
    {
	/* all realloc's succeeded, change size to reflect new size. */
	*old_p = new ;
	return (TRUE) ;
    }
}
