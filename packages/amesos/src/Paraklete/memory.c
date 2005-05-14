#include "pk.h"

/* ========================================================================== */
/* === memory handler ======================================================= */
/* ========================================================================== */

/* my_malloc2, my_calloc2, and my_realloc2 pretend to fail if my_tries goes to
 * zero, to test CHOLMOD's memory error handling.   No failure occurs if
 * my_tries is negative.
 *
 * PARAKLETE version 0.1: parallel sparse LU factorization.  May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

int my_tries = -1 ;

void *my_malloc2 (size_t size)
{
    void *p ;
    if (my_tries == 0)
    {
	/* pretend to fail */
	return (NULL) ;
    }
    if (my_tries > 0)
    {
	my_tries-- ;
    }
    p = malloc (size) ;
    return (p) ;
}

void *my_calloc2 (size_t n, size_t size)
{
    void *p ;
    if (my_tries == 0)
    {
	/* pretend to fail */
	return (NULL) ;
    }
    if (my_tries > 0)
    {
	my_tries-- ;
    }
    p = calloc (n, size) ;
    return (p) ;
}

void *my_realloc2 (void *p, size_t size)
{
    if (my_tries == 0)
    {
	/* pretend to fail */
	return (NULL) ;
    }
    if (my_tries > 0)
    {
	my_tries-- ;
    }
    return (realloc (p, size)) ;
}

void my_free2 (void *p)
{
    free (p) ;
}

void normal_memory_handler ( cholmod_common *cm )
{
    cm->malloc_memory = malloc ;
    cm->calloc_memory = calloc ;
    cm->realloc_memory = realloc ;
    cm->free_memory = free ;
    cm->error_handler = my_handler ;
    cholmod_free_work (cm) ;
}

void test_memory_handler ( cholmod_common *cm )
{
    cm->malloc_memory = my_malloc2 ;
    cm->calloc_memory = my_calloc2 ;
    cm->realloc_memory = my_realloc2 ;
    cm->free_memory = my_free2 ;
    cm->error_handler = NULL ;
    cholmod_free_work (cm) ;
    my_tries = 0 ;
}
