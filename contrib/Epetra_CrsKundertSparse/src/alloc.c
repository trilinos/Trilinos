/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/
#ifndef lint
static char *cvs_file_id = "$Id$";
#endif

/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
**********/

/*
 * Memory alloction functions
 */

#include "spice.h"
#include "stdio.h"
#include "misc.h"
#include "suffix.h"
#include "util.h"
#ifdef DEBUG_MALLOC
#include "sm.h"
#endif

#include <stdio.h>

/* Malloc num bytes and initialize to zero. Fatal error if the space can't
 * be malloc'd.   Return NULL for a request for 0 bytes.
 */
#undef SHARED_MEM

void bye_bye(i)
{
    printf ("inv = %d\n",1/i);
}
/*
*/
char *
tmalloc(num)
    int num;
{
    char *s;

    if (!num)
	return NULL;

#ifdef DEBUG_MALLOC
    s = sm_malloc((unsigned) num);
#else
    s = malloc((unsigned) num);
#endif
    if (!s) {
        fprintf(stderr, 
		"malloc: Internal Error: can't allocate %d bytes.\n", num);
        exit(EXIT_BAD);
    }

    bzero(s, num);

    return(s);
}

char *
trealloc(str, num)
    char *str;
    int num;
{
    char *s;

    if (!num) {
	if (str)
#ifdef SHARED_MEM
		FREE(str);
#else
#ifdef DEBUG_MALLOC
		sm_free(str);
#else
		free(str);
#endif
#endif
	return NULL;
    }
#ifdef SHARED_MEM
    if (!str)
	s = (char *) MALLOC_SM(num);
    else
        s = (char *) REALLOC_SM(str, (unsigned) num);
#else
    if (!str)
	s = tmalloc(num);
    else
#ifdef DEBUG_MALLOC
        s = sm_realloc(str, (unsigned) num);
#else
        s = realloc(str, (unsigned) num);
#endif
#endif
    if (!s) {
        fprintf(stderr, 
		"realloc: Internal Error: can't allocate %d bytes.\n", num);
        perror ("realloc");
        s = malloc((unsigned) num);
        bye_bye(0);
        fprintf (stderr, "From malloc of %d bytes: %lx\n",num,s);
        perror ("malloc");
        exit(EXIT_BAD);
    }
    return(s);
}

void
txfree(ptr)
	char	*ptr;
{
	if (ptr)
#ifdef SHARED_MEM
		FREE(ptr);
#else
#ifdef DEBUG_MALLOC
		sm_free(ptr);
#else
		free(ptr);
#endif
#endif
}
