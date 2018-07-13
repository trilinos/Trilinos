/* ========================================================================= */
/* === TRILINOS_AMD_dump ============================================================ */
/* ========================================================================= */

/* ------------------------------------------------------------------------- */
/* AMD, Copyright (c) Timothy A. Davis,					     */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README.txt for License.     */
/* email: davis at cise.ufl.edu    CISE Department, Univ. of Florida.        */
/* web: http://www.cise.ufl.edu/research/sparse/amd                          */
/* ------------------------------------------------------------------------- */

/* Debugging routines for AMD.  Not used if NDEBUG is not defined at compile-
 * time (the default).  See comments in trilinos_amd_internal.h on how to enable
 * debugging.  Not user-callable.
 */

#include "trilinos_amd_internal.h"

#ifndef NDEBUG

/* This global variable is present only when debugging */
GLOBAL Int TRILINOS_AMD_debug = -999 ;		/* default is no debug printing */

/* ========================================================================= */
/* === TRILINOS_AMD_debug_init ====================================================== */
/* ========================================================================= */

/* Sets the debug print level, by reading the file debug.amd (if it exists) */

GLOBAL void TRILINOS_AMD_debug_init ( char *s )
{
    FILE *f ;
    f = fopen ("debug.amd", "r") ;
    if (f == (FILE *) NULL)
    {
	TRILINOS_AMD_debug = -999 ;
    }
    else
    {
	fscanf (f, ID, &TRILINOS_AMD_debug) ;
	fclose (f) ;
    }
    if (TRILINOS_AMD_debug >= 0)
    {
	printf ("%s: TRILINOS_AMD_debug_init, D= "ID"\n", s, TRILINOS_AMD_debug) ;
    }
}

/* ========================================================================= */
/* === TRILINOS_AMD_dump ============================================================ */
/* ========================================================================= */

/* Dump AMD's data structure, except for the hash buckets.  This routine
 * cannot be called when the hash buckets are non-empty.
 */

GLOBAL void TRILINOS_AMD_dump (
    Int n,	    /* A is n-by-n */
    Int Pe [ ],	    /* pe [0..n-1]: index in iw of start of row i */
    Int Iw [ ],	    /* workspace of size iwlen, iwlen [0..pfree-1]
		     * holds the matrix on input */
    Int Len [ ],    /* len [0..n-1]: length for row i */
    Int iwlen,	    /* length of iw */
    Int pfree,	    /* iw [pfree ... iwlen-1] is empty on input */
    Int Nv [ ],	    /* nv [0..n-1] */
    Int Next [ ],   /* next [0..n-1] */
    Int Last [ ],   /* last [0..n-1] */
    Int Head [ ],   /* head [0..n-1] */
    Int Elen [ ],   /* size n */
    Int Degree [ ], /* size n */
    Int W [ ],	    /* size n */
    Int nel
)
{
    Int i, pe, elen, nv, len, e, p, k, j, deg, w, cnt, ilast ;

    if (TRILINOS_AMD_debug < 0) return ;
    ASSERT (pfree <= iwlen) ;
    TRILINOS_AMD_DEBUG3 (("\nAMD dump, pfree: "ID"\n", pfree)) ;
    for (i = 0 ; i < n ; i++)
    {
	pe = Pe [i] ;
	elen = Elen [i] ;
	nv = Nv [i] ;
	len = Len [i] ;
	w = W [i] ;

	if (elen >= EMPTY)
	{
	    if (nv == 0)
	    {
		TRILINOS_AMD_DEBUG3 (("\nI "ID": nonprincipal:    ", i)) ;
		ASSERT (elen == EMPTY) ;
		if (pe == EMPTY)
		{
		    TRILINOS_AMD_DEBUG3 ((" dense node\n")) ;
		    ASSERT (w == 1) ;
		}
		else
		{
		    ASSERT (pe < EMPTY) ;
		    TRILINOS_AMD_DEBUG3 ((" i "ID" -> parent "ID"\n", i, FLIP (Pe[i])));
		}
	    }
	    else
	    {
		TRILINOS_AMD_DEBUG3 (("\nI "ID": active principal supervariable:\n",i));
		TRILINOS_AMD_DEBUG3 (("   nv(i): "ID"  Flag: %d\n", nv, (nv < 0))) ;
		ASSERT (elen >= 0) ;
		ASSERT (nv > 0 && pe >= 0) ;
		p = pe ;
		TRILINOS_AMD_DEBUG3 (("   e/s: ")) ;
		if (elen == 0) TRILINOS_AMD_DEBUG3 ((" : ")) ;
		ASSERT (pe + len <= pfree) ;
		for (k = 0 ; k < len ; k++)
		{
		    j = Iw [p] ;
		    TRILINOS_AMD_DEBUG3 (("  "ID"", j)) ;
		    ASSERT (j >= 0 && j < n) ;
		    if (k == elen-1) TRILINOS_AMD_DEBUG3 ((" : ")) ;
		    p++ ;
		}
		TRILINOS_AMD_DEBUG3 (("\n")) ;
	    }
	}
	else
	{
	    e = i ;
	    if (w == 0)
	    {
		TRILINOS_AMD_DEBUG3 (("\nE "ID": absorbed element: w "ID"\n", e, w)) ;
		ASSERT (nv > 0 && pe < 0) ;
		TRILINOS_AMD_DEBUG3 ((" e "ID" -> parent "ID"\n", e, FLIP (Pe [e]))) ;
	    }
	    else
	    {
		TRILINOS_AMD_DEBUG3 (("\nE "ID": unabsorbed element: w "ID"\n", e, w)) ;
		ASSERT (nv > 0 && pe >= 0) ;
		p = pe ;
		TRILINOS_AMD_DEBUG3 ((" : ")) ;
		ASSERT (pe + len <= pfree) ;
		for (k = 0 ; k < len ; k++)
		{
		    j = Iw [p] ;
		    TRILINOS_AMD_DEBUG3 (("  "ID"", j)) ;
		    ASSERT (j >= 0 && j < n) ;
		    p++ ;
		}
		TRILINOS_AMD_DEBUG3 (("\n")) ;
	    }
	}
    }

    /* this routine cannot be called when the hash buckets are non-empty */
    TRILINOS_AMD_DEBUG3 (("\nDegree lists:\n")) ;
    if (nel >= 0)
    {
	cnt = 0 ;
	for (deg = 0 ; deg < n ; deg++)
	{
	    if (Head [deg] == EMPTY) continue ;
	    ilast = EMPTY ;
	    TRILINOS_AMD_DEBUG3 ((ID": \n", deg)) ;
	    for (i = Head [deg] ; i != EMPTY ; i = Next [i])
	    {
		TRILINOS_AMD_DEBUG3 (("   "ID" : next "ID" last "ID" deg "ID"\n",
		    i, Next [i], Last [i], Degree [i])) ;
		ASSERT (i >= 0 && i < n && ilast == Last [i] &&
		    deg == Degree [i]) ;
		cnt += Nv [i] ;
		ilast = i ;
	    }
	    TRILINOS_AMD_DEBUG3 (("\n")) ;
	}
	ASSERT (cnt == n - nel) ;
    }

}

#endif
