/* ========================================================================= */
/* === TRILINOS_CAMD_dump =========================================================== */
/* ========================================================================= */

/* ------------------------------------------------------------------------- */
/* CAMD, Copyright (c) Timothy A. Davis, Yanqing Chen,			     */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README.txt for License.     */
/* email: davis at cise.ufl.edu    CISE Department, Univ. of Florida.        */
/* web: http://www.cise.ufl.edu/research/sparse/camd                         */
/* ------------------------------------------------------------------------- */

/* Debugging routines for CAMD.  Not used if NDEBUG is not defined at compile-
 * time (the default).  See comments in camd_internal.h on how to enable
 * debugging.  Not user-callable.
 */

#include "trilinos_camd_internal.h"

#ifndef NDEBUG

/* This global variable is present only when debugging */
GLOBAL Int TRILINOS_CAMD_debug = -999 ;		/* default is no debug printing */

/* ========================================================================= */
/* === TRILINOS_CAMD_debug_init ===================================================== */
/* ========================================================================= */

/* Sets the debug print level, by reading the file debug.camd (if it exists) */

GLOBAL void TRILINOS_CAMD_debug_init ( char *s )
{
    FILE *f ;
    f = fopen ("debug.camd", "r") ;
    if (f == (FILE *) NULL)
    {
	TRILINOS_CAMD_debug = -999 ;
    }
    else
    {
	fscanf (f, ID, &TRILINOS_CAMD_debug) ;
	fclose (f) ;
    }
    if (TRILINOS_CAMD_debug >= 0)
    {
	printf ("%s: TRILINOS_CAMD_debug_init, D= "ID"\n", s, TRILINOS_CAMD_debug) ;
    }
}

/* ========================================================================= */
/* === TRILINOS_CAMD_dump =========================================================== */
/* ========================================================================= */

/* Dump CAMD's data structure, except for the hash buckets.  This routine
 * cannot be called when the hash buckets are non-empty.
 */

GLOBAL void TRILINOS_CAMD_dump (
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
    Int nel,
    Int BucketSet [ ],
    const Int C [ ],
    Int CurC
)
{
    Int i, pe, elen, nv, len, e, p, k, j, deg, w, cnt, ilast ;

    if (TRILINOS_CAMD_debug < 0) return ;
    ASSERT (pfree <= iwlen) ;
    TRILINOS_CAMD_DEBUG3 (("\nCAMD dump, pfree: "ID"\n", pfree)) ;
    for (i = 0 ; i < n ; i++)
    {
	pe = Pe [i] ;
	elen = Elen [i] ;
	nv = Nv [i] ;
	len = Len [i] ;
	w = W [i] ;

	if (elen >= TRILINOS_CAMD_EMPTY)
	{
	    if (nv == 0)
	    {
		TRILINOS_CAMD_DEBUG4 (("\nI "ID": nonprincipal:    ", i)) ;
		ASSERT (elen == TRILINOS_CAMD_EMPTY) ;
		if (pe == FLIP(n))
		{
		    TRILINOS_CAMD_DEBUG4 ((" dense node\n")) ;
		    ASSERT (w == 1) ;
		}
		else
		{
		    ASSERT (pe < TRILINOS_CAMD_EMPTY) ;
		    TRILINOS_CAMD_DEBUG4 ((" i "ID" -> parent "ID"\n", i, FLIP (Pe[i])));
		}
	    }
	    else
	    {
		TRILINOS_CAMD_DEBUG4 (("\nI "ID": active principal supervariable:\n",i));
		TRILINOS_CAMD_DEBUG4 (("   nv(i): "ID"  Flag: %d\n", nv, (nv < 0))) ;
		ASSERT (elen >= 0) ;
		ASSERT (nv > 0 && pe >= 0) ;
		p = pe ;
		TRILINOS_CAMD_DEBUG4 (("   e/s: ")) ;
		if (elen == 0) TRILINOS_CAMD_DEBUG4 ((" : ")) ;
		ASSERT (pe + len <= pfree) ;
		for (k = 0 ; k < len ; k++)
		{
		    j = Iw [p] ;
		    TRILINOS_CAMD_DEBUG4 (("  "ID"", j)) ;
		    ASSERT (j >= 0 && j < n) ;
		    if (k == elen-1) TRILINOS_CAMD_DEBUG4 ((" : ")) ;
		    p++ ;
		}
		TRILINOS_CAMD_DEBUG4 (("\n")) ;
	    }
	}
	else
	{
	    e = i ;
	    if (w == 0)
	    {
		TRILINOS_CAMD_DEBUG4 (("\nE "ID": absorbed element: w "ID"\n", e, w)) ;
		ASSERT (nv > 0 && pe < 0) ;
		TRILINOS_CAMD_DEBUG4 ((" e "ID" -> parent "ID"\n", e, FLIP (Pe [e]))) ;
	    }
	    else
	    {
		TRILINOS_CAMD_DEBUG4 (("\nE "ID": unabsorbed element: w "ID"\n", e, w)) ;
		ASSERT (nv > 0 && pe >= 0) ;
		p = pe ;
		TRILINOS_CAMD_DEBUG4 ((" : ")) ;
		ASSERT (pe + len <= pfree) ;
		for (k = 0 ; k < len ; k++)
		{
		    j = Iw [p] ;
		    TRILINOS_CAMD_DEBUG4 (("  "ID"", j)) ;
		    ASSERT (j >= 0 && j < n) ;
		    p++ ;
		}
		TRILINOS_CAMD_DEBUG4 (("\n")) ;
	    }
	}
	TRILINOS_CAMD_DEBUG4 (("C[i] is :"ID"\n", (C == NULL) ? 0 : C [i]));
    }

    /* this routine cannot be called when the hash buckets are non-empty */
    TRILINOS_CAMD_DEBUG4 (("\nDegree lists:\n")) ;
    if (nel >= 0)
    {
	cnt = 0 ;
	for (deg = 0 ; deg < n ; deg++)
	{
	    if (Head [deg] == TRILINOS_CAMD_EMPTY) continue ;
	    ilast = TRILINOS_CAMD_EMPTY ;
	    TRILINOS_CAMD_DEBUG4 ((ID": \n", deg)) ;
	    for (i = Head [deg] ; i != TRILINOS_CAMD_EMPTY ; i = Next [i])
	    {
		TRILINOS_CAMD_DEBUG4 (("   "ID" : next "ID" last "ID" deg "ID"\n",
		    i, Next [i], Last [i], Degree [i])) ;
		ASSERT (i >= 0 && i < n && ilast == Last [i] &&
		    deg == Degree [i]) ;
		cnt += Nv [i] ;
		ilast = i ;
	    }
	    TRILINOS_CAMD_DEBUG4 (("\n")) ;
	}
    }
    
    TRILINOS_CAMD_DEBUG4(("\nCurrent C[i] is "ID". current Buckets are:\n", CurC)) ;
    for (i = 0 ; i < n ; i++)
    {
	if ((C == NULL) ? 1 : (C [BucketSet [i]] <= CurC))
            TRILINOS_CAMD_DEBUG4((ID",",BucketSet [i]));
    }
    TRILINOS_CAMD_DEBUG4 (("\n")) ;
}

#endif
