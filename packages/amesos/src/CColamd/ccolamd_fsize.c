/* ========================================================================== */
/* === CCOLAMD_fsize ======================================================== */
/* ========================================================================== */

/*
 * CCOLAMD version 0.1, May 13, 2005
 * Copyright (C) 2005, Univ. of Florida.  Authors: Timothy A. Davis
 * and Sivasankaran Rajamanickam
 * See License.txt for the Version 2.1 of the GNU Lesser General Public License
 * http://www.cise.ufl.edu/research/sparse
 */

/* Determine the largest frontal matrix size for each subtree. 
 * Only required to sort the children of each
 * node prior to AMD_postorder. */

/* -------------------- */
/* modified for ccolamd */
#include "ccolamd.h"
#include "ccolamd_internal.h"
extern Int colamd_debug  ;
/* -------------------- */

GLOBAL void CCOLAMD_fsize
(
    Int nn,
    Int Fsize [ ],
    Int Fnrows [ ],
    Int Fncols [ ],
    Int Parent [ ],
    Int Npiv [ ]
)
{
    double dr, dc ;
    Int j, parent, frsize, r, c ;

    for (j = 0 ; j < nn ; j++)
    {
	Fsize [j] = EMPTY ;
    }

    /* ---------------------------------------------------------------------- */
    /* find max front size for tree rooted at node j, for each front j */
    /* ---------------------------------------------------------------------- */

    DEBUG1 (("\n\n========================================FRONTS:\n")) ;
    for (j = 0 ; j < nn ; j++)
    {
	if (Npiv [j] > 0)
	{
	    /* this is a frontal matrix */
	    parent = Parent [j] ;
	    r = Fnrows [j] ;
	    c = Fncols [j] ;
	    /* avoid integer overflow */
	    dr = (double) r ;
	    dc = (double) c ;
	    frsize = (INT_OVERFLOW (dr * dc)) ?  Int_MAX : (r * c) ;
	    DEBUG1 ((""ID" : npiv "ID" size "ID" parent "ID" ",
		j, Npiv [j], frsize, parent)) ;
	    Fsize [j] = MAX (Fsize [j], frsize) ;
	    DEBUG1 (("Fsize [j = "ID"] = "ID"\n", j, Fsize [j])) ;
	    if (parent != EMPTY)
	    {
		/* find the maximum frontsize of self and children */
		ASSERT (Npiv [parent] > 0) ;
		ASSERT (parent > j) ;
		Fsize [parent] = MAX (Fsize [parent], Fsize [j]) ;
		DEBUG1 (("Fsize [parent = "ID"] = "ID"\n",
		    parent, Fsize [parent]));
	    }
	}
    }
    DEBUG1 (("fsize done\n")) ;
}
