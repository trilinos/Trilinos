/* ========================================================================== */
/* === klu_dump ============================================================= */
/* ========================================================================== */

/* debug routines for klu / klu_btf */

#include "klu_btf_internal.h"

#ifndef NDEBUG

/* ========================================================================== */
/* === klu_valid ============================================================ */
/* ========================================================================== */

/* Check if a column-form matrix is valid or not.  The matrix A is
 * n-by-n.  The row indices of entries in column j are in
 * Ai [Ap [j] ... Ap [j+1]-1].  Required conditions are:
 *
 *	n >= 0
 *	nz = Ap [n_col] >= 0	    number of entries in the matrix
 *	Ap [0] == 0
 *	Ap [j] <= Ap [j+1] for all j in the range 0 to n_col.
 *	row indices in Ai [Ap [j] ... Ap [j+1]-1]
 *	    must be in the range 0 to n_row-1,
 *	    and no duplicate entries can exist (TODO not yet checked).
 *
 * Not user-callable.
 */

int klu_valid (int n, int Ap [ ], int Ai [ ], double Ax [ ])
{
    int nz, j, p1, p2, ilast, i, p ;
    PRINTF (("\ncolumn oriented matrix, n = %d\n", n)) ;
    if (n <= 0)
    {
	PRINTF (("n must be >= 0: %d\n", n)) ;
	return (FALSE) ;
    }
    nz = Ap [n] ;
    if (Ap [0] != 0 || nz < 0)
    {
	/* column pointers must start at Ap [0] = 0, and Ap [n] must be >= 0 */
	PRINTF (("column 0 pointer bad or nz < 0\n")) ;
	return (FALSE) ;
    }
    for (j = 0 ; j < n ; j++)
    {
	p1 = Ap [j] ;
	p2 = Ap [j+1] ;
	PRINTF (("\nColumn: %d p1: %d p2: %d\n", j, p1, p2)) ;
	if (p1 > p2)
	{
	    /* column pointers must be ascending */
	    PRINTF (("column %d pointer bad\n", j)) ;
	    return (FALSE) ;
	}
	for (p = p1 ; p < p2 ; p++)
	{
	    i = Ai [p] ;
	    PRINTF (("row: %d", i)) ;
	    if (i < 0 || i >= n)
	    {
		/* row index out of range */
		PRINTF (("index out of range, col %d row %d\n", j, i)) ;
		return (FALSE) ;
	    }
	    if (Ax != (double *) NULL)
	    {
		PRINTF ((" (%16.6g)", Ax [p])) ;
	    }
	    PRINTF (("\n")) ;
	}
    }
    return (TRUE) ;
}
#endif
