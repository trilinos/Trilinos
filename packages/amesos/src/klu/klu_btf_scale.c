/* ========================================================================== */
/* === klu_btf_sacle ======================================================== */
/* ========================================================================== */

/* Scale a matrix.  Can be called by the user, but not needed.
 * This is called by klu_btf_factor and klu_btf_refactor.
 */

#include "klu_btf.h"
#include "klu_kernel.h"
#include "klu_dump.h"

int klu_btf_scale
(
    /* inputs, not modified */
    int n,
    int Ap [ ],		/* size n+1, column pointers */
    int Ai [ ],		/* size nz, row indices */
    double Ax [ ],
    /* outputs, not defined on input */
    double Rs [ ]
)
{
    int row, col, p, pend ;

    for (row = 0 ; row < n ; row++)
    {
	Rs [row] = 0.0 ;
    }

    /* the scale factor for row i is sum (abs (A (i,:))) */
    for (col = 0 ; col < n ; col++)
    {
	pend = Ap [col+1] ;
	for (p = Ap [col] ; p < pend ; p++)
	{
	    Rs [Ai [p]] += ABS (Ax [p]) ;
	}
    }

    for (row = 0 ; row < n ; row++)
    {
	/* matrix is singular */
	PRINTF (("Rs [%d] = %g\n", row, Rs [row])) ;
	if (Rs [row] == 0.0)
	{
	    PRINTF (("Row %d of A is all zero\n", row)) ;
	    return (FALSE) ;
	}
    }
    return (TRUE) ;
}
