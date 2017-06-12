/* ========================================================================= */
/* === TRILINOS_AMD_info ============================================================ */
/* ========================================================================= */

/* ------------------------------------------------------------------------- */
/* AMD, Copyright (c) Timothy A. Davis,					     */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README.txt for License.     */
/* email: davis at cise.ufl.edu    CISE Department, Univ. of Florida.        */
/* web: http://www.cise.ufl.edu/research/sparse/amd                          */
/* ------------------------------------------------------------------------- */

/* User-callable.  Prints the output statistics for AMD.  See amd.h
 * for details.  If the Info array is not present, nothing is printed.
 */

#include "trilinos_amd_internal.h"

#define PRI(format,x) { if (x >= 0) { PRINTF ((format, x)) ; }}

GLOBAL void TRILINOS_AMD_info
(
    double Info [ ]
)
{
    double n, ndiv, nmultsubs_ldl, nmultsubs_lu, lnz, lnzd ;

    PRINTF (("\nAMD version %d.%d.%d, %s, results:\n",
	TRILINOS_AMD_MAIN_VERSION, TRILINOS_AMD_SUB_VERSION, TRILINOS_AMD_SUBSUB_VERSION, TRILINOS_AMD_DATE)) ;

    if (!Info)
    {
	return ;
    }

    n = Info [TRILINOS_AMD_N] ;
    ndiv = Info [TRILINOS_AMD_NDIV] ;
    nmultsubs_ldl = Info [TRILINOS_AMD_NMULTSUBS_LDL] ;
    nmultsubs_lu = Info [TRILINOS_AMD_NMULTSUBS_LU] ;
    lnz = Info [TRILINOS_AMD_LNZ] ;
    lnzd = (n >= 0 && lnz >= 0) ? (n + lnz) : (-1) ;

    /* AMD return status */
    PRINTF (("    status: ")) ;
    if (Info [TRILINOS_AMD_STATUS] == TRILINOS_AMD_OK)
    {
	PRINTF (("OK\n")) ;
    }
    else if (Info [TRILINOS_AMD_STATUS] == TRILINOS_AMD_OUT_OF_MEMORY)
    {
	PRINTF (("out of memory\n")) ;
    }
    else if (Info [TRILINOS_AMD_STATUS] == TRILINOS_AMD_INVALID)
    {
	PRINTF (("invalid matrix\n")) ;
    }
    else if (Info [TRILINOS_AMD_STATUS] == TRILINOS_AMD_OK_BUT_JUMBLED)
    {
	PRINTF (("OK, but jumbled\n")) ;
    }
    else
    {
	PRINTF (("unknown\n")) ;
    }

    /* statistics about the input matrix */
    PRI ("    n, dimension of A:                                  %.20g\n", n);
    PRI ("    nz, number of nonzeros in A:                        %.20g\n",
	Info [TRILINOS_AMD_NZ]) ;
    PRI ("    symmetry of A:                                      %.4f\n",
	Info [TRILINOS_AMD_SYMMETRY]) ;
    PRI ("    number of nonzeros on diagonal:                     %.20g\n",
	Info [TRILINOS_AMD_NZDIAG]) ;
    PRI ("    nonzeros in pattern of A+A' (excl. diagonal):       %.20g\n",
	Info [TRILINOS_AMD_NZ_A_PLUS_AT]) ;
    PRI ("    # dense rows/columns of A+A':                       %.20g\n",
	Info [TRILINOS_AMD_NDENSE]) ;

    /* statistics about AMD's behavior  */
    PRI ("    memory used, in bytes:                              %.20g\n",
	Info [TRILINOS_AMD_MEMORY]) ;
    PRI ("    # of memory compactions:                            %.20g\n",
	Info [TRILINOS_AMD_NCMPA]) ;

    /* statistics about the ordering quality */
    PRINTF (("\n"
	"    The following approximate statistics are for a subsequent\n"
	"    factorization of A(P,P) + A(P,P)'.  They are slight upper\n"
	"    bounds if there are no dense rows/columns in A+A', and become\n"
	"    looser if dense rows/columns exist.\n\n")) ;

    PRI ("    nonzeros in L (excluding diagonal):                 %.20g\n",
	lnz) ;
    PRI ("    nonzeros in L (including diagonal):                 %.20g\n",
	lnzd) ;
    PRI ("    # divide operations for LDL' or LU:                 %.20g\n",
	ndiv) ;
    PRI ("    # multiply-subtract operations for LDL':            %.20g\n",
	nmultsubs_ldl) ;
    PRI ("    # multiply-subtract operations for LU:              %.20g\n",
	nmultsubs_lu) ;
    PRI ("    max nz. in any column of L (incl. diagonal):        %.20g\n",
	Info [TRILINOS_AMD_DMAX]) ;

    /* total flop counts for various factorizations */

    if (n >= 0 && ndiv >= 0 && nmultsubs_ldl >= 0 && nmultsubs_lu >= 0)
    {
	PRINTF (("\n"
	"    chol flop count for real A, sqrt counted as 1 flop: %.20g\n"
	"    LDL' flop count for real A:                         %.20g\n"
	"    LDL' flop count for complex A:                      %.20g\n"
	"    LU flop count for real A (with no pivoting):        %.20g\n"
	"    LU flop count for complex A (with no pivoting):     %.20g\n\n",
	n + ndiv + 2*nmultsubs_ldl,
	    ndiv + 2*nmultsubs_ldl,
	  9*ndiv + 8*nmultsubs_ldl,
	    ndiv + 2*nmultsubs_lu,
	  9*ndiv + 8*nmultsubs_lu)) ;
    }
}
