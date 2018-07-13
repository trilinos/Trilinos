/* ========================================================================= */
/* === TRILINOS_AMD_control ========================================================= */
/* ========================================================================= */

/* ------------------------------------------------------------------------- */
/* AMD, Copyright (c) Timothy A. Davis,					     */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README.txt for License.     */
/* email: davis at cise.ufl.edu    CISE Department, Univ. of Florida.        */
/* web: http://www.cise.ufl.edu/research/sparse/amd                          */
/* ------------------------------------------------------------------------- */

/* User-callable.  Prints the control parameters for AMD.  See amd.h
 * for details.  If the Control array is not present, the defaults are
 * printed instead.
 */

#include "trilinos_amd_internal.h"

GLOBAL void TRILINOS_AMD_control
(
    double Control [ ]
)
{
    double alpha ;
    Int aggressive ;

    if (Control != (double *) NULL)
    {
	alpha = Control [TRILINOS_AMD_DENSE] ;
	aggressive = Control [TRILINOS_AMD_AGGRESSIVE] != 0 ;
    }
    else
    {
	alpha = TRILINOS_AMD_DEFAULT_DENSE ;
	aggressive = TRILINOS_AMD_DEFAULT_AGGRESSIVE ;
    }

    PRINTF (("\nAMD version %d.%d.%d, %s: approximate minimum degree ordering\n"
	"    dense row parameter: %g\n", TRILINOS_AMD_MAIN_VERSION, TRILINOS_AMD_SUB_VERSION,
	TRILINOS_AMD_SUBSUB_VERSION, TRILINOS_AMD_DATE, alpha)) ;

    if (alpha < 0)
    {
	PRINTF (("    no rows treated as dense\n")) ;
    }
    else
    {
	PRINTF ((
	"    (rows with more than max (%g * sqrt (n), 16) entries are\n"
	"    considered \"dense\", and placed last in output permutation)\n",
	alpha)) ;
    }

    if (aggressive)
    {
	PRINTF (("    aggressive absorption:  yes\n")) ;
    }
    else
    {
	PRINTF (("    aggressive absorption:  no\n")) ;
    }

    PRINTF (("    size of AMD integer: %d\n\n", sizeof (Int))) ;
}
