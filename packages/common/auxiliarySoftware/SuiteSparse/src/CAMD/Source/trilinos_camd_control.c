/* ========================================================================= */
/* === TRILINOS_CAMD_control ======================================================== */
/* ========================================================================= */

/* ------------------------------------------------------------------------- */
/* CAMD, Copyright (c) Timothy A. Davis, Yanqing Chen,			     */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README.txt for License.     */
/* email: davis at cise.ufl.edu    CISE Department, Univ. of Florida.        */
/* web: http://www.cise.ufl.edu/research/sparse/camd                         */
/* ------------------------------------------------------------------------- */

/* User-callable.  Prints the control parameters for CAMD.  See camd.h
 * for details.  If the Control array is not present, the defaults are
 * printed instead.
 */

#include "trilinos_camd_internal.h"

GLOBAL void TRILINOS_CAMD_control
(
    double Control [ ]
)
{
    double alpha ;
    Int aggressive ;

    if (Control != (double *) NULL)
    {
	alpha = Control [TRILINOS_CAMD_DENSE] ;
	aggressive = Control [TRILINOS_CAMD_AGGRESSIVE] != 0 ;
    }
    else
    {
	alpha = TRILINOS_CAMD_DEFAULT_DENSE ;
	aggressive = TRILINOS_CAMD_DEFAULT_AGGRESSIVE ;
    }

    PRINTF (("\ncamd version %d.%d, %s:  approximate minimum degree ordering:\n"
	"    dense row parameter: %g\n", TRILINOS_CAMD_MAIN_VERSION, TRILINOS_CAMD_SUB_VERSION,
	TRILINOS_CAMD_DATE, alpha)) ;

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
	PRINTF (("    aggressive absorption:  yes\n\n")) ;
    }
    else
    {
	PRINTF (("    aggressive absorption:  no\n\n")) ;
    }
}
