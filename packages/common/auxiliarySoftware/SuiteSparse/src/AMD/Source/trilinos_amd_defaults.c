/* ========================================================================= */
/* === TRILINOS_AMD_defaults ======================================================== */
/* ========================================================================= */

/* ------------------------------------------------------------------------- */
/* AMD, Copyright (c) Timothy A. Davis,					     */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README.txt for License.     */
/* email: davis at cise.ufl.edu    CISE Department, Univ. of Florida.        */
/* web: http://www.cise.ufl.edu/research/sparse/amd                          */
/* ------------------------------------------------------------------------- */

/* User-callable.  Sets default control parameters for AMD.  See amd.h
 * for details.
 */

#include "trilinos_amd_internal.h"

/* ========================================================================= */
/* === AMD defaults ======================================================== */
/* ========================================================================= */

GLOBAL void TRILINOS_AMD_defaults
(
    double Control [ ]
)
{
    Int i ;

    if (Control != (double *) NULL)
    {
	for (i = 0 ; i < TRILINOS_AMD_CONTROL ; i++)
	{
	    Control [i] = 0 ;
	}
	Control [TRILINOS_AMD_DENSE] = TRILINOS_AMD_DEFAULT_DENSE ;
	Control [TRILINOS_AMD_AGGRESSIVE] = TRILINOS_AMD_DEFAULT_AGGRESSIVE ;
    }
}
