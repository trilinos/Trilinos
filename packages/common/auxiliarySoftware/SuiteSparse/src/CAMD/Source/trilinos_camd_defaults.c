/* ========================================================================= */
/* === TRILINOS_CAMD_defaults ======================================================= */
/* ========================================================================= */

/* ------------------------------------------------------------------------- */
/* CAMD, Copyright (c) Timothy A. Davis, Yanqing Chen,			     */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README.txt for License.     */
/* email: davis at cise.ufl.edu    CISE Department, Univ. of Florida.        */
/* web: http://www.cise.ufl.edu/research/sparse/camd                         */
/* ------------------------------------------------------------------------- */

/* User-callable.  Sets default control parameters for CAMD.  See camd.h
 * for details.
 */

#include "trilinos_camd_internal.h"

/* ========================================================================= */
/* === CAMD defaults ======================================================= */
/* ========================================================================= */

GLOBAL void TRILINOS_CAMD_defaults
(
    double Control [ ]
)
{
    Int i ;
    if (Control != (double *) NULL)
    {
	for (i = 0 ; i < TRILINOS_CAMD_CONTROL ; i++)
	{
	    Control [i] = 0 ;
	}
	Control [TRILINOS_CAMD_DENSE] = TRILINOS_CAMD_DEFAULT_DENSE ;
	Control [TRILINOS_CAMD_AGGRESSIVE] = TRILINOS_CAMD_DEFAULT_AGGRESSIVE ;
    }
}
