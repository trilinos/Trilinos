/* ========================================================================== */
/* === TRILINOS_KLU_free_numeric ===================================================== */
/* ========================================================================== */

/* Free the KLU Numeric object. */

/* This file should make the long int version of KLU */
#define DLONG 1

#include "trilinos_klu_internal.h"

Int TRILINOS_KLU_free_numeric
(
    TRILINOS_KLU_numeric **NumericHandle,
    TRILINOS_KLU_common	*Common
)
{
    TRILINOS_KLU_numeric *Numeric ;
    Unit **LUbx ;
    size_t *LUsize ;
    Int block, n, nzoff, nblocks ;

    if (Common == NULL)
    {
	return (FALSE) ;
    }
    if (NumericHandle == NULL || *NumericHandle == NULL)
    {
	return (TRUE) ;
    }

    Numeric = *NumericHandle ;

    n = Numeric->n ;
    nzoff = Numeric->nzoff ;
    nblocks = Numeric->nblocks ;
    LUsize = Numeric->LUsize ;

    LUbx = (Unit **) Numeric->LUbx ;
    if (LUbx != NULL)
    {
	for (block = 0 ; block < nblocks ; block++)
	{
	    TRILINOS_KLU_free (LUbx [block], LUsize ? LUsize [block] : 0,
		sizeof (Unit), Common) ;
	}
    }

    TRILINOS_KLU_free (Numeric->Pnum, n, sizeof (Int), Common) ;
    TRILINOS_KLU_free (Numeric->Offp, n+1, sizeof (Int), Common) ;
    TRILINOS_KLU_free (Numeric->Offi, nzoff+1, sizeof (Int), Common) ;
    TRILINOS_KLU_free (Numeric->Offx, nzoff+1, sizeof (Entry), Common) ;

    TRILINOS_KLU_free (Numeric->Lip,  n, sizeof (Int), Common) ;
    TRILINOS_KLU_free (Numeric->Llen, n, sizeof (Int), Common) ;
    TRILINOS_KLU_free (Numeric->Uip,  n, sizeof (Int), Common) ;
    TRILINOS_KLU_free (Numeric->Ulen, n, sizeof (Int), Common) ;

    TRILINOS_KLU_free (Numeric->LUsize, nblocks, sizeof (size_t), Common) ;

    TRILINOS_KLU_free (Numeric->LUbx, nblocks, sizeof (Unit *), Common) ;

    TRILINOS_KLU_free (Numeric->Udiag, n, sizeof (Entry), Common) ;

    TRILINOS_KLU_free (Numeric->Rs,   n, sizeof (double), Common) ;
    TRILINOS_KLU_free (Numeric->Pinv, n, sizeof (Int), Common) ;

    TRILINOS_KLU_free (Numeric->Work, Numeric->worksize, 1, Common) ;

    TRILINOS_KLU_free (Numeric, 1, sizeof (TRILINOS_KLU_numeric), Common) ;

    *NumericHandle = NULL ;
    return (TRUE) ;
}
