/* ========================================================================== */
/* === TRILINOS_KLU_free_symbolic ==================================================== */
/* ========================================================================== */

/* Free the KLU Symbolic object. */

#include "trilinos_klu_internal.h"

Int TRILINOS_KLU_free_symbolic
(
    TRILINOS_KLU_symbolic **SymbolicHandle,
    TRILINOS_KLU_common	 *Common
)
{
    TRILINOS_KLU_symbolic *Symbolic ;
    Int n ;
    if (Common == NULL)
    {
	return (FALSE) ;
    }
    if (SymbolicHandle == NULL || *SymbolicHandle == NULL)
    {
	return (TRUE) ;
    }
    Symbolic = *SymbolicHandle ;
    n = Symbolic->n ;
    TRILINOS_KLU_free (Symbolic->P, n, sizeof (Int), Common) ;
    TRILINOS_KLU_free (Symbolic->Q, n, sizeof (Int), Common) ;
    TRILINOS_KLU_free (Symbolic->R, n+1, sizeof (Int), Common) ;
    TRILINOS_KLU_free (Symbolic->Lnz, n, sizeof (double), Common) ;
    TRILINOS_KLU_free (Symbolic, 1, sizeof (TRILINOS_KLU_symbolic), Common) ;
    *SymbolicHandle = NULL ;
    return (TRUE) ;
}
