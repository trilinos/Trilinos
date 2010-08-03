/* ========================================================================== */
/* === KLU_free_symbolic ==================================================== */
/* ========================================================================== */

/* Free the KLU Symbolic object. */

#ifndef KLU2_FREE_SYMBOLIC_HPP
#define KLU2_FREE_SYMBOLIC_HPP

#include "tklu_internal.h"
#include "klu2_memory.hpp"

template <typename Entry, typename Int>
Int KLU_free_symbolic
(
    KLU_symbolic<Entry, Int> **SymbolicHandle,
    KLU_common<Entry, Int>   *Common
)
{
    KLU_symbolic<Entry, Int> *Symbolic ;
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
    KLU_free (Symbolic->P, n, sizeof (Int), Common) ;
    KLU_free (Symbolic->Q, n, sizeof (Int), Common) ;
    KLU_free (Symbolic->R, n+1, sizeof (Int), Common) ;
    KLU_free (Symbolic->Lnz, n, sizeof (double), Common) ; /* TODO: Entry ?? */
    KLU_free (Symbolic, 1, sizeof (KLU_symbolic<Entry, Int>), Common) ;
    *SymbolicHandle = NULL ;
    return (TRUE) ;
}

#endif
