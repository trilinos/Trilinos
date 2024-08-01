/* ========================================================================== */
/* === KLU_free_symbolic ==================================================== */
/* ========================================================================== */
// @HEADER
// *****************************************************************************
//                   KLU2: A Direct Linear Solver package
//
// Copyright 2011 NTESS and the KLU2 contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

/* Free the KLU Symbolic object. */

#ifndef KLU2_FREE_SYMBOLIC_HPP
#define KLU2_FREE_SYMBOLIC_HPP

#include "klu2_internal.h"
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
