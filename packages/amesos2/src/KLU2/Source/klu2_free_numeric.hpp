/* ========================================================================== */
/* === KLU_free_numeric ===================================================== */
/* ========================================================================== */
// @HEADER
// *****************************************************************************
//                   KLU2: A Direct Linear Solver package
//
// Copyright 2011 NTESS and the KLU2 contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

/* Free the KLU Numeric object. */

#ifndef KLU2_FREE_NUMERIC_HPP
#define KLU2_FREE_NUMERIC_HPP

#include "klu2_internal.h"
#include "klu2_memory.hpp"

template <typename Entry, typename Int>
Int KLU_free_numeric
(
    KLU_numeric<Entry, Int> **NumericHandle,
    KLU_common<Entry, Int>  *Common
)
{
    KLU_numeric<Entry, Int> *Numeric ;
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
            KLU_free (LUbx [block], LUsize ? LUsize [block] : 0,
                sizeof (Unit), Common) ;
        }
    }

    KLU_free (Numeric->Pnum, n, sizeof (Int), Common) ;
    KLU_free (Numeric->Offp, n+1, sizeof (Int), Common) ;
    KLU_free (Numeric->Offi, nzoff+1, sizeof (Int), Common) ;
    KLU_free (Numeric->Offx, nzoff+1, sizeof (Entry), Common) ;

    KLU_free (Numeric->Lip,  n, sizeof (Int), Common) ;
    KLU_free (Numeric->Llen, n, sizeof (Int), Common) ;
    KLU_free (Numeric->Uip,  n, sizeof (Int), Common) ;
    KLU_free (Numeric->Ulen, n, sizeof (Int), Common) ;

    KLU_free (Numeric->LUsize, nblocks, sizeof (size_t), Common) ;

    KLU_free (Numeric->LUbx, nblocks, sizeof (Unit *), Common) ;

    KLU_free (Numeric->Udiag, n, sizeof (Entry), Common) ;

    KLU_free (Numeric->Rs,   n, sizeof (double), Common) ;
    KLU_free (Numeric->Pinv, n, sizeof (Int), Common) ;

    KLU_free (Numeric->Work, Numeric->worksize, 1, Common) ;

    KLU_free (Numeric, 1, sizeof (KLU_numeric<Entry, Int>), Common) ;

    *NumericHandle = NULL ;
    return (TRUE) ;
}

#endif
