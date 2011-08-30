/* ========================================================================== */
/* === KLU_sort ============================================================= */
/* ========================================================================== */
// @HEADER
// ***********************************************************************
//
//                   KLU2: A Direct Linear Solver package
//                    Copyright 2011 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, with Sandia Corporation, the 
// U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Mike A. Heroux (maherou@sandia.gov)
//
// KLU2 is derived work from KLU, licensed under LGPL, and copyrighted by
// University of Florida. The Authors of KLU are Timothy A. Davis and
// Eka Palamadai. See Doc/KLU_README.txt for the licensing and copyright
// information for KLU.
//
// ***********************************************************************
// @HEADER

/* sorts the columns of L and U so that the row indices appear in strictly
 * increasing order.
 */

#ifndef KLU2_SORT_HPP
#define KLU2_SORT_HPP

#include "klu2_internal.h"
#include "klu2_memory.hpp"

/* ========================================================================== */
/* === sort ================================================================= */
/* ========================================================================== */

/* Sort L or U using a double-transpose */

template <typename Entry, typename Int>
static void sort (Int n, Int *Xip, Int *Xlen, Unit *LU, Int *Tp, Int *Tj,
    Entry *Tx, Int *W)
{
    Int *Xi ;
    Entry *Xx ;
    Int p, i, j, len, nz, tp, xlen, pend ;

    ASSERT (KLU_valid_LU (n, FALSE, Xip, Xlen, LU)) ;

    /* count the number of entries in each row of L or U */ 
    for (i = 0 ; i < n ; i++)
    {
        W [i] = 0 ;
    }
    for (j = 0 ; j < n ; j++)
    {
        GET_POINTER (LU, Xip, Xlen, Xi, Xx, j, len) ;
        for (p = 0 ; p < len ; p++)
        {
            W [Xi [p]]++ ;
        }
    }

    /* construct the row pointers for T */
    nz = 0 ;
    for (i = 0 ; i < n ; i++)
    {
        Tp [i] = nz ;
        nz += W [i] ;
    }
    Tp [n] = nz ;
    for (i = 0 ; i < n ; i++)
    {
        W [i] = Tp [i] ;
    }

    /* transpose the matrix into Tp, Ti, Tx */
    for (j = 0 ; j < n ; j++)
    {
        GET_POINTER (LU, Xip, Xlen, Xi, Xx, j, len) ;
        for (p = 0 ; p < len ; p++)
        {
            tp = W [Xi [p]]++ ;
            Tj [tp] = j ;
            Tx [tp] = Xx [p] ;
        }
    }

    /* transpose the matrix back into Xip, Xlen, Xi, Xx */
    for (j = 0 ; j < n ; j++)
    {
        W [j] = 0 ;
    }
    for (i = 0 ; i < n ; i++)
    {
        pend = Tp [i+1] ;
        for (p = Tp [i] ; p < pend ; p++)
        {
            j = Tj [p] ;
            GET_POINTER (LU, Xip, Xlen, Xi, Xx, j, len) ;
            xlen = W [j]++ ;
            Xi [xlen] = i ;
            Xx [xlen] = Tx [p] ;
        }
    }

    ASSERT (KLU_valid_LU (n, FALSE, Xip, Xlen, LU)) ;
}


/* ========================================================================== */
/* === KLU_sort ============================================================= */
/* ========================================================================== */

template <typename Entry, typename Int>
Int KLU_sort
(
    KLU_symbolic<Entry, Int> *Symbolic,
    KLU_numeric<Entry, Int> *Numeric,
    KLU_common<Entry, Int> *Common
)
{
    Int *R, *W, *Tp, *Ti, *Lip, *Uip, *Llen, *Ulen ;
    Entry *Tx ;
    Unit **LUbx ;
    Int n, nk, nz, block, nblocks, maxblock, k1 ;
    size_t m1 ;

    if (Common == NULL)
    {
        return (FALSE) ;
    }
    Common->status = KLU_OK ;

    n = Symbolic->n ;
    R = Symbolic->R ;
    nblocks = Symbolic->nblocks ;
    maxblock = Symbolic->maxblock ;

    Lip  = Numeric->Lip ;
    Llen = Numeric->Llen ;
    Uip  = Numeric->Uip ;
    Ulen = Numeric->Ulen ;
    LUbx = (Unit **) Numeric->LUbx ;

    m1 = ((size_t) maxblock) + 1 ;

    /* allocate workspace */
    nz = MAX (Numeric->max_lnz_block, Numeric->max_unz_block) ;
    W  = (Int *) KLU_malloc (maxblock, sizeof (Int), Common) ;
    Tp = (Int *) KLU_malloc (m1, sizeof (Int), Common) ;
    Ti = (Int *) KLU_malloc (nz, sizeof (Int), Common) ;
    Tx = (Entry *) KLU_malloc (nz, sizeof (Entry), Common) ;

    PRINTF (("\n======================= Start sort:\n")) ;

    if (Common->status == KLU_OK)
    {
        /* sort each block of L and U */
        for (block = 0 ; block < nblocks ; block++)
        {
            k1 = R [block] ;
            nk = R [block+1] - k1 ;
            if (nk > 1)
            {
                PRINTF (("\n-------------------block: %d nk %d\n", block, nk)) ;
                sort (nk, Lip + k1, Llen + k1, LUbx [block], Tp, Ti, Tx, W) ;
                sort (nk, Uip + k1, Ulen + k1, LUbx [block], Tp, Ti, Tx, W) ;
            }
        }
    }

    PRINTF (("\n======================= sort done.\n")) ;

    /* free workspace */
    KLU_free (W, maxblock, sizeof (Int), Common) ;
    KLU_free (Tp, m1, sizeof (Int), Common) ;
    KLU_free (Ti, nz, sizeof (Int), Common) ;
    KLU_free (Tx, nz, sizeof (Entry), Common) ;
    return (Common->status == KLU_OK) ;
}

#endif
