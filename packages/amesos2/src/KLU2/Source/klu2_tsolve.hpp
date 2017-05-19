/* ========================================================================== */
/* === KLU_tsolve =========================================================== */
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
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

/* Solve A'x=b using the symbolic and numeric objects from KLU_analyze
 * (or KLU_analyze_given) and KLU_factor.  Note that no iterative refinement is
 * performed.  Uses Numeric->Xwork as workspace (undefined on input and output),
 * of size 4n Entry's (note that columns 2 to 4 of Xwork overlap with
 * Numeric->Iwork).
 */

#ifndef KLU2_TSOLVE_HPP
#define KLU2_TSOLVE_HPP

#include "klu2_internal.h"
#include "klu2.hpp"

template <typename Entry, typename Int>
Int KLU_tsolve
(
    /* inputs, not modified */
    KLU_symbolic<Entry, Int> *Symbolic,
    KLU_numeric<Entry, Int> *Numeric,
    Int d,                  /* leading dimension of B */
    Int nrhs,               /* number of right-hand-sides */

    /* right-hand-side on input, overwritten with solution to Ax=b on output */
    Entry B [ ],           /* size n*nrhs, in column-oriented form, with
                             * leading dimension d. */
#ifdef COMPLEX
    Int conj_solve,         /* TRUE for conjugate transpose solve, FALSE for
                             * array transpose solve.  Used for the complex
                             * case only. */
#endif
    /* --------------- */
    KLU_common<Entry, Int> *Common
)
{
    Entry x [4], offik, s ;
    double rs, *Rs ;
    Entry *Offx, *X, *Bz, *Udiag ;
    Int *Q, *R, *Pnum, *Offp, *Offi, *Lip, *Uip, *Llen, *Ulen ;
    Unit **LUbx ;
    Int k1, k2, nk, k, block, pend, n, p, nblocks, chunk, nr, i ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    if (Common == NULL)
    {
        return (FALSE) ;
    }
    if (Numeric == NULL || Symbolic == NULL || d < Symbolic->n || nrhs < 0 ||
        B == NULL)
    {
        Common->status = KLU_INVALID ;
        return (FALSE) ;
    }
    Common->status = KLU_OK ;

    /* ---------------------------------------------------------------------- */
    /* get the contents of the Symbolic object */
    /* ---------------------------------------------------------------------- */

    Bz = (Entry *) B ;
    n = Symbolic->n ;
    nblocks = Symbolic->nblocks ;
    Q = Symbolic->Q ;
    R = Symbolic->R ;

    /* ---------------------------------------------------------------------- */
    /* get the contents of the Numeric object */
    /* ---------------------------------------------------------------------- */

    ASSERT (nblocks == Numeric->nblocks) ;
    Pnum = Numeric->Pnum ;
    Offp = Numeric->Offp ;
    Offi = Numeric->Offi ;
    Offx = (Entry *) Numeric->Offx ;

    Lip  = Numeric->Lip ;
    Llen = Numeric->Llen ;
    Uip  = Numeric->Uip ;
    Ulen = Numeric->Ulen ;
    LUbx = (Unit **) Numeric->LUbx ;
    Udiag = (Entry *) Numeric->Udiag ;

    Rs = Numeric->Rs ;
    X = (Entry *) Numeric->Xwork ;
    ASSERT (KLU_valid (n, Offp, Offi, Offx)) ;

    /* ---------------------------------------------------------------------- */
    /* solve in chunks of 4 columns at a time */
    /* ---------------------------------------------------------------------- */

    for (chunk = 0 ; chunk < nrhs ; chunk += 4)
    {

        /* ------------------------------------------------------------------ */
        /* get the size of the current chunk */
        /* ------------------------------------------------------------------ */

        nr = MIN (nrhs - chunk, 4) ;

        /* ------------------------------------------------------------------ */
        /* permute the right hand side, X = Q'*B */
        /* ------------------------------------------------------------------ */

        switch (nr)
        {

            case 1:

                for (k = 0 ; k < n ; k++)
                {
                    X [k] = Bz  [Q [k]] ;
                }
                break ;

            case 2:

                for (k = 0 ; k < n ; k++)
                {
                    i = Q [k] ;
                    X [2*k    ] = Bz [i      ] ;
                    X [2*k + 1] = Bz [i + d  ] ;
                }
                break ;

            case 3:

                for (k = 0 ; k < n ; k++)
                {
                    i = Q [k] ;
                    X [3*k    ] = Bz [i      ] ;
                    X [3*k + 1] = Bz [i + d  ] ;
                    X [3*k + 2] = Bz [i + d*2] ;
                }
                break ;

            case 4:

                for (k = 0 ; k < n ; k++)
                {
                    i = Q [k] ;
                    X [4*k    ] = Bz [i      ] ;
                    X [4*k + 1] = Bz [i + d  ] ;
                    X [4*k + 2] = Bz [i + d*2] ;
                    X [4*k + 3] = Bz [i + d*3] ;
                }
                break ;

        }

        /* ------------------------------------------------------------------ */
        /* solve X = (L*U + Off)'\X */
        /* ------------------------------------------------------------------ */

        for (block = 0 ; block < nblocks ; block++)
        {

            /* -------------------------------------------------------------- */
            /* the block of size nk is from rows/columns k1 to k2-1 */
            /* -------------------------------------------------------------- */

            k1 = R [block] ;
            k2 = R [block+1] ;
            nk = k2 - k1 ;
            PRINTF (("tsolve %d, k1 %d k2-1 %d nk %d\n", block, k1,k2-1,nk)) ;

            /* -------------------------------------------------------------- */
            /* block back-substitution for the off-diagonal-block entries */
            /* -------------------------------------------------------------- */

            if (block > 0)
            {
                switch (nr)
                    {

                    case 1:

                        for (k = k1 ; k < k2 ; k++)
                        {
                            pend = Offp [k+1] ;
                            for (p = Offp [k] ; p < pend ; p++)
                            {
#ifdef COMPLEX
                                if (conj_solve)
                                {
                                    MULT_SUB_CONJ (X [k], X [Offi [p]],
                                            Offx [p]) ;
                                }
                                else
#endif
                                {
                                    MULT_SUB (X [k], Offx [p], X [Offi [p]]) ;
                                }
                            }
                        }
                        break ;

                    case 2:

                        for (k = k1 ; k < k2 ; k++)
                        {
                            pend = Offp [k+1] ;
                            x [0] = X [2*k    ] ;
                            x [1] = X [2*k + 1] ;
                            for (p = Offp [k] ; p < pend ; p++)
                            {
                                i = Offi [p] ;
#ifdef COMPLEX
                                if (conj_solve)
                                {
                                    CONJ (offik, Offx [p]) ;
                                }
                                else
#endif
                                {
                                    offik = Offx [p] ;
                                }
                                MULT_SUB (x [0], offik, X [2*i]) ;
                                MULT_SUB (x [1], offik, X [2*i + 1]) ;
                            }
                            X [2*k    ] = x [0] ;
                            X [2*k + 1] = x [1] ;
                        }
                        break ;

                    case 3:

                        for (k = k1 ; k < k2 ; k++)
                        {
                            pend = Offp [k+1] ;
                            x [0] = X [3*k    ] ;
                            x [1] = X [3*k + 1] ;
                            x [2] = X [3*k + 2] ;
                            for (p = Offp [k] ; p < pend ; p++)
                            {
                                i = Offi [p] ;
#ifdef COMPLEX
                                if (conj_solve)
                                {
                                    CONJ (offik, Offx [p]) ;
                                }
                                else
#endif
                                {
                                    offik = Offx [p] ;
                                }
                                MULT_SUB (x [0], offik, X [3*i]) ;
                                MULT_SUB (x [1], offik, X [3*i + 1]) ;
                                MULT_SUB (x [2], offik, X [3*i + 2]) ;
                            }
                            X [3*k    ] = x [0] ;
                            X [3*k + 1] = x [1] ;
                            X [3*k + 2] = x [2] ;
                        }
                        break ;

                    case 4:

                        for (k = k1 ; k < k2 ; k++)
                        {
                            pend = Offp [k+1] ;
                            x [0] = X [4*k    ] ;
                            x [1] = X [4*k + 1] ;
                            x [2] = X [4*k + 2] ;
                            x [3] = X [4*k + 3] ;
                            for (p = Offp [k] ; p < pend ; p++)
                            {
                                i = Offi [p] ;
#ifdef COMPLEX
                                if (conj_solve)
                                {
                                    CONJ(offik, Offx [p]) ;
                                }
                                else
#endif
                                {
                                    offik = Offx [p] ;
                                }
                                MULT_SUB (x [0], offik, X [4*i]) ;
                                MULT_SUB (x [1], offik, X [4*i + 1]) ;
                                MULT_SUB (x [2], offik, X [4*i + 2]) ;
                                MULT_SUB (x [3], offik, X [4*i + 3]) ;
                            }
                            X [4*k    ] = x [0] ;
                            X [4*k + 1] = x [1] ;
                            X [4*k + 2] = x [2] ;
                            X [4*k + 3] = x [3] ;
                        }
                        break ;
                    }
            }

            /* -------------------------------------------------------------- */
            /* solve the block system */
            /* -------------------------------------------------------------- */

            if (nk == 1)
            {
#ifdef COMPLEX
                if (conj_solve)
                {
                    CONJ (s, Udiag [k1]) ;
                }
                else
#endif
                {
                    s = Udiag [k1] ;
                }
                switch (nr)
                {

                    case 1:
                        DIV (X [k1], X [k1], s) ;
                        break ;

                    case 2:
                        DIV (X [2*k1], X [2*k1], s) ;
                        DIV (X [2*k1 + 1], X [2*k1 + 1], s) ;
                        break ;

                    case 3:
                        DIV (X [3*k1], X [3*k1], s) ;
                        DIV (X [3*k1 + 1], X [3*k1 + 1], s) ;
                        DIV (X [3*k1 + 2], X [3*k1 + 2], s) ;
                        break ;

                    case 4:
                        DIV (X [4*k1], X [4*k1], s) ;
                        DIV (X [4*k1 + 1], X [4*k1 + 1], s) ;
                        DIV (X [4*k1 + 2], X [4*k1 + 2], s) ;
                        DIV (X [4*k1 + 3], X [4*k1 + 3], s) ;
                        break ;

                }
            }
            else
            {
                KLU_utsolve (nk, Uip + k1, Ulen + k1, LUbx [block],
                        Udiag + k1, nr,
#ifdef COMPLEX
                        conj_solve,
#endif
                        X + nr*k1) ;
                KLU_ltsolve (nk, Lip + k1, Llen + k1, LUbx [block], nr,
#ifdef COMPLEX
                        conj_solve,
#endif
                        X + nr*k1) ;
            }
        }

        /* ------------------------------------------------------------------ */
        /* scale and permute the result, Bz  = P'(R\X) */
        /* ------------------------------------------------------------------ */

        if (Rs == NULL)
        {

            /* no scaling */
            switch (nr)
            {

                case 1:

                    for (k = 0 ; k < n ; k++)
                    {
                        Bz  [Pnum [k]] = X [k] ;
                    }
                    break ;

                case 2:

                    for (k = 0 ; k < n ; k++)
                    {
                        i = Pnum [k] ;
                        Bz  [i      ] = X [2*k    ] ;
                        Bz  [i + d  ] = X [2*k + 1] ;
                    }
                    break ;

                case 3:

                    for (k = 0 ; k < n ; k++)
                    {
                        i = Pnum [k] ;
                        Bz  [i      ] = X [3*k    ] ;
                        Bz  [i + d  ] = X [3*k + 1] ;
                        Bz  [i + d*2] = X [3*k + 2] ;
                    }
                    break ;

                case 4:

                    for (k = 0 ; k < n ; k++)
                    {
                        i = Pnum [k] ;
                        Bz  [i      ] = X [4*k    ] ;
                        Bz  [i + d  ] = X [4*k + 1] ;
                        Bz  [i + d*2] = X [4*k + 2] ;
                        Bz  [i + d*3] = X [4*k + 3] ;
                    }
                    break ;
            }

        }
        else
        {

            switch (nr)
            {

                case 1:

                    for (k = 0 ; k < n ; k++)
                    {
                        SCALE_DIV_ASSIGN (Bz [Pnum [k]], X [k], Rs [k]) ;
                    }
                    break ;

                case 2:

                    for (k = 0 ; k < n ; k++)
                    {
                        i = Pnum [k] ;
                        rs = Rs [k] ;
                        SCALE_DIV_ASSIGN (Bz [i], X [2*k], rs) ;
                        SCALE_DIV_ASSIGN (Bz [i + d], X [2*k + 1], rs) ;
                    }
                    break ;

                case 3:

                    for (k = 0 ; k < n ; k++)
                    {
                        i = Pnum [k] ;
                        rs = Rs [k] ;
                        SCALE_DIV_ASSIGN (Bz [i], X [3*k], rs) ;
                        SCALE_DIV_ASSIGN (Bz [i + d], X [3*k + 1], rs) ;
                        SCALE_DIV_ASSIGN (Bz [i + d*2], X [3*k + 2], rs) ;
                    }
                    break ;

                case 4:

                    for (k = 0 ; k < n ; k++)
                    {
                        i = Pnum [k] ;
                        rs = Rs [k] ;
                        SCALE_DIV_ASSIGN (Bz [i], X [4*k], rs) ;
                        SCALE_DIV_ASSIGN (Bz [i + d], X [4*k + 1], rs) ;
                        SCALE_DIV_ASSIGN (Bz [i + d*2], X [4*k + 2], rs) ;
                        SCALE_DIV_ASSIGN (Bz [i + d*3], X [4*k + 3], rs) ;
                    }
                    break ;
            }
        }

        /* ------------------------------------------------------------------ */
        /* go to the next chunk of B */
        /* ------------------------------------------------------------------ */

        Bz  += d*4 ;
    }
    return (TRUE) ;
}


// Require Xout and B to have equal number of entries, pre-allocated
template <typename Entry, typename Int>
Int KLU_tsolve2
(
    /* inputs, not modified */
    KLU_symbolic<Entry, Int> *Symbolic,
    KLU_numeric<Entry, Int> *Numeric,
    Int d,                  /* leading dimension of B */
    Int nrhs,               /* number of right-hand-sides */

    /* right-hand-side on input, overwritten with solution to Ax=b on output */
    Entry B [ ],            /* size n*1, in column-oriented form, with
                             * leading dimension d. */
    Entry Xout [ ],         /* size n*1, in column-oriented form, with
                             * leading dimension d. */
#ifdef COMPLEX
    Int conj_solve,         /* TRUE for conjugate transpose solve, FALSE for
                             * array transpose solve.  Used for the complex
                             * case only. */
#endif
    /* --------------- */
    KLU_common<Entry, Int> *Common
)
{
    Entry x [4], offik, s ;
    double rs, *Rs ;
    Entry *Offx, *X, *Bz, *Udiag ;
    Int *Q, *R, *Pnum, *Offp, *Offi, *Lip, *Uip, *Llen, *Ulen ;
    Unit **LUbx ;
    Int k1, k2, nk, k, block, pend, n, p, nblocks, chunk, nr, i ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    if (Common == NULL)
    {
        return (FALSE) ;
    }
    if (Numeric == NULL || Symbolic == NULL || d < Symbolic->n || nrhs < 0 ||
        B == NULL || Xout == NULL)
    {
        Common->status = KLU_INVALID ;
        return (FALSE) ;
    }
    Common->status = KLU_OK ;

    /* ---------------------------------------------------------------------- */
    /* get the contents of the Symbolic object */
    /* ---------------------------------------------------------------------- */

    Bz = (Entry *) B ;
    n = Symbolic->n ;
    nblocks = Symbolic->nblocks ;
    Q = Symbolic->Q ;
    R = Symbolic->R ;

    /* ---------------------------------------------------------------------- */
    /* get the contents of the Numeric object */
    /* ---------------------------------------------------------------------- */

    ASSERT (nblocks == Numeric->nblocks) ;
    Pnum = Numeric->Pnum ;
    Offp = Numeric->Offp ;
    Offi = Numeric->Offi ;
    Offx = (Entry *) Numeric->Offx ;

    Lip  = Numeric->Lip ;
    Llen = Numeric->Llen ;
    Uip  = Numeric->Uip ;
    Ulen = Numeric->Ulen ;
    LUbx = (Unit **) Numeric->LUbx ;
    Udiag = (Entry *) Numeric->Udiag ;

    Rs = Numeric->Rs ;
    X = (Entry *) Numeric->Xwork ;
    ASSERT (KLU_valid (n, Offp, Offi, Offx)) ;

    /* ---------------------------------------------------------------------- */
    /* solve in chunks of 4 columns at a time */
    /* ---------------------------------------------------------------------- */

    for (chunk = 0 ; chunk < nrhs ; chunk += 4)
    {

        /* ------------------------------------------------------------------ */
        /* get the size of the current chunk */
        /* ------------------------------------------------------------------ */

        nr = MIN (nrhs - chunk, 4) ;

        /* ------------------------------------------------------------------ */
        /* permute the right hand side, X = Q'*B */
        /* ------------------------------------------------------------------ */

        switch (nr)
        {

            case 1:

                for (k = 0 ; k < n ; k++)
                {
                    X [k] = Bz  [Q [k]] ;
                }
                break ;

            case 2:

                for (k = 0 ; k < n ; k++)
                {
                    i = Q [k] ;
                    X [2*k    ] = Bz [i      ] ;
                    X [2*k + 1] = Bz [i + d  ] ;
                }
                break ;

            case 3:

                for (k = 0 ; k < n ; k++)
                {
                    i = Q [k] ;
                    X [3*k    ] = Bz [i      ] ;
                    X [3*k + 1] = Bz [i + d  ] ;
                    X [3*k + 2] = Bz [i + d*2] ;
                }
                break ;

            case 4:

                for (k = 0 ; k < n ; k++)
                {
                    i = Q [k] ;
                    X [4*k    ] = Bz [i      ] ;
                    X [4*k + 1] = Bz [i + d  ] ;
                    X [4*k + 2] = Bz [i + d*2] ;
                    X [4*k + 3] = Bz [i + d*3] ;
                }
                break ;

        }

        /* ------------------------------------------------------------------ */
        /* solve X = (L*U + Off)'\X */
        /* ------------------------------------------------------------------ */

        for (block = 0 ; block < nblocks ; block++)
        {

            /* -------------------------------------------------------------- */
            /* the block of size nk is from rows/columns k1 to k2-1 */
            /* -------------------------------------------------------------- */

            k1 = R [block] ;
            k2 = R [block+1] ;
            nk = k2 - k1 ;
            PRINTF (("tsolve %d, k1 %d k2-1 %d nk %d\n", block, k1,k2-1,nk)) ;

            /* -------------------------------------------------------------- */
            /* block back-substitution for the off-diagonal-block entries */
            /* -------------------------------------------------------------- */

            if (block > 0)
            {
                switch (nr)
                    {

                    case 1:

                        for (k = k1 ; k < k2 ; k++)
                        {
                            pend = Offp [k+1] ;
                            for (p = Offp [k] ; p < pend ; p++)
                            {
#ifdef COMPLEX
                                if (conj_solve)
                                {
                                    MULT_SUB_CONJ (X [k], X [Offi [p]],
                                            Offx [p]) ;
                                }
                                else
#endif
                                {
                                    MULT_SUB (X [k], Offx [p], X [Offi [p]]) ;
                                }
                            }
                        }
                        break ;

                    case 2:

                        for (k = k1 ; k < k2 ; k++)
                        {
                            pend = Offp [k+1] ;
                            x [0] = X [2*k    ] ;
                            x [1] = X [2*k + 1] ;
                            for (p = Offp [k] ; p < pend ; p++)
                            {
                                i = Offi [p] ;
#ifdef COMPLEX
                                if (conj_solve)
                                {
                                    CONJ (offik, Offx [p]) ;
                                }
                                else
#endif
                                {
                                    offik = Offx [p] ;
                                }
                                MULT_SUB (x [0], offik, X [2*i]) ;
                                MULT_SUB (x [1], offik, X [2*i + 1]) ;
                            }
                            X [2*k    ] = x [0] ;
                            X [2*k + 1] = x [1] ;
                        }
                        break ;

                    case 3:

                        for (k = k1 ; k < k2 ; k++)
                        {
                            pend = Offp [k+1] ;
                            x [0] = X [3*k    ] ;
                            x [1] = X [3*k + 1] ;
                            x [2] = X [3*k + 2] ;
                            for (p = Offp [k] ; p < pend ; p++)
                            {
                                i = Offi [p] ;
#ifdef COMPLEX
                                if (conj_solve)
                                {
                                    CONJ (offik, Offx [p]) ;
                                }
                                else
#endif
                                {
                                    offik = Offx [p] ;
                                }
                                MULT_SUB (x [0], offik, X [3*i]) ;
                                MULT_SUB (x [1], offik, X [3*i + 1]) ;
                                MULT_SUB (x [2], offik, X [3*i + 2]) ;
                            }
                            X [3*k    ] = x [0] ;
                            X [3*k + 1] = x [1] ;
                            X [3*k + 2] = x [2] ;
                        }
                        break ;

                    case 4:

                        for (k = k1 ; k < k2 ; k++)
                        {
                            pend = Offp [k+1] ;
                            x [0] = X [4*k    ] ;
                            x [1] = X [4*k + 1] ;
                            x [2] = X [4*k + 2] ;
                            x [3] = X [4*k + 3] ;
                            for (p = Offp [k] ; p < pend ; p++)
                            {
                                i = Offi [p] ;
#ifdef COMPLEX
                                if (conj_solve)
                                {
                                    CONJ(offik, Offx [p]) ;
                                }
                                else
#endif
                                {
                                    offik = Offx [p] ;
                                }
                                MULT_SUB (x [0], offik, X [4*i]) ;
                                MULT_SUB (x [1], offik, X [4*i + 1]) ;
                                MULT_SUB (x [2], offik, X [4*i + 2]) ;
                                MULT_SUB (x [3], offik, X [4*i + 3]) ;
                            }
                            X [4*k    ] = x [0] ;
                            X [4*k + 1] = x [1] ;
                            X [4*k + 2] = x [2] ;
                            X [4*k + 3] = x [3] ;
                        }
                        break ;
                    }
            }

            /* -------------------------------------------------------------- */
            /* solve the block system */
            /* -------------------------------------------------------------- */

            if (nk == 1)
            {
#ifdef COMPLEX
                if (conj_solve)
                {
                    CONJ (s, Udiag [k1]) ;
                }
                else
#endif
                {
                    s = Udiag [k1] ;
                }
                switch (nr)
                {

                    case 1:
                        DIV (X [k1], X [k1], s) ;
                        break ;

                    case 2:
                        DIV (X [2*k1], X [2*k1], s) ;
                        DIV (X [2*k1 + 1], X [2*k1 + 1], s) ;
                        break ;

                    case 3:
                        DIV (X [3*k1], X [3*k1], s) ;
                        DIV (X [3*k1 + 1], X [3*k1 + 1], s) ;
                        DIV (X [3*k1 + 2], X [3*k1 + 2], s) ;
                        break ;

                    case 4:
                        DIV (X [4*k1], X [4*k1], s) ;
                        DIV (X [4*k1 + 1], X [4*k1 + 1], s) ;
                        DIV (X [4*k1 + 2], X [4*k1 + 2], s) ;
                        DIV (X [4*k1 + 3], X [4*k1 + 3], s) ;
                        break ;

                }
            }
            else
            {
                KLU_utsolve (nk, Uip + k1, Ulen + k1, LUbx [block],
                        Udiag + k1, nr,
#ifdef COMPLEX
                        conj_solve,
#endif
                        X + nr*k1) ;
                KLU_ltsolve (nk, Lip + k1, Llen + k1, LUbx [block], nr,
#ifdef COMPLEX
                        conj_solve,
#endif
                        X + nr*k1) ;
            }
        }

        /* ------------------------------------------------------------------ */
        /* scale and permute the result, Bz  = P'(R\X) */
        /* store result in Xout */
        /* ------------------------------------------------------------------ */

        if (Rs == NULL)
        {

            /* no scaling */
            switch (nr)
            {

                case 1:

                    for (k = 0 ; k < n ; k++)
                    {
                        Xout  [Pnum [k]] = X [k] ;
                    }
                    break ;

                case 2:

                    for (k = 0 ; k < n ; k++)
                    {
                        i = Pnum [k] ;
                        Xout  [i      ] = X [2*k    ] ;
                        Xout  [i + d  ] = X [2*k + 1] ;
                    }
                    break ;

                case 3:

                    for (k = 0 ; k < n ; k++)
                    {
                        i = Pnum [k] ;
                        Xout  [i      ] = X [3*k    ] ;
                        Xout  [i + d  ] = X [3*k + 1] ;
                        Xout  [i + d*2] = X [3*k + 2] ;
                    }
                    break ;

                case 4:

                    for (k = 0 ; k < n ; k++)
                    {
                        i = Pnum [k] ;
                        Xout  [i      ] = X [4*k    ] ;
                        Xout  [i + d  ] = X [4*k + 1] ;
                        Xout  [i + d*2] = X [4*k + 2] ;
                        Xout  [i + d*3] = X [4*k + 3] ;
                    }
                    break ;
            }

        }
        else
        {

            switch (nr)
            {

                case 1:

                    for (k = 0 ; k < n ; k++)
                    {
                        SCALE_DIV_ASSIGN (Xout [Pnum [k]], X [k], Rs [k]) ;
                    }
                    break ;

                case 2:

                    for (k = 0 ; k < n ; k++)
                    {
                        i = Pnum [k] ;
                        rs = Rs [k] ;
                        SCALE_DIV_ASSIGN (Xout [i], X [2*k], rs) ;
                        SCALE_DIV_ASSIGN (Xout [i + d], X [2*k + 1], rs) ;
                    }
                    break ;

                case 3:

                    for (k = 0 ; k < n ; k++)
                    {
                        i = Pnum [k] ;
                        rs = Rs [k] ;
                        SCALE_DIV_ASSIGN (Xout [i], X [3*k], rs) ;
                        SCALE_DIV_ASSIGN (Xout [i + d], X [3*k + 1], rs) ;
                        SCALE_DIV_ASSIGN (Xout [i + d*2], X [3*k + 2], rs) ;
                    }
                    break ;

                case 4:

                    for (k = 0 ; k < n ; k++)
                    {
                        i = Pnum [k] ;
                        rs = Rs [k] ;
                        SCALE_DIV_ASSIGN (Xout [i], X [4*k], rs) ;
                        SCALE_DIV_ASSIGN (Xout [i + d], X [4*k + 1], rs) ;
                        SCALE_DIV_ASSIGN (Xout [i + d*2], X [4*k + 2], rs) ;
                        SCALE_DIV_ASSIGN (Xout [i + d*3], X [4*k + 3], rs) ;
                    }
                    break ;
            }
        }

        /* ------------------------------------------------------------------ */
        /* go to the next chunk of B */
        /* ------------------------------------------------------------------ */

        Xout  += d*4 ;
    }
    return (TRUE) ;
}

#endif
