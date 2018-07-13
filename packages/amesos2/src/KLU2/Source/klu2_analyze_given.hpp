/* ========================================================================== */
/* === klu_analyze_given ==================================================== */
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

/* Given an input permutation P and Q, create the Symbolic object.  BTF can
 * be done to modify the user's P and Q (does not perform the max transversal;
 * just finds the strongly-connected components). */

#ifndef KLU2_ANALYZE_GIVEN_HPP
#define KLU2_ANALYZE_GIVEN_HPP

#include "klu2_internal.h"
#include "klu2_memory.hpp"

/* ========================================================================== */
/* === klu_alloc_symbolic =================================================== */
/* ========================================================================== */

/* Allocate Symbolic object, and check input matrix.  Not user callable. */

template <typename Entry, typename Int>
KLU_symbolic<Entry, Int> *KLU_alloc_symbolic
(
    Int n,
    Int *Ap,
    Int *Ai,
    KLU_common<Entry, Int> *Common
)
{
    KLU_symbolic<Entry, Int> *Symbolic ;
    Int *P, *Q, *R ;
    double *Lnz ;
    Int nz, i, j, p, pend ;

    if (Common == NULL)
    {
        return (NULL) ;
    }
    Common->status = KLU_OK ;

    /* A is n-by-n, with n > 0.  Ap [0] = 0 and nz = Ap [n] >= 0 required.
     * Ap [j] <= Ap [j+1] must hold for all j = 0 to n-1.  Row indices in Ai
     * must be in the range 0 to n-1, and no duplicate entries can be present.
     * The list of row indices in each column of A need not be sorted.
     */

    if (n <= 0 || Ap == NULL || Ai == NULL)
    {
        /* Ap and Ai must be present, and n must be > 0 */
        Common->status = KLU_INVALID ;
        return (NULL) ;
    }

    nz = Ap [n] ;
    if (Ap [0] != 0 || nz < 0)
    {
        /* nz must be >= 0 and Ap [0] must equal zero */
        Common->status = KLU_INVALID ;
        return (NULL) ;
    }

    for (j = 0 ; j < n ; j++)
    {
        if (Ap [j] > Ap [j+1])
        {
            /* column pointers must be non-decreasing */
            Common->status = KLU_INVALID ;
            return (NULL) ;
        }
    }
    P = (Int *) KLU_malloc (n, sizeof (Int), Common) ;
    if (Common->status < KLU_OK)
    {
        /* out of memory */
        Common->status = KLU_OUT_OF_MEMORY ;
        return (NULL) ;
    }
    for (i = 0 ; i < n ; i++)
    {
        P [i] = EMPTY ;
    }
    for (j = 0 ; j < n ; j++)
    {
        pend = Ap [j+1] ;
        for (p = Ap [j] ; p < pend ; p++)
        {
            i = Ai [p] ;
            if (i < 0 || i >= n || P [i] == j)
            {
                /* row index out of range, or duplicate entry */
                KLU_free (P, n, sizeof (Int), Common) ;
                Common->status = KLU_INVALID ;
                return (NULL) ;
            }
            /* flag row i as appearing in column j */
            P [i] = j ;
        }
    }

    /* ---------------------------------------------------------------------- */
    /* allocate the Symbolic object */
    /* ---------------------------------------------------------------------- */

    Symbolic = (KLU_symbolic<Entry, Int> *) KLU_malloc (sizeof (KLU_symbolic<Entry, Int>), 1, Common) ;
    if (Common->status < KLU_OK)
    {
        /* out of memory */
        KLU_free (P, n, sizeof (Int), Common) ;
        Common->status = KLU_OUT_OF_MEMORY ;
        return (NULL) ;
    }

    Q = (Int *) KLU_malloc (n, sizeof (Int), Common) ;
    R = (Int *) KLU_malloc (n+1, sizeof (Int), Common) ;
    Lnz = (double *) KLU_malloc (n, sizeof (double), Common) ;

    Symbolic->n = n ;
    Symbolic->nz = nz ;
    Symbolic->P = P ;
    Symbolic->Q = Q ;
    Symbolic->R = R ;
    Symbolic->Lnz = Lnz ;

    if (Common->status < KLU_OK)
    {
        /* out of memory */
        KLU_free_symbolic (&Symbolic, Common) ;
        Common->status = KLU_OUT_OF_MEMORY ;
        return (NULL) ;
    }

    return (Symbolic) ;
}


/* ========================================================================== */
/* === KLU_analyze_given ==================================================== */
/* ========================================================================== */

template <typename Entry, typename Int>
KLU_symbolic<Entry, Int> *KLU_analyze_given     /* returns NULL if error, or a valid
                                       KLU_symbolic object if successful */
(
    /* inputs, not modified */
    Int n,              /* A is n-by-n */
    Int Ap [ ],         /* size n+1, column pointers */
    Int Ai [ ],         /* size nz, row indices */
    Int Puser [ ],      /* size n, user's row permutation (may be NULL) */
    Int Quser [ ],      /* size n, user's column permutation (may be NULL) */
    /* -------------------- */
    KLU_common<Entry, Int> *Common
)
{
    KLU_symbolic<Entry, Int> *Symbolic ;
    double *Lnz ;
    Int nblocks, nz, block, maxblock, *P, *Q, *R, nzoff, p, pend, do_btf, k ;

    /* ---------------------------------------------------------------------- */
    /* determine if input matrix is valid, and get # of nonzeros */
    /* ---------------------------------------------------------------------- */

    Symbolic = KLU_alloc_symbolic (n, Ap, Ai, Common) ;
    if (Symbolic == NULL)
    {
        return (NULL) ;
    }
    P = Symbolic->P ;
    Q = Symbolic->Q ;
    R = Symbolic->R ;
    Lnz = Symbolic->Lnz ;
    nz = Symbolic->nz ;

    /* ---------------------------------------------------------------------- */
    /* Q = Quser, or identity if Quser is NULL */
    /* ---------------------------------------------------------------------- */

    if (Quser == (Int *) NULL)
    {
        for (k = 0 ; k < n ; k++)
        {
            Q [k] = k ;
        }
    }
    else
    {
        for (k = 0 ; k < n ; k++)
        {
            Q [k] = Quser [k] ;
        }
    }

    /* ---------------------------------------------------------------------- */
    /* get the control parameters for BTF and ordering method */
    /* ---------------------------------------------------------------------- */

    do_btf = Common->btf ;
    do_btf = (do_btf) ? TRUE : FALSE ;
    Symbolic->ordering = 2 ;
    Symbolic->do_btf = do_btf ;

    /* ---------------------------------------------------------------------- */
    /* find the block triangular form, if requested */
    /* ---------------------------------------------------------------------- */

    if (do_btf)
    {

        /* ------------------------------------------------------------------ */
        /* get workspace for BTF_strongcomp */
        /* ------------------------------------------------------------------ */

        Int *Pinv, *Work, *Bi, k1, k2, nk, oldcol ;

        Work = (Int *) KLU_malloc (4*n, sizeof (Int), Common) ;
        Pinv = (Int *) KLU_malloc (n, sizeof (Int), Common) ;
        if (Puser != (Int *) NULL)
        {
            Bi = (Int *) KLU_malloc (nz+1, sizeof (Int), Common) ;
        }
        else
        {
            Bi = Ai ;
        }

        if (Common->status < KLU_OK)
        {
            /* out of memory */
            KLU_free (Work, 4*n, sizeof (Int), Common) ;
            KLU_free (Pinv, n, sizeof (Int), Common) ;
            if (Puser != (Int *) NULL)
            {
                KLU_free (Bi, nz+1, sizeof (Int), Common) ;
            }
            KLU_free_symbolic (&Symbolic, Common) ;
            Common->status = KLU_OUT_OF_MEMORY ;
            return (NULL) ;
        }

        /* ------------------------------------------------------------------ */
        /* B = Puser * A */
        /* ------------------------------------------------------------------ */

        if (Puser != (Int *) NULL)
        {
            for (k = 0 ; k < n ; k++)
            {
                Pinv [Puser [k]] = k ;
            }
            for (p = 0 ; p < nz ; p++)
            {
                Bi [p] = Pinv [Ai [p]] ;
            }
        }

        /* ------------------------------------------------------------------ */
        /* find the strongly-connected components */
        /* ------------------------------------------------------------------ */

        /* TODO : Correct version of BTF */
        /* modifies Q, and determines P and R */
        /*nblocks = BTF_strongcomp (n, Ap, Bi, Q, P, R, Work) ;*/
        nblocks = KLU_OrdinalTraits<Int>::btf_strongcomp (n, Ap, Bi, Q, P, R, 
                    Work) ;

        /* ------------------------------------------------------------------ */
        /* P = P * Puser */
        /* ------------------------------------------------------------------ */

        if (Puser != (Int *) NULL)
        {
            for (k = 0 ; k < n ; k++)
            {
                Work [k] = Puser [P [k]] ;
            }
            for (k = 0 ; k < n ; k++)
            {
                P [k] = Work [k] ;
            }
        }

        /* ------------------------------------------------------------------ */
        /* Pinv = inverse of P */
        /* ------------------------------------------------------------------ */

        for (k = 0 ; k < n ; k++)
        {
            Pinv [P [k]] = k ;
        }

        /* ------------------------------------------------------------------ */
        /* analyze each block */
        /* ------------------------------------------------------------------ */

        nzoff = 0 ;         /* nz in off-diagonal part */
        maxblock = 1 ;      /* size of the largest block */

        for (block = 0 ; block < nblocks ; block++)
        {

            /* -------------------------------------------------------------- */
            /* the block is from rows/columns k1 to k2-1 */
            /* -------------------------------------------------------------- */

            k1 = R [block] ;
            k2 = R [block+1] ;
            nk = k2 - k1 ;
            PRINTF (("BLOCK %d, k1 %d k2-1 %d nk %d\n", block, k1, k2-1, nk)) ;
            maxblock = MAX (maxblock, nk) ;

            /* -------------------------------------------------------------- */
            /* scan the kth block, C */
            /* -------------------------------------------------------------- */

            for (k = k1 ; k < k2 ; k++)
            {
                oldcol = Q [k] ;
                pend = Ap [oldcol+1] ;
                for (p = Ap [oldcol] ; p < pend ; p++)
                {
                    if (Pinv [Ai [p]] < k1)
                    {
                        nzoff++ ;
                    }
                }
            }

            /* fill-in not estimated */
            Lnz [block] = EMPTY ;
        }

        /* ------------------------------------------------------------------ */
        /* free all workspace */
        /* ------------------------------------------------------------------ */

        KLU_free (Work, 4*n, sizeof (Int), Common) ;
        KLU_free (Pinv, n, sizeof (Int), Common) ;
        if (Puser != (Int *) NULL)
        {
            KLU_free (Bi, nz+1, sizeof (Int), Common) ;
        }

    }
    else
    {

        /* ------------------------------------------------------------------ */
        /* BTF not requested */
        /* ------------------------------------------------------------------ */

        nzoff = 0 ;
        nblocks = 1 ;
        maxblock = n ;
        R [0] = 0 ;
        R [1] = n ;
        Lnz [0] = EMPTY ;

        /* ------------------------------------------------------------------ */
        /* P = Puser, or identity if Puser is NULL */
        /* ------------------------------------------------------------------ */

        for (k = 0 ; k < n ; k++)
        {
            P [k] = (Puser == NULL) ? k : Puser [k] ;
        }
    }

    /* ---------------------------------------------------------------------- */
    /* return the symbolic object */
    /* ---------------------------------------------------------------------- */

    Symbolic->nblocks = nblocks ;
    Symbolic->maxblock = maxblock ;
    Symbolic->lnz = EMPTY ;
    Symbolic->unz = EMPTY ;
    Symbolic->nzoff = nzoff ;

    return (Symbolic) ;
}

#endif /* KLU2_ANALYZE_GIVEN_HPP */
