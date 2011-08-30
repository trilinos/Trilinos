/* ========================================================================== */
/* === KLU_defaults ========================================================= */
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

/* Sets default parameters for KLU */
#ifndef KLU2_DEFAULTS_HPP
#define KLU2_DEFAULTS_HPP

#include "klu2_internal.h"

template <typename Entry, typename Int>
Int KLU_defaults
(
    KLU_common<Entry, Int> *Common
)
{
    if (Common == NULL)
    {
        return (FALSE) ;
    }

    /* parameters */
    Common->tol = 0.001 ;       /* pivot tolerance for diagonal */
    Common->memgrow = 1.2;      /* realloc size ratio increase for LU factors */
    Common->initmem_amd = 1.2 ; /* init. mem with AMD:  c*nnz(L) + n */
    Common->initmem = 10 ;      /* init. mem otherwise: c*nnz(A) + n */
    Common->btf = TRUE ;        /* use BTF pre-ordering, or not */
    Common->maxwork = 0 ;       /* no limit to work done by btf_order */
    Common->ordering = 0 ;      /* 0: AMD, 1: COLAMD, 2: user-provided P and Q,
                                 * 3: user-provided function */
    Common->scale = 2 ;         /* scale: -1: none, and do not check for errors
                                 * in the input matrix in KLU_refactor.
                                 * 0: none, but check for errors,
                                 * 1: sum, 2: max */
    Common->halt_if_singular = TRUE ;   /* quick halt if matrix is singular */

    /* memory management routines */
    Common->malloc_memory  = malloc ;
    Common->calloc_memory  = calloc ;
    Common->free_memory    = free ;
    Common->realloc_memory = realloc ;

    /* user ordering function and optional argument */
    Common->user_order = NULL ;
    Common->user_data = NULL ;

    /* statistics */
    Common->status = KLU_OK ;
    Common->nrealloc = 0 ;
    Common->structural_rank = EMPTY ;
    Common->numerical_rank = EMPTY ;
    Common->noffdiag = EMPTY ;
    Common->flops = EMPTY ;
    Common->rcond = EMPTY ;
    Common->condest = EMPTY ;
    Common->rgrowth = EMPTY ;
    Common->work = 0 ;          /* work done by btf_order */

    Common->memusage = 0 ;
    Common->mempeak = 0 ;

    return (TRUE) ;
}

#endif
