// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLUBASKER_MATRIX_VIEW_DECL_HPP
#define SHYLUBASKER_MATRIX_VIEW_DECL_HPP

#define BASKER_DEBUG_MATRIX_VIEW

/*Basker Includes*/
#include "shylubasker_matrix_decl.hpp"
#include "shylubasker_types.hpp"

/*System Includes*/
#include <iostream>
#include <stdio.h>

#ifdef BASKER_KOKKOS
#include <Kokkos_Core.hpp>
#else
#include <omp.h>
#endif

using namespace std;

//--------------------This is all depracated at this time---------//

namespace BaskerNS
{
 
  template <class Int, class Entry, class Exe_Space>
  class BaskerMatrixView
  {
  public:
    BASKER_INLINE
    BaskerMatrixView();
    BASKER_INLINE
    BaskerMatrixView(BASKER_MATRIX *_base, Int _sr, Int _m,
		     Int _sc, Int _n);
    BASKER_INLINE
    void init(BASKER_MATRIX *_base, Int _sr, Int _m, Int _sc, Int _n);
    BASKER_INLINE
    void init_perm(INT_1DARRAY *pinv);
    BASKER_INLINE
    void init_offset(Int k, Int prev);

    /*
    BASKER_INLINE
    void 2D_convert(BASKER_MATRIX &M);
    */

    //NOTE ::  NOT GOING TO USE IN FUTURE
    BASKER_INLINE
    Int good(Int i);
    BASKER_INLINE
    Int good_extend(Int i);
    BASKER_INLINE
    Int col_ptr(Int k);
    BASKER_INLINE
    Int row_idx(Int i);
    BASKER_INLINE
    Entry val(Int i);
    BASKER_INLINE
    Int nnz();
    BASKER_INLINE
    void flip_base();
    BASKER_INLINE
    void flip_base(BASKER_MATRIX *_base);
    BASKER_INLINE
    void print();
    BASKER_INLINE
    void info();

    Int scol;
    Int srow;
    Int ncol;
    Int nrow;

    Int k_offset;

    Int offset;
    Int srow_offset;
    Int m_offset;

    BASKER_BOOL perm;
    BASKER_BOOL roll_back;

    BASKER_MATRIX *base;
    BASKER_MATRIX *base_backup;
    INT_1DARRAY *lpinv;

    //Offsets

  };//end BaskerMatrixView
  
}//end namespace Basker
#endif //endif basker_matrix_view_decl
