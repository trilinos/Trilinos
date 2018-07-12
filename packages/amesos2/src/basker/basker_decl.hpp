// @HEADER
// ***********************************************************************
//
//                   Basker: A Direct Linear Solver package
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
// ***********************************************************************
// @HEADER

#ifndef BASKER_DECL_HPP
#define BASKER_DECL_HPP

#include "basker_types.hpp"

namespace BaskerClassicNS{

  template <class Int, class Entry>
  class BaskerClassic
  {

  public:
    BaskerClassic();
    BaskerClassic(Int nnzL, Int nnzU);
    ~BaskerClassic();
    int preorder(Int *row_perm, Int *col_perm);
    int factor(Int nrow, Int ncol , Int nnz, Int *col_ptr, Int *row_idx, Entry *val);
    int returnL(Int *dim, Int *nnz, Int **col_ptr, Int **row_idx, Entry **val);
    int returnU(Int *dim, Int *nnz, Int **col_ptr, Int **row_idx, Entry **val);
    int returnP(Int **p);
    int solve( Entry* b, Entry* x);
    int solveMultiple(Int nrhs, Entry *b, Entry *x);

    Int get_NnzL();
    Int get_NnzU();
    Int get_NnzLU();
    //int solve();

  private:
    int basker_dfs(
                   Int n,
                   Int j,
                   Int *Li,
                   Int *Lp,
                   Int *color,
                   Int *pattern, /* o/p */
                   Int *top,       /* o/p */
                 
                   Int *tpinv,
                   Int *stack
                   );
    void free_factor();
    void free_perm_matrix();
    int low_tri_solve_csc(Int n, Int* col_ptr, Int *row_idx, Entry *val,  Entry *x,  Entry *b);
    int up_tri_solve_csc(Int n, Int* col_ptr, Int *row_idx, Entry *val, Entry *x, Entry *b);
    int permute_row(Int *p, basker_matrix<Int,Entry> *B);
    int permute_column(Int *p, basker_matrix<Int, Entry> *B);
    int sort_factors();
    Entry* entry_realloc(Entry *old, Int old_size, Int new_size);
    Int* int_realloc(Int *old, Int old_size, Int new_size);
    basker_matrix<Int, Entry> *A;
    basker_matrix<Int, Entry> *L;
    basker_matrix<Int, Entry> *U;
    Int *in_perm;
    Int *pinv;
    Int actual_lnnz;
    Int actual_unnz;
    bool been_fact;
    bool perm_flag;

  };

}/*End namespace*/
#endif
