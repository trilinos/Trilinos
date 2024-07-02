// @HEADER
// *****************************************************************************
//                   Basker: A Direct Linear Solver package
//
// Copyright 2011 NTESS and the Basker contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
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
