// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLUBASKER_DECL_HPP
#define SHYLUBASKER_DECL_HPP

//Basker Includes
#include "ShyLU_NodeBasker_config.h"
#include "shylubasker_matrix_decl.hpp"
#include "shylubasker_matrix_view_decl.hpp"
#include "shylubasker_stats.hpp"
#include "shylubasker_types.hpp"
#include "shylubasker_structs.hpp"
#include "shylubasker_thread.hpp"
#include "shylubasker_scalar_traits.hpp"

//Kokkos Includes
#ifdef BASKER_KOKKOS
#include "Kokkos_Core.hpp"
#else
#include <omp.h>
#endif

//System Includes
#include <iostream>
#include <stdio.h>

namespace BaskerNS
{
  template <class Int, class Entry, class Exe_Space>
  class Basker
  {
  public:

    #ifdef BASKER_KOKKOS
    typedef Kokkos::TeamPolicy<Exe_Space>    TeamPolicy;
    typedef typename TeamPolicy::member_type TeamMember;
    #endif

    typedef Basker_ScalarTraits<Entry>        EntryOP;

    //Options
    basker_options<Int,Entry,Exe_Space> Options;

    //basker_def.hpp
    BASKER_INLINE
    Basker();

    BASKER_INLINE
    ~Basker();

    BASKER_INLINE
    int InitMatrix(string filename);

    BASKER_INLINE
    int InitMatrix(Int nrow, Int ncol, Int nnz, Int *col_ptr, Int *row_idx, Entry *val);

    BASKER_INLINE
    int Order(Int option); //Future replacement for InitOrder

    BASKER_INLINE
    int InitOrder(Int option);

    BASKER_INLINE
    int InitOrder(Int *perm, Int nblks, Int parts, Int *row_tabs, Int *col_tabs, Int *tree_tabs);

    BASKER_INLINE
    int Symbolic(Int option);

    BASKER_INLINE
    int Symbolic(Int nrow, Int ncol, Int nnz, Int *col_ptr, Int *row_idx, Entry *val, bool transpose_needed = false);

    BASKER_INLINE
    int Factor(Int option);

    BASKER_INLINE
    int Factor(Int nrow, Int ncol, Int nnz, Int *col_ptr, Int *row_idx, Entry *val);

    BASKER_INLINE
    int Factor_Inc(Int option);

    BASKER_INLINE
    int Solve(Entry *b, Entry *x, bool transpose = false);

    BASKER_INLINE
    int Solve(Int nrhs, Entry *b, Entry *x, bool transpose = false);

    BASKER_INLINE
    int Solve(ENTRY_1DARRAY b, ENTRY_1DARRAY x, bool transpose = false);

    BASKER_INLINE
    int Solve(Int nrhs, Entry *b, Entry *x, Int option, bool transpose = false);

    BASKER_INLINE
    int SetThreads(Int nthreads);

    BASKER_INLINE
    int GetLnnz(Int &Lnnz);

    BASKER_INLINE
    int GetUnnz(Int &Unnz);

    BASKER_INLINE
    int GetL(Int &n, Int &nnz, Int **col_ptr, Int **row_idx, Entry **val);

    BASKER_INLINE
    int GetU(Int &n, Int &nnz, Int **col_ptr, Int **row_idx, Entry **val);

    BASKER_INLINE
    int GetPerm(Int **lp, Int **rp);

    BASKER_INLINE
    void PrintTime();
    
    BASKER_INLINE
    int Info();

    BASKER_INLINE
    int KokkosPlay();

    BASKER_INLINE
    void PRINT_C();

    BASKER_INLINE
    void DEBUG_PRINT();
    
    BASKER_INLINE
    void Finalize();

 
    BASKER_INLINE
    int t_nfactor_blk(Int kid);

    BASKER_INLINE
    int t_nfactor_blk_inc_lvl(Int kid);

    int t_nfactor_blk_old(Int kid);

    BASKER_INLINE
    void t_init_workspace(bool flag, Int kid);

    BASKER_INLINE
    void t_init_2DA(Int kid,
                    BASKER_BOOL _alloc = BASKER_FALSE,
                    BASKER_BOOL keep_zeros = BASKER_TRUE);

    BASKER_INLINE
    #ifdef BASKER_KOKKOS
    int t_nfactor_sep(Int kid, Int lvl, Int team_leader, const TeamMember &thread);
    #else
    int t_nfactor_sep(Int kid, Int lvl, Int team_leader);
    #endif

    BASKER_INLINE
    #ifdef BASKER_KOKKOS
    int t_nfactor_sep_old(Int kid, Int lvl, Int team_leader, const TeamMember &thread);
    #else
    int t_nfactor_sep_old(Int kid, Int lvl, Int team_leader);
    #endif

    int t_nfactor_sep2(const Int kid, const Int lvl, const Int team_leader, const TeamMember &thread);

    void t_nfactor_sep2_inc_lvl(const Int kid, const Int lvl, const Int team_leader, const TeamMember &thread);

    //BASKER_INLINE
    //int t_init_solve_rhs(Int kid, Entry *x , Entry *b);

    BASKER_INLINE
    void t_reset_BTF_factor(Int kid);
    void t_reset_BTF_factor_top(Int start, Int end);

    void t_reset_ND_factor(Int kid);

    BASKER_INLINE
    void t_init_factor(Int kid);

    inline
    Int t_get_kid(const TeamMember &thread);

    //BTF array
    int t_nfactor_diag(Int kid, Int schunk, Int nchunk);

    INT_1DARRAY   btf_tabs; // stores starting col id (global) of btf blocks
    Int           btf_tabs_offset; // stores offset of first btf block in BTF_C, after the nd blocks BTF_A
    Int           btf_nblks; // set during Symbolic() btf ordering step - total of BTF_C blocks and ND blocks?

    Int           btf_top_tabs_offset; // store offset of "first" btf block in BTF_A (BTF_A has only one block)
    Int           btf_top_nblks;       // store the number of diagonal blocks in BTF_D

    //These are temp arrys that are used for ordering and sfactor
    INT_1DARRAY btf_blk_work;
    INT_1DARRAY btf_blk_nnz;
    INT_1DARRAY btf_schedule;

    Int btf_total_work;
 
  private:

    //basker_tree
    BASKER_INLINE
    void malloc_tree(Int, Int, Int);

    BASKER_INLINE
    void init_tree_struc();

    BASKER_INLINE
    int init_tree(Int*, Int, Int, Int*, Int*, Int*, Int);

    BASKER_INLINE
    int init_tree_thread(Int*, Int, Int, Int*, Int*, Int*);
                         
    BASKER_INLINE
    int init_tree_thread();
  
    BASKER_INLINE
    int init_tree_old(Int*, Int, Int, Int*, Int*, Int*);

    BASKER_INLINE
    void rec_tabs(Int, Int, Int, Int, Int, Int*, Int* , Int*, Int *, INT_1DARRAY, INT_1DARRAY, INT_1DARRAY);

    BASKER_INLINE
    int init_tree_lvl(Int*, Int, Int, Int*, Int*, Int*, Int);

    BASKER_INLINE
    int update_lpinv_tree(Int,Int);

    BASKER_INLINE
    void update_pivot_tree(Int, Int );

    BASKER_INLINE
    void matrix_to_views(BASKER_MATRIX&, MATRIX_VIEW_2DARRAY&);

    BASKER_INLINE
    void matrix_to_views_2D(BASKER_MATRIX&);

    BASKER_INLINE
    void find_2D_convert(BASKER_MATRIX&);

    BASKER_INLINE
    int clean_2d();


    //basker_order.hpp
    BASKER_INLINE
    int default_order();

    BASKER_INLINE
    void user_order(Int *perm, Int nblks, Int parts, Int *row_tabs, Int *col_tabs, Int *tree_tabs);

    BASKER_INLINE
    int btf_order();

    BASKER_INLINE
    int btf_order2();

    BASKER_INLINE
    void order_incomplete();

    BASKER_INLINE
    int partition(int option);

    BASKER_INLINE
    int match_ordering(int option);

    BASKER_INLINE
    int apply_scotch_partition(BASKER_BOOL keep_zeros = BASKER_TRUE,
                               BASKER_BOOL compute_nd = BASKER_TRUE,
                               BASKER_BOOL apply_nd   = BASKER_TRUE);

    BASKER_INLINE
    int scotch_partition(BASKER_MATRIX &M, BASKER_BOOL apply_nd = BASKER_TRUE);

    BASKER_INLINE
    int scotch_partition(BASKER_MATRIX &M, BASKER_MATRIX &MMT, BASKER_BOOL apply_nd = BASKER_TRUE);

    BASKER_INLINE
    int permute_inv(INT_1DARRAY, INT_1DARRAY, Int);

    BASKER_INLINE
    int permute_inv(ENTRY_1DARRAY, INT_1DARRAY, Int);

    BASKER_INLINE
    int permute_inv(ENTRY_1DARRAY, INT_1DARRAY, Int, Int, Int);

    BASKER_INLINE
    int permute(ENTRY_1DARRAY, INT_1DARRAY, Int);

    BASKER_INLINE
    int permute(ENTRY_1DARRAY, INT_1DARRAY, Int, Int, Int);

    BASKER_INLINE
    int permute(BASKER_MATRIX &M, INT_1DARRAY row, INT_1DARRAY col);

  // NDE - added routines for solve performance improvement
    BASKER_INLINE
    int permute(Entry*, INT_1DARRAY, Int);

    BASKER_INLINE
    int permute(INT_1DARRAY, INT_1DARRAY, Int);

    BASKER_INLINE
    int permute_array_with_workspace(Entry * vec,
                                     INT_1DARRAY & p,
                                     Int n);

    BASKER_INLINE
    int permute_with_workspace(INT_1DARRAY & vec,
                               INT_1DARRAY & p,
                               Int n,
                               Int istart = 0,
                               Int offset = 0);

    BASKER_INLINE
    int permute_inv_with_workspace(INT_1DARRAY & vec,
                                   INT_1DARRAY & p,
                                   Int n,
                                   Int istart = 0,
                                   Int offset = 0);

    BASKER_INLINE
    int permute_with_workspace(ENTRY_1DARRAY&, INT_1DARRAY&, Int);

    BASKER_INLINE
    int permute_inv_with_workspace(ENTRY_1DARRAY&, INT_1DARRAY&, Int);

    BASKER_INLINE
    int permute_inv_array_with_workspace(Entry*, INT_1DARRAY&, Int);

    BASKER_INLINE
    int permute_inv_and_init_for_solve
    (
     Entry* ,
     ENTRY_1DARRAY &,
     ENTRY_1DARRAY &,
     INT_1DARRAY &,
     Int 
    );

    BASKER_INLINE
    int permute_and_init_for_solve
    (
     Entry* y,
     ENTRY_1DARRAY &xcon,
     ENTRY_1DARRAY &ycon,
     INT_1DARRAY  &p, 
     Int n
    );

    BASKER_INLINE
    int permute_inv_and_finalcopy_after_solve
    (
     Entry* x,
     ENTRY_1DARRAY &xconv,
     ENTRY_1DARRAY &yconv,
     INT_1DARRAY  &p,
     Int n
    );

    BASKER_INLINE
    int permute_and_finalcopy_after_solve
    ( 
     Entry* ,
     ENTRY_1DARRAY &,
     ENTRY_1DARRAY &,
     INT_1DARRAY &,
     Int
    );

    BASKER_INLINE
    void permute_composition_for_solve(Int);
  // end NDE

    BASKER_INLINE
    int permute_col(BASKER_MATRIX &M, INT_1DARRAY col, Int frow = 0);

    BASKER_INLINE
    int permute_row(BASKER_MATRIX &M, INT_1DARRAY row);

    BASKER_INLINE
    int permute_row(Int nnz, Int *row_idx, Int *row);

    BASKER_INLINE
    int sort_matrix(Int nnz, Int ncol, Int *col_ptr, Int *row_idx, Entry *val);

    BASKER_INLINE
    int sort_matrix(BASKER_MATRIX &M);

    //basker_order_match.hpp
    BASKER_INLINE
    int mwm(BASKER_MATRIX &M, INT_1DARRAY _perm);

    BASKER_INLINE
    int mc64(Int n_, Int nnz_, Int *colptr, Int *rowidx, Entry *val,
             Int _job, Int *_perm, Entry *_scale_row, Entry *_scale_col);

    BASKER_INLINE
    int mc64(BASKER_MATRIX &M, Int _job, INT_1DARRAY _perm, ENTRY_1DARRAY _scale_row, ENTRY_1DARRAY _scale_col);

    BASKER_INLINE
    int mc64(Int _job, INT_1DARRAY _perm, ENTRY_1DARRAY _scale_row, ENTRY_1DARRAY _scale_col);

    //basker_order_scotch.hpp
    BASKER_INLINE
    int AplusAT(BASKER_MATRIX &M, BASKER_MATRIX &C, BASKER_BOOL keep_zeros = BASKER_TRUE);

    int part_scotch(BASKER_MATRIX &M, BASKER_TREE &BT);

    BASKER_INLINE
    int part_scotch(BASKER_MATRIX &M, BASKER_TREE &BT, Int num_domains);

    void to_complete_tree(Int lvl, Int iblks, Int nblks, INT_1DARRAY tabs, INT_1DARRAY tree);

    void rec_build_tree(Int lvl, Int &lpos, Int &rpos, Int &mynum, INT_1DARRAY tree);

    
    //basker_order_btf.hpp
    BASKER_INLINE
    int find_btf(BASKER_MATRIX &M);

    BASKER_INLINE
    int find_btf2(BASKER_MATRIX &M);

    BASKER_INLINE
    int break_into_parts(BASKER_MATRIX &M, Int nblks, INT_1DARRAY btf_tabs);

    BASKER_INLINE
    int break_into_parts2(BASKER_MATRIX &M, Int nblks, INT_1DARRAY btf_tabs);
    
    BASKER_INLINE
    void find_btf_schedule(BASKER_MATRIX &M, Int nblks, INT_1DARRAY btf_tabs);

    /*
    BASKER_INLINE
    int strong_component(BASKER_MATRIX &M,
                         Int &nblks,
                         INT_1DARRAY &perm_in,
                         INT_1DARRAY &perm,
                         INT_1DARRAY &CC);
    */

    BASKER_INLINE
    int strong_component(BASKER_MATRIX &M,Int &nblks, INT_1DARRAY &perm, INT_1DARRAY &CC);

    //basker_sfactor.hpp
    BASKER_INLINE 
    int sfactor();

    BASKER_INLINE
    int symmetric_sfactor();

    BASKER_INLINE
    int unsymmetric_sfactor();

    BASKER_INLINE
    void e_tree(BASKER_MATRIX &MV, BASKER_SYMBOLIC_TREE &ST, Int ata_option);

    BASKER_INLINE
    int sfactor_copy();

    BASKER_INLINE
    int sfactor_copy2(bool alloc_BTFA = false, bool copy_BTFA = true);


    //old
    BASKER_INLINE
    void e_tree(BASKER_MATRIX_VIEW &MV, BASKER_SYMBOLIC_TREE &ST, Int ata_option);

    BASKER_INLINE
    void post_order(BASKER_MATRIX &MV, BASKER_SYMBOLIC_TREE &ST);

    //old
    BASKER_INLINE
    void post_order(BASKER_MATRIX_VIEW &MV, BASKER_SYMBOLIC_TREE &ST);
    BASKER_INLINE
    Int post_dfs
    (
     Int j, 
     Int k,
     Int *head,
     Int *next, 
     INT_1DARRAY post,
     Int *stack
    );

    BASKER_INLINE
    void col_count(BASKER_MATRIX &MV, BASKER_SYMBOLIC_TREE &ST);

    //old
    BASKER_INLINE
    void col_count(BASKER_MATRIX_VIEW &MV, BASKER_SYMBOLIC_TREE &ST);
 
    BASKER_INLINE
    Int least_common
    (
     Int i, 
     Int j, 
     Int* first,
     Int *mfirst, 
     Int *pleaf, 
     Int *past, 
     Int *jleaf
    );

    BASKER_INLINE
    void U_blk_sfactor
    (
     BASKER_MATRIX &MV,
     BASKER_SYMBOLIC_TREE &ST, 
     INT_1DARRAY gcol,
     INT_1DARRAY grow, 
     Int off_diag
    );

    //old
    BASKER_INLINE
    void U_blk_sfactor
    (
     BASKER_MATRIX_VIEW &MV,
     BASKER_SYMBOLIC_TREE &ST, 
     INT_1DARRAY gcol,
     INT_1DARRAY grow, 
     Int off_diag
    );

    BASKER_INLINE
    void L_blk_sfactor
    (
     BASKER_MATRIX &MV,
     BASKER_SYMBOLIC_TREE &ST, 
     INT_1DARRAY gcol, 
     INT_1DARRAY grow
    );

    //old
    BASKER_INLINE
    void L_blk_sfactor
    (
     BASKER_MATRIX_VIEW &MV,
     BASKER_SYMBOLIC_TREE &ST, 
     INT_1DARRAY gcol, 
     INT_1DARRAY grow
    );

    BASKER_INLINE
    void S_sfactor_reduce
    (
     BASKER_MATRIX &MV,
     BASKER_SYMBOLIC_TREE &ST,
     INT_1DARRAY gcol,
     INT_1DARRAY grow
    );

    //old
    BASKER_INLINE
    void S_sfactor_reduce
    (
     BASKER_MATRIX_VIEW &MV,
     BASKER_SYMBOLIC_TREE &ST,
     INT_1DARRAY gcol,
     INT_1DARRAY grow
    );

    BASKER_INLINE
    void S_blk_sfactor
    (
     BASKER_MATRIX &MV,
     BASKER_SYMBOLIC_TREE &ST, 
     INT_1DARRAY gcol,
     INT_1DARRAY grow
    );

    //old
    BASKER_INLINE
    void S_blk_sfactor
    (
     BASKER_MATRIX_VIEW &MV,
     BASKER_SYMBOLIC_TREE &ST, 
     INT_1DARRAY gcol,
     INT_1DARRAY grow
    );

    BASKER_INLINE
    void leaf_assign_nnz
    (
     BASKER_MATRIX &M,
     BASKER_SYMBOLIC_TREE &ST,
     Int option
    );

    BASKER_INLINE
    void U_assign_nnz
    (
     BASKER_MATRIX &M,
     BASKER_SYMBOLIC_TREE &ST,
     double fill_factor,
     Int option
    );

    BASKER_INLINE
    void L_assign_nnz
    (
     BASKER_MATRIX &M,
     BASKER_SYMBOLIC_TREE &ST,
     double fill_factor,
     Int option
    );

    BASKER_INLINE
    void S_assign_nnz
    (
     BASKER_MATRIX &M,
     BASKER_SYMBOLIC_TREE &ST,
     Int option
    );

    BASKER_INLINE
    void btf_last_dense(bool flag);

    BASKER_INLINE
    int factor_inc_lvl(Int Option);

    //basker_sfactor_inc.hpp
    BASKER_INLINE
    int sfactor_inc();

    BASKER_INLINE
    void sfactor_nd_estimate();

    BASKER_INLINE
    void sfactor_nd_dom_estimate(BASKER_MATRIX &M, BASKER_MATRIX &LM, BASKER_MATRIX &UM);

    BASKER_INLINE
    void sfactor_nd_lower_estimate(BASKER_MATRIX &M, BASKER_MATRIX &ML);

    BASKER_INLINE
    void sfactor_nd_upper_estimate(BASKER_MATRIX &M, BASKER_MATRIX &UM);

    BASKER_INLINE
    void sfactor_nd_sep_upper_estimate(BASKER_MATRIX &M, BASKER_MATRIX &UM);

    BASKER_INLINE
    void sfactor_nd_sep_lower_estimate(BASKER_MATRIX &M, BASKER_MATRIX &LM);

    BASKER_INLINE
    void sfactor_nd_sep_estimate(BASKER_MATRIX &M, BASKER_MATRIX &ML, BASKER_MATRIX &MU);

    //basker_nfactor.hpp
    BASKER_INLINE
    int factor_token(Int option);

    BASKER_INLINE
    int factor_notoken(Int option);

    BASKER_INLINE
    int t_factor_tree(Int kid);

    BASKER_INLINE
    int copy_schedule(INT_2DARRAY &s, INT_2DARRAY &ls, Int l, Int sl, Int t);

    BASKER_INLINE
    int nfactor_domain_error(INT_1DARRAY);
    
    BASKER_INLINE
    int nfactor_sep_error(INT_1DARRAY);

    BASKER_INLINE
    int nfactor_diag_error(INT_1DARRAY, INT_1DARRAY);
    
    BASKER_INLINE
    void reset_error();
    
    
    //BASKER_INLINE
    inline
    void t_prune(const Int, const Int, const Int, const Int, const Int);

    inline 
    void t_local_reach_short(const Int,const Int, const Int, const Int, Int&);

    inline
    void t_local_reach(const Int, const Int, const Int, Int, Int &);
    
    inline 
    void t_local_reach_short_inc_rlvl(const Int,const Int, const Int, const Int, Int&);

    inline
    void t_local_reach_inc_rlvl(const Int, const Int, const Int, Int, Int &);

    inline
    int t_local_reach_old(Int,Int,Int,Int,Int*);

    BASKER_INLINE
    int t_local_reach_old_old(Int,Int,Int,Int,Int*);

    BASKER_INLINE
    int t_local_reach_inc_lvl(Int,Int,Int,Int,Int*);

    inline
    int t_back_solve(Int,Int,Int,Int,Int,Int);

    BASKER_INLINE
    int t_back_solve_old(Int,Int,Int,Int,Int,Int);

    BASKER_INLINE
    int t_back_solve_inc_lvl(Int,Int,Int,Int,Int,Int);

    BASKER_INLINE
    int t_back_solve_inc_rlvl(Int,Int,Int,Int,Int,Int,Entry&);

    BASKER_INLINE
    int t_upper_col_factor(Int kid, Int team_leader, Int lvl, Int l, Int k, BASKER_BOOL);

    BASKER_INLINE
    int t_upper_col_factor_inc_lvl(Int kid, Int team_leader, Int lvl, Int l, Int k, BASKER_BOOL);

    BASKER_INLINE
    int t_upper_col_factor_old(Int kid, Int team_leader, Int lvl, Int l, Int k, BASKER_BOOL);

    BASKER_INLINE
    int t_upper_col_factor_offdiag(Int kid, Int lvl, Int l, Int k);

    BASKER_INLINE
    int t_upper_col_factor_offdiag_old(Int kid, Int lvl, Int l, Int k);

    BASKER_INLINE
    int t_back_solve_atomic(Int kid, Int team_leader, Int lvl, Int l, Int k, Int top, Int xnnz);

    BASKER_INLINE
    int t_lower_col_factor(Int kid, Int team_leader, Int lvl, Int l, Int k, Entry &opivot);

     BASKER_INLINE
    int t_lower_col_factor_inc_lvl(Int kid, Int team_leader, Int lvl, Int l, Int k, Entry &opivot);
   
    BASKER_INLINE
    int t_lower_col_factor_old(Int kid, Int team_leader, Int lvl, Int l, Int k, Entry &opivot);

    BASKER_INLINE
    int t_lower_col_factor_offdiag(Int kid, Int lvl, Int l, Int k, Entry pivot);

    BASKER_INLINE
    int t_lower_col_factor_offdiag_old(Int kid, Int lvl, Int l, Int k, Entry pivot);

    BASKER_INLINE
    int t_col_barrier(Int kid);

    BASKER_INLINE
    int t_dense_move_offdiag_L(Int kid, 
                         Int blkcol, Int blkrow,
                         Int X_col, Int X_row,
                         Int k , Entry pivot);

     BASKER_INLINE
    int t_dense_move_offdiag_L_inc_lvl(Int kid, 
                         Int blkcol, Int blkrow,
                         Int X_col, Int X_row,
                         Int k , Entry pivot);

    BASKER_INLINE
    int t_dense_move_offdiag_L_inc_lvl_old(Int kid, 
                         Int blkcol, Int blkrow,
                         Int X_col, Int X_row,
                         Int k , Entry pivot);


    BASKER_INLINE
    int t_move_offdiag_L(Int kid, 
                         Int blkcol, Int blkrow,
                         Int X_col, Int X_row,
                         Int k , Entry pivot);

    BASKER_INLINE
    int t_move_offdiag_L_inc_lvl(Int kid, 
                         Int blkcol, Int blkrow,
                         Int X_col, Int X_row,
                         Int k , Entry pivot);

    BASKER_INLINE
    int t_move_offdiag_L_old(Int kid, 
                         Int blkcol, Int blkrow,
                         Int X_col, Int X_row,
                         Int k , Entry pivot);

    BASKER_INLINE
    int t_dense_back_solve_offdiag(Int kid,
       Int blkcol, Int blkrow,
       Int X_col, Int X_row,
       Int k, Int &view_offset,
       ENTRY_1DARRAY x,
       INT_1DARRAY x_indx,
       Int x_size, Int x_offset,
       BASKER_BOOL A_option);

    BASKER_INLINE
    int t_dense_back_solve_offdiag_inc_lvl(Int kid,
       Int blkcol, Int blkrow,
       Int X_col, Int X_row,
       Int k, Int &view_offset,
       ENTRY_1DARRAY x,
       INT_1DARRAY x_indx,
       INT_1DARRAY x_fill,
       Int x_size, Int x_offset,
       BASKER_BOOL A_option);

    BASKER_INLINE
    void t_same_pattern_back_solve_offdiag_inc_lvl(Int kid,
       Int blkcol, Int blkrow,
       Int X_col, Int X_row,
       Int UP_col, Int UP_row,
       Int LP_col, Int LP_row,
       Int k, Int &view_offset,
       ENTRY_1DARRAY x,
       INT_1DARRAY x_indx,
       INT_1DARRAY x_fill,
       Int x_size, Int x_offset,
       BASKER_BOOL A_option);

    BASKER_INLINE
    int t_dense_back_solve_offdiag_inc_lvl_old(Int kid,
      Int blkcol, Int blkrow,
      Int X_col, Int X_row,
      Int k, Int &view_offset,
      ENTRY_1DARRAY x,
      INT_1DARRAY x_indx,
      INT_1DARRAY x_fill,
      Int x_size, Int x_offset,
      BASKER_BOOL A_option);
    
    BASKER_INLINE
    void t_dom_lower_col_offdiag_find_fill(const Int kid, const Int pbrow,
                                           const Int blkcol, const Int blkrow,
                                           const Int X_col, const Int X_row,
                                           const Int k,
                                           INT_1DARRAY x_idx,
                                           const Int x_size,
                                           const Int x_offset,
                                           const BASKER_BOOL A_option);

    BASKER_INLINE
    int t_lower_col_offdiag_find_fill(Int kid,
                                      Int blkcol, Int blkrow,
                                      Int X_col, Int X_row,
                                      Int k,
                                      ENTRY_1DARRAY x,
                                      INT_1DARRAY x_idx,
                                      INT_1DARRAY x_fill,
                                      Int x_size, Int x_offset);
    
    BASKER_INLINE
    int t_lower_col_offdiag_find_fill_rlvl(Int kid,
                                      Int blkcol, Int blkrow,
                                      Int X_col, Int X_row,
                                      Int k,
                                      ENTRY_1DARRAY x,
                                      INT_1DARRAY x_idx,
                                      INT_1DARRAY x_fill,
                                      Int x_size, Int x_offset);

    BASKER_INLINE
    void t_populate_col_fill(const Int kid,
                             const Int blkcol, const Int blkrow,
                             const Int X_col, const Int X_row,
                             const Int k, 
                             const BASKER_BOOL lower );

    BASKER_INLINE
    void t_reduce_col_fill(const Int kid, const Int lvl,
                           const Int sl, const Int l,
                           const Int k, const BASKER_BOOL lower);

    BASKER_INLINE
    int t_back_solve_offdiag(Int kid,
                             Int blkcol, Int blkrow,
                             Int X_col, Int X_row,
                             Int k, Int &view_offset,
                             ENTRY_1DARRAY x,
                             INT_1DARRAY x_indx,
                             Int x_size, Int x_offset,
                             BASKER_BOOL A_option);

    BASKER_INLINE
    int t_back_solve_offdiag_old(Int kid,
                             Int blkcol, Int blkrow,
                             Int X_col, Int X_row,
                             Int k, Int &view_offset,
                             ENTRY_1DARRAY x,
                             INT_1DARRAY x_indx,
                             Int x_size, Int x_offset,
                             BASKER_BOOL A_option);

    BASKER_INLINE
    int t_back_solve_offdiag_inc_lvl(Int kid, Int pbrow,
                             Int blkcol, Int blkrow,
                             Int X_col, Int X_row,
                             Int k, Int &view_offset,
                             ENTRY_1DARRAY x,
                             INT_1DARRAY x_indx,
                             Int x_size, Int x_offset,
                             BASKER_BOOL A_option);

     BASKER_INLINE
     void t_back_solve_offdiag_same_pattern_inc_lvl(Int kid, Int pbrow,
                             Int blkcol, Int blkrow,
                             Int X_col, Int X_row,
                             Int k, Int &view_offset,
                             ENTRY_1DARRAY x,
                             INT_1DARRAY x_indx,
                             Int x_size, Int x_offset,
                             BASKER_BOOL A_option);
   
    BASKER_INLINE
    int t_back_solve_offdiag_inc_lvl_old(Int kid, Int pbrow,
                             Int blkcol, Int blkrow,
                             Int X_col, Int X_row,
                             Int k, Int &view_offset,
                             ENTRY_1DARRAY x,
                             INT_1DARRAY x_indx,
                             Int x_size, Int x_offset,
                             BASKER_BOOL A_option);

    BASKER_INLINE
    int t_dense_blk_col_copy_atomic(Int kid, Int team_leader, 
           Int lvl, Int l, Int k);

    BASKER_INLINE
    int t_blk_col_copy_atomic(Int kid, Int team_leader,
                              Int lvl, Int l, Int k);

    BASKER_INLINE
    int t_dense_copy_update_matrix(Int kid, Int team_leader,
                             Int lvl, Int l, Int k);

    BASKER_INLINE
    int t_copy_update_matrix(Int kid, Int team_leader,
                             Int lvl, Int l, Int k);

    BASKER_INLINE
    int t_copy_update_matrix_old(Int kid, Int team_leader,
                             Int lvl, Int l, Int k);

    void t_add_extend(const TeamMember &thread,
                      const Int kid, 
                      const Int lvl, 
          const Int l,
                      const Int k,
                      const Int k_offset,
                      const BASKER_BOOL lower);

    void t_add_extend_inc_lvl(const TeamMember &thread,
                      const Int kid, 
                      const Int lvl, 
          const Int l,
                      const Int k,
                      const Int k_offset,
                      const BASKER_BOOL lower);

    void t_upper_col_factor_offdiag2(const Int kid,
                                  const Int lvl, 
                                  const Int sl,
                                  const Int l,
                                  const Int k, 
                                  const BASKER_BOOL lower);

    void t_upper_col_ffactor_offdiag2_inc_lvl(const Int kid,
                                  const Int lvl, 
                                  const Int sl,
                                  const Int l,
                                  const Int k, 
                                  const BASKER_BOOL lower);

    void t_upper_col_factor_offdiag2_same_pattern_inc_lvl(const Int kid,
                                  const Int lvl, 
                                  const Int sl,
                                  const Int l,
                                  const Int k, 
                                  const BASKER_BOOL lower);

    void t_upper_col_factor_offdiag2_inc_lvl(const Int kid,
                                  const Int lvl, 
                                  const Int sl,
                                  const Int l,
                                  const Int k, 
                                  const BASKER_BOOL lower);

    int t_lower_col_diag_find_fill(Int kid, 
                                  Int blkcol, 
          Int blkrow,
                                  Int X_col, 
          Int X_row,
                                  Int k,
                                  ENTRY_1DARRAY x,
                                  INT_1DARRAY x_idx,
                                  INT_1DARRAY x_fill,
                                  Int x_size, 
          Int x_offset);
   
    void t_dense_blk_col_copy_atomic2(const Int kid, 
                                  const Int team_leader,
                                  const Int lvl,
                                  const Int sl,
                                  const Int l,
                                  const Int k, 
                                  const BASKER_BOOL lower);

    void t_dense_blk_col_copy_atomic2_inc_lvl(const Int kid, 
                                  const Int team_leader,
                                  const Int lvl,
                                  const Int sl,
                                  const Int l,
                                  const Int k, 
                                  const BASKER_BOOL lower);

    void t_same_pattern_col_copy_inc_lvl(const Int kid,
                                  const Int lvl,
                                  const Int sl,
                                  const Int l,
                                  const Int k, 
                                  const BASKER_BOOL lower);

    void t_dense_copy_update_matrix2(const Int kid,
                                  const Int team_leader,
                                   const Int lvl, 
          const Int l,
                                  const Int k);
    
    void t_dense_copy_update_matrix2_inc_lvl(const Int kid,
                                  const Int team_leader,
                                  const Int lvl, 
          const Int l,
                                  const Int k);

    void t_same_pattern_update_matrix_inc_lvl(const Int kid,
                                        const Int team_leader,
                                        const Int lvl,
                                        const Int l, 
                                        const Int k);

    void t_lower_col_factor_offdiag2(const Int kid,
                                  const Int lvl,
                                  const Int l,
                                  const Int k,
                                  Entry pivot);

    void t_lower_col_factor_offdiag2_inc_lvl(const Int kid,
                                  const Int lvl,
                                  const Int l,
                                  const Int k,
                                  Entry pivot);

    void t_lower_col_factor_offdiag2_cleanup_inc_lvl(const Int kid,
                                  const Int lvl,
                                  const Int l,
                                  const Int k);
    
    void t_add_orig_fill(const Int kid, const Int lvl,
                            const Int l, 
          const Int k, 
          const BASKER_BOOL lower);

    BASKER_INLINE
    Int find_leader(Int kid, Int l);

    BASKER_INLINE
    Int find_leader_inc_lvl(Int kid, Int l);


    //basker_nfactor_diag
    BASKER_INLINE
    int t_single_nfactor(Int kid, Int c);

    BASKER_INLINE
    int t_blk_nfactor(Int kid, Int c);

    BASKER_FINLINE
    void t_local_reach_short_btf(const Int, const Int, Int &);

    BASKER_INLINE
    int t_local_reach_btf(Int, BASKER_MATRIX&,Int,Int&, Int,Int);

    void t_prune_btf(const Int, const BASKER_MATRIX &, const BASKER_MATRIX&, const Int, const Int);

    BASKER_INLINE
    int t_local_reach_old(Int,BASKER_MATRIX&,Int,Int,Int,Int*);

    BASKER_INLINE
    int t_back_solve(Int,Int,BASKER_MATRIX&,Int,Int,Int,Int,Int);

    //basker_thread.hpp
    //BASKER_INLINE
    inline
    void t_basker_barrier(const TeamMember &thread,
                          const Int my_kid,
                          const Int leader_kid, 
                          const Int size,
                          const Int function_n,
                          const Int k, 
                          const Int l);
                          
    inline
    void t_basker_barrier_inc_lvl(const TeamMember &thread,
                          const Int my_kid,
                          const Int leader_kid, 
                          const Int size,
                          const Int function_n,
                          const Int k, 
                          const Int l);

    BASKER_INLINE
    void t_basker_barrier_old(const TeamMember &thread,
                          const Int leader_kid,
                          const Int sublvl,
                          const Int function_n,
                          const Int size);

    //basker_util.hpp
    //Memory Util
    //On host
    BASKER_INLINE
    static
    void init_value(INT_1DARRAY, Int, Int);
    
    BASKER_INLINE
    static
    void init_value(INT_1DARRAY, Int, Int*);
    
    BASKER_INLINE
    static
    void init_value(ENTRY_1DARRAY, Int, Entry);

    BASKER_INLINE
    static
    void init_value(ENTRY_1DARRAY, Int, Entry*);

    BASKER_INLINE
    static
    void init_value(BOOL_1DARRAY, Int, BASKER_BOOL);

    BASKER_INLINE
    static
    void init_value(BOOL_1DARRAY, Int, BASKER_BOOL*);

    void init_value(INT_1DARRAY, Int, Int, Int);

    void init_value(ENTRY_1DARRAY, Int, Entry, Int);
    
    //Workspace Util
    
    //Print Util
    void print_factor(BASKER_MATRIX &L, BASKER_MATRIX &U);
    int printL();
    int printL2D();
    int printU();
    int printLMTX();
    int printUMTX();
    void printMTX(std::string fname, BASKER_MATRIX &M);
    void printMTX(std::string fname, BASKER_MATRIX &M, BASKER_BOOL  off);
    void readMTX(std::string fname, BASKER_MATRIX &M);
    int printRHS();
    int printSOL();
    void printTree();

    BASKER_INLINE
    int get_L(Int &n, Int &nnz, Int **col_ptr, Int **row_idx, Entry **val);

    BASKER_INLINE
    int get_U(Int &n, Int &nnz, Int **col_ptr, Int **row_idx, Entry **val);

    BASKER_INLINE
    int get_p(Int **p);

    BASKER_INLINE
    void printVec(INT_1DARRAY, Int);
    
    BASKER_INLINE
    void printVec(ENTRY_1DARRAY, Int);

    BASKER_INLINE
    void printVec(std::string, INT_1DARRAY, Int);

    BASKER_INLINE
    void printVec(std::string, ENTRY_1DARRAY, Int);

    BASKER_INLINE
    void printVec(std::string, BASKER_ENTRY*, Int);

    void get_total_perm(INT_1DARRAY, INT_1DARRAY);

    //inline
    //Int t_get_kid(const TeamMember &thread);

    void print_sep_bal();

    //Matrix helper
    BASKER_INLINE
    void matrix_transpose(BASKER_MATRIX &M,
                          BASKER_MATRIX &MT,
                          BASKER_BOOL keep_zeros = BASKER_TRUE);

    BASKER_INLINE
    void matrix_transpose(BASKER_MATRIX_VIEW &, BASKER_MATRIX &);

    BASKER_INLINE
    void matrix_transpose(
        const Int sm_, 
        const Int m_,
                          const Int sn_, 
        const Int n_,
                          const Int nnz_,
                          Int *col_ptr,
                          Int *row_idx,
                          Entry *val,
                          BASKER_MATRIX &AT);

    BASKER_INLINE
    void matrix_transpose(
        const Int sm_, 
        const Int m_,
                          const Int sn_, 
        const Int n_,
                          const Int nnz_,
                          Int *col_ptr,
                          Int *row_idx,
                          Entry *val,
                          BASKER_MATRIX &AT,
        INT_1DARRAY &vals_transpose_local);

    //basker_solve_rhs.hpp
    BASKER_INLINE
    int solve_interface(Entry *, Entry*);

    BASKER_INLINE
    int solve_interface(Int, Entry *, Entry*);

    BASKER_INLINE
    int solve_interface(ENTRY_1DARRAY &, ENTRY_1DARRAY &);

    BASKER_INLINE
    int serial_solve(ENTRY_1DARRAY &, ENTRY_1DARRAY &);

    BASKER_INLINE
    int serial_forward_solve(ENTRY_1DARRAY &, ENTRY_1DARRAY &);

    BASKER_INLINE
    int serial_backward_solve(ENTRY_1DARRAY &, ENTRY_1DARRAY &);

    BASKER_INLINE
    int serial_btf_solve(ENTRY_1DARRAY &, ENTRY_1DARRAY &);

    BASKER_INLINE
    int spmv(BASKER_MATRIX &, ENTRY_1DARRAY, ENTRY_1DARRAY);

    BASKER_INLINE
    int neg_spmv(BASKER_MATRIX &M,
                 ENTRY_1DARRAY x,
                 ENTRY_1DARRAY y,
                 Int offset = 0);

    BASKER_INLINE
    int neg_spmv_perm(BASKER_MATRIX &M,
                      ENTRY_1DARRAY &y,
                      ENTRY_1DARRAY &x,
                      Int offset = 0);

    BASKER_INLINE
    int lower_tri_solve(BASKER_MATRIX &M,
                        ENTRY_1DARRAY &x,
                        ENTRY_1DARRAY &y,
                        Int offset = 0);

    BASKER_INLINE
    int upper_tri_solve(BASKER_MATRIX &M,
                        ENTRY_1DARRAY &x,
                        ENTRY_1DARRAY &y,
                        Int offset = 0);

    BASKER_INLINE
    int spmv_BTF(Int tab,
                 BASKER_MATRIX &M,
                 ENTRY_1DARRAY &x, // modified rhs
                 ENTRY_1DARRAY &y,
                 bool full = true);


    BASKER_INLINE
    int solve_interfacetr(Entry *, Entry*);

    //BASKER_INLINE
    int solve_interfacetr(Int, Entry *, Entry*);

    BASKER_INLINE
    int solve_interfacetr(ENTRY_1DARRAY &, ENTRY_1DARRAY &);

    BASKER_INLINE
    int serial_btf_solve_tr(ENTRY_1DARRAY &, ENTRY_1DARRAY &);

    BASKER_INLINE
    int l_tran_brfa_solve(ENTRY_1DARRAY &, ENTRY_1DARRAY &);

    BASKER_INLINE
    int u_tran_btfa_solve(ENTRY_1DARRAY &, ENTRY_1DARRAY &);

    BASKER_INLINE
    int neg_spmv_tr(BASKER_MATRIX &M,
                 ENTRY_1DARRAY x,
                 ENTRY_1DARRAY y,
                 Int offset = 0);

    BASKER_INLINE
    int neg_spmv_perm_tr(BASKER_MATRIX &M,
                      ENTRY_1DARRAY &y,
                      ENTRY_1DARRAY &x,
                      Int offset = 0);

    BASKER_INLINE
    int lower_tri_solve_tr(BASKER_MATRIX &M,
                        ENTRY_1DARRAY &x,
                        ENTRY_1DARRAY &y,
                        Int offset = 0);

    BASKER_INLINE
    int upper_tri_solve_tr(BASKER_MATRIX &M,
                        ENTRY_1DARRAY &x,
                        ENTRY_1DARRAY &y,
                        Int offset = 0);

    BASKER_INLINE
    int spmv_BTF_tr(Int tab,
                 BASKER_MATRIX &M,
                 ENTRY_1DARRAY &x, // modified rhs
                 ENTRY_1DARRAY &y,
                 bool full = true);



    //basker_stats.hpp
    BASKER_INLINE
    void print_local_time_stats();

    BASKER_INLINE
    Int get_Lnnz();

    BASKER_INLINE
    Int get_Unnz();

    /*----------------Package Variables----------------*/

    /* ----------------TYPEDEF TYPES-------------------*/

    bool crs_transpose_needed;
    //OLD
    //For BTF Option
    // [BTF_A BTF_B]
    // [0     BTF_C]

    //NEW
    //For BTF Option
    // [BTF_D  BTF_E]
    // [0      BTF_A BTF_B]
    // [0      0     BTF_C]

    BASKER_MATRIX A;
    BASKER_MATRIX BTF_A;
    BASKER_MATRIX BTF_C;
    BASKER_MATRIX BTF_B;
    BASKER_MATRIX BTF_D;
    BASKER_MATRIX BTF_E;

    Int gn;
    Int gm;

    // view of views of 2D blocks; stores CCS of the BTF_A upper matrix 2D blocks; btf_tabs_offset blocks in BTF_A
    MATRIX_2DARRAY       AVM;
    MATRIX_2DARRAY       ALM;

    BASKER_MATRIX At;
    
    MATRIX_2DARRAY LL;   // view of views of 2D blocks; stores CCS factored ALM
    MATRIX_2DARRAY LU;   // view of views of 2D blocks; stores CCS factored AVM
    INT_1DARRAY LL_size; // tracks the number of 2D blocks ('rows') in a given 'column'
    INT_1DARRAY LU_size;

    //Used for BTF
#define BASKER_SPLIT_A
#if defined(BASKER_SPLIT_A)    
    MATRIX_1DARRAY L_D; //lower blocks for BTF_D; total of btf_top_nblks
    MATRIX_1DARRAY U_D; //upper blocks for BTF_D
#endif
    MATRIX_1DARRAY LBTF; //lower blocks for BTF_C; total of btf_nblks - btf_tabs_offset
    MATRIX_1DARRAY UBTF; //upper blocks for BTF_C
    
    //Thread Array 
    //2D-1D Format, stores workspace and token
    //2D Format, stores only token
    //Could be used for light weight task scheduling in future
    THREAD_1DARRAY  thread_array;

    //INT_2DARRAY lvl_task;
    INT_2DARRAY  S; //schedule
                    // S maps a tree level and thread id to 2D col id

    //Made 2d Inorder to be able to use ref in Kokkos
    INT_1DARRAY gperm; //global perm due to pivoting
    INT_1DARRAY gpermi;
    INT_1DARRAY gperm_same;
    INT_1DARRAY gperm_array;  // used to map new rows (after pivot) to original rows
    INT_1DARRAY gpermi_array; // i.e., gperm_array(i) = k means that the current i-th row is the k-th row before pivot


    // NDE
    ENTRY_1DARRAY x_view_ptr_scale;
    ENTRY_1DARRAY y_view_ptr_scale;

    ENTRY_1DARRAY x_view_ptr_copy;
    ENTRY_1DARRAY y_view_ptr_copy;
    INT_1DARRAY perm_inv_comp_array;
    INT_1DARRAY perm_comp_array;

    INT_1DARRAY perm_comp_iworkspace_array;
    ENTRY_1DARRAY perm_comp_fworkspace_array;

    // Matrix dims stored within Symbolic
    Int sym_gn;
    Int sym_gm;

    // sfactor_copy2 mapping of input vals to reordered vals
    INT_1DARRAY vals_perm_composition; //this will store the btf permutation+sorts of val (for use in Factor)
    INT_1DARRAY_PAIRS vals_block_map_perm_pair; //this will map perm(val) indices to BTF_A, BTF_B, and BTF_C 

    // These store the permutation of indices of val during permute_col and sort calls
    INT_1DARRAY vals_order_btf_array;
    INT_1DARRAY vals_order_blk_amd_array;
    INT_1DARRAY vals_order_scotch_array;
    INT_1DARRAY vals_order_csym_array;

    // These store the permutation indices of the block vals during ND permute_col and sort calls
    INT_1DARRAY vals_order_ndbtfd_array;
    INT_1DARRAY vals_order_ndbtfe_array;
    //
    INT_1DARRAY vals_order_ndbtfa_array;
    INT_1DARRAY vals_order_ndbtfb_array;
    INT_1DARRAY vals_order_ndbtfc_array;

    INT_1DARRAY inv_vals_order_ndbtfd_array;
    INT_1DARRAY inv_vals_order_ndbtfe_array;
    //
    INT_1DARRAY inv_vals_order_ndbtfa_array;
    INT_1DARRAY inv_vals_order_ndbtfb_array;
    INT_1DARRAY inv_vals_order_ndbtfc_array;

    //
    INT_1DARRAY symbolic_row_iperm_array;
    INT_1DARRAY symbolic_col_iperm_array;

    INT_1DARRAY symbolic_row_perm_array;
    INT_1DARRAY symbolic_col_perm_array;

    INT_1DARRAY numeric_row_iperm_array;
    INT_1DARRAY numeric_col_iperm_array;

    // For transpose
    INT_1DARRAY vals_crs_transpose; //this will store shuffling and sort of vals due to transpose

    // To hold the nnz and avoid some compiler errors if BTF_A.nnz undefined, for example
    Int btfd_nnz; 
    Int btfe_nnz;
    //
    Int btfa_nnz; 
    Int btfb_nnz;
    Int btfc_nnz;

    ENTRY_1DARRAY input_vals_unordered; //may be needed if a copy of val input is needed to be stored


    BASKER_INLINE
    int permute_col_store_valperms
    (
     BASKER_MATRIX &M,
     INT_1DARRAY &col,
     INT_1DARRAY &vals_order_perm
    );

    BASKER_INLINE
    int sort_matrix_store_valperms( BASKER_MATRIX &M, INT_1DARRAY &order_vals_perms );

    BASKER_INLINE
    int ndsort_matrix_store_valperms(BASKER_MATRIX &M,
                                     INT_1DARRAY &order_vals_perms,
                                     BASKER_BOOL track_perm = BASKER_TRUE);

    //end NDE


    //RHS and solutions (These are not used anymore)
    ENTRY_2DARRAY rhs;
    ENTRY_2DARRAY sol;
    Int nrhs;

    
    BASKER_TREE   part_tree;
    BASKER_TREE   tree;
    BASKER_SYMBOLIC_TREE stree;
    //std::vector<BASKER_SYMBOLIC_TREE> stree_list;
    BASKER_BOOL   part_tree_saved;
    BASKER_TREE   part_tree_orig;
    

    BASKER_STATS stats;
    
    /*flags*/
    BASKER_BOOL matrix_flag;
    BASKER_BOOL order_flag;
    BASKER_BOOL tree_flag;
    BASKER_BOOL symb_flag;
    BASKER_BOOL workspace_flag;
    BASKER_BOOL factor_flag;
    BASKER_BOOL rhs_flag;
    BASKER_BOOL solve_flag;
    BASKER_BOOL match_flag;
    BASKER_BOOL btf_flag;
    BASKER_BOOL nd_flag;
    BASKER_BOOL amd_flag;
    BASKER_BOOL same_pattern_flag;
   
    Int num_threads;
    Int global_nnz;


    //Don't think we use this anymore
    //Post ordering for symmetric
    //INT_1DARRAY perm_post_order;

    BaskerPointBarrier<Int,Entry,Exe_Space> basker_barrier;

    /*Incomplete Factorization Arrays*/
    //#ifdef BASKER_INC_LVL
    INT_1DARRAY INC_LVL_TEMP;
    INT_1DARRAY INC_LVL_ARRAY_CNT;
    INT_1DARRAY INC_LVL_ARRAY;
    //#endif

    //ordering perms
    //This should be all compounded in the future
    INT_1DARRAY order_match_array;
    INT_1DARRAY order_btf_array;
    INT_1DARRAY order_scotch_array;
    INT_1DARRAY order_csym_array;
    INT_1DARRAY order_c_csym_array;
    INT_1DARRAY order_blk_mwm_array;
    INT_1DARRAY order_blk_mwm_inv;
    //invert ordering
    INT_1DARRAY order_csym_inv;
    INT_1DARRAY order_nd_inv;

    // row/col scaling
    ENTRY_1DARRAY scale_row_array;
    ENTRY_1DARRAY scale_col_array;
    //for experimental 
    INT_1DARRAY order_blk_amd_array;
    INT_1DARRAY order_blk_amd_inv;

    void blk_amd(BASKER_MATRIX &M, INT_1DARRAY p);

    int btf_blk_mwm_amd(Int b_start, Int b_num, BASKER_MATRIX &M,
                        INT_1DARRAY p_mwm, INT_1DARRAY p_amd,
                        INT_1DARRAY btf_nnz, INT_1DARRAY btf_work);
    int btf_blk_mwm_amd(BASKER_MATRIX &M, INT_1DARRAY p_mwm, INT_1DARRAY p_amd,
                        INT_1DARRAY btf_nnz, INT_1DARRAY btf_work);

    //basker_order_amd
    void amd_order(BASKER_MATRIX &M,INT_1DARRAY p);
    
    void csymamd_order(BASKER_MATRIX &M, INT_1DARRAY p, INT_1DARRAY cmember);
  };

}//End namespace Basker

#endif //End ifndef basker_def_hpp
