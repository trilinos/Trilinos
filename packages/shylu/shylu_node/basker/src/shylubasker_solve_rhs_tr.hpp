// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLUBASKER_SOLVE_RHS_TR_HPP
#define SHYLUBASKER_SOLVE_RHS_TR_HPP

/*Basker Includes*/
#include "shylubasker_matrix_decl.hpp"
#include "shylubasker_matrix_view_decl.hpp"
#include "shylubasker_types.hpp"
#include "shylubasker_util.hpp"

/*Kokkos Includes*/
#ifdef BASKER_KOKKOS
#include <Kokkos_Core.hpp>
#include <Kokkos_Timer.hpp>
#else
#include <omp.h>
#endif

/*System Includes*/
#include <iostream>
#include <string>

namespace BaskerNS
{

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::solve_interfacetr
  (
   Int _nrhs,
   Entry *_x, // Solution
   Entry *_y  // rhs
  )
  {
    for(Int r = 0; r < _nrhs; r++)
    {
      solve_interfacetr(&(_x[r*gm]), &(_y[r*gm]));
    }

    return 0;
  }//end solve_interface(_nrhs,x,y);

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::solve_interfacetr
  (
   Entry *_x, // Solution (len = gn)
   Entry *_y  // rhs
  )
  {

    if (Options.blk_matching != 0) {
      for(Int i = 0; i < gn; i++) {
        Int col = symbolic_col_iperm_array(i);
        _y[i] = scale_col_array(col) * _y[i];
      }
    }
    // rhs content from _y has been permuted and copied to x_view_ptr_copy, and will act as the rhs-to-update during solve calls; 
    // y_view_ptr_copy is initialized to be zeros, and will store the pivots (i.e. solutions)
    permute_inv_and_init_for_solve(_y, x_view_ptr_copy, y_view_ptr_copy, perm_comp_array, gn);

    if (Options.matrix_scaling != 0) {
      for (Int i = 0; i < gn; i++) {
        x_view_ptr_copy(i) = x_view_ptr_copy(i) * scale_col_array(i);
      }
    }
    if (Options.blk_matching != 0 || Options.static_delayed_pivot != 0) {
      permute_inv_with_workspace(x_view_ptr_copy, numeric_col_iperm_array, gn);
    }

    // ==================
    // transpose-solve
    // after the solve: solution is stored in [x_view_ptr_copy(1:poffset); y_view_ptr_copy(1:(gn-poffset))]
    solve_interfacetr(x_view_ptr_copy, y_view_ptr_copy);

    if (Options.no_pivot == BASKER_FALSE) {
      // apply partial pivoting from numeric
      permute_and_finalcopy_after_solve(&(y_view_ptr_scale(0)), x_view_ptr_copy, y_view_ptr_copy, gperm, gn);
      // y_view_ptr_scale contains the solution after germi is applied

      // put it back into [x_view_ptr_copy; y_view_ptr_copy]
      const Int poffset = btf_tabs(btf_tabs_offset);
      for (Int i = 0; i < poffset; i++) {
        x_view_ptr_copy(i) = y_view_ptr_scale(i);
      }
      for (Int i = poffset; i < gn; i++) {
        y_view_ptr_copy(i) = y_view_ptr_scale(i);
      }
    }

    if (Options.blk_matching != 0 || Options.static_delayed_pivot != 0) {
      permute_inv_and_finalcopy_after_solve(_x, x_view_ptr_copy, y_view_ptr_copy, numeric_row_iperm_array, gn);
      permute_array_with_workspace(_x, perm_inv_comp_array, gn);

      for(Int i = 0; i < gn; i++) {
        //Int row = order_blk_mwm_array(symbolic_row_iperm_array(i));
        Int row = symbolic_row_perm_array(order_blk_mwm_inv(i));
        _x[i] = scale_row_array(row) * _x[i];
      }
    } else {
      if (Options.matrix_scaling != 0) {
        const Int poffset = btf_tabs(btf_tabs_offset);
        for (Int i = 0; i < poffset; i++) {
          x_view_ptr_copy(i) = x_view_ptr_copy(i) * scale_row_array(i);
        }
        for (Int i = poffset; i < gn; i++) {
          y_view_ptr_copy(i) = y_view_ptr_copy(i) * scale_row_array(i);
        }
      }
      permute_and_finalcopy_after_solve(_x, x_view_ptr_copy, y_view_ptr_copy, perm_inv_comp_array, gn);
    }
    // final solution is permuted back to original ordering and copied to _x

    return 0;
  }

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::solve_interfacetr
  ( 
   ENTRY_1DARRAY & x, // x is permuted rhs at input
   ENTRY_1DARRAY & y  // y is 0 at input 
  )
  {

    serial_btf_solve_tr(x,y);

    return 0;
  }

  // This results in the "solution" of the diagonal block transpose solve - input Matrix interpretted as "upper triangular" via transpose
  // Called after upper_tri_solve_tr
  // x should be mod. rhs, though x == y at input in range of M
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::lower_tri_solve_tr
  (
   BASKER_MATRIX &M,
   ENTRY_1DARRAY &x, // mod rhs at input; output to be garbage
   ENTRY_1DARRAY &y, // y != x at input likely (U^T likely not unit diag); output to be soln
   Int offset
  )
  {
    // block diagonal offset stays the same as non-transpose
    const Int bcol = M.scol + offset;
    const Int brow = M.scol + offset;

    // k is a col of L CCS; for transpose solve, treat k as a row of L^T
    for(Int k = M.ncol-1; k >= 0; --k)
    {
      //const Int istart = M.col_ptr(k)+1; // start one after diagonal; pivot solved after removing all non-diag contributions from rhs
      // these will be the offsets for col k, which we are treating as transposed matrix row k
      const Int istart = M.col_ptr(k);
      const Int iend = M.col_ptr(k+1);
      // Start 1 past the diagonal, skip for row updates before final pivot; 
      // For L, the diagonal is the first entry of the column (or transposed row), as opposed to being final entry with U
      // For the first iteration, this loop need not be entered because the only entry is the diagonal
      for(Int i = istart+1; i < iend; ++i)
      {
        // global row id for dense vectors
        const Int j = (Options.no_pivot == BASKER_FALSE) ? 
                        gperm(M.row_idx(i)+brow) :
                             (M.row_idx(i)+brow) ;

        x(k+brow) -= M.val(i)*y(j);
      } //over all nnz in a column
      // Complete the solution and store in rhs x 
      y(k+brow) = x(k+bcol) / M.val(M.col_ptr(k));
    } //over each column

    return 0;
  } //end lower_tri_solve_tr


  // Input matrix interpretted as "lower tri" matrix via transpose
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::upper_tri_solve_tr
  (
   BASKER_MATRIX &M, // U^T i.e. lower tri matrix with diagonal as last entry or indices pointer
   ENTRY_1DARRAY &x,
   ENTRY_1DARRAY &y,
   Int offset
  )
  {
    // block diagonal offset stays the same as non-transpose
    const Int bcol = M.scol + offset;
    const Int brow = M.srow + offset;

    y(brow) = x(bcol) / M.val(M.col_ptr(1)-1);

    // k is a col of U CCS; for transpose solve, treat k as a row of U^T
    for(Int k = 1; k < M.ncol; k++) // k == 0 already handled above
    {
      const Int istart = M.col_ptr(k);
      const Int iend = M.col_ptr(k+1);

      // skip the diagonal during row updates
      // for U, the diagonal should be stored as last entry of a column (or row for U^T) (not first, like L)
      for(Int i = istart; i < iend-1; ++i)
      {
        const Int j = M.row_idx(i) + brow;
        x(k+brow) -= M.val(i)*y(j);
      }
      // finish the diag 
      y(k+brow) = x(k+bcol) / M.val(iend-1); // y == x in M range assumed true at end of this routine, but not automatic as with non-transpose lower_tri_solve since U^T diagonal is not necessarily 1's
      x(k+bcol) = y(k+brow); // enforce x == y at end to avoid issues with lower_tri_solve_tr
    }//end over all columns

    x(bcol) = y(brow);  // set k == 0 values equal after updates complete

    return 0;
  } //end upper_tri_solve_tr


  // non-transpose A*x=b CCS
  // x1*col1(A) + ... + xn*coln(A) = b

  // transpose spmv
  // row1(AT) CRS == col1(A) CCS
  // transpose A^T*x=b, interpretting A^T as CRS to reuse A pointers since A transpose is not computed nor stored;
  // row1(A^T).*x + row2(A^T).*x ... = b, i.e.
  // col1(A).*x + col2(A).*x ... = b
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::spmv_BTF_tr
  (
   Int tab, // treat this as info to determine where the spmv terminates (i.e. final col of transpose matrix)
   BASKER_MATRIX &M,
   ENTRY_1DARRAY &x, // modified rhs
   ENTRY_1DARRAY &y, // intermediate solution
   bool full
  )
  {
    // What identifying block info needed before starting transpose spmv rhs update?
    //Tab = block in    
    Int bcol = btf_tabs(tab)- M.scol;
    Int ecol = btf_tabs(tab+1) - M.scol;
    if (ecol > M.ncol) {
      // for D block, btf_tabs(tab+1) > ncol.
      ecol = M.ncol;
    }

    //loop over each row of transpose (column of original)
    for(Int k = bcol; k < ecol; ++k) {
      const Int istart = M.col_ptr(k);
      const Int iend = M.col_ptr(k+1);
      const Int gk = k+M.scol;

      for(Int i = istart; i < iend; ++i) {
        const Int gj = (Options.no_pivot == BASKER_FALSE) ? 
                       gperm(M.row_idx(i) + M.srow) :
                       (M.row_idx(i) + M.srow) ;

        if (full || gj < bcol+M.scol) // bcol will act as a "starting bcol offset" for portion of the matrix to skip in the spmv^T update; we use the local col id j (row id of non-transpose)
        {
          x(gk) -= M.val(i)*y(gj);
        }
      } //over all nnz in row
    } // end for over col

    return 0;
  } //end spmv_BTF_tr;

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::neg_spmv_tr
  (
   BASKER_MATRIX &M,
   ENTRY_1DARRAY x, 
   ENTRY_1DARRAY y,
   Int offset
  )
  {
    const Int bcol  = M.scol + offset;
    const Int msrow = M.srow + offset;
    for(Int k=0; k < M.ncol; ++k)
    {
      const Int istart = M.col_ptr(k);
      const Int iend = M.col_ptr(k+1);

      const Int gk = k+bcol;
      for(Int i = istart; i < iend; ++i)
      {
        const Int gj = M.row_idx(i) + msrow;
        x(gk) -= M.val(i)*y(gj);
      }
    }

    return 0;
  } //neg_spmv_tr

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::neg_spmv_perm_tr
  (
   BASKER_MATRIX &M,
   ENTRY_1DARRAY &y, // input
   ENTRY_1DARRAY &x, // output
   Int offset
  )
  {
    const Int bcol  = M.scol + offset;
    const Int msrow = M.srow + offset;

    for(Int k=0; k < M.ncol; ++k)
    {
      const Int istart = M.col_ptr(k);
      const Int iend   = M.col_ptr(k+1);

      const Int gk = k+bcol;
      for(Int i = istart; i < iend; ++i)
      {
        const Int gj = (Options.no_pivot == BASKER_FALSE) ? 
                       gperm(M.row_idx(i) + msrow) :
                       (M.row_idx(i) + msrow) ;
        x(gk) -= M.val(i)*y(gj);
      }
    }

    return 0;
  } //neg_spmv_perm_tr


  // solve L^T*x=y - transpose means pattern of "backward" solve
  // this is the final stage of the solve
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::l_tran_brfa_solve
  (
   ENTRY_1DARRAY & y, // in: partial solution 
   ENTRY_1DARRAY & x  // out: solution in btfa range
  )
  {
    Int scol_top = btf_tabs[btf_top_tabs_offset]; // the first column index of A
    for(int b = tree.nblks-1; b >= 0; b--)
    {
      // Update off-diag in the block-row before the diag solve
      for(int bb = LL_size(b)-1; bb > 0; bb--)
      {
        BASKER_MATRIX &LD = LL(b)(bb);
        neg_spmv_perm_tr(LD, x, y, scol_top); // update y as mod. rhs, x as solution
      }
      BASKER_MATRIX &L = LL(b)(0);
      if (L.nrow != 0 && L.ncol != 0) // Avoid degenerate case e.g. empty block following nd-partitioning
        lower_tri_solve_tr(L, y, x, scol_top); // x and y should be equal after in M range...
    }

    return 0;
  }//end l_tran_brfa_solve()


  // U^T*y = x, transpose implies "forward" solve pattern
  template<class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::u_tran_btfa_solve
  (
   ENTRY_1DARRAY & x,
   ENTRY_1DARRAY & y
  )
  {
    Int scol_top = btf_tabs[btf_top_tabs_offset]; // the first column index of A
    for(Int b = 0; b < tree.nblks; b++)
    {
      for(Int bb = 0; bb <  LU_size(b)-1; bb++)
      {
        // update offdiag corresponding to the block-row
        BASKER_MATRIX &UB = LU(b)(bb);
        neg_spmv_tr(UB, x, y, scol_top);
      }
      BASKER_MATRIX &U = LU(b)(LU_size(b)-1);
      if (U.nrow != 0 && U.ncol != 0) // Avoid degenerate case
        upper_tri_solve_tr(U, x, y, scol_top);
    }


    if (BTF_A.ncol > 0) {
      for (Int i = 0; i < BTF_A.ncol; i++) {
        x(scol_top+i) = y(scol_top+i);
      }
    }

    return 0;
  }//end u_tran_btfa_solve()


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::serial_btf_solve_tr
  (
   ENTRY_1DARRAY & x, // Permuted rhs at input
   ENTRY_1DARRAY & y  // 0 at input
  )
  {
    // P1 T
    // BTF_D T solve
    //  Input : X
    //  Output: Y
    int offset = 0;
if (Options.verbose) std::cout << "BTF_D^T begin: from 0 to " << btf_top_tabs_offset << std::endl;
    if(btf_top_tabs_offset >  0)
    {
      for(Int b = 0; b < btf_top_tabs_offset; b++)
      {
        BASKER_MATRIX &UC = U_D(b);

        if ( b > 0 )
          spmv_BTF_tr(b, BTF_D, x, y, false);

        if (UC.nrow != 0 && UC.ncol != 0) // Avoid degenerate case
          upper_tri_solve_tr(UC, x, y);

        BASKER_MATRIX &LC = L_D(b);

        if (LC.nrow != 0 && LC.ncol != 0) // Avoid degenerate case
          lower_tri_solve_tr(LC, x, y);

        offset += UC.nrow;
      }
      // Checkpoint: x is mod rhs (and has some garbage in BTF_D range), y stores solution from BTF_D range

      // Update for offdiag BTF_E T
      //  Update X with Y
      if (Options.verbose) std::cout << "BTF_E^T update begin:" << std::endl;

      neg_spmv_perm_tr(BTF_E, y, x);

      // Checkpoint: x is mod. rhs (and has some garbage in BTF_D range), y stores solution from BTF_D range
      Int srow_d = BTF_D.srow;
      Int erow_d = BTF_D.nrow + srow_d;
      for (Int i = srow_d; i < erow_d; i++) {
        x(i) = y(i);
      }
    }


    // P2 T
    // BTF_A T solves
    if (BTF_A.nrow > 0 || BTF_A.ncol > 0) {
      u_tran_btfa_solve(x,y); // U^T*y=x
      l_tran_brfa_solve(y,x); // L^T*x=y
      // Checkpoint: in BTF_A range, y is mod. rhs (and has some garbage in BTF_A range), x stores solution
    }

    // Update for offdiag BTF_B T
    if(btf_tabs_offset !=  0)
    {
      neg_spmv_perm_tr(BTF_B,x,x); // x is updated rhs in BTF_C range, solution in BTF_D and BTF_A range
    }

    // P3 T
    Int nblks_c = btf_nblks-btf_tabs_offset;

    if (nblks_c > 0) {
      Int offset = 0;
      for(Int b = 0;  b < nblks_c; b++) {
        BASKER_MATRIX &UC = UBTF(b);

        // Update off-diag
        //  Update X with Y
        if ( b > 0 )
          spmv_BTF_tr(b+btf_tabs_offset, BTF_C, x, y, false);

        if (UC.nrow != 0 && UC.ncol != 0) // Avoid degenerate case
          upper_tri_solve_tr(UC,x,y);

        BASKER_MATRIX &LC = LBTF(b);

        if (LC.nrow != 0 && LC.ncol != 0) // Avoid degenerate case
          lower_tri_solve_tr(LC,x,y);

        offset += UC.nrow;
      }

      // copy the solution for C from y to x
      Int srow_c = BTF_C.srow;
      Int erow_c = BTF_C.nrow + srow_c;
      for (Int i = srow_c; i < erow_c; i++) {
        x(i) = y(i);
      }
    }
    return 0;
  }

} // namespace BaskerNS

#endif
