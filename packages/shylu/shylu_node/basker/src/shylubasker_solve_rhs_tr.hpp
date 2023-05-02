#ifndef SHYLUBASKER_SOLVE_RHS_TR_HPP
#define SHYLUBASKER_SOLVE_RHS_TR_HPP

/*Basker Includes*/
//#include "shylubasker_decl.hpp"
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

//#define BASKER_DEBUG_SOLVE_RHS_TR
//#define BASKER_DEBUG_SOLVE_RHS_TR_PRINT_MTX
//#define BASKER_DEBUG_SOLVE_RHS_TR_PRINT_BLOCKS


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
#ifdef BASKER_DEBUG_SOLVE_RHS_TR
    printf("---- Pre-Solve printOptions: ----\n");
    // TODO: Add other permutation options
    //Options.printOptions();
    printf( " >> gperm; print non-identity gperm output\n" );
    for (Int i = 0; i < gn; i++) {
      if (gperm(i) != i) {
        printf( "  >> gperm(%d) = %d\n", (int)i, (int)gperm(i) );
      }
    }
    printf( " >> gpermi; print non-identity gpermi output\n" );
    for (Int i = 0; i < gn; i++) {
      if (gpermi(i) != i) {
        printf( "  >> gpermi(%d) = %d\n", (int)i, (int)gpermi(i) );
      }
    }
#endif

    // Transpose: Swap permutation order - only handles case without 
    // TODO: determine which case+options the below routine supports - i.e. WHAT to enable and WHAT to disable
#if 1//def BASKER_DEBUG_SOLVE_RHS_TR
    //for (Int i = 0; i < gn; i++) printf( "  perm_comp_array(%d) = %d\n",i,perm_comp_array[i] );
    //for (Int i = 0; i < gn; i++) printf( "  (unused) perm_inv_comp_array(%d) = %d\n",i,perm_inv_comp_array[i] );
    if (Options.verbose) {
      printf( "perm=[\n");
      for (Int i = 0; i < gn; i++) printf( " %d %d %d\n",i,perm_inv_comp_array[i],perm_comp_array[i] );
      printf( "];\n");
      printf("b_0=[\n");
      for (Int i = 0; i < gn; i++) printf("%d %.16e\n",i,_y[i]);
      printf("];\n");
    }
#endif

    if (Options.blk_matching != 0) {
if (Options.verbose)
printf( " +++ Matching +++\n" );
      for(Int i = 0; i < gn; i++) {
        Int col = symbolic_col_iperm_array(i);
if (Options.verbose && scale_col_array(col) != 1.0) printf( "scal_col(%d) = %e\n",col,scale_col_array(col) );
        _y[i] = scale_col_array(col) * _y[i];
      }
    }
    // rhs content from _y has been permuted and copied to x_view_ptr_copy, and will act as the rhs-to-update during solve calls; 
    // y_view_ptr_copy is initialized to be zeros, and will store the pivots (i.e. solutions)
    permute_inv_and_init_for_solve(_y, x_view_ptr_copy, y_view_ptr_copy, perm_comp_array, gn);

    if (Options.matrix_scaling != 0) {
if (Options.verbose) 
printf( " +++ Scaling +++\n");
      for (Int i = 0; i < gn; i++) {
        x_view_ptr_copy(i) = x_view_ptr_copy(i) * scale_col_array(i);
      }
    }
    if (Options.blk_matching != 0 || Options.static_delayed_pivot != 0) {
if (Options.verbose) 
{
  printf( " +++ Numeric +++\n");
  for (Int i=0; i < gn; i++) if (numeric_col_iperm_array(i) != i) printf( " numeric_col_iperm_array(%d) = %d\n",i,numeric_col_iperm_array(i) );
}
      permute_inv_with_workspace(x_view_ptr_copy, numeric_col_iperm_array, gn);
    }
if (Options.verbose) {
  printf("b=[\n");
  for (Int i = 0; i < gn; i++) printf("%d %.16e\n",i,x_view_ptr_copy(i));
  printf("];\n");
}
    // ==================
    // transpose-solve
    // after the solve: solution is stored in [x_view_ptr_copy(1:poffset); y_view_ptr_copy(1:(gn-poffset))]
    solve_interfacetr(x_view_ptr_copy, y_view_ptr_copy);
if (Options.verbose) {
  printf("x=[\n");
  for (Int i = 0; i < gn; i++) printf("%.16e %.16e\n",x_view_ptr_copy(i),y_view_ptr_copy(i));
  printf("];\n");
}

    if (Options.no_pivot == BASKER_FALSE) {
if (Options.verbose) 
printf( " +++ Pivot +++\n" );
      // apply partial pivoting from numeric
      //for (Int i = 0; i < gn; i++) printf( " gperm(%d) = %d\n",i,gperm(i) );
      //std::cout << "  Permute no_pivot == false: permute_with_workspace" << std::endl;
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

#ifdef BASKER_DEBUG_SOLVE_RHS_TR
    for (Int i = 0; i < gn; i++) printf( "  perm_inv_comp_array(%d) = %d\n",i,perm_inv_comp_array[i] );
#endif

    if (Options.blk_matching != 0 || Options.static_delayed_pivot != 0) {
if (Options.verbose) {
  printf( " +++ Matching (out) +++\n" );
  for (Int i=0; i < gn; i++) if (numeric_row_iperm_array(i) != i) printf( " numeric_row_iperm_array(%d) = %d\n",i,numeric_row_iperm_array(i) );
}
      permute_inv_and_finalcopy_after_solve(_x, x_view_ptr_copy, y_view_ptr_copy, numeric_row_iperm_array, gn);
      permute_array_with_workspace(_x, perm_inv_comp_array, gn);

      for(Int i = 0; i < gn; i++) {
        //Int row = order_blk_mwm_array(symbolic_row_iperm_array(i));
        Int row = symbolic_row_perm_array(order_blk_mwm_inv(i));
        _x[i] = scale_row_array(row) * _x[i];
if (Options.verbose && scale_row_array(row) != 1.0) printf( "scal(%d) = %e\n",row,scale_row_array(row) );
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
if (Options.verbose) {
  printf("x_0=[\n");
  for (Int i = 0; i < gn; i++) printf("%.16e\n",_x[i]);
  printf("];\n");
}

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

    // TODO Add other options
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

#if 1//def BASKER_DEBUG_SOLVE_RHS_TR_PRINT_MTX
if (Options.verbose) {
    std::cout << "  L^T begin: bcol = " << bcol << "  brow = " << brow << "  offset = " << offset << std::endl;
    printf("b=[\n");
    for(Int k = 0; k < M.ncol; ++k) printf("%.16e\n",x(k+brow));
    printf("];\n");
    printf("L=[\n");
    for(Int k = 0; k < M.ncol; ++k) {
      for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); ++i)
      {
        printf("%d %d %.16e\n",M.row_idx(i),k,M.val(i));
      }
    }
    printf("];\n");
    //if (M.nrow != 0 && M.ncol != 0)
    //  M.print();
}
#endif
    // k is a col of L CCS; for transpose solve, treat k as a row of L^T
    for(Int k = M.ncol-1; k >= 0; --k)
    {
#ifdef BASKER_DEBUG_SOLVE_RHS_TR
      //Test if zero pivot value
      BASKER_ASSERT(M.val[M.col_ptr[k]]!=0.0, "LOWER PIVOT 0");

      std::cout << "  LT: k = " << k << "  bcol = " << bcol << "  brow = " << brow << std::endl;
#endif

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

#ifdef BASKER_DEBUG_SOLVE_RHS_TR
        BASKER_ASSERT(j != BASKER_MAX_IDX,"Using nonperm\n");
#endif

#ifdef BASKER_DEBUG_SOLVE_RHS_TR
        std::cout << "    inner loop: i = " << i << "  j = " << j << "  brow = " << brow << std::endl;
        std::cout << "    Before update to k+brow: x(k+brow) = " << x(k+brow) << "  y(j) = " << y(j) << "  M.val(i) = " << M.val(i) << std::endl;
#endif
        x(k+brow) -= M.val(i)*y(j);
#ifdef BASKER_DEBUG_SOLVE_RHS_TR
        std::cout << "    After update to k+brow: x(k+brow) = " << x(k+brow) << std::endl;
#endif
      } //over all nnz in a column
      // Complete the solution and store in rhs x 
#ifdef BASKER_DEBUG_SOLVE_RHS_TR
      std::cout << "  LT Pre-row k solve: y(k+brow) = " << y(k+brow) << " x(k+bcol) = " << x(k+bcol) << " M.val(istart) (diag entry) = " << M.val(istart) << std::endl;
#endif
      y(k+brow) = x(k+bcol) / M.val(M.col_ptr(k));
#ifdef BASKER_DEBUG_SOLVE_RHS_TR
      std::cout << "  After row k solve: y(k+brow) = " << y(k+brow) << std::endl;
#endif
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

#if 1//def BASKER_DEBUG_SOLVE_RHS_TR_PRINT_MTX
if (Options.verbose) {
    std::cout << "  U^T begin: bcol = " << bcol << "  brow = " << brow << "  offset = " << offset << std::endl;
    printf("b=[\n");
    for(Int k = 0; k < M.ncol; ++k) printf("%.16e\n",x(bcol+k));
    printf("];\n");
    printf("U=[\n");
    for(Int k = 0; k < M.ncol; ++k) {
      for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); ++i)
      {
        printf("%d %d %.16e\n",M.row_idx(i),k,M.val(i));
      }
    }
    printf("];\n");
    //if (M.nrow != 0 && M.ncol != 0)
    //  M.print();
}
#endif
#if defined(BASKER_DEBUG_SOLVE_RHS_TR_PRINT_BLOCKS) || defined(BASKER_DEBUG_SOLVE_RHS_TR)
    // Solve initial row 0 transpose (i.e. col 0) before inner loop
    std::cout << "\n  UT Pre-row 0 solve: y(brow) = " << y(brow) << " x(bcol) = " << x(bcol) << " M.val(M.col_ptr(1)-1) (diag entry) = " << M.val(M.col_ptr(1)-1) << std::endl;
#endif
    y(brow) = x(bcol) / M.val(M.col_ptr(1)-1);
#if defined(BASKER_DEBUG_SOLVE_RHS_TR_PRINT_BLOCKS) || defined(BASKER_DEBUG_SOLVE_RHS_TR)
    std::cout << "  After row 0 solve: y(brow) = " << y(brow) << std::endl;
#endif

    // k is a col of U CCS; for transpose solve, treat k as a row of U^T
    for(Int k = 1; k < M.ncol; k++) // k == 0 already handled above
    {
#if defined(BASKER_DEBUG_SOLVE_RHS_TR_PRINT_BLOCKS) || defined(BASKER_DEBUG_SOLVE_RHS_TR)
      BASKER_ASSERT(M.val[M.col_ptr[k]-1]!=0.0, "UPPER PIVOT == 0\n");
        /*
        printf("Upper Tri Solve, scol: %d ncol: %d \n",
          M.scol, M.ncol);
        */
      std::cout << "  UT: k = " << k << "  bcol = " << bcol << "  brow = " << brow << std::endl;
#endif

      const Int istart = M.col_ptr(k);
      const Int iend = M.col_ptr(k+1);

      // skip the diagonal during row updates
      // for U, the diagonal should be stored as last entry of a column (or row for U^T) (not first, like L)
      for(Int i = istart; i < iend-1; ++i)
      {
        const Int j = M.row_idx(i) + brow;
#if defined(BASKER_DEBUG_SOLVE_RHS_TR_PRINT_BLOCKS) || defined(BASKER_DEBUG_SOLVE_RHS_TR)
        std::cout << "    inner loop: i = " << i << "  j = " << j << "  brow = " << brow << std::endl;
        std::cout << "    Before update to k+brow: x(k+brow) = " << x(k+brow) << "  y(j) = " << y(j) << "  M.val(i) = " << M.val(i) << std::endl;
#endif
        x(k+brow) -= M.val(i)*y(j);
#if defined(BASKER_DEBUG_SOLVE_RHS_TR_PRINT_BLOCKS) || defined(BASKER_DEBUG_SOLVE_RHS_TR)
        std::cout << "    After update to k+brow: x(k+brow) = " << x(k+brow) << std::endl;
#endif
      }
      // finish the diag 
#if defined(BASKER_DEBUG_SOLVE_RHS_TR_PRINT_BLOCKS) || defined(BASKER_DEBUG_SOLVE_RHS_TR)
      std::cout << "  UT Pre-row k solve: y(k+brow) = " << y(k+brow) << " x(k+bcol) = " << x(k+bcol) << " M.val(iend-1) (diag entry) = " << M.val(iend-1) << std::endl;
#endif
      y(k+brow) = x(k+bcol) / M.val(iend-1); // y == x in M range assumed true at end of this routine, but not automatic as with non-transpose lower_tri_solve since U^T diagonal is not necessarily 1's
#if defined(BASKER_DEBUG_SOLVE_RHS_TR_PRINT_BLOCKS) || defined(BASKER_DEBUG_SOLVE_RHS_TR)
      std::cout << "  After row k solve: y(k+brow) = " << y(k+brow) << std::endl;
      std::cout << "    about to update x: k+bcol = " << k+bcol << " k+brow = " << k+brow << std::endl;
#endif
      x(k+bcol) = y(k+brow); // enforce x == y at end to avoid issues with lower_tri_solve_tr
#if defined(BASKER_DEBUG_SOLVE_RHS_TR_PRINT_BLOCKS) || defined(BASKER_DEBUG_SOLVE_RHS_TR)
      std::cout << "    x update: x(k+bcol) = " << x(k+bcol) << "  y(k+brow) = " << y(k+brow) << std::endl;
#endif

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
#ifdef BASKER_DEBUG_SOLVE_RHS_TR_PRINT_MTX
    if (M.nrow != 0 && M.ncol != 0)
      M.print();
#endif
    Int bcol = btf_tabs(tab)- M.scol;
    Int ecol = btf_tabs(tab+1) - M.scol;
    if (ecol > M.ncol) {
      // for D block, btf_tabs(tab+1) > ncol.
      ecol = M.ncol;
    }

#if defined(BASKER_DEBUG_SOLVE_RHS_TR_PRINT_BLOCKS) || defined(BASKER_DEBUG_SOLVE_RHS_TR)
    Int brow = M.srow;
    std::cout << "\n  spmv_BTF_tr tab = " << tab << "  M.scol = " << M.scol << "  bcol = " << bcol << "  ecol = " << ecol << std::endl;
    Int erow = 0;
    if(tab > 0)
    {
      erow = btf_tabs(tab);
    }
    else
    {
      erow = brow-1;
    }
    printf("  BTF_UPDATE, TAB: %d [%d %d] [%d %d] \n",
        tab, brow, erow, bcol, ecol);
#endif

    //loop over each row of transpose (column of original)
if (Options.verbose) {
  printf( " SPMV(%d:%d), tab=%d(%d), sroc=%d,scol=%d\n",bcol,ecol-1,btf_tabs(tab),tab,M.srow,M.scol );
  /*printf("D=[\n");
  for(Int k = bcol; k < ecol; ++k) {
    const Int istart = M.col_ptr(k);
    const Int iend = M.col_ptr(k+1);
    for(Int i = istart+1; i < iend; ++i)
    {
      printf("%d %d %.16e\n",M.row_idx(i),k,M.val(i));
    }
  }
  printf("];\n");*/
}
    for(Int k = bcol; k < ecol; ++k) {
      const Int istart = M.col_ptr(k);
      const Int iend = M.col_ptr(k+1);
#if defined(BASKER_DEBUG_SOLVE_RHS_TR_PRINT_BLOCKS) || defined(BASKER_DEBUG_SOLVE_RHS_TR)
      std::cout << "  spmv_BTF_tr: k = " << k << "  istart = " << istart << "  iend = " << iend << std::endl;
#endif

      const Int gk = k+M.scol;

      for(Int i = istart; i < iend; ++i) {
        const Int gj = (Options.no_pivot == BASKER_FALSE) ? 
                       gperm(M.row_idx(i) + M.srow) :
                       (M.row_idx(i) + M.srow) ;

#if defined(BASKER_DEBUG_SOLVE_RHS_TR_PRINT_BLOCKS) || defined(BASKER_DEBUG_SOLVE_RHS_TR)
        const auto j = M.row_idx(i);
        printf("   BTF_UPDATE-val, j: %d x: %f y: %f, val: %f \n",
            gj, x[gj], y[k+M.scol], M.val[i]);
        std::cout << "    inner loop : gk = " << gk << "  j = " << j << "  gj = " << gj << "  bcol = " << bcol << std::endl;
#endif

        if (full || gj < bcol+M.scol) // bcol will act as a "starting bcol offset" for portion of the matrix to skip in the spmv^T update; we use the local col id j (row id of non-transpose)
        {
#if defined(BASKER_DEBUG_SOLVE_RHS_TR_PRINT_BLOCKS) || defined(BASKER_DEBUG_SOLVE_RHS_TR)
          std::cout << "     Pre-update j < bcol: x(gk) = " << x(gk) << "  y(gj) = " << y(gj) << "  M.val(i) = " << M.val(i) << std::endl;
#endif
if (Options.verbose) printf("x(%d) (%d -> %d) : %e - %e * %e = %e (D(%d,%d)=%e)\n",gk,M.row_idx(i),gj, x(gk),M.val(i),y(gj),x(gk)-M.val(i)*y(gj), M.row_idx(i),k,M.val(i));
          x(gk) -= M.val(i)*y(gj);
#if defined(BASKER_DEBUG_SOLVE_RHS_TR_PRINT_BLOCKS) || defined(BASKER_DEBUG_SOLVE_RHS_TR)
          std::cout << "     After update j < bcol: x(gk) = " << x(gk) << std::endl;
#endif
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
#ifdef BASKER_DEBUG_SOLVE_RHS_TR
    printf("neg_spmv_tr. scol: %d ncol: %d \n", M.scol, M.ncol);
#endif

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
        //y(gk) -= M.val(i)*x(gj);
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
#ifdef BASKER_DEBUG_SOLVE_RHS_TR
    printf("neg_spmv_perm_tr. scol: %d ncol: %d, srow: %d nrow: %d \n", M.scol, M.ncol, M.srow, M.nrow);
    if (Options.no_pivot == BASKER_FALSE) {
      printf("P=[\n");
      for(Int k=0; k < M.nrow; ++k) printf("%d+%d, %d\n",msrow,k,gperm(k+msrow));
      printf("];\n");
    }
    printf("M=[\n");
    for(Int k=0; k < M.ncol; ++k) {
      for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); ++i)
        printf( "%d %d %d %.16e\n",gperm(M.row_idx(i) + msrow)-msrow, M.row_idx(i),k,M.val(i) );
    }
    printf("];\n");
#endif

    for(Int k=0; k < M.ncol; ++k)
    {
      const Int istart = M.col_ptr(k);
      const Int iend   = M.col_ptr(k+1);

      const Int gk = k+bcol;
      for(Int i = istart; i < iend; ++i) //NDE retest with const vars, scope tightly
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
#ifdef BASKER_DEBUG_SOLVE_RHS_TR
    printf("Called l_tran_brfa_solve\n");

    printf(" Vector y before l_tran_btfa_solve \n");
    printVec(y, gn);
    printf(" Vector x before l_tran_btfa_solve \n");
    printVec(x, gn);
#endif

    Int scol_top = btf_tabs[btf_top_tabs_offset]; // the first column index of A
    for(int b = tree.nblks-1; b >= 0; b--)
    {
#ifdef BASKER_DEBUG_SOLVE_RHS_TR
      std::cout << "\nLT: Block Row b = " << b << "\n\n" << std::endl;
      std::cout << "    LL_size(b) = " << LL_size(b) << std::endl;
#endif
      // Update off-diag in the block-row before the diag solve
      for(int bb = LL_size(b)-1; bb > 0; bb--)
      {
        BASKER_MATRIX &LD = LL(b)(bb);
#ifdef BASKER_DEBUG_SOLVE_RHS_TR_PRINT_MTX
        std::cout << "    LT rhs update (BTF_A). bb = " << bb << std::endl;
        printf("LT update blk (%d, %d): size=(%dx%d) srow=%d, scol=%d\n",b,bb, (int)LD.nrow,(int)LD.ncol, (int)LD.srow,(int)LD.scol);
        if (LD.nrow != 0 && LD.ncol != 0)
          LD.print();
#endif
        neg_spmv_perm_tr(LD, x, y, scol_top); // update y as mod. rhs, x as solution
#ifdef BASKER_DEBUG_SOLVE_RHS_TR
        printf(" Vector y after neg_spmv_perm_tr \n");
        printVec(y, gn);
        printf(" Vector x after neg_spmv_perm_tr \n");
        printVec(x, gn);
#endif
      }
      BASKER_MATRIX &L = LL(b)(0);
#ifdef BASKER_DEBUG_SOLVE_RHS_TR
      std::cout << "  LT solve (BTF_A). b = " << b << std::endl;
#endif
      if (L.nrow != 0 && L.ncol != 0) // Avoid degenerate case e.g. empty block following nd-partitioning
        lower_tri_solve_tr(L, y, x, scol_top); // x and y should be equal after in M range...
#ifdef BASKER_DEBUG_SOLVE_RHS_TR
        printf("LT diag blk (%d, %d): size=(%dx%d) srow=%d, scol=%d\n",b,0, (int)L.nrow,(int)L.ncol, (int)L.srow,(int)L.scol);
        printf(" Vector y after lower_tri_solve_tr \n");
        printVec(y, gn);
        printf(" Vector x after lower_tri_solve_tr \n");
        printVec(x, gn);
#endif
    }

#ifdef BASKER_DEBUG_SOLVE_RHS_TR
        printf(" Vector y after l_tran_btfa_solve \n");
        printVec(y, gn);
        printf(" Vector x after l_tran_btfa_solve \n");
        printVec(x, gn);
#endif

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
#ifdef BASKER_DEBUG_SOLVE_RHS_TR
    printf("called u_tran_btfa_solve\n");
    printf(" Vector y before u_tran_btfa_solve \n");
    printVec(y, gn);
    printf(" Vector x before u_tran_btfa_solve \n");
    printVec(x, gn);
#endif

    Int scol_top = btf_tabs[btf_top_tabs_offset]; // the first column index of A
    for(Int b = 0; b < tree.nblks; b++)
    {
#ifdef BASKER_DEBUG_SOLVE_RHS_TR
      std::cout << "\nUT: Block Row b = " << b << "\n\n" << std::endl;
      std::cout << "      LU_size(b) = " << LU_size(b) << std::endl;
#endif
      for(Int bb = 0; bb <  LU_size(b)-1; bb++)
      {
        // update offdiag corresponding to the block-row
        BASKER_MATRIX &UB = LU(b)(bb);
#ifdef BASKER_DEBUG_SOLVE_RHS_TR_PRINT_MTX
        std::cout << "    UT update rhs (BTF_A). bb = " << bb << std::endl;
        printf("UT update blk (%d, %d): size=(%dx%d) srow=%d, scol=%d\n",b,bb, (int)UB.nrow,(int)UB.ncol, (int)UB.srow,(int)UB.scol);
        if (UB.nrow != 0 && UB.ncol != 0)
          UB.print();
#endif
        neg_spmv_tr(UB, x, y, scol_top);
#ifdef BASKER_DEBUG_SOLVE_RHS_TR
        printf(" Vector y after neg_spmv_tr \n");
        printVec(y, gn);
        printf(" Vector x after neg_spmv_tr \n");
        printVec(x, gn);
#endif
      }
      BASKER_MATRIX &U = LU(b)(LU_size(b)-1);
#ifdef BASKER_DEBUG_SOLVE_RHS_TR
      std::cout << "    UT solve (BTF_A). b = " << b << std::endl;
#endif
      if (U.nrow != 0 && U.ncol != 0) // Avoid degenerate case
        upper_tri_solve_tr(U, x, y, scol_top);
#ifdef BASKER_DEBUG_SOLVE_RHS_TR
        printf("UT diag blk (%d, %d): size=(%dx%d) srow=%d, scol=%d\n",b,LU_size(b)-1, (int)U.nrow,(int)U.ncol, (int)U.srow,(int)U.scol);
        printf("  Right after upper_tri_solve_tr call\n");
        printf(" Vector y after upper_tri_solve_tr \n");
        printVec(y, gn);
        printf(" Vector x after upper_tri_solve_tr \n");
        printVec(x, gn);
#endif
    }


    if (BTF_A.ncol > 0) {
      for (Int i = 0; i < BTF_A.ncol; i++) {
        x(scol_top+i) = y(scol_top+i);
      }
    }

#ifdef BASKER_DEBUG_SOLVE_RHS_TR
    printf(" Vector y after u_tran_btfa_solve \n");
    printVec(y, gn);
    printf(" Vector x after u_tran_btfa_solve \n");
    printVec(x, gn);
#endif

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
#ifdef BASKER_DEBUG_SOLVE_RHS_TR
    std::cout << " BTF Stats:\n";
    std::cout << "    btf_nblks = " << btf_nblks << std::endl;
    std::cout << "    btf_tabs_offset = " << btf_tabs_offset << std::endl;
    std::cout << "    btf_top_tabs_offset = " << btf_top_tabs_offset << std::endl;
    std::cout << "    btf_top_nblks = " << btf_top_nblks << std::endl;

    printf(" Vector y (before solves)\n");
    printVec(y, gn);
    printf(" Vector x (before solves)\n");
    printVec(x, gn);
#endif

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
#if 1//def BASKER_DEBUG_SOLVE_RHS_TR
        if (Options.verbose) std::cout << "  diag block b = " << b << "(" << UC.nrow << "x" << UC.ncol << ")" << std::endl;
        //for (int i=offset; i<offset+UC.nrow; i++) printf("%d %.16e %.16e\n",i,x(i),y(i));
        //printf(" Vector y after spmv_BTF_tr\n");
        //printVec(y, gn);
        //printf(" Vector x after spmv_BTF_tr\n");
        //printVec(x, gn);
#endif
        if (UC.nrow != 0 && UC.ncol != 0) // Avoid degenerate case
          upper_tri_solve_tr(UC, x, y);

#ifdef BASKER_DEBUG_SOLVE_RHS_TR
        printf( "\n After UT solve (b=%d): x(i),y(i)\n",b ); fflush(stdout);
        for (Int i=offset; i<offset+UC.nrow; i++) printf("%d: %e %e\n",i,x(i),y(i));
        //for (Int i = 0; i < gn; i++) printf( " %e %e\n",x(i),y(i));
        printf( "\n");
#endif

        BASKER_MATRIX &LC = L_D(b);
        if (LC.nrow != 0 && LC.ncol != 0) // Avoid degenerate case
          lower_tri_solve_tr(LC, x, y); // TODO Want x == mod. rhs, y == soln  after this op

#ifdef BASKER_DEBUG_SOLVE_RHS_TR
        printf( "\n After LT solve (b=%d): x(i),y(i)\n",b ); fflush(stdout);
        for (Int i=offset; i<offset+UC.nrow; i++) printf("%d: %e %e\n",i,x(i),y(i));
        //for (Int i = 0; i < gn; i++) printf( " %e %e\n",x(i),y(i));
        printf( "\n");
#endif
        offset += UC.nrow;
      }
      // Checkpoint: x is mod rhs (and has some garbage in BTF_D range), y stores solution from BTF_D range
#if 1//def BASKER_DEBUG_SOLVE_RHS_TR
      if (Options.verbose) {
        printf(" Vector y after BTF_D^T solve\n");
        printVec(y, gn);
        //printf(" Vector x after BTF_D^T solve\n");
        //printVec(x, gn);
      }
#endif

      // Update for offdiag BTF_E T
      //  Update X with Y
      if (Options.verbose) std::cout << "BTF_E^T update begin:" << std::endl;
      neg_spmv_perm_tr(BTF_E, y, x);
#ifdef BASKER_DEBUG_SOLVE_RHS_TR
      printf(" Vector y after neg_spmv_perm_tr\n");
      printVec(y, gn);
      printf(" Vector x after neg_spmv_perm_tr\n");
      printVec(x, gn);
#endif
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
#ifdef BASKER_DEBUG_SOLVE_RHS_TR_PRINT_MTX
      std::cout << "BTF_A matrix(" << BTF_A.nrow << "x" << BTF_A.ncol << ")" << std::endl;
      BTF_A.print();
#endif
#ifdef BASKER_DEBUG_SOLVE_RHS_TR
      std::cout << "BTF_A serial_backward_solve_tr" << std::endl;
#endif
      u_tran_btfa_solve(x,y); // U^T*y=x
                            //
#ifdef BASKER_DEBUG_SOLVE_RHS_TR
      std::cout << "BTF_A serial_forward_solve_tr" << std::endl;
#endif
      l_tran_brfa_solve(y,x); // L^T*x=y

      // Checkpoint: in BTF_A range, y is mod. rhs (and has some garbage in BTF_A range), x stores solution
#if 1//def BASKER_DEBUG_SOLVE_RHS_TR
      if (Options.verbose) {
        printf(" Vector y after BTF_A solve \n");
        printVec(y, gn);
        printf(" Vector x after BTF_A solve \n");
        printVec(x, gn);
      }
#endif
    }

    // Update for offdiag BTF_B T
if (Options.verbose) std::cout << "BTF_B^T update with offset=" << btf_tabs_offset << std::endl;
    if(btf_tabs_offset !=  0)
    {
      if (Options.verbose) std::cout << "BTF_B^T update begin:" << std::endl;
      neg_spmv_perm_tr(BTF_B,x,x); // x is updated rhs in BTF_C range, solution in BTF_D and BTF_A range
#ifdef BASKER_DEBUG_SOLVE_RHS_TR
      printf(" Vector x after BTF_B update \n");
      printVec(x, gn);
      printf(" Vector y after BTF_B update \n");
      printVec(y, gn);
#endif
    }

    // P3 T
    Int nblks_c = btf_nblks-btf_tabs_offset;
#ifdef BASKER_DEBUG_SOLVE_RHS_TR
    std::cout << "\nBTF_C^T region:" << std::endl;
#endif

#ifdef BASKER_DEBUG_SOLVE_RHS_TR_PRINT_BLOCKS
    if (BTF_C.nrow != 0 && BTF_C.ncol != 0) {
      //BTF_C.print();
      printMTX("BTF_C.mtx", BTF_C);
    }
#endif

#if 1//defined(BASKER_DEBUG_SOLVE_RHS_TR) || defined(BASKER_DEBUG_SOLVE_RHS_TR_PRINT_BLOCKS)
    if (Options.verbose) {
      std::cout << "\n\nBTF_C^T with " << nblks_c << " blocks" << std::endl;
      printf(" Vector x\n");
      printVec(x, gn);
      //printf(" Vector y\n");
      //printVec(y, gn);
      printf("(%dx%d)\n",BTF_D.nrow,BTF_D.ncol);
      printf("C=[\n");
      for (Int j=0; j<BTF_C.ncol; j++) {
        for(Int k=BTF_C.col_ptr(j); k<BTF_C.col_ptr(j+1); k++) {
          printf("%d %d %.16e\n",BTF_C.row_idx(k),j,BTF_C.val(k));
        }
      }
      printf("]\n");
    }
#endif
    if (nblks_c > 0) {
      Int offset = 0;
      for(Int b = 0;  b < nblks_c; b++) {
        BASKER_MATRIX &UC = UBTF(b);

#if 1//defined(BASKER_DEBUG_SOLVE_RHS_TR) || defined(BASKER_DEBUG_SOLVE_RHS_TR_PRINT_BLOCKS)
        if (Options.verbose) {
          std::cout << "\n\nBTF_C^T diag block b = " << b << "(" << UC.nrow << "x" << UC.ncol << ")" << std::endl;
          //for (Int i=0; i<offset; i++) printf("- %d(%d): %e %e\n",i,i,x(i),y(i));
          //for (Int i=offset; i<offset+UC.nrow; i++) printf("+ %d(%d): %e %e\n",i,i-offset,x(i),y(i));
          for (Int i=offset; i<offset+UC.nrow; i++) printf("%d %d %.16e %.16e\n",i,i-offset,x(i),y(i));
        }
#endif

        // Update off-diag
        //  Update X with Y
        if ( b > 0 )
          spmv_BTF_tr(b+btf_tabs_offset, BTF_C, x, y, false);

#if 1//def BASKER_DEBUG_SOLVE_RHS_TR
        if (Options.verbose) {
          printf(" Vector y print (after spmv_BTF_tr update )\n");
          //printVec(y, gn);
          //printf(" Vector x print (after spmv_BTF_tr update )\n");
          //printVec(x, gn);
          //for (Int i=0; i<offset; i++) printf("- %d(%d): %e %e\n",i,i,x(i),y(i));
          //for (Int i=offset; i<offset+UC.nrow; i++) printf("+ %d(%d): %e %e\n",i,i-offset,x(i),y(i));
          for (Int i=offset; i<offset+UC.nrow; i++) printf("%d %d %.16e %.16e\n",i,i-offset,x(i),y(i));
          printf("\n");
       }
#endif
#ifdef BASKER_DEBUG_SOLVE_RHS_TR_PRINT_MTX
        std::cout << "UC  b: " << b << std::endl;
        if (UC.nrow != 0 && UC.ncol != 0)
          UC.print();
#endif

        if (UC.nrow != 0 && UC.ncol != 0) // Avoid degenerate case
          upper_tri_solve_tr(UC,x,y);
#if 1//def BASKER_DEBUG_SOLVE_RHS_TR
        if (Options.verbose) {
          printf(" Vector y (after upper_tri_tr solves)\n");
          printf("U=[\n");
          for (Int j=0; j<UC.ncol; j++) {
            for(Int k=UC.col_ptr(j); k<UC.col_ptr(j+1); k++) {
              printf("%d %d %.16e\n",UC.row_idx(k),j,UC.val(k));
            }
          }
          printf("]\n");
          //printVec(y, gn);
          //printf(" Vector x (after upper_tri_tr solves)\n");
          //printVec(x, gn);
          //for (Int i=0; i<offset; i++) printf("- %d(%d): %e %e\n",i,i,x(i),y(i));
          //for (Int i=offset; i<offset+UC.nrow; i++) printf("+ %d(%d): %e %e\n",i,i-offset,x(i),y(i));
          for (Int i=offset; i<offset+UC.nrow; i++) printf("%d %d %.16e %.16e\n",i,i-offset,x(i),y(i));
          printf("\n");
        }
#endif

        BASKER_MATRIX &LC = LBTF(b);
#ifdef BASKER_DEBUG_SOLVE_RHS_TR_PRINT_MTX
        std::cout << "LC  b: " << b << std::endl;
        if (LC.nrow != 0 && LC.ncol != 0)
          LC.print();
#endif

        if (LC.nrow != 0 && LC.ncol != 0) // Avoid degenerate case
          lower_tri_solve_tr(LC,x,y);
#if 1//def BASKER_DEBUG_SOLVE_RHS_TR
        if (Options.verbose) {
          printf(" Vector y print (after lower_tri_tr solves)\n");
          printf("L=[\n");
          for (Int j=0; j<LC.ncol; j++) {
            for(Int k=LC.col_ptr(j); k<LC.col_ptr(j+1); k++) {
              printf("%d %d %.16e\n",LC.row_idx(k),j,LC.val(k));
            }
          }
          printf("]\n");
          //printVec(y, gn);
          //printf(" Vector x print (after lower_tri_tr solves)\n");
          //printVec(x, gn);
          //for (Int i=0; i<offset; i++) printf("- %d(%d): %e %e\n",i,i,x(i),y(i));
          //for (Int i=offset; i<offset+UC.nrow; i++) printf("+ %d(%d): %e %e\n",i,i-offset,x(i),y(i));
          for (Int i=offset; i<offset+UC.nrow; i++) printf("%d %d %.16e %.16e\n",i,i-offset,x(i),y(i));
          printf("\n");
        }
#endif
        offset += UC.nrow;
      }
#if 1//def BASKER_DEBUG_SOLVE_RHS_TR
      if (Options.verbose) {
        printf(" Vector y print (after C solves)\n");
        printVec(y, gn);
        //printf(" Vector x print (after C solves)\n");
        //printVec(x, gn);
      }
#endif

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
