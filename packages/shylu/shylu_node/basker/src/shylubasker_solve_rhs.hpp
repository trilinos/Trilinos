// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLUBASKER_SOLVE_RHS_HPP
#define SHYLUBASKER_SOLVE_RHS_HPP

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

//#define BASKER_DEBUG_SOLVE_RHS

using namespace std;

namespace BaskerNS
{


  //Note: we will want to come back and make
  //a much better multivector solve interface
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::solve_interface
  (
   Int _nrhs,
   Entry *_x, // Solution
   Entry *_y  // rhs
  )
  {
    for(Int r = 0; r < _nrhs; r++)
    {
      solve_interface(&(_x[r*gm]), &(_y[r*gm]));
    }

    return 0;
  }//end solve_interface(_nrhs,x,y);


  // _x will be solution (properly permuted)
  // _y is originally the rhs
  // In this function, rhs _y is copied and permuted to x_view_ptr_copy
  // In subsequent solver calls, x_view_ptr_copy (initially permuted rhs)
  // is updated/modified during block solves
  // y_view_ptr_copy stores the solution pivots
  // After solve is complete, the y_view_ptr_copy results are permuted
  // and copied to the raw pointer _x
  //
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::solve_interface
  (
   Entry *_x, // Solution (len = gn)
   Entry *_y  // rhs
  )
  {
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf( "\n -- solve_interface --\n" );
    //for (Int i = 0; i < gn; i++) printf( " input: x(%d) = %e\n",i,_y[i] );
    //printf( "\n" );
    #endif
    if (Options.blk_matching != 0 || Options.static_delayed_pivot != 0) {
      // apply mwm+amd row scaling from numeric
      for(Int i = 0; i < gn; i++) {
        Int row = order_blk_mwm_array(symbolic_row_iperm_array(i));
        y_view_ptr_scale(i) = scale_row_array(row) * _y[i];
        //printf( " symbolic_row_iperm(%d) = %d\n",i,symbolic_row_iperm_array(i) );
        //printf( " scale_row(%d) = %e\n",row,scale_row_array(row) );
      }
      //printf( " > after scale:\n" );
      //for (Int i = 0; i < gn; i++) printf( " > y(%d) = %.16e\n",i,y_view_ptr_scale(i) );

      // apply mwm row-perm from symbolic
      //for (Int i = 0; i < gn; i++) printf( " > iperm(%d) = %d\n",i,perm_inv_comp_array(i) );
      permute_inv_and_init_for_solve(&(y_view_ptr_scale(0)), x_view_ptr_copy, y_view_ptr_copy, perm_inv_comp_array, gn);
      //printf( " > after symbolic-perm:\n" );
      //for (Int i = 0; i < gn; i++) printf( " %d %.16e %.16e\n",i, x_view_ptr_copy(i),y_view_ptr_copy(i) );

      // apply row-perm from setup at numeric phase
      //for (Int i = 0; i < gn; i++) printf( " > iperm(%d) = %d\n",i,numeric_row_iperm_array(i) );
      permute_with_workspace(x_view_ptr_copy, numeric_row_iperm_array, gn);
    } else {
      //for (Int i = 0; i < gn; i++) printf( " > iperm(%d) = %d\n",i,perm_inv_comp_array(i) );
      // apply matrix_ordering from symbolic
      permute_inv_and_init_for_solve(_y, x_view_ptr_copy, y_view_ptr_copy, perm_inv_comp_array, gn);
      if (Options.matrix_scaling != 0) {
        for(Int i = 0; i < gn; i++) {
          x_view_ptr_copy(i) = x_view_ptr_copy(i) * scale_row_array(i);
        }
      }
    }
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf( " > after perm:\n" );
    for (Int i = 0; i < gn; i++) printf( " > perm_inv(%d) = %d\n",i, perm_inv_comp_array(i) );
    for (Int i = 0; i < gn; i++) printf( " %d %.16e %.16e\n",i, x_view_ptr_copy(i),y_view_ptr_copy(i) );
    printf( "\n" );
    #endif

    if (Options.no_pivot == BASKER_FALSE) {
      // apply partial pivoting from numeric
      //for (Int i = 0; i < gn; i++) printf( " gperm(%d) = %d\n",i,gperm(i) );
      permute_inv_with_workspace(x_view_ptr_copy, gperm, gn);
      //printf( " > after partial-pivot:\n" );
      //for (Int i = 0; i < gn; i++) printf( " %d %.16e %.16e\n",i, x_view_ptr_copy(i),y_view_ptr_copy(i) );
      //printf( "\n" );
    }

    // solve
    #ifdef BASKER_DEBUG_SOLVE_RHS
    //printf( "\n before solver-interface\n" );
    printf( " c = [\n" );
    for (Int i = 0; i < gn; i++) printf( " %d %.16e\n",i,x_view_ptr_copy(i) );
    printf( "];\n" );
    #endif
    solve_interface(x_view_ptr_copy, y_view_ptr_copy); //x is now permuted rhs; y is 0 
    #ifdef BASKER_DEBUG_SOLVE_RHS
    //printf( "\n after solver-interface\n" );
    printf( " y = [\n" );
    for (Int i = 0; i < gn; i++) printf( " %d %.16e %.16e\n",i,x_view_ptr_copy(i),y_view_ptr_copy(i) );
    printf( "];\n\n" );
    #endif

    if (Options.blk_matching != 0 || Options.static_delayed_pivot != 0) {
      // apply amd col-permutation from numeric
      permute_and_finalcopy_after_solve(&(y_view_ptr_scale(0)), x_view_ptr_copy, y_view_ptr_copy, numeric_col_iperm_array, gn);
      //for (Int i = 0; i < gn; i++) printf( " > %d:%d: %.16e %.16e -> %.16e\n",i,numeric_col_iperm_array(i),x_view_ptr_copy(i),y_view_ptr_copy(i), y_view_ptr_scale(i));

      /*const Int scol_top = btf_tabs[btf_top_tabs_offset]; // the first column index of A
      for (Int i = 0; i < scol_top; i++) {
        x_view_ptr_copy(i) = y_view_ptr_scale(i);
      }*/
      const Int poffset = btf_tabs(btf_tabs_offset);
      for (Int i = 0; i < poffset; i++) {
        x_view_ptr_copy(i) = y_view_ptr_scale(i);
      }
      for (Int i = poffset; i < gn; i++) {
        y_view_ptr_copy(i) = y_view_ptr_scale(i);
      }
    }
    if (Options.matrix_scaling != 0) {
      const Int poffset = btf_tabs(btf_tabs_offset);
      for (Int i = 0; i < poffset; i++) {
        x_view_ptr_copy(i) = x_view_ptr_copy(i) * scale_col_array(i);
      }
      for (Int i = poffset; i < gn; i++) {
        y_view_ptr_copy(i) = y_view_ptr_copy(i) * scale_col_array(i);
      }
    }

    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf( " before calling permute_and_finalcopy_after_solve\n" );
    for (Int i = 0; i < gn; i++) printf( " %d:%d; %e %e\n",i,perm_comp_array(i),x_view_ptr_copy(i),y_view_ptr_copy(i));
    #endif
    permute_and_finalcopy_after_solve(_x, x_view_ptr_copy, y_view_ptr_copy, perm_comp_array, gn);
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf( " after calling permute_and_finalcopy_after_solve\n" );
    for (Int i = 0; i < gn; i++) printf( " %d:%d; %e %e\n",i,symbolic_col_iperm_array(i), _x[i],scale_col_array(i));
    #endif

    if (Options.blk_matching != 0) {
      for(Int i = 0; i < gn; i++) {
        Int col = symbolic_col_iperm_array(i);
        _x[i] = scale_col_array(col) * _x[i];
      }
    }
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf( " final\n" );
    for (Int i = 0; i < gn; i++) printf( " %d %e\n",i,_x[i]);
    #endif
    return 0;
  }


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::solve_interface
  ( 
   ENTRY_1DARRAY & x, // x is permuted rhs at input
   ENTRY_1DARRAY & y  // y is 0 at input 
  )
  {
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("\n\n");
    printf("X: \n");
    for(Int i = 0; i < gn; i++)
    {
      printf("%f, " , x(i));
    }
    printf("\n\n");
    printf("RHS: \n");
    for(Int i =0; i < gm; i++)
    {
      printf("%f, ", y(i)); 
    }
    printf("\n\n");
    #endif

    if(Options.btf == BASKER_FALSE)
    {
      if(btf_tabs_offset != 0)
      {
        serial_solve(x,y);
      }
    }
    else
    {
      //A\y -> y
      serial_btf_solve(x,y);
    }

    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("\n\n");
    printf("X: \n");
    for(Int i = 0; i < gn; i++)
    {
      printf("%f, " , x(i));
    }
    printf("\n\n");
    printf("RHS: \n");
    for(Int i = 0; i < gm; i++)
    {
      printf("%f, ", y(i)); 
    }
    printf("\n\n");
    #endif

    return 0;
  }
  

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::serial_solve
  (
   ENTRY_1DARRAY & x, // Permuted rhs at input
   ENTRY_1DARRAY & y  // 0 at input
  )
  {
    //L\x -> y
    serial_forward_solve(x,y);

    //printVec(y,gn);

    for(Int i =0; i<gn; ++i)
    {
      x(i) = 0;
    }
    //U\y -> x
    serial_backward_solve(y,x);

    return 0;
  }//end serial solve()


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::serial_btf_solve
  (
   ENTRY_1DARRAY & x, // Permuted rhs at input
   ENTRY_1DARRAY & y  // 0 at input
  )
  {
    //Start in C and go backwards
    //In first level, only do U\L\x->y
    /*printf(" C = [\n" );
    for(Int j = 0; j < BTF_C.ncol; j++) {
      for(Int k = BTF_C.col_ptr[j]; k < BTF_C.col_ptr[j+1]; k++) {
        printf("%d %d %.16e\n", BTF_C.row_idx[k], j, BTF_C.val[k]);
      }
    }
    printf("];\n");*/
    Int nblks_c = btf_nblks-btf_tabs_offset;
    for(Int b = nblks_c-1; b>= 0; b--) {

      //---Lower solve
      BASKER_MATRIX &LC = LBTF(b);
      #ifdef BASKER_DEBUG_SOLVE_RHS
      printf("\n\n btf b=%ld (%d x %d), LBTF(%d)\n", (long)b, (int)LC.nrow, (int)LC.ncol, (int)b);
      #endif

      //L(C)\x -> y (x = y at output, with unit diagonal L)
      lower_tri_solve(LC,x,y);

      //printVec(y,gn);

      BASKER_MATRIX &UC = UBTF(b);
      //U(C)\x -> y
      upper_tri_solve(UC,x,y);

      #ifdef BASKER_DEBUG_SOLVE_RHS
      printf("Before spmv\n"); fflush(stdout);
      printf("Inner Vector y print\n");
      printVec(y, gn);
      printf("Inner Vector x print\n");
      printVec(x, gn);
      printf("\n");
      #endif

      //-----Update
      //if(b > btf_tabs_offset)
      {
        //x = BTF_C*y; 
        //  x(b) gets updated, though not needed (y(b) contains the solution)
        spmv_BTF(b+btf_tabs_offset, BTF_C, x, y);
      }

      #ifdef BASKER_DEBUG_SOLVE_RHS
      printf("After spmv\n"); fflush(stdout);
      printf("Inner Vector y print\n");
      printVec(y, gn);
      printf("Inner Vector x print\n");
      printVec(x, gn);
      #endif
    }
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf( " t1 = [\n" );
    for (Int i = 0; i < gn; i++) printf( " %d %.16e %.16e\n",i,x(i),y(i) );
    printf( "];\n\n" );
    #endif

    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("Done, BTF-C Solve \n"); fflush(stdout);
    for (Int i = 0; i < gn; i++) printf( " %e %e\n",x(i),y(i));
    printf( "\n");
    #endif

    //Update B
    //BTF_B*y -> x
    if(btf_tabs_offset !=  0)
    {
      neg_spmv_perm(BTF_B,y,x);
    }
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf( " t2 = [\n" );
    for (Int i = 0; i < gn; i++) printf( " %d %.16e %.16e\n",i,x(i),y(i) );
    printf( "];\n\n" );
    #endif

    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("Done, SPMV BTF_B UPDATE \n"); fflush(stdout);
    for (Int i = 0; i < gn; i++) printf( " %e %e\n",x(i),y(i));
    printf( "\n");
    #endif

    //now do the forward backward solve
    //L\x ->y
    serial_forward_solve(x,y);
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("Done, serial_forward \n"); fflush(stdout);
    for (Int i = 0; i < gn; i++) printf( "%d  %e %e\n",i,x(i),y(i));
    printf( "\n");
    #endif
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf( " t3 = [\n" );
    for (Int i = 0; i < gn; i++) printf( " %d %.16e %.16e\n",i,x(i),y(i) );
    printf( "];\n\n" );
    #endif

    //U\y->x
    serial_backward_solve(y,x);
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("Done, serial_backward \n"); fflush(stdout);
    for (Int i = 0; i < gn; i++) printf( " %d %e %e\n",i,x(i),y(i));
    printf( "\n");
    #endif
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf( " t4 = [\n" );
    for (Int i = 0; i < gn; i++) printf( " %d %.16e %.16e\n",i,x(i),y(i) );
    printf( "];\n\n" );
    #endif

    if(btf_top_tabs_offset >  0)
    {
      // copy the solution for C from y to x
      Int srow_c = BTF_C.srow;
      Int erow_c = BTF_C.nrow + srow_c;
      for (Int i = srow_c; i < erow_c; i++) {
        x(i) = y(i);
      }

      #ifdef BASKER_DEBUG_SOLVE_RHS
      printf( " \n D solve..\n" );
      for (Int i = 0; i < gn; i++) printf( " %d %e\n",i,x(i));
      printf( "\n");
      printf( " E = [\n" );
      for(Int j = 0; j < BTF_E.ncol; ++j) {
        for(Int k = BTF_E.col_ptr(j); k < BTF_E.col_ptr(j+1); k++) {
          printf( "%d %d %e\n",BTF_E.row_idx(k),j,BTF_E.val(k) );
        }
      }
      printf("];\n");
      #endif
      neg_spmv_perm(BTF_E, x, x);
      #ifdef BASKER_DEBUG_SOLVE_RHS
      printf( "\n after update\n" );
      for (Int i = 0; i < gn; i++) printf( " %e\n",x(i));
      printf( "\n");
      #endif

      for(Int b = btf_top_tabs_offset-1; b>= 0; b--)
      {
        //L(C)\x -> y 
        BASKER_MATRIX &LC = L_D(b);
        lower_tri_solve(LC, x, y);
        #ifdef BASKER_DEBUG_SOLVE_RHS
        printf( "\n after L solve (b=%d)\n",b ); fflush(stdout);
        for (Int i = 0; i < gn; i++) printf( " %e %e\n",x(i),y(i));
        printf( "\n");
        #endif

        //U(C)\y -> x
        BASKER_MATRIX &UC = U_D(b);
        upper_tri_solve(UC, y, x);
        #ifdef BASKER_DEBUG_SOLVE_RHS
        printf( "\n after U solve\n" ); fflush(stdout);
        for (Int i = 0; i < gn; i++) printf( " %e %e\n",x(i),y(i));
        printf( "\n");
        #endif
        {
          //x = BTF_C*y;
          spmv_BTF(b, BTF_D, x, x, false);
        }
        #ifdef BASKER_DEBUG_SOLVE_RHS
        printf( "\n after spmv\n" ); fflush(stdout);
        for (Int i = 0; i < gn; i++) printf( " %e %e\n",x(i),y(i));
        printf( "\n");
        #endif
      }
      #ifdef BASKER_DEBUG_SOLVE_RHS
      printf("Done, D solve\n");
      printf("\n x \n");
      printVec(x, gn);
      printf("\n y \n");
      printVec(y, gn);
      printf("\n\n");
      #endif
    }

    return 0;
  }//end serial_btf_solve


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::serial_forward_solve
  (
   ENTRY_1DARRAY & x, // modified rhs
   ENTRY_1DARRAY & y  // partial solution
  )
  {
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("Called serial forward solve \n");
    #endif

    Int scol_top = btf_tabs[btf_top_tabs_offset]; // the first column index of A
    //Forward solve on A
    for(Int b = 0; b < tree.nblks; ++b)
    {
      BASKER_MATRIX &L = LL(b)(0);

      //L\x -> y 
      lower_tri_solve(L, x, y, scol_top);
      #ifdef BASKER_DEBUG_SOLVE_RHS
      printf("Lower Solve blk (%d, 0): size=(%dx%d) srow=%d, scol=%d\n",b,(int)L.nrow,(int)L.ncol, (int)L.srow,(int)L.scol);
      printf("[\n");
      for (Int i = 0; i < gn; i++) printf( "%d %.16e %.16e\n",i,x(i),y(i));
      printf( "];\n");
      if (b == 3) {
        const Int bcol = L.scol + scol_top;
        const Int brow = L.scol + scol_top;
        printf( " P = [\n" );
        for (Int k = 0; k < L.ncol; k++) printf( "%d %d\n",brow+k,gperm(brow+k) );
        printf("];\n");
        char filename[200];
        sprintf(filename,"LL_%d_%d.dat",b,0);
        L.print_matrix(filename);
      }
      #endif

      //Update offdiag
      for(Int bb = 1; bb < LL_size(b); ++bb)
      {
        BASKER_MATRIX &LD = LL(b)(bb);
        //x = LD*y;
        #ifdef BASKER_DEBUG_SOLVE_RHS
        char filename[200];
        sprintf(filename,"LL_%d_%d.dat",b,bb);
        LD.print_matrix(filename);
        //printf( "];\n" );
        //std::cout << " ++ neg_spmv_perm ( LL(" << b << "," << bb << ") )" << std::endl;
        #endif
        neg_spmv_perm(LD, y, x, scol_top);
        #ifdef BASKER_DEBUG_SOLVE_RHS
        printf("Lower Solver Update blk = (%d %d): size=(%dx%d) srow=%d, scol=%d \n",
               b, bb, (int)LD.nrow,(int)LD.ncol, (int)LD.srow,(int)LD.scol);
        printf("[\n");
        for (Int i = 0; i < gn; i++) printf( "%d  %e %e\n",i,x(i),y(i));
        printf( "];\n");
        #endif
      }
      //printVec(y,gn);
    }

    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("Done forward solve A \n");
    printVec(y, gn);
    #endif

    return 0;
  }//end serial_forward_solve()

  template<class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::serial_backward_solve
  (
   ENTRY_1DARRAY & y,
   ENTRY_1DARRAY & x
  )
  {
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("called serial backward solve \n");
    #endif

    Int scol_top = btf_tabs[btf_top_tabs_offset]; // the first column index of A
    for(Int b = tree.nblks-1; b >=0; b--)
    {
      #ifdef BASKER_DEBUG_SOLVE_RHS
      printf("Upper solve blk: %d, %d \n", b,LU_size(b)-1);
      #endif

      //U\y -> x
      BASKER_MATRIX &U = LU(b)(LU_size(b)-1);
      upper_tri_solve(U, y, x, scol_top); // NDE: y , x positions swapped...
                                          //      seems role of x and y changed...
      #ifdef BASKER_DEBUG_SOLVE_RHS
      {
        char filename[200];
        sprintf(filename,"LU_%d_%d.dat",b,LU_size(b)-1);
        U.print_matrix(filename);
      }
      #endif

      for(Int bb = LU_size(b)-2; bb >= 0; bb--)
      {
        #ifdef BASKER_DEBUG_SOLVE_RHS
        printf("Upper solver spmv: %d %d \n",
            b, bb);
        #endif

        //y = UB*x;
        BASKER_MATRIX &UB = LU(b)(bb);
        neg_spmv(UB, x, y, scol_top);

        #ifdef BASKER_DEBUG_SOLVE_RHS
        {
          char filename[200];
          sprintf(filename,"LU_%d_%d.dat",b,bb);
          UB.print_matrix(filename);
        }
        #endif
      }
    }//end over all blks

    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("Done with Upper Solve: \n");
    printVec(x, gn);
    #endif

    return 0;
  }//end serial_backward_solve()


  //Horrible, cheap spmv
  //y = M*x
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::spmv
  (
    BASKER_MATRIX &M,
    ENTRY_1DARRAY x,
    ENTRY_1DARRAY y
  )
  {
    //Add checks
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("SPMV. scol: %d ncol: %d nnz: %d \n",
        M.scol, M.ncol, M.nnz);
    M.info();
    #endif

    const Int bcol = M.scol;
    const Int brow = M.srow;
    //for(Int k=M.scol; k < (M.scol+M.ncol); k++)
    for(Int k = 0; k < M.ncol; ++k)
    {
      const auto xkbcol = x(k+bcol);
      const Int istart = M.col_ptr(k);
      const Int iend = M.col_ptr(k+1);

      //for(Int i = M.col_ptr(k); i<M.col_ptr(k+1); ++i)
      for(Int i = istart; i<iend; ++i)
      {
        const Int j = M.row_idx(i);

        //y(j+brow) += M.val(i)*x(k+bcol);
        y(j+brow) += M.val(i)*xkbcol;

      }
    }
    return 0;
  }//spmv


  //Horrible, cheap spmv
  //y = M*x
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::neg_spmv
  (
   BASKER_MATRIX &M,
   ENTRY_1DARRAY x, 
   ENTRY_1DARRAY y,
   Int offset
  )
  {
    //Add checks
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("SPMV. scol: %d ncol: %d \n", M.scol, M.ncol);
    #endif
    //M.print_matrix("M.dat");

    const Int bcol  = M.scol + offset;
    const Int msrow = M.srow + offset;
    //const Int brow = M.srow;
    for(Int k=0; k < M.ncol; ++k)
    {
      const auto xkbcol = x(k+bcol);
      const Int istart = M.col_ptr(k);
      const Int iend = M.col_ptr(k+1);
      for(Int i = istart; i < iend; ++i)
      {
        const Int j = M.row_idx(i) + msrow;

        //y(j) -= M.val(i)*x(k+bcol);
        y(j) -= M.val(i)*xkbcol;
      }
    }

    return 0;
  }//neg_spmv


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::neg_spmv_perm
  (
   BASKER_MATRIX &M,
   ENTRY_1DARRAY &y, 
   ENTRY_1DARRAY &x,
   Int offset
  )
  {
    //Add checks
    const Int bcol  = M.scol + offset;
    const Int msrow = M.srow + offset;
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("SPMV. offset = %d, srow = %d, scol = %d: nrow = %d x ncol = %d\n", offset, M.srow, M.scol, M.nrow, M.ncol);

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

    //std::cout << " + neg_spmv_perm ( " << M.nrow << " x " << M.ncol << " )" << std::endl;
    for(Int k=0; k < M.ncol; ++k)
    {
      const Int istart = M.col_ptr(k);
      const Int iend   = M.col_ptr(k+1);
      const auto ykbcol = y(k+bcol);

      for(Int i = istart; i < iend; ++i) //NDE retest with const vars, scope tightly
      {
        const Int j = (Options.no_pivot == BASKER_FALSE) ? 
                       gperm(M.row_idx(i) + msrow) :
                       (M.row_idx(i) + msrow) ;

        //if (M.nrow == 289 && M.ncol == 100) printf( " x(%d) = %e - %e * %e (%d, %d->%d->%d)\n",j,x(j),M.val(i),ykbcol,k, M.row_idx(i),M.row_idx(i) + msrow,gperm(M.row_idx(i) + msrow) );
        x(j) -= M.val(i)*ykbcol;
      }
    }

    return 0;
  }//neg_spmv


  //M\x = y
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::lower_tri_solve
  (
   BASKER_MATRIX &M,
   ENTRY_1DARRAY &x, 
   ENTRY_1DARRAY &y,
   Int offset
  )
  {
    const Int bcol = M.scol + offset;
    const Int brow = M.scol + offset;

    /*printf( " P = [\n" );
    for (Int k = 0; k < M.ncol; k++) printf( "%d %d\n",brow+k,gperm(brow+k) );
    printf("];\n");
    M.print_matrix("L.dat"); */

    for(Int k = 0; k < M.ncol; ++k)
    {
      //Test if zero pivot value
      #ifdef BASKER_DEBUG_SOLVE_RHS
      BASKER_ASSERT(M.val[M.col_ptr[k]]!=0.0, "LOWER PIVOT 0");
      //printf("Lower tri.  k: %d out: %f in: %f piv: %f \n",
      //   k+bcol, y[k+bcol], x[k+bcol], M.val[M.col_ptr[k]]);
      #endif

      // TODO NDE: Need to make sure this is properly checked in numeric factorization
      /*
      if(M.val[M.col_ptr[k]] == 0.0) 
      {
        printf("Lower Pivot: %d %f \n", 
            M.row_idx[M.col_ptr[k]],
            M.val[M.col_ptr[k]]);
        return -1;
      }
      */

      //Replace with Entry divide in future
      const Int istart = M.col_ptr(k);
      const Int iend = M.col_ptr(k+1);

      //printf( " %d %d %e\n",M.row_idx(M.col_ptr(k)),k,M.val(M.col_ptr(k)));
      //printf( " -> %e %e (%d, %d,%d)\n",y(k+brow),x(k+bcol),k,brow,bcol );
      y(k+brow) = x(k+bcol) / M.val(M.col_ptr(k));
      //printf( " y(%d) = %e / %e = %e\n",k+brow,x(k+bcol), M.val(M.col_ptr(k)), y(k+bcol));

      const auto ykbcol = y(k+bcol);
      //for(Int i = M.col_ptr(k)+1; i < M.col_ptr(k+1); ++i)
      for(Int i = istart+1; i < iend; ++i)
      {
        const Int j = (Options.no_pivot == BASKER_FALSE) ? 
                        gperm(M.row_idx(i)+brow) :
                             (M.row_idx(i)+brow) ;

        #ifdef BASKER_DEBUG_SOLVE_RHS
        BASKER_ASSERT(j != BASKER_MAX_IDX,"Using nonperm\n");
        #endif

        //x(j) -= M.val(i)*y(k+bcol);
        x(j) -= M.val(i)*ykbcol;
      } //over all nnz in a column

    } //over each column

    return 0;
  } //end lower_tri_solve


  //U\x = y
  // Note: In serial_backward_solve usage, the vars do not match up
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::upper_tri_solve
  (
   BASKER_MATRIX &M,
   ENTRY_1DARRAY &x,
   ENTRY_1DARRAY &y,
   Int offset
  )
  {
    const Int bcol = M.scol + offset;
    const Int brow = M.srow + offset;
    // M.print_matrix("U.dat);

    for(Int k = M.ncol; k >= 1; k--)
    {
      //printf( " -- k = %d --\n",k );
      #ifdef BASKER_DEBUG_SOLVE_RHS
      BASKER_ASSERT(M.val[M.col_ptr[k]-1]!=0.0,"UpperPivot\n");
      printf("Upper Tri Solve, scol: %d ncol: %d \n",
        M.scol, M.ncol);

      #endif

      // TODO NDE: Need to make sure this is properly checked in numeric factorization
      /*
      if(M.val(M.col_ptr(k)-1)==0)
      {
        printf("Upper pivot: %d %f \n",
            M.row_idx[M.col_ptr[k]-1],
            M.val[M.col_ptr[k]-1]);
        return -1;
      }
      */

      //Comeback and do with and entry divide
      const Int istart = M.col_ptr(k);
      const Int iend = M.col_ptr(k-1);

      //printf( " > y(%d + %d -1 = %d) = x(%d + %d -1 = %d) -> %e / %e \n",k,brow,k-brow-1, k,bcol, k+bcol-1, x(k+bcol-1), M.val(M.col_ptr(k)-1) );
      y(k+brow-1)  =  x(k+bcol-1) / M.val(M.col_ptr(k)-1);

      const auto ykbcol = y(k+bcol-1);
      //for(Int i = M.col_ptr(k)-2; i >= M.col_ptr(k-1); --i) 
      for(Int i = istart-2; i >= iend; --i)
      {
        const Int j = M.row_idx(i) + brow; //NDE: why isn't gperm here like above?

        //x(j) -= M.val(i) * y(k+bcol-1);
        //printf( " > x (%d) -= %e %e\n",j,M.val(i),ykbcol );
        x(j) -= M.val(i) * ykbcol;
      }

    }//end over all columns

    return 0;
  } //end upper_tri_solve


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::spmv_BTF
  (
   Int tab,
   BASKER_MATRIX &M,
   ENTRY_1DARRAY &x, // modified rhs
   ENTRY_1DARRAY &y, // intermediate solution
   bool full
  )
  {
    //Tab = block in    
    Int bcol = btf_tabs(tab)- M.scol;
    Int mscol = M.scol;
    Int brow = M.srow;
    Int ecol = btf_tabs(tab+1) - M.scol;
    if (ecol > M.ncol) {
      // for D block, btf_tabs(tab+1) > ncol.
      ecol = M.ncol;
    }

    #ifdef BASKER_DEBUG_SOLVE_RHS
    Int erow = 0;
    if(tab > 0)
    {
      erow = btf_tabs(tab);
    }
    else
    {
      erow = brow-1;
    }
    printf("BTF_UPDATE, TAB: %d [%d %d] [%d %d] \n",
        tab, brow, erow, bcol, ecol);
    #endif

    //loop over each column
    //printf( "spmv_BTF(bcol = %d, ecol = %d) tab = %d, scol = %d\n",bcol,ecol,tab,M.scol );
    for(Int k = bcol; k < ecol; ++k) {
      //const Int kcol = k+M.scol;
      const Int istart = M.col_ptr(k);
      const Int iend = M.col_ptr(k+1);

      //printf( " M.col_ptr(%d) = %d, M.col_ptr(%d) = %d\n",k,istart,k+1,iend );
      const auto ykmcol = y(k+mscol);
      //for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); ++i) 
      for(Int i = istart; i < iend; ++i) {
        const Int j = (Options.no_pivot == BASKER_FALSE) ? 
                        gperm(M.row_idx(i)+brow) :
                        (M.row_idx(i)+brow) ;

       #ifdef BASKER_DEBUG_SOLVE_RHS
        printf("BTF_UPDATE-val, j: %d x: %f y: %f, val: %f \n",
            j, x[j], y[k+M.scol], M.val[i]);
       #endif

        if (full || j < bcol) {
          //printf( " x(%d) -= M(%d,%d) * x(%d)\n",j,j,i,k+mscol );
          //x(j) -= M.val(i)*y(k+M.scol);
          x(j) -= M.val(i)*ykmcol;
        }
      } //over all nnz in row
    } // end for over col

    return 0;
  } //end spmv_BTF();
  
} //end namespace BaskerNS
#endif //end ifndef basker_solver_rhs
