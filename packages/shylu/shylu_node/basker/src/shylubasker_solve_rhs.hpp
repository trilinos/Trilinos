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
#include <impl/Kokkos_Timer.hpp>
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
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::test_solve()
  {
    ENTRY_1DARRAY  x_known;
    ENTRY_1DARRAY  x;
    ENTRY_1DARRAY  y;


    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("test_solve called \n");
    printf("Global pivot permuation\n");
    printVec(gperm, gn);
    printf("\n");
    printf("Global pivot permutation inverse\n");
    printVec(gpermi, gn);
    printf("\n");
    #endif


    BASKER_ASSERT(gn > 0, "solve testsolve gn");
    MALLOC_ENTRY_1DARRAY(x_known, gn);
    init_value(x_known, gn , (Entry)1.0);


    //temp
    for(Int i = 0; i < gn; i++)
    {
      x_known(i) = (Entry) 1.0;
    }
    //JDB: used for other test
    //permute(x_known, order_csym_array, gn);

    MALLOC_ENTRY_1DARRAY(x, gn);
    init_value(x, gn, (Entry) 0.0);
    BASKER_ASSERT(gm > 0, "solve testsolve gm");
    MALLOC_ENTRY_1DARRAY(y, gm);
    init_value(y, gm, (Entry) 0.0);

    if(btf_nblks > 0)
    {
      sort_matrix(BTF_C);
      //printMTX("C_BEFORE_SOLVE.mtx", BTF_C);
    }

    if(Options.btf == BASKER_TRUE)
    {
      if(btf_tabs_offset != 0)
      {
        spmv(BTF_A, x_known,y);
        if(btf_nblks> 1)
        {
          spmv(BTF_B, x_known, y);
        }
      }
      if(btf_nblks > 1)
      {
        spmv(BTF_C, x_known, y);
      }
    }
    else
    {
      //printf("other\n");
      //spmv(BTF_A, x_known,y);
    }

    for(Int i = 0; i < gn; i++)
    {
      if(gperm(i) < 0)
      {
        printf("Basker test_solve error: %ld %ld \n",
            (long)i, (long)gperm(i));
      }
      if(gperm(i) > gn)
      {
        printf("Basker test_solve serror: %ld %ld \n",
            (long)i, (long)gperm(i));
      }
      x(gperm(i)) = y(i);
    }
    for(Int i = 0; i < gn; i++)
    {
      y(i) = x(i);
      x(i) = 0;
    }

    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("\n\n");
    printf("RHS: \n");
    for(Int i =0; i < gm; i++)
    {
      printf("%ld %e,\n ", (long)i, y(i)); 
    }
    printf("\n\n");
    #endif

    if(Options.btf == BASKER_FALSE)
    {

      if(btf_tabs_offset != 0)
      {
        serial_solve(y,x);
      }
    }
    else
    {
      //A\y -> y
      //serial_btf_solve(y,x);
      //printf("before btf serial solve\n");
      serial_btf_solve(y,x);
    }


    Entry diff =0.0;

    for(Int i = 0; i < gn; i++)
    {
      diff += (x_known(i) - x(i));
    }
    diff = diff/(Entry) gn;

    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("\n\n");
    printf("Solve Compare: \n");
    for(Int i = 0; i < gn; i++)
    {
      printf("%ld %f %f \n", 
          (long)i, x_known(i), x(i));

    }
    printf("\n\n");
    #endif

    printf("\n");
    printf("TEST_SOLVE: ||x-x||/||x| = %e", diff);
    printf("\n");

    if((diff > -1e-2) && (diff < 1e-2))
    {
      printf("TEST PASSED \n");
    }  

    return 0;
  }//end test_solve


  //Note: we will want to come back and make
  //a much better multivector solve interface
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::solve_interface
  (
   Int nrhs,
   Entry *_x, // Solution
   Entry *_y  // rhs
  )
  {

    for(Int r = 0; r < nrhs; r++)
    {
      solve_interface(&(_x[r*gm]), &(_y[r*gm]));
    }

    return 0;
  }//end solve_interface(nrhs,x,y);


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
    permute_inv_and_init_for_solve(_y, x_view_ptr_copy, y_view_ptr_copy, perm_inv_comp_array , gn);
    if (Options.no_pivot == BASKER_FALSE) {
      permute_inv_with_workspace(x_view_ptr_copy, gperm, gn);
    }

    solve_interface(x_view_ptr_copy,y_view_ptr_copy); //x is now permuted rhs; y is 0 

    permute_and_finalcopy_after_solve(_x, x_view_ptr_copy, y_view_ptr_copy, perm_comp_array, gn);

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
    for(Int b = (btf_nblks-btf_tabs_offset)-1; b>= 0; b--)
    {

    #ifdef BASKER_DEBUG_SOLVE_RHS
      printf("\n\n btf b: %ld \n", (long)b);
    #endif

      //---Lower solve
      BASKER_MATRIX &LC = LBTF(b);
      //L\x -> y 
      lower_tri_solve(LC,x,y);

      //printVec(y,gn);

      BASKER_MATRIX &UC = UBTF(b);
      //U\x -> y
      upper_tri_solve(UC,x,y);

      #ifdef BASKER_DEBUG_SOLVE_RHS
      printf("Before spmv\n");
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
        spmv_BTF(b+btf_tabs_offset, BTF_C, x, y);
      }

      #ifdef BASKER_DEBUG_SOLVE_RHS
      printf("After spmv\n");
      printf("Inner Vector y print\n");
      printVec(y, gn);
      printf("Inner Vector x print\n");
      printVec(x, gn);
      #endif

      //BASKER_MATRIX &UC = UBTF[b];
      //U\x -> y
      //upper_tri_solve(UC,x,y);
    }

    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("Done, BTF-C Solve \n");
    printf("\n x \n");
    printVec(x, gn);
    printf("\n y \n");
    printVec(y, gn);
    printf("\n\n");
    #endif

    //Update B
    //BTF_B*y -> x
    if(btf_tabs_offset !=  0)
    {
      neg_spmv_perm(BTF_B,y,x);
    }

    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("Done, SPMV BTF_B UPDATE \n");
    printf("\n x \n");
    printVec(x, gn);
    printf("\n y \n");
    printVec(y, gn);
    printf("\n\n");
    #endif

    //now do the forward backward solve
    //L\x ->y
    serial_forward_solve(x,y);
    //U\y->x
    serial_backward_solve(y,x);

    //copy lower part down
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("copying lower starting: %ld \n",
        (long)btf_tabs[btf_tabs_offset]);
    #endif

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

    //Forward solve on A
    for(Int b = 0; b < tree.nblks; ++b)
    {
    #ifdef BASKER_DEBUG_SOLVE_RHS
      printf("Lower Solve blk: %d \n",b);
    #endif

      BASKER_MATRIX &L = LL(b)(0);
      //L\x -> y 
      lower_tri_solve(L, x, y);

      //Update offdiag
      for(Int bb = 1; bb < LL_size(b); ++bb)
      {
        #ifdef BASKER_DEBUG_SOLVE_RHS
        printf("Lower Solver Update blk: %d %d \n",
            b, bb);
        #endif

        BASKER_MATRIX &LD = LL(b)(bb);
        //x = LD*y;
        neg_spmv_perm(LD, y, x);
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

    for(Int b = tree.nblks-1; b >=0; b--)
    {
      #ifdef BASKER_DEBUG_SOLVE_RHS
      printf("Upper solve blk: %d \n", b);
      #endif

      BASKER_MATRIX &U = LU(b)(LU_size(b)-1);

      //U\y -> x
      upper_tri_solve(U,y,x); // NDE: y , x positions swapped...
                              //      seems role of x and y changed...

      for(Int bb = LU_size(b)-2; bb >= 0; bb--)
      {
        #ifdef BASKER_DEBUG_SOLVE_RHS
        printf("Upper solver spmv: %d %d \n",
            b, bb);
        #endif

        BASKER_MATRIX &UB = LU(b)(bb);
        //y = UB*x;
        neg_spmv(UB,x,y);
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
   ENTRY_1DARRAY y  
  )
  {
    //Add checks
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("SPMV. scol: %d ncol: %d \n", M.scol, M.ncol);
    #endif

    const Int bcol = M.scol;
    const Int msrow = M.srow;
    //const Int brow = M.srow;
    for(Int k=0; k < M.ncol; ++k)
    {
      const auto xkbcol = x(k+bcol);
      const Int istart = M.col_ptr(k);
      const Int iend = M.col_ptr(k+1);
      //for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); ++i)
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
   ENTRY_1DARRAY &x  
  )
  {
    //Add checks
    #ifdef BASKER_DEBUG_SOLVE_RHS
    printf("SPMV. scol: %d ncol: %d \n", M.scol, M.ncol);
    #endif

    const Int bcol = M.scol;
    const Int msrow = M.srow;

    //for(Int k=M.scol; k < (M.scol+M.ncol); k++)
    for(Int k=0; k < M.ncol; ++k)
    {
      const Int istart = M.col_ptr(k);
      const Int iend   = M.col_ptr(k+1);
      const auto ykbcol = y(k+bcol);

      //for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); ++i) 
      for(Int i = istart; i < iend; ++i) //NDE retest with const vars, scope tightly
      {
        // const Int j = M.row_idx(i) + msrow;

        const Int j = (Options.no_pivot == BASKER_FALSE) ? 
                       gperm(M.row_idx(i) + msrow) :
                       (M.row_idx(i) + msrow) ;

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
   ENTRY_1DARRAY &y  
  )
  {
    const Int bcol = M.scol;
    const Int brow = M.scol;

    //M.info();

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

      y(k+brow) = x(k+bcol) / M.val(M.col_ptr(k));

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
   ENTRY_1DARRAY &y 
  )
  {
    const Int bcol = M.scol;
    const Int brow = M.srow;

    for(Int k = M.ncol; k >= 1; k--)
    {

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

      y(k+brow-1)  =  x(k+bcol-1) / M.val(M.col_ptr(k)-1);

      const auto ykbcol = y(k+bcol-1);
      //for(Int i = M.col_ptr(k)-2; i >= M.col_ptr(k-1); --i) 
      for(Int i = istart-2; i >= iend; --i)
      {
        const Int j = M.row_idx(i) + brow; //NDE: why isn't gperm here like above?

        //x(j) -= M.val(i) * y(k+bcol-1);
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
   ENTRY_1DARRAY &y  // intermediate solution
  )
  {
    //Tab = block in    
    const Int bcol = btf_tabs(tab)- M.scol;
    const Int mscol = M.scol;
    const Int brow = M.srow;
    const Int ecol = btf_tabs(tab+1) - M.scol;

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
    for(Int k = bcol; k < ecol; ++k)
    {
      //const Int kcol = k+M.scol;
      const Int istart = M.col_ptr(k);
      const Int iend = M.col_ptr(k+1);

      const auto ykmcol = y(k+mscol);
      //for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); ++i) 
      for(Int i = istart; i < iend; ++i)
      {
        const Int j = (Options.no_pivot == BASKER_FALSE) ? 
                        gperm(M.row_idx(i)+brow) :
                        (M.row_idx(i)+brow) ;

       #ifdef BASKER_DEBUG_SOLVE_RHS
        printf("BTF_UPDATE-val, j: %d x: %f y: %f, val: %f \n",
            j, x[j], y[k+M.scol], M.val[i]);
       #endif

        //x(j) -= M.val(i)*y(k+M.scol);
        x(j) -= M.val(i)*ykmcol;
      } //over all nnz in row

    } // end for over col

    return 0;
  } //end spmv_BTF();
  
} //end namespace BaskerNS
#endif //end ifndef basker_solver_rhs
