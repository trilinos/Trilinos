// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLUBASKER_NFACTOR_COL_HPP
#define SHYLUBASKER_NFACTOR_COL_HPP

//#include "shylubasker_decl.hpp"
#include "shylubasker_matrix_decl.hpp"
#include "shylubasker_matrix_view_decl.hpp"
#include "shylubasker_matrix_view_def.hpp"
#include "shylubasker_types.hpp"
#include "shylubasker_stats.hpp"
#include "shylubasker_thread.hpp"

#include "shylubasker_nfactor_blk.hpp"
#include "shylubasker_nfactor_blk_inc.hpp"
#include "shylubasker_nfactor_col_inc.hpp"

#include <Kokkos_Core.hpp>
#include <Kokkos_Timer.hpp>

#ifdef BASKER_DEBUG
#include <assert.h>
#endif

#ifdef HAVE_VTUNE
#include "ittnotify.h"
#endif


//#define MY_DEBUG_BASKER
//#define BASKER_DEBUG_NFACTOR_COL
//#define BASKER_DEBUG_TIME
//#define BASKER_COUNT_OPS

namespace BaskerNS
{
  //local idx for local blks
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_upper_col_factor
  (
   Int kid,
   Int team_leader,
   Int lvl,
   Int l,
   Int k, 
   BASKER_BOOL sep_flg
  )
  {
    const Entry zero (0.0);

    //Get needed variables
    const Int L_col = S(l)(kid);
    const Int U_col = S(lvl)(kid);

    Int my_row_leader = S(0)(find_leader(kid,lvl-1));
    //Int my_new_row = 
    // L_col - my_row_leader;
    Int U_row = L_col - my_row_leader;

    const Int X_col = S(0)(kid);
    const Int X_row = l; //X_row = lower(L)
    //const Int col_idx_offset = 0; //we might be able to remove
  
    #ifdef MY_DEBUG_BASKER
    const Int L_row = 0; //NDE - warning: unused 
    if(kid >= 0)
      printf("kid %d, upper using L(%d, %d)  U(%d, %d)  X(%d, %d)\n",
             (int)kid, (int)L_col, (int)L_row, (int)U_col, (int)U_row, (int)X_col, (int)X_row);
    #endif
    //end get needed variables//

    BASKER_MATRIX  &U = LU(U_col)(U_row);
    
    //Ask C++ guru if this is ok
    BASKER_MATRIX        *Bp;
    //Int                  bbcol = k;
    //if(sep_flg == BASKER_FALSE)
    if(l == 0)
    {
      Bp = &(AVM(U_col)(U_row));
      //bbcol = Bp->scol;
    }
    else
    {
      Bp = &(thread_array(kid).C);
      //printf("Using temp matrix, kid: %d\n", kid);
      //Bp->print();
    }
    BASKER_MATRIX    &B = *Bp;
    //B.print();
    //printf("POINT TEST, kid: %d X: %d %d \n",
    //           kid, X_col, X_row);


    INT_1DARRAY ws     = LL(X_col)(X_row).iws;
    const Int ws_size  = LL(X_col)(X_row).iws_size;
    ENTRY_1DARRAY X    = LL(X_col)(X_row).ews;

    const Int scol_top = btf_tabs[btf_top_tabs_offset]; // the first column index of A
    const Int brow_a = U.srow; // offset within A
    const Int brow_g = brow_a + scol_top; // global offset
    //const Int bcol = U.scol;

    Int *color     = &(ws(0));
    Int *pattern   = &(color[ws_size]);
    
    Int j, t, xnnz;
    Int top = ws_size;

    Int unnz = 0 ; 
    //if(k != U.scol)
    if(k!=0)
    { unnz = U.col_ptr(k); }
    //{unnz = U.col_ptr(k-U.scol);}

    Int uunnz = U.nnz;
   
    Int lcnt = 0;
    Int ucnt = 0; 
    
    #ifdef MY_DEBUG_BASKER
    if(kid == 0) {
      printf("kid: %d col_start Unnz: %d %d \n", kid, unnz, uunnz);
    }
    #endif


    #ifdef BASKER_DEBUG_NFACTOR_COL
    printf("---------Checking Workspace, kid: %d-------\n", kid);
    BASKER_ASSERT(top == ws_size);
    for(Int i = 0 ; i < ws_size; i++)
    {
      if(ws[i] != 0)
      {
        printf("--------------ERROR---------");
        BASKER_MATRIX  &L = LL(L_col)(L_row); //NDE - warning: unused L
        printf("kid: %d k: %d i: %d ws[i]=%d L.scol: %d \n",
               kid, k, i, ws[i], L.scol);
        BASKER_ASSERT(ws[i] == 0);
      }
    }
    #endif

    Int k_offset = 0;
    if(l != 0)
    {
      //printf("sep flag true, offset: %d \n",
      //     B.scol);
      k_offset = k;
    }

    int info = BASKER_SUCCESS;
    for(Int i = B.col_ptr(k-k_offset); i < B.col_ptr(k-k_offset+1); ++i)
    {
      j = B.row_idx(i);

      //if (B.val(i) != zero) 
      {
        X(j) = B.val(i);

        #ifdef BASKER_DEBUG_NFACTOR_COL
        //if(kid>=0)
        printf("t_upper_col_factor(kid=%d, l=%d, k=%d): X(%d) = %e\n", 
               kid,l,k, j, B.val(i));
        #endif

        #ifdef MY_DEBUG_BASKER
         #ifdef BASKER_2DL
         if(kid == 0)
           printf(" > kid=%d: X(%d) = %g color = %d \n",
                  (int)kid, (int)j, X[j], (int)ws[j] );
         #else
         if(kid == 0)
           printf(" > kid=%d: X(%d) = %g color = %d \n",
                  kid, j, X[j], ws[0 + j] );
         #endif
        #endif

        #ifdef BASKER_2DL
        //if(color[j-brow] == 0)
        if(color[j] == 0)
        #else
        if(color[j] == 0)
        #endif
        {
          #ifdef BASKER_INC_LVL
          //t_local_reach_selective(kid, l, l, j, &top);
          #else
          //t_local_reach(kid, l, l, j, &top); //Note: comeback
          t_local_reach(kid,l,l,j,top);
          #endif
        }//if not colored
      }// if not zero
    }//end over each nnz in column
    xnnz = ws_size - top;
    
    #ifdef BASKER_DEBUG_NFACTOR_COL
    if(kid == 0)
      printf("xnnz: %d ws_size: %d top: %d , kid: %d\n",
          xnnz, ws_size, top, kid);
    #endif

    //Count ops to show imbalance
    #ifdef BASKER_COUNT_OPS
    thread_array(kid).ops_counts[0][l] += xnnz;
    #endif

    //WE SHOUD DO A UNNZ COUNT
    //count number of nnz
    for(Int i=top; i < ws_size; ++i)
    {
      j = pattern[i];
      t = gperm(j+brow_g);
      if(t == BASKER_MAX_IDX)
      {
        lcnt++;
      }
    }
    //Note: This +1 causes some trouble
    ucnt = ws_size - top - lcnt +1;
     
    #ifdef BASKER_DEBUG_NFACTOR_COL
     if(kid>=0)
       printf("lcnt: %d ucnt: %d , kid: %d \n", (int)lcnt, (int)ucnt, (int)kid);
    #endif


    //#ifdef BASKER_DEBUG_NFACTOR_COL
    if(unnz+ucnt-1 > uunnz)
    {
      if (Options.verbose == BASKER_TRUE)
      {
        printf("kid: %ld col: %ld need to realloc, unnz: %ld ucnt: %ld uunnz: %ld U_col: %ld U_row: %ld \n",
               (long)kid, (long)k, (long)unnz, (long)ucnt, (long)uunnz, (long)U_col, (long)U_row);
      }
      //Note: commented out.. Does this work?
      //BASKER_ASSERT(0==1, "USIZE\n");

      Int newsize = (unnz+U.nrow) * 1.2  ;

      thread_array(kid).error_blk    = U_col;
      thread_array(kid).error_subblk = U_row;
      if(Options.realloc == BASKER_FALSE)
      {
        thread_array(kid).error_type = BASKER_ERROR_NOMALLOC;
        return BASKER_ERROR;
      }
      else
      {
        thread_array(kid).error_type = BASKER_ERROR_REMALLOC;
        thread_array(kid).error_info = newsize;
        return BASKER_ERROR;
      }//if/else realloc
    }
    //#endif


    // --------------
    // solve with L
    #ifdef BASKER_INC_LVL
    //Note we have not made a t_back_solve_atomic_selective
    t_back_solve_selective(kid, l,l, k, top, xnnz);
    #else //BASKER_INC_LVL
    #ifdef BASKER_ATOMIC_2
    t_back_solve(kid, l, l, k, top, xnnz);
    #else
    t_back_solve_atomic(kid, team_leader,
                 lvl, l,k , top,
                 xnnz);
    #endif
    #endif //BASKER_INC_LVL


    // --------------
    //move over nnz to U 
    for(Int i = top; i < ws_size && info == BASKER_SUCCESS; ++i)
    {
      j = pattern[i];
      t = gperm(j+brow_g);

      #ifdef BASKER_DEBUG_NFACTOR_COL
      if(kid == 0)
      {
        printf( " %d: gperm(%d+%d = %d) = %d\n",kid,j,brow_g,j+brow_g,t );
        printf( " %d: considering j: %d t:%d val: %e\n",kid,j, t, X[j]);
      }
      #endif

      #ifdef BASKER_2DL
      //if(X[j-brow] !=0)
      if (X(j) != zero)
      #else
      if (X[j] != zero)
      //kkos_nfactor_sep2
      #endif
      {
        //SHOULD:: REMOVE, check if this fine
        //if( (t != L.max_idx) && (t >= L.scol) && 
        //(t<(L.scol+L.ncol)))
        //if(t!=L.max_idx)

        //Note, if we remove this test, 
        //we might get this to unroll!!!
        if(t != BASKER_MAX_IDX)
        {
          #ifdef BASKER_DEBUG_NFACTOR_COL
          if(kid == 0)
            printf("kid: %d adding x[%d] to U(row=%d, val=%e)\n", kid, j, t, X(j)); 
          #endif

          //U.row_idx[unnz] = gperm[j];
          //printf("kid: %d addU: %d %d \n",
          //     kid, unnz, U.nnz);
          U.row_idx(unnz) = t-brow_g;
          #ifdef BASKER_2DL
          U.val(unnz) = X(j);
          #else
          U.val[unnz] = X[j];
          #endif
          unnz++;

          #ifdef BASKER_2DL
          //X[j-brow] = 0;
          X(j) = 0;
          #else
          X[j] = 0;
          #endif
        }//if in U
        else
        {
          #ifdef BASKER_2DL
          std::cout << "----Error--- kid = " << kid << ": extra L[" << j << "]="
                    << X[j] << " with gperm( " << brow_g << " + " << j << " ) = " << t
                    << std::endl;
          thread_array(kid).error_type = BASKER_ERROR_OTHER;
          thread_array(kid).error_blk    = lvl;
          thread_array(kid).error_subblk = l;
          thread_array(kid).error_info = k;
          info = BASKER_ERROR;
          //BASKER_ASSERT(t != BASKER_MAX_IDX, "lower entry in U");
          #endif
        }//lower
      }//end not 0
    }//over all x

    //U.col_ptr[k+1-bcol] = unnz;
    U.col_ptr(k+1) = unnz;

    //DEBUG
    //printf("TEST SIZE: %d %d \n", 
    //            U.col_ptr(k), U.col_ptr(k+1));

    #ifdef BASKER_2DL
  
    #ifdef BASKER_DEBUG_NFACTOR_COL
    printf("JB TEST kid: %d lvl: %d l: %d L_col: %d size: %d Xrow: %d\n",
           kid, lvl, l, L_col, LL_size[X_col], X_row);
    #endif

    #ifndef BASKER_MULTIPLE_UPPER
    //----------------------Update offdiag-----------------//
    X_row++;
    L_row++;
    //for(; X_row < LL_size[X_col]; X_row++, L_row++)
    for(; X_row < LL_size(X_col); X_row++, L_row++)
    {

      //printf("xrow: %d \n", X_row);
      const BASKER_BOOL A_option = BASKER_FALSE;
      //if((kid == team_leader) && blk_row == 1)
      //  {A_option = BASKER_TRUE;}

      //back_solve
      //come back to the inc case
      #ifdef BASKER_INC_LVL
       t_back_solve_offdiag_selective(kid,
           L_col, L_row,
           X_col, X_row,
           k, col_idx_offset,
           U.val,
           U.row_idx,
           U.col_ptr[k-bcol+1]-U.col_ptr[k-bcol],
           U.col_ptr[k-bcol],
           A_option);

      #else
       /*
          printf("t_bsolve_d test, kid: %d xsize: %d\n",
          kid, 
          U.col_ptr[k-bcol+1]-U.col_ptr[k-bcol]);
          */
       /*        
          t_back_solve_offdiag(kid,
          L_col, L_row,
          X_col, X_row,
          k, col_idx_offset,
          U.val,
          U.row_idx,
       //U.col_ptr[k-bcol+1]-U.col_ptr[k-bcol],
       U.col_ptr(k+1)-U.col_ptr(k),
       //U.col_ptr[k-bcol],
       U.col_ptr(k),
       A_option);
       */

       t_dense_back_solve_offdiag(kid,
           L_col, L_row,
           X_col, X_row,
           k, col_idx_offset,
           U.val,
           U.row_idx,
           //U.col_ptr[k-bcol+1]-U.col_ptr[k-bcol],
           U.col_ptr(k+1)-U.col_ptr(k),
           //U.col_ptr[k-bcol],
           U.col_ptr(k),
           A_option);

      #endif
     }//end for over all offdiag
     #endif // end of BASKER_MULTIPLE_UPPER
     
     #endif // end of BASKER_2DL

     //Bgood(removed)
     //B.flip_base();
     
     return info;
  }//end t_upper_col_factor()


  //Did not update for 2DL since we will not use
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_back_solve_atomic
  (
   Int kid,
   Int team_leader,
   Int lvl, 
   Int l,
   Int k, 
   Int top,
   Int xnnz
  )
  {

    Int            b = S(l)(kid);
    BASKER_MATRIX &L = LL(b)(0);
    INT_1DARRAY   ws = thread_array(kid).iws;
    ENTRY_1DARRAY  X = thread_array(team_leader).ews;
    Int      ws_size = thread_array(kid).iws_size;
    Int     ews_size = thread_array(team_leader).ews_size;
  
    #ifdef BASKER_DEBUG_NFACTOR_COL
    if(kid>3)
    printf("Back solve, kid: %d b: %d lvl: %d l: %d\n",
           kid, b, lvl,l);
    #endif

    Int *color    = &(ws[0]);
    Int *pattern  = &(color[ws_size]);

    Int top1 = top;
    Int j,t,pp, p, p2;
    Entry xj = 0;

    Int Lend = L.scol + L.ncol;

    for(pp = 0; pp < xnnz; pp++)
    {
      j = pattern[top1];
      top1++;
      color[j] = 0;
      t = gperm[j];

#ifdef BASKER_DEBUG_NFACTOR_COL
      if(kid>3)
        printf("Backsolve, j: %d t: %d \n",
            j, t);
#endif

      if((t >= L.scol) && (t < (Lend)))
      {
        //ATOMIC : NOT NEEDED
        xj = X[j];
        Int local_offset = t-L.scol;
        p2 = L.col_ptr[local_offset+1];
        p =  L.col_ptr[local_offset]+1;

#ifdef BASKER_DEBUG_NFACTOR_COL
        if(kid>3)
          printf("Atomic_tri_solve: Updating col: %d %d with val: %f \n",
              j, t, xj);
#endif

        for( ; p < p2; p++)
        {
          Int row_idx = L.row_idx[p];

          if(gperm[row_idx] < Lend)
          {  
#ifdef BASKER_DEBUG_NFACTOR_COL
            if(kid>3)
              printf("kid: %d, bs  row_idx: %d x: %f val: %f xj: %f \n",
                  kid, row_idx, X[row_idx], L.val[p], xj); 
#endif
            X[row_idx] -= L.val[p]*xj;
          }
          else
          {
#ifdef BASKER_DEBUG_NFACTOR_COL
            if(kid>3)
              printf("kid: %d, abs  row_idx: %d x: %f val: %f xj: %f \n",
                  kid, row_idx, X[row_idx], L.val[p], xj); 
#endif
            Kokkos::atomic_fetch_sub(&(X[L.row_idx[p]]),
                L.val[p]*xj);
          }
        }//end for() over all nnz in a column

      }//end if() not permuted

    }//end for() over all nnz in LHS
    return 0;
  }//end t_back_solve_atomic


  //uses local idx for local blks
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::t_lower_col_factor
  (
   Int kid,
   Int team_leader,
   Int lvl,
   Int l, 
   Int k, 
   Entry &opivot
  )
  {
    using STS = Teuchos::ScalarTraits<Entry>;
    using Mag = typename STS::magnitudeType;
    const Entry zero (0.0);
    const Entry one  (1.0);
    const Mag eps = STS::eps ();
    const Mag normA = BTF_A.gnorm;
    const Mag normA_blk = BTF_A.anorm;

    //Get needed variables
    const Int L_col = S(lvl)(kid);
    const Int L_row = 0;
    const Int U_col = S(lvl)(kid);
    const Int U_row = LU_size(U_col)-1;
    const Int X_col = S(0)(kid);
    //Int col_idx_offset = 0; //can we get rid of now?

    #ifdef BASKER_DEBUG_NFACTOR_COL
    printf("LOWER_COL_FACTOR kid: %d \n", kid);
    printf("kid %d using L (%d %d),  U (%d %d),  and X (%d) \n",
           kid, L_col, L_row, U_col, U_row, X_col);
    #endif
    //end get needed variables

    BASKER_MATRIX        &L = LL(L_col)(L_row);
    BASKER_MATRIX        &U = LU(U_col)(U_row); 
   
    BASKER_MATRIX        &B = thread_array(kid).C;
    
    #ifdef BASKER_DEBUG_NFACTOR_COL
    if(kid >= 0)
    {
      printf("---------------LOWER FACTOR KID: %d ------\n",
          kid);
      B.info();
      B.print();
      printf("After matrix print \n");
    }
    #endif
    //B.print();


    INT_1DARRAY   ws      = LL(X_col)(l+1).iws;
    const Int     ws_size = LL(X_col)(l+1).iws_size;
    ENTRY_1DARRAY X       = LL(X_col)(l+1).ews;

    Int   scol_top = btf_tabs[btf_top_tabs_offset]; // the first column index of A
    const Int brow_a  = U.srow; // offset within A
    const Int brow_g  = brow_a + scol_top; // global offset
    const Int lval  = L.col_ptr(k);
    const Int uval  = U.col_ptr(k);
    
    Int *color     = &(ws(0));
    Int *pattern   = &(color[ws_size]);

    Int i,j;
    Int top, maxindex, t;
    Int lnnz, unnz, xnnz, lcnt, ucnt;
    Int cu_ltop, cu_utop;

    Int newsize;
    Entry pivot, value;
    Mag absv, maxv, digv;
    Int llnnz = L.mnnz;
    Int uunnz = U.mnnz;
    
    cu_ltop = lval;
    cu_utop = uval;
    top = ws_size;
    
    lnnz = lval;
    unnz = uval;

    #ifdef BASKER_DEBUG_NFACTOR_COL
    printf("-------------lower_col_factor-------\n");
    printf("JB ecol: %d L.nnz: %d U.nnz: %d brow: %d ws_size: %d  kid: %d\n",
           U.scol+U.ncol, llnnz, uunnz, brow, ws_size, kid);
    #endif

    #ifdef MY_DEBUG_BASKER
    printf("\n----------------- t_lower_col_factor(K = %d, L = LL(%d,%d) with mnnz=%d, U = LU(%d,%d) with mnnz=%d) --------------\n", 
           k+U.scol,L_col,L_row,llnnz, U_col,U_row,uunnz);
    #endif
              
    value = zero;
    pivot = zero;

    maxindex = BASKER_MAX_IDX;
    lcnt = 0;
    ucnt = 0;
    
    #ifdef BASKER_DEBUG_NFACTOR_COL
    BASKER_ASSERT(top == ws_size);
    for(i = 0 ; i < ws_size; i++){
      if(X(i) !=0)
      {
        printf("--error, kid: %d x[%d] = %f \n", kid, i,X(i)); 
      }
      assert(x[i] == 0);
    }
    for(i = 0; i <  ws_size; i++)
    {assert(ws[i] == 0 );}
    #endif

    #ifdef MY_DEBUG_BASKER
    printf(" k=%d: B.col(k): %d  B.col(k+1): %d \n", 
            k,(int)B.col_ptr(0), (int)B.col_ptr(1));
    #endif

    //printf( "\n >> t_lower_col_factor (%d) <<\n",k );
    //for(i = B.col_ptr(0); i < B.col_ptr(1); ++i) printf( "%d %d %e\n",B.row_idx(i),k,B.val(i) );

    for(i = B.col_ptr(0); i < B.col_ptr(1); ++i)
    {
      j = B.row_idx(i);

      //Dont think we need this anymore ... do we?
      if(j > U.nrow)
      {
        printf("j continue -- do be need? \n");
        break;
      }

      //X[j-brow] = B.val(i);
      X(j) = B.val(i);
      #ifdef MY_DEBUG_BASKER
      printf(" t_lower_col_factor: Input(kid=%d)t: X(%d) = LL(%d, %d)(%d) = %e \n", kid, j, X_col,l+1,j, X(j));
      #endif

#ifdef BASKER_DEBUG_NFACTOR_COL
      if(kid>=0)
      {
        printf("i: %ld  j: %ld %ld  val: %g  top: %d \n", 
            i, j, gperm(j+brow_g), B.val(i), top);
      }
      #ifdef BASKER_2DL
      if(kid>=0)
        printf("Nxk in Ak %d %g color = %d \n",
            j, X[j], ws[j ] );
      #else
      if(kid>=0)
        printf("Nx in Ak %d %g color = %d \n",
            j, X[j], ws[j ] );
      #endif
#endif

      #ifdef BASKER_2DL
      if(color[j] == 0)
      #else
      if(color[j] == 0)
      #endif
      {
        //printf("doing reach: %d \n", j);
        //#ifdef BASKER_INC_LVL
        //t_local_reach_selective(kid,lvl,l+1,j, &top);
        //#else
        //t_local_reach(kid,lvl,l+1, j, &top);
        if(gperm(j+brow_g) != BASKER_MAX_IDX)
        {
          //printf("COL LOCAL REACH\n");
          t_local_reach(kid,lvl,l+1,j,top);
        }
        else
        {
          //printf("COL LOCAL SHORT\n");
          t_local_reach_short(kid,lvl,l+1,j,top);
        }
        //#endif
      }
    }//over each nnz in the column
    xnnz = ws_size - top;
   
    #ifdef MY_DEBUG_BASKER
    if(kid>=0)
    {
      printf("xnnz: %d ws_size: %d top: %d \n",
             (int)xnnz, (int)ws_size, (int)top);
    }
    #endif

    #ifdef BASKER_OPS_COUNT
    thread_array(kid).ops_counts[0][l] += xnnz;
    #endif
   
    t_back_solve(kid, lvl,l+1, k, top, xnnz); // note: l not lvl given

    //find pivot
    maxv = abs(zero);
    digv = abs(zero);
    for(i = top; i < ws_size; i++)
    { 
      j = pattern[i];
      //t = gperm[j];
      t = gperm(j+brow_g);
      //printf( " pattern[%d] = %d , gperm(%d+%d) = %d\n",i,j,brow_g,j,t );

      value = X(j);

      absv = EntryOP::approxABS(value);
      if(t == BASKER_MAX_IDX)
      {
        lcnt++;
        if(EntryOP::gt(absv,maxv) || maxindex == BASKER_MAX_IDX) // \gt, or the first time
        {
          maxv     = absv;
          pivot    = value;
          maxindex = j;
        }
        if(j == k) {
          digv = absv;
        }
      }
      #ifdef MY_DEBUG_BASKER
      {
        printf(" kid=%d: considering: j=%d value=%f, absv=%f, maxv=%f, pivot=%f, digv=%f (t=%d)\n", (int)kid, (int)j, value,absv,maxv,pivot,digv, (int)t);
      }
      #endif
    } //for (i = top; i < ws_size)          
    ucnt = ws_size - top - lcnt +1;
    #ifdef MY_DEBUG_BASKER
    {
      printf(" kid=%d:   ==> maxv=%f, pivot=%f, digv=%f\n", (int)kid, maxv,pivot,digv);
    }
    #endif
   
    //----------------------------Sym-----
    //SYM
    if(Options.no_pivot == BASKER_TRUE || digv > maxv * Options.pivot_tol)
    {
      maxindex = k;
      pivot = X(k);
      #ifdef MY_DEBUG_BASKER
      {
        printf( " kid=%d: using diagonal for pivot, pivot=%f\n",(int)kid,pivot);
      }
      #endif
    }

    #ifdef BASKER_DEBUG_NFACTOR_COL
    //if(kid>=0)
    printf("pivot found: %f , kid: %d (normA=%e -> tol=%e)\n", pivot, kid, normA_blk, normA_blk * sqrt(eps));
    #endif

    if((maxindex == BASKER_MAX_IDX) || (pivot == zero) )
    {
      if(Options.verbose == BASKER_TRUE)
      {
        std::cout << std::endl << std::endl;
        std::cout << "---------------------------" << std::endl;
        std::cout << "Error: Col Matrix is singular, col k = " << k
                  << ", lvl = " << lvl << ", l = " << l << std::endl;
        std::cout << "MaxIndex: " << maxindex << " pivot " << pivot << std::endl;
        std::cout << " norm(A)   = " << normA     << " (global)" << std::endl
                  << " norm(A)   = " << normA_blk << " (block)"  << std::endl
                  << " replace_zero_pivot = " << Options.replace_zero_pivot << std::endl;
        if (Options.replace_tiny_pivot && normA_blk > abs(zero) && maxindex != BASKER_MAX_IDX) {
          std::cout << "  + replace zero pivot with " << normA_blk * sqrt(eps) << std::endl;
        } else if (Options.replace_zero_pivot && normA_blk > abs(zero) && maxindex != BASKER_MAX_IDX) {
          std::cout << "  - replace zero pivot with " << normA_blk * eps << std::endl;
        }
        std::cout << "---------------------------" << std::endl;
      }
       
      if (Options.replace_tiny_pivot && normA_blk > abs(zero) && maxindex != BASKER_MAX_IDX) {
          pivot = normA_blk * sqrt(eps);
          X(maxindex) = pivot;
      } else if (Options.replace_zero_pivot && normA_blk > abs(zero) && maxindex != BASKER_MAX_IDX) {
        pivot = normA_blk * eps;
        X(maxindex) = pivot;
      } else {
        // replace-tiny-pivot not requested, or the current column is structurally empty after elimination
        thread_array(kid).error_type   = BASKER_ERROR_SINGULAR;
        thread_array(kid).error_blk    = L_col;
        thread_array(kid).error_subblk = -1;
        thread_array(kid).error_info   = k;
        return BASKER_ERROR;
      }
    } else if (Options.replace_tiny_pivot && normA_blk > abs(zero) && abs(pivot) < normA_blk * sqrt(eps)) {
      if (Options.verbose == BASKER_TRUE)
      {
        std::cout << std::endl << std::endl;
        std::cout << "---------------------------" << std::endl;
        std::cout << "Col Matrix : replace tiny pivot col k = " << k
                  << ", lvl = " << lvl << ", l = " << l << std::endl;
        std::cout << " pivot " << pivot << " -> "
                  << (STS::real(pivot) >= abs(zero) ? normA_blk * sqrt(eps) : -normA_blk * sqrt(eps))
                  << std::endl;
        std::cout << "---------------------------" << std::endl;
      }
      if (STS::real(pivot) >= abs(zero)) {
        pivot = normA_blk * sqrt(eps);
      } else {
        pivot = -normA_blk * sqrt(eps);
      }
      X(maxindex) = pivot;
    }
    U.tpivot = pivot;
    //printf("lower pivot: %e, k: %d, kid: %d \n", U.tpivot, k, kid);

  
    //gperm[maxindex] = k;
    gperm(maxindex+brow_g) = k+brow_g; 
    //gpermi[k] = maxindex;
    gpermi(k+brow_g) = maxindex+brow_g;
   
    if(lnnz + lcnt > llnnz)
    {
      //Note: comeback
      newsize = lnnz * 1.1 + 2 *A.nrow + 1;
      if (Options.verbose == BASKER_TRUE)
      {
        std::cout << "Lower Col Reallocing L oldsize: " << llnnz 
                  << " newsize: " << newsize << " kid = " << kid
                  << std::endl;
        //cout << " > k = " << k << " lnnz = " << lnnz << " lcnt = " << lcnt << endl;
        //cout << " > L_col = " << L_col << " L_row = " << L_row << endl;
      }

      thread_array(kid).error_blk    = L_col;
      thread_array(kid).error_subblk = -1;
      if(Options.realloc == BASKER_FALSE)
      {
        thread_array(kid).error_type = BASKER_ERROR_NOMALLOC;
        return BASKER_ERROR;
      }
      else
      {
        thread_array(kid).error_type = BASKER_ERROR_REMALLOC;
        thread_array(kid).error_info = newsize;
        return BASKER_ERROR;
      }
    }
    if(unnz+ucnt > uunnz)
    {
      //Note: comeback
      newsize = uunnz*1.1 + 2*A.nrow+1;
      if (Options.verbose == BASKER_TRUE)
      {
        std::cout << "Lower Col Reallocing U oldsize: " << uunnz 
                  << " newsize " << newsize << " kid = " << kid
                  << std::endl;
      }

      thread_array(kid).error_blk    = U_col;
      thread_array(kid).error_subblk = U_row;
      if(Options.realloc == BASKER_FALSE)
      {
        thread_array(kid).error_type = BASKER_ERROR_NOMALLOC;
        return BASKER_ERROR;
      }
      else
      {
        thread_array(kid).error_type = BASKER_ERROR_REMALLOC;
        thread_array(kid).error_info = newsize;
        return BASKER_ERROR;
      }
    }

    //printf("kid: %d lnnz: %d llnz: %d \n", 
    //kid, lnnz, llnnz);
    L.row_idx(lnnz) = maxindex;
    L.val(lnnz) = one;
    lnnz++;

    Entry lastU = zero;
  
    //For every nnz in LHS
    for( i = top; i < ws_size; i++)
    {
      j = pattern[i];
      //t = gperm[j];
      t = gperm(j+brow_g);

#ifdef BASKER_DEBUG_NFACTOR_COL
      if(k>=0)
        printf("j %d  t %d , kid: %d \n", j, t, kid);
#endif

      //if not zero
      #ifdef BASKER_2DL
      if(X(j) != zero)
      #else
      if(X[j] != zero)
      #endif
      {
#ifdef BASKER_DEBUG_NFACTOR_COL
         #ifdef BASKER_2DL
         if(kid>=0)
         {
           printf("found value: %f at %d, kid: %d \n",
               X[j], j, kid);
         }
         #else
         if(kid>=0)
           printf("found value: %f at %d, kid: %d \n",
               X[j], j, kid);
         #endif
#endif

         if(t != BASKER_MAX_IDX)
         {
           if(t < k+brow_g)
           {
#ifdef BASKER_DEBUG_NFACTOR_COL
             //if(kid>=0)
             if (L_col == 2 && kid == 0)
             {
             #ifdef BASKER_2DL
               printf("U insert: %e %d at %d \n",
                   X[j], t-brow_g, unnz);
             #else
               printf("U insert: %e %d at %d \n",
                   X[j], t-brow_g, unnz);
             #endif
             }
#endif
             //U.row_idx[unnz] = gperm[j];
             //can't we reuse, this seems
             //stupid
             U.row_idx(unnz) = t-brow_g;
             #ifdef BASKER_2DL
             U.val(unnz) = X(j);
             #else
             U.val[unnz] = X[j];
             #endif
             unnz++;
           }
           else
           {
             lastU = X(j);
           }
         }
         else if (t == BASKER_MAX_IDX)
         {
#ifdef BASKER_DEBUG_NFACTOR_COL
           #ifdef BASKER_2DL
           if(kid>=0)
           {
             // printf("inserting %f at %d into %d \n", 
             //          X[j-brow]/pivot, j, lnnz );
             printf("inserting %f at %d into %d \n", 
                 X[j]/pivot, j, lnnz );
           }
           #else
           if(kid>=0)
           {
             //printf("inserting %f at %d into %d \n", 
             //          X[j]/pivot, j, lnnz );
             printf("inserting %f at %d into %d \n", 
                 X[j]/pivot, j, lnnz );
           }
           #endif
#endif
           L.row_idx(lnnz) = j;
           #ifdef BASKER_2DL
           L.val(lnnz) = EntryOP::divide(X(j), pivot);
           #else
           L.val(lnnz) = EntryOP::divide(X(j), pivot);
           #endif
           lnnz++;
         }
      }//end if() not zero             
      //move inside if not zero..... extra set ops not needed

      #ifdef MY_DEBUG_BASKER
      printf(" t_lower_col_factor(kid=%d)t: X(%d) = ZERO\n",kid,j);
      #endif
      #ifdef BASKER_2DL
      X(j) = zero;
      #else
      X[j] = zero;
      #endif
    }//if(x[j-brow] != 0)
   
    //Fill in last element of U
    U.row_idx(unnz) = k;
    U.val(unnz) = lastU;
#ifdef BASKER_DEBUG_NFACTOR_COL
    if (L_col == 2 && kid == 0)
    {
      printf( " Insert U(last) = %e %d at %d\n",lastU,k,unnz );
    }
#endif
    unnz++;
   
    xnnz = 0;
    top = ws_size;
   
    #ifdef BASKER_DEBUG_NFACTOR_COL
    printf("setting col: k=%d, cu_ltop=%d, lnnz=%d\n",
            k, cu_ltop, lnnz);
    #endif
    L.col_ptr(k) = cu_ltop;
    L.col_ptr(k+1) = lnnz;
    cu_ltop = lnnz;
   
    U.col_ptr(k) = cu_utop;
    U.col_ptr(k+1) = unnz;
    cu_utop = unnz;

    #ifdef BASKER_DEBUG_NFACTOR_COL
    if(kid>=0)
    printf("col_fact k: %d Unnz: %d   Lnnz: %d \n",
           k, unnz, lnnz);
    #endif
  
    #ifndef BASKER_MULTIPLE_LOWER
    #ifdef BASKER_2DL
    //----Update offdiag (this can be done in parallel l>0)---//
    //-Need to fix
    #ifdef BASKER_DEBUG_NFACTOR_COL
    printf("k: %d -----TTTEST(blk_row-- %d %d \n",
           kid, lvl+1, LL_size[L_col]);
    printf("k: %d -----TTTEST(x_row)--- %d %d \n",
           kid, l+2, LL_size[X_col]);
    #endif

    #ifdef MY_DEBUG_BASKER
    printf(" t_lower_col_factor(kid=%d)t: calling t_dense_back_solve_offdiag(blk_row = %d:%d)\n", kid,L_row+1,LL_size(L_col)-1 );
    #endif
    for(Int blk_row = L_row+1, x_row = l+2; blk_row < LL_size(L_col); blk_row++, x_row++)
    { 
      /*old
       t_back_solve_offdiag(kid,
       L_col, blk_row,
       X_col, x_row,
       k, col_idx_offset,
       U.val, U.row_idx,
       //U.col_ptr[k-bcol+1]-U.col_ptr[k-bcol],
       U.col_ptr(k+1)-U.col_ptr(k),
       //U.col_ptr[k-bcol],
       U.col_ptr(k),
       BASKER_TRUE);
      */
      t_dense_back_solve_offdiag(kid,
          L_col, blk_row,
          X_col, x_row,
          k, col_idx_offset,
          U.val, U.row_idx,
          //U.col_ptr[k-bcol+1]-U.col_ptr[k-bcol],
          U.col_ptr(k+1)-U.col_ptr(k),
          //U.col_ptr[k-bcol],
          U.col_ptr(k),
          BASKER_TRUE);

      t_dense_move_offdiag_L(kid, 
          L_col, blk_row,
          X_col, x_row,
          k, pivot);
    }//end for over all offdiag blks
    #endif
    #endif

    t_prune(kid, lvl, l+1, k, maxindex);

    #ifdef BASKER_DEBUG_NFACTOR_DEBUG
    print_factor(L,U);
    #endif

   //Bgood(remove)
   //B.flip_base();
   
    return 0;
  }//end t_lower_col_fact()


  //This has been replace with t_basker_barrier in thread.hpp
  template <class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::t_col_barrier(Int kid)
  {
    return 0;
  }//end t_col_barrier()


  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::t_basker_barrier
  (
   const TeamMember &thread,
   const Int my_kid,
   const Int leader_kid, 
   const Int size,
   const Int function_n,
   const Int k, 
   const Int l
  )
  {
    //printf("before call. lkid: %d kid: %d task: %d size: %d k: %d \n",
    //               leader_kid, my_kid, function_n, size, k);

    #ifdef HAVE_VTUNE
    //__itt_pause();
    #endif

    //printf( " t_basker_barrier(size = %d, thread.team_size = %d: my_id = %d, leader_id = %d)\n",size,thread.team_size(), my_kid,leader_kid ); fflush(stdout);
    //if(size < 0)
    if(size <= thread.team_size())
    {
      thread.team_barrier();
    }
    else
    {
      basker_barrier.BarrierDomain(leader_kid,
         my_kid, 
         function_n,
         size, 
         k, 
         l );
    }

    #ifdef HAVE_VTUNE
    //__itt_resume();
    #endif

  }//end t_basker_barrier


  template <class Int, class Entry,class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::t_basker_barrier_old
  (
   const TeamMember &thread,
   const Int leader_kid,
   const Int sublvl,
   const Int function_n, 
   const Int size
  )
  {

    if(size <= thread.team_size())
    {
      thread.team_barrier();
    }
    else
    {

      //basker_barrier.Barrier(leader_kid, 

      /* Old Atomic Barrier
         BaskerBarrier<Int,Entry,Exe_Space> BB;
         BB.Barrier(thread_array(leader_kid).token[sublvl][function_n],
         thread_array(leader_kid).token[sublvl][1],
         size);
         */
    }
  }//end t_basker_barrier()

}//end namespace BaskerNS--functions

#undef BASKER_DEBUG_NFACTOR_COL
#endif //end ifndef backer_nfactor_col
