// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLUBASKER_NFACTOR_COL2_HPP
#define SHYLUBASKER_NFACTOR_COL2_HPP

//#include "shylubasker_decl.hpp"
#include "shylubasker_matrix_decl.hpp"
#include "shylubasker_matrix_view_decl.hpp"
#include "shylubasker_matrix_view_def.hpp"
#include "shylubasker_types.hpp"
#include "shylubasker_stats.hpp"
#include "shylubasker_thread.hpp"

#include "shylubasker_nfactor_blk.hpp"
#include "shylubasker_nfactor_col.hpp"

#ifdef BASKER_KOKKOS
#include <Kokkos_Core.hpp>
#include <Kokkos_Timer.hpp>
#endif 

#ifdef BASKER_DEBUG
#include <assert.h>
#endif

#ifdef HAVE_VTUNE
#include "ittnotify.h"
#endif

//#define BASKER_DEBUG_NFACTOR_COL2
//#define BASKER_DEBUG_TIME
//#define BASKER_COUNT_OPS
//#define BASKER_TIMER

namespace BaskerNS
{
  template <class Int, class Entry, class Exe_Space>
  struct kokkos_nfactor_sep2
  {
    #ifdef BASKER_KOKKOS
    typedef Exe_Space                         execution_space;
    typedef Kokkos::TeamPolicy<Exe_Space>     TeamPolicy;
    typedef typename TeamPolicy::member_type  TeamMember;
    #endif
    
    Basker<Int,Entry,Exe_Space> *basker;
    Int lvl;
    
    kokkos_nfactor_sep2()
    {}
    
    kokkos_nfactor_sep2(
                Basker<Int,Entry,Exe_Space>* _basker, Int _lvl)
    {
      basker = _basker;
      lvl = _lvl;
    }

    BASKER_INLINE
    #ifdef BASKER_KOKKOS
    void operator()(const TeamMember &thread) const
    #else
    void operator()(Int kid) const  
    #endif
    {
      #ifdef BASKER_KOKKOS
      //Int kid = (Int)(thread.league_rank()*thread.team_size()+
      //                thread.team_rank());
      Int kid = basker->t_get_kid(thread);

      //team leader is not using
      Int team_leader = (Int)(thread.league_rank()*thread.team_size());
      #else
      Int team_leader = 0; //Note: come back and fix
      #endif

      #ifdef HAVE_VTUNE
      __itt_pause();
      #endif

      {
      #ifdef BASKER_KOKKOS
        basker->t_nfactor_sep2(kid, lvl, team_leader, thread);
      #else
      
      #endif
      }

      #ifdef HAVE_VTUNE
      __itt_resume();
      #endif
    }//end operator ()
  };//end col_factor_funct


  template <class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::t_nfactor_sep2
  (
   const Int kid,
   const Int lvl,
   const Int team_leader,
   const TeamMember &thread
   )
  {
    // this routine performs left-looking factorization of the diagonabl block LU(U_col)(U_row)
    //
    // before the call to nfactor_sep: off-diagonal blocks of L have been already factored
    // as part of nfactor_blk
    //
    // so, this routine:
    // 1) compute the previous blocks U(1:U_col-1)(U_row)
    // 2) update LU(U_col)(U_row) using the previous diagonal blocks (l < lvl)
    // 3) factor LU(U_col)(U_row)
    // 4) apply U(U_col)(U_row)^{-1} to the blocks below A(U_col)(U_row)
    //
    // For instance, we factor A(7,7):
    //
    // LU(1,1)         U(1,3)                        A(1,7)
    //        LU(2,2)  U(2,3)                        A(2,7)
    //  L(3,1) L(3,2) LU(3,3)                        A(3,7)
    //                        LU(4,4)         U(4,6) A(4,7)
    //                               LU(5,5)  U(5,6) A(5,7)
    //                         L(6,4) L(6,5) LU(6,6) A(6,7)
    //  L(7,1) L(7,2)  L(7,3)  L(7,4) L(7,5)  L(7,6) A(7,7)
    //
    // 1) t_upper_col_factor          : apply L-solve to A(:,7)
    // 2) t_add_extend                : accumulate update into A(7,7)
    // 3) t_lower_col_factor          : factor A(7,7), sequential
    // 4) t_lower_col_factor_offdiag2 : compute L(8:end, 7)

    const Int U_col = S(lvl)(kid);
    const Int U_row = 0;
    Int ncol = LU(U_col)(U_row).ncol;
    Int my_leader = find_leader(kid, 0);
    if(Options.verbose == BASKER_TRUE && kid == my_leader)
    {
      printf(" > kid = %ld(%ld): factoring_col current_chunk: lvl = %ld size=%ld\n",
            (long)kid, (long)my_leader, (long)lvl, (long)ncol);
    }

    #ifdef BASKER_TIMER
    Kokkos::Timer timer;
    #endif

    #ifdef BASKER_DEBUG_NFACTOR_COL2
    printf("\n\n  LVL=%d  ----- kid: %d -----\n\n", lvl, kid);
    #endif
    // > Apply lower-triangular solve with L(1,1)
    //   to compute the first block U(1)(U_row) in the upper-triangular part
    int info = BASKER_SUCCESS;
    for(Int k = 0; k < ncol && info == BASKER_SUCCESS; ++k)
    {
      #ifdef BASKER_DEBUG_NFACTOR_COL2
      printf("UPPER, kid: %d k: %d \n", kid, k);
      #endif
      info = t_upper_col_factor(kid, team_leader,
                                lvl, 0,
                                k,
                                BASKER_FALSE);
    }//over all columns / domains / old sublevel 0

    #ifdef BASKER_TIMER
    printf("Time Upper-Col1(%d): %lf \n", (int)kid, timer.seconds()); fflush(stdout);
    timer.reset();
    #endif
    //------Need because extend does not 
    //-------------Barrier--Between Domains-------------
    //barrier k = 0 usedl1
    //#define USE_TEAM_BARRIER_NFACTOR_COL2
    #ifdef USE_TEAM_BARRIER_NFACTOR_COL2
    thread.team_barrier();
    #else
    Int b_size = pow(2, 1);
    t_basker_barrier(thread, kid, my_leader,
                     b_size, 0, LU(U_col)(U_row).scol, 0);
    for(Int tid = 0; tid < num_threads; tid++) {
      if (thread_array(tid).error_type != BASKER_SUCCESS) {
        info = BASKER_ERROR;
      }
    }
    #endif
    #ifdef BASKER_DEBUG_NFACTOR_COL2
    printf("\n done with 1st UPPER, kid: %d \n\n", kid); fflush(stdout);
    #endif


    //----------------Sep level upper tri-------------
    for(Int l = 1; l < (lvl) && info == BASKER_SUCCESS; ++l)
    {
      for(Int k = 0; k < ncol; ++k)
      {
        #ifdef BASKER_DEBUG_NFACTOR_COL2
        printf("\nSep, upper update, kid: %d k=%d \n",
               kid, k+LU(U_col)(U_row).scol);
        #endif

        // > Accumulate the update from (l-1)th level:
        //    LU(U_col)(U_row) -= L(U_col)(l-1) * U(l-1)(U_row)
        t_add_extend(thread, kid, lvl, l-1, k,
                     LU(U_col)(U_row).scol,
                     BASKER_FALSE);

        if(kid%((Int)pow(2, l)) == 0)
        {
          //my_leader = find_leader(kid,l);
          //b_size    = pow(2,l+1);
          #ifdef BASKER_DEBUG_NFACTOR_COL2
          printf("\n SEP UPPER, kid: %d \n",
                 kid);
          #endif

          // > Apply lower-triangular solve with next diagonal L(l,l)
          //   to compute next U(l)(U_row) in upper-triangular part
          info = t_upper_col_factor(kid, team_leader, 
                                    lvl, l, 
                                    k,
                                    BASKER_FALSE);
        }//if correct kid to do this sublevels upperfactor
      }//over all columns
      //Trick to just add 1 in case of sep size == 1
      //t_basker_barrier(thread, kid, kid, 
      //                 1, 1, LU(U_col)(U_row).ncol+1, l-1);
      //t_basker_barrier(thread, kid, kid,
      //                 1, 2, LU(U_col)(U_row).ncol+1, l-1);

    }//for - over all sublevel 1...lvl-2
    #ifdef BASKER_TIMER
    printf("Time Upper-Col(%d): %lf \n", (int)kid, timer.seconds());
    timer.reset();
    #endif

    //---------Lower Factor (old sublevel lvl-1)-------
    
    //printf("\n\n");
    //printf("lower team size: %d \n", thread.team_size());
    #ifdef USE_TEAM_BARRIER_NFACTOR_COL2
    thread.team_barrier();
    #else
    my_leader = find_leader(kid, lvl-1);
    b_size    = pow(2, lvl);
    // printf("[3] barrier test, kid: %d leader: %d b_size: %d lvl: %d \n",
    //        kid,  my_leader, b_size, lvl);
    t_basker_barrier(thread, kid, my_leader,
                     b_size, 3, LU(U_col)(U_row).scol, 0);
    for(Int ti = 0; ti < num_threads; ti++) {
      if (thread_array(kid).error_type != BASKER_SUCCESS) {
        info = BASKER_ERROR;
      }
    }
    #endif
    #ifdef BASKER_DEBUG_NFACTOR_COL2
    printf("\n done with UPPER, kid: %d \n\n", kid);
    #endif

    //printf("\n\n======= LOWER, KID: %d ======= \n\n", kid);
    //return;
    // > accumulate the last update
    // > factor the diagonal block LU(U_col)(U_row)
    // > apply U-solve with U(U_col)(U_row)
    //   to compute off-diagonal L(:)(U_row)
    {
      #ifdef BASKER_TIMER
      double time_extend = 0.0;
      double time_faccol = 0.0;
      double time_facoff = 0.0;
      Kokkos::Timer timer_extend;
      Kokkos::Timer timer_faccol;
      Kokkos::Timer timer_facoff;
      #endif
      for(Int k = 0; k < ncol && info == BASKER_SUCCESS; ++k)
      {
        // ------------------------------------------------------- //
        // > accumulate the last update into k-th column of LU(U_col)(U_row)
        #ifdef BASKER_TIMER
        timer_extend.reset();
        #endif
        if (info == BASKER_SUCCESS)
        {
          #ifdef BASKER_DEBUG_NFACTOR_COL2
          printf( " kid=%d: calling t_add_extend(k=%d/%d)\n",kid,k,ncol ); fflush(stdout);
          #endif
          t_add_extend(thread, kid,lvl,lvl-1, k,
                       LU(U_col)(U_row).scol,
                       BASKER_TRUE);
        }
        #ifdef BASKER_TIMER
        time_extend += timer_extend.seconds();
        #endif
        #ifdef BASKER_DEBUG_NFACTOR_COL2
        printf( " kid=%d: done calling t_add_extend(k=%d/%d)\n",kid,k,ncol ); fflush(stdout);
        #endif

        // ------------------------------------------------------- //
        // > factor the k-th column of LU(U_col)(U_row)
        #ifdef BASKER_TIMER
        timer_faccol.reset();
        #endif
        Entry pivot (0.0);
        if((kid%(Int)(pow(2,lvl))) == 0)
        {
          #ifdef BASKER_DEBUG_NFACTOR_COL2
          printf(" > calling lower factor, kid: %d k: %d \n", kid, k); fflush(stdout);
          #endif
          if (info == BASKER_SUCCESS) {
            // factor the kth column of LU(U_col)(U_row)
            info = t_lower_col_factor(kid, team_leader, 
                                      lvl, lvl-1, 
                                      k, pivot);
          }
        }
        #ifdef BASKER_DEBUG_NFACTOR_COL2
        printf(" > done calling lower factor, kid: %d k: %d info=%d\n", kid, k, info); fflush(stdout);
        #endif
        #ifdef BASKER_DEBUG_NFACTOR_COL2
        else {
          printf(" + skipping lower factor, kid: %d k: %d \n", kid, k); fflush(stdout);
        }
        #endif
        //need barrier if multiple thread uppdate
        #ifdef USE_TEAM_BARRIER_NFACTOR_COL2
        thread.team_barrier();
        #else
        my_leader = find_leader(kid, lvl-1);
        b_size    = pow(2, lvl);
        #ifdef BASKER_DEBUG_NFACTOR_COL2
        printf("barrier test-4: kid = %d, k = %d, leader = %d, b_size = %d, lvl = %d \n",
               kid, k, my_leader, b_size, lvl); fflush(stdout);
        #endif
        t_basker_barrier(thread, kid, my_leader,
                         b_size, 4, k, lvl-1);
        for(Int tid = 0; tid < num_threads; tid++) {
          if (thread_array(tid).error_type != BASKER_SUCCESS) {
            info = BASKER_ERROR;
          }
        }
        #ifdef BASKER_DEBUG_NFACTOR_COL2
        printf("barrier test-4 done: kid = %d, k = %d, leader = %d, b_size = %d, lvl = %d, info = %d \n",
               kid, k, my_leader, b_size, lvl, info); fflush(stdout);
        #endif
        #endif
        #ifdef BASKER_TIMER
        time_faccol += timer_faccol.seconds();
        #endif

        // ------------------------------------------------------- //
        // > factor the k-th column of the off-diagonal blocks
        if (info == BASKER_SUCCESS) {
          #ifdef BASKER_TIMER
          timer_facoff.reset();
          #endif
          #ifdef BASKER_DEBUG_NFACTOR_COL2
          printf(" calling lower diag factor, kid: %d k: %d \n",
                 kid, k); fflush(stdout);
          #endif
          t_lower_col_factor_offdiag2(kid, lvl, lvl-1, k, pivot);
          #ifdef BASKER_DEBUG_NFACTOR_COL2
          printf(" done lower diag factor, kid: %d k: %d \n",
                 kid, k); fflush(stdout);
          #endif
        }

        #ifdef USE_TEAM_BARRIER_NFACTOR_COL2
        thread.team_barrier();
        #else
        my_leader = find_leader(kid, lvl-1);
        b_size    = pow(2, lvl);
        #ifdef BASKER_DEBUG_NFACTOR_COL2
        printf("barrier test-5: kid = %d, k = %d, leader = %d, b_size = %d, lvl = %d \n",
               kid, k, my_leader, b_size, lvl); fflush(stdout);
        #endif
        t_basker_barrier(thread, kid, my_leader,
                         b_size, 5, k, lvl-1);
        #ifdef BASKER_DEBUG_NFACTOR_COL2
        printf("barrier test-5 done: kid = %d, k = %d, leader = %d, b_size = %d, lvl = %d \n",
               kid, k, my_leader, b_size, lvl); fflush(stdout);
        #endif
        #endif
        #ifdef BASKER_TIMER
        time_facoff += timer_facoff.seconds();
        #endif
      }

      //Trick for sep = 1
      //t_basker_barrier(thread, kid, kid,
      //                 1, 1, LU(U_col)(U_row).ncol+1, lvl-1);
      //t_basker_barrier(thread, kid, kid,
      //                 1, 1, LU(U_col)(U_row).ncol+1, lvl-1);
      #ifdef BASKER_TIMER
      double time_factot = timer.seconds();
      if((kid%(Int)(pow(2,lvl))) == 0) {
        const Int L_col = S(lvl)(kid);
        const Int L_row = LU_size(U_col)-1;

        printf("Time Lower-Col(%d): %lf, n = %d, nnz(L) = %d, nnz(U) = %d \n", (int)kid, time_factot,
               (int)ncol, (int)LL(U_col)(U_row).col_ptr(ncol), (int)LU(L_col)(L_row).col_ptr(ncol));
        printf(" > Time Lower-Col(%d):extend-add: %lf \n", (int)kid, time_extend);
        printf(" > Time Lower-Col(%d):fac-col   : %lf \n", (int)kid, time_faccol);
        printf(" > Time Lower-Col(%d):fac-off   : %lf \n", (int)kid, time_facoff);
      }
      #endif
    }

    //printf( " >> kid=%d: returning with t_nfactor_sep2 with info = %d\n",kid,info );
    return info;
  }//end t_nfactor_sep2

  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::t_add_extend
  (
   const TeamMember &thread,
   const Int kid, 
   const Int lvl,
   const Int l,
   const Int k,
   const Int k_offset,
   const BASKER_BOOL lower
   )
  {
    Int my_leader = find_leader(kid,l);
    #ifndef USE_TEAM_BARRIER_NFACTOR_COL2
    Int b_size    = pow(2, l+1);
    #endif

    //loop over sublevels to perform 
    //off-dig multiple, reduce, etc
    for(Int sl = 0; sl <= l; ++sl)
    {
      #ifdef BASKER_DEBUG_NFACTOR_COL2
      if(lower == BASKER_FALSE)
      {
        printf("extend up, kid: %d sl: %d l: %d lvl: %d \n",
               kid, sl, l, lvl);
      }
      else
      {
        printf("extend low, kid: %d sl: %d l: %d lvl: %d \n", 
               kid, sl, l, lvl);
      }
      #endif

      //This will do the correct spmv
      if(thread_array(kid).error_type == BASKER_ERROR_NOERROR) {
        t_upper_col_factor_offdiag2(kid, lvl, sl,l, k, lower);
      }
      //Barrier--Start
      #ifdef USE_TEAM_BARRIER_NFACTOR_COL2
      thread.team_barrier();
      #else
      my_leader = find_leader(kid,sl);
      b_size    = pow(2,sl+1);
      t_basker_barrier(thread, kid, my_leader,
                       b_size, 1, k+k_offset, sl);
      #endif
      //Barrier--End

      if(kid%((Int)pow(2,sl)) == 0 &&
         thread_array(kid).error_type == BASKER_ERROR_NOERROR) {
        t_dense_blk_col_copy_atomic2(kid, my_leader,
                                     lvl, sl, l, k, lower);
      }

      //Barrier--Start
      //printf("[2] Barrier test, kid: %d leader: %d k: %d sl: %d \n",
      //     kid, my_leader, k, sl);
      #ifdef USE_TEAM_BARRIER_NFACTOR_COL2
      thread.team_barrier();
      #else
      t_basker_barrier(thread, kid, my_leader,
                       b_size, 2, k+k_offset, sl);
      #endif
    }//over all sublevels

    if(thread_array(kid).error_type == BASKER_ERROR_NOERROR) {
      t_dense_copy_update_matrix2(kid, my_leader, lvl, l, k);
    }
  }//end t_add_add

  //uses local idx for local blks X
  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::t_upper_col_factor_offdiag2
  (
   const Int kid,
   const Int lvl,
   const Int sl,
   const Int l,
   const Int k, 
   const BASKER_BOOL lower
   )
  {
    const Int my_leader = (sl==0)?(kid):find_leader(kid,sl-1);
    if(kid != my_leader)
    {
      /*
      if(lower == BASKER_TRUE)
      {
        printf("offd, kid: %d my_leader: %d \n",
               kid, my_leader);
      }
      */
      return;
    }

    Int my_row_leader = S(0)(find_leader(kid,lvl-1));
    const Int L_col = S(sl)(my_leader);
    const Int U_col = S(lvl)(kid);
    const Int X_col = S(0)(my_leader);
    Int L_row = l-sl+1; //Might have to think about th
    Int U_row = L_col-my_row_leader;
    Int X_row = l+1; //this will change for us 

    BASKER_MATRIX &U = LU(U_col)(U_row);
    #ifdef BASKER_DEBUG_NFACTOR_COL2
    if(L_row >= LL_size(L_col))
    {
      printf("assert error. Lrow:%d L_col: %d size: %d kid: %d\n",
             L_row, L_col, LL_size(L_col), kid);
    }
    BASKER_ASSERT(L_row < LL_size(L_col), "upper-off, Lrow >= size");
    BASKER_ASSERT(X_row < LL_size(X_col), "upper-off, Xrow >=size"); 
    
    //if(lower == BASKER_TRUE)
    {
      printf("Upper_fact_offdiag, kid: %d leader: %d l: %d lvl: %d works_size: %d X: %d %d L: %d %d U: %d %d k: %d \n",
             kid, my_leader, l, lvl, LL_size[X_col], X_col, X_row, L_col, L_row, U_col, U_row,  k+U.scol);
      printf("OFF-DIAG, kid: %d, l: %d  X: %d %d L: %d %d \n",
             kid, l,X_col, X_row, L_col, L_row);
    }
    #endif

    // backward-solve on first block
    Int col_idx_offset  = 0;
    t_dense_back_solve_offdiag(kid,
                               L_col, L_row,
                               X_col, X_row,
                               k, col_idx_offset,
                               U.val,
                               U.row_idx,
                               U.col_ptr(k+1)-U.col_ptr(k),
                               U.col_ptr(k),
                               BASKER_FALSE);
    
    //if lower, finish off the updates
    if(lower == BASKER_TRUE)
    {
      X_row++;
      L_row++;
      for(; X_row < LL_size(X_col); ++X_row, ++L_row)
      {
        #ifdef BASKER_DEBUG_NFACTOR_COL2
        printf("LLL OFF-DIAG,kid:%d, l: %d X: %d %d L: %d %d U: %d %d \n",
               kid, l,X_col, X_row, L_col, L_row, U_col, U_row);
        #endif

        t_dense_back_solve_offdiag(kid,
                                   L_col, L_row,
                                   X_col, X_row,
                                   k, col_idx_offset,
                                   U.val,
                                   U.row_idx,
                                   U.col_ptr(k+1)-U.col_ptr(k),
                                   U.col_ptr(k),
                                   BASKER_FALSE);
      }//for --over remainder blks
    }//if --lower
  
  }//end t_upper_col_factor_offdiag2()

 
  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry, Exe_Space>::t_dense_blk_col_copy_atomic2
  (
   const Int kid,
   const Int NOT_USED,
   const Int lvl,
   const Int sl,
   const Int l,
   const Int k,
   const BASKER_BOOL lower
   )
  {
    //Setup
    //printf("DEBUG, kid: %d k: %d A_col: %d A_row: %d \n", 
    //       kid, k, A_col, A_row);
    const Int my_idx     = S(0)(kid);
    //should remove either as a paramter or here
    Int team_leader      = find_leader(kid, sl);
    const Int leader_idx = S(0)(team_leader);
    #ifdef BASKER_DEBUG_NFACTOR_COL2
    if(lower == BASKER_TRUE)
    {
      printf("Called t_blk_col_copy_atomic kid: %d \n " , kid);
      printf("Copying col, kid: %d  k: %d lvl: %d l: %d \n", 
             kid,k,lvl, l);
      //printf("Copying Col, kid: %d k:%d  A: %d %d to tl: %d li: %d\n",
      //        kid, k, A_col, A_row, team_leader, leader_idx);
    }
    #endif
   
    //If I an not a leader, then need to copy over
    if(kid != team_leader)
    {
      Int endblk = (lower)?(LL_size(my_idx)):(l+2);
      for(Int blk = l+1; blk < endblk; ++blk)
      {
        ENTRY_1DARRAY &XL = LL(leader_idx)(blk).ews;
        Int      p_sizeL  = LL(leader_idx)(blk).p_size;
        ENTRY_1DARRAY &X  = LL(my_idx)(blk).ews;
        INT_1DARRAY   &ws = LL(my_idx)(blk).iws;
        Int       *color  = &(ws[0]);
        //printf( " + t_dense_blk_col_copy_atomic2(kid=%d: LL(%d)(%d) += LL(%d)(%d)\n",kid,leader_idx, blk,my_idx,blk);

        #ifdef BASKER_DEBUG_NFACTOR_COL2
        if(lower == BASKER_TRUE)
        {
          printf("kid: %d  COPY INDEX %d %d to %d %d \n",
                 kid, my_idx, blk, leader_idx, blk);
          printf("t_b_col_copy, kid: %d wsize: %d \n", 
                 kid, ws_size);
          printf("t_b_col_copy,kid:%d ps: %d XL: %d %d X: %d %d\n",
                 kid, p_size, leader_idx, blk, my_idx, blk);
        }
        #endif

        //over all nnnz found
        for(Int jj = 0; jj < LL(my_idx)(blk).nrow; ++jj)
        {
          color[jj] = 0;
          #ifdef BASKER_DEBUG_NFACTOR_COL2
          if(lower == BASKER_TRUE)
          {
            printf("kid: %d jj: %d nrows: %d xsize: %d %d \n",
                   kid, jj, LL(my_idx)(blk).nrow, 
                   LL(my_idx)(blk).iws_size,
                   LL(leader_idx)(blk).iws_size);

            printf("Adding X(%d)%f to XL(%d) %f, leader: %d kid: %d\n",
                   jj+brow, X[jj], jj, XL[jj], team_leader, kid);

          }
          #endif

          //if(X(jj) != (Entry)(0) )
          {
            #ifdef BASKER_DEBUG_NFACTOR_COL2
            if(lower == BASKER_TRUE)
            {
              printf("Atomic Adding X(%d) %f to XL(%d) %f, kid: %d\n",
                     jj+brow, X[jj], jj, XL[jj], kid);
            }
            #endif

            XL(jj) += X(jj);
            X(jj)   = 0;
          }//if X(j) != 0
        }//for - over all nnz
   
        #ifdef BASKER_DEBUG_NFACTOR_COL2
        if(lower == BASKER_TRUE)
        {
          printf("----DEBUG--- X %d %d brow: %d kid: %d ws_size: %d\n",
                 my_idx, blk, brow, kid, ws_size);
          for(Int j = 0; j< ws_size; ++j)
          {
            printf("X[%d,%d] = %f , ptern: %d color: %d kid: %d \n",
                   j, j+brow, X[j], pattern[j], color[j],kid);
          }
        }
        #endif

        //This can be removed in the future
        if(kid != team_leader)
        {
          LL(my_idx)(blk).p_size = 0;
        }
        else
        {
          #ifdef BASKER_DEBUG_NFACTOR_COL22
          printf("SETTING PS: %d L:%d %d kid: %d\n",
                 p_sizeL, leader_idx, blk, kid);
          #endif
          LL(leader_idx)(blk).p_size = p_sizeL;
          //p_size = 0; //not needed
        }//over all blks
      }
    }//if not team_leader
  }//end t_dense_blk_col_copy_atomic2()


  //local idx local blk
  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::t_dense_copy_update_matrix2
  (
   const Int kid,
   const Int team_leader,
   const Int lvl,
   const Int l, 
   const Int k
   )
  {
    //printf("\n\n\n\n");
    //printf("-----------------copy_update_matrx----------");
    //printf("\n\n\n\n");

    const Entry zero (0.0);
    const Int leader_idx = S(0)(kid);
    BASKER_MATRIX     &C = thread_array(kid).C;  
    Int nnz = 0;

    //Over each blk    
    /*
    if(lvl ==(l+1))
    {
      last_blk = LL_size(leader_idx);
    }
    */
    
    // X += B(:, k)
    {
      Int bl = l+1;
      Int A_col = S(lvl)(kid);

      Int my_row_leader = S(0)(find_leader(kid,lvl-1));
      Int A_row = S(bl)(kid) - my_row_leader;

      BASKER_MATRIX *Bp;
      if(A_row != (LU_size(A_col)-1))
      {
        //printf("upper picked, kid: %d \n", kid);
        //printf("up: %d %d kid: %d \n",
        //       A_col, A_row, kid);
        Bp = &(AVM(A_col)(A_row));
      }
      else
      {
        //printf("lower picked, kid: %d\n", kid);
        Bp = &(ALM(A_col)(0));
      }
      #ifdef BASKER_DEBUG_NFACTOR_COL2
      printf("copy, kid: %d bl: %d  A: %d %d \n", 
             kid, bl, A_col, A_row);
      #endif

      // X += B(:, k)
      BASKER_MATRIX &B = *Bp;
      ENTRY_1DARRAY  X = LL(leader_idx)(bl).ews;
      //printf( " -- t_dense_copy_update_matrix2(kid=%d: LL(%d)(%d) += B)\n",kid,leader_idx,bl );
      //printf("ADDING UPDATES TO B\n");
      //B.info();
      //B.print();
      for(Int i = B.col_ptr(k); i < B.col_ptr(k+1); ++i)
      {
        Int B_row = B.row_idx(i);

        #ifdef BASKER_DEBUG_NFACTOR_COL22
        printf("Scanning_2 A: %d %d lvl %d l: %d bl:%d brow: % d %d K: %d \n",
               B_row, j, lvl, l, bl, brow, B.srow, kid);
        printf("Adding Aval: %f to xval: %f \n", 
               X[B_row], B.val(i));
        #endif
 
        X(B_row) += B.val(i);
      }//end for over all nnz
    }

    //-------------move into C------------------- 
    //(Right now have to do dense but like to do sparse)

    //last_blk = LL_size(leader_idx); //NDE - warning: unused 
    //printf("------maybe l:%d lastblk: %d kid: %d\n",
    //   l, last_blk, kid);
    
    // copy C = X
    {
      Int bl = l+1;
  
      /*
      Int A_row = (lvl==1)?(2):S(bl)(kid)%(LU_size(A_col));
      //maybe no???
      if((S(bl)(kid) > 14) &&
         (S(bl)(kid) > LU_size(A_col)) &&
         (lvl != 1))
      {
        //printf("test cm %d %d %d \n",
        //     kid, S(bl)(kid), LU_size(A_col));

        Int tm = (S(bl)(kid)+1)/16;
        A_row  = ((S(bl)(kid)+1) - (tm*16))%LU_size(A_col);
      }

      // printf("kid: %d leader_idx: %d bl: %d \n",
      //        kid, leader_idx, bl);
      */

      //For recounting patterns in dense blk
      //Need better sparse update
      ENTRY_1DARRAY   X   = LL(leader_idx)(bl).ews;
      INT_1DARRAY    ws   = LL(leader_idx)(bl).iws;
      const Int      nrow = LL(leader_idx)(bl).nrow;
      Int *color   = &(ws(0));
      #ifdef BASKER_DEBUG_NFACTOR_COL2
      printf("moving, kid: %d  A: %d %d %d %d p_size: %d \n", 
             kid, A_col, A_row, team_leader, bl,p_size);
      #endif

      //over all dim(S)
      for(Int jj=0; jj < nrow; ++jj)
      {
        //Int j = pattern[jj];
        Int j = jj;
        #ifdef BASKER_DEBUG_NFACTOR_COL22
        printf("considering: %d %d %f, kid: %d\n",
               j,brow,X[j], kid);
        #endif

        if (X(j) != zero)
        {
          if(bl == l+1)
          {
            #ifdef BASKER_DEBUG_NFACTOR_COL22
            printf("moving X[%d] %f kid: %d nnz: %d Csize: %d \n",
                   j+brow, X(j), kid, nnz, C.nrow);
            #endif
            C.row_idx(nnz) = j;
            C.val(nnz)     = X(j);
            nnz++;
            X(j)           = zero;
            color[j] = 0;
          }
          else
          {
            #ifdef BASKER_DEBUG_NFACTOR_COL22
            printf("counting [%d] %f kid: %d \n",
                   j+brow, X[j], kid);
            #endif
            BASKER_ASSERT("1==0", "copy, should never call");
          }
        }//if X(j) != 0
      }//for -- length of x
    }
    C.col_ptr(0) = 0;
    C.col_ptr(1) = nnz;

    //printf("Done with move, kid: %d found nnz: %d \n",
    //   kid, nnz);
    
    //C.info();
    //C.print();

  }//end t_dense_copy_update_matrix2()


  //local idx and lock blk
  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::t_lower_col_factor_offdiag2
  (
   const Int kid,
   const Int lvl,
   const Int l,
   const Int k,
   Entry pivot
   )
  {
    #ifndef BASKER_MULTIPLE_LOWER
    BASKER_ASSERT(0==1,"BASKER_MULTIPLE_LOWER ERROR");
    return 0;
    #endif

    const Int leader_id   = find_leader(kid, l);
    const Int lteam_size  = pow(2,l+1);

    const Int L_col       = S(lvl)(leader_id);
    const Int U_col       = S(lvl)(leader_id);

    Int L_row             = 0;
    Int U_row             = LU_size(U_col)-1;

    Int X_col             = S(0)(leader_id);
    Int X_row             = l+1;

    Int col_idx_offset    = 0;  //can get rid of?
   
    BASKER_MATRIX        &U = LU(U_col)(U_row); 
    pivot = U.tpivot;
    
    //BASKER_MATRIX        &L = LL(L_col)(L_row); //NDE - warning: unused L
    //const Int  ws_size = LL(X_col)(X_row).iws_size;
    //INT_1DARRAY     ws = LL(X_col)(X_row).iws;
    //ENTRY_1DARRAY    X = LL(X_col)(X_row).ews;

    //const Int brow     = U.srow;
    //const Int bcol     = U.scol;
    
    //printf("OFF_DIAG_LOWER, kid: %d leaderid: %d t_size: %d \n",
    //       kid, leader_id, lteam_size);
    
    L_row += (kid-leader_id)+1;
    X_row += (kid-leader_id)+1;
    for( ; 
         L_row < LL_size(L_col);
         X_row+=(lteam_size), L_row+=(lteam_size))
    {
      //printf("OFF_DIAG_LOWER. kid: %d k: %d  U: %d %d L: %d %d X: %d %d pivot: %f \n", kid, k, U_col, U_row, L_col, L_row, X_col, X_row, pivot);
      /*old
       t_back_solve_offdiag(leader_id,
                            L_col, L_row,
                            X_col, X_row,
                            k, col_idx_offset,
                            U.val, U.row_idx,
                            //U.col_ptr[k-bcol+1]-U.col_ptr[k-bcol],
                            U.col_ptr(k+1)-U.col_ptr(k),
                            //U.col_ptr[k-bcol],
                            U.col_ptr(k),
                            BASKER_TRUE);
      */

      //We might still have to do sparse here
      t_dense_back_solve_offdiag(leader_id,
                                 L_col, L_row,
                                 X_col, X_row,
                                 k, col_idx_offset,
                                 U.val, U.row_idx,
                                 U.col_ptr(k+1)-U.col_ptr(k),
                                 U.col_ptr(k),
                                 BASKER_TRUE);

      //printf("MOVING OFF, kid: %d k: %d L: %d %d X %d %d \n",
      //     kid, k, L_col, L_row, X_col, X_row);
      t_dense_move_offdiag_L(leader_id, 
                             L_col, L_row,
                             X_col, X_row,
                             k, pivot);
    }//end for over all offdiag blks
  }//end t_lower_col_factor_offdiag2()
}//end namespace BaskerNS

#undef BASKER_TIMER
#endif //end ifndef
