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


//#define MY_DEBUG_BASKER
//#define BASKER_DEBUG_NFACTOR_COL
//#define BASKER_DEBUG_TIME
//#define BASKER_COUNT_OPS

namespace BaskerNS
{
  template <class Int, class Entry, class Exe_Space>
  struct kokkos_nfactor_sep
  {
    #ifdef BASKER_KOKKOS
    typedef Exe_Space                         execution_space;
    typedef Kokkos::TeamPolicy<Exe_Space>     TeamPolicy;
    typedef typename TeamPolicy::member_type  TeamMember;
    #endif
    
    Basker<Int,Entry,Exe_Space> *basker;
    Int lvl;
    
    kokkos_nfactor_sep()
    {}
    
    kokkos_nfactor_sep(Basker<Int,Entry,Exe_Space>* _basker, Int _lvl)
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
      //              thread.team_rank());
      Int kid = basker->t_get_kid(thread);

      //team leader is not using
      Int team_leader = (Int)(thread.league_rank()*thread.team_size());
      #else
      Int team_leader = 0; //Note: come back and fix
      #endif


      #ifdef HAVE_VTUNE
      __itt_pause();
      #endif


      //if(kid ==0 || kid ==1)
      {
      #ifdef BASKER_KOKKOS
      basker->t_nfactor_sep(kid, lvl, team_leader, thread);
      #else
      basker->t_nfactor_sep(kid, lvl, team_leader);
      #endif
      }

      #ifdef HAVE_VTUNE
      __itt_resume();
      #endif
            
    }//end operator ()
  };//end col_factor_funct

   //old using global index for local blk
  #ifdef BASKER_KOKKOS
  template <class Int, class Entry, class Exe_Space>
  int Basker<Int, Entry,Exe_Space>::t_nfactor_sep
  (
   Int kid,
   Int lvl, 
   Int team_leader,
   const TeamMember &thread
  )
  #else
  template <class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::t_nfactor_sep
  (
   Int kid,
   Int lvl,
   Int team_leader
  )
  #endif
  {
    #ifdef BASKER_DEBUG_TIME
    double upper_time = 0; 
    double upper_diag_time = 0;
    double lower_time = 0;
    double lower_diag_time = 0;
    double reduce_time = 0;
    double copy_time = 0;
    double barrier_time = 0;
    #endif
    
    Int U_col = S(lvl)(kid);
    Int U_row = 0;
    
    const Int scol = LU(U_col)(U_row).scol;
    const Int ecol = LU(U_col)(U_row).ecol;
    const Int ncol = LU(U_col)(U_row).ncol;

    //for(Int k = scol; k < ecol; k++)
    //might have to use k+scol for barrier
    for(Int k = 0; k < ncol; ++k)
    {
      #ifdef BASKER_DEBUG_NFACTOR_COL
      printf("-------------------k=%d--------------------\n",
          k+LU(U_col)(U_row).scol);
      #endif

      for(Int l = 0; l < lvl; ++l) //sublevel
      {
        #ifdef BASKER_DEBUG_NFACTOR_COL
        printf("\n\n--------sl=%d-------------\n\n", l);
        #endif

        //Used for barrier
        Int my_leader = find_leader(kid, l);
        Int b_size = pow(2, l+1);

        //Can remove in future, don't need anymore
        BASKER_BOOL sep_lvl_flg = BASKER_FALSE;

        #ifdef BASKER_DEBUG_TIME
        Kokkos::Timer timer;
        #endif

        //-----------Upper Col Factor------------//
        if(kid%((Int)pow(2,l)) == 0)
        {
          #ifdef BASKER_DEBUG_NFACTOR_COL
          printf("------STARTING TRI-SOLVE: %d %d -----\n", kid, l);
          #endif


          #ifdef HAVE_VTUNE
          __itt_resume();
          #endif

          t_upper_col_factor(kid, team_leader, lvl, l, k, sep_lvl_flg);

          #ifdef HAVE_VTUNE
          __itt_pause();
          #endif
        }//end if kid ... I should do an upper_col_factor
        #ifdef BASKER_DEBUG_TIME
        upper_time += timer.seconds();
        timer.reset();
        #endif


        #ifdef BASKER_MULTIPLE_UPPER
         //Upper_col_factor_updates
         my_leader = (l == 0 ? kid : find_leader(kid,l-1));
         b_size = pow(2,l);

         //printf("t_upper barrier, kid: %d b_size: %d \n",
         //           kid, b_size);

         //printf("Barrier0. kid: %d lkid: %d bsize: %d l: %d\n",
         //           kid, my_leader, b_size, l);

         //t_basker_barrier(thread,my_leader,l,0,b_size);
         t_basker_barrier(thread, kid, my_leader, b_size, 0, k, l);

         //printf("Leave Barrier0 kid: %d \n", kid);

         #ifdef BASKER_DEBUG_TIME
         barrier_time += timer.seconds();
         timer.reset();
         #endif


         #ifdef HAVE_VTUNE
         __itt_resume();
         #endif

         t_upper_col_factor_offdiag(kid,lvl,l,k);
         my_leader = find_leader(kid,l);
         b_size = pow(2,l+1);

         #ifdef HAVE_VTUNE
         __itt_pause();
         #endif


         #ifdef BASKER_DEBUG_TIME
         upper_diag_time += timer.seconds();
         timer.reset();
         #endif

        #else // else not defined BASKER_MULTIPLE_UPPER

         my_leader = find_leader(kid,l);
         b_size = pow(2,l+1);
        #endif // end of BASKER_MULTIPLE_UPPER


        #ifdef BASKER_OLD_BARRIER
        thread.team_barrier();
        #else
        //Int my_leader = find_leader(kid, l);
        //Int b_size = pow(2,l+1);
        //Debug

        //printf("Barrier1, kid: %d lkid: %d bsize: %d l: %d \n",
        //           kid, my_leader, b_size, l);
        //t_basker_barrier(thread, my_leader, l, 1, b_size);
        t_basker_barrier(thread, kid, my_leader, b_size, 1, k, l);

        //printf("Leave Barrier1: %d \n", kid);
        #endif


        #ifdef BASKER_DEBUG_TIME
        barrier_time += timer.seconds();
        timer.reset();
        #endif


#ifdef BASKER_2DL
        if(kid%((Int)pow(2,l))==0)
        {
          #ifdef HAVE_VTUNE
          __itt_resume();
          #endif

          //Rename---not really atomic anymore
          //if((kid==0)||(kid==1)||(kid==2)||(kid==3))

          //t_blk_col_copy_atomic(kid, team_leader, lvl, 
          //                      l, k);

          t_dense_blk_col_copy_atomic(kid, team_leader, lvl, l, k);


          #ifdef HAVE_VTUNE
          __itt_pause();
          #endif
        }
#endif

        #ifdef BASKER_DEBUG_TIME
        reduce_time += timer.seconds();
        timer.reset();
        #endif

        #ifdef BASKER_OLD_BARRIER
        thread.team_barrier();
        #else
        //printf("Barrier2, kid: %d lkid: %d bsize: %d l: %d\n",
        //           kid, my_leader, b_size, l);
        //t_basker_barrier(thread, my_leader, l, 2,b_size);
        t_basker_barrier(thread, kid, my_leader, b_size, 2, k,l);

        //printf("Leave Barrier2 kid: %d \n", kid);
        #endif


        #ifdef BASKER_DEBUG_TIME
        barrier_time += timer.seconds();
        timer.reset();
        #endif


        #ifdef BASKER_ATOMIC_2 //O(dim(S)*p)
        //t_n_col_copy_atomic(kid, team_leader, lvl, l, k);
        #endif

        //#ifdef BASKER_KOKKOS
        //thread.team_barrier();//might have to be added back
        //#else
        //t_col_barrier(kid,lvl,l);
        //#endif

        if(kid%((Int)pow((double)2,(double)l+1)) == 0)
        {  
          #ifdef HAVE_VTUNE
          __itt_resume();
          #endif

          #ifdef BASKER_ATOMIC
          //t_col_copy_atomic(kid, team_leader, lvl, l, k);
          //if((kid==0)||(kid==1)||(kid==2)||(kid==3))
          //t_copy_update_matrix(kid, team_leader,lvl, l, k);
          t_dense_copy_update_matrix(kid, team_leader,lvl, l, k);
          #endif


          #ifdef HAVE_VTUNE
          __itt_pause();
          #endif
        }//end if(kid) should I do reduce


        #ifdef BASKER_DEBUG_TIME
        copy_time += timer.seconds();
        timer.reset();
        #endif

        //Used if O(dim(S)*num_threads) Atomic
        #ifdef BASKER_ATOMIC_2
        //thread.team_barrier();
        #endif

        //-----------------lower factor-----------------//  
        Entry pivot = 0;
        if(((l+1)==lvl)&&
            (kid%((Int)pow((double)2,(double)l+1)) == 0))
        {
          #ifdef HAVE_VTUNE
          __itt_resume();
          #endif

          //if((kid==0)||(kid==1)||(kid==2)||(kid==3))
          t_lower_col_factor(kid, team_leader,
              lvl, l, k, pivot);

          #ifdef HAVE_VTUNE
          __itt_pause();
          #endif

        }//end if I should perfom factor
        #ifdef BASKER_DEBUG_TIME
        lower_time += timer.seconds();
        timer.reset();
        #endif

        //printf("Barrier3, kid: %d lkid: %d bsize: %d l: %d\n", kid, my_leader, b_size, l);

        t_basker_barrier(thread,kid, my_leader, b_size, 3, k, l);

        //printf("Leave Barrier3, kid: %d \n", kid);

#ifdef BASKER_MULTIPLE_LOWER
        if((l+1) == lvl)
        {
          //need a barrier to make sure upper_lower done
          //printf("Barrier: lower_diag off-diag \n");
          //t_basker_barrier(thread, my_leader, l, 3,
          //                     b_size);

          #ifdef BASKER_DEBUG_TIME
          barrier_time += timer.seconds();
          timer.reset();
          #endif

          #ifdef HAVE_VTUNE
          __itt_resume();
          #endif

          t_lower_col_factor_offdiag(kid,lvl,l,k, pivot);

          #ifdef HAVE_VTUNE
          __itt_pause();
          #endif

          #ifdef BASKER_DEBUG_TIME
          lower_diag_time += timer.seconds();
          timer.reset();
          #endif
        }
#endif

        //We might be able to remove this in the future
        #ifdef BASKER_OLD_BARRIER
        thread.team_barrier();
        #else
        //my_leader = 0;
        //b_size = num_threads;
        //printf("Barrier4, kid: %d lkid: %d bsize: %d l: %d\n",
        //           kid, my_leader, b_size, l);
        //t_basker_barrier(thread, my_leader, l,4, b_size);
        t_basker_barrier(thread, kid, my_leader, b_size, 4, k, l);

        //printf("Leave Barrier4 kid: %d \n", kid);
        #endif

        #ifdef BASKER_DEBUG_TIME
        barrier_time += timer.seconds();
        timer.reset();
        #endif
      }// end over all sublvl
    }//over each column


    #ifdef BASKER_DEBUG_TIME
    printf("COL-TIME KID: %d  UPPER: %f UPPER-DIAG: %f LOWER: %f LOWER_DIAG: %f REDUCE: %f COPY: %f BARRIER: %f \n", 
           kid, upper_time, upper_diag_time, lower_time, lower_diag_time, reduce_time, copy_time, barrier_time);
    #endif

    #ifdef BASKER_OPS_COUNT
    for(Int l = 0; l < lvl; l++)
    {
      printf("OPS.   KID : %d  LVL: %d OPS : %d \n",
          kid, l, thread_array[kid].ops_counts[l][0]);
      thread_array[kid].ops_count[1][0] = 0;
    }
    #endif

    return 0;
  }//end t_nfactor_sep()


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
      Bp = &(thread_array[kid].C);
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
    thread_array[kid].ops_counts[0][l] += xnnz;
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


  //uses local idx for local blks
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_upper_col_factor_offdiag
  (
   Int kid,
   Int lvl,
   Int l,
   Int k
  )
  {
    //printf("t_upper_col_factor_offdiag called kid: %d \n", 
    //           kid);
   
    //Note: We can remove this and pass it
    //Might be faster
    Int my_leader = kid ;
    if(l>0)
    {
      my_leader = find_leader(kid, l-1);
    }
    //Note: We can remove this and pass it
    //Might be faster
    //hard coding for 2 right now
    int lteam_size = pow(2, l);
    
    #ifdef BASKER_2DL
    Int L_col = S(l)(my_leader);
    Int L_row = 0;
    Int U_col = S(lvl)(kid);
    Int U_row = (lvl==1)?(kid%2):S(l)(kid)%LU_size(U_col);
    Int X_col = S(0)(my_leader);
    Int X_row = l; //this will change for us 
    Int col_idx_offset = 0;
    BASKER_MATRIX        &U = LU(U_col)(U_row);
    const Int bcol = U.scol;
    #else
    BASKER_ASSERT(0==1, "t_upper_col_factor_offdiag, only work with with 2D layout");
    #endif

    #ifdef BASKER_2DL
    INT_1DARRAY ws     = LL[X_col][X_row].iws;
    const Int ws_size  = LL[X_col][X_row].iws_size;
    ENTRY_1DARRAY X    = LL[X_col][X_row].ews;
    #else
    BASKER_ASSERT(0==1, "t_upper_col_factor_offdiag, only works with 2D layout");
    #endif
    

    #ifdef BASKER_DEBUG_NFACTOR_COL
    printf("Upper_fact_offdiag, kid: %d leader: %d l: %d lvl: %d ltsize: %d works_size: %d X: %d %d L_row: %d %d \n",
           kid, my_leader, l, lvl,lteam_size, LL_size[X_col], X_col, X_row, L_col, L_row);
    #endif

     //----------------------Update offdiag-----------------//
    //Want to assign in a round-robin fashion

    //X_row++;
    //L_row++;
    X_row += (kid-my_leader)+1;
    L_row += (kid-my_leader)+1;
    for(; X_row < LL_size(X_col); X_row+=lteam_size, L_row+=lteam_size)
    {

      #ifdef BASKER_DEBUG_NFACTOR_COL
      printf("OFF-DIAG, kid: %d, l: %d  X: %d %d L: %d %d \n",
          kid, l,X_col, X_row, L_col, L_row);
      #endif

      const BASKER_BOOL A_option = BASKER_FALSE;

      //printf("\n\n size: %d \n\n",
      //U.col_ptr(k+1)-U.col_ptr(k));

      //back_solve 
      /*old
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


    }//end for over all offdiag     

    return BASKER_SUCCESS;
  }//end t_upper_col_factor_offdiag()


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

    Int            b = S[l][kid];
    BASKER_MATRIX &L = LL[b][0];
    INT_1DARRAY   ws = thread_array[kid].iws;
    ENTRY_1DARRAY  X = thread_array[team_leader].ews;
    Int      ws_size = thread_array[kid].iws_size;
    Int     ews_size = thread_array[team_leader].ews_size;
  
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
#ifdef BASKER_KOKKOS
            Kokkos::atomic_fetch_sub(&(X[L.row_idx[p]]),
                L.val[p]*xj);
#else
            //Note: come back
#endif         
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
    thread_array[kid].ops_counts[0][l] += xnnz;
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
        cout << endl << endl;
        cout << "---------------------------" << endl;
        cout << "Error: Col Matrix is singular, col k = " << k
             << ", lvl = " << lvl << ", l = " << l << endl;
        cout << "MaxIndex: " << maxindex << " pivot " << pivot << endl;
        cout << " norm(A)   = " << normA     << " (global)" << endl
             << " norm(A)   = " << normA_blk << " (block)"  << endl
             << " replace_zero_pivot = " << Options.replace_zero_pivot << endl;
        if (Options.replace_tiny_pivot && normA_blk > abs(zero) && maxindex != BASKER_MAX_IDX) {
          cout << "  + replace zero pivot with " << normA_blk * sqrt(eps) << endl;
        } else if (Options.replace_zero_pivot && normA_blk > abs(zero) && maxindex != BASKER_MAX_IDX) {
          cout << "  - replace zero pivot with " << normA_blk * eps << endl;
        }
        cout << "---------------------------" << endl;
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
        cout << endl << endl;
        cout << "---------------------------" << endl;
        cout << "Col Matrix : replace tiny pivot col k = " << k
             << ", lvl = " << lvl << ", l = " << l << endl;
        cout << " pivot " << pivot << " -> "
             << (STS::real(pivot) >= abs(zero) ? normA_blk * sqrt(eps) : -normA_blk * sqrt(eps))
             << endl;
        cout << "---------------------------" << endl;
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
        cout << "Lower Col Reallocing L oldsize: " << llnnz 
             << " newsize: " << newsize << " kid = " << kid
             << endl;
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
        cout << "Lower Col Reallocing U oldsize: " << uunnz 
             << " newsize " << newsize << " kid = " << kid
             << endl;
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


  //local idx and lock blk
  template <class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::t_lower_col_factor_offdiag
  (
   Int kid,
   Int lvl,
   Int l,
   Int k,
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
    Int L_row             = 0;
    const Int U_col       = S(lvl)(leader_id);
    Int U_row             = LU_size(U_col)-1;
    Int X_col             = S(0)(leader_id);
    Int X_row             = l+1;
    Int col_idx_offset    = 0;  //can get rid of?
   
    BASKER_MATRIX        &L = LL(L_col)(L_row);
    BASKER_MATRIX        &U = LU(U_col)(U_row); //U.fill();
    
    INT_1DARRAY     ws = LL(X_col)(X_row).iws;
    const Int  ws_size = LL(X_col)(X_row).iws_size;
    ENTRY_1DARRAY    X = LL(X_col)(X_row).ews;

    const Int bcol     = U.scol;

    pivot = U.tpivot;

    //printf("\n\n\n");
    printf(" t_lower_col_factor_offdiag( LL(%d)(%d) and LU(%d)(%d)) \n",L_col,L_row, U_col,U_row );
    //printf("lower_off, kid: %d leader_id: %d lsize: %d X: %d %d U %d %d L: %d %d \n",
    //           kid, leader_id, lteam_size, 
    //           X_col, X_row, U_col, U_row, L_col, L_row);
    //printf("\n\n\n");

    L_row += (kid-leader_id)+1;
    X_row += (kid-leader_id)+1;
    for( ; 
        L_row < LL_size(L_col);
        X_row+=(lteam_size), L_row+=(lteam_size)
       )
    { 
      // printf("OFF_DIAG_LOWER. kid: %d U: %d %d L: %d %d X: %d %d pivot: %f \n", kid, U_col, U_row, L_col, L_row, X_col, X_row, pivot);

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
      t_dense_back_solve_offdiag(leader_id,
          L_col, L_row,
          X_col, X_row,
          k, col_idx_offset,
          U.val, U.row_idx,
          //U.col_ptr[k-bcol+1]-U.col_ptr[k-bcol],
          U.col_ptr(k+1)-U.col_ptr(k),
          //U.col_ptr[k-bcol],
          U.col_ptr(k),
          BASKER_TRUE);


      /*
         t_move_offdiag_L(leader_id, 
         L_col, L_row,
         X_col, X_row,
         k, pivot);
         */
      t_dense_move_offdiag_L(leader_id, 
          L_col, L_row,
          X_col, X_row,
          k, pivot);

    }//end for over all offdiag blks

    return 0;
  }//end t_lower_col_factor_offdiag()


  //This has been replace with t_basker_barrier in thread.hpp
  template <class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::t_col_barrier(Int kid)
  {
    return 0;
  }//end t_col_barrier()


  //Used for 2DL O(dim(S)*p) (No more atomics)
  //Now done in ||
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry, Exe_Space>::t_dense_blk_col_copy_atomic
  (
   Int kid,
   Int team_leader,
   Int lvl,
   Int l,
   Int k
  )
  {
    //printf("\n\n\n\n");
    //printf("---------------blk_col_copy_atomic----------");
    //printf("\n\n");


    //Setup
    Int A_col = S(lvl)(kid);
    Int A_row = (lvl==1)?(2):S(l+1)(kid)%(LU_size(A_col));
   
    BASKER_MATRIX &B = AVM(A_col)(A_col);

    const Int my_idx     = S(0)(kid);
    team_leader = find_leader(kid, l);
    const Int leader_idx = S(0)(team_leader);
    Int loop_col_idx = S(l)(kid);

    #ifdef BASKER_DEBUG_NFACTOR_COL
    printf("Called t_blk_col_copy_atomic kid: %d " , kid);
    printf("Copying col, kid: %d  lvl: %d l: %d \n", 
           kid,lvl, l);
    printf("Copying Col, kid: %d  A: %d %d to tl: %d li: %d\n",
           kid, A_col, A_row, team_leader, leader_idx);
    #endif
   
#ifdef BASKER_2DL
    //If I an not a leader, then need to copy over
    if(kid != team_leader)
    {
      //over all blks
      //Split over threads (leader and nonleader)
      for(Int blk=l+1; blk<LL_size(my_idx); blk++)
      {
        ENTRY_1DARRAY &XL = LL(leader_idx)(blk).ews;
        INT_1DARRAY  &wsL = LL(leader_idx)(blk).iws;
        Int      p_sizeL  = LL(leader_idx)(blk).p_size;
        ENTRY_1DARRAY &X  = LL(my_idx)(blk).ews;
        INT_1DARRAY   &ws = LL(my_idx)(blk).iws;
        const Int ws_size = LL(my_idx)(blk).iws_size;
        Int       p_size  = LL(my_idx)(blk).p_size;
        Int       *color  = &(ws[0]);
        Int     *pattern  = &(color[ws_size]); 
        Int      brow     = LL(my_idx)(blk).srow;
        Int      browL    = LL(leader_idx)(blk).srow;

#ifdef BASKER_DEBUG_NFACTOR_COL
        printf("kid: %d  COPY INDEX %d %d to %d %d \n",
            kid, my_idx, blk, leader_idx, blk);
#endif


        Int *colorL   = &(wsL(0));
        Int *patternL = &(colorL[ws_size]);

#ifdef BASKER_DEBUG_NFACTOR_COL
        printf("t_b_col_copy, kid: %d wsize: %d \n", 
            kid, ws_size);
        printf("t_b_col_copy, kid: %d ps: %d XL: %d %d X: %d %d\n",
            kid, p_size, leader_idx, blk, my_idx, blk);
#endif


        //over all nnnz found

        //for(Int j=0; j < p_size; j++)
        for(Int jj = 0; jj < LL(my_idx)(blk).nrow; ++jj)
        {
          //Int jj = pattern[j];
          //color[jj-brow] = 0;
          color[jj] = 0;

#ifdef BASKER_DEBUG_NFACTOR_COL
          //printf("Adding j: %d X(%d) %f to XL(%d) %f, leader: %d kid: %d\n",
          //      j, jj, X[jj-brow], jj, XL[jj-browL], team_leader, kid);

          printf("Adding j: %d X(%d) %f to XL(%d) %f, leader: %d kid: %d\n",
              jj, X[jj], jj, XL[jj], team_leader, kid);

#endif

          //Note: may want to change to vary access pattern
          //if(kid != team_leader)
          {
            //if(X[jj-brow]!=0)
            if(X(jj) !=0)
            {
#ifdef BASKER_DEBUG_NFACTOR_COL

              printf("Atomic Adding X(%d) %f to XL(%d) %f, kid: %d\n",
                  jj, X[jj], jj, XL[jj], kid);
#endif


              XL(jj) += X(jj);
              X(jj)   = 0;

              /*
              //if(colorL[jj-browL] != 1)
              if(color[jj] != 1)
              {
              //colorL[jj-browL] = 1;
              colorL[jj] = 1;
              //patternL[p_sizeL++] = jj;
              patternL[p_sizeL++] = jj;
              }
              */

            }//if X(j) != 0


          }//no usec
        }//end over all nnz


#ifdef BASKER_DEBUG_NFACTOR_COL
        printf("----DEBUG--- X %d %d brow: %d kid: %d\n",
            my_idx, blk, brow, kid);
        for(Int j = 0; j< ws_size; j++)
        {

          //printf("X[%d] = %f , pttern: %d color: %d kid: %d \n", 
          //     j+brow, X[j], pattern[j], color[j],kid);
          printf("X[%d] = %f , pttern: %d color: %d kid: %d \n", 
              j+brow, X[j], pattern[j], color[j],kid);

        }

#endif


        if(kid != team_leader)
        {
          //LL[my_idx][blk].p_size = 0;
          LL(my_idx)(blk).p_size = 0;
        }
        else
        {
#ifdef BASKER_DEBUG_NFACTOR_COL
          printf("SETTING PS: %d L:%d %d kid: %d\n",
              p_sizeL, leader_idx, blk, kid);
#endif
          //LL[leader_idx][blk].p_size = p_sizeL;
          LL(leader_idx)(blk).p_size = p_sizeL;
        }
        p_size = 0;
      }//over all blks
    }//if not team_leader
   
#endif

    return 0;
  }//end t_dense_blk_col_copy_atomic()


  //Used for 2DL O(dim(S)*p) (No more atomics)
  //Now done in ||
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry, Exe_Space>::t_blk_col_copy_atomic
  (
   Int kid,
   Int team_leader,
   Int lvl,
   Int l,
   Int k
  )
  {
    //printf("\n\n\n\n");
    //printf("---------------blk_col_copy_atomic----------");
    //printf("\n\n");


    //Setup
    Int A_col = S(lvl)(kid);
    Int A_row = (lvl==1)?(2):S(l+1)(kid)%(LU_size(A_col));

    BASKER_MATRIX &B = AVM(A_col)(A_col);

    const Int my_idx     = S(0)(kid);
    team_leader = find_leader(kid, l);
    const Int leader_idx = S(0)(team_leader);
    Int loop_col_idx = S(l)(kid);

#ifdef BASKER_DEBUG_NFACTOR_COL
    printf("Called t_blk_col_copy_atomic kid: %d " , kid);
    printf("Copying col, kid: %d  lvl: %d l: %d \n", 
        kid,lvl, l);
    printf("Copying Col, kid: %d  A: %d %d to tl: %d li: %d\n",
        kid, A_col, A_row, team_leader, leader_idx);
#endif

#ifdef BASKER_2DL
    //If I an not a leader, then need to copy over
    if(kid != team_leader)
    {
      //over all blks
      //Split over threads (leader and nonleader)
      for(Int blk=l+1; blk<LL_size(my_idx); blk++)
      {
        ENTRY_1DARRAY &XL = LL(leader_idx)(blk).ews;
        INT_1DARRAY  &wsL = LL(leader_idx)(blk).iws;
        Int      p_sizeL  = LL(leader_idx)(blk).p_size;
        ENTRY_1DARRAY &X  = LL(my_idx)(blk).ews;
        INT_1DARRAY   &ws = LL(my_idx)(blk).iws;
        const Int ws_size = LL(my_idx)(blk).iws_size;
        Int       p_size  = LL(my_idx)(blk).p_size;
        Int       *color  = &(ws[0]);
        Int     *pattern  = &(color[ws_size]); 
        Int      brow     = LL(my_idx)(blk).srow;
        Int      browL    = LL(leader_idx)(blk).srow;

#ifdef BASKER_DEBUG_NFACTOR_COL
        printf("kid: %d  COPY INDEX %d %d to %d %d \n",
            kid, my_idx, blk, leader_idx, blk);
#endif

        Int *colorL   = &(wsL(0));
        Int *patternL = &(colorL[ws_size]);

#ifdef BASKER_DEBUG_NFACTOR_COL
        printf("t_b_col_copy, kid: %d wsize: %d \n", 
            kid, ws_size);
        printf("t_b_col_copy, kid: %d ps: %d XL: %d %d X: %d %d\n",
            kid, p_size, leader_idx, blk, my_idx, blk);
#endif

        //over all nnnz found
        for(Int j=0; j < p_size; j++)
        {
          Int jj = pattern[j];
          //color[jj-brow] = 0;
          color[jj] = 0;

#ifdef BASKER_DEBUG_NFACTOR_COL
          //printf("Adding j: %d X(%d) %f to XL(%d) %f, leader: %d kid: %d\n",
          //      j, jj, X[jj-brow], jj, XL[jj-browL], team_leader, kid);

          printf("Adding j: %d X(%d) %f to XL(%d) %f, leader: %d kid: %d\n",
              j, jj, X[jj], jj, XL[jj], team_leader, kid);
#endif
          //Note: may want to change to vary access pattern
          if(kid != team_leader)
          {
            //if(X[jj-brow]!=0)
            if(X(jj) !=0)
            {
#ifdef BASKER_DEBUG_NFACTOR_COL
              //printf("Atomic Adding X(%d) %f to XL(%d) %f, kid: %d\n",
              // jj, X[jj-brow], jj, XL[jj-browL], kid);
              printf("Atomic Adding X(%d) %f to XL(%d) %f, kid: %d\n",
                  jj, X[jj], jj, XL[jj], kid);
#endif

              //------------HERE-------//
              //DO DENSE in PLACE
              //Need to figure out a sparse way

              //____COMEBACK
              //Kokkos::atomic_fetch_add(&(XL[jj-browL]),
              //                     X[jj-brow]);


              //XL[jj-brow] += X[jj-brow];
              XL(jj) += X(jj);


              //X[jj-brow]=0;
              X(jj) = 0;
              //color[jj-brow] =0;

              //if(colorL[jj-browL] != 1)
              if(color[jj] != 1)
              {
                //colorL[jj-browL] = 1;
                colorL[jj] = 1;
                //patternL[p_sizeL++] = jj;
                patternL[p_sizeL++] = jj;

#ifdef BASKER_DEBUG_NFACTOR_COL
                printf("Filling: %d p_size: %d kid: %d \n", jj, p_size, kid);
#endif
              }
            }
          }
        }//end over all nnz

#ifdef BASKER_DEBUG_NFACTOR_COL
        printf("----DEBUG--- X %d %d brow: %d kid: %d\n", my_idx, blk, brow, kid);
        for(Int j = 0; j< ws_size; j++)
        {
          //printf("X[%d] = %f , pttern: %d color: %d kid: %d \n", 
          //     j+brow, X[j], pattern[j], color[j],kid);
          printf("X[%d] = %f , pttern: %d color: %d kid: %d \n", 
              j+brow, X[j], pattern[j], color[j],kid);
        }
#endif

        if(kid != team_leader)
        {
          //LL[my_idx][blk].p_size = 0;
          LL(my_idx)(blk).p_size = 0;
        }
        else
        {
#ifdef BASKER_DEBUG_NFACTOR_COL
          printf("SETTING PS: %d L:%d %d kid: %d\n",
              p_sizeL, leader_idx, blk, kid);
#endif
          //LL[leader_idx][blk].p_size = p_sizeL;
          LL(leader_idx)(blk).p_size = p_sizeL;
        }
        p_size = 0;
      }//over all blks
    }//if not team_leader

#endif

    return 0;
  }//end t_blk_col_copy_atomic()


  //local idx local blk
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_dense_copy_update_matrix
  (
   Int kid,
   Int team_leader,
   Int lvl,
   Int l, 
   Int k
  )
  {
    //printf("\n\n\n\n");
    //printf("-----------------copy_update_matrx----------");
    //printf("\n\n\n\n");

    Int       leader_idx = S(0)(kid);
    BASKER_MATRIX     &C = thread_array(kid).C;  
    Int nnz = 0;
    //COME BACK HERE

    //team_leader = find_leader(kid, l);
    //----------------Add A------------------
    //Over each blk
    Int last_blk = l+2;   
    if(lvl ==(l+1))
    {
      last_blk = LL_size(leader_idx);
    }
    // for(Int bl = l+1; bl < last_blk; bl++)
    {
      Int bl = l+1;
      Int A_col = S(lvl)(kid);
      Int A_row = (lvl==1)?(2):S(bl)(kid)%(LU_size(A_col));
      Int CM_idx = kid;

      BASKER_MATRIX  *Bp;
      if(A_row != (LU_size(A_col)-1))
      {
        //printf("upper picked, kid: %d \n", kid);
        //printf("up: %d %d kid: %d \n",
        //           A_col, A_row, kid);
        Bp = &(AVM(A_col)(A_row));
      }
      else
      {
        //printf("lower picked, kid: %d\n", kid);
        Bp = &(ALM[A_col][0]);
      }

      BASKER_MATRIX   &B  = *Bp;

      //printf("ADDING UPDATES TO B\n");
      //B.info();
      //B.print();

      team_leader       = find_leader(kid, l);
      ENTRY_1DARRAY   X = LL(leader_idx)(bl).ews;
      INT_1DARRAY    ws = LL(leader_idx)(bl).iws;
      const Int brow    = LL(leader_idx)(bl).srow;
      const Int nrow    = LL(leader_idx)(bl).nrow;
      Int p_size        = LL(leader_idx)(bl).p_size;
      const Int ws_size = LL(leader_idx)(bl).iws_size;
      Int *color        = &(ws(0));
      Int *pattern      = &(color[ws_size]);


#ifdef BASKER_DEBUG_NFACTOR_COL
      printf("copy, kid: %d bl: %d  A: %d %d \n", 
          kid, bl, A_col, A_row);
#endif


      const Int bbcol = B.scol;
      //Int
      //printf("k: %d bbcol: %d \n", k, bbcol);
      //for(Int i = B.col_ptr[k-bbcol]; 
      //  i < B.col_ptr[k-bbcol+1]; ++i)
      for(Int i = B.col_ptr(k);
          i < B.col_ptr(k+1); ++i)
      {
        Int B_row = B.row_idx(i);

        Int j = gperm(B_row+B.srow);


#ifdef BASKER_DEBUG_NFACTOR_COL
        printf("Scanning_2 A: %d %d lvl %d l: %d bl:%d brow: % d %d K: %d \n",
            B_row, j, lvl, l, bl, brow, B.srow, kid);
        // printf("Adding Aval: %f to xval: %f \n", X[B_row-brow], B.val(i));
        printf("Adding Aval: %f to xval: %f \n", X[B_row], B.val(i));
#endif


        X(B_row) += B.val(i);

        /*
           if(color[B_row] == 0)
           {
        //color[B_row-brow] = 1;
        color[B_row] = 1;
        pattern[p_size++] = B_row;
        }
        */

      }//end for over all nnz
    }//end over all blks

    //-------------move into C 
    //(Right now have to do dense but like to do sparse)

    last_blk = LL_size(leader_idx);
    //printf("------maybe l:%d lastblk: %d kid: %d\n",
    //   l, last_blk, kid);
    for(Int bl=l+1; bl<last_blk; ++bl)
    {
      Int A_col = S(lvl)(kid);
      Int A_row = (lvl==1)?(2):S(bl)(kid)%(LU_size(A_col));
      Int CM_idx = kid;
      ENTRY_1DARRAY   X   = LL(leader_idx)(bl).ews;
      INT_1DARRAY    ws   = LL(leader_idx)(bl).iws;
      const Int   ws_size = LL(leader_idx)(bl).ews_size;
      const Int      brow = LL(leader_idx)(bl).srow;
      const Int      nrow = LL(leader_idx)(bl).nrow;
      Int p_size          = LL[leader_idx][bl].p_size;

      //For recounting patterns in dense blk
      //Need better sparse update
      Int p_count  =0 ; 

      Int *color   = &(ws(0));
      Int *pattern = &(color[ws_size]);

#ifdef BASKER_DEBUG_NFACTOR_COL
      printf("moving, kid: %d  A: %d %d %d %d p_size: %d \n", 
          kid, A_col, A_row, team_leader, bl,p_size);
#endif

      //over all dim(S)
      //for(Int jj=brow; jj < (brow+nrow); jj++)
      for(Int jj=0; jj < nrow; ++jj)
      {
        //Int j = pattern[jj];
        Int j = jj;
#ifdef BASKER_DEBUG_NFACTOR_COL
        //printf("considering: %d %d %f, kid: %d\n",
        // j,brow,X[j-brow], kid);
        printf("considering: %d %d %f, kid: %d\n",
            j,brow,X[j], kid);
#endif


        if(X(j) != 0)
        {
          if(bl == l+1)
          {
#ifdef BASKER_DEBUG_NFACTOR_COL
            //printf("moving X[%d] %f kid: %d \n",
            //j, X[j-brow], kid);
            printf("moving X[%d] %f kid: %d \n",
                j, X(j), kid);
#endif
            C.row_idx(nnz) = j;
            //C.val(nnz) = X[j-brow];
            C.val(nnz) = X(j);
            nnz++;
            //X[j-brow] = 0;
            X(j) = 0;
            //color[j-brow]  = 0;
            color[j] = 0;
          }
          else
          {
#ifdef BASKER_DEBUG_NFACTOR_COL

            printf("counting [%d] %f kid: %d \n",
                j, X[j], kid);
#endif
            //pattern[p_count++] = j;
            //color[j] = 1;
          }
        }//if not empyt
      }//over all dim(S)

      if(bl == l+1)
      {

#ifdef BASKER_DEBUG_NFACTOR_COL
        printf("SETTING move_over set 0, L: %d %d kid: %d \n",
            leader_idx, bl, kid);
#endif
        //LL[leader_idx][bl].p_size = 0;
        LL(leader_idx)(bl).p_size = 0;
        p_count =0;
      }
      else
      {

        //printf("--------------SHOULD NOT BE CALLED----------\n");

#ifdef BASKER_DEBUG_NFACTOR_COL
        printf("SETTING Re-pop pattern: %d %d size: %d \n",
            leader_idx, bl, p_count);
#endif
        //LL[leader_idx][bl].p_size = p_count;
        LL(leader_idx)(bl).p_size = p_count;
      }

    }//over all blks

    //printf("kid: %d col_ptr: %d \n",
    //   kid, nnz);

    C.col_ptr(0) = 0;
    C.col_ptr(1) = nnz;

    //printf("Done with move, kid: %d found nnz: %d \n",
    //   kid, nnz);

    //C.info();
    //C.print();

    return 0;
  }//end t_dense_copy_update_matrix()


  //local idx local blk
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_copy_update_matrix
  (
   Int kid,
   Int team_leader,
   Int lvl,
   Int l, 
   Int k
  )
  {
    //printf("\n\n\n\n");
    //printf("-----------------copy_update_matrx----------");
    //printf("\n\n\n\n");

    Int       leader_idx = S(0)(kid);
    BASKER_MATRIX     &C = thread_array(kid).C;  
    Int nnz = 0;
    //COME BACK HERE

    //team_leader = find_leader(kid, l);
    //----------------Add A------------------
    //Over each blk
    Int last_blk = l+2;   
    if(lvl ==(l+1))
    {
      last_blk = LL_size(leader_idx);
    }
    // for(Int bl = l+1; bl < last_blk; bl++)
    {
      Int bl = l+1;
      Int A_col = S(lvl)(kid);
      Int A_row = (lvl==1)?(2):S(bl)(kid)%(LU_size(A_col));
      Int CM_idx = kid;

      //Bgood(remove)
      //BASKER_MATRIX_VIEW &B = AV[A_col][A_row];
      //B.init_perm(&gperm);
      //B.init_offset(k, 0); //could be made faster
      BASKER_MATRIX  *Bp;
      if(A_row != (LU_size(A_col)-1))
      {
        //printf("upper picked, kid: %d \n", kid);
        //printf("up: %d %d kid: %d \n",
        //           A_col, A_row, kid);
        Bp = &(AVM(A_col)(A_row));
      }
      else
      {
        //printf("lower picked, kid: %d\n", kid);
        Bp = &(ALM[A_col][0]);
      }

      BASKER_MATRIX   &B  = *Bp;


      //printf("ADDING UPDATES TO B\n");
      //B.info();
      //B.print();

      team_leader = find_leader(kid, l);
      //ENTRY_1DARRAY   X = LL[team_leader][bl].ews;
      ENTRY_1DARRAY   X = LL(leader_idx)(bl).ews;
      //INT_1DARRAY    ws = LL[team_leader][bl].iws;
      INT_1DARRAY    ws = LL(leader_idx)(bl).iws;
      //Int brow = LL[team_leader][bl].srow;
      //Int nrow = LL[team_leader][bl].nrow;
      const Int brow = LL(leader_idx)(bl).srow;
      const Int nrow = LL(leader_idx)(bl).nrow;
      //Int p_size = LL[team_leader][bl].p_size;
      Int p_size = LL(leader_idx)(bl).p_size;
      //Int ws_size = LL[team_leader][bl].iws_size;
      const Int ws_size = LL(leader_idx)(bl).iws_size;
      Int *color = &(ws(0));
      Int *pattern = &(color[ws_size]);


#ifdef BASKER_DEBUG_NFACTOR_COL
      printf("copy, kid: %d bl: %d  A: %d %d \n", 
          kid, bl, A_col, A_row);
#endif


      //Bgood(remove)
      //Note: note where called
      //for(Int i = B.offset; i < B.m_offset; i++)
      const Int bbcol = B.scol;
      //Int
      //printf("k: %d bbcol: %d \n", k, bbcol);
      //for(Int i = B.col_ptr[k-bbcol]; 
      //  i < B.col_ptr[k-bbcol+1]; ++i)
      for(Int i = B.col_ptr(k); i < B.col_ptr(k+1); ++i)
      {
        Int B_row = B.row_idx(i);
        Int j = gperm(B_row+B.srow);

        //Note: Do we need this anymore?
        //This is being used, might want to check out why
        //HERE!!!! (Check out G2)
        //if(B_row < brow)
        //  {
        //        
        //        printf("continue at : %d %d %d, kid: %d \n",
        //               B_row, B.srow, brow, kid);
        //        continue;
        //}

        //Note: Do we need this anymore?
        //if(B_row > (brow+nrow))
        //{
        // printf("breaking at %d %d \n", B_row, j);
        // break;
        //}

#ifdef BASKER_DEBUG_NFACTOR_COL
        printf("Scanning_2 A: %d %d lvl %d l: %d bl:%d brow: % d %d K: %d \n",
            B_row, j, lvl, l, bl, brow, B.srow, kid);
        // printf("Adding Aval: %f to xval: %f \n", X[B_row-brow], B.val(i));
        printf("Adding Aval: %f to xval: %f \n", X[B_row], B.val(i));
#endif

        //X[B_row-brow] += B.val(i);
        X(B_row) += B.val(i);
        //if(color[B_row-brow] == 0)
        if(color[B_row] == 0)
        {
          //color[B_row-brow] = 1;
          color[B_row] = 1;
          pattern[p_size++] = B_row;
        }

      }//end for over all nnz
    }//end over all blks

    //-------------move into C 
    //(Right now have to do dense but like to do sparse)


    //last_blk = LL_size[leader_idx];
    last_blk = LL_size(leader_idx);
    //printf("------maybe l:%d lastblk: %d kid: %d\n",
    //   l, last_blk, kid);
    for(Int bl=l+1; bl<last_blk; ++bl)
    {
      Int A_col = S(lvl)(kid);
      Int A_row = (lvl==1)?(2):S(bl)(kid)%(LU_size(A_col));
      Int CM_idx = kid;
      //ENTRY_1DARRAY   X = LL[team_leader][bl].ews;
      ENTRY_1DARRAY   X = LL(leader_idx)(bl).ews;
      //INT_1DARRAY    ws = LL[team_leader][bl].iws;
      INT_1DARRAY    ws = LL(leader_idx)(bl).iws;
      //Int        ws_size =LL[team_leader][bl].ews_size;
      const Int   ws_size =LL(leader_idx)(bl).ews_size;
      //Int brow = LL[team_leader][bl].srow;
      const Int brow = LL(leader_idx)(bl).srow;
      //Int nrow = LL[team_leader][bl].nrow;
      const Int nrow = LL(leader_idx)(bl).nrow;
      //Int p_size = LL[team_leader][bl].p_size;
      Int p_size = LL[leader_idx][bl].p_size;

      //For recounting patterns in dense blk
      //Need better sparse update
      Int p_count  =0 ; 

      Int *color   = &(ws(0));
      Int *pattern = &(color[ws_size]);

#ifdef BASKER_DEBUG_NFACTOR_COL
      printf("moving, kid: %d  A: %d %d %d %d p_size: %d \n", 
          kid, A_col, A_row, team_leader, bl,p_size);
#endif

      //over all dim(S)
      //for(Int jj=brow; jj < (brow+nrow); jj++)
      for(Int jj=0; jj < nrow; ++jj)
      {
        //Int j = pattern[jj];
        Int j = jj;
#ifdef BASKER_DEBUG_NFACTOR_COL
        //printf("considering: %d %d %f, kid: %d\n",
        // j,brow,X[j-brow], kid);
        printf("considering: %d %d %f, kid: %d\n",
            j,brow,X[j], kid);
#endif

        //if(X[j-brow] != 0)
        if(X(j) != 0)
        {
          if(bl == l+1)
          {
#ifdef BASKER_DEBUG_NFACTOR_COL
            //printf("moving X[%d] %f kid: %d \n",
            //j, X[j-brow], kid);
            printf("moving X[%d] %f kid: %d \n",
                j, X(j), kid);
#endif
            C.row_idx(nnz) = j;
            //C.val(nnz) = X[j-brow];
            C.val(nnz) = X(j);
            nnz++;
            //X[j-brow] = 0;
            X(j) = 0;
            //color[j-brow]  = 0;
            color[j] = 0;
          }
          else
          {
#ifdef BASKER_DEBUG_NFACTOR_COL
            //printf("counting [%d] %f kid: %d \n",
            // j, X[j-brow], kid);
            printf("counting [%d] %f kid: %d \n",
                j, X[j], kid);
#endif
            pattern[p_count++] = j;
            //color[j-brow] = 1;
            color[j] = 1;
          }
        }//if not empyt
      }//over all dim(S)

      if(bl == l+1)
      {

#ifdef BASKER_DEBUG_NFACTOR_COL
        printf("SETTING move_over set 0, L: %d %d kid: %d \n",
            leader_idx, bl, kid);
#endif
        //LL[leader_idx][bl].p_size = 0;
        LL(leader_idx)(bl).p_size = 0;
        p_count =0;
      }
      else
      {
#ifdef BASKER_DEBUG_NFACTOR_COL
        printf("SETTING Re-pop pattern: %d %d size: %d \n",
            leader_idx, bl, p_count);
#endif
        //LL[leader_idx][bl].p_size = p_count;
        LL(leader_idx)(bl).p_size = p_count;
      }

    }//over all blks

    //printf("kid: %d col_ptr: %d \n",
    //   kid, nnz);

    C.col_ptr(0) = 0;
    C.col_ptr(1) = nnz;

    //printf("Done with move, kid: %d found nnz: %d \n",
    //   kid, nnz);

    //C.info();
    //C.print();

    //set pointer 
    //Bgood(remove)
    /*
       for(Int bl = l+1; bl < l+2; bl++)
       {
       Int A_col = S[lvl][kid];
       Int A_row = (lvl==1)?(2):S[bl][kid]%(LU_size[A_col]);
       Int CM_idx = kid;
       BASKER_MATRIX_VIEW &B = AV[A_col][A_row];

       B.flip_base(&(thread_array[kid].C));
       B.k_offset = k;

       if(kid == 0)
       {
    //            printf("\n\n----------------TEST, kid: %d-----------\n\n", kid);
    // B.base->print();
    }
    }
    */

    return 0;
  }//end t_copy_update_matrix()


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
         BB.Barrier(thread_array[leader_kid].token[sublvl][function_n],
         thread_array[leader_kid].token[sublvl][1],
         size);
         */
    }
  }//end t_basker_barrier()

}//end namespace BaskerNS--functions

#undef BASKER_DEBUG_NFACTOR_COL
#endif //end ifndef backer_nfactor_col
