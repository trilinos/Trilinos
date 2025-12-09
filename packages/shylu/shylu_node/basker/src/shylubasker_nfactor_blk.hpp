// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLUBASKER_NFACTOR_BLK_HPP
#define SHYLUBASKER_NFACTOR_BLK_HPP

#include "shylubasker_matrix_decl.hpp"
#include "shylubasker_matrix_view_decl.hpp"
#include "shylubasker_matrix_view_def.hpp"
#include "shylubasker_types.hpp"
#include "shylubasker_stats.hpp"

#include <string>

#include "Kokkos_Core.hpp"
#include "Kokkos_Timer.hpp"

#include "Teuchos_ScalarTraits.hpp"

//#define BASKER_DEBUG_NFACTOR_BLK
//#define BASKER_TIMER
//#define BASKER_DETAILED_TIMER

namespace BaskerNS
{

  template <class Int, class Entry, class Exe_Space>
  struct kokkos_nfactor_domain
  {
    typedef Exe_Space                        execution_space;
    typedef Kokkos::TeamPolicy<Exe_Space>    TeamPolicy;
    typedef typename TeamPolicy::member_type TeamMember;

    Basker<Int,Entry,Exe_Space> *basker;

    kokkos_nfactor_domain()
    {}

    kokkos_nfactor_domain(Basker<Int,Entry,Exe_Space> *_basker)
    { basker = _basker;}

    BASKER_INLINE
    void operator()(const TeamMember &thread) const
    {
#if 0
      //Int kid = (Int)(thread.league_rank()*thread.team_size()+
      //              thread.team_rank());
      Int kid = basker->t_get_kid(thread);
      #ifdef BASKER_DEBUG_NFACTOR_BLK
      printf("\n-----------BLK---Kid: %d -------------\n", (int)kid);
      #endif

      // Factor kid-th leaf block
      basker->t_nfactor_blk(kid);
#else
      // Factor kid-th leaf block
      basker->t_nfactor_blk(thread);
#endif
    }//end operator
  };//end kokkos_nfactor_domain struct


  template <class Int, class Entry, class Exe_Space>
  struct kokkos_nfactor_domain_remalloc
  {
    typedef Exe_Space                        execution_space;
    typedef Kokkos::TeamPolicy<Exe_Space>    TeamPolicy;
    typedef typename TeamPolicy::member_type TeamMember;

    Basker<Int,Entry,Exe_Space> *basker;
    INT_1DARRAY                 thread_start;

    kokkos_nfactor_domain_remalloc()
    {}

    kokkos_nfactor_domain_remalloc
    (
     Basker<Int,Entry,Exe_Space> *_basker, 
     INT_1DARRAY                 _thread_start
     )
    {
      basker       = _basker;
      thread_start = _thread_start;
    }

    BASKER_INLINE
    void operator()(const TeamMember &thread) const
    {
      Int kid = thread.league_rank();
      if(thread_start(kid) != BASKER_MAX_IDX)
      {
        #ifdef BASKER_DEBUG_NFACTOR_BLK
        printf("\n-----------BLK---Kid: %d -------------\n", (int)kid);
        #endif

        basker->t_nfactor_blk(thread);
      }

    }//end operator
  };//end kokkos_nfactor_domain_remalloc struct
 

  //use local number on local blks (Crazy idea)
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_nfactor_blk(const TeamMember &thread)
  {
    using STS = Teuchos::ScalarTraits<Entry>;
    using Mag = typename STS::magnitudeType;
    const Entry zero (0.0);
    const Entry one (1.0);
    const Mag eps = STS::eps ();
    const Mag normA     = BTF_A.gnorm;
    const Mag normA_blk = BTF_A.anorm;

    Int league_rank = thread.league_rank();
    Int team_rank = thread.team_rank();
    Int team_size = thread.team_size();
    Int kid = league_rank;

    Int b = S(0)(kid); //Which blk from schedule
    BASKER_MATRIX &L   = LL(b)(0);
    BASKER_MATRIX &U   = LU(b)(LU_size(b)-1);
    BASKER_MATRIX &M   = ALM(b)(0); //A->blk
#ifdef BASKER_2DL
    //printf("Accessing blk: %d kid: %d  \n", b, kid);
    INT_1DARRAY   ws   = LL(b)(0).iws;
    ENTRY_1DARRAY X    = LL(b)(0).ews;
    Int        ws_size = LL(b)(0).iws_size;
#else  //else if BASKER_2DL
    INT_1DARRAY   ws   = thread_array(kid).iws;
    ENTRY_1DARRAY X    = thread_array(kid).ews;
    Int       ws_size  = thread_array(kid).iws_size;
#endif
    //Int          bcol  = L.scol;  //begining col //NOT UD
    Int          scol_top = btf_tabs[btf_top_tabs_offset]; // the first column index of A
    Int          brow_g   = L.srow + scol_top;   // global offset

    // integer workspaces
    Int *color    = &(ws(0));
    Int *pattern  = &(color[ws_size]);

    Int          lval  = 0;
    Int          uval  = 0;
    #ifdef BASKER_DEBUG_NFACTOR_COL
    if (Options.verbose == BASKER_TRUE) {
      {
        int  pid;       //process id
        char tcomm[256];//filename of the executable
        char state[2];  //state (R is running, S is sleeping, D is sleeping in an
                        //uninterruptible wait, Z is zombie, T is traced or stopped)
        int ppid;//          process id of the parent process
        int pgrp;//          pgrp of the process
        int sid;//           session id
        int tty_nr;//        tty the process uses
        int tty_pgrp;//      pgrp of the tty
        int flags;//         task flags
        int min_flt;//       number of minor faults
        int cmin_flt;//      number of minor faults with child's
        int maj_flt;//       number of major faults
        int cmaj_flt;//      number of major faults with child's
        int utime;//         user mode jiffies
        int stime;//         kernel mode jiffies
        int cutime;//        user mode jiffies with child's
        int cstime;//        kernel mode jiffies with child's
        int priority;//      priority level
        int nice;//          nice level
        int num_threads;//   number of threads
        int it_real_value;//  (obsolete, always 0)
        int start_time;//    time the process started after system boot
        int vsize;//         virtual memory size
        int rss;//           resident set memory size
        int rsslim;//        current limit in bytes on the rss
        int start_code;//    address above which program text can run
        int end_code;//      address below which program text can run
        int start_stack;//   address of the start of the stack
        int esp;//           current value of ESP
        int eip;//           current value of EIP
        int pending;//       bitmap of pending signals
        int blocked;//       bitmap of blocked signals
        int sigign;//        bitmap of ignored signals
        int sigcatch;//      bitmap of catched signals
        int wchan;//         address where process went to sleep
        int i0;//             (place holder)
        int i1;//             (place holder)
        int exit_signal;//   signal to send to parent thread on exit
        int task_cpu;//      which CPU the task is scheduled on
        int rt_priority;//   realtime priority
        int policy;//        scheduling policy (man sched_setscheduler)
        int blkio_ticks;//   time spent waiting for block IO
        int gtime;//         guest time of the task in jiffies
        int cgtime;//        guest time of the task children in jiffies

        FILE* f = fopen("/proc/self/stat", "r");
        if (f) {
          fscanf(f,  "%d%s%s%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d",
                 &pid, tcomm, state, &ppid, &pgrp, &sid, &tty_nr, &tty_pgrp, &flags,
                 &min_flt, &cmin_flt, &maj_flt, &cmaj_flt, &utime, &stime,  &cutime, &cstime,
                 &priority, &nice, &num_threads, &it_real_value, &start_time,  &vsize, &rss,
                 &rsslim, &start_code, &end_code, &start_stack, &esp, &eip,  &pending, &blocked,
                 &sigign, &sigcatch, &wchan, &i0, &i1, &exit_signal,  &task_cpu, &rt_priority, &policy,
                 &blkio_ticks, &gtime, &cgtime);
          printf(" pid=%d, tcomm=%s, state=%s, ppid=%d, pgrp=%d, sid=%d, tty_nr=%d, tty_pgrp=%d, flags=%d, min_flt=%d, cmin_flt=%d, maj_flt=%d, cmaj_flt%d, utime=%d, stime=%d, cutime=%d, cstime=%d, priority=%d, nice=%d, num_threads=%d, it_real_value=%d, start_time=%d, vsize=%d, rss=%d, rsslim=%d, start_code=%d, end_code=%d, start_stack=%d, esp=%d, eip=%d, pending=%d, blocked=%d, sigign=%d, sigcatch=%d, wchan=%d, i0=%d, i1=%d, exit_signal=%d, task_cpu=%d, rt_priority=%d, policy=%d, blkio_ticks=%d, gtime=%d, cgtime=%d\n",
                 pid, tcomm, state, ppid, pgrp, sid, tty_nr, tty_pgrp, flags,
                 min_flt, cmin_flt, maj_flt, cmaj_flt, utime, stime, cutime, cstime,
                 priority, nice, num_threads, it_real_value, start_time, vsize, rss,
                 rsslim, start_code, end_code, start_stack, esp, eip, pending, blocked,
                 sigign, sigcatch, wchan, i0, i1, exit_signal,  task_cpu, rt_priority, policy,
                 blkio_ticks, gtime, cgtime);
          printf("CPU : pic = %d, task_cpu = %d\n", pid, task_cpu);
        } else {
          printf( " Failed to open /proc/self/stat file\n" );
        }
      }
      printf(" thread-%d(%d): t_nfactor_blk: M(%d, 0) using LL(%d, 0) & LU(%d, %d)\n",
             kid,team_rank, b, b, b, LU_size(b)-1); fflush(stdout);
    }
    #endif

    #ifdef BASKER_DETAILED_TIMER
    Kokkos::Timer timer1;
    double time1 = 0.0, time2 = 0.0, time3 = 0.0, time4 = 0.0, time5 = 0.0, time6 = 0.0;
    #endif
    #if defined(BASKER_TIMER) | defined(BASKER_DETAILED_TIMER)
    Kokkos::Timer timer;
    #endif
    Int flops = 0;
    Int npivots = 0;

    Int newsize;
    Entry pivot (0.0);
    Mag absv (0.0);
    Mag maxv (0.0);
    Mag digv (0.0);

    Int llnnz = L.mnnz;
    Int uunnz = U.mnnz;

    Int top = ws_size;
    Int lnnz = lval;
    Int unnz = uval;
    Int xnnz;

    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("b: %d  ls: %d us: %d llnzz: %d uunzz: %d \n", 
        b, lnnz, unnz, L.nnz, U.nnz);
    printf("b: %d gperm: %d \n", b, gperm(L.srow));
    #endif

    //---TEMP--- DEBUG ---
    //ecol = 5;

    //for each column
    //for(k = scol; k < ecol; k++)
    //Note: might want to add const trick for vectorize,
    //though this loop should really not vectorize

    #ifdef BASKER_DEBUG_NFACTOR_BLK
    for(i = 0; i < M.ncol; i++)
    {
      if(L.pend(i) != BASKER_MAX_IDX)
      {
        printf("pend error: i: %d b: %d p: %d \n",
            i, b, L.pend(i));

      }
      BASKER_ASSERT(L.pend(i) == BASKER_MAX_IDX, "pend");
    }
    #endif

    if(Options.verbose == BASKER_TRUE)
    {
      printf(" thread-%ld: >  factoring_blk : b = %d, size = %ld, brow = %ld\n",
          (long)kid, (int)b, (long)M.ncol, (long)brow_g); fflush(stdout);
    }
    //#define MY_DEBUG_BASKER
    #ifdef MY_DEBUG_BASKER
    #define debug_kid 0
    {
      //printf( " t_nfactor_blk(kid = %d, %dx%d): wsize=%d\n",kid,M.nrow,M.ncol,ws_size );
      char filename[200];
      sprintf(filename,"D_%d.dat",kid);
      M.print_matrix(filename);
    }
    char filename[200];
    sprintf(filename,"t_nfactor_blk_%d_%d.dat",kid,M.ncol);
    FILE *fp = fopen(filename,"w");
    #endif

    int main_thread_id = 0;
    int worker_thread_id = (team_size == 1 ? 0 : 1);

    // reset counter used for cordinating with worker thread
    Kokkos::View<Int*, BASKER_EXE_SPACE, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::Atomic>> atomic_check(
          &(U.col_ptr(M.ncol)), 1);
    if(team_size > 1) {
      U.col_ptr(M.ncol) = 0;
      thread.team_barrier();
    }

    // factor each column
    Int maxindex = BASKER_MAX_IDX;
    for(Int k = 0; k < M.ncol; ++k)
    {
      #ifdef BASKER_DEBUG_NFACTOR_BLK
      //if (kid == debug_kid)
      {
        printf("\n---------------- K=%d/%d (%d:%d) --------------\n", 
               k,M.ncol, league_rank,team_rank); fflush(stdout);
      }
      #endif

      #ifdef BASKER_DEBUG_NFACTOR_BLK
      ASSERT(top == ws_size);
      //ASSERT entry workspace is clean
      for(i = 0 ; i < ws_size; i++)
      {
        BASKER_ASSERT(X(i) == 0, "Xerror");
      }
      //ASSERT int workspace is clean
      for(i = 0; i <  ws_size; i++)
      {
        BASKER_ASSERT(ws(i) == 0, "wserror");
      }
      #endif

      if (team_rank == main_thread_id) {

        Int lcnt = 0;
        Int ucnt = 0;

        // --------------------------------------------
        // Find sparsity/fill pattern of the current column
        #ifdef BASKER_DETAILED_TIMER
        timer1.reset();
        #endif
        //  for each nnz in column
        //  Want to change this to local blk anyway
        for(Int i = M.col_ptr(k); i < M.col_ptr(k+1); ++i)
        {
          Int j = M.row_idx(i);

          #ifdef BASKER_2D
          //Do we need this anymore ?? Don't think
          if(j >= ecol)
          {
            #ifdef BASKER_DEBUG_NFACTOR_BLK
            printf("col_break, kid: %d idx: %d \n", kid, i);
            #endif
            break;
          }
          #endif

          if (M.val(i) != zero ||
             (i+1 == M.col_ptr(k+1) && top == ws_size)) // the last element, and have not found non-zero entry, then go ahead and process this zero entry
          {
            #ifdef BASKER_2DL
            X(j) = M.val(i);
            #else
            X[j] = M.val[i];
            #endif
            #ifdef MY_DEBUG_BASKER
            if (kid == debug_kid) {
              printf( " > original: X(%d) = %e (color = %d)\n",j, X(j), color[j] );
            }
            #endif

            //NOTE:  Need a quick skip of dfs if
            //j i not pivotal (KLU)
            if(color[j] == 0)
            {
              //we want to skip the call if we can
              if(gperm(j+brow_g) != BASKER_MAX_IDX)
              {
                //printf("local_reach\n");
                t_local_reach(kid, 0, 0, j, top);
              }
              else
              {
                //printf(" %d: short (%d,%d -> %d)\n",kid, k,j, gperm(j+brow_g));
                t_local_reach_short(kid, 0, 0, j, top);
              }
            }
          }
        }//end for() each nnz in column
        xnnz = ws_size - top;
        #ifdef BASKER_DETAILED_TIMER
        time1 += timer1.seconds();
        timer1.reset();
        #endif

        #ifdef BASKER_DEBUG_NFACTOR_BLK
        if (kid == 0) {
          printf("xnnz: %d ws_size: %d top: %d \n", xnnz, ws_size, top);
        }
        #endif

        // --------------------------------------------
        // back-solve to compute the k-th column
        //t_back_solve_selective(kid, 0, 0, k, top, xnnz);
        flops += t_back_solve(kid, 0,0,  k, top, xnnz);
        #ifdef BASKER_DETAILED_TIMER
        time2 += timer1.seconds(); // time to update k-th column with previous columns
        timer1.reset();
        #endif

        #ifdef BASKER_DEBUG_NFACTOR_BLK
        if (kid == 0) {
          printf("xnnz: %d ws_size: %d top: %d (after solve)\n", xnnz, ws_size, top);
        }
        #endif

        // --------------------------------------------
        //find pivot
        pivot = zero;
        maxv = abs(zero);
        digv = abs(zero);
        maxindex = BASKER_MAX_IDX;
        Int digj = BASKER_MAX_IDX;
        for(Int i = top; i < ws_size; i++)
        {
          Int j = pattern[i];
          Int t = gperm(j+brow_g);

          #ifdef BASKER_2DL
          Entry value = X(j);
          #else
          Entry value = X[j];
          #endif

          #ifdef MY_DEBUG_BASKER
          if (kid == debug_kid)
          {
            printf("\n thread-%d: k=%d consider: j=%d t=%d value=%e maxv=%e pivot=%e\n",
                   kid,k, j+brow_g, t, value, maxv, pivot); fflush(stdout);
          }
          #endif

          absv = EntryOP::approxABS(value);

          if(t == BASKER_MAX_IDX)
          {
            lcnt++;

            if(EntryOP::gt(absv,maxv) || maxindex == BASKER_MAX_IDX) // \gt, or the first time
            {
              maxv     = absv;
              pivot    = value;
              maxindex = j;
              #ifdef MY_DEBUG_BASKER
              if (kid == debug_kid) {
                printf( " thread-%d: -> new pivot %e\n",kid, absv ); fflush(stdout);
              }
              #endif
            }

            #ifdef MY_DEBUG_BASKER
            if (kid == debug_kid) {
              printf( " thread-%d: check: gperm_array(%d + %d) = %d, germi_array(%d + %d) = %d vs %d+%d\n",kid, j,brow_g,gperm_array(j+brow_g), j,brow_g,gpermi_array(j+brow_g), k,brow_g);
              fflush(stdout);
            }
            #endif
            if (j == k)
            {
              digv = absv;
              digj = j;
              #ifdef MY_DEBUG_BASKER
              if (kid == debug_kid) {
                printf( " thread-%d -> diag %e\n",kid,absv ); fflush(stdout);
              }
              #endif
            }
          }
        }//for (i = top; i < ws_size)
        #ifdef BASKER_DETAILED_TIMER
        time3 += timer1.seconds(); // time to find pivot
        timer1.reset();
        #endif
        #ifdef MY_DEBUG_BASKER
        if (kid == debug_kid)
        {
          printf(" thread-%d > k=%d maxindex=%d pivot=%e maxv=%e, diag=%e diagj=%d tol=%e eps*normA=%e*%e=%e (nopivot=%d)\n", 
                 kid, k, maxindex, pivot, maxv, digv,digj, Options.pivot_tol,eps,normA_blk,eps*normA_blk,Options.no_pivot);
        }
        {
          const Mag eps = STS::eps ();
          const Mag normA_blk = BTF_A.anorm;
          fprintf(fp, " thread-%d > k=%d maxindex=%d pivot=%e maxv=%e, diag=%e diagj=%d tol=%e eps*normA=%e*%e=%e (nopivot=%d)\n", 
                  kid, k, maxindex, pivot, maxv, digv,digj, Options.pivot_tol,eps,normA_blk,eps*normA_blk,Options.no_pivot);
          fflush(stdout);
        }
        #endif

        //Need a BIAS towards the diagonl
        if(digj == BASKER_MAX_IDX) {
          #ifdef MY_DEBUG_BASKER // diagonal may be zero (and not in CSC) in the original matrix
          if (Options.verbose == BASKER_TRUE)
          {
            cout << "---------------------------------" << std::endl;
            cout << "  thread-" << kid << "  Failed to find diagonal" << std::endl;
            cout << "---------------------------------" << std::endl;
          }
          #endif
        } else {
          if(Options.no_pivot == BASKER_TRUE || digv > maxv * Options.pivot_tol)
          {
            #ifdef MY_DEBUG_BASKER
            if (kid == debug_kid) {
              printf(" thread-%d: using diag: %e\n",kid, X(k)); fflush(stdout);
            }
            #endif
            pivot    = X(digj);
            maxindex = digj;
          }
        }

        bool explicit_pivot = false;
        Entry lastU = zero;
        ucnt = ws_size - top - lcnt +1;
        if((maxindex == BASKER_MAX_IDX) || (pivot == zero) )
        {
          if (Options.verbose == BASKER_TRUE)
          {
            std::cout << std::endl << std::endl;
            std::cout << "---------------------------" << std::endl;
            std::cout << "  thread-" << kid 
                      << " Error: Dom Matrix(" << b << "), k = " << k
                      << " ( " << M.nrow << " x " << M.ncol << " )"
                      << " with nnz = " << M.nnz
                      << " is singular"
                      << std::endl;
            std::cout << "  norm(A)   = " << normA     << " (global)" << std::endl
                      << "  norm(A)   = " << normA_blk << " (block)"  << std::endl
                      << "  replace_tiny_pivot = " << (Options.replace_tiny_pivot ? " true " : "false" ) << std::endl
                      << "  replace_zero_pivot = " << (Options.replace_zero_pivot ? " true " : "false" ) << std::endl;
            if (Options.replace_tiny_pivot && normA_blk > abs(zero) && maxindex != BASKER_MAX_IDX) {
              std::cout << "  + replace tiny pivot with " << normA_blk * sqrt(eps) << std::endl;
            } else if (Options.replace_zero_pivot && normA_blk > abs(zero) && maxindex != BASKER_MAX_IDX) {
              std::cout << "  - replace zero pivot with " << normA_blk * eps << std::endl;
            }
            std::cout << "  Ptr       = " << M.col_ptr(k) 
                              << " " << M.col_ptr(k+1)-1 << std::endl;
            std::cout << "  Top       = " << top         << std::endl
                      << "  WS_size   = " << ws_size     << std::endl
                      << "  MaxIndex  = " << maxindex    << std::endl
                      << "  Pivot     = " << pivot       << std::endl;
            if (maxindex != BASKER_MAX_IDX) {
              std::cout << "  x(MaxInd) = " << X(maxindex) << std::endl;
            } else {
              std::cout << "  x(MaxInd) = Empty Colum"     << std::endl;
            }
            std::cout << " Lcount    = " << lcnt       << std::endl
                      << "---------------------------" << std::endl;
            /*if (kid == 0)
            {
              M.print_matrix("D.dat");
            }*/
          }

          if (Options.replace_tiny_pivot && normA_blk > abs(zero) && maxindex != BASKER_MAX_IDX) {
            pivot = normA_blk * sqrt(eps);
            X(maxindex) = pivot;
            npivots ++;
          } else if (Options.replace_zero_pivot && normA_blk > abs(zero) && maxindex != BASKER_MAX_IDX) {
            pivot = normA_blk * eps;
            X(maxindex) = pivot;
            npivots ++;
          } else {
            // replace-tiny-pivot not requested, or the current column is structurally empty after elimination
            if (Options.replace_tiny_pivot && normA_blk > abs(zero)) {
              // just insert tiny pivot on diagonal
              maxindex = k;
              while (gperm(maxindex+brow_g) != BASKER_MAX_IDX && maxindex < M.ncol) {
                maxindex ++;
              }
              if (maxindex < M.ncol) {
                if (Options.verbose == BASKER_TRUE)
                {
                  std::cout << "  thread-" << kid << " Explicit tiny pivot for maxind = " << maxindex << std::endl;
                }
                pivot = normA_blk * sqrt(eps);
                lastU = pivot;
                npivots ++;
                explicit_pivot = true;
              }
            } else if (Options.replace_zero_pivot && normA_blk > abs(zero)) {
              // just insert tiny pivot on diagonal
              maxindex = k;
              while (gperm(maxindex+brow_g) != BASKER_MAX_IDX && maxindex < M.ncol-1) {
                maxindex ++;
              }
              if (maxindex < M.ncol) {
                if (Options.verbose == BASKER_TRUE)
                {
                  std::cout << "  thread-" << kid << " Explicit nonzero pivot for maxind = " << maxindex
                            << "(" << gperm(maxindex+brow_g) << ")" << std::endl;
                }
                pivot = normA_blk * eps;
                lastU = pivot;
                npivots ++;
                explicit_pivot = true;
              }
            }
            if (!explicit_pivot) {
              thread_array(kid).error_type =
                BASKER_ERROR_SINGULAR;
              thread_array(kid).error_blk    = b;
              thread_array(kid).error_subblk = 0; 
              thread_array(kid).error_info   = k;
              // let worker thread know something went wrong
              atomic_check(0) = k+1;
              return BASKER_ERROR;
            }
          }
        } else if (Options.replace_tiny_pivot && normA_blk > abs(zero) && abs(pivot) < normA_blk * sqrt(eps)) {
          if (Options.verbose == BASKER_TRUE)
          {
            std::cout << std::endl << std::endl;
            std::cout << "---------------------------" << std::endl;
            std::cout << "  thread-" << kid 
                      << " Dom Matrix(" << b << "), k = " << k
                      << " ( " << M.nrow << " x " << M.ncol << " )"
                      << " with nnz = " << M.nnz
                      << " : replace tiny pivot( " << pivot << " -> "
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
          npivots ++;
        }

        // store pivot
        gperm(maxindex+brow_g) = k+brow_g;
        gpermi(k+brow_g) = maxindex + brow_g;
        #ifdef MY_DEBUG_BASKER
        if (kid == debug_kid)
        {
          printf(" + %d: gperm(%d + %d) = %d\n",kid,maxindex,brow_g,k+brow_g );
        }
        {
          fprintf(fp, " + %d: gperm(%d + %d) = %d\n",kid,maxindex,brow_g,k+brow_g );
        }
        #endif

        #ifdef MY_DEBUG_BASKER
        if (kid == debug_kid) {
          for (int ii=brow_g; ii<brow_g+M.ncol; ii++) {
            if (gperm_array(ii) != ii) printf( " thread-%d: gperm_array(%d) = %d\n",kid,ii,gperm_array(ii) );
          }
          printf( "\n" );
          fflush(stdout);
        }
        #endif

        #ifdef BASKER_DEBUG_NFACTOR
        //if(maxindex != k)
        //  {
        //    cout << "Permuting Pivot: " << k << " as row " 
        //         << maxindex << endl;
        //  }
        #endif

        //Note: Come back to this!!!!
        if(lnnz + lcnt > llnnz)
        {
          newsize = lnnz * 1.1 + 2 *M.nrow + 1;

          if (Options.verbose == BASKER_TRUE)
          {
            printf("kid = %d, b = %ld: Reallocing L oldsize: %ld current: %ld count: %ld (new size = %d*1.1 = 2*%d =  %ld)\n",
                   (int)kid, (long)b, (long)llnnz, (long)lnnz, (long)lcnt, (int)lnnz, (int)M.nrow, (long)newsize);
          }

          thread_array(kid).error_blk = b;
          thread_array(kid).error_subblk = 0;
          if(Options.realloc == BASKER_FALSE)
          {
            thread_array(kid).error_type = BASKER_ERROR_NOMALLOC;
            // let worker thread know something went wrong
            atomic_check(0) = k+1;
            return BASKER_ERROR;
          }
          else
          {
            thread_array(kid).error_type = BASKER_ERROR_REMALLOC;
            thread_array(kid).error_info = newsize;
            // let worker thread know something went wrong
            atomic_check(0) = k+1;
            return BASKER_ERROR;
          }
        }
        if(unnz+ucnt > uunnz)
        {
          newsize = uunnz*1.1 + 2*M.nrow+1;

          if (Options.verbose == BASKER_TRUE)
          {
            printf("kid = %d, b = %ld: Reallocing U oldsize: %ld newsize: %ld  k: %ld (new size = %d*1.1 + 2*%d = %d)\n",
                   (int)kid, (long)b, (long)uunnz, (long)unnz+ucnt, (long)k, (int)uunnz, (int)M.nrow, (int)newsize);
          }

          thread_array(kid).error_blk = b;
          thread_array(kid).error_subblk = -1;
          if(Options.realloc == BASKER_FALSE)
          {
            thread_array(kid).error_type = BASKER_ERROR_NOMALLOC;
            // let worker thread know something went wrong
            atomic_check(0) = k+1;
            return BASKER_ERROR;
          }
          else
          {
            thread_array(kid).error_type = BASKER_ERROR_REMALLOC;
            thread_array(kid).error_info = newsize;
            // let worker thread know something went wrong
            atomic_check(0) = k+1;
            return BASKER_ERROR;
          }
        }

        L.row_idx(lnnz) = maxindex;
        L.val(lnnz)     = one;
        #ifdef MY_DEBUG_BASKER
        if (kid == debug_kid) {
          printf(" thread-%d: > L(%d,%d): %e \n",
                  kid, maxindex, k, L.val(lnnz));
          fflush(stdout);
        }
        #endif
        lnnz++;

        #ifdef MY_DEBUG_BASKER
        if (kid == debug_kid) {
          printf(" > for: ws_size: %d top: %d \n", ws_size, top);
        }
        #endif
        #ifdef BASKER_DETAILED_TIMER
        timer1.reset();
        #endif
        for(Int i = top; i < ws_size; i++)
        {
          Int j = pattern[i];
          Int t = gperm(j+brow_g);

          #ifdef MY_DEBUG_BASKER
          if (kid == debug_kid) {
            printf("> insert(pattern[%d] => j = %d): t = %d x = %e with diag = %d \n", i, j+brow_g, t, X(j), k+brow_g);
          }
          #endif            

          //Note can not exclude numeric cancel for prune
          //if fill-in
          #ifdef BASKER_2DL
          //if(X(j) != 0)
          #else
          //if(X[j] != 0)
          #endif
          {
            if(t != BASKER_MAX_IDX)
            {
              if(t < (k+brow_g))
              {
                U.row_idx(unnz) = t-brow_g;
                #ifdef BASKER_2DL
                U.val(unnz) = X(j);
                #else
                U.val[unnz] = X[j];
                #endif

                #ifdef MY_DEBUG_BASKER
                if (kid == debug_kid) {
                  printf(" thread-%d: U(%d,%d): %e \n",
                         kid, t-brow_g, k, X(j));
                  fflush(stdout);
                }
                #endif
                unnz++;
              }
              else
              {
                #ifdef BASKER_2DL
                lastU = X(j);
                #else
                lastU = X[j];
                #endif
                #ifdef MY_DEBUG_BASKER
                if (kid == debug_kid) {
                  printf(" thread-%d: lastU(%d): %e \n",
                         kid, k, X(j));
                  fflush(stdout);
                }
                #endif
              }
            }
            else if (t == BASKER_MAX_IDX)
            {
              L.row_idx(lnnz) = j;
              #ifdef BASKER_2DL
              L.val(lnnz) = EntryOP::divide(X(j),pivot);
              #else
              L.val(lnnz) = EntryOP::divde(X(j),pivot);
              #endif

              #ifdef MY_DEBUG_BASKER
              if (kid == debug_kid) {
                printf(" thread-%d: L(%d,%d): %e/%e = %e \n",
                       kid, j, k, X(j),pivot, L.val(lnnz));
                fflush(stdout);
              }
              #endif
              lnnz++;
            }
          }//end if() not 0

          //Note: move x[j] inside of if() not 0....extra ops this way
          #ifdef BASKER_DEBUG_NFACTOR_BLK
          printf("Zeroing element: %d \n", j);
          #endif
#ifdef SHYLU_BASKER_WORKER_THREAD_VERSION_1
          // X(maxindex) store pivot val,
          // which is also stored as the last entry of U
          if (j != maxindex)
#endif
          {
            #ifdef BASKER_2DL
            X(j) = zero;
            #else
            X[j] = zero;
            #endif
          }
        }//end if(x[i] != 0)
        #ifdef BASKER_DETAILED_TIMER
        time4 += timer1.seconds(); // time to scale with pivot
        #endif

        //Fill in last element of U
        U.row_idx(unnz) = k;
        U.val(unnz) = lastU;
        #ifdef MY_DEBUG_BASKER
        //if (kid == debug_kid)
        {
          printf(" thread-%d: > last U(%d,%d): %e %s, %d\n",
                  kid, k, k, U.val(unnz), (digv == EntryOP::approxABS(lastU) ? "(dig)" : 
                                          (maxv == EntryOP::approxABS(lastU) ? "(piv)" : "(warn)")),unnz);
          fflush(stdout);
        }
        #endif
        unnz++;

        xnnz = 0;
        top = ws_size;
        L.col_ptr(k+1) = lnnz;

        if (team_size > 1) {
          if (k+1 < M.ncol) {
            U.col_ptr(k+1) = unnz; // set U.col_ptr(k+1) before setting check and lettng worker thread use this column
            atomic_check(0) = k+1; // counter to cordinate with worker threads
            //printf( " > %d:%d: U.col_ptr(%d) = %d, U.col_ptr(%d) = %d\n",league_rank,team_rank,M.ncol,U.col_ptr(M.ncol), k+1,U.col_ptr(k+1) ); fflush(stdout);
          } else {
            // The last column (just set col_ptr(M.ncol) = unnz)
            atomic_check(0) = unnz; // counter to cordinate with worker threads
            //printf( " > %d:%d: U.col_ptr(%d) = %d\n",league_rank,team_rank, k+1,U.col_ptr(k+1) ); fflush(stdout);
          }
        } else {
          U.col_ptr(k+1) = unnz;
        }

        #ifdef MY_DEBUG_BASKER
        if (kid == debug_kid) {
          BASKER_MATRIX &U = LU(b)(LU_size(b)-1);
          printf("L.col_ptr(%d) = %d:%d with L = LL(%d)(0)\n",  (int)k, (int)L.col_ptr(k), (int)L.col_ptr(k+1), (int)b);
          printf("U.col_ptr(%d) = %d:%d with U = LU(%d)(%d)\n", (int)k, (int)U.col_ptr(k), (int)U.col_ptr(k+1), (int)b, (int)LU_size(b)-1);
        }
        #endif

        t_prune(kid, 0, 0, k, maxindex);
      }

#ifdef SHYLU_BASKER_WORKER_THREAD_VERSION_1
      // !! Synch threads !!
      thread.team_barrier();
      if (team_rank == worker_thread_id) {
        // worker thread read pivot, and set X(maxindex) to zero
        maxindex = gpermi(k+brow_g) - brow_g;
        pivot = X(maxindex);
        X(maxindex) = zero;
      }
      // Make sure X(maxindex) is zeroed out
      thread.team_barrier();
#endif


      if (team_rank == worker_thread_id)
      {
        #ifdef BASKER_2DL
        //-----------------------Update offdiag-------------//
        #ifdef BASKER_DETAILED_TIMER
        timer1.reset();
        #endif

#ifndef SHYLU_BASKER_WORKER_THREAD_VERSION_1
        if (team_size > 1) {
          // busy-wait for main thread to finish factoring the column
          Int check_k = atomic_check(0);
          while (check_k <= k) {
            check_k = atomic_check(0); //U.col_ptr(M.ncol);
          }
          // worker thread read pivot, and set X(maxindex) to zero
          pivot = U.val(U.col_ptr(k+1)-1);
          if (thread_array(kid).error_type != BASKER_ERROR_NOERROR) {
            return BASKER_ERROR;
          }
        }
        //printf( " -> %d:%d: pivot = %e (col_ptr(%d)=%d, %d)\n",league_rank,team_rank, pivot, M.ncol,U.col_ptr(M.ncol),k); fflush(stdout);
#endif

        //if (kid == 0) printf( "\n blk_row = %d:%d\n",1,LL_size(b)-1 );
        // for each row of remainng L, apply the inverse with the k-th column of U
        for(Int blk_row = 1; blk_row < LL_size(b); ++blk_row)
        {
          //Do back solve of off-diag blocks
          #ifndef BASKER_INC_LVL
          t_back_solve_offdiag(kid,
              b, blk_row,
              b, blk_row,
              //Note, different meaning
              k,
              U.val, U.row_idx,
              U.col_ptr(k+1)-U.col_ptr(k),
              U.col_ptr(k),
              BASKER_TRUE);
          #endif
          #ifdef BASKER_DETAILED_TIMER
          time5 += timer1.seconds();
          timer1.reset();
          #endif

          //Move these factors into Local Ls
          Int move_error = 
            t_move_offdiag_L(kid,
                b, blk_row,
                b, blk_row,
                k, pivot);

          if(move_error == BASKER_ERROR)
          {
            return BASKER_ERROR;
          }
          #ifdef BASKER_DETAILED_TIMER
          time6 += timer1.seconds();
          timer1.reset();
          #endif
        }//end over all diag
        #endif
      }

      #if 0//def BASKER_DETAILED_TIMER
      if (k%1000 == 0) {
        double time_facto_k = timer.seconds();
        printf( " > %d: Time : %lf (%lf %lf %lf %lf %lf), %d %d, %d %d\n",k,time_facto_k,time1,time2,time3,time4,time5,lnnz,unnz,npivots,flops );
      }
      #endif
    }//end for() over all columns
    #ifdef BASKER_DETAILED_TIMER
    {
      double time_facto_k = timer.seconds();
      printf( " > %d:%d: Time : %lf (%lf %lf %lf %lf %lf %lf), %d %d, %d %d\n", int(kid),int(team_rank), time_facto_k,time1,time2,time3,time4,time5,time6, int(lnnz),int(unnz),int(npivots),int(flops) );
    }
    #endif

    if (team_rank == main_thread_id) {
      L.nnz = lnnz;
      U.nnz = unnz;
    }
    if(Options.verbose == BASKER_TRUE)
    {
      printf(" thread-%ld: >  factoring_blk : nnzL = %ld, nnzU = %ld (%ld x %ld)\n",
          (long)kid, (long)lnnz, (long)unnz, (long)M.nrow, (long)M.ncol ); fflush(stdout);
    }

    #ifdef BASKER_TIMER
    double time_facto = timer.seconds();
    printf("Time Dom Facto(%d): %lf, n = %d, nnz(L) = %d, nnz(U) = %d \n", (int)kid, time_facto,
           (int)L.ncol, (int)L.col_ptr(L.ncol), (int)U.col_ptr(U.ncol)); fflush(stdout);
    #endif

    #ifdef MY_DEBUG_BASKER
    fclose(fp);
    {
      char filename[200];
      sprintf(filename,"P_%d.dat",kid);
      FILE *fp = fopen(filename, "w");
      for (int j = 0; j < L.ncol; j++) {
        fprintf(fp,"%d %d\n",gperm(j+brow_g)-brow_g,gpermi(j+brow_g)-brow_g);
      }
      fclose(fp);

      // D(P(:,2), :) = L(P(:,2), :) * U
      sprintf(filename,"L_%d.dat",kid);
      L.print_matrix(filename);

      sprintf(filename,"U_%d.dat",kid);
      U.print_matrix(filename);
    }
    #endif

    return 0;
  }//end t_nfactor_blk()


 
  template<class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  void Basker<Int,Entry,Exe_Space>::t_local_reach_short
  (
   const Int kid, 
   const Int lvl,
   const Int l,
   const Int j,
   Int &top
   )
  {
    //Setup variables
    const Int wsb    = S(0)(kid);

    INT_1DARRAY  ws   = LL(wsb)(l).iws;
    const Int ws_size = LL(wsb)(l).iws_size;

    Int *color       = &(ws(0));
    Int *pattern     = &(ws(ws_size));

    color[j]       = 2;
    pattern[--top] = j;
    #ifdef MY_DEBUG_BASKER
    if (kid == debug_kid) {
      printf( " > t_local_reach_short(top = %d)\n",(int)top );
    }
    #endif
  }//end t_local_reach short()
  

  template <class Int, class Entry, class Exe_Space>
  inline
  void Basker<Int,Entry,Exe_Space>::t_prune
  (
   const Int kid,
   const Int lvl, 
   const Int l, 
   const Int k, 
   const Int pivotrow 
   )
  {
    const Int  scol_top = btf_tabs[btf_top_tabs_offset]; // the first column index of A
    const Int  b      = S(lvl)(kid);

    //const Int wsb    = S(0)(kid);
    BASKER_MATRIX &L = LL(b)(0);
    const Int U_col  = S(lvl)(kid);
    Int U_row        = LU_size(U_col)-1;
    BASKER_MATRIX &U = LU(U_col)(U_row);

    //if (kid == 0 && k == 0) {
    //  printf(" %d: t_prune,k = %d  L = LL(%d)(%d) U = LU(%d)(%d) pivotrow = %d\n", 
    //         kid, k,b,0, U_col, U_row, pivotrow);
    //}
    //Scan U and check if we any columns we can prune :)
    //don't include the last entry
    for(Int ui = U.col_ptr(k); ui < U.col_ptr(k+1)-1; ++ui)
    {
      Int j = U.row_idx(ui);
      //printf("prune, j = %d k = %d \n", j, k);

      //if(j>=k)
      {
        //printf("Error: j: %d k: %d U: %d %d ui: %d %d \n",
        //j, k, U_col, U_row, U.col_ptr(k), U.col_ptr(k+1));
        BASKER_ASSERT(j < k, "Pruning, j not less than k");
      }

      //Check if this column has not been pruned
      if(L.pend(j)==BASKER_MAX_IDX)
      {
        //printf("This column can be pruned\n");

        for(Int li = L.col_ptr(j); li < L.col_ptr(j+1); ++li)
        {
          //printf("considering row: %d %d \n",
          //     L.row_idx(li),
          //     L.row_idx(li)+L.srow);

          //if we can prune ?
          if(pivotrow == L.row_idx(li))
          {
            //printf("This can be pruned\n");

            //order L 
            //partition those [above pivotrow, below]
            Int phead = L.col_ptr(j);
            Int ptail = L.col_ptr(j+1);

            //printf("phead: %d ptail: %d \n",
            //           phead, ptail);

            //partion and bubble up
            while(phead < ptail)
            {
              Int i = L.row_idx(phead);
              if(gperm(i+L.srow+scol_top) >= 0)
              {
                //advance top ptr
                //printf("adv head(%d) %d %d \n",
                //           phead+1,
                //           i+L.srow,
                //           gperm(i+L.srow));
                phead++;
              }
              else
              {
                //flip head and tail
                //printf("flipping head:%d tail:%d \n",
                //           phead, ptail);
                //printf("head: %d %d %f \n",
                //           phead, 
                //           L.row_idx(phead)+L.srow,
                //           L.val(phead));
                //printf("tail: %d %d %f \n",
                //           ptail-1,
                //           L.row_idx(ptail-1)+L.srow,
                //           L.val(ptail-1));

                ptail--;
                L.row_idx(phead) = L.row_idx(ptail);
                L.row_idx(ptail) = i;
                Entry x = L.val(phead);
                L.val(phead) = L.val(ptail);
                L.val(ptail) = x;
              }//End Flip
            }//end -- while(for flip)

            //printf("pruned, L.pend(%d) = %d\n",j, ptail);
            L.pend(j) = ptail;
            break;
          }//if pivotrow == L.row_idx(..)
        }//end scan over all Li
      }//if L.pend == enmpty
    }//over all Ui
  }//end t_prune()
 

  //updated local reach based
  //Note that this function needs to be tight
  template <class Int, class Entry, class Exe_Space>
  inline
  void Basker<Int,Entry,Exe_Space>::t_local_reach
  (
   const Int kid, 
   const Int lvl,
   const Int l,
   Int j, 
   Int &top
  )
  {
    
    //Setup variables
    const Int  b        = S(lvl)(kid);
    const Int  wsb      = S(0)(kid);
    BASKER_MATRIX  &L   = LL(b)(0);
    const Int  scol_top = btf_tabs[btf_top_tabs_offset]; // the first column index of A
    const Int  brow_g   = L.srow + scol_top;   // global offset

    INT_1DARRAY    ws  = LL(wsb)(l).iws;
    const Int  ws_size = LL(wsb)(l).iws_size;
   
    Int *pattern     = &(ws(ws_size));
    Int *stack       = &(pattern[ws_size]);
    Int *store       = &(stack[ws_size]);
    

    Int i, t, head, i1;
    Int start, end;

    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("local_reach, L: %d %d  X: %d %d j: %d, kid: %d \n",
           b, 0, wsb, l, j, kid);
    #endif

    start    = -1;
    head     = 0;
    stack[0] = j;

    while(head != BASKER_MAX_IDX)
    { 
      #ifdef BASKER_DEBUG_LOCAL_REACH
      if (kid == 0) {
        printf(" > head: %d \n", head);
      }
      #endif

      j = stack[head];
      t = gperm(j+brow_g);
        
      #ifdef BASKER_DEBUG_LOCAL_REACH
      printf("--------DFS: %d %d %d -------------\n", j, j+brow, t);
      #endif

      if(ws(j) == 0)
      {
          ws(j) = 1;

          if(t!=BASKER_MAX_IDX)
          {
            #ifdef BASKER_DEBUG_LOCAL_REACH
            printf("reach.... j: %d t:%d L.scol %d \n",
                j, t, L.scol);
            printf("colend: %d pend(%d): %d %d \n",
                L.col_ptr(t+1-L.scol),
                t-L.scol,
                L.pend(j),
                L.pend(t-L.scol));
            #endif
            start = 
              (L.pend(t-brow_g) == BASKER_MAX_IDX) ? L.col_ptr(t+1-brow_g) : L.pend(t-brow_g); 
          }
      }
      else
      {
        start = store[j];
      }

      //We want to go backwards through this
      //This little insight can save time
      end = L.col_ptr(t-brow_g);
      #ifdef BASKER_DEBUG_NFACTOR_BLK
      if (kid == 0) {
        printf("t: %d start: %d end: %d \n",t, start, end);
      }
      #endif
      for(i1 = --start; i1 >= end; --i1)
      {
        i = L.row_idx(i1);

        //printf("Search i1: %d  i: %d %d %d \n",
        //           i1, i, i+L.scol, gperm(i+L.scol));

        if(ws(i) != 0)
        {
            //printf("continue called\n");
            continue;
        }
        else
        {
          #ifdef BASKER_DEBUG_NFACTOR_BLK
          if (kid == 0) {
            printf( " gperm(%d + %d) = %d\n",i,brow_g,gperm(i+brow_g) );
          }
          #endif
          if(gperm(i+brow_g) >= 0)
          {
            store[j] = i1;
            stack[++head] = i;
            break;
          }
          else
          {
            ws(i) = 2;
            pattern[--top] = i;
            //printf("Adding idx: %d %d  to pattern at location: %d \n",i, i+L.scol, top);
          }
        }
      }// end for i1

      //if we have used up the whole column
      #ifdef BASKER_DEBUG_NFACTOR_BLK
      if (kid == 0) {
        printf( " i1 = %d, end = %d\n",i1,end );
      }
      #endif
      if(i1 < end)
      {
        head--;
        ws(j) = 2;
        pattern[--top] = j;
      }
      //printf("head: %d \n", head);
      /*
        if(head == 0) {head = BASKER_MAX_IDX;}
        else {--head;}
      */
      //printf("Adding idx: %d to pattern at location: %d \n",j, *top);
      //printf("head2: %d \n", head);
    }//end while 

  }//end t_local_reach (new)

   //Note: need to fix indexing for color
  //Uses local idx for local blks
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_local_reach_old
  (Int kid, Int lvl, Int l,
   Int j, Int *top)
  {

    //Setup variables
    const Int      b   = S(lvl)(kid);
    const Int     wsb  = S(0)(kid);
    BASKER_MATRIX  &L  = LL(b)(0);
    #ifdef BASKER_2DL
    INT_1DARRAY    ws  = LL(wsb)(l).iws;
    const Int  ws_size = LL(wsb)(l).iws_size;
    #else
    INT_1DARRAY    ws  = thread_array(kid).iws;
    Int        ws_size = thread_array(kid).iws_size;
    #endif

    const Int  scol_top = btf_tabs[btf_top_tabs_offset]; // the first column index of A
    const Int  brow_g   = L.srow + scol_top;   // global offset

    Int *color       = &(ws[0]);
    Int *pattern     = &(ws[ws_size]);
    Int *stack       = &(pattern[ws_size]);
    Int *store       = &(stack[ws_size]);
    
    /*
    auto color = Kokkos::subview(ws, 
                                 std::make_pair((Int)0, (Int)ws_size+1));
    auto pattern = Kokkos::subview(ws,
                                   std::make_pair((Int)ws_size, (Int)(2*ws_size)-1));
    auto stack   = Kokkos::subview(ws,
                             std::make_pair((Int)2*ws_size, (Int)(3*ws_size)-1));
    auto store   = Kokkos::subview(ws,
                                 std::make_pair((Int)3*ws_size, (Int)(4*ws_size)-1));
    */

    Int i, t, head, i1;
    Int start, end, done;
    
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("local_reach, L: %d %d  X: %d %d, kid: %d \n", b, 0, wsb, l, kid);
    #endif

    start = -1;
    head = 0;
    stack[head] = j;

    //while(head < L.max_idx)
    while(head != BASKER_MAX_IDX)
    { 
      #ifdef BASKER_DEBUG_LOCAL_REACH
      printf("stack_offset: %d head: %d \n", stack_offset , head);
      ASSERT(head > -1);
      #endif

      j = stack[head];
      t = gperm(j+brow_g);

      #ifdef BASKER_DEBUG_LOCAL_REACH
      printf("----------DFS: %d %d -------------\n", j, t);
      #endif

      if(color[j] == 0)
      {            
        color[j] = 1;

        //if not permuted, might be able to make this smaller now
        //if((t!=L.max_idx) && (t>=L.scol) && (t<(L.scol+L.ncol)))
        //if(t!=L.max_idx)
        if(t!=BASKER_MAX_IDX)
        {
          #ifdef BASKER_DEBUG_LOCAL_REACH
          printf("reach.... j: %d t:%d L.scol %d \n", j, t, L.scol);
          #endif
          // continue from top of jth col
          start = L.col_ptr(t-L.scol);
        }
        else
        {
          #ifdef BASKER_DEBUG_LOCAL_REACH
          printf("L.scol: %d  L.ncol: %d t: %d \n", L.scol, L.ncol,t);
          #endif
        }
      }
      else
      {
        // continue from left-out
        start = store[j];
      }
      done = 1;
  
      //if  permuted
      //if(t != L.max_idx), might be able to make this smaller now
      //if((t!=L.max_idx) && (t>=L.scol) && (t<(L.scol+L.ncol)))
      //if(t!=L.max_idx)
      if(t!=BASKER_MAX_IDX)
      {
        // shift into A
        t -= scol_top;

        //We want to go through this backward
        //That way our first search will end fast
        //our later searchs will use our first searchs
        //end = L.col_ptr[t+1-L.scol];
        end = L.col_ptr(t+1-L.scol);
        for(i1 = start; i1 < end; i1++)
        {
          i = L.row_idx(i1);
          if(color[i] == 0)
          {
            // ith row has not been processed
            //  -> push into stack and break
            #ifdef BASKER_INC_LVL
            inc_lvl++;
            #endif
            head++;
            stack[head] = i;
            store[j] = i1+1;
            done = 0;
            break;
          }
        }
      }
      if(done)
      {
        pattern[--*top] = j;
        color[j] = 2;
        if(head == 0) {
          head = BASKER_MAX_IDX;
        }
        else
        {
          head--;
        }
        //printf("Adding idx: %d to pattern at location: %d \n",j, *top);
      }
    }//end while 
    return 0;
  }//end t_local_reach()



  //local idx for local blks
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry,Exe_Space>::t_back_solve
  (
   Int kid, 
   Int lvl,
   Int l,
   Int k, Int top,
   Int xnnz)
  {
    const Int      b = S(lvl)(kid);
    const Int    wsb = S(0)(kid);
    BASKER_MATRIX &L = LL(b)(0);
    #ifdef BASKER_2DL
    INT_1DARRAY   ws = LL(wsb)(l).iws;
    ENTRY_1DARRAY X  = LL(wsb)(l).ews;
    Int      ws_size = LL(wsb)(l).iws_size;
    #else
    INT_1DARRAY   ws = thread_array(kid).iws;
    ENTRY_1DARRAY  X = thread_array(kid).ews;
    Int      ws_size = thread_array(kid).iws_size;
    #endif
    
    const Entry zero (0.0);
    Int *color   = &(ws[0]);
    Int *pattern = &(color[ws_size]);
    
    const Int  scol_top = btf_tabs[btf_top_tabs_offset]; // the first column index of A
    const Int  brow_g   = L.srow + scol_top;   // global offset
    const Int  local_offset = L.scol; //L.srow

    Int flops = 0;
    Int top1 = top;
    for(Int pp = 0; pp < xnnz; ++pp)
    {
      Int j = pattern[top1];
      Int t = gperm(j+brow_g);
      if(t != BASKER_MAX_IDX)
      {
        const Entry xj = X(j);
        if (xj != zero) {
          //We could use kokkos local memory,
          //to get rid of this pragma
          //
          //What we want to think about is how to 
          //do this in a dense way
          //#pragma ivdep (may be slower)
          //
          // shift into A
          const Int col = (t-scol_top)-local_offset;
          for(Int p = L.col_ptr(col)+1; p < L.col_ptr(col+1); ++p)
          {
            X(L.row_idx(p)) -= L.val(p)*xj;
            flops += 2;
          }//end for() over each nnz in the column
        }
      }//end if() not permuted
      color[j] = 0;
      top1++;
    }//end for() over all nnz in LHS

    return flops;
  }//end t_back_solve()

  //uses local idx for local blks
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_dense_move_offdiag_L
  (Int kid,
   Int blkcol, Int blkrow,
   Int X_col, Int X_row,
   Int k, Entry pivot)
  {
    BASKER_MATRIX &L    = LL(blkcol)(blkrow);
   
    INT_1DARRAY   ws    = LL(X_col)(X_row).iws;
    ENTRY_1DARRAY X     = LL(X_col)(X_row).ews;
   

    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("t_dense_move_offdiag_L, kid=%d, k=%d:  L (%d %d) X (%d %d), nnz=%d\n",
           kid, k, blkcol,blkrow, X_col, X_row, L.nnz);
    #endif

   
    Int *color   = &(ws(0));
    Int    lnnz  = L.col_ptr(k);
   
    /*
    if((p_size) > (llnnz-lnnz))
      {
        printf("-Warning, Need to remalloc L: %d %d kid: %d current size: %d used_size: %d  addition: %d \n",
               blkcol, blkrow, kid, llnnz,lnnz,p_size  );
        
      }
    */

    for(Int j = 0; j < L.nrow; ++j)
    {
      if(X(j) != (Entry)(0) )
      {
        //Int t = gperm(j+brow);
        if (lnnz >= L.nnz) { // this should not happen since allocated as dense separator blocks
          if (Options.verbose == BASKER_TRUE)
          {
            printf("Move Off-diag L failed with insufficient storage L(%d,%d).nnz = %d\n",
                   (int)blkcol, (int)blkrow, (int)L.nnz );
          }
          BASKER_ASSERT(true, "\n Not enough memory allocated for off-diagonal L\n");
          return BASKER_ERROR;
        }
#ifdef BASKER_DEBUG_NFACTOR_BLK
        printf("L-Moving, kid: %d j: %d val: %f lnnz: %d \n",
            kid, j, X[j]/pivot, lnnz);
#endif

        color[j] = 0;

        L.row_idx(lnnz) = j;

        //L.val(lnnz) = X(j)/pivot;
        L.val(lnnz) = EntryOP::divide(X(j),pivot);
        #ifdef BASKER_DEBUG_NFACTOR_BLK
        if (blkcol == 2 && blkrow == 1) printf( " L.val(%d) = (%e %d) for %d)\n",lnnz,L.val(lnnz),j, k );
        #endif
        X(j) = 0;

        //Note come back to fix
#ifdef BASKER_INC_LVL
        L.inc_lvl[lnnz] = INC_LVL_TEMP[j];
#endif
        lnnz++;
      }
    }
    //printf("L-Moving, kid: %d col_ptr(%d): %d \n",
    //           kid, k-bcol+1, lnnz);
    
    //L.col_ptr[k-bcol+1] = lnnz;
    L.col_ptr(k+1) = lnnz;
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    if (blkcol == 2 && blkrow == 1) printf( " L.colptr(%d) = %d\n",k+1,lnnz );
    #endif

    //LL(X_col)(X_row).p_size = 0;
    LL(X_col)(X_row).p_size = 0;

    return 0;
  }//end t_dense_offdiag_mov_L()


  //uses local idx for local blks
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::t_move_offdiag_L
  (Int kid,
   Int blkcol, Int blkrow,
   Int X_col, Int X_row,
   Int k, Entry pivot)
  {
    BASKER_MATRIX &L    = LL(blkcol)(blkrow);
   
    INT_1DARRAY   ws    = LL(X_col)(X_row).iws;
    ENTRY_1DARRAY X     = LL(X_col)(X_row).ews;
    const Int   ws_size = LL(X_col)(X_row).iws_size;
    const Int   p_size  = LL(X_col)(X_row).p_size;
   

    #ifdef BASKER_DEBUG_NFACTOR_BLK
    if(kid == 0 && k == 0)
    {
      printf(" %d: t_move_offdiag_L( L = LL(%d %d) X = LL(%d %d) )\n",
           kid, blkcol,blkrow, X_col, blkrow);
    }
    #endif

    Int *color   = &(ws(0));
    Int *pattern = &(color[ws_size]);

    const Int    llnnz = L.nnz;
          Int    lnnz  = L.col_ptr(k);

    if((p_size) > (llnnz-lnnz))
    {

      Int newsize = llnnz*1.2 + L.ncol;

      if (Options.verbose == BASKER_TRUE)
      {
        printf("-Warning, Need to remalloc L: %ld %ld kid: %ld current size: %ld used_size: %ld  addition: %ld \n",
            (long)blkcol, (long)blkrow, (long)kid, (long)llnnz, (long)lnnz, (long)p_size  );
      }

      thread_array(kid).error_blk    = blkcol;
      thread_array(kid).error_subblk = blkrow;
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
      //BASKER_ASSERT(0==1, "REALLOC LOWER BLOCK\n");
    }

    for(Int i = 0; i < p_size; i++)
    {
      Int j = pattern[i];

#ifdef BASKER_DEBUG_NFACTOR_BLK
      printf("L-Moving, kid: %d j: %d val: %f lnnz: %d \n",
          kid, j, X[j]/pivot, lnnz);
#endif

      color[j] = 0;
      L.row_idx(lnnz) = j;
      L.val(lnnz) = EntryOP::divide(X(j),pivot);
      X(j) = 0;

      #ifdef MY_DEBUG_BASKER
      if (kid == debug_kid) {
        printf( " L.val(%d) = %e at (%d,%d)\n",lnnz, L.val(lnnz),L.row_idx(lnnz),k );
      }
      #endif

      //Note come back to fix
#ifdef BASKER_INC_LVL
      L.inc_lvl[lnnz] = INC_LVL_TEMP[j];
#endif

      lnnz++;
    }
    //printf("L-Moving, kid: %d col_ptr(%d): %d \n",
    //           kid, k-bcol+1, lnnz);
    
    L.col_ptr(k+1) = lnnz;
    #ifdef MY_DEBUG_BASKER
    if (kid == debug_kid) {
      printf( " > L.colptr(%d) = %d\n\n",k+1,lnnz );
    }
    #endif

    LL(X_col)(X_row).p_size = 0;

    return 0;
  }//end t_offdiag_mov_L()

  
  template <class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::t_dense_back_solve_offdiag
  (
   Int kid, 
   Int blkcol, Int blkrow,
   Int X_col,  Int X_row,
   Int k,      Int &view_offset,
   ENTRY_1DARRAY x, 
   INT_1DARRAY   x_idx,
   Int x_size, Int x_offset,
   BASKER_BOOL A_option)
  {
    //Note:  need to add support for offdiag permuation
    BASKER_MATRIX &L            = LL(blkcol)(blkrow);
    BASKER_MATRIX &B            = ALM(blkcol)(blkrow);

    INT_1DARRAY   ws            = LL(X_col)(X_row).iws;
    ENTRY_1DARRAY X             = LL(X_col)(X_row).ews;
    
    Int nnz = LL(X_col)(X_row).p_size;
  
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    if (kid == 0 && X_row == 1) {
      printf( " %d: t_dense_back_solve_offdiag ( L = LL(%d,%d), X = LL(%d,%d), and B = ALM(%d,%d)\n", kid, blkcol,blkrow,X_col,X_row, blkcol,blkrow );
    }

    Int         ws_size         = LL(X_col)(X_row).iws_size;
    const Int    brow           = L.srow;
    const Int    bcol           = L.scol;
    printf("\n\n");
    printf("t_back_solve_diag, kid: %d, blkcol: %d blkrow: %d k: %d \n",
           kid, blkcol, blkrow, k);
    printf("t_back_solve_diag, kid: %d, Xcol: %d Xrow: %d \n",
           kid, X_col, X_row);
    printf("t_back_solve_diag, kid: %d, brow = %d bcol = %d\n",
           kid, brow, bcol);
    printf("t_back_solve_diag, kid: %d, ws: %d starting psize: %d \n",
           kid, ws_size, nnz);
    printf("t_back_solve_diag, kid: %d, ALM(%d)(%d): %dx%d\n",kid,blkcol,blkrow,B.nrow,B.ncol );
    printf("t_back_solve_diag, kid: %d,  LL(%d)(%d): %dx%d, nnz=%d, X.nnz=%d\n",kid,blkcol,blkrow,L.nrow,L.ncol,LL(blkcol)(blkrow).nnz,X.extent(0) );
    printf("\n\n"); fflush(stdout);
    #endif
    //B.info();
    //B.print();

    //Preload with A
    //printf(" t_back_solve(kid=%d): A_OPTION = %d, accumlate update into X = LL(col=%d,row=%d).ems, k=%d\n",kid,A_option,X_col,X_row,k);
    if(A_option == BASKER_TRUE)
    {
      for(Int i = B.col_ptr(k); i < B.col_ptr(k+1); ++i)
      {
        const Int j = B.row_idx(i);
        #ifdef BASKER_DEBUG_NFACTOR_BLK
        //Bgood(remove)
        //printf("t_back_solve_diag, kid: %d i: %d g: %d\n",
        //           kid, i, B.good(i));
        printf("t_back_solve_d, add A, kid: %d psize:%d \n",
            kid, nnz);
        printf("t_back_solve_diag, kid: %d A(%d) %f \n",
            kid, B.row_idx(i), B.val(i));
        printf("t_back_solve_diag, kid: %d x: %f %f \n",
            kid, X(j), X(j)+B.val(i));
        #endif

        //if (blkcol == 2 && blkrow == 1 && j == 87) {
        //  printf( " t_dense_back_solve_offdiag:: load X(%d) = %e + %e = %e\n",j,X(j),B.val(i),X(j)+B.val(i) );
        //}
        X(j) = X(j) + B.val(i);

      }//over all nnz in subview
    }//end if preload
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("t_back_solve_d, kid: %d xsize: %ld \n",
           kid, x_size);
    #endif

    //SPMV
    //printf("x_offset: %d \n", x_offset);
    const Entry zero (0.0);
    for(Int i = 0; i < x_size; ++i)
    {
      const Int colid = x_idx[i+x_offset];
      const Entry xj  = x(i+x_offset);
      //if (blkcol == 2 && blkrow == 1) {
      //  printf("t_back_solve_diag kid = %d: i = %d / %d k = %d, xj = %e \n", kid, i,x_size, k, xj); fflush(stdout);
      //}
      #ifdef BASKER_DEBUG_NFACTOR_BLK
      printf("t_back_solve_diag, kid: %d  k: %d %g  x_size: %d [%d %d] \n",
             kid, colid, xj, x_size,  L.col_ptr[k], L.col_ptr[k+1]);
      #endif
      if (xj != zero)
      {
        //if (blkcol == 2 && blkrow == 1) {
        //  printf( " > L.colptr(%d) = %d : %d\n",colid,L.col_ptr(colid),L.col_ptr(colid+1)-1 );
        //}
        for(Int j = L.col_ptr(colid); j < L.col_ptr(colid+1); j++)
        {
          const Int jj = L.row_idx(j);
          //if (blkcol == 2 && blkrow == 1 && jj == 87) {
          //  printf( " > kid=%d: t_dense_back_solve_offdiag:: update X(%d) = %e - %e * %e = %e\n",kid,jj,X(jj),L.val(j),xj,X(jj)-L.val(j)*xj );
          //}
          X(jj) -= L.val(j)*xj; 

          #ifdef BASKER_DEBUG_NFACTOR_BLK
          printf("t_b_solve_d, kid: %d j: %d color: %d \n",
                 kid, jj, color[jj]);
          printf("t_back_solve_d,id:%d  row_idx: %d b4: %f mult: %f %f\n",
                 kid, jj,X[jj], L.val[j], xj);
          #endif
        }
      }
    }//over all nonzero in left

    #ifdef BASKER_2DL
    //LL(X_col)(X_row).p_size = nnz;
    LL(X_col)(X_row).p_size = nnz;
    #endif

    //Debug
    #ifdef BASKER_DEBUG_NFACTOR_BLK
     #ifdef BASKER_2DL
     printf("---PATTERN End test: kid: %d nnz: %d pattern: %d \n",
            kid, nnz, pattern[nnz-1]); 
     printf("SETTING dig PS: %d kid: %d L: %d %d\n",
            nnz, kid, X_col, X_row);
     printf("kid %d Ending nnz: %d \n",kid, nnz);
     #endif
     printf("kid: %d all values: \n", kid);
     for(Int i = 0; i < L.nrow; ++i)
     {
       Int jj = i;
       //Int jj = pattern[i];
       //printf("X[%d] = %f , kid: %d  \n",
       //     jj, X[jj-brow], kid);
       printf("X[%d](%d) = %f , kid: %d  \n",
           jj,jj+B.srow, X[jj], kid);

     }
     printf("\n\n");
    #endif

    return 0;
  }//end t_offdiag_dense_back_solve();


  //local idx for local blks
  //Bgood done.
  template <class Int, class Entry, class Exe_Space>
  int Basker<Int,Entry,Exe_Space>::t_back_solve_offdiag
  (
   Int kid, 
   Int blkcol, Int blkrow,
   Int X_col,  Int X_row,
   Int col,
   ENTRY_1DARRAY  x, 
   INT_1DARRAY    x_idx,
   Int x_size,
   Int x_offset,
   BASKER_BOOL A_option)
  {
    //Note:  need to add support for offdiag permuation
    BASKER_MATRIX &L  = LL(blkcol)(blkrow);
    BASKER_MATRIX &B  = ALM(blkcol)(blkrow);

    INT_1DARRAY   ws  = LL(X_col)(X_row).iws;
    ENTRY_1DARRAY  X  = LL(X_col)(X_row).ews;

    Int      ws_size  = LL(X_col)(X_row).iws_size;
    Int          nnz  = LL(X_col)(X_row).p_size;
  
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    if(kid == 0 and col == 0)
    {
      printf( " t_back_solve_offdiag: L = LL(%d)(%d) with %d x %d, X = LL(%d)(%d)\n",blkcol,blkrow,L.nrow,L.ncol, X_col,X_row );
      //printf( "[\n" );
      //for (Int j = 0; j < L.ncol; j++) {
      //  for (Int k = L.col_ptr(j); k < L.col_ptr(j+1); k++) printf( " %d %d %e\n",L.row_idx(k),j,L.val(k) );
      //}
      //printf( "];\n" );
    }
    #endif
    // B.info();
    //B.print();

    Int *color =   &(ws(0));
    Int *pattern = &(color[ws_size]);
    
    //Preload with A
    #ifdef MY_DEBUG_BASKER
    //if (kid == debug_kid) 
    //if (blkcol == 2 && blkrow == 1)
    {
      printf( " > t_back_solve_offdiag with LL(%d,%d) and X(%d,%d) and x_size=%d\n",blkcol,blkrow, X_col,X_row, x_size);
    }
    #endif
    if(A_option == BASKER_TRUE)
    {
      for(Int i = B.col_ptr(col); i < B.col_ptr(col+1); ++i)
      {
#ifdef BASKER_DEBUG_NFACTOR_BLK
        if(kid == 8)
        {
          printf("t_back_solve_d, add A, kid: %d psize:%d \n",
              kid, nnz);
          printf("t_back_solve_diag, kid: %d A(%d) %f \n",
              kid, B.row_idx(i), B.val(i));
        }
#endif
        const Int j = B.row_idx(i);
        color[j] = 1;
        X(j) = B.val(i);
        #ifdef MY_DEBUG_BASKER
        if (blkcol == 2 && blkrow == 1)
        {
          printf( " load: X(%d) = %e\n",j,X(j) );
        }
        #endif
        pattern[nnz++] = j;
      }//over all nnz in subview
    }//end if preload
  
    //SPMV
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    if(kid == 0)
    {
      printf("t_back_solve_d, kid: %d xsize: %d \n",
               (int)kid, (int)x_size);
    }
    #endif
    const Entry zero (0.0);
    for(Int i = 0; i < x_size; ++i)
    {
      const Int   k  = x_idx[i+x_offset];
      const Entry xj = x(i+x_offset);
      if(k >= col || xj == zero)
      {
        // skip diagonal element
        continue;
      }

      #ifdef BASKER_DEBUG_NFACTOR_BLK
      if(kid == 8)
      {
        printf("t_back_solve_diag, kid: %d k: %d [%d %d] \n",
            kid, k, L.col_ptr[k], L.col_ptr[k+1]);
      }
      #endif

      #ifdef MY_DEBUG_BASKER
      if (blkcol == 2 && blkrow == 1)
      {
        printf( " L.col_ptr(i=%d -> k=%d): %d:%d\n", (int)i,(int)k,(int)L.col_ptr(k),(int)L.col_ptr(k+1));
      }
      #endif
      for(Int j = L.col_ptr(k); j < L.col_ptr(k+1); j++)
      {
        const Int jj = L.row_idx(j);
        #ifdef BASKER_DEBUG_NFACTOR_BLK
        if (blkcol == 2 && blkrow == 1)
        {
          printf("t_b_solve_d, kid: %d j: %d color: %d \n",
              kid, jj, color[jj]);
        }
        #endif
        if(color[jj] != 1)
        {
          color[jj] = 1;
          pattern[nnz++] = jj;
          #ifdef BASKER_DEBUG_NFACTOR_BLK
          if(kid == 8)
          {
            printf("pattern index: %d kid: %d \n",
                nnz, kid);
          }
          #endif
        }

        #ifdef BASKER_DEBUG_NFACTOR_BLK
        //if(kid == 8)
        if (blkcol == 2 && blkrow == 1)
        {
          printf("t_back_solve_d,id:%d  row_idx: %d b4: %f mult: %f %f\n",
              kid, jj,X[jj], L.val[j], xj);
        }
        #endif 

        X(jj) -= L.val(j)*xj;
        #ifdef MY_DEBUG_BASKER
        if (kid == debug_kid) {
          printf( " -> X(%d) -= %e * %e\n",jj, L.val(j),xj );
        }
        #endif
      }
    }//over all nonzero in left

    #ifdef BASKER_2DL
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    printf("---PATTERN End test: kid: %d nnz: %d pattern: %d \n",
           kid, nnz, pattern[nnz-1]); 
    printf("SETTING dig PS: %d kid: %d L: %d %d\n",
           nnz, kid, X_col, X_row);
    printf("kid %d Ending nnz: %d \n",kid, nnz);
    #endif
    //LL(X_col)(X_row).p_size = nnz;
    LL(X_col)(X_row).p_size = nnz;
    #endif

    //Debug
    #ifdef BASKER_DEBUG_NFACTOR_BLK
    if(kid == 8)
    {
      printf("kid: %d all values: \n", kid);
      for(Int i = 0; i < nnz; ++i)
      {
        Int jj = pattern[i];
        //printf("X[%d] = %f , kid: %d  \n",
        //     jj, X[jj-brow], kid);
        printf("X[%d](%d) = %f , kid: %d  \n",
            jj,jj+B.srow, X[jj], kid);

      }
      printf("\n\n");
    }
    #endif
      
    return 0;
  }//end t_offdiag_back_solve();

}//end namespace baskerNS -- contains functions

#undef BASKER_DETAILED_TIMER
#undef BASKER_TIMER
#undef MY_DEBUG_BASKER
#endif//end ifndef basker_nfactor_blk
