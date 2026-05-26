// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLUBASKER_NFACTOR_HPP
#define SHYLUBASKER_NFACTOR_HPP

/*Basker Includes*/
#include "shylubasker_types.hpp"
#include "shylubasker_util.hpp"
#include "shylubasker_structs.hpp"
#include "shylubasker_matrix_view_decl.hpp"
#include "shylubasker_matrix_view_def.hpp"

#include "shylubasker_nfactor_dom.hpp"
#include "shylubasker_nfactor_sep2.hpp"
#include "shylubasker_nfactor_btf.hpp"

#include "shylubasker_error_manager.hpp"

/*Kokkos Includes*/
#include <Kokkos_Core.hpp>
#include <Kokkos_Timer.hpp>

/*System Includes*/
#include <iostream>
#include <string>

//#define BASKER_DEBUG_NFACTOR 
//#define BASKER_TIMER

namespace BaskerNS
{
  
  template <class Int, class Entry, class Exe_Space>
  struct kokkos_nfactor_funct
  {
    typedef Exe_Space                         execution_space;
    typedef Kokkos::TeamPolicy<Exe_Space>     TeamPolicy;
    typedef typename TeamPolicy::member_type  TeamMember;
    
    Basker<Int,Entry,Exe_Space> *basker;
    
    kokkos_nfactor_funct()
    {}

    kokkos_nfactor_funct(Basker<Int,Entry,Exe_Space> *_basker)
    {basker = _basker;}

    BASKER_INLINE
    void operator()(const TeamMember &thread) const
    {}//end operator()

  };//end basker_nfactor_funct

  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::factor_notoken(Int option)
  {

    //printf("factor no token called \n");

    gn = A.ncol;
    gm = A.nrow;
    BASKER_MATRIX ATEMP;

    //Kokkos::Timer tza;
    int info = BASKER_SUCCESS;
    if(Options.btf == BASKER_TRUE)
    {
      //JDB: We can change this for the new inteface
      gn = A.ncol;
      gm = A.nrow;
      ATEMP = A;
      A = BTF_A;
    }
    //printf("Switch time: %f \n", tza.seconds());

    // --------------------------------------------------------------- //
    // ----------------------- Big A block --------------------------- //
    Kokkos::Timer timer;
    if(btf_tabs_offset != 0)
    {
      //Spit into Domain and Sep

      // -------------------------------------------------------- //
      // -----------------------Domain--------------------------- //
      Int nworker_threads = (Options.worker_threads ? 2 : 1);
      Int num_doms = num_threads;
      if(Options.verbose == BASKER_TRUE)
      {
        printf("Factoring Dom(# teams = %ld with # threads = %ld)\n", (long)num_doms,(long)nworker_threads);
        //for (int tid = 0; tid < num_doms; tid ++) {
        //  printf( " error[%d] = %d\n",tid,thread_array(tid).error_type );
        //}
        fflush(stdout);
        timer.reset();
      }
      //if(Options.partial_facto != 0) {
      //  if(Options.verbose == BASKER_TRUE) printf( " > halving # teams for partial factorization (%d -> %d)\n",num_doms,num_doms/2 );
      //  num_doms /= 2;
      //}

      Int domain_restart = 0;
      kokkos_nfactor_domain <Int,Entry,Exe_Space>
        domain_nfactor(this);
      Kokkos::parallel_for(TeamPolicy(num_doms,nworker_threads),
          domain_nfactor);
      Kokkos::fence();

      //=====Check for error======
      while(info != BASKER_ERROR)
      {
        INT_1DARRAY thread_start;
        MALLOC_INT_1DARRAY(thread_start, num_doms+1);
        init_value(thread_start, num_doms+1, (Int) BASKER_MAX_IDX);

        info = nfactor_domain_error(thread_start);
        if(Options.verbose == BASKER_TRUE) {
          printf( " nfactor_domain: info = %d\n",(int)info );
        }
        if(info == BASKER_SUCCESS)
        {
          break;
        }
        else if(domain_restart > BASKER_RESTART)
        {
          if(Options.verbose == BASKER_TRUE)
          {
            printf(" nfactor_domain_error reports info=%d and max restartt reached (%d)\n", (int)info, (int)domain_restart);
          }
          break;
        }
        else if (info == BASKER_ERROR)
        {
          if(Options.verbose == BASKER_TRUE)
          {
            printf(" nfactor_domain_error reports BASKER_ERROR - numeric factorization failed\n");
          }
          break;
        }
        else
        {
          domain_restart++;
          if(Options.verbose == BASKER_TRUE)
          {
            printf("\n restart factorization\n");
          }
          kokkos_nfactor_domain_remalloc <Int, Entry, Exe_Space>
            nfactor_dom_remalloc(this, thread_start);
          Kokkos::parallel_for(TeamPolicy(num_doms,nworker_threads),
              nfactor_dom_remalloc);
          Kokkos::fence();
        }
      }//end while
      if(Options.verbose == BASKER_TRUE)
      {
        printf("Time DOMAIN: %lf \n\n", timer.seconds());
        timer.reset();
      }
      //printVec("domperm.csc", gpermi, A.nrow);
      //printf( " End Dom: info = %d (%d, %d)\n",info,BASKER_SUCCESS,BASKER_ERROR );
      //-------------------End--Domian--------------------------//


      // -------------------------------------------------------- //
      // ---------------------------Sep-------------------------- //
      if(info == BASKER_SUCCESS) {
        if(Options.partial_facto == 2) {
          // zero out schur complement
          for (Int i=0; i<schur_size*schur_size; i++) schur_out_ptr[i] = Entry(0.0);
        }
        #ifdef BASKER_NO_LAMBDA
        for(Int l=1; l <= tree.nlvls; l++)
        {
          #ifdef BASKER_TIMER
          Kokkos::Timer timer_inner_sep;
          #endif
          //#ifdef BASKER_OLD_BARRIER
          //Int lthreads = pow(2,l);
          //Int lnteams = num_threads/lthreads;
          //#else
          Int lthreads = 1;
          Int lnteams = num_threads/lthreads;
          //#endif

          Int sep_restart = 0;
          if(Options.verbose == BASKER_TRUE)
          {
            printf("Factoring Sep(# teams = %ld with # threads = %ld) at level = %d\n",
                   (long)lnteams, (long)lthreads, (int)l); fflush(stdout);
          }
          //if(Options.partial_facto != 0) {
          //  if(Options.verbose == BASKER_TRUE) printf( " > halving # teams for partial factorization (%d -> %d)\n",lnteams,lnteams/2 );
          //  lnteams /= 2;
          //}

          kokkos_nfactor_sep2 <Int, Entry, Exe_Space> sep_nfactor(this, l, tree.nlvls);
          Kokkos::parallel_for(TeamPolicy(lnteams,lthreads),
                               sep_nfactor);
          Kokkos::fence();

          //======Check for error=====
          while(info != BASKER_ERROR)
          {
            INT_1DARRAY thread_start;
            MALLOC_INT_1DARRAY(thread_start, num_threads+1);
            init_value(thread_start, num_threads+1, (Int) BASKER_MAX_IDX);

            info = nfactor_sep_error(thread_start);
            if(info == BASKER_SUCCESS)
            {
              FREE_INT_1DARRAY(thread_start);
              break;
            }
            else if (sep_restart > BASKER_RESTART)
            {
              if(Options.verbose == BASKER_TRUE)
              {
                printf("%s: nfactor_separator_error reports info=%d and max restartt reached (%d)\n",__FILE__, info, (int)sep_restart);
              }
              break;
            }
            else if (info == BASKER_ERROR)
            {
              if(Options.verbose == BASKER_TRUE)
              {
                printf("%s: nfactor_separator_error reports BASKER_ERROR - numeric factorization failed\n",__FILE__);
              }
              break;
            }
            else
            {
              sep_restart++;
              if (Options.verbose == BASKER_TRUE)
              {
                printf("\n restart factorization\n\n");
              }
              Kokkos::parallel_for(TeamPolicy(lnteams,lthreads), sep_nfactor);
              //printf( "\n ***** done sep_nfactor *****\n\n" ); fflush(stdout);
              Kokkos::fence();
            }
          }//end while-true
          #ifdef BASKER_TIMER
          printf(" > Time INNERSEP %ld: %lf \n", (long int)l, timer_inner_sep.seconds());
          #endif
        }//end over each level
        if(Options.partial_facto == 2 && schur_out_ptr != nullptr) {
          // copy out schur complement
          for (Int i=0; i<schur_size; i++) {
            for (Int j=0; j<schur_size; j++) schur_out_ptr[i+j*schur_size] = schur_out(i,j);
          }
        }
        #else //ELSE BASKER_NO_LAMBDA
        //Note: to be added
        #endif //end BASKER_NO_LAMBDA
        //-------------------------End Sep----------------//
      }// info != BASKER_ERROR

      //printf( " End Sep: info = %d (%d, %d)\n",info,BASKER_SUCCESS,BASKER_ERROR );
      if(Options.verbose == BASKER_TRUE)
      {
        printf("Time SEP: %lf \n\n", timer.seconds());
        timer.reset();
      }
    }

    // ---------------------------------------------------------------------------------------- //
    // ----------------------- Small blocks in bottom C (and top D) --------------------------- //
    if(info == BASKER_SUCCESS && Options.btf == BASKER_TRUE)
    {
      Int btf_restart = 0;

      if(Options.verbose == BASKER_TRUE)
      {
        printf("Factoring BLKs num_threads: %ld \n",
            (long)num_threads);
        timer.reset();
      }

      //======Call BTF factor====
      kokkos_nfactor_btf <Int, Entry, Exe_Space> 
        nfactor_btf(this);
      Kokkos::parallel_for(TeamPolicy(num_threads,1),
        nfactor_btf);
      Kokkos::fence();

      //=====Check for error======
      while(info != BASKER_ERROR)
      {
        INT_1DARRAY thread_start;
        MALLOC_INT_1DARRAY(thread_start, num_threads+1);
        init_value(thread_start, num_threads+1, (Int) BASKER_MAX_IDX);

        INT_1DARRAY thread_start_top;
        MALLOC_INT_1DARRAY(thread_start_top, num_threads+1);
        init_value(thread_start_top, num_threads+1, (Int) BASKER_MAX_IDX);

        info = nfactor_btf_error(thread_start_top, thread_start);
        //printf( " nfactor_btf: info = %d\n\n",(int)info );
        //printf("RETURNED: %d (success=%d, error=%d)\n", info, BASKER_SUCCESS, BASKER_ERROR);
        if(info == BASKER_SUCCESS)
        {
          break;
        }
        else if (btf_restart > BASKER_RESTART)
        {
          if(Options.verbose == BASKER_TRUE)
          {
            printf("%s: nfactor_btfonal_error reports info=%d and max restartt reached (%d)\n",__FILE__, info, (int)btf_restart);
          }
          break;
        }
        else if (info == BASKER_ERROR)
        {
          if(Options.verbose == BASKER_TRUE)
          {
            printf("%s: nfactor_btfonal_error reports BASKER_ERROR - numeric factorization failed\n",__FILE__);
          }
          break;
        }
        else
        {
          btf_restart++;
          if (Options.verbose == BASKER_TRUE)
          {
            printf("\n restart factorization\n");
          }
          kokkos_nfactor_btf_remalloc <Int, Entry, Exe_Space>
            nfactor_btf_remalloc(this, thread_start_top, thread_start);
          Kokkos::parallel_for(TeamPolicy(num_threads,1),
            nfactor_btf_remalloc);
          Kokkos::fence();
          //printf( " nfactor_btf_remalloc done\n\n" );
        }
      }//end while

      if(Options.verbose == BASKER_TRUE)
      {
        printf("Time BTF: %lf \n\n", timer.seconds());
      }
    }//end btf call

    Kokkos::Timer tzback;
    if(Options.btf == BASKER_TRUE)
    {
      A = ATEMP;
    }
    //printf("Switch back: %f \n",
    //            tzback.seconds());

    return info;
  }//end factor_notoken()
  
}//end namespace baskerNS

#undef BASKER_TIMER
#endif //end ifndef basker_nfactor_hpp
