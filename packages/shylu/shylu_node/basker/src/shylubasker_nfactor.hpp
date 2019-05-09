#ifndef SHYLUBASKER_NFACTOR_HPP
#define SHYLUBASKER_NFACTOR_HPP

//#define BASKER_DEBUG_NFACTOR 

//#define BASKER_TIME

/*Basker Includes*/
#include "shylubasker_types.hpp"
#include "shylubasker_util.hpp"
#include "shylubasker_structs.hpp"
#include "shylubasker_matrix_view_decl.hpp"
#include "shylubasker_matrix_view_def.hpp"

#include "shylubasker_nfactor_blk.hpp"
#include "shylubasker_nfactor_col.hpp"
#include "shylubasker_nfactor_col2.hpp"
#include "shylubasker_nfactor_diag.hpp"

#include "shylubasker_error_manager.hpp"

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

namespace BaskerNS
{
  
  template <class Int, class Entry, class Exe_Space>
  struct kokkos_nfactor_funct
  {
    #ifdef BASKER_KOKKOS
    typedef Exe_Space                         execution_space;
    typedef Kokkos::TeamPolicy<Exe_Space>     TeamPolicy;
    typedef typename TeamPolicy::member_type  TeamMember;
    #endif 
    
    Basker<Int,Entry,Exe_Space> *basker;
    
    kokkos_nfactor_funct()
    {}

    kokkos_nfactor_funct(Basker<Int,Entry,Exe_Space> *_basker)
    {basker = _basker;}

    BASKER_INLINE
    #ifdef BASKER_KOKKOS
    void operator()(const TeamMember &thread) const
    #else
    void operator()(Int kid) const
    #endif
    {

    }//end operator()

  };//end basker_nfactor_funct

  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::factor_notoken(Int option)
  {

    //printf("factor no token called \n");

    gn = A.ncol;
    gm = A.nrow;
    BASKER_MATRIX ATEMP;

    //Kokkos::Impl::Timer tza;
    if(Options.btf == BASKER_TRUE)
    {
      //JDB: We can change this for the new inteface
      gn = A.ncol;
      gm = A.nrow;
      ATEMP = A;
      A = BTF_A; 
    }
    //printf("Switch time: %f \n", tza.seconds());



    //Spit into Domain and Sep
    //----------------------Domain-------------------------//
#ifdef BASKER_KOKKOS
    //====TIMER==
#ifdef BASKER_TIME
    Kokkos::Impl::Timer       timer;
#endif
    //===TIMER===

    typedef Kokkos::TeamPolicy<Exe_Space>        TeamPolicy;

    if(btf_tabs_offset != 0)
    {

      if(Options.verbose == BASKER_TRUE)
      {
        printf("Factoring Dom num_threads: %ld \n",
            (long)num_threads);
      }

      Int domain_restart = 0;
      kokkos_nfactor_domain <Int,Entry,Exe_Space>
        domain_nfactor(this);
      Kokkos::parallel_for(TeamPolicy(num_threads,1), 
          domain_nfactor);
      Kokkos::fence();

      //=====Check for error======
      while(true)
      {
        INT_1DARRAY thread_start;
        MALLOC_INT_1DARRAY(thread_start, num_threads+1);
        init_value(thread_start, num_threads+1, 
            (Int) BASKER_MAX_IDX);
        int nt = nfactor_domain_error(thread_start);
        if((nt == BASKER_SUCCESS) ||
            (nt == BASKER_ERROR) ||
            (domain_restart > BASKER_RESTART))
        {
          break;
        }
        else
        {
          domain_restart++;
          if(Options.verbose == BASKER_TRUE)
          {
            printf("restart \n");
          }
          kokkos_nfactor_domain_remalloc <Int, Entry, Exe_Space>
            diag_nfactor_remalloc(this, thread_start);
          Kokkos::parallel_for(TeamPolicy(num_threads,1),
              diag_nfactor_remalloc);
          Kokkos::fence();
        }
      }//end while

      //====TIMER===
#ifdef BASKER_TIME
      printf("Time DOMAIN: %lf \n", timer.seconds());
      timer.reset();
#endif
      //====TIMER====


#else// else basker_kokkos
#pragma omp parallel
      {

      }//end omp parallel
#endif //end basker_kokkos

    }
    //-------------------End--Domian--------------------------//

    //printVec("domperm.csc", gpermi, A.nrow);

    //---------------------------Sep--------------------------//

    if(btf_tabs_offset != 0)
    {
      //for(Int l=1; l<=4; l++)
      for(Int l=1; l <= tree.nlvls; l++)
      {

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
          printf("Factoring Sep num_threads: %ld %ld \n",
              (long)lnteams, (long)lthreads);
        }

#ifdef BASKER_KOKKOS
        Kokkos::Impl::Timer  timer_inner_sep;
  #ifdef BASKER_NO_LAMBDA

        //kokkos_nfactor_sep <Int, Entry, Exe_Space> 
        //sep_nfactor(this, l);

        kokkos_nfactor_sep2 <Int, Entry, Exe_Space>
          sep_nfactor(this,l);

        Kokkos::parallel_for(TeamPolicy(lnteams,lthreads),
            sep_nfactor);
        Kokkos::fence();

        //======Check for error=====
        while(true)
        {
          INT_1DARRAY thread_start;
          MALLOC_INT_1DARRAY(thread_start, num_threads+1);
          init_value(thread_start, num_threads+1,
              (Int) BASKER_MAX_IDX);
          int nt = nfactor_sep_error(thread_start);
          if((nt == BASKER_SUCCESS)||
              (nt == BASKER_ERROR) ||
              (sep_restart > BASKER_RESTART))
          {
            FREE_INT_1DARRAY(thread_start);
            break;
          }
          else
          {
            sep_restart++;
            if (Options.verbose == BASKER_TRUE)
            {
              printf("restart \n");
            }
            Kokkos::parallel_for(TeamPolicy(lnteams,lthreads),  sep_nfactor);
            Kokkos::fence();

          }
        }//end while-true


    #ifdef BASKER_TIME
        printf("Time INNERSEP: %ld %lf \n", 
            (long)l, timer_inner_sep.seconds());
    #endif
  #else //ELSE BASKER_NO_LAMBDA
        //Note: to be added
  #endif //end BASKER_NO_LAMBDA

#else

#pragma omp parallel
        {

        }//end omp parallel
#endif
      }//end over each level

#ifdef BASKER_TIME
      printf("Time SEP: %lf \n", timer.seconds());
#endif
    }

    //-------------------------End Sep----------------//


    //-------------------IF BTF-----------------------//
    if(Options.btf == BASKER_TRUE)
    {
      Int btf_restart = 0;

      if(Options.verbose == BASKER_TRUE)
      {
        printf("Factoring BLKs num_threads: %ld \n",
            (long)num_threads);
      }

      //=====Timer
#ifdef BASKER_TIME
      Kokkos::Impl::Timer  timer_btf;
#endif
      //====Timer

      //======Call diag factor====
      kokkos_nfactor_diag <Int, Entry, Exe_Space> 
        diag_nfactor(this);
      Kokkos::parallel_for(TeamPolicy(num_threads,1),
        diag_nfactor);
      Kokkos::fence();

      //=====Check for error======
      while(true)
      {
        INT_1DARRAY thread_start;
        MALLOC_INT_1DARRAY(thread_start, num_threads+1);
        init_value(thread_start, num_threads+1, 
            (Int) BASKER_MAX_IDX);
        int nt = nfactor_diag_error(thread_start);
        //printf("RETURNED: %d \n", nt);
        if((nt == BASKER_SUCCESS) || 
            (nt == BASKER_ERROR) ||
            (btf_restart > BASKER_RESTART))
        {
          break;
        }
        else
        {
          btf_restart++;
          if (Options.verbose == BASKER_TRUE)
          {
            printf("restart \n");
          }
          kokkos_nfactor_diag_remalloc <Int, Entry, Exe_Space>
            diag_nfactor_remalloc(this, thread_start);
          Kokkos::parallel_for(TeamPolicy(num_threads,1),
            diag_nfactor_remalloc);
          Kokkos::fence();
        }
      }//end while

      //====TIMER
#ifdef BASKER_TIME
      printf("Time BTF: %lf \n", 
          timer_btf.seconds());
#endif
      //===TIMER

    }//end btf call

    Kokkos::Impl::Timer tzback;
    if(Options.btf == BASKER_TRUE)
    {
      A = ATEMP;
    }
    //printf("Switch back: %f \n",
    //	    tzback.seconds());

    return 0;
  }//end factor_notoken()
  
}//end namespace baskerNS
#endif //end ifndef basker_nfactor_hpp
