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

#include "shylubasker_nfactor_blk.hpp"
#include "shylubasker_nfactor_col.hpp"
#include "shylubasker_nfactor_col2.hpp"
#include "shylubasker_nfactor_diag.hpp"

#include "shylubasker_error_manager.hpp"

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

//#define BASKER_DEBUG_NFACTOR 
//#define BASKER_TIMER

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

#ifdef BASKER_KOKKOS
    Kokkos::Timer timer;

    typedef Kokkos::TeamPolicy<Exe_Space>        TeamPolicy;
    // --------------------------------------------------------------- //
    // ----------------------- Big A block --------------------------- //
    if(btf_tabs_offset != 0)
    {
      //Spit into Domain and Sep

      // -------------------------------------------------------- //
      // -----------------------Domain--------------------------- //
      if(Options.verbose == BASKER_TRUE)
      {
        printf("Factoring Dom num_threads: %ld \n", (long)num_threads);
        //for (int tid = 0; tid < num_threads; tid ++) {
        //  printf( " error[%d] = %d\n",tid,thread_array(tid).error_type );
        //}
        fflush(stdout);
        timer.reset();
      }

      Int domain_restart = 0;
      kokkos_nfactor_domain <Int,Entry,Exe_Space>
        domain_nfactor(this);
      Kokkos::parallel_for(TeamPolicy(num_threads,1), 
          domain_nfactor);
      Kokkos::fence();

      //=====Check for error======
      while(info != BASKER_ERROR)
      {
        INT_1DARRAY thread_start;
        MALLOC_INT_1DARRAY(thread_start, num_threads+1);
        init_value(thread_start, num_threads+1, (Int) BASKER_MAX_IDX);

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
            diag_nfactor_remalloc(this, thread_start);
          Kokkos::parallel_for(TeamPolicy(num_threads,1),
              diag_nfactor_remalloc);
          Kokkos::fence();
        }
      }//end while
      if(Options.verbose == BASKER_TRUE)
      {
        printf("Time DOMAIN: %lf \n", timer.seconds());
        timer.reset();
      }
      #ifdef BASKER_TIMER
      printf("Time DOMAIN: %lf \n", timer.seconds());
      timer.reset();
      #endif

#else// else basker_kokkos
      #pragma omp parallel
      {

      }//end omp parallel
#endif //end basker_kokkos
      //printVec("domperm.csc", gpermi, A.nrow);
      //printf( " End Dom: info = %d (%d, %d)\n",info,BASKER_SUCCESS,BASKER_ERROR );
      //-------------------End--Domian--------------------------//


      // -------------------------------------------------------- //
      // ---------------------------Sep-------------------------- //
      if(info == BASKER_SUCCESS) {
#ifdef BASKER_KOKKOS
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
            printf("Factoring Sep(lnteams = %ld, lthreads = %ld)\n",
                   (long)lnteams, (long)lthreads);
          }

          kokkos_nfactor_sep2 <Int, Entry, Exe_Space> sep_nfactor(this, l);
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
            //printf( "\n ***** nfactor_separator: info = %d *****\n",(int)info );
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
        #else //ELSE BASKER_NO_LAMBDA
        //Note: to be added
        #endif //end BASKER_NO_LAMBDA
#else

        #pragma omp parallel
        {

        }//end omp parallel
#endif
        //-------------------------End Sep----------------//
      }// info != BASKER_ERROR
      //printf( " End Sep: info = %d (%d, %d)\n",info,BASKER_SUCCESS,BASKER_ERROR );
      if(Options.verbose == BASKER_TRUE)
      {
        printf("Time SEP: %lf \n", timer.seconds());
        timer.reset();
      }
      #ifdef BASKER_TIMER
      printf("Time SEP: %lf \n", timer.seconds());
      timer.reset();
      #endif
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

      //======Call diag factor====
      kokkos_nfactor_diag <Int, Entry, Exe_Space> 
        diag_nfactor(this);
      Kokkos::parallel_for(TeamPolicy(num_threads,1),
        diag_nfactor);
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

        info = nfactor_diag_error(thread_start_top, thread_start);
        //printf( " nfactor_diag: info = %d\n\n",(int)info );
        //printf("RETURNED: %d (success=%d, error=%d)\n", info, BASKER_SUCCESS, BASKER_ERROR);
        if(info == BASKER_SUCCESS)
        {
          break;
        }
        else if (btf_restart > BASKER_RESTART)
        {
          if(Options.verbose == BASKER_TRUE)
          {
            printf("%s: nfactor_diagonal_error reports info=%d and max restartt reached (%d)\n",__FILE__, info, (int)btf_restart);
          }
          break;
        }
        else if (info == BASKER_ERROR)
        {
          if(Options.verbose == BASKER_TRUE)
          {
            printf("%s: nfactor_diagonal_error reports BASKER_ERROR - numeric factorization failed\n",__FILE__);
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
          kokkos_nfactor_diag_remalloc <Int, Entry, Exe_Space>
            diag_nfactor_remalloc(this, thread_start_top, thread_start);
          Kokkos::parallel_for(TeamPolicy(num_threads,1),
            diag_nfactor_remalloc);
          Kokkos::fence();
          //printf( " diag_nfactor_remalloc done\n\n" );
        }
      }//end while

      if(Options.verbose == BASKER_TRUE)
      {
        printf("Time BTF: %lf \n", timer.seconds());
      }
      #ifdef BASKER_TIMER
      printf("Time BTF: %lf \n", timer.seconds());
      #endif
    }//end btf call

    Kokkos::Timer tzback;
    if(Options.btf == BASKER_TRUE)
    {
      A = ATEMP;
    }
    //printf("Switch back: %f \n",
    //	    tzback.seconds());

    return info;
  }//end factor_notoken()
  
}//end namespace baskerNS

#undef BASKER_TIMER
#endif //end ifndef basker_nfactor_hpp
