// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLUBASKER_NFACTOR_INC_HPP
#define SHYLUBASKER_NFACTOR_INC_HPP

//#define BASKER_TIME

/*Basker Includes*/
#include "shylubasker_types.hpp"
#include "shylubasker_util.hpp"
#include "shylubasker_structs.hpp"
#include "shylubasker_matrix_view_decl.hpp"
#include "shylubasker_matrix_view_def.hpp"

#include "shylubasker_nfactor_blk_inc.hpp"
#include "shylubasker_nfactor_col_inc.hpp"

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

namespace BaskerNS
{  
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::factor_inc_lvl(Int option)
  {

    //printf("Factor Inc Level Called \n");
    BASKER_BOOL fatal_error = BASKER_FALSE;
    
    gn = A.ncol;
    gm = A.nrow;

    if(Options.btf == BASKER_TRUE)
      {
	//call reference copy constructor
	//JDB check if reference constructor or deep-copy is called
	gn = A.ncol;
	gm = A.nrow;
	A = BTF_A; 
      }
   
  
    //----------------------Domain NFactor-------------------------//
    #ifdef BASKER_KOKKOS

    //====TIMER==
    #ifdef BASKER_TIME
    //Kokkos::Timer       timer;
    #endif
    //===TIMER===

    typedef Kokkos::TeamPolicy<Exe_Space>        TeamPolicy;
    if(btf_tabs_offset != 0)
      {
	//printf("domain\n");
	Int domain_restart = 0;
	kokkos_nfactor_domain_inc_lvl <Int,Entry,Exe_Space>
	  domain_nfactor(this);
	Kokkos::parallel_for(TeamPolicy(num_threads,1), 
			     domain_nfactor);
	Kokkos::fence();
        //printf("done domain\n");
    
	//=====Check for error======
	while(true)
	  {
	INT_1DARRAY thread_start;
	MALLOC_INT_1DARRAY(thread_start, num_threads+1);
	init_value(thread_start, num_threads+1, 
		   (Int) BASKER_MAX_IDX);
	int nt = nfactor_domain_error(thread_start);
	if((nt == BASKER_SUCCESS))
	  {
	    FREE_INT_1DARRAY(thread_start);
	    break;
	  }
	else if((nt == BASKER_ERROR) || 
		(domain_restart > BASKER_RESTART))
	  {
	    //This error is not recoverable
	    //Therefore, we should not try any other work
	    FREE_INT_1DARRAY(thread_start);
	    if(Options.verbose == BASKER_TRUE)
	      {
	    printf("DOM FATAL ERROR\n");
	      }
	    fatal_error = BASKER_TRUE;
	    break;
	  }
	else
	  {
	    domain_restart++;
	    if(Options.verbose == BASKER_TRUE)
	      {
	    printf("dom restart \n");
	      }
	    kokkos_nfactor_domain_remalloc_inc_lvl <Int, Entry, Exe_Space>
	      diag_nfactor_remalloc(this, thread_start);
	    Kokkos::parallel_for(TeamPolicy(num_threads,1),
				 diag_nfactor_remalloc);
	    Kokkos::fence();
	  }
      }//end while
   
    //====TIMER===
    #ifdef BASKER_TIME
	//printf("Time DOMAIN: %f \n", timer.seconds());
	//timer.reset();
    #endif
    //====TIMER====
    

    #else// else basker_kokkos
    #pragma omp parallel
    {}//end omp parallel
    #endif //end basker_kokkos

      }
    //-------------------End--Domian NFactor------------------------//

   
    //---------------------------Sep  NFactor--------------------------//

    
    //if(false)
    if((btf_tabs_offset != 0) && (fatal_error == BASKER_FALSE))
      {
	
	//for(Int l=1; l<=2; l++)
	for(Int l=1; l <= tree.nlvls; l++)
      {

        //Come back for syncs
        
	//#ifdef BASKER_OLD_BARRIER
	Int lthreads = pow(2,l);
	Int lnteams = num_threads/lthreads;
	//#else
	//Int lthreads = 1;
	//Int lnteams = num_threads/lthreads;
	//#endif

	Int sep_restart = 0;

	if(Options.verbose == BASKER_TRUE)
	  {
	printf("\n\n   ============ SEP: %ld ======\n\n",(long)l);
	  }

	#ifdef BASKER_KOKKOS
	//Kokkos::Timer  timer_inner_sep;
	#ifdef BASKER_NO_LAMBDA
	
	kokkos_nfactor_sep2_inc_lvl <Int, Entry, Exe_Space>
	  sep_nfactor(this,l);
	
	Kokkos::parallel_for(TeamPolicy(lnteams,lthreads),
			     sep_nfactor);
	Kokkos::fence();

	//printf("AFTER SEP \n");
       
	//======Check for error=====
	while(true)
	  {
	    INT_1DARRAY thread_start;
	    MALLOC_INT_1DARRAY(thread_start, num_threads+1);
	    init_value(thread_start, num_threads+1,
		      (Int) BASKER_MAX_IDX);
	    int nt = nfactor_sep_error(thread_start);
	    //printf("AFTER SEP ERROR %d \n", nt);
	    if((nt == BASKER_SUCCESS))
	      {
		FREE_INT_1DARRAY(thread_start);
		break;
	      }
	    else if((nt == BASKER_ERROR) ||
		    (sep_restart > BASKER_RESTART))
	      {
		FREE_INT_1DARRAY(thread_start);
		fatal_error = BASKER_TRUE;
		if(Options.verbose == BASKER_TRUE)
		  {
		printf("SEP FATAL ERROR\n");
		  }
		break;
	      }
	    else
	      {
		sep_restart++;
		if(Options.verbose == BASKER_TRUE)
		  {
		printf("sep restart l: %ld \n", (long)l);
		  }
		//exit(0);
		Kokkos::parallel_for(TeamPolicy(lnteams,lthreads),  sep_nfactor);
		Kokkos::fence();

	      }
	  }//end while-true

	
	#ifdef BASKER_TIME
	//printf("Time INNERSEP: %d %f \n", 
	//     l, timer_inner_sep.seconds());
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
       //printf("Time SEP: %f \n", timer.seconds());
    #endif
      }

   
    //-------------------------End Sep----------------//


    //-------------------IF BTF-----------------------//
    if(false)
    //if(Options.btf == BASKER_TRUE)
      {
	//=====Timer
	#ifdef BASKER_TIME
	//Kokkos::Timer  timer_btf;
	#endif
	//====Timer
	
	//======Call diag factor====
	/*
	kokkos_nfactor_diag <Int, Entry, Exe_Space> 
	  diag_nfactor(this);
	Kokkos::parallel_for(TeamPolicy(num_threads,1),
			     diag_nfactor);
	Kokkos::fence();
	*/
	//=====Check for error======
	//while(true)
	// {
	    //INT_1DARRAY thread_start;
	    // MALLOC_INT_1DARRAY(thread_start, num_threads+1);
	    //init_value(thread_start, num_threads+1, 
	    //	       (Int) BASKER_MAX_IDX);
	    //int nt = nfactor_diag_error(thread_start);
	    // if(nt == BASKER_SUCCESS)
	    //  {
	    ///		break;
	    // }
	    //else
	    // {
		/*
		break;
		printf("restart \n");
		kokkos_nfactor_diag_remalloc <Int, Entry, Exe_Space>
		  diag_nfactor_remalloc(this, thread_start);
		Kokkos::parallel_for(TeamPolicy(num_threads,1),
				     diag_nfactor);
		Kokkos::fence();
		*/
	    //}
	    // }//end while

	//====TIMER
	#ifdef BASKER_TIME
	//printf("Time BTF: %f \n", 
	//     timer_btf.seconds());
	#endif
	//===TIMER

      }//end btf call

    
    return 0;
  }//end factor_lvl_inc()
  
}//end namespace BaskerNS


#endif //End BASKER_NFACTOR_INC_HPP

