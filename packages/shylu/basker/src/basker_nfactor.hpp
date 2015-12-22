#ifndef BASKER_NFACTOR_HPP
#define BASKER_NFACTOR_HPP

//#define BASKER_DEBUG_NFACTOR 

#define BASKER_TIME

/*Basker Includes*/
#include "basker_types.hpp"
#include "basker_util.hpp"
#include "basker_structs.hpp"
#include "basker_matrix_view_decl.hpp"
#include "basker_matrix_view_def.hpp"

#include "basker_nfactor_blk.hpp"
#include "basker_nfactor_col.hpp"
#include "basker_nfactor_col2.hpp"
#include "basker_nfactor_diag.hpp"

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

    if(Options.btf == BASKER_TRUE)
      {

	//JDB: We can change this for the new inteface

	//call reference copy constructor
	gn = A.ncol;
	gm = A.nrow;
	A = BTF_A; 
	//printf("\n\n Switching A, newsize: %d \n",
        //     A.ncol);
      }
   

    //Spit into Domain and Sep
    //----------------------Domain-------------------------//
    #ifdef BASKER_KOKKOS
    Kokkos::Impl::Timer       timer;


    typedef Kokkos::TeamPolicy<Exe_Space>        TeamPolicy;

    if(btf_tabs_offset != 0)
      {
    #ifdef BASKER_NO_LAMBDA
    printf("Nfactor, NoToken, Kokkos, NoLambda \n");
    kokkos_nfactor_domain <Int,Entry,Exe_Space>
      domain_nfactor(this);
    Kokkos::parallel_for(TeamPolicy(num_threads,1), domain_nfactor);
    Kokkos::fence();
    #else //else no_lambda
    typedef typename TeamPolicy::member_type     MemberType;
    //Note: To be added
    #endif //end basker_no_lambda
    
    #ifdef BASKER_TIME
    printf("Time DOMAIN: %f \n", timer.seconds());
    #endif
    timer.reset();

    #else// else basker_kokkos
    #pragma omp parallel
    {


    }//end omp parallel
    #endif //end basker_kokkos

      }
    //-------------------End--Domian--------------------------//

   
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

	#ifdef BASKER_KOKKOS
	Kokkos::Impl::Timer  timer_inner_sep;
	#ifdef BASKER_NO_LAMBDA
	/*
	kokkos_nfactor_sep <Int, Entry, Exe_Space> 
	  sep_nfactor(this, l);
	*/
	kokkos_nfactor_sep2 <Int, Entry, Exe_Space>
	  sep_nfactor(this,l);
	
	Kokkos::parallel_for(TeamPolicy(lnteams,lthreads),
			     sep_nfactor);
	Kokkos::fence();
	#ifdef BASKER_TIME
	printf("Time INNERSEP: %d %f \n", 
	       l, timer_inner_sep.seconds());
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
    printf("Time SEP: %f \n", timer.seconds());
    #endif
      }

    //-------------------------End Sep----------------//


    //-------------------IF BTF-----------------------//
    if(Options.btf == BASKER_TRUE)
      {
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
	    if(nt == BASKER_SUCCESS)
	      {
		break;
	      }
	    else
	      {
		printf("restart \n");
		kokkos_nfactor_diag_remalloc <Int, Entry, Exe_Space>
		  diag_nfactor_remalloc(this, thread_start);
		Kokkos::parallel_for(TeamPolicy(num_threads,1),
				     diag_nfactor);
		Kokkos::fence();
	      }
	  }//end while

	//====TIMER
	#ifdef BASKER_TIME
	printf("Time BTF: %f \n", 
	       timer_btf.seconds());
	#endif
	//===TIMER

      }//end btf call

    
    return 0;
  }//end factor_notoken()

  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int,Entry,Exe_Space>::nfactor_diag_error
  (
   INT_1DARRAY threads_start
   )
  {
    Int nthread_remalloc = 0;
    for(Int ti = 0; ti < num_threads; ti++)
      {

	//Note: jdb we can make this into a switch
	if(thread_array(ti).error_type ==
	   BASKER_ERROR_NOERROR)
	  {
	    threads_start(ti) = BASKER_MAX_IDX;
	    continue;
	  }//end if NOERROR

	if(thread_array(ti).error_type ==
	   BASKER_ERROR_SINGULAR)
	  {
	    printf("ERROR THREAD: %d DIAGBLK SINGULAR: %d\n",
		   ti,
		   thread_array(ti).error_blk);
	    return BASKER_ERROR;
	  }//end if SINGULAR
	
	if(thread_array(ti).error_type ==
	   BASKER_ERROR_NOMALLOC)
	  {
	    printf("ERROR THREADS: %d DIAGBLK NOMALLOC: %d\n",
		   ti,
		   thread_array(ti).error_blk);
	    return BASKER_ERROR;
	  }//end if NOMALLOC
	
	if(thread_array(ti).error_type ==
	   BASKER_ERROR_REMALLOC)
	  {

	    BASKER_ASSERT(thread_array(ti).error_blk > 0,
			  "nfactor_diag_error error_blk");

	    printf("ERROR THREADS: %d DIAGBLK MALLOC: %d \n",
		   ti,
		   thread_array(ti).error_blk);
	    

	    //Resize L
	    BASKER_MATRIX &L = LBTF(thread_array(ti).error_blk);
	    RESIZE_INT_1DARRAY(L.row_idx,
			       L.nnz,
			       thread_array(ti).error_info);
	    RESIZE_ENTRY_1DARRAY(L.val,
				 L.nnz,
				 thread_array(ti).error_info);
	    L.nnz = thread_array(ti).error_info;
	    
	    //Resize U
	    BASKER_MATRIX &U = UBTF(thread_array(ti).error_blk);
	    RESIZE_INT_1DARRAY(U.row_idx,
			       U.nnz,
			       thread_array(ti).error_info);
	    RESIZE_ENTRY_1DARRAY(U.val,
				 U.nnz,
				 thread_array(ti).error_info);
	    U.nnz = thread_array(ti).error_info;
	    
	    threads_start(ti) = thread_array(ti).error_blk;
	

	    //Reset 
	    thread_array(ti).error_type = BASKER_ERROR_NOERROR;
	    thread_array(ti).error_blk  = BASKER_MAX_IDX;
	    thread_array(ti).error_info = BASKER_MAX_IDX;

	    nthread_remalloc++;
    
	  }//if REMALLOC

      }//for all threads
    
    if(nthread_remalloc == 0)
      {
	return BASKER_SUCCESS;
      }
    else
      {
	return nthread_remalloc;
      }

    //Should never be here
    BASKER_ASSERT(0==1, "nfactor_diag_error, should never");
    return BASKER_SUCCESS;
  }//end nfactor_diag_error
  

  
}//end namespace baskerNS
#endif //end ifndef basker_nfactor_hpp
