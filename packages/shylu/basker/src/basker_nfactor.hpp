#ifndef BASKER_NFACTOR_HPP
#define BASKER_NFACTOR_HPP

#define BASKER_DEBUG_NFACTOR 

#define BASKER_TIME

/*Basker Includes*/
//#include "basker_decl.hpp"
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

//using namespace std;

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

	//printf("Now do after btf stuff \n");
	Kokkos::Impl::Timer  timer_btf;
	kokkos_nfactor_diag <Int, Entry, Exe_Space> 
	  diag_nfactor(this);
	Kokkos::parallel_for(TeamPolicy(num_threads,1),
			     diag_nfactor);
	Kokkos::fence();
	#ifdef BASKER_TIME
	printf("Time BTF: %f \n", 
	       timer_btf.seconds());
	#endif

      }//end btf call

    
    return 0;
  }//end factor_notoken()
  

  /*
  template <class Int, class Entry, class Exe_Space>
  BASKER_INLINE
  int Basker<Int, Entry, Exe_Space>::factor_token(Int option)
  {
    #ifdef BASKER_KOKKOS
    //Note: May want to turn the parallel_fors into parallel_reduces ???
    typedef Kokkos::TeamPolicy<Exe_Space>       TeamPolicy;
    typedef typename TeamPolicy::member_type    MemberType;
    #ifdef BASKER_NO_LAMBDA
    kokkos_nfactor_funct<Int,Entry,Exe_Space>  
      nfactor(this);
    Kokkos::parallel_for(TeamPolicy(1,num_threads), nfactor);
    Kokkos::fence();
    #else
    //Note: comebask
    //Kokkos::parallel_for(TeamPolicy(1,num_threads),);
    #endif //end BASKER_NO_LAMBDA
    #else //BASKER_KOKKOS
    #pragma omp parallel
    {
      //call work here
    }
    #endif //BASKER_KOKKOS
    return 0;
  }//end factor()
  */

  /*
  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry,Exe_Space>::t_factor_tree(Int kid)
  {
    //Init workspace
    t_init_workspace(kid);
    
    for(Int l = 0; l < tree.nlvls+1; l++) //all lvls in factor tree
      {
	if(l == 0)
	  {
	    t_nfactor_blk(kid);
	  }
	else
	  {
	    //sep calls

	  }
      }//end over all lvls in tree
  }//end t_factor_tree(Int kid)
  */

}//end namespace baskerNS
#endif //end ifndef basker_nfactor_hpp
