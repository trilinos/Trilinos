// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLUBASKER_STATS_HPP
#define SHYLUBASKER_STATS_HPP

/*Basker Includes*/
//#include "shylubasker_decl.hpp"
#include "shylubasker_matrix_decl.hpp"
#include "shylubasker_matrix_view_decl.hpp"
#include "shylubasker_types.hpp"

/*Kokkos Includes*/
#ifdef BASKER_KOKKOS
#include <Kokkos_Core.hpp>
#else
#include <omp.h>
#endif

/*System Includes*/
#include <iostream>
#include <string>

using namespace std;

namespace BaskerNS
{
  template <class Int, class Entry, class Exe_Space>
  class BaskerStats
  {
  public:
    BaskerStats()
    {
      Unnz   = 0;
      Lnnz   = 0;
      memory = 0;
      time_setup = 0;
      time_tree = 0;
      time_sfactor = 0;
      time_nfactor = 0;
      time_upper_solve = 0;
      time_lower_solve = 0;
    }//BaskerStats()
    ~BaskerStats()
    {
      Finalize();
    }//end ~BaskerStats()
    void Finalize()
    {
      //printf("baskerstats finalize todo \n");
    }//end Finalize()
    

    Int Unnz;
    Int Lnnz;
    Int memory;
    
    double time_setup;
    double time_tree;
    double time_sfactor;
    double time_nfactor;

    double time_upper_solve;
    double time_lower_solve;


  private:

    
  };//end of Basker_stats
  

  /*
  //print stats
  template <class Int, class Entry, class Exe_Space>
  void Basker<Int,Entry, Exe_Space>::print_local_time_stats()
  {

    printf("\n\n");
    printf("-----------------------STATS----------------------");
    printf("\n\n");

    printf("---------------------Domains---------------------");
    printf("\n");

    

    for(Int i = 0; i < num_threads; i++)
      {

        printf("kid: %d Total Domain Time: %f \n", 
               i, stats.nfactor_domain_time[i]);

        printf("kid: %d Domain Solve Time: %f  Flop: %d Flops: %e \n",
               i, stats.nfactor_domain_solve_time[i], 
               stats.nfactor_domain_flop[i],
               stats.nfactor_domain_flop[i]/stats.nfactor_domain_solve_time[i]);
        
        printf("kid: %d Domain Reach Time: %f  nodes: %d  node/time: %e \n",
               i, stats.nfactor_domain_reach_time[i], 
               stats.nfactor_domain_reach_flop[i],
               stats.nfactor_domain_reach_flop[i]/stats.nfactor_domain_reach_time[i]);
           
        printf("\n");

      }

    printf("-----------------Sep---------------------------");
    printf("\n");
    for(Int i=0; i < num_threads; i++)
      {

        printf("kid: %d Total Sep Time: %f \n",
               i, stats.nfactor_sep_time[i]);

        printf("kid: %d Sep Solve Time: %f  NonAtomicFlop: %d AtomicFlop: %d Flops: %e \n", 
               i, stats.nfactor_sep_solve_time[i], 
               stats.nfactor_sep_nonatomic_flop[i], 
               stats.nfactor_sep_atomic_flop[i], 
               (stats.nfactor_sep_nonatomic_flop[i]+stats.nfactor_sep_atomic_flop[i])/stats.nfactor_sep_solve_time[i]);

        printf("kid: %d Sep Reach Time: %f nodes: %d node/time: %e \n",
               i, stats.nfactor_sep_reach_time[i],
               stats.nfactor_sep_reach_flop[i],
               stats.nfactor_sep_reach_flop[i]/stats.nfactor_sep_reach_time[i]);
        printf("\n");

      }


  }

  //returns the total number of Lnnz
  template <class Int, class Entry, class Exe_Space>
  Int Basker<Int,Entry,Exe_Space>::get_Lnnz()
  {
    if(factor_flag == -1)
      {return 0;}
  
    if(stats.Lnnz != 0)
      {return stats.Lnnz;}

    for(Int l = 0; l < tree.nblks; l++)
      {
        MATRIX &myL = LL[l][0];
        stats.Lnnz += LL[l][0].nnz;
      }//over all Ls

    return stats.Lnnz;
  }//end get_Lnnz();

  template <class Int, class Entry, class Exe_Space>
  Int Basker<Int,Entry,Exe_Space>::get_Unnz()
  {
    if(factor_flag == false)
      {return 0;}

    if(stats.Unnz != 0)
      {return stats.Unnz;}

    for(Int l = 0; l < tree.nblks; l++)
      {
        for(Int r=0; r<LU[l].size(); r++)
          {
            Int k = LU[l][r].ncol;
            stats.Unnz += LU[l][r].col_ptr[k];
            //know problem with .nnz 9too much... need to fix
          }
      }
    return stats.Unnz;
  }//end get_Unnz()


  */
  //To be continued....

}//end namespace basker

#endif //end IFNDEF basker_stats_hpp
