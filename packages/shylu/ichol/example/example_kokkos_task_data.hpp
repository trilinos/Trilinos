#pragma once
#ifndef __EXAMPLE_KOKKOS_TASK_DATA_HPP__
#define __EXAMPLE_KOKKOS_TASK_DATA_HPP__

#include <iostream>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "ShyLUIChol_config.h"

#include "util.hpp"

#include <Kokkos_Threads.hpp>
#include <Threads/Kokkos_Threads_TaskPolicy.hpp>

#define MAXTASKS 100000

#define BIG        1000
#define SMALL      1000
#define TINY         10

namespace Example {

  using namespace std;

  template<typename PolicyType>
  class SimpleTask {
  public:
    SimpleTask(const bool use_barrier) : _use_barrier(use_barrier) { }

    typedef int value_type; 
    typedef typename PolicyType::member_type member_type;

    bool _use_barrier;
    double _dummy[SMALL];
    
    // task only interface 
    void apply(value_type &r_val) {

      for (int iter=0;iter<BIG;++iter) {
        for (long i=0;i<SMALL;++i) {
          double tmp = 0.0;
          for (long j=0;j<TINY;++j) {
            tmp += j;
          }
          _dummy[i] += (tmp + 1);
        }
      }
      r_val = 0;
    }
    
    // task team interface
    void apply(const member_type &member, value_type &r_val) {
      for (int iter=0;iter<BIG;++iter) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, SMALL),
                             [&](const long i) {
                               double tmp = 0.0;
                               for (long j=0;j<TINY;++j)
                                 tmp += j;
                               _dummy[i] += (tmp + 1);
                             });
        if (_use_barrier)
          member.team_barrier();
      }
    }
  };

  template<typename SpaceType>
  KOKKOS_INLINE_FUNCTION
  int exampleKokkosTaskData(const int ntasks,
                            const int max_task_dependence,
                            const int team_size, 
                            const bool verbose) {

    typedef Kokkos::Experimental::TaskPolicy<SpaceType> policy_type ;
    typedef SimpleTask<policy_type> simple_task_type;

    typedef Kokkos::Experimental::Future<typename simple_task_type::value_type,SpaceType> future_type ;
    
    policy_type policy;

    Kokkos::Impl::Timer timer;

    for (int use_barrier=0;use_barrier<2;++use_barrier) {
      cout << "KokkosTaskData:: use barrier " << (use_barrier ? "yes" : "no") << endl;
      {
        timer.reset();
        
        future_type f = policy.create(simple_task_type(use_barrier), max_task_dependence);
        policy.spawn(f);
        
        Kokkos::Experimental::wait( policy );
        
        const double t = timer.seconds();
        cout << "KokkosTaskData:: single task is spawned :: time = " << t << endl;
      }
      
      {
        timer.reset();
        
        future_type f = policy.create_team(simple_task_type(use_barrier), max_task_dependence);
        policy.spawn(f);
        
        Kokkos::Experimental::wait( policy );
        
        const double t = timer.seconds();
        cout << "KokkosTaskData:: single team task is spawned :: time = " << t << endl;
      }
      
      {
        timer.reset();
        
        future_type f[MAXTASKS];
        for (int i=0;i<ntasks;++i) {
          f[i] = policy.create(simple_task_type(use_barrier), max_task_dependence);
          policy.spawn(f[i]);
        }
        Kokkos::Experimental::wait( policy );
      
        const double t = timer.seconds();
        cout << "KokkosTaskData:: " << ntasks << " tasks are spawned :: time = " << t << endl;
      }
      
      {
        timer.reset();
        
        future_type f[MAXTASKS];
        for (int i=0;i<ntasks;++i) {
          f[i] = policy.create_team(simple_task_type(use_barrier), max_task_dependence);
          policy.spawn(f[i]);
        }
        
        Kokkos::Experimental::wait( policy );
        
        const double t = timer.seconds();
        cout << "KokkosTaskData:: " << ntasks << " team tasks are spawned :: time = " << t << endl;
      }
    }
    return 0;
  }
}

#endif
