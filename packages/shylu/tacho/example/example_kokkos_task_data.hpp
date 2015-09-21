#pragma once
#ifndef __EXAMPLE_KOKKOS_TASK_DATA_HPP__
#define __EXAMPLE_KOKKOS_TASK_DATA_HPP__

#include <iostream>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "ShyLUTacho_config.h"

#include "util.hpp"

#include <Kokkos_Threads.hpp>
#include <Threads/Kokkos_Threads_TaskPolicy.hpp>

#define MAXTASKS 100000

#define BIG        100
#define SMALL      10000
#define TINY         1

namespace Tacho {

  using namespace std;

  template<typename PolicyType>
  class SimpleTask {
  public:
    SimpleTask(Kokkos::View<double*,typename PolicyType::execution_space> &a,
               Kokkos::View<double*,typename PolicyType::execution_space> &b,
               const long position = 0) 
      : _a(a), 
        _b(b),
        _position(position) { }

    typedef int value_type; 
    typedef typename PolicyType::member_type member_type;

    const Kokkos::View<double*,typename PolicyType::execution_space,Kokkos::MemoryUnmanaged> _a;
    const Kokkos::View<double*,typename PolicyType::execution_space,Kokkos::MemoryUnmanaged> _b;
    const long _position;

    // task only interface 
    KOKKOS_INLINE_FUNCTION
    void apply(value_type &r_val) {

      for (int iter=0;iter<BIG;++iter) {
        for (long i=0;i<SMALL;++i) {
          double tmp = 0.0;
          for (long j=0;j<TINY;++j) 
            tmp += j;
          _a[i+_position*SMALL] = _b[i+_position*SMALL] + (tmp + 1);
        }
      }
      r_val = 0;
    }
    
    // task team interface
    KOKKOS_INLINE_FUNCTION
    void apply(const member_type &member, value_type &r_val) {
      for (int iter=0;iter<BIG;++iter) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, SMALL),
                             [&](const long i) {
                               double tmp = 0.0;
                               for (long j=0;j<TINY;++j)
                                 tmp += j;
                               _a[i+_position*SMALL] = _b[i+_position*SMALL] + (tmp + 1);
                             });
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
    
    policy_type policy(4, team_size);

    Kokkos::Impl::Timer timer;
    Kokkos::View<double*,SpaceType> a("TaskDataExample::a", ntasks*SMALL);
    Kokkos::View<double*,SpaceType> b("TaskDataExample::b", ntasks*SMALL);

    {
      //timer.reset();
      
      for (int i=0;i<10;++i) {
        future_type f = policy.create(simple_task_type(a,b), max_task_dependence);
        policy.spawn(f);
        
        Kokkos::Experimental::wait( policy );
      }
      //const double t = timer.seconds();
      //cout << "KokkosTaskData:: single task pre-run :: time = " << t << endl;
    }
    {
      timer.reset();
      
      future_type f = policy.create(simple_task_type(a,b), max_task_dependence);
      policy.spawn(f);
      
      Kokkos::Experimental::wait( policy );
      
      const double t = timer.seconds();
      cout << "KokkosTaskData:: single task is spawned :: time = " << t << endl;
    }
    {
      timer.reset();
      
      future_type f = policy.create_team(simple_task_type(a,b), max_task_dependence);
      policy.spawn(f);
      
      Kokkos::Experimental::wait( policy );
      
      const double t = timer.seconds();
      cout << "KokkosTaskData:: single team task is spawned :: time = " << t << endl;
    }
    {
      timer.reset();
      
      future_type f[MAXTASKS];
      for (int i=0;i<ntasks;++i) {
        f[i] = policy.create(simple_task_type(a,b,i), max_task_dependence);
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
        f[i] = policy.create_team(simple_task_type(a,b,i), max_task_dependence);
        policy.spawn(f[i]);
      }
      
      Kokkos::Experimental::wait( policy );
      
      const double t = timer.seconds();
      cout << "KokkosTaskData:: " << ntasks << " team tasks are spawned :: time = " << t << endl;
    }
    return 0;
  }
}

#endif
