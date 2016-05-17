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

  template<typename PolicyType, typename ValueType>
  class SimpleTask {
  public:
    SimpleTask(Kokkos::View<ValueType*,typename PolicyType::execution_space> &a,
               Kokkos::View<ValueType*,typename PolicyType::execution_space> &b,
               const int itask = 0) 
      : _a(a), _b(b), _itask(itask) { }

    typedef int value_type; 
    typedef typename PolicyType::member_type member_type;

    const Kokkos::View<ValueType*,typename PolicyType::execution_space,Kokkos::MemoryUnmanaged> _a;
    const Kokkos::View<ValueType*,typename PolicyType::execution_space,Kokkos::MemoryUnmanaged> _b;
    const int _itask;

    // task only interface 
    KOKKOS_INLINE_FUNCTION
    void apply(value_type &r_val) {
      const int offset = _itask*SMALL;
      auto a = &_a[offset];
      auto b = &_b[offset];

      for (int iter=0;iter<BIG;++iter) {
        for (int i=0;i<SMALL;++i) {
          ValueType tmp = 0.0;
          for (int j=0;j<TINY;++j) 
            tmp += j;
          a[i] = b[i] + (tmp + 1);
        }
      }
      r_val = 0;
    }
    
    // task team interface
    KOKKOS_INLINE_FUNCTION
    void apply(const member_type &member, value_type &r_val) {
      const int offset = _itask*SMALL;
      auto a = &_a[offset];
      auto b = &_b[offset];

      for (int iter=0;iter<BIG;++iter) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, SMALL),
                             [&](const int i) {
                               ValueType tmp = 0.0;
                               for (int  j=0;j<TINY;++j)
                                 tmp += j;
                               a[i] = b[i] + (tmp + 1);
                             });
        // For testing
        member.team_barrier();
      }
    }
  };

  template<typename SpaceType, typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int exampleKokkosTaskData(const int ntasks,
                            const int max_task_dependence,
                            const int team_size,
                            const bool verbose) {
    typedef SpaceType space_type;
    typedef ValueType value_type;

    typedef Kokkos::View<value_type*,space_type> array_type;

    typedef Kokkos::Experimental::TaskPolicy<space_type> policy_type ;
    typedef SimpleTask<policy_type,value_type> simple_task_type;

    typedef Kokkos::Experimental::Future<typename simple_task_type::value_type,space_type> future_type ;

    Kokkos::Impl::Timer timer;    

    policy_type policy(ntasks,
                       sizeof(simple_task_type),
                       max_task_dependence, 
                       team_size);

    const int nsize = ntasks*SMALL;
    array_type a("TaskDataExample::a", nsize);
    array_type b("TaskDataExample::b", nsize);
    array_type c("TaskDataExample::c", nsize);

    value_type aa = 0.0, bb = 10.0;

    {
      Kokkos::Impl::ViewFill<array_type>(a, aa);
      Kokkos::Impl::ViewFill<array_type>(b, bb);
    }
    {
      future_type f[MAXTASKS];
      for (int i=0;i<ntasks;++i) {
        f[i] = policy.create(simple_task_type(a,b,i), max_task_dependence);
        policy.spawn(f[i]);
      }
      Kokkos::Experimental::wait( policy );

      Kokkos::deep_copy(c, a);
    }
    {
      Kokkos::Impl::ViewFill<array_type>(a, aa);
      Kokkos::Impl::ViewFill<array_type>(b, bb);
    }
    {
      timer.reset();
      
      future_type f = policy.create(simple_task_type(a,b), max_task_dependence);
      policy.spawn(f);
      
      Kokkos::Experimental::wait( policy );
      
      const double t = timer.seconds();

      value_type diff = 0.0;
      for (int i=0;i<SMALL;++i)
        diff += abs(a[i]-c[i]);
      cout << "KokkosTaskData:: single task is spawned :: time = " << t << "  diff = " << diff << endl;
    }
    {
      Kokkos::Impl::ViewFill<array_type>(a, aa);
      Kokkos::Impl::ViewFill<array_type>(b, bb);
    }
    {
      timer.reset();
      
      future_type f = policy.create_team(simple_task_type(a,b), max_task_dependence);
      policy.spawn(f);
      
      Kokkos::Experimental::wait( policy );
      
      const double t = timer.seconds();

      value_type diff = 0.0;
      for (int i=0;i<SMALL;++i) 
        diff += abs(a[i]-c[i]);
      cout << "KokkosTaskData:: single team task is spawned :: time = " << t << "  diff = " << diff << endl;
    }
    {
      Kokkos::Impl::ViewFill<array_type>(a, aa);
      Kokkos::Impl::ViewFill<array_type>(b, bb);
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

      value_type diff = 0.0;
      for (int i=0;i<SMALL;++i)
        diff += abs(a[i]-c[i]);
      cout << "KokkosTaskData:: " << ntasks << " tasks are spawned :: time = " << t << "  diff = " << diff << endl;
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

      value_type diff = 0.0;
      for (int i=0;i<SMALL;++i)
        diff += abs(a[i]-c[i]);
      cout << "KokkosTaskData:: " << ntasks << " team tasks are spawned :: time = " << t << "  diff = " << diff << endl;
    }
    {
      try {
        // For Carter
        for (int i=0;i<ntasks;++i) {
          future_type f = policy.create(simple_task_type(a,b,i), max_task_dependence);
          policy.spawn(f);
        }
        Kokkos::Experimental::wait( policy );
      } catch (exception &e) {
        cout << "Exception is caught " << endl << e.what() << endl;
      }
    }

    return 0;
  }
}

#endif
