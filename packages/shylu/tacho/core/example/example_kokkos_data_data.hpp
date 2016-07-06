#pragma once
#ifndef __EXAMPLE_KOKKOS_DATA_DATA_HPP__
#define __EXAMPLE_KOKKOS_DATA_DATA_HPP__

#include <iostream>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "ShyLUTacho_config.h"

#include "util.hpp"

#include <Kokkos_Threads.hpp>

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
               Kokkos::View<ValueType*,typename PolicyType::execution_space> &b)
      : _a(a), _b(b) { }
    const Kokkos::View<ValueType*,typename PolicyType::execution_space,Kokkos::MemoryUnmanaged> _a;
    const Kokkos::View<ValueType*,typename PolicyType::execution_space,Kokkos::MemoryUnmanaged> _b;

    void operator()(const typename PolicyType::member_type i) const {
      for (int iter=0;iter<BIG;++iter) {
        ValueType tmp = 0.0;
        for (int j=0;j<TINY;++j) 
          tmp += j;
        _a[i] = _b[i] + (tmp + 1);
      }
    }
  };

  template<typename PolicyType, typename ValueType>
  class SimpleTaskTeam {
  public:
    SimpleTaskTeam(Kokkos::View<ValueType*,typename PolicyType::execution_space> &a,
                   Kokkos::View<ValueType*,typename PolicyType::execution_space> &b,
                   int ntasks)
      : _a(a), _b(b), _ntasks(ntasks) { }
    const Kokkos::View<ValueType*,typename PolicyType::execution_space,Kokkos::MemoryUnmanaged> _a;
    const Kokkos::View<ValueType*,typename PolicyType::execution_space,Kokkos::MemoryUnmanaged> _b;
    int _ntasks;

    void operator()(const typename PolicyType::member_type &member) const {
      const int lrank = member.league_rank();
      const int lsize = member.league_size();

      for (int k=lrank;k<_ntasks;k+=lsize) {
        const int offset = k*SMALL;        
        auto a = &_a[offset];
        auto b = &_b[offset];

        for (int iter=0;iter<BIG;++iter) {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(member, SMALL), [&](const int i) {
              ValueType tmp = 0.0;
              for (int j=0;j<TINY;++j) 
                tmp += j;
              a[i] = b[i] + (tmp + 1);
            });
        }
      }
    }
  };
  
  template<typename SpaceType, typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int exampleKokkosDataData(const int ntasks,
                            const int league_size,
                            const int team_size,
                            const bool verbose) {
    typedef SpaceType space_type;
    typedef ValueType value_type;

    typedef Kokkos::View<value_type*,space_type> array_type;

    Kokkos::Impl::Timer timer;
    
    const int nsize = ntasks*SMALL;
    array_type a("DataDataExample::a", nsize);
    array_type b("DataDataExample::b", nsize);
    array_type c("DataDataExample::c", nsize);

    value_type aa = 0.0, bb = 10.0;

    {
      Kokkos::Impl::ViewFill<array_type>(a, aa);
      Kokkos::Impl::ViewFill<array_type>(b, bb);
    }
    {
      typedef Kokkos::RangePolicy<space_type> policy_type;
      typedef SimpleTask<policy_type,value_type> simple_task_type;
      
      Kokkos::RangePolicy<space_type> policy(0,nsize);
      Kokkos::parallel_for( policy, simple_task_type(a,b) );
      
      Kokkos::deep_copy(c, a);
    }
    {
      Kokkos::Impl::ViewFill<array_type>(a, aa);
      Kokkos::Impl::ViewFill<array_type>(b, bb);
    }
    {
      typedef Kokkos::RangePolicy<space_type> policy_type;
      typedef SimpleTask<policy_type,value_type> simple_task_type;
      
      timer.reset();

      policy_type policy(0, SMALL);
      Kokkos::parallel_for( policy, simple_task_type(a,b) );

      const double t = timer.seconds();

      value_type diff = 0.0;
      for (int i=0;i<SMALL;++i) 
        diff += abs(a[i]-c[i]);
      cout << "KokkosDataData:: single task is spawned :: time = " << t << "  diff = " << diff << endl;
    }
    {
      Kokkos::Impl::ViewFill<array_type>(a, aa);
      Kokkos::Impl::ViewFill<array_type>(b, bb);
    }
    {
      typedef Kokkos::TeamPolicy<space_type> policy_type;
      typedef SimpleTaskTeam<policy_type,value_type> simple_task_type;

      const int league_size_single = 1;

      timer.reset();

      policy_type policy(league_size_single, team_size);
      Kokkos::parallel_for( policy, simple_task_type(a,b,1) );
        
      const double t = timer.seconds();

      value_type diff = 0.0;
      for (int i=0;i<SMALL;++i) 
        diff += abs(a[i]-c[i]);
      cout << "KokkosDataData:: single team task is spawned :: time = " << t << "  diff = " << diff << endl;
    }
    {
      Kokkos::Impl::ViewFill<array_type>(a, aa);
      Kokkos::Impl::ViewFill<array_type>(b, bb);
    }
    {
      typedef Kokkos::RangePolicy<space_type> policy_type;
      typedef SimpleTask<policy_type,value_type> simple_task_type;
      
      timer.reset();

      policy_type policy(0, SMALL*ntasks);
      Kokkos::parallel_for( policy, simple_task_type(a,b) );
      
      const double t = timer.seconds();

      value_type diff = 0.0;
      for (int i=0;i<nsize;++i) 
        diff += abs(a[i]-c[i]);
      cout << "KokkosDataData:: " << ntasks << " tasks are spawned :: time = " << t << "  diff = " << diff << endl;
    }
    {
      Kokkos::Impl::ViewFill<array_type>(a, aa);
      Kokkos::Impl::ViewFill<array_type>(b, bb);
    }
    {
      typedef Kokkos::TeamPolicy<space_type> policy_type;
      typedef SimpleTaskTeam<policy_type,value_type> simple_task_type;

      timer.reset();

      policy_type policy(league_size, team_size);
      Kokkos::parallel_for( policy, simple_task_type(a,b,ntasks) );
      
      const double t = timer.seconds();
      value_type diff = 0.0;
      for (int i=0;i<nsize;++i) 
        diff += abs(a[i]-c[i]);
      cout << "KokkosDataData:: " << ntasks << " team tasks are spawned :: time = " << t << "  diff = " << diff << endl;
    }
    return 0;
  }
}

#endif
