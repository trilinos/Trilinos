#ifndef __TACHO_TEST_TASKFUNCTOR_HPP__
#define __TACHO_TEST_TASKFUNCTOR_HPP__

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "Tacho_Util.hpp"
#include "Tacho_TaskFunctor_MemoryPool.hpp"

using namespace Tacho;

TEST( TaskFunctor, MemoryPool ) {
  TEST_BEGIN;
  typedef Kokkos::TaskScheduler<typename HostDeviceType::execution_space> scheduler_type;
  typedef Kokkos::MemoryPool<HostDeviceType> memory_pool_type;
  
  typedef typename HostDeviceType::memory_space memory_space;

  typedef TaskFunctor_MemoryPool_Allocate<scheduler_type> functor_allocate;
  typedef TaskFunctor_MemoryPool_Deallocate<scheduler_type> functor_deallocate;
  typedef TaskFunctor_MemoryPool_TestView<scheduler_type> functor_test_view;
  typedef TaskFunctor_MemoryPool_TestViewSee<scheduler_type> functor_test_view_see;

  // enum { MemoryCapacity = 4000 }; // Triggers infinite loop in memory pool.
  enum { MemoryCapacity = 16000 };
  enum { MinBlockSize   =   64 };
  enum { MaxBlockSize   = 1024 };
  enum { SuperBlockSize = 1u << 12 };
  
  scheduler_type sched( typename scheduler_type::memory_space()
                    , MemoryCapacity
                    , MinBlockSize
                    , MaxBlockSize
                    , SuperBlockSize );

  memory_pool_type pool(memory_space()
                    , MemoryCapacity
                    , MinBlockSize
                    , MaxBlockSize
                    , SuperBlockSize );
  
  const ordinal_type bufsize = 10*10*sizeof(double);
  auto f_alloc    = Kokkos::host_spawn(Kokkos::TaskSingle(sched), functor_allocate(pool, bufsize));
  auto f_view     = Kokkos::host_spawn(Kokkos::TaskSingle(f_alloc), functor_test_view(pool, f_alloc, 10, 10));
  auto f_view_see = Kokkos::host_spawn(Kokkos::TaskSingle(f_view), functor_test_view_see(pool, f_view));
  Kokkos::host_spawn( Kokkos::TaskSingle(f_view_see), functor_deallocate(pool, f_alloc, bufsize) );

  Kokkos::wait( sched );  

  TEST_END;
}

#endif

