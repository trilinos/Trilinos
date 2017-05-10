#ifndef __TACHO_TEST_TASKFUNCTOR_HPP__
#define __TACHO_TEST_TASKFUNCTOR_HPP__

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "TachoExp_Util.hpp"
#include "TachoExp_TaskFunctor_MemoryPool.hpp"

using namespace Tacho::Experimental;

TEST( TaskFunctor, MemoryPool ) {
  typedef Kokkos::TaskScheduler<DeviceSpaceType> sched_type;
  typedef Kokkos::MemoryPool<DeviceSpaceType> memory_pool_type;
  
  typedef typename sched_type::memory_space memory_space;

  typedef TaskFunctor_MemoryPool_Allocate<DeviceSpaceType> functor_allocate;
  typedef TaskFunctor_MemoryPool_Deallocate<DeviceSpaceType> functor_deallocate;

  // enum { MemoryCapacity = 4000 }; // Triggers infinite loop in memory pool.
  enum { MemoryCapacity = 16000 };
  enum { MinBlockSize   =   64 };
  enum { MaxBlockSize   = 1024 };
  enum { SuperBlockSize = 1u << 12 };
  
  sched_type sched( memory_space()
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
  auto f = Kokkos::host_spawn( Kokkos::TaskSingle( sched ), functor_allocate( pool, bufsize) );
  Kokkos::host_spawn( Kokkos::TaskSingle(f), functor_deallocate(pool, f, bufsize) );
  
  Kokkos::wait( sched );
}

#endif

