#ifndef __TACHO_TEST_TASKFUNCTOR_HPP__
#define __TACHO_TEST_TASKFUNCTOR_HPP__

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "TachoExp_Util.hpp"
#include "TachoExp_TaskFunctor_MemoryPool.hpp"

using namespace Tacho::Experimental;

TEST( TaskFunctor, MemoryPool ) {

  typedef Kokkos::TaskScheduler<HostSpaceType> sched_type;
  typedef TaskFunctor_MemoryPool_Allocate<HostSpaceType> functor_allocate;
  typedef TaskFunctor_MemoryPool_Deallocate<HostSpaceType> functor_deallocate;
  typedef TaskFunctor_MemoryPool_TestView<HostSpaceType> functor_testview;
  typedef TaskFunctor_MemoryPool_TestViewSee<HostSpaceType> functor_testview_see;

  typedef typename functor_allocate::memory_pool_type memory_pool_type;
  typedef typename sched_type::memory_space memory_space;

  const size_type task_queue_capacity = 1000*sizeof(double);
  const ordinal_type
    min_block_alloc_size = 32,
    max_block_alloc_size = 512,
    min_superblock_size = 1024;

#if defined (__KK__)
  sched_type sched(memory_space(),
                   task_queue_capacity,
                   min_block_alloc_size,
                   max_block_alloc_size,
                   min_superblock_size);
#else
  sched_type sched(memory_space(),
                   task_queue_capacity);
#endif


  const size_type min_total_alloc_size = 4096*10;
#if defined (__KK__)
  memory_pool_type pool(memory_space(),
                        min_total_alloc_size,
                        min_block_alloc_size,
                        max_block_alloc_size,
                        min_superblock_size);
#else
  memory_pool_type pool(memory_space(),
                        min_total_alloc_size);
#endif

  const ordinal_type bufsize = 10*10*sizeof(double);
  auto f_allocate = Kokkos::host_spawn(Kokkos::TaskSingle(sched),
                                       functor_allocate(pool, bufsize));
 
  auto f_view = Kokkos::host_spawn(Kokkos::TaskSingle(sched),// Kokkos::TaskSingle(f_allocate),
                                   functor_testview(pool, 10, 10, f_allocate));

  auto f_view_see = Kokkos::host_spawn(Kokkos::TaskSingle(sched),//Kokkos::TaskSingle(f_view),
                                       functor_testview_see(pool, f_view));

  // auto f_deallocate = Kokkos::host_spawn(Kokkos::TaskSingle(sched),//Kokkos::TaskSingle(f_view_see),
  //                                        functor_deallocate(pool, bufsize, f_allocate));

  Kokkos::wait(sched);
}

#endif

