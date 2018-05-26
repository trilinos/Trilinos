#ifndef __TACHO_TEST_TASKFUNCTOR_HPP__
#define __TACHO_TEST_TASKFUNCTOR_HPP__

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>

#include "TachoExp_Util.hpp"
#include "TachoExp_TaskFunctor_MemoryPool.hpp"

using namespace Tacho::Experimental;




template< class Space, typename future_type >
struct TestTaskSpawn {
  typedef Kokkos::TaskScheduler< Space >  sched_type;
  typedef int                            value_type;

  sched_type   m_sched ;
  future_type  m_future ;

  KOKKOS_INLINE_FUNCTION
  TestTaskSpawn( const sched_type & arg_sched
                 , const future_type & arg_future
                 )
    : m_sched( arg_sched )
    , m_future( arg_future )
  {}

  KOKKOS_INLINE_FUNCTION
  void operator()( typename sched_type::member_type &, value_type &r_val )
  {
    if ( ! m_future.is_null() ) {
      printf("I am task spawn\n");
      Kokkos::task_spawn( Kokkos::TaskSingle( m_sched ) , TestTaskSpawn( m_sched , future_type() ) );
      r_val = 0;
    } else {
      printf("I have future and I won't spawn\n");
      r_val = 1;
    }

  }
};

TEST( TaskFunctor, MemoryPool ) {
  typedef Kokkos::TaskScheduler<HostSpaceType> sched_type;
  typedef Kokkos::MemoryPool<HostSpaceType> memory_pool_type;
  typedef TaskFunctor_MemoryPool_Allocate<HostSpaceType> functor_allocate;
  typedef TaskFunctor_MemoryPool_Deallocate<HostSpaceType> functor_deallocate;

  typedef Kokkos::Future< int, HostSpaceType >         future_type;
  typedef typename sched_type::memory_space memory_space;

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
  Kokkos::host_spawn( Kokkos::TaskSingle( f ), functor_deallocate( pool, bufsize) );
  
  Kokkos::wait( sched );






//   typedef TaskFunctor_MemoryPool_TestView<HostSpaceType> functor_testview;
//   typedef TaskFunctor_MemoryPool_TestViewSee<HostSpaceType> functor_testview_see;

//   typedef typename functor_allocate::memory_pool_type memory_pool_type;
//   typedef typename sched_type::memory_space memory_space;

//   const size_type task_queue_span = 1000*sizeof(double);
//   const ordinal_type
//     min_block_alloc_size = 32,
//     max_block_alloc_size = 512,
//     min_superblock_size = 1024;

// #if defined (__KK__)
//   sched_type sched(memory_space(),
//                    task_queue_span,
//                    min_block_alloc_size,
//                    max_block_alloc_size,
//                    min_superblock_size);
// #else
//   sched_type sched(memory_space(),
//                    task_queue_span);
// #endif


//   const size_type min_total_alloc_size = 4096*10;
// #if defined (__KK__)
//   memory_pool_type pool(memory_space(),
//                         min_total_alloc_size,
//                         min_block_alloc_size,
//                         max_block_alloc_size,
//                         min_superblock_size);
// #else
//   memory_pool_type pool(memory_space(),
//                         min_total_alloc_size);
// #endif

//   const ordinal_type bufsize = 10*10*sizeof(double);
// #if defined (__KK__)
//   auto f_allocate = Kokkos::host_spawn(Kokkos::TaskSingle(sched), 
//                                        functor_allocate(pool, bufsize));
//   Kokkos::host_spawn(Kokkos::TaskSingle(f_allocate, Kokkos::TaskPriority::High), 
//                      functor_deallocate(pool, bufsize, f_allocate));
  
//   // auto f_allocate = Kokkos::host_spawn(Kokkos::TaskSingle(sched),
//   //                                      functor_allocate(pool, bufsize));
 
//   // // auto f_view = Kokkos::host_spawn(Kokkos::TaskSingle(f_allocate),
//   // //                                  functor_testview(pool, 10, 10, f_allocate));

//   // // auto f_view_see = Kokkos::host_spawn(Kokkos::TaskSingle(f_view),
//   // //                                      functor_testview_see(pool, f_view));

//   // auto f_deallocate = Kokkos::host_spawn(Kokkos::TaskSingle(f_allocate),
//   //                                        functor_deallocate(pool, bufsize, f_allocate));
// #endif

//   Kokkos::wait(sched);
}

#endif

