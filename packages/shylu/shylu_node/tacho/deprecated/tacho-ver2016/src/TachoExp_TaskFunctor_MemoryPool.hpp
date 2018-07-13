#ifndef __TACHOEXP_TASKFUNCTOR_MEMORYPOOL_HPP__
#define __TACHOEXP_TASKFUNCTOR_MEMORYPOOL_HPP__

/// \file TachoExp_TaskFunctor_MemoryPool.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Util.hpp"

namespace Tacho {

  namespace Experimental {

    template<typename ExecSpace>
    struct TaskFunctor_MemoryPool_Allocate {
      typedef ExecSpace exec_space;
      
      // task scheduler/future
      typedef Kokkos::TaskScheduler<exec_space> sched_type;
      typedef typename sched_type::member_type member_type;
      typedef Kokkos::Future<void*,exec_space> future_type;

      // memory pool
#if defined(__KK__)
      typedef Kokkos::MemoryPool<exec_space> memory_pool_type;
#else
      typedef Kokkos::Experimental::MemoryPool<exec_space> memory_pool_type;
#endif      

      typedef void* value_type; // functor return type

    private:
      memory_pool_type _pool;
      size_type _bufsize;

    public:

      KOKKOS_INLINE_FUNCTION
      TaskFunctor_MemoryPool_Allocate() = delete;

      KOKKOS_INLINE_FUNCTION
      TaskFunctor_MemoryPool_Allocate(const memory_pool_type &pool,
                                      const size_type bufsize)
        : _pool(pool),
          _bufsize(bufsize) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(member_type &member, value_type &r_val) {
        printf("task pool allocate\n");
        if (get_team_rank(member) == 0) {
          if (_bufsize)
            r_val = (void*)_pool.allocate(_bufsize);
          else
            r_val = NULL;
        }
      }
    };

    template<typename ExecSpace>
    struct TaskFunctor_MemoryPool_Deallocate {
      typedef ExecSpace exec_space;
      
      // task scheduler/future
      typedef Kokkos::TaskScheduler<exec_space> sched_type;
      typedef typename sched_type::member_type member_type;
      typedef Kokkos::Future<int,exec_space> future_type;

      typedef int value_type;

      // memory pool
#if defined(__KK__)
      typedef Kokkos::MemoryPool<exec_space> memory_pool_type;
#else
      typedef Kokkos::Experimental::MemoryPool<exec_space> memory_pool_type;
#endif      

      typedef Kokkos::Future<void*,exec_space> future_ptr_type;

    private:
      memory_pool_type _pool;
      size_type _bufsize;
      
      future_ptr_type _future_ptr;

    public:

      KOKKOS_INLINE_FUNCTION
      TaskFunctor_MemoryPool_Deallocate() = delete;

      KOKKOS_INLINE_FUNCTION
      TaskFunctor_MemoryPool_Deallocate(const memory_pool_type &pool,
                                        const size_type bufsize) {
        //                                        const future_ptr_type &future_ptr) 
        : _pool(pool),
          _bufsize(bufsize) {}
        //          _future_ptr(future_ptr) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(member_type &member, int &r_val) {
        printf("dealloc\n");
        r_val = 1;
        // if (get_team_rank(member) == 0) {
        //   if (_bufsize)
        //     _pool.deallocate((void*)future_ptr.get(), _bufsize);
        // }
      }
    };

    template<typename ExecSpace>
    struct TaskFunctor_MemoryPool_TestView {
      typedef ExecSpace exec_space;
      
      // task scheduler/future
      typedef Kokkos::TaskScheduler<exec_space> sched_type;
      typedef typename sched_type::member_type member_type;
      typedef Kokkos::Future<exec_space> future_type;

      typedef Kokkos::Future<void*,exec_space> future_ptr_type;
      typedef Kokkos::View<double**,Kokkos::LayoutLeft,exec_space> value_type;

      // memory pool
#if defined(__KK__)
      typedef Kokkos::MemoryPool<exec_space> memory_pool_type;
#else
      typedef Kokkos::Experimental::MemoryPool<exec_space> memory_pool_type;
#endif      

    private:
      memory_pool_type _pool;
      ordinal_type _m, _n;
      
      Kokkos::Future<void*,exec_space> _future_ptr;

    public:

      inline
      TaskFunctor_MemoryPool_TestView() = delete;

      inline
      TaskFunctor_MemoryPool_TestView(const memory_pool_type &pool,
                                      const ordinal_type m,
                                      const ordinal_type n,
                                      const Kokkos::Future<void*,exec_space> &future_ptr)
        : _pool(pool),
          _m(m),
          _n(n),
          _future_ptr(future_ptr) {}
      
      inline
      void operator()(member_type &member, value_type &r_val) {
        if (get_team_rank(member) == 0) {
          if (_m && _n) {
            Kokkos::View<double**,Kokkos::LayoutLeft,exec_space> A((double*)_future_ptr.get(), _m, _n);
            ordinal_type cnt = 0;
            for (ordinal_type i=0;i<_m;++i)
              for (ordinal_type j=0;j<_n;++j)
                A(i,j) = cnt++;
            r_val = A;
          } else {
            r_val = value_type();
          }
        }
      }
    };

    template<typename ExecSpace>
    struct TaskFunctor_MemoryPool_TestViewSee {
      typedef ExecSpace exec_space;
      
      // task scheduler/future
      typedef Kokkos::TaskScheduler<exec_space> sched_type;
      typedef typename sched_type::member_type member_type;
      typedef Kokkos::Future<exec_space> future_type;

      typedef Kokkos::View<double**,Kokkos::LayoutLeft,exec_space> view_type;
      typedef Kokkos::Future<view_type,exec_space> future_view_type;
      
      typedef int value_type;

      // memory pool
#if defined(__KK__)
      typedef Kokkos::MemoryPool<exec_space> memory_pool_type;
#else
      typedef Kokkos::Experimental::MemoryPool<exec_space> memory_pool_type;
#endif      

    private:
      memory_pool_type _pool;
      ordinal_type _m, _n;

      future_view_type _A;

    public:

      inline
      TaskFunctor_MemoryPool_TestViewSee() = delete;

      inline
      TaskFunctor_MemoryPool_TestViewSee(const memory_pool_type &pool,
                                         const future_view_type &A)
        : _pool(pool),
          _A(A) {}
      
      inline
      void operator()(member_type &member, value_type &r_val) {
        if (get_team_rank(member) == 0) {
          const auto A = _A.get();

          const ordinal_type m = A.extent(0);
          const ordinal_type n = A.extent(1);

          printf("A in TestViewSee\n");
          for (ordinal_type i=0;i<m;++i) {
            for (ordinal_type j=0;j<n;++j)
              printf(" %4d ", int(A(i,j)));
            printf("\n");
          }
          r_val = 1;
        }
      }
    };

  }
}

#endif
