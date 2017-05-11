#ifndef __TACHOEXP_TASKFUNCTOR_CHOLESKY_SUPERNODES_HPP__
#define __TACHOEXP_TASKFUNCTOR_CHOLESKY_SUPERNODES_HPP__

/// \file TachoExp_TaskFunctor_FactorizeSupernodes.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Util.hpp"

#include "TachoExp_CholSupernodes.hpp"
#include "TachoExp_CholSupernodes_Serial.hpp"

namespace Tacho {

  namespace Experimental {

    template<typename MatValueType, typename ExecSpace>
    struct TaskFunctor_CholSupernodes {
    public:
      typedef ExecSpace exec_space;

      typedef Kokkos::TaskScheduler<exec_space> sched_type;
      typedef typename sched_type::member_type member_type;
      typedef int value_type; // functor return type
      typedef Kokkos::Future<int,exec_space> future_type;

      typedef Kokkos::MemoryPool<exec_space> memory_pool_type;
      typedef MatValueType mat_value_type; // matrix value type

      typedef SupernodeInfo<mat_value_type,exec_space> supernode_info_type;
      typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;

    private:
      sched_type _sched;
      memory_pool_type _pool;
      
      supernode_info_type _info;
      ordinal_type _sid, _sidpar;
      
      // respawn control
      ordinal_type _state;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_CholSupernodes() = delete;

      KOKKOS_INLINE_FUNCTION
      TaskFunctor_CholSupernodes(const sched_type &sched,
                                 const memory_pool_type &pool,
                                 const supernode_info_type &info,
                                 const ordinal_type sid,
                                 const ordinal_type sidpar)                                     
        : _sched(sched),
          _pool(pool),
          _info(info),
          _sid(sid),
          _sidpar(sidpar),
          _state(0) {}
      
      KOKKOS_INLINE_FUNCTION
      void operator()(member_type &member, value_type &r_val) {
        if (get_team_rank(member) == 0) {
          TACHO_TEST_FOR_ABORT(_state > 1, "dead lock");
          bool is_execute = (_state == 1);
          if (_state == 0) {
            ///
            ///   have children   : invoke children tasks recursively
            ///   leaf (no child) : compute factorization
            ///
            const ordinal_type 
              ibeg = _info.stree_ptr(_sid), 
              iend = _info.stree_ptr(_sid+1),
              isize = iend - ibeg;

            if (isize > 0) {
              future_type depbuf[MaxDependenceSize], *dep = &depbuf[0];
              const size_type depsize = (isize > MaxDependenceSize ? isize*sizeof(future_type) : 0);
              if (depsize) {
                dep = (future_type*)_pool.allocate(depsize);
                TACHO_TEST_FOR_ABORT(dep == NULL, "pool allocation fails");
              }
              
              // spawn child tasks
              for (ordinal_type i=0;i<isize;++i) {
                // the first child has a higher priority
                const auto priority = i ? Kokkos::TaskPriority::Low : Kokkos::TaskPriority::High;

                const ordinal_type parent = _sid, child = _info.stree_children(i+ibeg);
                auto f = Kokkos::task_spawn(Kokkos::TaskSingle(_sched, priority),
                                            TaskFunctor_CholSupernodes(_sched, _pool, _info, child, parent));
                TACHO_TEST_FOR_ABORT(f.is_null(), "task allocation fails");
                dep[i] = f;
              }

              // respawn with updating state
              ++_state;
              Kokkos::respawn(this, Kokkos::when_all(dep, isize), Kokkos::TaskPriority::High);
              
              // if memory pool is used, deallocate the memory
              if (depsize) {
                _pool.deallocate((void*)dep, depsize);
              }
            } else {
              is_execute = true;
            }
          } 
          
          if (is_execute) {
            // for (ordinal_type i=0;i<isize;++i) {
            //   const ordinal_type child = _info.stree_children(i+ibeg);
              
            //   CholSupernodes<Algo::Workflow::Serial>
            //     ::update(_sched, _pool, member, 
            //              _info, ABR, child,
            //              0, NULL);
            //   _pool.deallocate((void*)ABR.data(), ABR.dimension_0()*ABR.dimension_1());
            // }
            CholSupernodes<Algo::Workflow::Serial>
              ::factorize(_sched, _pool, member, 
                          _info, _sid, _sidpar,
                          0, NULL);
          }          
        }
      }
    };
  }
}

#endif
