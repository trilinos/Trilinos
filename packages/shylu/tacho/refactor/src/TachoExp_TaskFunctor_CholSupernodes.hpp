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
      // TODO:: for now, we use three memory pool
      // - sched.memory() - tasking allocation (small/uniform superblocks)
      // - workpool - misc allocation (small/many superblocks)
      // - bufpool - large/small number of superblocks
      // later merge workpool + sched.memory() but keep bufpool
      sched_type _sched;
      memory_pool_type _workpool, _bufpool;
      
      supernode_info_type _info;
      ordinal_type _sid, _sidpar;
      
      // respawn control
      ordinal_type _state;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_CholSupernodes() = delete;

      KOKKOS_INLINE_FUNCTION
      TaskFunctor_CholSupernodes(const sched_type &sched,
                                 const memory_pool_type &workpool,
                                 const supernode_info_type &info,
                                 const ordinal_type sid,
                                 const ordinal_type sidpar)                                     
        : _sched(sched),
          _workpool(workpool),
          _info(info),
          _sid(sid),
          _sidpar(sidpar),
          _state(0) {}
      
      KOKKOS_INLINE_FUNCTION
      void operator()(member_type &member, value_type &r_val) {
        if (get_team_rank(member) == 0) {
          TACHO_TEST_FOR_ABORT(_state > 1, "dead lock");

          // children information
          const ordinal_type 
            ibeg = _info.stree_ptr(_sid), 
            iend = _info.stree_ptr(_sid+1),
            isize = iend - ibeg;

          bool is_execute = (_state == 1 /* respawned */ || isize == 0 /* leaf */);
          if (is_execute) { 
            typedef typename supernode_info_type::value_type_matrix value_type_matrix;

            const ordinal_type 
              bbeg = _info.super_schur_ptr[_sid],
              bend = _info.super_schur_ptr[_sid+1],
              bufsize = (bend - bbeg)*sizeof(mat_value_type);
            
            void *buf = (void*)&_info.super_schur_buf[bbeg];
            CholSupernodes<Algo::Workflow::Serial>
              ::factorize(_sched, member,
                          _info, _sid, _sidpar,
                          bufsize, buf);

            ordinal_type m, n; _info.getSuperPanelSize(_sid, m, n);
            UnmanagedViewType<value_type_matrix> ABR((mat_value_type*)buf, n-m, n-m);

            CholSupernodes<Algo::Workflow::Serial>
              ::update(_sched, member,
                       _info, ABR, _sid,
                       bufsize - ABR.span()*sizeof(mat_value_type),
                       (void*)((mat_value_type*)buf + ABR.span()));
            
            // this is done
            _state = 2;
          } else {
            // allocate dependence array to handle variable number of children schur contributions
            future_type depbuf[MaxDependenceSize] /* 3 */, *dep = &depbuf[0];
            const size_type depsize = (isize > MaxDependenceSize ? isize*sizeof(future_type) : 0);
            if (depsize) {
              dep = (future_type*)_workpool.allocate(depsize);
              TACHO_TEST_FOR_ABORT(dep == NULL, "pool allocation fails");
            }

            // spawn children tasks and this (their parent) depends on the children tasks
            for (ordinal_type i=0;i<isize;++i) {
              // the first child has a higher priority
              const auto priority = (i ? Kokkos::TaskPriority::Low : Kokkos::TaskPriority::High);
              
              const ordinal_type child = _info.stree_children(i+ibeg);
              auto f = Kokkos::task_spawn(Kokkos::TaskSingle(_sched, priority),
                                          TaskFunctor_CholSupernodes(_sched, _workpool, _info, child, _sid));
              TACHO_TEST_FOR_ABORT(f.is_null(), "task allocation fails");
              dep[i] = f;
            }
            
            // respawn with updating state
            ++_state;
            Kokkos::respawn(this, Kokkos::when_all(dep, isize), Kokkos::TaskPriority::Regular);
              
            // deallocate dependence array
            if (depsize) {
              // manually reset future to decrease the reference count
              for (ordinal_type i=0;i<isize;++i) dep[i] = future_type();
              _workpool.deallocate((void*)dep, depsize);
            }
          } 

        }
      }
    };
  }
}

#endif

            
            

