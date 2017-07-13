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
    struct TaskFunctor_UpdateSupernodes {
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
      ordinal_type _sid, _srow;
      
    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_UpdateSupernodes() = delete;
      
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_UpdateSupernodes(const sched_type &sched,
                                   const memory_pool_type &pool,
                                   const supernode_info_type &info,
                                   const ordinal_type sid,
                                   const ordinal_type srow)
        : _sched(sched),
          _pool(pool),
          _info(info),
          _sid(sid),
          _srow(srow) {}
        
      KOKKOS_INLINE_FUNCTION
      void operator()(member_type &member, value_type &r_val) {
        if (get_team_rank(member) == 0) {
          const bool is_deallocate = (_srow == -2);
          if (_srow == -1) {
            const size_type
              sbeg = _info.sid_super_panel_ptr(_sid)+1,
              send = _info.sid_super_panel_ptr(_sid+1)-1,
              srows = send - sbeg;
            
            future_type depbuf[MaxDependenceSize], *dep = &depbuf[0];
            const size_type depsize = (srows > MaxDependenceSize ? srows*sizeof(future_type) : 0);
            if (depsize) {
              dep = (future_type*)_pool.allocate(depsize);
              TACHO_TEST_FOR_ABORT(dep == NULL, "pool allocation fails");
            }
            
            for (ordinal_type i=sbeg;i<send;++i) {
              const ordinal_type tgt = _info.sid_super_panel_colidx(i);
              auto f_tgt = _info.supernodes_future(tgt);
              if (f_tgt.is_null()) {
                auto f = Kokkos::task_spawn
                  (Kokkos::TaskSingle(_sched, Kokkos::TaskPriority::High),
                   TaskFunctor_UpdateSupernodes(_sched, _pool, _info, _sid, i));
                _info.supernodes_future(tgt) = f;
                dep[i-sbeg] = f;
              } else {
                auto f = Kokkos::task_spawn
                  (Kokkos::TaskSingle(f_tgt, Kokkos::TaskPriority::High),
                   TaskFunctor_UpdateSupernodes(_sched, _pool, _info, _sid, i));
                _info.supernodes_future(tgt) = f;
                dep[i-sbeg] = f;
              }
            }

            --_srow;
            Kokkos::respawn(this, Kokkos::when_all(dep, srows), Kokkos::TaskPriority::High);

            if (depsize) {
              _pool.deallocate((void*)dep, depsize);
            }
          }
          
          if (_srow >= 0) {
            const Kokkos::pair<size_type,size_type> srow(_srow, _srow+1);
            CholSupernodes<Algo::Workflow::Serial>
              ::update(_sched, _pool, 
                       member, 
                       _info, _info.supernodes_abr(_sid),
                       srow, _sid,
                       0, NULL);
          }

          if (is_deallocate) {
            auto ABR = _info.supernodes_abr(_sid);
            //_pool.deallocate((void*)ABR.data(), ABR.span()*sizeof(value_type));
            _info.supernodes_abr(_sid) = typename supernode_info_type::value_type_matrix();
          }
        }
      }
    };

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
              future_type depbuf[MaxDependenceSize+1], *dep = &depbuf[0];
              const size_type depsize = (isize > MaxDependenceSize ? (isize+1)*sizeof(future_type) : 0);
              if (depsize) {
                dep = (future_type*)_pool.allocate(depsize);
                TACHO_TEST_FOR_ABORT(dep == NULL, "pool allocation fails");
              }
              
              // spawn child tasks
              for (ordinal_type i=0;i<isize;++i) {
                // the first child has a higher priority
                const auto priority = i ? Kokkos::TaskPriority::Low : Kokkos::TaskPriority::High;

                const ordinal_type child = _info.stree_children(i+ibeg);
                auto f = Kokkos::task_spawn(Kokkos::TaskSingle(_sched, priority),
                                            TaskFunctor_CholSupernodes(_sched, _pool, _info, child, _sid));
                TACHO_TEST_FOR_ABORT(f.is_null(), "task allocation fails");
                dep[i] = f;
              }
              dep[isize] = _info.supernodes_future(_sid);
              
              // respawn with updating state
              ++_state;
              Kokkos::respawn(this, Kokkos::when_all(dep, isize+1), Kokkos::TaskPriority::Regular);
              
              // if memory pool is used, deallocate the memory
              if (depsize) {
                _pool.deallocate((void*)dep, depsize);
              }
            } else {
              is_execute = true;
            }
          } 
          
          if (is_execute) {
            CholSupernodes<Algo::Workflow::Serial>
              ::factorize(_sched, _pool, 
                          member, 
                          _info, _sid, _sidpar,
                          _info.supernodes_abr(_sid),
                          0, NULL);

            _info.supernodes_future(_sid) = future_type();

            Kokkos::task_spawn(Kokkos::TaskSingle(_sched, Kokkos::TaskPriority::High),
                               TaskFunctor_UpdateSupernodes<mat_value_type,exec_space>(_sched, _pool, _info, _sid, -1));
          }
        }
      }
    };
  }
}

#endif

            
            

