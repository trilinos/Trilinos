#ifndef __TACHOEXP_TASKFUNCTOR_SOLVE_UPPER_CHOL_HPP__
#define __TACHOEXP_TASKFUNCTOR_SOLVE_UPPER_CHOL_HPP__

/// \file TachoExp_TaskFunctor_FactorizeChol.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Util.hpp"

#include "TachoExp_CholSupernodes.hpp"
#include "TachoExp_CholSupernodes_Serial.hpp"

namespace Tacho {
  
  namespace Experimental {
    
    template<typename MatValueType, typename ExecSpace>
    struct TaskFunctor_SolveUpperChol {
    public:
      typedef ExecSpace exec_space;

      typedef Kokkos::TaskScheduler<exec_space> sched_type;
      typedef typename sched_type::member_type member_type;
      
      typedef Kokkos::MemoryPool<exec_space> memory_pool_type;
      
      typedef int value_type; // functor return type
      typedef Kokkos::Future<int,exec_space> future_type;

      typedef MatValueType mat_value_type; // matrix value type

      typedef SupernodeInfo<mat_value_type,exec_space> supernode_info_type;
      typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
      
    private:
      sched_type _sched;
      memory_pool_type _bufpool;

      supernode_info_type _info;
      ordinal_type _sid;
      
    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_SolveUpperChol() = delete;
      
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_SolveUpperChol(const sched_type &sched,
                                 const memory_pool_type &bufpool,
                                 const supernode_info_type &info,
                                 const ordinal_type sid)                                     
        : _sched(sched),
          _bufpool(bufpool),
          _info(info),
          _sid(sid) {}
      
      KOKKOS_INLINE_FUNCTION
      void operator()(member_type &member, value_type &r_val) {
        if (get_team_rank(member) == 0) {

          const auto &s = _info.supernodes(_sid);

          if (_info.serial_thres_size > s.max_decendant_supernode_size) {
            const ordinal_type n = s.max_decendant_schur_size, nrhs = _info.x.dimension_1();
            const size_type bufsize = n*nrhs*sizeof(mat_value_type);
            
            mat_value_type *buf = bufsize > 0 ? (mat_value_type*)_bufpool.allocate(bufsize) : NULL;
            TACHO_TEST_FOR_ABORT(buf == NULL && bufsize != 0, "bufmemory pool allocation fails");   
            
            CholSupernodes<Algo::Workflow::Serial>
              ::solve_upper_recursive_serial(_sched, member, _info, _sid, true, buf, bufsize);
            
            if (bufsize) 
              _bufpool.deallocate(buf, bufsize);

          } else {
            const ordinal_type n = s.n - s.m, nrhs = _info.x.dimension_1();
            const size_type bufsize = n*nrhs*sizeof(mat_value_type);
            
            mat_value_type *buf = bufsize > 0 ? (mat_value_type*)_bufpool.allocate(bufsize) : NULL;
            TACHO_TEST_FOR_ABORT(buf == NULL && bufsize != 0, "bufmemory pool allocation fails");   
            
            CholSupernodes<Algo::Workflow::Serial>
              ::solve_upper_recursive_serial(_sched, member, _info, _sid, false, buf, bufsize);
            
            if (bufsize > 0)
              _bufpool.deallocate(buf, bufsize);
            
            // allocate dependence array to handle variable number of children schur contributions
            future_type dep[MaxDependenceSize]; /* 4 */
            
            // spawn children tasks and this (their parent) depends on the children tasks
            for (ordinal_type i=0;i<s.nchildren;++i) {
              auto f = Kokkos::task_spawn(Kokkos::TaskSingle(_sched, Kokkos::TaskPriority::Regular),
                                          TaskFunctor_SolveUpperChol(_sched, _bufpool, _info, s.children[i]));
              TACHO_TEST_FOR_ABORT(f.is_null(), "task allocation fails");
              dep[i] = f;
            }            
          }
        }
      }
    };
  }
}

#endif

            
            

