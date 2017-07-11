#ifndef __TACHOEXP_TASKFUNCTOR_FACTORIZE_CHOL_HPP__
#define __TACHOEXP_TASKFUNCTOR_FACTORIZE_CHOL_HPP__

/// \file TachoExp_TaskFunctor_FactorizeChol.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Util.hpp"

#include "TachoExp_CholSupernodes.hpp"
#include "TachoExp_CholSupernodes_Serial.hpp"

namespace Tacho {
  
  namespace Experimental {
    
    template<typename MatValueType, typename ExecSpace>
    struct TaskFunctor_FactorizeChol {
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
      
      // respawn control
      ordinal_type _state;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_FactorizeChol() = delete;

      KOKKOS_INLINE_FUNCTION
      TaskFunctor_FactorizeChol(const sched_type &sched,
                                 const memory_pool_type &bufpool,
                                 const supernode_info_type &info,
                                 const ordinal_type sid)                                     
        : _sched(sched),
          _bufpool(bufpool),
          _info(info),
          _sid(sid),
          _state(0) {}
      
      KOKKOS_INLINE_FUNCTION
      void operator()(member_type &member, value_type &r_val) {
        if (get_team_rank(member) == 0) {
          TACHO_TEST_FOR_ABORT(_state > 2, "dead lock");

          // children information
          const ordinal_type 
            ibeg = _info.stree_ptr(_sid), 
            iend = _info.stree_ptr(_sid+1),
            isize = iend - ibeg;

          switch (_state) {
          case 0: { // tree parallelsim
            if (_info.serial_thres_size > _info.max_decendant_supernode_size(_sid)) {
              _state = 1; 
              Kokkos::respawn(this, future_type(), Kokkos::TaskPriority::Regular);
            } else {
              // allocate dependence array to handle variable number of children schur contributions
              future_type depbuf[MaxDependenceSize] /* 3 */, *dep = &depbuf[0];
              const size_type depsize = (isize > MaxDependenceSize ? isize*sizeof(future_type) : 0);
              if (depsize) {
                dep = (future_type*)_sched.memory()->allocate(depsize);
                TACHO_TEST_FOR_ABORT(dep == NULL, "sched memory pool allocation fails");
                clear((char*)dep, depsize);
              }
              
              // spawn children tasks and this (their parent) depends on the children tasks
              for (ordinal_type i=0;i<isize;++i) {
                // the first child has a higher priority
                const auto priority = (i ? Kokkos::TaskPriority::Low : Kokkos::TaskPriority::High);
                
                const ordinal_type child = _info.stree_children(i+ibeg);
                auto f = Kokkos::task_spawn(Kokkos::TaskSingle(_sched, priority),
                                            TaskFunctor_FactorizeChol(_sched, _bufpool, _info, child));
                TACHO_TEST_FOR_ABORT(f.is_null(), "task allocation fails");
                dep[i] = f;
              }
              
              // respawn with updating state
              _state = 2;
              Kokkos::respawn(this, Kokkos::when_all(dep, isize), Kokkos::TaskPriority::Regular);
              
              // deallocate dependence array
              if (depsize) {
                // manually reset future to decrease the reference count
                for (ordinal_type i=0;i<isize;++i) (dep+i)->~future_type();
                _sched.memory()->deallocate((void*)dep, depsize);
              }
            }
            break;
          }
          case 1: { 
            const ordinal_type n = _info.max_decendant_schur_size(_sid);
            const size_type bufsize = (n*n + _info.max_schur_size)*sizeof(mat_value_type);

            mat_value_type *buf = bufsize > 0 ? (mat_value_type*)_bufpool.allocate(bufsize) : NULL;
            TACHO_TEST_FOR_ABORT(buf == NULL && bufsize != 0, "bufmemory pool allocation fails");   

            CholSupernodes<Algo::Workflow::Serial>
              ::factorize_recursive_serial(_sched, member, _info, _sid, true, buf, bufsize);

            _bufpool.deallocate(buf, bufsize);

            _state = 3; // done
            break;
          }
          case 2: {
            ordinal_type pm, pn; _info.getSuperPanelSize(_sid, pm, pn);
            const ordinal_type n = pn - pm;
            const size_type bufsize = (n*n + _info.max_schur_size)*sizeof(mat_value_type);
            mat_value_type *buf = bufsize > 0 ? (mat_value_type*)_bufpool.allocate(bufsize) : NULL;
            TACHO_TEST_FOR_ABORT(buf == NULL && bufsize != 0, "bufmemory pool allocation fails");   

            CholSupernodes<Algo::Workflow::Serial>
              ::factorize_recursive_serial(_sched, member, _info, _sid, false, buf, bufsize);

            _bufpool.deallocate(buf, bufsize);

            _state = 3; // done 
            break;
          }
          }         
        }
      }
    };
  }
}

#endif

            
            

