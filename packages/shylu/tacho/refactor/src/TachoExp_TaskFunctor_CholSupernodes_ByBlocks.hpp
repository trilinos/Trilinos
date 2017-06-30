#ifndef __TACHOEXP_TASKFUNCTOR_CHOLESKY_SUPERNODES_BYBLOCKS_HPP__
#define __TACHOEXP_TASKFUNCTOR_CHOLESKY_SUPERNODES_BYBLOCKS_HPP__

/// \file TachoExp_TaskFunctor_FactorizeSupernodes.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Util.hpp"

#include "TachoExp_CholSupernodes.hpp"
#include "TachoExp_CholSupernodes_ByBlocks.hpp"

namespace Tacho {
  
  namespace Experimental {
    
    template<typename MatValueType, typename ExecSpace>
    struct TaskFunctor_CholSupernodes_ByBlocks {
    public:
      typedef ExecSpace exec_space;

      typedef Kokkos::TaskScheduler<exec_space> sched_type;
      typedef typename sched_type::member_type member_type;
      typedef int value_type; // functor return type
      typedef Kokkos::Future<int,exec_space> future_type;

      typedef Kokkos::MemoryPool<exec_space> memory_pool_type;
      typedef MatValueType mat_value_type; // matrix value type
      typedef Kokkos::View<mat_value_type**,Kokkos::LayoutLeft,exec_space> mat_value_type_matrix;

      typedef SupernodeInfo<mat_value_type,exec_space> supernode_info_type;

      typedef DenseBlock<mat_value_type,exec_space> dense_block_type;
      typedef Kokkos::View<dense_block_type**,Kokkos::LayoutLeft,exec_space> dense_block_type_matrix;

    private:
      sched_type _sched;
      memory_pool_type _bufpool;
      
      supernode_info_type _info;
      ordinal_type _sid;

      // respawn control
      ordinal_type _state;

      // blocksize
      ordinal_type _mb;

      // matrix of blocks
      UnmanagedViewType<dense_block_type_matrix> _abr;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_CholSupernodes_ByBlocks() = delete;

      KOKKOS_INLINE_FUNCTION
      TaskFunctor_CholSupernodes_ByBlocks(const sched_type &sched,
                                          const memory_pool_type &bufpool,
                                          const supernode_info_type &info,
                                          const ordinal_type sid,
                                          const ordinal_type mb)                                     
        : _sched(sched),
          _bufpool(bufpool),
          _info(info),
          _sid(sid),
          _state(0),
          _mb(mb) {}
      
      KOKKOS_INLINE_FUNCTION
      void operator()(member_type &member, value_type &r_val) {
        if (get_team_rank(member) == 0) {
          // children information
          const ordinal_type 
            ibeg = _info.stree_ptr(_sid), 
            iend = _info.stree_ptr(_sid+1),
            isize = iend - ibeg;
          
          // no children, skip spawning
          if (isize == 0 && _state == 0) _state = 1;
          
          switch (_state) {
          case 2: {
            ///
            /// atomic update to its parent 
            ///
            if (_abr.span() > 0) {
              const ordinal_type bm = _abr.dimension_0(), bn = _abr.dimension_1();
              deallocateStorageByBlocks(_abr, _bufpool);

              for (ordinal_type j=0;j<bn;++j)
                for (ordinal_type i=0;i<bm;++i)
                  _abr(i,j).future.~future_type();
              
              _sched.memory()->deallocate(_abr.data(), bm*bn*sizeof(dense_block_type));            
            }
            _state = 3; // done
            break;
          }
          case 1: {
            ///
            /// matrix parallelism
            ///

            // get supernode panel pointer
            mat_value_type *ptr = info.getSuperPanelPtr(sid);

            // characterize the panel size
            ordinal_type pm, pn;
            info.getSuperPanelSize(sid, pm, pn);

            // panel is divided into diagonal and interface block (i.e., ATL and ATR)
            const ordinal_type m = pm, n = pn - pm;
            
            // block matrix size
            const ordinal_type bm = m/_mb + (m%_mb > 0), bn = n/_mb + (n%_mb > 0);

            // m and n are available, then factorize the supernode block
            if (m > 0) {
              dense_block_type
                *atl_buf = (dense_block_type*)_sched.memory()->allocate(bm*bm*sizeof(dense_block_type));
              TACHO_TEST_FOR_ABORT(atl_buf == NULL, "sched memory pool allocation fails");

              _atl = UnmanagedViewType<dense_block_type_matrix>(atl_buf, bm, bm);
              setMatrixOfBlocks(_atl, m, m, _mb);
              attachBaseBuffer(_atl, ptr, m, m, _mb); ptr += m*m;

              if (n > 0) {
                dense_block_type
                  *atr_buf = (dense_block_type*)_sched.memory()->allocate(bm*bn*sizeof(dense_block_type));
                TACHO_TEST_FOR_ABORT(atr_buf == NULL, "sched memory pool allocation fails");

                UnmanagedViewType<dense_block_type_matrix> atr(atr_buf, bm, bn);
                setMatrixOfBlocks(atr, m, n, _mb);
                attachBaseBuffer(atr, ptr, m, n, _mb);
              
                dense_block_type
                  *abr_buf = (dense_block_type*)_sched.memory()->allocate(bn*bn*sizeof(dense_block_type));

                _abr = UnmanagedViewType<dense_block_type_matrix>(abr_buf, bn, bn);
                setMatrixOfBlocks(_abr, n, n, _mb);
                allocateStorageByBlocks(_abr, _bufpool);

                // in place destruction of atr
                for (ordinal_type j=0;j<bn;++j)
                  for (ordinal_type i=0;i<bm;++i)
                    atr(i,j).future.~future_type();

                _sched.memory()->deallocate(atr_buf, bm*bn*sizeof(dense_block_type));
              }

              // in place destruction of atl
              for (ordinal_type j=0;j<bm;++j)
                for (ordinal_type i=0;i<bm;++i)
                  atl(i,j).future.~future_type();

              _sched.memory()->deallocate(atl_buf, bm*bm*sizeof(dense_block_type));
            }

            if (_abr.span() > 0) {
              const size_type depsize = _abr.span()*sizeof(future_type);
              future_type *dep = (future_type*)_sched.memory()->allocate(depsize);
              TACHO_TEST_FOR_ABORT(dep == NULL, "sched memory pool allocation fails");
              
              ordinal_type cnt = 0;
              for (ordinal_type j=0;j<_abr.dimension_1();++j)
                for (ordinal_type i=0;i<_abr.dimension_0();++i)
                  dep[cnt++] = _abr(i,j).future;
              
              _state = 2;
              Kokkos::respawn(this, Kokkos::when_all(dep, cnt), Kokkos::TaskPriority::Regular);
              
              // depedences
              for (ordinal_type k=0;k<cnt;++k) (dep+k)->~future_type();
              _sched.memory()->deallocate((void*)dep, depsize);
            } else {
              _state = 3; // done
            }
            break;
          }
          case 0: {
            ///
            /// tree parallelism
            ///
            future_type depbuf[MaxDependenceSize] /* 3 */, *dep = &depbuf[0];
            const size_type depsize = (isize > MaxDependenceSize ? isize*sizeof(future_type) : 0);
            if (depsize) {
              dep = (future_type*)_sched.memory()->allocate(depsize);
              TACHO_TEST_FOR_ABORT(dep == NULL, "sched memory pool allocation fails");
            }
            
            // spawn children tasks and this (their parent) depends on the children tasks
            for (ordinal_type i=0;i<isize;++i) {
              // the first child has a higher priority
              const auto priority = (i ? Kokkos::TaskPriority::Low : Kokkos::TaskPriority::High);
              
              const ordinal_type child = _info.stree_children(i+ibeg);
              auto f = Kokkos::task_spawn(Kokkos::TaskSingle(_sched, priority),
                                          TaskFunctor_CholSupernodes_ByBlocks(_sched, _bufpool, _info, child, _sid, _mb));
              TACHO_TEST_FOR_ABORT(f.is_null(), "task allocation fails");
              dep[i] = f;
            }
            
            // respawn with updating state
            _state = 1;
            Kokkos::respawn(this, Kokkos::when_all(dep, isize), Kokkos::TaskPriority::Regular);
              
            // deallocate dependence array
            if (depsize) {
              // manually reset future to decrease the reference count
              for (ordinal_type i=0;i<isize;++i) dep[i] = future_type();
              _sched.memory()->deallocate((void*)dep, depsize);
            }            
            break;
          }
          }
        }
      }
    };
  }
}

#endif

            
            

