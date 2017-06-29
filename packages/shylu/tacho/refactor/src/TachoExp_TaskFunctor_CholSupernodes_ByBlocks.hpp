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

      typedef SupernodeInfo<mat_value_type,exec_space> supernode_info_type;

      typedef DenseBlock<mat_value_type,exec_space> dense_block_type;
      typedef Kokkos::View<dense_block_type**,Kokkos::LayoutLeft,exec_space> dense_block_type_matrix;

    private:
      sched_type _sched;
      memory_pool_type _bufpool;
      
      supernode_info_type _info;
      ordinal_type _sid, _sidpar;

      // respawn control
      ordinal_type _state;

      // temporary future matrix
      ordinal_type _mb;
      UnmanagedViewType<future_type_matrix> _atl, _atr, _abr;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_CholSupernodes_ByBlocks() = delete;

      KOKKOS_INLINE_FUNCTION
      TaskFunctor_CholSupernodes_ByBlocks(const sched_type &sched,
                                          const memory_pool_type &bufpool,
                                          const supernode_info_type &info,
                                          const ordinal_type sid,
                                          const ordinal_type sidpar,
                                          const ordinal_type mb)                                     
        : _sched(sched),
          _bufpool(bufpool),
          _info(info),
          _sid(sid),
          _sidpar(sidpar),
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
            break;
          }
          case 1: {
            ///
            /// matrix parallelism
            ///
            std::cout << " state 1, sid = " << _sid << " Allocate matrix of blocks \n";
            ordinal_type m, tmp; _info.getSuperPanelSize(_sid, m, tmp);
            
            const ordinal_type 
              n = tmp - m,
              bm = m/_mb + (m%_mb > 0),
              bn = n/_mb + (n%_mb > 0);
            
            dense_block_type
              *atl_buf = (dense_block_type*)_sched.memory()->allocate(bm*bm*sizeof(dense_block_type)),
              *atr_buf = (dense_block_type*)_sched.memory()->allocate(bm*bn*sizeof(dense_block_type)),
              *abr_buf = (dense_block_type*)_sched.memory()->allocate(bn*bn*sizeof(dense_block_type));

            TACHO_TEST_FOR_ABORT(atl_buf == NULL || atr_buf == NULL || abr_buf == NULL, 
                                 "sched memory pool allocation fails");            

            _atl = UnmanagedViewType<dense_block_type_matrix>(atl_buf, bm, bm);
            _atr = UnmanagedViewType<dense_block_type_matrix>(atr_buf, bm, bn);
            _abr = UnmanagedViewType<dense_block_type_matrix>(abr_buf, bn, bn);

            auto set_matrix_of_blocks = [&](const UnmanagedViewType<dense_block_type_matrix> &mat) {
              dense_block_type blk; blk.rs = 1;
              const ordinal_type mm = mat.dimension_0(), nn = mat.dimension_1();              
              for (ordinal_type j=0;j<nn;++j) {
                const ordinal_type 
                  col_beg = j*_mb, col_tmp = col_beg + _mb, 
                  ncols = col_tmp > m ? m - col_beg : _mb; 
                blk.n = ncols;
                for (ordinal_type i=0;i<mm;++i) {
                  const ordinal_type 
                    row_beg = i*_mb, row_tmp = row_beg + _mb, 
                    nrows = row_tmp > m ? m - row_beg : _mb; 
                  blk.cs = blk.m = nrows;
                  blk.buf = (mat_value_type*)_bufpool.allocate(blk.m*blk.n*sizeof(mat_value_type));
                  TACHO_TEST_FOR_ABORT(blk.buf == NULL, "buf memory pool allocation fails");

                  mat(i,j) = blk;
                }
              }
            };
            set_matrix_of_blocks(_atl);
            set_matrix_of_blocks(_atr);
            set_matrix_of_blocks(_abr);

            std::cout << " state 1, sid = " << _sid << " Spawn algorithm by blocks  \n";              
            std::cout << " state 1, sid = " << _sid << " Respawn this with dependences  \n";              

            const size_type depsize = _abr.span()*sizeof(future_type);
            future_type *dep = (future_type*)_sched.memory()->allocate(depsize);
            TACHO_TEST_FOR_ABORT(dep == NULL, "sched memory pool allocation fails");
            
            ordinal_type cnt = 0;
            for (ordinal_type j=0;j<_abr.dimension_1();++j)
              for (ordinal_type i=0;i<_abr.dimension_0();++i)
                dep[cnt++] = _abr(i,j).future;

            _state = 2;
            Kokkos::respawn(this, Kokkos::when_all(dep, cnt), Kokkos::TaskPriority::Regular);

            for (ordinal_type k=0;k<cnd;++k) dep[k] = future_type();
            _sched.memory()->deallocate((void*)dep, depsize);
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

            
            

