#ifndef __TACHOEXP_TASKFUNCTOR_FACTORIZE_CHOLBYBLOCKS_PANEL_HPP__
#define __TACHOEXP_TASKFUNCTOR_FACTORIZE_CHOLBYBLOCKS_PANEL_HPP__

/// \file TachoExp_TaskFunctor_FactorizeCholByBlocks.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Util.hpp"

#include "TachoExp_Chol_ByBlocks.hpp"
#include "TachoExp_Trsm_ByBlocks.hpp"
#include "TachoExp_Herk_ByBlocks.hpp"
#include "TachoExp_Gemm_ByBlocks.hpp"

#include "TachoExp_CholSupernodes.hpp"
#include "TachoExp_CholSupernodes_SerialPanel.hpp"

#include "TachoExp_TaskFunctor_FactorizeCholPanel.hpp"

namespace Tacho {

  namespace Experimental {

    template<typename MatValueType, typename ExecSpace>
    struct TaskFunctor_FactorizeCholByBlocksPanel {
    public:
      typedef ExecSpace exec_space;

      typedef Kokkos::TaskScheduler<exec_space> sched_type;
      typedef typename sched_type::member_type member_type;

      typedef Kokkos::MemoryPool<exec_space> memory_pool_type;

      typedef int value_type; // functor return type
      typedef Kokkos::Future<int,exec_space> future_type;

      typedef MatValueType mat_value_type; // matrix value type

      typedef SupernodeInfo<mat_value_type,exec_space> supernode_info_type;
      typedef typename supernode_info_type::supernode_type supernode_type;

      typedef typename supernode_info_type::value_type_matrix mat_value_type_matrix;
      typedef typename supernode_info_type::dense_block_type dense_block_type;
      typedef typename supernode_info_type::dense_matrix_of_blocks_type dense_matrix_of_blocks_type;

    private:
      sched_type _sched;
      memory_pool_type _bufpool;

      supernode_info_type _info;
      ordinal_type _sid;

      supernode_type _s;

      ordinal_type _state;

      ordinal_type _mb, _nb, _offn;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_FactorizeCholByBlocksPanel() = delete;

      KOKKOS_INLINE_FUNCTION
      TaskFunctor_FactorizeCholByBlocksPanel(const sched_type &sched,
                                             const memory_pool_type &bufpool,
                                             const supernode_info_type &info,
                                             const ordinal_type sid,
                                             const ordinal_type mb,
                                             const ordinal_type nb,
                                             const ordinal_type state = 0,
                                             const ordinal_type offn = 0)
        : _sched(sched),
          _bufpool(bufpool),
          _info(info),
          _sid(sid),
          _s(info.supernodes(sid)),
          _state(state),
          _mb(mb),
          _nb(nb),
          _offn(offn) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(member_type &member, value_type &r_val) {
        if (get_team_rank(member) == 0) {
          constexpr ordinal_type done = 5;
          TACHO_TEST_FOR_ABORT(_state == done, "dead lock");

          if (_s.nchildren == 0 && _state == 0) _state = 1;

          switch (_state) {
          case 3: {
            ///
            /// partial update and its parent
            ///
            const ordinal_type n = _s.n - _s.m, nb = min(_nb, n - _offn), nn = _offn + nb;
            const size_type bufsize = n*sizeof(ordinal_type) + nn*nb*sizeof(mat_value_type);  
            mat_value_type *buf = bufsize > 0 ? (mat_value_type*)_bufpool.allocate(bufsize) : NULL;  
            if (buf == NULL && bufsize) {
              Kokkos::respawn(this, _sched, Kokkos::TaskPriority::Low);   
            } else { 
              CholSupernodes<Algo::Workflow::SerialPanel>
                ::update(_sched, member, _info, _offn, nb, _sid, bufsize, (void*)buf);
              _bufpool.deallocate(buf, bufsize); 
            }      
            break;
          }
          case 2: {
            ///
            /// update to its parent
            ///
            const ordinal_type n = _s.n - _s.m, bn = n/_nb + (n%_nb > 0), depsize = bn*sizeof(future_type);
            future_type *dep = (future_type*)_sched.memory()->allocate(depsize);  
            TACHO_TEST_FOR_ABORT(dep == NULL, "sched memory pool allocation fails"); 
            clear((char*)dep, depsize);  

            const ordinal_type state = 3;
            for (ordinal_type i=0;i<bn;++i) {
              auto f = Kokkos::task_spawn(Kokkos::TaskSingle(_sched, Kokkos::TaskPriority::Regular),
                                          TaskFunctor_FactorizeCholByBlocksPanel
                                          (_sched, _bufpool, _info, _sid, _mb, _nb, state, i*_nb));
              TACHO_TEST_FOR_ABORT(f.is_null(), "task allocation fails");
              dep[i] = f;
            }

            _state = 4;
            Kokkos::respawn(this, Kokkos::when_all(dep, bn), Kokkos::TaskPriority::Regular);

            for (ordinal_type k=0;k<bn;++k) (dep+k)->~future_type();
            if (depsize) _sched.memory()->deallocate(dep, depsize);

            break;
          }
          case 1: {
            ///
            /// matrix parallelism
            ///

            // get supernode panel pointer
            mat_value_type *ptr = _s.buf;

            // panel is divided into diagonal and interface block (i.e., ATL and ATR)
            const ordinal_type m = _s.m, n = _s.n - _s.m;

            // block matrix size
            const ordinal_type bm = m/_mb + (m%_mb > 0), bn = n/_mb + (n%_mb > 0);

            // allocation for matrix of blocks
            const size_t bufsize = (bm*bm + bm*bn)*sizeof(dense_block_type);
            char *buf = (char*)_bufpool.allocate(bufsize);

            dense_block_type *buf_ptr = (dense_block_type*)buf;

            if (buf == NULL && bufsize) {
              Kokkos::respawn(this, _sched, Kokkos::TaskPriority::Low);  
            } else {
              clear((char*)buf, bufsize);

              // m and n are available, then factorize the supernode block
              dense_matrix_of_blocks_type htr;
              if (m > 0) {
                dense_matrix_of_blocks_type htl(buf_ptr, bm, bm);
                buf_ptr += bm*bm;
                
                setMatrixOfBlocks(htl, m, m, _mb);
                attachBaseBuffer(htl, ptr, 1, m); ptr += m*m;
                
                // chol
                Chol<Uplo::Upper,Algo::ByBlocks>::invoke(_sched, member, htl);
                
                if (n > 0) {
                  htr = dense_matrix_of_blocks_type(buf_ptr, bm, bn);
                  buf_ptr += bm*bn;
                  
                  setMatrixOfBlocks(htr, m, n, _mb);
                  attachBaseBuffer(htr, ptr, 1, m);
                  
                  // trsm
                  Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,Algo::ByBlocks>
                    ::invoke(_sched, member, Diag::NonUnit(), 1.0, htl, htr);                  
                }
                for (ordinal_type k=0;k<(bm*bm);++k) htl[k].set_future();
              }

              if (bm > 0 && bn > 0) {
                const size_t bmn = bm*bn, depsize = bm*bn*sizeof(future_type);
                future_type *dep = (future_type*)_sched.memory()->allocate(depsize);  
                TACHO_TEST_FOR_ABORT(dep == NULL, "sched memory pool allocation fails"); 
                clear((char*)dep, depsize);  

                for (ordinal_type j=0,k=0;j<bn;++j)
                  for (ordinal_type i=0;i<bm;++i,++k) {
                    dep[k] = htr(i,j).future();
                    htr(i,j).set_future();
                  }

                _state = 2;
                Kokkos::respawn(this, Kokkos::when_all(dep, bmn), Kokkos::TaskPriority::Regular); 

                for (ordinal_type k=0;k<static_cast<ordinal_type>(bmn);++k) (dep+k)->~future_type();
                _sched.memory()->deallocate(dep, depsize);   
              } else {
                _state = 4;
              }
            }
            if (bufsize) _bufpool.deallocate(buf, bufsize);
            break;
          }
          case 0: {
            ///
            /// tree parallelism
            ///
            const bool use_byblocks = (_mb*1.5 < _s.max_decendant_supernode_size);

            // spawn children tasks and this (their parent) depends on the children tasks
            future_type dep[MaxDependenceSize];
            if (use_byblocks) {
              for (ordinal_type i=0;i<_s.nchildren;++i) {
                auto f = Kokkos::task_spawn(Kokkos::TaskSingle(_sched, Kokkos::TaskPriority::Regular),
                                            TaskFunctor_FactorizeCholByBlocksPanel
                                            (_sched, _bufpool, _info, _s.children[i], _mb, _nb));
                TACHO_TEST_FOR_ABORT(f.is_null(), "task allocation fails");
                dep[i] = f;
              }
            } else {
              for (ordinal_type i=0;i<_s.nchildren;++i) {
                auto f = Kokkos::task_spawn(Kokkos::TaskSingle(_sched, Kokkos::TaskPriority::Regular),
                                            TaskFunctor_FactorizeCholPanel<mat_value_type,exec_space>
                                            (_sched, _bufpool, _info, _s.children[i], _nb));
                TACHO_TEST_FOR_ABORT(f.is_null(), "task allocation fails");
                dep[i] = f;
              }
            }

            // respawn with updating state
            _state = 1;
            Kokkos::respawn(this, Kokkos::when_all(dep, _s.nchildren), Kokkos::TaskPriority::Regular);
            break;
          }
          }
        }
      }
    };
  }
}

#endif
