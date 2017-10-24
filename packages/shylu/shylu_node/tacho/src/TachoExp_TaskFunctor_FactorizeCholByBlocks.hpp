#ifndef __TACHOEXP_TASKFUNCTOR_FACTORIZE_CHOLBYBLOCKS_HPP__
#define __TACHOEXP_TASKFUNCTOR_FACTORIZE_CHOLBYBLOCKS_HPP__

/// \file TachoExp_TaskFunctor_FactorizeCholByBlocks.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Util.hpp"

#include "TachoExp_Chol_ByBlocks.hpp"
#include "TachoExp_Trsm_ByBlocks.hpp"
#include "TachoExp_Herk_ByBlocks.hpp"
#include "TachoExp_Gemm_ByBlocks.hpp"

#include "TachoExp_CholSupernodes.hpp"
#include "TachoExp_CholSupernodes_Serial.hpp"

#include "TachoExp_TaskFunctor_FactorizeChol.hpp"

namespace Tacho {

  namespace Experimental {

    template<typename MatValueType, typename ExecSpace>
    struct TaskFunctor_FactorizeCholByBlocks {
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

      ordinal_type _mb;

      // abr and matrix of blocks (deallocated after tasks are completed)
      char *_buf;
      size_t _bufsize;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_FactorizeCholByBlocks() = delete;

      KOKKOS_INLINE_FUNCTION
      TaskFunctor_FactorizeCholByBlocks(const sched_type &sched,
                                        const memory_pool_type &bufpool,
                                        const supernode_info_type &info,
                                        const ordinal_type sid,
                                        const ordinal_type mb)
        : _sched(sched),
          _bufpool(bufpool),
          _info(info),
          _sid(sid),
          _s(info.supernodes(sid)),
          _state(0),
          _mb(mb),
          _buf(NULL),
          _bufsize(0) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(member_type &member, value_type &r_val) {
        if (get_team_rank(member) == 0) {
          constexpr ordinal_type done = 3;
          TACHO_TEST_FOR_ABORT(_state == done, "dead lock");

          if (_s.nchildren == 0 && _state == 0) _state = 1;

          switch (_state) {
          case 2: {
            ///
            /// update to its parent
            ///
            const ordinal_type
              m = _s.m, n = _s.n - _s.m,
              bm = m/_mb + (m%_mb > 0), bn = n/_mb + (n%_mb > 0);

            const size_t
              matrix_of_blocks_bufsize = (bm*bm + bm*bn + bn*bn)*sizeof(dense_block_type);

            char *buf = _buf + matrix_of_blocks_bufsize;
            const size_t bufsize = _bufsize - matrix_of_blocks_bufsize;

            UnmanagedViewType<mat_value_type_matrix> ABR((mat_value_type*)buf, n, n);
            CholSupernodes<Algo::Workflow::Serial>
              ::update(_sched, member,
                       _info, ABR, _sid,
                       (bufsize - ABR.span()*sizeof(mat_value_type)),
                       (void*)((mat_value_type*)buf + ABR.span()));

            if (_bufsize)
              _bufpool.deallocate(_buf, _bufsize);

            _state = done;
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
            const size_t
              matrix_of_blocks_bufsize = (bm*bm + bm*bn + bn*bn)*sizeof(dense_block_type),
              abr_bufsize = (n*n + _info.max_schur_size)*sizeof(mat_value_type),
              bufsize = matrix_of_blocks_bufsize + abr_bufsize;

            char *buf = (char*)_bufpool.allocate(bufsize);
            if (buf == NULL && bufsize) {
              Kokkos::respawn(this, _sched, Kokkos::TaskPriority::Low);  
            } else {
              clear((char*)buf, matrix_of_blocks_bufsize);
              
              // assign buf information for deallocation
              _buf = buf; _bufsize = bufsize;
              
              dense_matrix_of_blocks_type hbr;
              
              // m and n are available, then factorize the supernode block
              if (m > 0) {
                dense_matrix_of_blocks_type htl((dense_block_type*)buf, bm, bm);
                buf += (bm*bm)*sizeof(dense_block_type);
                
                setMatrixOfBlocks(htl, m, m, _mb);
                attachBaseBuffer(htl, ptr, 1, m); ptr += m*m;
                
                // chol
                Chol<Uplo::Upper,Algo::ByBlocks>::invoke(_sched, member, htl);
                
                if (n > 0) {
                  dense_matrix_of_blocks_type htr((dense_block_type*)buf, bm, bn);
                  buf += (bm*bn)*sizeof(dense_block_type);
                  
                  setMatrixOfBlocks(htr, m, n, _mb);
                  attachBaseBuffer(htr, ptr, 1, m);
                  
                  // trsm
                  Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,Algo::ByBlocks>
                    ::invoke(_sched, member, Diag::NonUnit(), 1.0, htl, htr);
                  
                  hbr = dense_matrix_of_blocks_type((dense_block_type*)buf, bn, bn);
                  buf += (bn*bn)*sizeof(dense_block_type);
                  
                  setMatrixOfBlocks(hbr, n, n, _mb);
                  attachBaseBuffer(hbr, (mat_value_type*)buf, 1, n);
                  
                  // herk
                  Herk<Uplo::Upper,Trans::ConjTranspose,Algo::ByBlocks>
                    ::invoke(_sched, member, -1.0, htr, 0.0, hbr);
                  
                  for (ordinal_type k=0;k<(bm*bn);++k) htr[k].set_future();
                }
                for (ordinal_type k=0;k<(bm*bm);++k) htl[k].set_future();
              }
              
              _state = 2;
              if (bn > 0) {
                const size_t bn2 = bn*(bn+1)/2, depsize = bn2*sizeof(future_type);
                future_type *dep = (future_type*)_sched.memory()->allocate(depsize);
                TACHO_TEST_FOR_ABORT(dep == NULL, "sched memory pool allocation fails");
                clear((char*)dep, depsize);
                
                for (ordinal_type j=0,k=0;j<bn;++j)
                  for (ordinal_type i=0;i<=j;++i,++k) {
                    dep[k] = hbr(i,j).future();
                    hbr(i,j).set_future();
                  }
                Kokkos::respawn(this, Kokkos::when_all(dep, bn2), Kokkos::TaskPriority::Regular);
                
                for (ordinal_type k=0;k<static_cast<ordinal_type>(bn2);++k) (dep+k)->~future_type();
                _sched.memory()->deallocate((void*)dep, depsize);
              } else {
                if (_bufsize)
                  _bufpool.deallocate(_buf, _bufsize);
              }
            }
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
                                            TaskFunctor_FactorizeCholByBlocks
                                            (_sched, _bufpool, _info, _s.children[i], _mb));
                TACHO_TEST_FOR_ABORT(f.is_null(), "task allocation fails");
                dep[i] = f;
              }
            } else {
              for (ordinal_type i=0;i<_s.nchildren;++i) {
                auto f = Kokkos::task_spawn(Kokkos::TaskSingle(_sched, Kokkos::TaskPriority::Regular),
                                            TaskFunctor_FactorizeChol<mat_value_type,exec_space>
                                            (_sched, _bufpool, _info, _s.children[i]));
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
