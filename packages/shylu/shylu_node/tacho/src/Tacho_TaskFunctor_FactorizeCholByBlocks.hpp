#ifndef __TACHO_TASKFUNCTOR_FACTORIZE_CHOLBYBLOCKS_HPP__
#define __TACHO_TASKFUNCTOR_FACTORIZE_CHOLBYBLOCKS_HPP__

/// \file Tacho_TaskFunctor_FactorizeCholByBlocks.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

#include "Tacho_Chol_ByBlocks.hpp"
#include "Tacho_Trsm_ByBlocks.hpp"
#include "Tacho_Herk_ByBlocks.hpp"
#include "Tacho_Gemm_ByBlocks.hpp"

#include "Tacho_CholSupernodes.hpp"
#include "Tacho_CholSupernodes_Serial.hpp"

#include "Tacho_TaskFunctor_FactorizeChol.hpp"

namespace Tacho {

    template<typename MatValueType, typename ExecSpace>
    struct TaskFunctor_FactorizeCholByBlocks {
    public:
      typedef ExecSpace exec_space;

      typedef Kokkos::TaskScheduler<exec_space> scheduler_type;
      typedef typename scheduler_type::member_type member_type;

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
      scheduler_type _sched;
      memory_pool_type _bufpool;

      supernode_info_type _info;
      ordinal_type _sid;

      //supernode_type _s;

      ordinal_type _state;

      ordinal_type _mb;

      // abr and matrix of blocks (deallocated after tasks are completed)
      char *_buf;
      size_t _bufsize;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_FactorizeCholByBlocks() = delete;

      KOKKOS_INLINE_FUNCTION
      TaskFunctor_FactorizeCholByBlocks(const scheduler_type &sched,
                                        const memory_pool_type &bufpool,
                                        const supernode_info_type &info,
                                        const ordinal_type sid,
                                        const ordinal_type mb)
        : _sched(sched),
          _bufpool(bufpool),
          _info(info),
          _sid(sid),
          //_s(info.supernodes(sid)),
          _state(0),
          _mb(mb),
          _buf(NULL),
          _bufsize(0) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(member_type &member, value_type &r_val) {
        const auto &_s = _info.supernodes(_sid);
        constexpr ordinal_type done = 3;
        TACHO_TEST_FOR_ABORT(_state == done, "dead lock");
            
        Kokkos::single(Kokkos::PerTeam(member), [&] () {
            if (_s.nchildren == 0 && _state == 0) _state = 1;
          });

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
          
          Kokkos::single(Kokkos::PerTeam(member), [&]() {
              if (_bufsize) {
                _bufpool.deallocate(_buf, _bufsize);
                _buf = NULL; _bufsize = 0;
              }
              _state = done;
            });
          break;
        }
        case 1: {
          ///
          /// matrix parallelism
          ///
          
          // get supernode panel pointer
          mat_value_type* ptr = _s.buf;

          // panel is divided into diagonal and interface block (i.e., ATL and ATR)
          const ordinal_type m = _s.m, n = _s.n - _s.m;

          // block matrix size
          const ordinal_type bm = m/_mb + (m%_mb > 0), bn = n/_mb + (n%_mb > 0);

          // allocation for matrix of blocks
          const size_t
            team_size = member.team_size(),
            matrix_of_blocks_bufsize = (bm*bm + bm*bn + bn*bn)*sizeof(dense_block_type),
            abr_bufsize = (n*n + _info.max_schur_size*team_size)*sizeof(mat_value_type),
            bufsize = matrix_of_blocks_bufsize + abr_bufsize;
          
          char *buf = NULL;
          Kokkos::single(Kokkos::PerTeam(member), [&](char *&val) {
              val = (char*)_bufpool.allocate(bufsize);
              if (bufsize) {
                if (val == NULL) {
                  Kokkos::respawn(this, _sched, Kokkos::TaskPriority::Low);                    
                } else {
                  // assign buf information for deallocation
                  _buf = buf; _bufsize = bufsize;
                }
              } else {
                _buf = NULL; _bufsize = 0;
              }
            }, buf);
          
          if (buf != NULL && bufsize > 0) {
            clear(member, (char*)buf, matrix_of_blocks_bufsize);              
            
            dense_matrix_of_blocks_type hbr;
            
            // m and n are available, then factorize the supernode block
            if (m > 0) {
              dense_matrix_of_blocks_type htl((dense_block_type*)buf, bm, bm);
              buf += (bm*bm)*sizeof(dense_block_type);
                
              setMatrixOfBlocks(member, htl, m, m, _mb);
              attachBaseBuffer(member, htl, ptr, 1, m); ptr += m*m;
                
              // chol
              Chol<Uplo::Upper,Algo::ByBlocks>::invoke(_sched, member, htl);
                
              if (n > 0) {
                dense_matrix_of_blocks_type htr((dense_block_type*)buf, bm, bn);
                buf += (bm*bn)*sizeof(dense_block_type);
                  
                setMatrixOfBlocks(member, htr, m, n, _mb);
                attachBaseBuffer(member, htr, ptr, 1, m);
                  
                // trsm
                Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,Algo::ByBlocks>
                  ::invoke(_sched, member, Diag::NonUnit(), 1.0, htl, htr);
                
                hbr = dense_matrix_of_blocks_type((dense_block_type*)buf, bn, bn);
                buf += (bn*bn)*sizeof(dense_block_type);
                  
                setMatrixOfBlocks(member, hbr, n, n, _mb);
                attachBaseBuffer(member, hbr, (mat_value_type*)buf, 1, n);
                
                // herk
                Herk<Uplo::Upper,Trans::ConjTranspose,Algo::ByBlocks>
                  ::invoke(_sched, member, -1.0, htr, 0.0, hbr);

                clearFutureOfBlocks(member, htr);
              }            
              clearFutureOfBlocks(member, htl);
            }
              
            _state = 2;
            if (bn > 0) {
              const size_t bn2 = bn*(bn+1)/2, depsize = bn2*sizeof(future_type);

              future_type *dep = NULL;
              Kokkos::single(Kokkos::PerTeam(member), [&](future_type *&val) {                  
                  val = (future_type*)_sched.memory()->allocate(depsize);
                }, dep);
              TACHO_TEST_FOR_ABORT(dep == NULL, "sched memory pool allocation fails");
              clear(member, (char*)dep, depsize);

              Kokkos::parallel_for(Kokkos::TeamThreadRange(member,bn),[&](const ordinal_type &j) {
                  Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,j+1),[&](const ordinal_type &i) {
                      dep[j*(j+1)/2 + i] = hbr(i,j).future();
                      hbr(i,j).set_future();
                    });
                });
              Kokkos::single(Kokkos::PerTeam(member), [&]() {
                  Kokkos::respawn(this, Kokkos::when_all(dep, bn2), Kokkos::TaskPriority::Regular);
                  for (ordinal_type k=0;k<static_cast<ordinal_type>(bn2);++k) (dep+k)->~future_type();
                  _sched.memory()->deallocate((void*)dep, depsize);
                });
            } else {
              Kokkos::single(Kokkos::PerTeam(member), [&]() {
                  if (_bufsize) {
                    _bufpool.deallocate(_buf, _bufsize);
                    _buf = NULL; _bufsize = 0;
                  }
                });
            }
          }
          break;
        }
        case 0: {
          ///
          /// tree parallelism
          ///          
          // spawn children tasks and this (their parent) depends on the children tasks
          Kokkos::single(Kokkos::PerTeam(member), [&]() {          
              const bool use_byblocks = (_mb*1.5 < _s.max_decendant_supernode_size);
              future_type *dep = NULL, depbuf[MaxDependenceSize];
              size_t depbuf_size = _s.nchildren > MaxDependenceSize ? _s.nchildren*sizeof(future_type) : 0;
              if (depbuf_size) {
                dep = (future_type*)_sched.memory()->allocate(depbuf_size);
                clear((char*)dep, depbuf_size);
              } else {
                dep = &depbuf[0];
              }

              if (dep == NULL) {
                Kokkos::respawn(this, _sched, Kokkos::TaskPriority::Regular);
              } else {
                if (use_byblocks) {
                  for (ordinal_type i=0;i<_s.nchildren;++i) {
                    auto f = Kokkos::task_spawn(Kokkos::TaskTeam(_sched, Kokkos::TaskPriority::Regular),
                                                TaskFunctor_FactorizeCholByBlocks
                                                (_sched, _bufpool, _info, _s.children[i], _mb));
                    TACHO_TEST_FOR_ABORT(f.is_null(), "task allocation fails");
                    dep[i] = f;
                  }
                } else {
                  for (ordinal_type i=0;i<_s.nchildren;++i) {
                    auto f = Kokkos::task_spawn(Kokkos::TaskTeam(_sched, Kokkos::TaskPriority::Regular),
                                                TaskFunctor_FactorizeChol<mat_value_type,exec_space>
                                                (_sched, _bufpool, _info, _s.children[i]));
                    TACHO_TEST_FOR_ABORT(f.is_null(), "task allocation fails");
                    dep[i] = f;
                  }
                }
                
                // respawn with updating state
                _state = 1;
                Kokkos::respawn(this, Kokkos::when_all(dep, _s.nchildren), Kokkos::TaskPriority::Regular);

                if (depbuf_size) {
                  for (ordinal_type i=0;i<_s.nchildren;++i) (dep+i)->~future_type();
                  _sched.memory()->deallocate(dep, depbuf_size);
                }
              }
            });
          break;
        }
        }
      }
    };

}

#endif
