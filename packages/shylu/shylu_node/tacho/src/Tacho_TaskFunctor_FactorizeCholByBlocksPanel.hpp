#ifndef __TACHO_TASKFUNCTOR_FACTORIZE_CHOLBYBLOCKS_PANEL_HPP__
#define __TACHO_TASKFUNCTOR_FACTORIZE_CHOLBYBLOCKS_PANEL_HPP__

/// \file Tacho_TaskFunctor_FactorizeCholByBlocks.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

#include "Tacho_Chol_ByBlocks.hpp"
#include "Tacho_Trsm_ByBlocks.hpp"
#include "Tacho_Herk_ByBlocks.hpp"
#include "Tacho_Gemm_ByBlocks.hpp"

#include "Tacho_CholSupernodes.hpp"
#include "Tacho_CholSupernodes_SerialPanel.hpp"

#include "Tacho_TaskFunctor_FactorizeCholPanel.hpp"

namespace Tacho {

    template<typename MatValueType, typename ExecSpace>
    struct TaskFunctor_FactorizeCholByBlocksPanel {
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

      ordinal_type _mb, _nb, _offn;

    public:
      KOKKOS_INLINE_FUNCTION
      TaskFunctor_FactorizeCholByBlocksPanel() = delete;

      KOKKOS_INLINE_FUNCTION
      TaskFunctor_FactorizeCholByBlocksPanel(const scheduler_type &sched,
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
      //_s(info.supernodes(sid)),
        _state(state),
        _mb(mb),
        _nb(nb),
        _offn(offn) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(member_type &member, value_type &r_val) {
        const auto &_s = _info.supernodes(_sid);
        constexpr ordinal_type done = 5;
        TACHO_TEST_FOR_ABORT(_state == done, "dead lock");
        
        Kokkos::single(Kokkos::PerTeam(member), [&]() {
            if (_s.nchildren == 0 && _state == 0) _state = 1;
          });

        switch (_state) {
        case 3: {
          ///
          /// partial update and its parent
          ///
          const ordinal_type n = _s.n - _s.m;
          const ordinal_type nb = min(_nb, n - _offn);
          const ordinal_type nn = _offn + nb;
          const size_type 
            team_size = member.team_size(),
            bufsize = n*team_size*sizeof(mat_value_type) + nn*nb*sizeof(mat_value_type);  

          mat_value_type *buf = NULL; 
          Kokkos::single(Kokkos::PerTeam(member), [&](mat_value_type *&val) {
              val = bufsize > 0 ? (mat_value_type*)_bufpool.allocate(bufsize) : NULL;  
              if (val == NULL && bufsize) 
                Kokkos::respawn(this, _sched, Kokkos::TaskPriority::Low);
            }, buf);
          
          if (buf != NULL && bufsize > 0) {
            CholSupernodes<Algo::Workflow::SerialPanel>
              ::update(_sched, member, 
                       _info, _offn, nb, _sid, bufsize, (void*)buf);
            Kokkos::single(Kokkos::PerTeam(member), [&]() {
                _bufpool.deallocate(buf, bufsize); 
              });
          }
          break;
        }
        case 2: {
          ///
          /// update to its parent
          ///
          const ordinal_type n = _s.n - _s.m;
          const ordinal_type nb = _nb;
          const ordinal_type bn = n/nb + (n%nb > 0);
          const ordinal_type depsize = bn*sizeof(future_type);

          future_type *dep = NULL;
          Kokkos::single(Kokkos::PerTeam(member), [&](future_type *&val) {                  
              val = (future_type*)_sched.memory()->allocate(depsize);  
            }, dep);
          TACHO_TEST_FOR_ABORT(dep == NULL, "sched memory pool allocation fails"); 
          clear(member, (char*)dep, depsize);  
          
          Kokkos::single(Kokkos::PerTeam(member), [&]() {
              const ordinal_type state = 3;
              for (ordinal_type i=0;i<bn;++i) {
                auto f = Kokkos::task_spawn(Kokkos::TaskTeam(_sched, Kokkos::TaskPriority::Regular),
                                            TaskFunctor_FactorizeCholByBlocksPanel
                                            (_sched, _bufpool, _info, _sid, _mb, nb, state, i*nb));
                TACHO_TEST_FOR_ABORT(f.is_null(), "task allocation fails");
                dep[i] = f;
              }
          
              _state = 4;
              Kokkos::respawn(this, Kokkos::when_all(dep, bn), Kokkos::TaskPriority::Regular);
              for (ordinal_type k=0;k<bn;++k) (dep+k)->~future_type();
              if (depsize) _sched.memory()->deallocate(dep, depsize);
            });
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

#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
          const ordinal_type mb = _mb;
#else
          const ordinal_type mnmax = max(m,n);
          const ordinal_type mb = mnmax < 160 ? 40 : (mnmax < 320 ? 80 : _mb);
#endif

          // block matrix size
          const ordinal_type bm = m/mb + (m%mb > 0), bn = n/mb + (n%mb > 0);
          
          // allocation for matrix of blocks
          const size_t bufsize = (bm*bm + bm*bn)*sizeof(dense_block_type);
          char *buf = NULL;
          Kokkos::single(Kokkos::PerTeam(member), [&](char *&val) {          
              val = (char*)_bufpool.allocate(bufsize);
              if (val == NULL && bufsize) 
                Kokkos::respawn(this, _sched, Kokkos::TaskPriority::Low);                  
            }, buf);

          if (buf != NULL && bufsize) {
            clear(member, (char*)buf, bufsize);

            dense_block_type *buf_ptr = (dense_block_type*)buf;

            // m and n are available, then factorize the supernode block
            dense_matrix_of_blocks_type htr;                
            if (m > 0) {
              dense_matrix_of_blocks_type htl(buf_ptr, bm, bm);
              buf_ptr += bm*bm;
              
              setMatrixOfBlocks(member, htl, m, m, mb);
              attachBaseBuffer(member, htl, ptr, 1, m); ptr += m*m;
              member.team_barrier();
              // chol
              Chol<Uplo::Upper,Algo::ByBlocks>::invoke(_sched, member, htl);

              if (n > 0) {
                htr = dense_matrix_of_blocks_type(buf_ptr, bm, bn);
                buf_ptr += bm*bn;
                
                setMatrixOfBlocks(member, htr, m, n, mb);
                attachBaseBuffer(member, htr, ptr, 1, m);
                  
                // trsm
                Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,Algo::ByBlocks>
                  ::invoke(_sched, member, Diag::NonUnit(), 1.0, htl, htr);                  
                member.team_barrier();
              }
              clearFutureOfBlocks(member, htl);
            }

            if (m > 0 && n > 0) {
              const size_t bmn = bm*bn, depsize = bm*bn*sizeof(future_type);
              future_type *dep = NULL;
              Kokkos::single(Kokkos::PerTeam(member), [&](future_type *&val) {                  
                  val = (future_type*)_sched.memory()->allocate(depsize);  
                }, dep);
              TACHO_TEST_FOR_ABORT(dep == NULL, "sched memory pool allocation fails"); 
              clear(member, (char*)dep, depsize);  
              // GCC compiler bug workaround              
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
              for (ordinal_type j=0,l=0;j<bn;++j)
                for (ordinal_type i=0;i<bm;++i,++l) {
                  dep[l] = htr(i,j).future();
                  htr(i,j).set_future();                  
                }
#else              
              Kokkos::parallel_for(Kokkos::TeamThreadRange(member,bn),[&](const ordinal_type &j) {
                  Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,bm),[&](const ordinal_type &i) {
                      dep[j*bm+i] = htr(i,j).future();
                      htr(i,j).set_future();
                    });
                });
#endif
              member.team_barrier();
              Kokkos::single(Kokkos::PerTeam(member), [&]() {              
                  _state = 2;
                  Kokkos::respawn(this, Kokkos::when_all(dep, bmn), Kokkos::TaskPriority::Regular); 

                  for (ordinal_type k=0;k<static_cast<ordinal_type>(bmn);++k) (dep+k)->~future_type();
                  _sched.memory()->deallocate(dep, depsize);   
                });
            } else {
              Kokkos::single(Kokkos::PerTeam(member), [&]() {                            
                  _state = 4;
                });
            }
          }

          Kokkos::single(Kokkos::PerTeam(member), [&]() {                                      
              if (bufsize) _bufpool.deallocate(buf, bufsize);
            });
          break;
        }
        case 0: {
          ///
          /// tree parallelism
          ///
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
                                                TaskFunctor_FactorizeCholByBlocksPanel
                                                (_sched, _bufpool, _info, _s.children[i], _mb, _nb));
                    TACHO_TEST_FOR_ABORT(f.is_null(), "task allocation fails");
                    dep[i] = f;
                  }
                } else {
                  for (ordinal_type i=0;i<_s.nchildren;++i) {
                    auto f = Kokkos::task_spawn(Kokkos::TaskTeam(_sched, Kokkos::TaskPriority::Regular),
                                                TaskFunctor_FactorizeCholPanel<mat_value_type,exec_space>
                                                (_sched, _bufpool, _info, _s.children[i], _nb));
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
