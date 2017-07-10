#ifndef __TACHOEXP_TASKFUNCTOR_CHOLESKY_SUPERNODES_BYBLOCKS_HPP__
#define __TACHOEXP_TASKFUNCTOR_CHOLESKY_SUPERNODES_BYBLOCKS_HPP__

/// \file TachoExp_TaskFunctor_FactorizeSupernodes.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Util.hpp"

#include "TachoExp_Chol_ByBlocks.hpp"
#include "TachoExp_Trsm_ByBlocks.hpp"
#include "TachoExp_Herk_ByBlocks.hpp"
#include "TachoExp_Gemm_ByBlocks.hpp"

#include "TachoExp_CholSupernodes.hpp"
#include "TachoExp_CholSupernodes_Serial.hpp"

namespace Tacho {
  
  namespace Experimental {
    
    template<typename MatValueType, typename ExecSpace>
    struct TaskFunctor_CholSupernodes_ByBlocks {
    public:
      typedef ExecSpace exec_space;

      typedef Kokkos::TaskScheduler<exec_space> sched_type;
      typedef typename sched_type::member_type member_type;

      typedef Kokkos::MemoryPool<exec_space> memory_pool_type;

      typedef int value_type; // functor return type
      typedef Kokkos::Future<int,exec_space> future_type;

      typedef MatValueType mat_value_type; // matrix value type

      typedef SupernodeInfo<mat_value_type,exec_space> supernode_info_type;
      typedef typename supernode_info_type::value_type_matrix mat_value_type_matrix;
      typedef typename supernode_info_type::dense_block_type dense_block_type;
      typedef typename supernode_info_type::dense_matrix_of_blocks_type dense_matrix_of_blocks_type;

    private:
      sched_type _sched;
      memory_pool_type _bufpool;
      
      supernode_info_type _info;
      ordinal_type _sid;

      // respawn control
      ordinal_type _state;

      // blocksize
      ordinal_type _mb;

      // abr and matrix of blocks
      char *_buf;
      size_type _bufsize;

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
          _mb(mb),
          _buf(NULL),
          _bufsize(0) {}
      
      KOKKOS_INLINE_FUNCTION
      void operator()(member_type &member, value_type &r_val) {
        if (get_team_rank(member) == 0) {
          TACHO_TEST_FOR_ABORT(_state > 2, "dead lock");

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
            ordinal_type pm, pn; _info.getSuperPanelSize(_sid, pm, pn);

            const ordinal_type 
              m = pm, n = pn - pm,
              bm = m/_mb + (m%_mb > 0), bn = n/_mb + (n%_mb > 0);
            
            const size_type 
              matrix_of_blocks_bufsize = (bm*bm + bm*bn + bn*bn)*sizeof(dense_block_type);

            char *buf = _buf + matrix_of_blocks_bufsize;
            const size_type bufsize = _bufsize - matrix_of_blocks_bufsize;

            UnmanagedViewType<mat_value_type_matrix> ABR((mat_value_type*)buf, n, n);
            CholSupernodes<Algo::Workflow::Serial>
              ::update(_sched, member,
                       _info, ABR, _sid,
                       (bufsize - ABR.span()*sizeof(mat_value_type)),
                       (void*)((mat_value_type*)buf + ABR.span()));

            if (_bufsize) 
              _bufpool.deallocate(_buf, _bufsize);

            _state = 3; // done
            break;
          }
          case 1: {
            ///
            /// matrix parallelism
            ///

            // get supernode panel pointer
            mat_value_type *ptr = _info.getSuperPanelPtr(_sid);

            // characterize the panel size
            ordinal_type pm, pn;
            _info.getSuperPanelSize(_sid, pm, pn);

            // panel is divided into diagonal and interface block (i.e., ATL and ATR)
            const ordinal_type m = pm, n = pn - pm;
            
            // block matrix size
            const ordinal_type bm = m/_mb + (m%_mb > 0), bn = n/_mb + (n%_mb > 0);

            // allocation for matrix of blocks
            const size_type 
              matrix_of_blocks_bufsize = (bm*bm + bm*bn + bn*bn)*sizeof(dense_block_type),
              abr_bufsize = (n*n + _info.max_schur_size)*sizeof(mat_value_type),
              bufsize = matrix_of_blocks_bufsize + abr_bufsize;
            
            char *buf = (char*)_bufpool.allocate(bufsize); 
            TACHO_TEST_FOR_ABORT(buf == NULL && bufsize != 0, "bufpool allocation fails");
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
              for (ordinal_type k=0;k<(bm*bm);++k) htl[k].set_future(future_type());
            }

            _state = 2;              
            if (bn > 0) {
              const size_type bn2 = bn*(bn+1)/2, depsize = bn2*sizeof(future_type);
              future_type *dep = (future_type*)_sched.memory()->allocate(depsize);
              TACHO_TEST_FOR_ABORT(dep == NULL, "sched memory pool allocation fails");
              clear((char*)dep, depsize);
              
              ordinal_type k = 0;
              for (ordinal_type j=0;j<bn;++j)
                for (ordinal_type i=0;i<=j;++i,++k) {
                  dep[k] = hbr(i,j).future();
                  hbr(i,j).set_future();
              }
              Kokkos::respawn(this, Kokkos::when_all(dep, bn2), Kokkos::TaskPriority::Regular);
              
              for (ordinal_type k=0;k<bn2;++k) (dep+k)->~future_type();
              _sched.memory()->deallocate((void*)dep, depsize);
            } else {
              Kokkos::respawn(this, future_type(), Kokkos::TaskPriority::Regular);
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
              clear((char*)dep, depsize);
            }
            
            // spawn children tasks and this (their parent) depends on the children tasks
            for (ordinal_type i=0;i<isize;++i) {
              // the first child has a higher priority
              const auto priority = (i ? Kokkos::TaskPriority::Low : Kokkos::TaskPriority::High);
              
              const ordinal_type child = _info.stree_children(i+ibeg);
              auto f = Kokkos::task_spawn(Kokkos::TaskSingle(_sched, priority),
                                          TaskFunctor_CholSupernodes_ByBlocks(_sched, _bufpool, _info, child, _mb));
              TACHO_TEST_FOR_ABORT(f.is_null(), "task allocation fails");
              dep[i] = f;
            }
            
            // respawn with updating state
            _state = 1;
            Kokkos::respawn(this, Kokkos::when_all(dep, isize), Kokkos::TaskPriority::Regular);
              
            // deallocate dependence array
            if (depsize) {
              // manually reset future to decrease the reference count
              for (ordinal_type i=0;i<isize;++i) (dep+i)->~future_type();
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

            
            

