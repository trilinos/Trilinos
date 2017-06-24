#ifndef __TACHOEXP_NUMERIC_TOOLS_HPP__
#define __TACHOEXP_NUMERIC_TOOLS_HPP__

/// \file TachoExp_NumericTools.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Util.hpp"

#include "TachoExp_Chol.hpp"
#include "TachoExp_Chol_External.hpp"

#include "TachoExp_Trsm.hpp"
#include "TachoExp_Trsm_External.hpp"

#include "TachoExp_Herk.hpp"
#include "TachoExp_Herk_External.hpp"

#include "TachoExp_SupernodeInfo.hpp"

#include "TachoExp_CholSupernodes.hpp"
#include "TachoExp_CholSupernodes_Serial.hpp"

#include "TachoExp_TaskFunctor_CholSupernodes.hpp"

namespace Tacho {

  namespace Experimental {

    template<typename ValueType, typename ExecSpace>
    class NumericTools {
    public:
      typedef ValueType value_type;

      typedef Kokkos::DefaultHostExecutionSpace host_exec_space;
      typedef Kokkos::View<ordinal_type*,host_exec_space> ordinal_type_array_host;
      typedef Kokkos::View<size_type*,   host_exec_space> size_type_array_host;
      typedef Kokkos::View<value_type*,  host_exec_space> value_type_array_host;

      typedef ExecSpace device_exec_space;
      typedef Kokkos::View<ordinal_type*,device_exec_space> ordinal_type_array_device;
      typedef Kokkos::View<size_type*,   device_exec_space> size_type_array_device;
      typedef Kokkos::View<value_type*,  device_exec_space> value_type_array_device;

      typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;

      typedef SupernodeInfo<value_type,host_exec_space> supernode_info_type_host;
      typedef SupernodeInfo<value_type,device_exec_space> supernode_info_type_device;

      typedef Kokkos::TaskScheduler<host_exec_space> sched_type_host;
      typedef Kokkos::MemoryPool<host_exec_space> memory_pool_type_host;

    private:
      // matrix input
      ordinal_type _m;
      size_type_array_host _ap;
      ordinal_type_array_host _aj;
      value_type_array_host _ax;

      // graph ordering input
      ordinal_type_array_host _perm, _peri;

      // supernodes input
      ordinal_type _nsupernodes;
      ordinal_type_array_host _supernodes;
      
      // dof mapping to sparse matrix
      size_type_array_host _gid_super_panel_ptr;
      ordinal_type_array_host _gid_super_panel_colidx;

      // supernode map and panel size configuration
      size_type_array_host _sid_super_panel_ptr;
      ordinal_type_array_host _sid_super_panel_colidx, _blk_super_panel_colidx;

      // supernode tree
      ordinal_type_array_host _stree_parent;
      size_type_array_host _stree_ptr;
      ordinal_type_array_host _stree_children, _stree_roots;

      // output : factors
      size_type_array_host _super_panel_ptr;
      value_type_array_host _super_panel_buf;

    private:

      inline
      void
      recursiveSerialChol(const sched_type_host &sched,
                          const memory_pool_type_host &pool,
                          const supernode_info_type_host &info,
                          const ordinal_type sid, 
                          const ordinal_type sidpar,
                          const size_type bufsize,
                          /* */ void *buf) {
        // recursion
        const ordinal_type 
          ibeg = info.stree_ptr(sid), 
          iend = info.stree_ptr(sid+1);

        for (ordinal_type i=ibeg;i<iend;++i)
          recursiveSerialChol(sched, pool, 
                              info, info.stree_children(i), sid,
                              bufsize, buf);

        {
          const ordinal_type member = 0;
          UnmanagedViewType<typename supernode_info_type_host::value_type_matrix> ABR;
          CholSupernodes<Algo::Workflow::Serial>
            ::factorize(sched, pool, 
                        member,
                        info, sid, sidpar,
                        ABR,
                        bufsize, buf);
          
          TACHO_TEST_FOR_ABORT(ABR.span() && ABR.data() != buf, 
                               "memory pool is used in serial factorization");
          
          Kokkos::pair<size_type,size_type> srows(0,0);
          CholSupernodes<Algo::Workflow::Serial>          
            ::update(sched, pool, 
                     member,
                     info, ABR, srows, sid, 
                     bufsize - ABR.span()*sizeof(value_type), 
                     (void*)((value_type*)buf + ABR.span()));
        }
      }
      
    public:
      NumericTools() = default;
      NumericTools(const NumericTools &b) = default;

      ///
      /// construction (assume input matrix and symbolic are from host)
      ///
      NumericTools(// input matrix A
                   const ordinal_type m,
                   const size_type_array_host &ap,
                   const ordinal_type_array_host &aj,
                   const value_type_array_host &ax,
                   // input permutation
                   const ordinal_type_array_host &perm,
                   const ordinal_type_array_host &peri,
                   // supernodes
                   const ordinal_type nsupernodes,
                   const ordinal_type_array_host &supernodes,
                   const size_type_array_host &gid_super_panel_ptr,
                   const ordinal_type_array_host &gid_super_panel_colidx,
                   const size_type_array_host &sid_super_panel_ptr,
                   const ordinal_type_array_host &sid_super_panel_colidx,
                   const ordinal_type_array_host &blk_super_panel_colidx,
                   const ordinal_type_array_host &stree_parent,
                   const size_type_array_host &stree_ptr,
                   const ordinal_type_array_host &stree_children,
                   const ordinal_type_array_host &stree_roots)
        : _m(m),
          _ap(ap),
          _aj(aj),
          _ax(ax),
          _perm(perm),
          _peri(peri),
          _nsupernodes(nsupernodes),
          _supernodes(supernodes),
          _gid_super_panel_ptr(gid_super_panel_ptr),
          _gid_super_panel_colidx(gid_super_panel_colidx),
          _sid_super_panel_ptr(sid_super_panel_ptr),
          _sid_super_panel_colidx(sid_super_panel_colidx),
          _blk_super_panel_colidx(blk_super_panel_colidx),
          _stree_parent(stree_parent),
          _stree_ptr(stree_ptr),
          _stree_children(stree_children),
          _stree_roots(stree_roots) {}

      ///
      /// host only input (value array can be rewritten in the same sparse structure)
      ///
      inline
      void
      factorizeCholesky_Serial() {
        /// supernode info
        supernode_info_type_host info;
        {
          /// symbolic input
          info.supernodes             = _supernodes;
          info.gid_super_panel_ptr    = _gid_super_panel_ptr;
          info.gid_super_panel_colidx = _gid_super_panel_colidx;
          
          info.sid_super_panel_ptr    = _sid_super_panel_ptr;
          info.sid_super_panel_colidx = _sid_super_panel_colidx;
          info.blk_super_panel_colidx = _blk_super_panel_colidx;
          
          info.stree_ptr              = _stree_ptr;
          info.stree_children         = _stree_children;
        }

        {
          /// factor allocation and copy the matrix
          ordinal_type_array_host iwork("work", _m+1);

          info.allocateSuperPanels(_super_panel_ptr, _super_panel_buf, iwork);
          
          info.super_panel_ptr = _super_panel_ptr;
          info.super_panel_buf = _super_panel_buf;

          info.copySparseToSuperPanels(_ap, _aj, _ax, _perm, _peri, iwork);
        }

        {
          /// maximum workspace size for serial run
          const size_type worksize = info.computeWorkspaceSerialChol() + _m + 1;
          value_type_array_host work("work", worksize);

          sched_type_host sched;
          memory_pool_type_host pool;

          /// recursive tree traversal
          const ordinal_type nroots = _stree_roots.dimension_0();
          for (ordinal_type i=0;i<nroots;++i)
            recursiveSerialChol(sched, pool, 
                                info, _stree_roots(i), -1,
                                work.span()*sizeof(value_type), work.data());
        }
      }

      ///
      /// host only input (value array can be rewritten in the same sparse structure)
      ///
      inline
      void
      factorizeCholesky_Parallel() {
        /// supernode info
        supernode_info_type_host info;
        typename supernode_info_type_host::value_type_matrix_array abrs("abrs", _nsupernodes); 

        //typename supernode_info_type_host::future_type_array futures(Kokkos::ViewAllocateWithoutInitializing("futures"), _nsupernodes); 
        typename supernode_info_type_host::future_type_array futures("futures", _nsupernodes); 


        {
          /// symbolic input
          info.supernodes             = _supernodes;
          info.gid_super_panel_ptr    = _gid_super_panel_ptr;
          info.gid_super_panel_colidx = _gid_super_panel_colidx;
          
          info.sid_super_panel_ptr    = _sid_super_panel_ptr;
          info.sid_super_panel_colidx = _sid_super_panel_colidx;
          info.blk_super_panel_colidx = _blk_super_panel_colidx;
          
          info.stree_ptr              = _stree_ptr;
          info.stree_children         = _stree_children;

          info.supernodes_abr         = abrs;
          info.supernodes_future      = futures;
        }
        
        {
          /// factor allocation and copy the matrix
          ordinal_type_array_host iwork("work", _m+1);

          info.allocateSuperPanels(_super_panel_ptr, _super_panel_buf, iwork);
          
          info.super_panel_ptr = _super_panel_ptr;
          info.super_panel_buf = _super_panel_buf;

          info.copySparseToSuperPanels(_ap, _aj, _ax, _perm, _peri, iwork);
        }

        {
          typedef typename sched_type_host::memory_space memory_space;

          sched_type_host sched;
          {
            const size_type max_functor_size = ( sizeof(supernode_info_type_host) + 
                                                 sizeof(sched_type_host) + 
                                                 sizeof(memory_pool_type_host) + 128 );
            const size_type estimate_max_numtasks = _blk_super_panel_colidx.dimension_0();
            const size_type task_queue_capacity = max(estimate_max_numtasks*max_functor_size*1000,
                                                      2048*400);
            
            const ordinal_type 
              min_block_size  = max_functor_size*8,
              max_block_size  = max_functor_size*8,
              superblock_size = max_block_size*16;
            
            sched = sched_type_host(memory_space(),
                                    task_queue_capacity,
                                    min_block_size,
                                    max_block_size,
                                    superblock_size);
          }

          memory_pool_type_host pool;
          {
            const size_type worksize = info.computeWorkspaceSerialChol() + _m + 1;
            const size_type pool_memory_capacity = max(worksize*sizeof(value_type)*8, 4096);

            ordinal_type alpha = 10000;
            const size_type
              min_block_size  = 64*alpha,
              max_block_size  = 1024*alpha,
              superblock_size = 4096*alpha;
            
            pool = memory_pool_type_host(memory_space(), 
                                         pool_memory_capacity,
                                         min_block_size,
                                         max_block_size,
                                         superblock_size);
          }
      
          typedef TaskFunctor_CholSupernodes<value_type,host_exec_space> chol_supernode_functor_type;
          const ordinal_type nroots = _stree_roots.dimension_0();
          for (ordinal_type i=0;i<nroots;++i)
            Kokkos::host_spawn(Kokkos::TaskSingle(sched, Kokkos::TaskPriority::High),
                               chol_supernode_functor_type(sched, pool, info, _stree_roots(i), -1));
          Kokkos::wait(sched);          
        }
      }

      
      // template<typename ArgUplo>
      // inline
      // CrsMatrixBase<value_type,host_exec_space> 
      // Factors() {
      //   CrsMatrixBase<value_type,host_exec_space> r_val;
        
      //   // for (ordinal_type i=0;i<_m;++i)
      //   //   std::cout << _perm(i) << "  ";
      //   // std::cout << "\n";

      //   const value_type eps = std::numeric_limits<value_type>::epsilon() * 100;
      //   if (std::is_same<ArgUplo,Uplo::Upper>::value) {
      //     size_type_array_host count("count", _m);
      //     for (ordinal_type sid=0;sid<_nsupernodes;++sid) {
      //       ordinal_type m, n;
      //       getSuperPanelSize(sid, _supernodes, _sid_super_panel_ptr, _blk_super_panel_colidx, m, n);
      //       Kokkos::View<value_type**,Kokkos::LayoutLeft,
      //         typename value_type_array_host::execution_space,
      //         Kokkos::MemoryUnmanaged> tgt(&_super_panel_buf(_super_panel_ptr(sid)), m, n);          
            
      //       const ordinal_type soffset = _supernodes(sid);            
      //       for (ordinal_type i=0;i<m;++i) {
      //         const ordinal_type ii = i+soffset;
      //         for (ordinal_type j=i;j<n;++j) 
      //           if (std::abs(tgt(i,j)) > eps)
      //             ++count(ii);
      //       }
      //     }
          
      //     size_type_array_host ap("ap", _m+1);          
      //     for (ordinal_type i=0;i<_m;++i) 
      //       ap(i+1) = ap(i) + count(i);

      //     ordinal_type_array_host aj("aj", ap(_m));
      //     value_type_array_host ax("ax", ap(_m));

      //     Kokkos::deep_copy(count, 0);
          
      //     for (ordinal_type sid=0;sid<_nsupernodes;++sid) {
      //       ordinal_type m, n;
      //       getSuperPanelSize(sid, _supernodes, _sid_super_panel_ptr, _blk_super_panel_colidx, m, n);
      //       Kokkos::View<value_type**,Kokkos::LayoutLeft,
      //         typename value_type_array_host::execution_space,
      //         Kokkos::MemoryUnmanaged> tgt(&_super_panel_buf(_super_panel_ptr(sid)), m, n);          
            
      //       const ordinal_type soffset = _supernodes(sid);
      //       const ordinal_type goffset = _gid_super_panel_ptr(sid);
      //       for (ordinal_type i=0;i<m;++i) {
      //         const ordinal_type ii = (i+soffset);
      //         for (ordinal_type j=i;j<n;++j) { 
      //           if (std::abs(tgt(i,j)) > eps) {
      //             aj(ap(ii) + count(ii)) = (_gid_super_panel_colidx(j+goffset));
      //             ax(ap(ii) + count(ii)) = tgt(i,j);
      //             ++count(ii);
      //           }
      //         }
      //       }
      //     }
      //     r_val.setExternalMatrix(_m, _m, ap(_m), ap, aj, ax);
      //   } 
      //   else if (std::is_same<ArgUplo,Uplo::Lower>::value) {
      //   }

      //   return r_val;
      // }
      
    };

  }
}
#endif
