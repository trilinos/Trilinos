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

#include "TachoExp_Gemm.hpp"
#include "TachoExp_Gemm_External.hpp"

#include "TachoExp_SupernodeInfo.hpp"

#include "TachoExp_CholSupernodes.hpp"
#include "TachoExp_CholSupernodes_Serial.hpp"
//#include "TachoExp_CholSupernodes_ByBlocks.hpp"

#include "TachoExp_TaskFunctor_CholSupernodes.hpp"
//#include "TachoExp_TaskFunctor_CholSupernodes_ByBlocks.hpp"

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

      typedef typename supernode_info_type_host::value_type_matrix value_type_matrix_host;
      typedef typename supernode_info_type_device::value_type_matrix value_type_matrix_device;

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

      // temp : schur arrays (will be replaced with a memory pool)
      size_type_array_host _super_schur_ptr;
      value_type_array_host _super_schur_buf;

    private:

      inline
      void
      recursiveSerialChol(const sched_type_host &sched,
                          const supernode_info_type_host &info,
                          const ordinal_type sid, 
                          const size_type bufsize,
                          /* */ void *buf) {
        // recursion
        const ordinal_type 
          ibeg = info.stree_ptr(sid), 
          iend = info.stree_ptr(sid+1);

        for (ordinal_type i=ibeg;i<iend;++i)
          recursiveSerialChol(sched, 
                              info, info.stree_children(i), 
                              bufsize, buf);

        {
          // dummy member
          const ordinal_type member = 0;

          CholSupernodes<Algo::Workflow::Serial>
            ::factorize(sched, member,
                        info, sid, 
                        bufsize, buf);

          ordinal_type m, n; info.getSuperPanelSize(sid, m, n);          
          UnmanagedViewType<value_type_matrix_host> ABR((value_type*)buf, n-m, n-m);

          CholSupernodes<Algo::Workflow::Serial>          
            ::update(sched, member,
                     info, ABR, sid, 
                     bufsize - ABR.span()*sizeof(value_type), 
                     (void*)((value_type*)buf + ABR.span()));
        }
      }

      inline
      void
      recursiveSerialSolveLower(const sched_type_host &sched,
                                const supernode_info_type_host &info,
                                const ordinal_type sid, 
                                const size_type bufsize,
                                /* */ void *buf) {
        // recursion
        const ordinal_type 
          ibeg = info.stree_ptr(sid), 
          iend = info.stree_ptr(sid+1);
        
        for (ordinal_type i=ibeg;i<iend;++i)
          recursiveSerialSolveLower(sched, 
                                    info, info.stree_children(i), 
                                    bufsize, buf);

        {
          // dummy member
          const ordinal_type member = 0;

          CholSupernodes<Algo::Workflow::Serial>
            ::solve_lower(sched, member,
                          info, sid, 
                          bufsize, buf);
          
          const ordinal_type nrhs = info.x.dimension_1();
          ordinal_type m, n; info.getSuperPanelSize(sid, m, n);          
          UnmanagedViewType<value_type_matrix_host> xB((value_type*)buf, n-m, nrhs);

          CholSupernodes<Algo::Workflow::Serial>          
            ::update_solve_lower(sched, member,
                                 info, xB, sid);
        }
      }

      inline
      void
      recursiveSerialSolveUpper(const sched_type_host &sched,
                                const supernode_info_type_host &info,
                                const ordinal_type sid, 
                                const size_type bufsize,
                                /* */ void *buf) {
        {
          // dummy member
          const ordinal_type member = 0;

          const ordinal_type nrhs = info.x.dimension_1();
          ordinal_type m, n; info.getSuperPanelSize(sid, m, n);          

          UnmanagedViewType<value_type_matrix_host> xB((value_type*)buf, n-m, nrhs);
          TACHO_TEST_FOR_ABORT(xB.span()*sizeof(value_type) > bufsize,
                               "xB uses greater than bufsize");

          CholSupernodes<Algo::Workflow::Serial>          
            ::update_solve_upper(sched, member,
                                 info, xB, sid);

          CholSupernodes<Algo::Workflow::Serial>
            ::solve_upper(sched, member,
                          info, xB, sid);          
        }

        // recursion
        const ordinal_type 
          ibeg = info.stree_ptr(sid), 
          iend = info.stree_ptr(sid+1);
        
        for (ordinal_type i=ibeg;i<iend;++i)
          recursiveSerialSolveUpper(sched, 
                                    info, info.stree_children(i), 
                                    bufsize, buf);
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

          /// recursive tree traversal
          const ordinal_type nroots = _stree_roots.dimension_0();
          for (ordinal_type i=0;i<nroots;++i)
            recursiveSerialChol(sched, 
                                info, _stree_roots(i), 
                                work.span()*sizeof(value_type), work.data());
        }
      }

      inline
      void
      solveCholesky_Serial(const value_type_matrix_host &x,
                           const value_type_matrix_host &b) {
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

          info.super_panel_ptr        = _super_panel_ptr;
          info.super_panel_buf        = _super_panel_buf;

          info.x                      = x;
        }
        
        // copy b -> x
        TACHO_TEST_FOR_EXCEPTION(x.dimension_0() != b.dimension_0() ||
                                 x.dimension_1() != b.dimension_1(), std::logic_error,
                                 "supernode data structure is not allocated");
        if (x.data() != b.data())
          Kokkos::deep_copy(x, b);
        
        {
          /// maximum workspace size for serial run
          const size_type worksize = sqrt(info.computeWorkspaceSerialChol())*x.dimension_1() + _m;
          value_type_array_host work("work", worksize);

          sched_type_host sched;

          /// recursive tree traversal
          const ordinal_type nroots = _stree_roots.dimension_0();
          for (ordinal_type i=0;i<nroots;++i)
            recursiveSerialSolveLower(sched, 
                                      info, _stree_roots(i), 
                                      work.span()*sizeof(value_type), work.data());

          for (ordinal_type i=0;i<nroots;++i)
            recursiveSerialSolveUpper(sched, 
                                      info, _stree_roots(i), 
                                      work.span()*sizeof(value_type), work.data());
        }
      }

      inline
      void
      factorizeCholesky_Parallel() {
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

          /// allocate schur complement workspace (this uses maximum preallocated workspace)
          info.allocateWorkspacePerSupernode(_super_schur_ptr, _super_schur_buf, iwork);

          info.super_schur_ptr = _super_schur_ptr;
          info.super_schur_buf = _super_schur_buf;
        }

        {
          typedef typename sched_type_host::memory_space memory_space;
          typedef TaskFunctor_CholSupernodes<value_type,host_exec_space> chol_supernode_functor_type;
          
          sched_type_host sched;
          {
            int nchildren = 0;
            const int nsupernodes = _stree_ptr.dimension_0() - 1;
            for (ordinal_type sid=0;sid<nsupernodes;++sid)
              nchildren = max(nchildren, _stree_ptr(sid+1) - _stree_ptr(sid));

            const size_type max_dep_future_size 
              = nchildren > MaxDependenceSize ? sizeof(Kokkos::Future<int,host_exec_space>)*nchildren : 16;
            
            const size_type max_functor_size = sizeof(chol_supernode_functor_type);
            const size_type estimate_max_numtasks = _blk_super_panel_colidx.dimension_0();
            
            const ordinal_type
              task_queue_capacity = max(estimate_max_numtasks,128)*max_functor_size, 
              min_block_size  = min(max_dep_future_size,max_functor_size),
              max_block_size  = max(max_dep_future_size,max_functor_size),
              num_superblock  = 32,
              superblock_size = task_queue_capacity/num_superblock;
            
            sched = sched_type_host(memory_space(),
                                    task_queue_capacity,
                                    min_block_size,
                                    max_block_size,
                                    superblock_size);
          }

          const ordinal_type nroots = _stree_roots.dimension_0();
          for (ordinal_type i=0;i<nroots;++i)
            Kokkos::host_spawn(Kokkos::TaskSingle(sched, Kokkos::TaskPriority::High),
                               chol_supernode_functor_type(sched, info, _stree_roots(i)));
          Kokkos::wait(sched);
        }
      }

      ///
      /// utility
      ///
      inline
      CrsMatrixBase<value_type,host_exec_space>
      exportFactorsToCrsMatrixBase(const bool replace_value_with_one = false) {
        /// this only avail after factorization is done
        TACHO_TEST_FOR_EXCEPTION(_super_panel_ptr.dimension_0() == 0, std::logic_error,
                                 "supernode data structure is not allocated");
        
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
          
          info.super_panel_ptr        = _super_panel_ptr;
          info.super_panel_buf        = _super_panel_buf;
        }
        return info.createCrsMatrixBase(replace_value_with_one);
      }

      ///
      /// host only input (value array can be rewritten in the same sparse structure)
      ///
      // inline
      // void
      // factorizeCholesky_ParallelByBlocks(const ordinal_type mb) {
      //   /// supernode info
      //   supernode_info_type_host info;
      //   {
      //     /// symbolic input
      //     info.supernodes             = _supernodes;
      //     info.gid_super_panel_ptr    = _gid_super_panel_ptr;
      //     info.gid_super_panel_colidx = _gid_super_panel_colidx;
          
      //     info.sid_super_panel_ptr    = _sid_super_panel_ptr;
      //     info.sid_super_panel_colidx = _sid_super_panel_colidx;
      //     info.blk_super_panel_colidx = _blk_super_panel_colidx;
          
      //     info.stree_ptr              = _stree_ptr;
      //     info.stree_children         = _stree_children;
      //   }
        
      //   {
      //     /// factor allocation and copy the matrix
      //     ordinal_type_array_host iwork("work", _m+1);

      //     info.allocateSuperPanels(_super_panel_ptr, _super_panel_buf, iwork);
          
      //     info.super_panel_ptr = _super_panel_ptr;
      //     info.super_panel_buf = _super_panel_buf;

      //     info.copySparseToSuperPanels(_ap, _aj, _ax, _perm, _peri, iwork);
      //   }

      //   {
      //     typedef typename sched_type_host::memory_space memory_space;
      //     typedef TaskFunctor_CholSupernodes_ByBlocks<value_type,host_exec_space> chol_supernode_byblocks_functor_type;

      //     const ordinal_type 
      //       small = 32, 
      //       blksize = (mb < small ? small : mb), 
      //       max_future_array_size = info.computeWorkspaceCholByBlocks(blksize)*sizeof(Kokkos::Future<int,exec_space>);

      //     sched_type_host sched;
      //     {
      //       const size_type max_functor_size = sizeof(chol_supernode_byblocks_functor_type);
      //       const size_type estimate_max_numtasks = _blk_super_panel_colidx.dimension_0();
            
      //       const ordinal_type
      //         task_queue_capacity = max(estimate_max_numtasks,128)*max_functor_size, 
      //         min_block_size  = min(max_future_array_size,max_functor_size),
      //         max_block_size  = max(max_future_array_size,max_functor_size),
      //         num_superblock  = 32,
      //         superblock_size = task_queue_capacity/num_superblock;
            
      //       sched = sched_type_host(memory_space(),
      //                               task_queue_capacity,
      //                               min_block_size,
      //                               max_block_size,
      //                               superblock_size);
      //     }

      //     memory_pool_type_host bufpool;
      //     {
      //       const size_type
      //         min_block_size  = small*small*sizeof(value_type),
      //         max_block_size  = blksize*blksize*sizeof(value_type),
      //         superblock_size = info.computeWorkspaceSerialChol()*sizeof(value_type),
      //         num_superblock  = 4,
      //         pool_memory_capacity = num_superblock*superblock_size;
            
      //       bufpool = memory_pool_type_host(memory_space(), 
      //                                       pool_memory_capacity,
      //                                       min_block_size,
      //                                       max_block_size,
      //                                       superblock_size);
      //     }
          
      //     const ordinal_type nroots = _stree_roots.dimension_0();
      //     for (ordinal_type i=0;i<nroots;++i)
      //       Kokkos::host_spawn(Kokkos::TaskSingle(sched, Kokkos::TaskPriority::High),
      //                          chol_supernode_byblocks_functor_type(sched, bufpool, info, _stree_roots(i), -1));
      //     Kokkos::wait(sched);

      //     std::cout << " -- sched memory pool -- \n";
      //     sched.memory()->print_state(std::cout);

      //     std::cout << " -- buffer memory pool -- \n";
      //     bufpool.memory()->print_state(std::cout);          
      //   }
      // }

      
    };

  }
}
#endif
