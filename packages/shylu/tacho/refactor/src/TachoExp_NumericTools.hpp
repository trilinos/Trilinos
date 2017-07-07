#ifndef __TACHOEXP_NUMERIC_TOOLS_HPP__
#define __TACHOEXP_NUMERIC_TOOLS_HPP__

/// \file TachoExp_NumericTools.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Util.hpp"

#include "TachoExp_CrsMatrixBase.hpp"
#include "TachoExp_DenseMatrixView.hpp"

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
      typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;

      ///
      /// host space typedefs
      ///
      typedef Kokkos::DefaultHostExecutionSpace host_exec_space;
      typedef SupernodeInfo<value_type,host_exec_space> supernode_info_type_host;
      typedef typename supernode_info_type_host::crs_matrix_type crs_matrix_type_host;
      typedef typename supernode_info_type_host::value_type_matrix value_type_matrix_host;

      typedef typename crs_matrix_type_host::ordinal_type_array ordinal_type_array_host;
      typedef typename crs_matrix_type_host::size_type_array size_type_array_host;
      typedef typename crs_matrix_type_host::value_type_array value_type_array_host;

      typedef Kokkos::TaskScheduler<host_exec_space> sched_type_host;
      typedef Kokkos::MemoryPool<host_exec_space> memory_pool_type_host;

      ///
      /// device space typedefs
      ///
      typedef ExecSpace device_exec_space;
      typedef SupernodeInfo<value_type,device_exec_space> supernode_info_type_device;
      typedef typename supernode_info_type_device::crs_matrix_type crs_matrix_type_device;
      typedef typename supernode_info_type_device::value_type_matrix value_type_matrix_device;

      typedef typename crs_matrix_type_device::ordinal_type_array ordinal_type_array_device;
      typedef typename crs_matrix_type_device::size_type_array size_type_array_device;
      typedef typename crs_matrix_type_device::value_type_array value_type_array_device;

      typedef Kokkos::TaskScheduler<device_exec_space> sched_type_device;
      typedef Kokkos::MemoryPool<device_exec_space> memory_pool_type_device;
    
    private:

      ///
      /// supernode data structure memory "managed"
      /// this holds all necessary connectivity data 
      ///

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

      ///
      /// supernode info: supernode data structure with "unamanged" view
      /// this is passed into computation algorithm without reference counting
      ///
      supernode_info_type_host _info;

      ///
      ///
      ///
      struct {
        double t_factor, t_solve, t_copy, t_extra;
        double m_used, m_peak;
      } stat;
      
    private:

      inline 
      void
      track_alloc(const double in) {
        stat.m_used += in;
        stat.m_peak  = max(stat.m_peak, stat.m_used);
      }

      inline 
      void
      track_free(const double out) {
        stat.m_used -= out;
      }

      inline 
      void
      reset_stat() {
        stat.t_factor = 0;
        stat.t_solve = 0;
        stat.t_copy = 0;
        stat.t_extra = 0;
        stat.m_used = 0;
        stat.m_peak = 0;
      }

      inline 
      void 
      print_stat_factor() {
        printf("  Time\n");
        printf("             time for copying A into U:         %10.6f s\n", stat.t_copy);
        printf("             time for numeric factorization:    %10.6f s\n", stat.t_factor);
        printf("             total time spent:                  %10.6f s\n", (stat.t_copy+stat.t_factor));
        printf("\n");
        printf("  Memory\n");
        printf("             memory used in factorization:      %10.2f MB\n", stat.m_used/1024/1024);
        printf("             peak memory used in factorization: %10.2f MB\n", stat.m_peak/1024/1024);
      }

      inline 
      void 
      print_stat_solve() {
        printf("  Time\n");
        printf("             time for extra work e.g.,copy rhs: %10.6f s\n", stat.t_extra);
        printf("             time for numeric solve:            %10.6f s\n", stat.t_solve);
        printf("             total time spent:                  %10.6f s\n", (stat.t_solve+stat.t_extra));
        printf("\n");
      }

      inline
      void
      recursiveSerialChol(const ordinal_type sid, 
                          const size_type bufsize,
                          /* */ void *buf) {
        // recursion
        {
          const ordinal_type 
            ibeg = _info.stree_ptr(sid), 
            iend = _info.stree_ptr(sid+1);
          
          for (ordinal_type i=ibeg;i<iend;++i)
            recursiveSerialChol(_info.stree_children(i), bufsize, buf);
        }

        {
          // dummy member
          const ordinal_type sched = 0, member = 0;

          ordinal_type m, n; _info.getSuperPanelSize(sid, m, n);          
          UnmanagedViewType<value_type_matrix_host> ABR((value_type*)buf, n-m, n-m);

          CholSupernodes<Algo::Workflow::Serial>
            ::factorize(sched, member, _info, ABR, sid);

          CholSupernodes<Algo::Workflow::Serial>          
            ::update(sched, member, _info, ABR, sid, 
                     bufsize - ABR.span()*sizeof(value_type), 
                     (void*)((value_type*)buf + ABR.span()));
        }
      }

      inline
      void
      recursiveSerialSolveLower(const ordinal_type sid, 
                                const size_type bufsize,
                                /* */ void *buf,
                                const bool final = false) {
        // recursion
        if (final) {
          // do nothing
        } else {
          const ordinal_type 
            ibeg = _info.stree_ptr(sid), 
            iend = _info.stree_ptr(sid+1);
          
          for (ordinal_type i=ibeg;i<iend;++i)
            recursiveSerialSolveLower(_info.stree_children(i), bufsize, buf);
        }

        {
          // dummy member
          const ordinal_type sched = 0, member = 0;

          const ordinal_type nrhs = _info.x.dimension_1();
          ordinal_type m, n; _info.getSuperPanelSize(sid, m, n);          
          UnmanagedViewType<value_type_matrix_host> xB((value_type*)buf, n-m, nrhs);

          CholSupernodes<Algo::Workflow::Serial>
            ::solve_lower(sched, member, _info, xB, sid);
          
          CholSupernodes<Algo::Workflow::Serial>          
            ::update_solve_lower(sched, member, _info, xB, sid);
        }
      }

      inline
      void
      recursiveSerialSolveUpper(const ordinal_type sid, 
                                const size_type bufsize,
                                /* */ void *buf,
                                const bool final = false) {
        {
          // dummy member
          const ordinal_type sched = 0, member = 0;

          const ordinal_type nrhs = _info.x.dimension_1();
          ordinal_type m, n; _info.getSuperPanelSize(sid, m, n);          

          UnmanagedViewType<value_type_matrix_host> xB((value_type*)buf, n-m, nrhs);
          TACHO_TEST_FOR_ABORT(xB.span()*sizeof(value_type) > bufsize,
                               "xB uses greater than bufsize");

          CholSupernodes<Algo::Workflow::Serial>          
            ::update_solve_upper(sched, member, _info, xB, sid);

          CholSupernodes<Algo::Workflow::Serial>
            ::solve_upper(sched, member, _info, xB, sid);          
        }

        // recursion
        if (final) {
          // do nothing
        } else {
          const ordinal_type 
            ibeg = _info.stree_ptr(sid), 
            iend = _info.stree_ptr(sid+1);
          
          for (ordinal_type i=ibeg;i<iend;++i)
            recursiveSerialSolveUpper(_info.stree_children(i), bufsize, buf);
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
                   //const value_type_array_host &ax, // ax is given in factorization
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
          //_ax(ax),
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
          _stree_roots(stree_roots) {
        ///
        /// symbolic input
        /// 
        _info.supernodes             = _supernodes;
        _info.gid_super_panel_ptr    = _gid_super_panel_ptr;
        _info.gid_super_panel_colidx = _gid_super_panel_colidx;
        
        _info.sid_super_panel_ptr    = _sid_super_panel_ptr;
        _info.sid_super_panel_colidx = _sid_super_panel_colidx;
        _info.blk_super_panel_colidx = _blk_super_panel_colidx;
        
        _info.stree_ptr              = _stree_ptr;
        _info.stree_children         = _stree_children;
      }

      inline
      void
      factorizeCholesky_Serial(const value_type_array_host &ax,
                               const ordinal_type verbose = 0) {
        Kokkos::Impl::Timer timer;

        reset_stat();

        timer.reset();
        {
          /// matrix values
          _ax = ax;

          /// allocate super panels 
          ordinal_type_array_host iwork("work", _m+1); 
          _info.allocateSuperPanels(_super_panel_ptr, _super_panel_buf, iwork);

          track_alloc(iwork.span()*sizeof(ordinal_type));
          track_alloc(_super_panel_ptr.span()*sizeof(size_type));
          track_alloc(_super_panel_buf.span()*sizeof(value_type));
        
          /// assign data structure into info
          _info.super_panel_ptr = _super_panel_ptr;
          _info.super_panel_buf = _super_panel_buf;
          
          /// copy the input matrix into super panels
          _info.copySparseToSuperPanels(_ap, _aj, _ax, _perm, _peri, iwork);

          track_free(iwork.span()*sizeof(ordinal_type));
        }
        stat.t_copy += timer.seconds();
                
        timer.reset();
        {
          /// maximum workspace size for serial run
          const size_type worksize = _info.computeWorkspaceSerialChol() + _m + 1;
          value_type_array_host work("work", worksize);

          track_alloc(work.span()*sizeof(value_type));

          /// recursive tree traversal
          const ordinal_type nroots = _stree_roots.dimension_0();
          for (ordinal_type i=0;i<nroots;++i)
            recursiveSerialChol(_stree_roots(i), work.span()*sizeof(value_type), work.data());

          track_free(work.span()*sizeof(value_type));
        }
        stat.t_factor += timer.seconds();
        
        if (verbose) {
          printf("Summary: NumericTools (SerialFactorization)\n");
          printf("===========================================\n");

          print_stat_factor();
        }
      }

      inline
      void
      solveCholesky_Serial(const value_type_matrix_host &x,   // solution
                           const value_type_matrix_host &b,   // right hand side
                           const value_type_matrix_host &t,
                           const ordinal_type verbose = 0) { // temporary workspace (store permuted vectors)
        TACHO_TEST_FOR_EXCEPTION(x.dimension_0() != b.dimension_0() ||
                                 x.dimension_1() != b.dimension_1() || 
                                 x.dimension_0() != t.dimension_0() ||
                                 x.dimension_1() != t.dimension_1(), std::logic_error,
                                 "supernode data structure is not allocated");

        TACHO_TEST_FOR_EXCEPTION(x.data() == b.data() || 
                                 x.data() == t.data() ||
                                 t.data() == b.data(), std::logic_error,
                                 "x, b and t have the same data pointer");

        TACHO_TEST_FOR_EXCEPTION(_info.super_panel_ptr.data() == NULL || 
                                 _info.super_panel_buf.data() == NULL, std::logic_error,
                                 "info's super_panel_ptr/buf is not allocated (factorization is not performed)");

        Kokkos::Impl::Timer timer;

        _info.x = t;

        // copy b -> t
        timer.reset();
        applyRowPermutation(t, b, _peri);
        stat.t_extra += timer.seconds();

        timer.reset();
        {
          /// maximum workspace size for serial run
          const size_type worksize = sqrt(_info.computeWorkspaceSerialChol())*x.dimension_1() + _m;
          value_type_array_host work("work", worksize);

          /// recursive tree traversal
          constexpr bool final = false;
          if (final) {
            const ordinal_type bufsize = work.span()*sizeof(value_type);
            for (ordinal_type sid=0;sid<_nsupernodes;++sid)
              recursiveSerialSolveLower(sid, bufsize, work.data(), true);
            for (ordinal_type sid=(_nsupernodes-1);sid>=0;--sid)
              recursiveSerialSolveUpper(sid, bufsize, work.data(), true);
          } else {
            const ordinal_type nroots = _stree_roots.dimension_0();
            for (ordinal_type i=0;i<nroots;++i)
              recursiveSerialSolveLower(_stree_roots(i), work.span()*sizeof(value_type), work.data());            
            for (ordinal_type i=0;i<nroots;++i)
              recursiveSerialSolveUpper(_stree_roots(i), work.span()*sizeof(value_type), work.data());
          }
        }
        stat.t_solve += timer.seconds();
        
        // copy t -> x
        timer.reset();
        applyRowPermutation(x, t, _perm);
        stat.t_extra += timer.seconds();

        if (verbose) {
          printf("Summary: NumericTools (SerialSolve)\n");
          printf("===================================\n");

          print_stat_solve();
        }
      }

      // inline
      // void
      // factorizeCholesky_Parallel() {
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

      //     /// allocate schur complement workspace (this uses maximum preallocated workspace)
      //     info.allocateWorkspacePerSupernode(_super_schur_ptr, _super_schur_buf, iwork);

      //     info.super_schur_ptr = _super_schur_ptr;
      //     info.super_schur_buf = _super_schur_buf;
      //   }

      //   {
      //     typedef typename sched_type_host::memory_space memory_space;
      //     typedef TaskFunctor_CholSupernodes<value_type,host_exec_space> chol_supernode_functor_type;
          
      //     sched_type_host sched;
      //     {
      //       int nchildren = 0;
      //       const int nsupernodes = _stree_ptr.dimension_0() - 1;
      //       for (ordinal_type sid=0;sid<nsupernodes;++sid)
      //         nchildren = max(nchildren, _stree_ptr(sid+1) - _stree_ptr(sid));

      //       const size_type max_dep_future_size 
      //         = nchildren > MaxDependenceSize ? sizeof(Kokkos::Future<int,host_exec_space>)*nchildren : 16;
            
      //       const size_type max_functor_size = sizeof(chol_supernode_functor_type);
      //       const size_type estimate_max_numtasks = _blk_super_panel_colidx.dimension_0();
            
      //       const ordinal_type
      //         task_queue_capacity = max(estimate_max_numtasks,128)*max_functor_size, 
      //         min_block_size  = min(max_dep_future_size,max_functor_size),
      //         max_block_size  = max(max_dep_future_size,max_functor_size),
      //         num_superblock  = 32,
      //         superblock_size = task_queue_capacity/num_superblock;
            
      //       sched = sched_type_host(memory_space(),
      //                               task_queue_capacity,
      //                               min_block_size,
      //                               max_block_size,
      //                               superblock_size);
      //     }

      //     const ordinal_type nroots = _stree_roots.dimension_0();
      //     for (ordinal_type i=0;i<nroots;++i)
      //       Kokkos::host_spawn(Kokkos::TaskSingle(sched, Kokkos::TaskPriority::High),
      //                          chol_supernode_functor_type(sched, info, _stree_roots(i)));
      //     Kokkos::wait(sched);
      //   }
      // }

      ///
      /// utility
      ///
      static 
      inline
      double 
      computeResidual(const crs_matrix_type_host &A,
                      const value_type_matrix_host &x,
                      const value_type_matrix_host &b) {
        TACHO_TEST_FOR_EXCEPTION(A.NumRows() != A.NumCols() ||
                                 A.NumRows() != b.dimension_0() ||
                                 x.dimension_0() != b.dimension_0() ||
                                 x.dimension_1() != b.dimension_1(), std::logic_error,
                                 "A,x and b dimensions are not compatible");

        const ordinal_type m = A.NumRows(), k = b.dimension_1();        
        double diff = 0, norm = 0;
        for (ordinal_type p=0;p<k;++p) {
          for (ordinal_type i=0;i<m;++i) {
            value_type s = 0;
            const ordinal_type jbeg = A.RowPtrBegin(i), jend = A.RowPtrEnd(i);
            for (ordinal_type j=jbeg;j<jend;++j) {
              const ordinal_type col = A.Col(j);
              s += A.Value(j)*x(col,p);
            }
            norm += b(i,p)*b(i,p);
            diff += (b(i,p) - s)*(b(i,p) - s);
          }
        }
        return sqrt(diff/norm);
      }

      inline
      double 
      computeResidual(const value_type_matrix_host &x,
                      const value_type_matrix_host &b) {
        crs_matrix_type_host A;
        A.setExternalMatrix(_m, _m, _ap(_m),
                            _ap, _aj, _ax);
                            
        return computeResidual(A, x, b);
      }

      inline
      crs_matrix_type_host
      exportFactorsToCrsMatrix(const bool replace_value_with_one = false) {
        /// this only avail after factorization is done
        TACHO_TEST_FOR_EXCEPTION(_info.super_panel_ptr.data() == NULL || 
                                 _info.super_panel_buf.data() == NULL, std::logic_error,
                                 "info's super_panel_ptr/buf is not allocated (factorization is not performed)");

        return _info.createCrsMatrix(replace_value_with_one);
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
