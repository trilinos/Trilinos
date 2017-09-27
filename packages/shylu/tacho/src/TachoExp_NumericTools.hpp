#ifndef __TACHOEXP_NUMERIC_TOOLS_HPP__
#define __TACHOEXP_NUMERIC_TOOLS_HPP__

/// \file TachoExp_NumericTools.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "TachoExp_Util.hpp"
#include "TachoExp_DenseFlopCount.hpp"

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

#include "TachoExp_Trsv.hpp"
#include "TachoExp_Trsv_External.hpp"

#include "TachoExp_Gemv.hpp"
#include "TachoExp_Gemv_External.hpp"

#include "TachoExp_SupernodeInfo.hpp"

#include "TachoExp_CholSupernodes.hpp"
#include "TachoExp_CholSupernodes_Serial.hpp"

#include "TachoExp_TaskFunctor_FactorizeChol.hpp"
#include "TachoExp_TaskFunctor_FactorizeCholByBlocks.hpp"

#include "TachoExp_TaskFunctor_SolveLowerChol.hpp"
#include "TachoExp_TaskFunctor_SolveUpperChol.hpp"

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

      typedef typename supernode_info_type_host::ordinal_type_array ordinal_type_array_host;
      typedef typename supernode_info_type_host::size_type_array size_type_array_host;
      typedef typename supernode_info_type_host::value_type_array value_type_array_host;

      typedef typename supernode_info_type_host::ordinal_pair_type_array ordinal_pair_type_array_host;
      typedef typename supernode_info_type_host::value_type_matrix value_type_matrix_host;
      typedef typename supernode_info_type_host::supernode_type_array supernode_type_array_host;

      typedef typename supernode_info_type_host::dense_block_type dense_block_type_host;
      typedef typename supernode_info_type_host::dense_matrix_of_blocks_type dense_matrix_of_blocks_type_host;

      typedef Kokkos::TaskScheduler<host_exec_space> sched_type_host;
      typedef Kokkos::MemoryPool<host_exec_space> memory_pool_type_host;

      ///
      /// device space typedefs
      ///
      typedef ExecSpace device_exec_space;
      typedef SupernodeInfo<value_type,device_exec_space> supernode_info_type_device;
      typedef typename supernode_info_type_device::crs_matrix_type crs_matrix_type_device;

      typedef typename supernode_info_type_device::ordinal_type_array ordinal_type_array_device;
      typedef typename supernode_info_type_device::size_type_array size_type_array_device;
      typedef typename supernode_info_type_device::value_type_array value_type_array_device;

      typedef typename supernode_info_type_device::ordinal_pair_type_array ordinal_pair_type_array_device;
      typedef typename supernode_info_type_device::value_type_matrix value_type_matrix_device;
      typedef typename supernode_info_type_device::supernode_type_array supernode_type_array_device;

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

      // supernodes       
      ordinal_type _nsupernodes;
      supernode_type_array_host _supernodes;

      // dof mapping to sparse matrix
      ordinal_type_array_host _gid_colidx;

      // supernode map and panel size configuration (sid and column blksize)
      ordinal_pair_type_array_host _sid_block_colidx;

      // supernode tree
      ordinal_type_array_host _stree_roots;

      // output : factors
      value_type_array_host _superpanel_buf;

      ///
      /// supernode info: supernode data structure with "unamanged" view
      /// this is passed into computation algorithm without reference counting
      ///
      supernode_info_type_host _info;

      ///
      /// bufpool # of superblocks
      ///
      ordinal_type _max_num_superblocks;

      /// 
      /// solve phase memoyr pool (reused when it repeat solve)
      ///   - pool capacity returns garbage.
      sched_type_host _sched_solve; size_type _sched_solve_capacity;
      memory_pool_type_host _bufpool_solve; size_type _bufpool_solve_capacity;
      
      ///
      /// statistics
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
        double flop = 0;
        for (ordinal_type sid=0;sid<_nsupernodes;++sid) {
          auto &s = _supernodes(sid);
          const ordinal_type m = s.m, n = s.n - s.m;
          flop += DenseFlopCount<value_type>::Chol(m);
          flop += DenseFlopCount<value_type>::Trsm(true,  m, n);
          flop += DenseFlopCount<value_type>::Syrk(m, n);
        }
        printf("  Time\n");
        printf("             time for copying A into U:                       %10.6f s\n", stat.t_copy);
        printf("             time for numeric factorization:                  %10.6f s\n", stat.t_factor);
        printf("             total time spent:                                %10.6f s\n", (stat.t_copy+stat.t_factor));
        printf("\n");
        printf("  Memory\n");
        printf("             memory used in factorization:                    %10.2f MB\n", stat.m_used/1024/1024);
        printf("             peak memory used in factorization:               %10.2f MB\n", stat.m_peak/1024/1024);
        printf("\n");
        printf("  FLOPs\n");
        printf("             gflop   for numeric factorization:               %10.2f GFLOP\n", flop/1024/1024/1024);
        printf("             gflop/s for numeric factorization:               %10.2f GFLOP/s\n", flop/stat.t_factor/1024/1024/1024);
        printf("\n");
      }
      
      inline
      void
      print_stat_solve() {
        printf("  Time\n");
        printf("             time for extra work e.g.,copy rhs:               %10.6f s\n", stat.t_extra);
        printf("             time for numeric solve:                          %10.6f s\n", stat.t_solve);
        printf("             total time spent:                                %10.6f s\n", (stat.t_solve+stat.t_extra));
        printf("\n");
      }

      inline
      void
      print_stat_memory() {
        printf("  Memory\n"); // better get zero leak
        printf("             leak (or not tracked) memory:                    %10.2f MB\n", stat.m_used/1024/1024);
      }
      
    public:
      NumericTools() 
        : _m(0),
          _sched_solve_capacity(0), 
          _bufpool_solve_capacity(0),
          stat() {}

      NumericTools(const NumericTools &b) = default;
      
      ///
      /// construction (assume input matrix and symbolic are from host)
      ///
      NumericTools(// input matrix A
                   const ordinal_type m,
                   const size_type_array_host &ap,
                   const ordinal_type_array_host &aj,
                   // input permutation
                   const ordinal_type_array_host &perm,
                   const ordinal_type_array_host &peri,
                   // supernodes
                   const ordinal_type nsupernodes,
                   const ordinal_type_array_host &supernodes,
                   const size_type_array_host &gid_ptr,
                   const ordinal_type_array_host &gid_colidx,
                   const size_type_array_host &sid_ptr,
                   const ordinal_type_array_host &sid_colidx,
                   const ordinal_type_array_host &blk_colidx,
                   const ordinal_type_array_host &stree_parent,
                   const size_type_array_host &stree_ptr,
                   const ordinal_type_array_host &stree_children,
                   const ordinal_type_array_host &stree_roots)
      : _m(m), _ap(ap), _aj(aj),
        _perm(perm), _peri(peri),
        _nsupernodes(nsupernodes),
        _gid_colidx(gid_colidx),
        _stree_roots(stree_roots) {        

        reset_stat();

        ///
        /// symbolic input
        ///
        _info.initialize(_supernodes,      
                         _sid_block_colidx,
                         _superpanel_buf,
                         supernodes,
                         gid_ptr,
                         gid_colidx,
                         sid_ptr,
                         sid_colidx,
                         blk_colidx,
                         stree_parent,
                         stree_ptr,
                         stree_children);
        track_alloc(_superpanel_buf.span()*sizeof(value_type));

        ///
        /// max number of superblocks
        ///
        _max_num_superblocks = 4;

        ///
        /// initialize solve scheduler
        ///
        _sched_solve_capacity = 0;
        _bufpool_solve_capacity = 0;
      }

      inline
      void
      release(const ordinal_type verbose = 0) {
        // release bufpool for solve
        if (_bufpool_solve_capacity) {
          track_free(_bufpool_solve_capacity);
          _bufpool_solve = memory_pool_type_host();
        }

        // release scheduler for solve
        if (_sched_solve_capacity) {
          track_free(_sched_solve_capacity);
          _sched_solve = sched_type_host();
        }

        // release supernode buffer
        track_free(_superpanel_buf.span()*sizeof(value_type));
        _superpanel_buf = value_type_array_host();

        if (verbose) {
          printf("Summary: NumericTools (Release)\n");
          printf("===============================\n");
          
          // this should be zero
          print_stat_memory();
        }
      }

      inline
      void
      setSerialThresholdSize(const ordinal_type serial_thres_size) {
        _info.serial_thres_size = serial_thres_size;
      }

      inline
      void
      setMaxNumberOfSuperblocks(const ordinal_type max_num_superblocks) {
        _max_num_superblocks = max_num_superblocks;
      }

      inline
      ordinal_type
      getMaxSupernodeSize() const {
        return _info.max_supernode_size;
      }

      inline
      ordinal_type
      getMaxSchurSize() const {
        return _info.max_schur_size;
      }

      ///
      /// Serial
      /// ------

      inline
      void
      factorizeCholesky_Serial(const value_type_array_host &ax,
                               const ordinal_type verbose = 0) {
        Kokkos::Impl::Timer timer;
        {
          memory_pool_type_host bufpool;
          timer.reset();
          {
            const size_t
              max_block_size = _m*sizeof(ordinal_type);
            
            typedef typename host_exec_space::memory_space memory_space;
            bufpool = memory_pool_type_host(memory_space(),
                                            max_block_size,
                                            max_block_size,
                                            max_block_size,
                                            max_block_size);
            
            track_alloc(bufpool.capacity());
          }
          stat.t_extra = timer.seconds();
          
          timer.reset();
          {
            /// matrix values
            _ax = ax;
            
            /// copy the input matrix into super panels
            _info.copySparseToSuperpanels(_ap, _aj, _ax, _perm, _peri, bufpool);
          }
          stat.t_copy = timer.seconds();
          track_free(bufpool.capacity());
        }

        timer.reset();
        {
          value_type_array_host buf("buf", _info.max_schur_size*(_info.max_schur_size + 1));
          const size_t bufsize = buf.span()*sizeof(value_type);
          track_alloc(bufsize);
          
          /// recursive tree traversal
          const ordinal_type sched = 0, member = 0, nroots = _stree_roots.dimension_0();
          for (ordinal_type i=0;i<nroots;++i)
            CholSupernodes<Algo::Workflow::Serial>
              ::factorize_recursive_serial(sched, member, _info, _stree_roots(i), true, buf.data(), bufsize);
          
          track_free(bufsize);
        }
        stat.t_factor = timer.seconds();
        
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

        Kokkos::Impl::Timer timer;

        _info.x = t;

        // copy b -> t
        timer.reset();
        applyRowPermutation(t, b, _peri);
        stat.t_extra = timer.seconds();
        
        timer.reset();
        {
          value_type_array_host buf("buf", _info.max_schur_size*x.dimension_1());
          const size_t bufsize = buf.span()*sizeof(value_type);
          track_alloc(bufsize);
          
          /// recursive tree traversal
          const ordinal_type sched = 0, member = 0, nroots = _stree_roots.dimension_0();
          for (ordinal_type i=0;i<nroots;++i)
            CholSupernodes<Algo::Workflow::Serial>
              ::solve_lower_recursive_serial(sched, member, _info, _stree_roots(i), true, buf.data(), bufsize);
          for (ordinal_type i=0;i<nroots;++i)
            CholSupernodes<Algo::Workflow::Serial>
              ::solve_upper_recursive_serial(sched, member, _info, _stree_roots(i), true, buf.data(), bufsize);
          
          track_free(bufsize);
        }
        stat.t_solve = timer.seconds();
        
        // copy t -> x
        timer.reset();
        applyRowPermutation(x, t, _perm);
        stat.t_extra += timer.seconds();

        if (verbose) {
          printf("Summary: NumericTools (SerialSolve: %3d)\n", ordinal_type(x.dimension_1()));
          printf("========================================\n");

          print_stat_solve();
        }
      }

      ///
      /// Kokkos Tasking
      /// --------------

      inline
      void
      factorizeCholesky_Parallel(const value_type_array_host &ax,
                                 const ordinal_type verbose = 0) {
        Kokkos::Impl::Timer timer;

        timer.reset();
        typedef typename sched_type_host::memory_space memory_space;
        typedef TaskFunctor_FactorizeChol<value_type,host_exec_space> functor_type;
        //typedef Kokkos::Future<int,host_exec_space> future_type;
        
        sched_type_host sched;
        {
          const size_t max_functor_size = sizeof(functor_type);
          const size_t estimate_max_numtasks = _sid_block_colidx.dimension_0();
          
          const size_t
            task_queue_capacity = max(estimate_max_numtasks,128)*max_functor_size,
            min_block_size  = 1,
            max_block_size  = max_functor_size,
            num_superblock  = 32, // various small size blocks
            superblock_size = task_queue_capacity/num_superblock;
          
          sched = sched_type_host(memory_space(),
                                  task_queue_capacity,
                                  min_block_size,
                                  max_block_size,
                                  superblock_size);
          
          track_alloc(sched.memory()->capacity());
        }
        
        memory_pool_type_host bufpool;
        {
          const size_t
            min_block_size  = 1,
            max_block_size  = max((_info.max_schur_size*_info.max_schur_size +
                                   _info.max_schur_size)*sizeof(value_type),
                                  _m*sizeof(ordinal_type));
          
          size_t superblock_size = 1;
          for ( ;superblock_size<max_block_size;superblock_size*=2);
          
          const size_t
            //num_superblock  = host_exec_space::thread_pool_size(0), // # of threads is safe number
            num_superblock  = min(host_exec_space::thread_pool_size(0), _max_num_superblocks),
            memory_capacity = num_superblock*superblock_size;
          
          bufpool = memory_pool_type_host(memory_space(),
                                          memory_capacity,
                                          min_block_size,
                                          max_block_size,
                                          superblock_size);
          
          track_alloc(bufpool.capacity());
        }
        stat.t_extra = timer.seconds();
        
        timer.reset();
        {
          /// matrix values
          _ax = ax;
          
          /// copy the input matrix into super panels
          _info.copySparseToSuperpanels(_ap, _aj, _ax, _perm, _peri, bufpool);
        }
        stat.t_copy = timer.seconds();

        timer.reset();
        const ordinal_type nroots = _stree_roots.dimension_0();
        for (ordinal_type i=0;i<nroots;++i)
          Kokkos::host_spawn(Kokkos::TaskSingle(sched, Kokkos::TaskPriority::High),
                             functor_type(sched, bufpool, _info, _stree_roots(i)));
        Kokkos::wait(sched);
        stat.t_factor = timer.seconds();
        
        track_free(bufpool.capacity());
        track_free(sched.memory()->capacity());

        // reset solve scheduler and bufpool
        _sched_solve_capacity = 0;
        _bufpool_solve_capacity = 0;

        if (verbose) {
          printf("Summary: NumericTools (ParallelFactorization)\n");
          printf("=============================================\n");

          print_stat_factor();
        }
      }

      inline
      void
      solveCholesky_Parallel(const value_type_matrix_host &x,   // solution
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

        Kokkos::Impl::Timer timer;

        _info.x = t;

        // copy b -> t
        timer.reset();
        applyRowPermutation_Device(t, b, _peri);
        stat.t_extra = timer.seconds();

        {
          timer.reset();
          typedef typename sched_type_host::memory_space memory_space;
          typedef TaskFunctor_SolveLowerChol<value_type,host_exec_space> functor_lower_type;
          typedef TaskFunctor_SolveUpperChol<value_type,host_exec_space> functor_upper_type;
          //typedef Kokkos::Future<int,host_exec_space> future_type;
          
          {
            const size_t max_functor_size = max(sizeof(functor_lower_type), sizeof(functor_upper_type));
            const size_t estimate_max_numtasks = _sid_block_colidx.dimension_0();
            
            const size_t
              task_queue_capacity = max(estimate_max_numtasks,128)*max_functor_size,
              min_block_size  = 1,
              max_block_size  = max_functor_size,
              num_superblock  = 32, // various small size blocks
              superblock_size = task_queue_capacity/num_superblock;

            if (_sched_solve_capacity < task_queue_capacity) {
              _sched_solve = sched_type_host(memory_space(),
                                             task_queue_capacity,
                                             min_block_size,
                                             max_block_size,
                                             superblock_size);
              _sched_solve_capacity = _sched_solve.memory()->capacity();
              track_alloc(_sched_solve_capacity);
            }
          }
          
          {
            const size_t
              min_block_size  = 1,
              max_block_size  = 2*_info.max_schur_size*sizeof(value_type)*x.dimension_1();
            
            size_t superblock_size = 1;
            for ( ;superblock_size<max_block_size;superblock_size*=2);
            
            const size_t
              //num_superblock  = host_exec_space::thread_pool_size(0), // # of threads is safe number
              num_superblock  = min(host_exec_space::thread_pool_size(0), _max_num_superblocks),
              memory_capacity = num_superblock*superblock_size;
            
            if (_bufpool_solve_capacity < memory_capacity) {
              _bufpool_solve = memory_pool_type_host(memory_space(),
                                                     memory_capacity,
                                                     min_block_size,
                                                     max_block_size,
                                                     superblock_size);
              _bufpool_solve_capacity = _bufpool_solve.capacity();
              track_alloc(_bufpool_solve_capacity);
            }
          }
          stat.t_extra += timer.seconds();
          
          timer.reset();
          const ordinal_type nroots = _stree_roots.dimension_0();
          for (ordinal_type i=0;i<nroots;++i) {
            auto fl = Kokkos::host_spawn(Kokkos::TaskSingle(_sched_solve, Kokkos::TaskPriority::High),
                                         functor_lower_type(_sched_solve, _bufpool_solve, _info, _stree_roots(i)));
            auto fu = Kokkos::host_spawn(Kokkos::TaskSingle(fl,           Kokkos::TaskPriority::High),
                                         functor_upper_type(_sched_solve, _bufpool_solve, _info, _stree_roots(i)));
          }
          Kokkos::wait(_sched_solve);
          stat.t_solve = timer.seconds();
        }

        // copy t -> x
        timer.reset();
        applyRowPermutation_Device(x, t, _perm);
        stat.t_extra += timer.seconds();

        if (verbose) {
          printf("Summary: NumericTools (ParallelSolve: %3d)\n", ordinal_type(x.dimension_1()));
          printf("==========================================\n");

          print_stat_solve();
        }
      }

      ///
      /// Kokkos Tasking ByBlocks
      /// -----------------------

      inline
      void
      factorizeCholesky_ParallelByBlocks(const value_type_array_host &ax,
                                         const ordinal_type blksize,
                                         const ordinal_type verbose = 0) {
        Kokkos::Impl::Timer timer;

        timer.reset();
        typedef typename sched_type_host::memory_space memory_space;
        typedef TaskFunctor_FactorizeCholByBlocks<value_type,host_exec_space> functor_type;
        typedef Kokkos::Future<int,host_exec_space> future_type;
        
        const size_t 
          max_nrows_of_blocks = _info.max_supernode_size/blksize + 1,
          max_ncols_of_blocks = _info.max_schur_size/blksize + 1;
        
        sched_type_host sched;
        {
          const size_t max_dep_future_size = max_ncols_of_blocks*max_ncols_of_blocks*sizeof(future_type);
          const size_t max_functor_size = sizeof(functor_type);
          const size_t estimate_max_numtasks = _sid_block_colidx.dimension_0();
          
          const size_t
            task_queue_capacity = max(estimate_max_numtasks,128)*max_functor_size,
            min_block_size  = 1,
            max_block_size  = ( max_dep_future_size + max_functor_size ),
            num_superblock  = 32, // various small size blocks
            superblock_size = task_queue_capacity/num_superblock;
          
          sched = sched_type_host(memory_space(),
                                  task_queue_capacity,
                                  min_block_size,
                                  max_block_size,
                                  superblock_size);
          
          track_alloc(sched.memory()->capacity());
        }
        
        memory_pool_type_host bufpool;
        {
          const size_t
            min_block_size  = 1,
            max_block_size  = max(( (_info.max_schur_size*_info.max_schur_size +
                                     _info.max_schur_size)*sizeof(value_type) +
                                    (max_nrows_of_blocks*max_nrows_of_blocks +
                                     max_nrows_of_blocks*max_ncols_of_blocks +
                                     max_ncols_of_blocks*max_ncols_of_blocks)*sizeof(dense_block_type_host) ),
                                  _m*sizeof(ordinal_type));
                                  
          size_t superblock_size = 1;
          for ( ;superblock_size<max_block_size;superblock_size*=2);
          
          const size_t
            //num_superblock  = host_exec_space::thread_pool_size(0), // # of threads is safe number
            num_superblock  = min(host_exec_space::thread_pool_size(0), _max_num_superblocks),
            memory_capacity = num_superblock*superblock_size;
          
          bufpool = memory_pool_type_host(memory_space(),
                                          memory_capacity,
                                          min_block_size,
                                          max_block_size,
                                          superblock_size);
          
          track_alloc(bufpool.capacity());
        }
        stat.t_extra += timer.seconds();
        
        timer.reset();
        {
          /// matrix values
          _ax = ax;
          
          /// copy the input matrix into super panels
          _info.copySparseToSuperpanels(_ap, _aj, _ax, _perm, _peri, bufpool);
        }
        stat.t_copy += timer.seconds();

        timer.reset();
        const ordinal_type nroots = _stree_roots.dimension_0();
        for (ordinal_type i=0;i<nroots;++i)
          Kokkos::host_spawn(Kokkos::TaskSingle(sched, Kokkos::TaskPriority::High),
                             functor_type(sched, bufpool, _info, _stree_roots(i), blksize));
        Kokkos::wait(sched);
        stat.t_factor = timer.seconds();
        
        track_free(bufpool.capacity());
        track_free(sched.memory()->capacity());

        // reset solve scheduler and bufpool
        _sched_solve_capacity = 0;
        _bufpool_solve_capacity = 0;

        if (verbose) {
          printf("Summary: NumericTools (ParallelFactorizationByBlocks: %3d)\n", blksize);
          printf("==========================================================\n");

          print_stat_factor();
        }
      }


      ///
      /// utility
      ///
      static
      inline
      double
      computeRelativeResidual(const crs_matrix_type_host &A,
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
            norm += real(b(i,p)*conj(b(i,p)));
            diff += real((b(i,p) - s)*conj(b(i,p) - s));
          }
        }
        return sqrt(diff/norm);
      }

      inline
      double
      computeRelativeResidual(const value_type_matrix_host &x,
                              const value_type_matrix_host &b) {
        crs_matrix_type_host A;
        A.setExternalMatrix(_m, _m, _ap(_m),
                            _ap, _aj, _ax);

        return computeRelativeResidual(A, x, b);
      }

      inline
      crs_matrix_type_host
      exportFactorsToCrsMatrix(const bool replace_value_with_one = false) {
        /// this only avail after factorization is done
        // TACHO_TEST_FOR_EXCEPTION(_info.super_panel_ptr.data() == NULL ||
        //                          _info.super_panel_buf.data() == NULL, std::logic_error,
        //                          "info's super_panel_ptr/buf is not allocated (factorization is not performed)");

        return _info.createCrsMatrix(replace_value_with_one);
      }

    };

  }
}
#endif
