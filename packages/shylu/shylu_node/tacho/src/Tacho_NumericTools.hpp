#ifndef __TACHO_NUMERIC_TOOLS_HPP__
#define __TACHO_NUMERIC_TOOLS_HPP__

/// \file Tacho_NumericTools.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"
#include "Tacho_DenseFlopCount.hpp"

#include "Tacho_CrsMatrixBase.hpp"
#include "Tacho_DenseMatrixView.hpp"

#include "Tacho_Chol.hpp"
#include "Tacho_Chol_External.hpp"

#include "Tacho_Trsm.hpp"
#include "Tacho_Trsm_External.hpp"

#include "Tacho_Herk.hpp"
#include "Tacho_Herk_External.hpp"

#include "Tacho_Gemm.hpp"
#include "Tacho_Gemm_External.hpp"

#include "Tacho_Trsv.hpp"
#include "Tacho_Trsv_External.hpp"

#include "Tacho_Gemv.hpp"
#include "Tacho_Gemv_External.hpp"

#include "Tacho_SupernodeInfo.hpp"

#include "Tacho_CholSupernodes.hpp"
#include "Tacho_CholSupernodes_Serial.hpp"
#include "Tacho_CholSupernodes_SerialPanel.hpp"

#include "Tacho_TaskFunctor_FactorizeChol.hpp"
#include "Tacho_TaskFunctor_FactorizeCholPanel.hpp"
#include "Tacho_TaskFunctor_FactorizeCholByBlocks.hpp"
#include "Tacho_TaskFunctor_FactorizeCholByBlocksPanel.hpp"

#include "Tacho_TaskFunctor_SolveLowerChol.hpp"
#include "Tacho_TaskFunctor_SolveUpperChol.hpp"

namespace Tacho {

    template<typename ValueType, typename ExecSpace>
    class NumericTools {
    public:
      typedef ValueType value_type;
      typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;

      ///
      /// device space typedefs
      ///
      typedef ExecSpace exec_space;
      typedef SupernodeInfo<value_type,exec_space> supernode_info_type;
      typedef typename supernode_info_type::crs_matrix_type crs_matrix_type;

      typedef typename supernode_info_type::ordinal_type_array ordinal_type_array;
      typedef typename supernode_info_type::size_type_array size_type_array;
      typedef typename supernode_info_type::value_type_array value_type_array;

      typedef typename supernode_info_type::ordinal_pair_type_array ordinal_pair_type_array;
      typedef typename supernode_info_type::value_type_matrix value_type_matrix;
      typedef typename supernode_info_type::supernode_type_array supernode_type_array;

      typedef typename supernode_info_type::dense_block_type dense_block_type;
      typedef typename supernode_info_type::dense_matrix_of_blocks_type dense_matrix_of_blocks_type;

      typedef Kokkos::TaskScheduler<exec_space> scheduler_type;
      typedef Kokkos::MemoryPool<exec_space> memory_pool_type;

      typedef Kokkos::DefaultHostExecutionSpace host_space;
      typedef Kokkos::View<ordinal_type*,host_space> ordinal_type_array_host;

    private:

      ///
      /// supernode data structure memory "managed"
      /// this holds all necessary connectivity data
      ///

      // matrix input
      ordinal_type _m;
      size_type_array _ap;
      ordinal_type_array _aj;
      value_type_array _ax;

      // graph ordering input
      ordinal_type_array _perm, _peri;
      
      // supernodes       
      ordinal_type _nsupernodes;
      supernode_type_array _supernodes;

      // dof mapping to sparse matrix
      ordinal_type_array _gid_colidx;

      // supernode map and panel size configuration (sid and column blksize)
      ordinal_pair_type_array _sid_block_colidx;

      // supernode tree
      ordinal_type_array_host _stree_roots;

      // output : factors
      value_type_array _superpanel_buf;

      ///
      /// supernode info: supernode data structure with "unamanged" view
      /// this is passed into computation algorithm without reference counting
      ///
      supernode_info_type _info;

      ///
      /// bufpool # of superblocks
      ///
      ordinal_type _max_num_superblocks;

      /// 
      /// solve phase memory pool (reused when it repeat solve)
      ///   - pool capacity returns garbage.
      scheduler_type _sched_solve; size_t _sched_solve_capacity;
      memory_pool_type _bufpool_solve; size_t _bufpool_solve_capacity;
      
      ///
      /// statistics
      ///
      struct {
        double t_factor, t_solve, t_copy, t_extra;
        double m_used, m_peak;
        size_t b_min_block_size, b_max_block_size, b_capacity, b_num_superblocks;
        size_t s_min_block_size, s_max_block_size, s_capacity, s_num_superblocks;
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

        stat.b_min_block_size = 0;
        stat.b_max_block_size = 0;
        stat.b_capacity = 0;
        stat.b_num_superblocks = 0;

        stat.s_min_block_size = 0;
        stat.s_max_block_size = 0;
        stat.s_capacity = 0;
        stat.s_num_superblocks = 0;
      }

      inline
      void
      print_stat_factor() {
        double flop = 0;
        auto h_supernodes = Kokkos::create_mirror_view(_supernodes);                    
        Kokkos::deep_copy(h_supernodes, _supernodes);  
        
        for (ordinal_type sid=0;sid<_nsupernodes;++sid) {
          auto &s = h_supernodes(sid);
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
        printf("  Buffer Pool\n");
        printf("             number of superblocks:                           %10d\n", int(stat.b_num_superblocks));
        printf("             min and max blocksize:                           %e, %e Bytes\n", double(stat.b_min_block_size), double(stat.b_max_block_size));
        printf("             pool capacity:                                   %10.2f MB\n", double(stat.b_capacity)/1024/1024);
        printf("\n");
        printf("  Sched Memory Pool\n");
        printf("             number of superblocks:                           %10d\n", int(stat.s_num_superblocks));
        printf("             min and max blocksize:                           %e, %e Bytes\n", double(stat.s_min_block_size), double(stat.s_max_block_size));
        printf("             pool capacity:                                   %10.2f MB\n", double(stat.s_capacity)/1024/1024);
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
                   const    size_type_array &ap,
                   const ordinal_type_array &aj,
                   // input permutation
                   const ordinal_type_array &perm,
                   const ordinal_type_array &peri,
                   // supernodes
                   const ordinal_type nsupernodes,
                   const ordinal_type_array &supernodes,
                   const    size_type_array &gid_ptr,
                   const ordinal_type_array &gid_colidx,
                   const    size_type_array &sid_ptr,
                   const ordinal_type_array &sid_colidx,
                   const ordinal_type_array &blk_colidx,
                   const ordinal_type_array &stree_parent,
                   const    size_type_array &stree_ptr,
                   const ordinal_type_array &stree_children,
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
          _bufpool_solve = memory_pool_type();
        }

        // release scheduler for solve
        if (_sched_solve_capacity) {
          track_free(_sched_solve_capacity);
          _sched_solve = scheduler_type();
        }

        // release supernode buffer
        track_free(_superpanel_buf.span()*sizeof(value_type));
        _superpanel_buf = value_type_array();

        if (verbose) {
          printf("Summary: NumericTools (Release)\n");
          printf("===============================\n");
          
          // this should be zero
          print_stat_memory();
        }
      }

      inline
      void
      setFrontUpdateMode(const ordinal_type front_update_mode) {
        _info.front_update_mode = front_update_mode;
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

      ///
      /// This cannot be called for both host and gpus as nvcc compiles twice and
      /// it should match the interface
      /// 

#if !defined (KOKKOS_ENABLE_CUDA)
      inline
      void
      factorizeCholesky_Serial(const value_type_array &ax,
                               const ordinal_type verbose = 0) {
        if (!std::is_same<typename exec_space::memory_space,Kokkos::HostSpace>::value) {
          printf("factorizeCholesky_Serial cannot be executed on non-host devices\n");          
          return;
        }

        Kokkos::Impl::Timer timer;
        {
          timer.reset();
          {
            /// matrix values
            _ax = ax;
            
            /// copy the input matrix into super panels
            _info.copySparseToSuperpanels(_ap, _aj, _ax, _perm, _peri);
          }
          stat.t_copy = timer.seconds();
        }
          
        timer.reset();
        {
          value_type_array buf("buf", _info.max_schur_size*(_info.max_schur_size + 1));
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
      factorizeCholesky_SerialPanel(const value_type_array &ax,
                                    const ordinal_type panelsize,
                                    const ordinal_type verbose = 0) {
        if (!std::is_same<typename exec_space::memory_space,Kokkos::HostSpace>::value) {
          printf("factorizeCholesky_SerialPanel cannot be executed on non-host devices\n");          
          return;
        }

        Kokkos::Impl::Timer timer;
        {
          timer.reset();
          {
            /// matrix values
            _ax = ax;
            
            /// copy the input matrix into super panels
            _info.copySparseToSuperpanels(_ap, _aj, _ax, _perm, _peri);
          }
          stat.t_copy = timer.seconds();
        }

        const ordinal_type nb = panelsize > 0 ? panelsize : _info.max_schur_size;
        timer.reset();
        {
          value_type_array buf("buf", _info.max_schur_size*(nb + 1));
          const size_t bufsize = buf.span()*sizeof(value_type);
          track_alloc(bufsize);
          
          /// recursive tree traversal
          const ordinal_type sched = 0, member = 0, nroots = _stree_roots.dimension_0();
          for (ordinal_type i=0;i<nroots;++i)
            CholSupernodes<Algo::Workflow::SerialPanel>
              ::factorize_recursive_serial(sched, member, 
                                           _info, _stree_roots(i), 
                                           true, buf.data(), bufsize, nb);
          
          track_free(bufsize);
        }
        stat.t_factor = timer.seconds();
        
        if (verbose) {
          printf("Summary: NumericTools (SerialPanelFactorization: %3d)\n", nb);
          printf("=====================================================\n");

          print_stat_factor();
        }
      }

      inline
      void
      solveCholesky_Serial(const value_type_matrix &x,   // solution
                           const value_type_matrix &b,   // right hand side
                           const value_type_matrix &t,   // temporary workspace (store permuted vectors)
                           const ordinal_type verbose = 0) {
        if (!std::is_same<typename exec_space::memory_space,Kokkos::HostSpace>::value) {
          printf("solveCholesky_Serial cannot be executed on non-host devices\n");          
          return;
        }

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
          value_type_array buf("buf", _info.max_schur_size*x.dimension_1());
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
#endif

      ///
      /// Kokkos Tasking
      /// --------------

      inline
      void
      factorizeCholesky_Parallel(const value_type_array &ax,
                                 const ordinal_type verbose = 0) {
        Kokkos::Impl::Timer timer;

        timer.reset();

        typedef TaskFunctor_FactorizeChol<value_type,exec_space> functor_type;
        typedef Kokkos::Future<int,exec_space> future_type;

        scheduler_type sched;
        {
          const size_t max_functor_size = 2*sizeof(functor_type);
          const size_t max_dep_future_size = _info.max_nchildren*sizeof(future_type);          
          const size_t estimate_max_numtasks = _sid_block_colidx.dimension_0() >> 3;
          
          const size_t
            task_queue_capacity = max(estimate_max_numtasks,128)*max_functor_size,
            min_block_size  = 16,
            max_block_size  = (max_functor_size + max_dep_future_size),
            num_superblock  = 32, // various small size blocks
            superblock_size = task_queue_capacity/num_superblock;
          
          sched = scheduler_type(typename scheduler_type::memory_space(),
                             task_queue_capacity,
                             min_block_size,
                             max_block_size,
                             superblock_size);
          
          track_alloc(sched.memory()->capacity());
        }
        
        stat.s_min_block_size  = sched.memory()->min_block_size();
        stat.s_max_block_size  = sched.memory()->max_block_size();
        stat.s_capacity        = sched.memory()->capacity();
        stat.s_num_superblocks = sched.memory()->capacity()/sched.memory()->max_block_size();
        
        memory_pool_type bufpool;
        {
          const size_t cuda_max_team_size 
            = std::is_same<typename exec_space::memory_space,Kokkos::HostSpace>::value 
            ? 1 : 32;
          const size_t
            min_block_size  = 16,
            max_block_size  = max((_info.max_schur_size*_info.max_schur_size +
                                   _info.max_schur_size*cuda_max_team_size)*sizeof(value_type),
                                  _m*sizeof(ordinal_type));
          
          ordinal_type ishift = 0;
          size_t superblock_size = 1;
          for ( ;superblock_size<max_block_size;superblock_size <<= 1,++ishift);

          size_t num_superblock_device = 0;
          if (std::is_same<typename exec_space::memory_space,Kokkos::HostSpace>::value) {
            // host space 
            num_superblock_device = host_space::thread_pool_size(0);
          } else {
            // cuda space (what would be best number)
            num_superblock_device = _max_num_superblocks;
          }

          const size_t  
            max_num_superblocks = _max_num_superblocks, //min(1ul << (ishift > 31 ? 0 : 31 - ishift), _max_num_superblocks),
            num_superblock  = min(num_superblock_device, max_num_superblocks),
            memory_capacity = num_superblock*superblock_size;
          
          bufpool = memory_pool_type(typename exec_space::memory_space(),
                                     memory_capacity,
                                     min_block_size,
                                     max_block_size,
                                     superblock_size);
          
          track_alloc(bufpool.capacity());
        }
        stat.t_extra = timer.seconds();

        stat.b_min_block_size  = bufpool.min_block_size();
        stat.b_max_block_size  = bufpool.max_block_size();
        stat.b_capacity        = bufpool.capacity();
        stat.b_num_superblocks = bufpool.capacity()/bufpool.max_block_size();

        timer.reset();
        {
          /// matrix values
          _ax = ax;
          
          /// copy the input matrix into super panels
          _info.copySparseToSuperpanels(_ap, _aj, _ax, _perm, _peri);
        }
        stat.t_copy = timer.seconds();

        timer.reset();
        const ordinal_type nroots = _stree_roots.dimension_0();
        for (ordinal_type i=0;i<nroots;++i)
          Kokkos::host_spawn(Kokkos::TaskTeam(sched, Kokkos::TaskPriority::High),
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
      factorizeCholesky_ParallelPanel(const value_type_array &ax,
                                      const ordinal_type panelsize, 
                                      const ordinal_type verbose = 0) {
        Kokkos::Impl::Timer timer;

        timer.reset();

        typedef TaskFunctor_FactorizeCholPanel<value_type,exec_space> functor_type;
        typedef Kokkos::Future<int,exec_space> future_type;
        
        scheduler_type sched;
        {
          const size_t max_functor_size = 2*sizeof(functor_type);
          const size_t max_dep_future_size = _info.max_nchildren*sizeof(future_type);          
          const size_t estimate_max_numtasks = _sid_block_colidx.dimension_0() >> 3;
          
          const size_t
            task_queue_capacity = max(estimate_max_numtasks,128)*max_functor_size,
            min_block_size  = 16,
            max_block_size  = (max_functor_size + max_dep_future_size),
            num_superblock  = 32, // various small size blocks
            superblock_size = task_queue_capacity/num_superblock;
          
          sched = scheduler_type(typename scheduler_type::memory_space(),
                             task_queue_capacity,
                             min_block_size,
                             max_block_size,
                             superblock_size);
          
          track_alloc(sched.memory()->capacity());
        }

        stat.s_min_block_size  = sched.memory()->min_block_size();
        stat.s_max_block_size  = sched.memory()->max_block_size();
        stat.s_capacity        = sched.memory()->capacity();
        stat.s_num_superblocks = sched.memory()->capacity()/sched.memory()->max_block_size();
        
        memory_pool_type bufpool;
        const ordinal_type nb = panelsize > 0 ? panelsize : _info.max_schur_size;
        {
          const size_t cuda_max_team_size 
            = std::is_same<typename exec_space::memory_space,Kokkos::HostSpace>::value 
            ? 1 : 32;
          const size_t
            min_block_size       = 16,
            max_block_size_left  = (nb*_info.max_schur_size + 
                                    _info.max_schur_size*cuda_max_team_size)*sizeof(value_type),
            max_block_size_right = _m*sizeof(ordinal_type),
            max_block_size       = max(max_block_size_left, max_block_size_right);
          
          ordinal_type ishift = 0;
          size_t superblock_size = 1;
          for ( ;superblock_size<max_block_size;superblock_size <<= 1,++ishift);

          size_t num_superblock_device = 0;
          if (std::is_same<typename exec_space::memory_space,Kokkos::HostSpace>::value) {
            // host space 
            num_superblock_device = host_space::thread_pool_size(0);
          } else {
            // cuda space (what would be best number)
            num_superblock_device = _max_num_superblocks;
          }

          const size_t  
            max_num_superblocks = _max_num_superblocks, 
            num_superblock  = min(num_superblock_device, max_num_superblocks),
            memory_capacity = num_superblock*superblock_size;
          
          bufpool = memory_pool_type(typename exec_space::memory_space(),
                                     memory_capacity,
                                     min_block_size,
                                     max_block_size,
                                     superblock_size);
          
          track_alloc(bufpool.capacity());
        }
        stat.t_extra = timer.seconds();

        stat.b_min_block_size  = bufpool.min_block_size();
        stat.b_max_block_size  = bufpool.max_block_size();
        stat.b_capacity        = bufpool.capacity();
        stat.b_num_superblocks = bufpool.capacity()/bufpool.max_block_size();

        timer.reset();
        {
          /// matrix values
          _ax = ax;
          
          /// copy the input matrix into super panels
          _info.copySparseToSuperpanels(_ap, _aj, _ax, _perm, _peri);
        }
        stat.t_copy = timer.seconds();
        
        timer.reset();
        const ordinal_type nroots = _stree_roots.dimension_0();
        for (ordinal_type i=0;i<nroots;++i)
          Kokkos::host_spawn(Kokkos::TaskTeam(sched, Kokkos::TaskPriority::High),
                             functor_type(sched, bufpool, _info, _stree_roots(i), nb));
        Kokkos::wait(sched);
        stat.t_factor = timer.seconds();
        
        track_free(bufpool.capacity());
        track_free(sched.memory()->capacity());

        // reset solve scheduler and bufpool
        _sched_solve_capacity = 0;
        _bufpool_solve_capacity = 0;

        if (verbose) {
          printf("Summary: NumericTools (ParallelPanelFactorization: %3d)\n", nb);
          printf("=======================================================\n");

          print_stat_factor();
        }
      }

      inline
      void
      solveCholesky_Parallel(const value_type_matrix &x,   // solution
                             const value_type_matrix &b,   // right hand side
                             const value_type_matrix &t,
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

        {
          timer.reset();

          typedef TaskFunctor_SolveLowerChol<value_type,exec_space> functor_lower_type;
          typedef TaskFunctor_SolveUpperChol<value_type,exec_space> functor_upper_type;
          typedef Kokkos::Future<int,exec_space> future_type;
          
          {
            const size_t max_functor_size = 2*max(sizeof(functor_lower_type), sizeof(functor_upper_type));
            const size_t max_dep_future_size = _info.max_nchildren*sizeof(future_type);
            const size_t estimate_max_numtasks = _sid_block_colidx.dimension_0()/2;
            
            const size_t
              task_queue_capacity = max(estimate_max_numtasks,128)*max_functor_size,
              min_block_size  = 16,
              max_block_size  = (max_functor_size + max_dep_future_size),
              num_superblock  = 32, // various small size blocks
              superblock_size = task_queue_capacity/num_superblock;

            if (_sched_solve_capacity < task_queue_capacity) {
              _sched_solve = scheduler_type(typename scheduler_type::memory_space(),
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
              min_block_size  = 16,
              max_block_size  = 2*_info.max_schur_size*sizeof(value_type)*x.dimension_1();
            
            size_t superblock_size = 1;
            for ( ;superblock_size<max_block_size;superblock_size*=2);

            size_t num_superblock_device = 0;
            if (std::is_same<typename exec_space::memory_space,Kokkos::HostSpace>::value) {
              // host space 
              num_superblock_device = host_space::thread_pool_size(0);
            } else {
              // cuda space (what would be best number)
              num_superblock_device = _max_num_superblocks;
            }
            
            const size_t
              num_superblock  = min(num_superblock_device, _max_num_superblocks),
              memory_capacity = num_superblock*superblock_size;
            
            if (_bufpool_solve_capacity < memory_capacity) {
              _bufpool_solve = memory_pool_type(typename exec_space::memory_space(),
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
            auto fl = Kokkos::host_spawn(Kokkos::TaskTeam(_sched_solve, Kokkos::TaskPriority::High),
                                         functor_lower_type(_sched_solve, _bufpool_solve, _info, _stree_roots(i)));
            auto fu = Kokkos::host_spawn(Kokkos::TaskTeam(          fl, Kokkos::TaskPriority::High),
                                         functor_upper_type(_sched_solve, _bufpool_solve, _info, _stree_roots(i)));
          }
          Kokkos::wait(_sched_solve);
          stat.t_solve = timer.seconds();
        }

        // copy t -> x
        timer.reset();
        applyRowPermutation(x, t, _perm);
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
      factorizeCholesky_ParallelByBlocks(const value_type_array &ax,
                                         const ordinal_type blksize,
                                         const ordinal_type verbose = 0) {
        Kokkos::Impl::Timer timer;

        timer.reset();

        typedef TaskFunctor_FactorizeCholByBlocks<value_type,exec_space> functor_type;
        typedef Kokkos::Future<int,exec_space> future_type;
        
        const size_t 
          max_nrows_of_blocks = _info.max_supernode_size/blksize + 1,
          max_ncols_of_blocks = _info.max_schur_size/blksize + 1;
        
        scheduler_type sched;
        {
          const size_t max_dep_future_size 
            = (_info.max_nchildren + max_ncols_of_blocks*max_ncols_of_blocks)*sizeof(future_type);
          const size_t max_functor_size = 2*sizeof(functor_type);
          const size_t estimate_max_numtasks = _sid_block_colidx.dimension_0() >> 3;
          
          const size_t
            task_queue_capacity = max(estimate_max_numtasks,128)*max_functor_size,
            min_block_size  = 16,
            max_block_size  = ( max_dep_future_size + max_functor_size ),
            num_superblock  = 32, // various small size blocks
            superblock_size = task_queue_capacity/num_superblock;
          
          sched = scheduler_type(typename scheduler_type::memory_space(),
                             task_queue_capacity,
                             min_block_size,
                             max_block_size,
                             superblock_size);
          
          track_alloc(sched.memory()->capacity());
        }
        
        stat.s_min_block_size  = sched.memory()->min_block_size();
        stat.s_max_block_size  = sched.memory()->max_block_size();
        stat.s_capacity        = sched.memory()->capacity();
        stat.s_num_superblocks = sched.memory()->capacity()/sched.memory()->max_block_size();
        
        memory_pool_type bufpool;
        {
          const size_t cuda_max_team_size 
            = std::is_same<typename exec_space::memory_space,Kokkos::HostSpace>::value 
            ? 1 : 32;
          const size_t
            min_block_size  = 16,
            max_block_size  = max(( (_info.max_schur_size*_info.max_schur_size +
                                     _info.max_schur_size*cuda_max_team_size)*sizeof(value_type) +
                                    (max_nrows_of_blocks*max_nrows_of_blocks +
                                     max_nrows_of_blocks*max_ncols_of_blocks +
                                     max_ncols_of_blocks*max_ncols_of_blocks)*sizeof(dense_block_type) ),
                                  _m*sizeof(ordinal_type));
                                  
          ordinal_type ishift = 0;
          size_t superblock_size = 1;
          for ( ;superblock_size<max_block_size;superblock_size <<= 1,++ishift);

          size_t num_superblock_device = 0;
          if (std::is_same<typename exec_space::memory_space,Kokkos::HostSpace>::value) {
            // host space 
            num_superblock_device = host_space::thread_pool_size(0);
          } else {
            // cuda space (what would be best number)
            num_superblock_device = _max_num_superblocks;
          }
          
          const size_t 
            max_num_superblocks = _max_num_superblocks, //min(1ul << (ishift > 31 ? 0 : 31 - ishift), _max_num_superblocks),
            num_superblock  = min(num_superblock_device, max_num_superblocks),
            memory_capacity = num_superblock*superblock_size;

          bufpool = memory_pool_type(typename exec_space::memory_space(),
                                     memory_capacity,
                                     min_block_size,
                                     max_block_size,
                                     superblock_size);
          
          track_alloc(bufpool.capacity());
        }
        stat.t_extra += timer.seconds();
        
        stat.b_min_block_size = bufpool.min_block_size();
        stat.b_max_block_size = bufpool.max_block_size();
        stat.b_capacity = bufpool.capacity();
        stat.b_num_superblocks = bufpool.capacity()/bufpool.max_block_size();

        timer.reset();
        {
          /// matrix values
          _ax = ax;
          
          /// copy the input matrix into super panels
          _info.copySparseToSuperpanels(_ap, _aj, _ax, _perm, _peri);
        }
        stat.t_copy += timer.seconds();

        timer.reset();
        const ordinal_type nroots = _stree_roots.dimension_0();
        for (ordinal_type i=0;i<nroots;++i)
          Kokkos::host_spawn(Kokkos::TaskTeam(sched, Kokkos::TaskPriority::High),
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

      inline
      void
      factorizeCholesky_ParallelByBlocksPanel(const value_type_array &ax,
                                              const ordinal_type blksize,
                                              const ordinal_type panelsize,
                                              const ordinal_type verbose = 0) {
        Kokkos::Impl::Timer timer;

        timer.reset();

        typedef TaskFunctor_FactorizeCholByBlocksPanel<value_type,exec_space> functor_type;
        typedef Kokkos::Future<int,exec_space> future_type;
        
        const size_t 
          max_nrows_of_blocks = _info.max_supernode_size/blksize + 1,
          max_ncols_of_blocks = _info.max_schur_size/blksize + 1;
        
        scheduler_type sched;
        {
          const size_t max_dep_future_size 
            = (_info.max_nchildren + max_nrows_of_blocks*max_ncols_of_blocks)*sizeof(future_type);
          const size_t max_functor_size = 2*sizeof(functor_type);
          const size_t estimate_max_numtasks = _sid_block_colidx.dimension_0() >> 3;
          
          const size_t
            task_queue_capacity = max(estimate_max_numtasks,128)*max_functor_size,
            min_block_size  = 16,
            max_block_size  = ( max_dep_future_size + max_functor_size ),
            num_superblock  = 32, // various small size blocks
            superblock_size = task_queue_capacity/num_superblock;
          
          sched = scheduler_type(typename scheduler_type::memory_space(),
                             task_queue_capacity,
                             min_block_size,
                             max_block_size,
                             superblock_size);
          
          track_alloc(sched.memory()->capacity());
        }
        
        stat.s_min_block_size  = sched.memory()->min_block_size();
        stat.s_max_block_size  = sched.memory()->max_block_size();
        stat.s_capacity        = sched.memory()->capacity();
        stat.s_num_superblocks = sched.memory()->capacity()/sched.memory()->max_block_size();
        
        memory_pool_type bufpool;
        {
          const size_t cuda_max_team_size 
            = std::is_same<typename exec_space::memory_space,Kokkos::HostSpace>::value 
            ? 1 : 32;
          const size_t
            min_block_size  = 16,
            max_block_size  = max(( (_info.max_schur_size*min(panelsize,_info.max_schur_size) +
                                     _info.max_schur_size*cuda_max_team_size)*sizeof(value_type) +
                                    (max_nrows_of_blocks*max_nrows_of_blocks +
                                     max_nrows_of_blocks*max_ncols_of_blocks + 
                                     max_nrows_of_blocks)*sizeof(dense_block_type) ),
                                  _m*sizeof(ordinal_type));
          
          ordinal_type ishift = 0;
          size_t superblock_size = 1;
          for ( ;superblock_size<max_block_size;superblock_size <<= 1,++ishift);

          size_t num_superblock_device = 0;
          if (std::is_same<typename exec_space::memory_space,Kokkos::HostSpace>::value) {
            // host space 
            num_superblock_device = host_space::thread_pool_size(0);
          } else {
            // cuda space (what would be best number)
            num_superblock_device = _max_num_superblocks;
          }

          const size_t // max 2 GB allows
            max_num_superblocks = _max_num_superblocks, 
            num_superblock  = min(num_superblock_device, max_num_superblocks),
            memory_capacity = num_superblock*superblock_size;

          bufpool = memory_pool_type(typename exec_space::memory_space(),
                                     memory_capacity,
                                     min_block_size,
                                     max_block_size,
                                     superblock_size);
          
          track_alloc(bufpool.capacity());
        }
        stat.t_extra += timer.seconds();
        
        stat.b_min_block_size = bufpool.min_block_size();
        stat.b_max_block_size = bufpool.max_block_size();
        stat.b_capacity = bufpool.capacity();
        stat.b_num_superblocks = bufpool.capacity()/bufpool.max_block_size();

        timer.reset();
        {
          /// matrix values
          _ax = ax;
          
          /// copy the input matrix into super panels
          _info.copySparseToSuperpanels(_ap, _aj, _ax, _perm, _peri);
        }
        stat.t_copy += timer.seconds();
        
        timer.reset();
        const ordinal_type nroots = _stree_roots.dimension_0();
        for (ordinal_type i=0;i<nroots;++i)
          Kokkos::host_spawn(Kokkos::TaskTeam(sched, Kokkos::TaskPriority::High),
                             functor_type(sched, bufpool, _info, _stree_roots(i), blksize, panelsize));
        Kokkos::wait(sched);
        stat.t_factor = timer.seconds();
        
        track_free(bufpool.capacity());
        track_free(sched.memory()->capacity());

        // reset solve scheduler and bufpool
        _sched_solve_capacity = 0;
        _bufpool_solve_capacity = 0;

        if (verbose) {
          printf("Summary: NumericTools (ParallelFactorizationByBlocksPanel: %3d, %3d)\n", blksize, panelsize);
          printf("====================================================================\n");

          print_stat_factor();
        }
      }


      ///
      /// Utility on device
      ///
      inline
      void
      exportFactorsToCrsMatrix(crs_matrix_type &A,
                               const bool replace_value_with_one = false) {
        /// this only avail after factorization is done
        // TACHO_TEST_FOR_EXCEPTION(_info.super_panel_ptr.data() == NULL ||
        //                          _info.super_panel_buf.data() == NULL, std::logic_error,
        //                          "info's super_panel_ptr/buf is not allocated (factorization is not performed)");
        _info.createCrsMatrix(A, replace_value_with_one);
      }

      ///
      /// Utility on host with deep copy
      ///
      static
      inline
      double
      computeRelativeResidual(const crs_matrix_type &A,
                              const value_type_matrix &x,
                              const value_type_matrix &b) {
        TACHO_TEST_FOR_EXCEPTION(size_t(A.NumRows()) != size_t(A.NumCols()) ||
                                 size_t(A.NumRows()) != size_t(b.dimension_0()) ||
                                 size_t(x.dimension_0()) != size_t(b.dimension_0()) ||
                                 size_t(x.dimension_1()) != size_t(b.dimension_1()), std::logic_error,
                                 "A,x and b dimensions are not compatible");
        CrsMatrixBase<value_type,host_space> h_A;
        h_A.createMirror(A); h_A.copy(A);

        auto h_x = Kokkos::create_mirror_view(x); Kokkos::deep_copy(h_x, x);
        auto h_b = Kokkos::create_mirror_view(b); Kokkos::deep_copy(h_b, b);

        typedef ArithTraits<value_type> arith_traits;
        const ordinal_type m = h_A.NumRows(), k = h_b.dimension_1();
        double diff = 0, norm = 0;
        for (ordinal_type i=0;i<m;++i) {
          for (ordinal_type p=0;p<k;++p) {
            value_type s = 0;
            const ordinal_type jbeg = h_A.RowPtrBegin(i), jend = h_A.RowPtrEnd(i);
            for (ordinal_type j=jbeg;j<jend;++j) {
              const ordinal_type col = h_A.Col(j);
              s += h_A.Value(j)*h_x(col,p);
            }
            norm += arith_traits::real(h_b(i,p)*arith_traits::conj(h_b(i,p)));
            diff += arith_traits::real((h_b(i,p) - s)*arith_traits::conj(h_b(i,p) - s));
          }
        }
        return sqrt(diff/norm);
      }

      inline
      double
      computeRelativeResidual(const value_type_matrix &x,
                              const value_type_matrix &b) {
        crs_matrix_type A;
        auto d_last = Kokkos::subview(_ap, _m);
        auto h_last = Kokkos::create_mirror_view(d_last);
        Kokkos::deep_copy(h_last, d_last);
        A.setExternalMatrix(_m, _m, h_last(), //_ap(_m),
                            _ap, _aj, _ax);

        return computeRelativeResidual(A, x, b);
      }

    };

}
#endif
