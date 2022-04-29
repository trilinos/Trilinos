#ifndef __TACHO_NUMERIC_TOOLS_LEVELSET_HPP__
#define __TACHO_NUMERIC_TOOLS_LEVELSET_HPP__

/// \file Tacho_NumericTools_LevelSet.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_NumericTools_Base.hpp"

#include "Tacho_DenseFlopCount.hpp"
#include "Tacho_SupernodeInfo.hpp"

#include "Tacho_Copy.hpp"
#include "Tacho_Copy_OnDevice.hpp"

#include "Tacho_SetIdentity.hpp"
#include "Tacho_SetIdentity_OnDevice.hpp"

#include "Tacho_Symmetrize.hpp"
#include "Tacho_Symmetrize_OnDevice.hpp"

#include "Tacho_ApplyPivots.hpp"
#include "Tacho_ApplyPivots_OnDevice.hpp"

#include "Tacho_Scale2x2_BlockInverseDiagonals.hpp"
#include "Tacho_Scale2x2_BlockInverseDiagonals_OnDevice.hpp"

#include "Tacho_Chol.hpp"
#include "Tacho_Chol_OnDevice.hpp"

#include "Tacho_LDL.hpp"
#include "Tacho_LDL_OnDevice.hpp"

#include "Tacho_Trsm.hpp"
#include "Tacho_Trsm_OnDevice.hpp"

#include "Tacho_Herk.hpp"
#include "Tacho_Herk_OnDevice.hpp"

#include "Tacho_Trsv.hpp"
#include "Tacho_Trsv_OnDevice.hpp"

#include "Tacho_Gemv.hpp"
#include "Tacho_Gemv_OnDevice.hpp"

#include "Tacho_Gemm.hpp"
#include "Tacho_Gemm_OnDevice.hpp"

#include "Tacho_GemmTriangular.hpp"
#include "Tacho_GemmTriangular_OnDevice.hpp"

#include "Tacho_SupernodeInfo.hpp"

#include "Tacho_TeamFunctor_FactorizeChol.hpp"
#include "Tacho_TeamFunctor_SolveLowerChol.hpp"
#include "Tacho_TeamFunctor_SolveUpperChol.hpp"

#include "Tacho_TeamFunctor_FactorizeLDL.hpp"
#include "Tacho_TeamFunctor_SolveLowerLDL.hpp"
#include "Tacho_TeamFunctor_SolveUpperLDL.hpp"

//#define TACHO_TEST_LEVELSET_TOOLS_KERNEL_OVERHEAD
//#define TACHO_ENABLE_LEVELSET_TOOLS_USE_LIGHT_KERNEL

namespace Tacho {

  template<typename ValueType, typename DeviceType, int Var>
  class NumericToolsLevelSet : public NumericToolsBase<ValueType,DeviceType> {
  public:
    enum { variant = Var,
           max_factor_team_size = 64 };

    ///
    /// types
    ///
    using base_type = NumericToolsBase<ValueType,DeviceType>;
    using typename base_type::value_type;
    using typename base_type::device_type;
    using typename base_type::exec_space;
    using typename base_type::exec_memory_space;
    using typename base_type::host_space;
    using typename base_type::host_memory_space;
    using typename base_type::range_type;
    using typename base_type::supernode_info_type;
    using typename base_type::ordinal_type_array;    
    using typename base_type::size_type_array;    
    using typename base_type::value_type_array;    
    using typename base_type::value_type_matrix;    
    using typename base_type::ordinal_type_array_host;    
    using typename base_type::size_type_array_host;    
    using typename base_type::supernode_type_array_host;    

    using base_type::base_type;

  private:
    using base_type::_m;
    using base_type::_ap;
    using base_type::_aj;
    using base_type::_ax;
    using base_type::_perm;
    using base_type::_peri;
    using base_type::_stree_roots;
    using base_type::_superpanel_buf;
    using base_type::_piv;
    using base_type::_diag; 
    using base_type::_info;
    using base_type::_nsupernodes;
    using base_type::_stree_level;

    using base_type::stat;
    using base_type::track_alloc;
    using base_type::track_free;
    using base_type::reset_stat;
    using base_type::print_stat_memory;
    
    // supernode host information for level kernels launching
    supernode_type_array_host _h_supernodes;

    ///
    /// level set infrastructure and tuning parameters
    ///

    // 0: device level function, 1: team policy, 2: team policy recursive
    ordinal_type _device_factorize_thres, _device_solve_thres;
    ordinal_type _device_level_cut, _team_serial_level_cut;

    ordinal_type_array_host _h_factorize_mode, _h_solve_mode;
    ordinal_type_array        _factorize_mode,   _solve_mode;

    // level details on host
    ordinal_type _nlevel;
    size_type_array_host _h_level_ptr;
    ordinal_type_array_host _h_level_sids;

    // level sids on device
    ordinal_type_array _level_sids;

    // buf level pointer
    ordinal_type_array_host _h_buf_level_ptr;

    // workspace metadata for factorization; 
    size_type_array_host _h_buf_factor_ptr;
    size_type_array _buf_factor_ptr;

    // workspace meta data for solve
    size_type_array_host _h_buf_solve_ptr, _h_buf_solve_nrhs_ptr;
    size_type_array _buf_solve_ptr, _buf_solve_nrhs_ptr;

    // workspace
    size_type _bufsize_factorize, _bufsize_solve;
    value_type_array _buf;

    // common for host and cuda
    int _status;

    // cuda stream
    int _nstreams;
#if defined(KOKKOS_ENABLE_CUDA)
    bool _is_cublas_created, _is_cusolver_dn_created;
    cublasHandle_t _handle_blas;
    cusolverDnHandle_t _handle_lapack;
    using cuda_stream_array_host = std::vector<cudaStream_t>;
    cuda_stream_array_host _cuda_streams;

    using exec_instance_array_host = std::vector<exec_space>;
    exec_instance_array_host _exec_instances;
#else 
    int _handle_blas, _handle_lapack; // dummy handle for convenience
#endif

    ///
    /// statistics
    ///
    struct {
      int n_device_factorize, n_team_factorize, n_kernel_launching_factorize;
      int n_device_solve, n_team_solve, n_kernel_launching_solve;
      int n_kernel_launching;
    } stat_level;

    ///
    /// error check for cuda things
    ///
    inline 
    void 
    checkStatus(const char *func, const char *lib) {
      if (_status != 0) {
        printf("Error: %s, %s returns non-zero status %d\n", 
               lib, func, _status);
        std::runtime_error("checkStatus failed");
      }
    }
    inline void checkDeviceLapackStatus(const char *func) { 
#if defined(KOKKOS_ENABLE_CUDA)
      constexpr bool is_host = std::is_same<exec_memory_space,Kokkos::HostSpace>::value;
      checkStatus(func, is_host ? "HostLapack" : "CuSolverDn");
#else
      checkStatus(func, "HostLapack");      
#endif  
    }
    inline void checkDeviceBlasStatus(const char *func) { 
#if defined(KOKKOS_ENABLE_CUDA)
      constexpr bool is_host = std::is_same<exec_memory_space,Kokkos::HostSpace>::value;
      checkStatus(func, is_host ? "HostBlas" : "CuBlas");
#else
      checkStatus(func, "HostBlas");      
#endif  
    }
    inline void checkDeviceStatus(const char *func) {
#if defined(KOKKOS_ENABLE_CUDA)
      constexpr bool is_host = std::is_same<exec_memory_space,Kokkos::HostSpace>::value;
      checkStatus(func, is_host ? "Host" : "Cuda");
#else
      checkStatus(func, "Host");  
#endif
    }

    inline
    void
    print_stat_init() override {
      base_type::print_stat_init();
      const double kilo(1024);
      printf("  Time\n");
      printf("             time for initialization:                         %10.6f s\n", stat.t_init);
      printf("             time for compute mode classification:            %10.6f s\n", stat.t_mode_classification);
      printf("             total time spent:                                %10.6f s\n", (stat.t_init+stat.t_mode_classification));
      printf("\n");
      printf("  Memory\n");
      printf("             workspace allocated for solve:                   %10.3f MB\n", stat.m_used/kilo/kilo);
      printf("             peak memory used:                                %10.3f MB\n", stat.m_peak/kilo/kilo);
      printf("\n");
      printf("  Compute Mode in Factorize with a Threshold(%d)\n", _device_factorize_thres);
      printf("             # of subproblems using device functions:         %6d\n", stat_level.n_device_factorize);
      printf("             # of subproblems using team functions:           %6d\n", stat_level.n_team_factorize);
      printf("             total # of subproblems:                          %6d\n", (stat_level.n_device_factorize+stat_level.n_team_factorize));
      printf("\n");
      printf("  Compute Mode in Solve with a Threshold(%d)\n", _device_solve_thres);
      printf("             # of subproblems using device functions:         %6d\n", stat_level.n_device_solve);
      printf("             # of subproblems using team functions:           %6d\n", stat_level.n_team_solve);
      printf("             total # of subproblems:                          %6d\n", (stat_level.n_device_solve+stat_level.n_team_solve));
      printf("\n");
    }

    inline
    void
    print_stat_factor() override {
      base_type::print_stat_factor();
      double flop = 0;
      switch (this->getSolutionMethod()) {
      case 1: {
        for (ordinal_type sid=0;sid<_nsupernodes;++sid) {
          auto &s = _h_supernodes(sid);
          const ordinal_type m = s.m, n = s.n - s.m;
          flop += DenseFlopCount<value_type>::Chol(m);
          if (variant == 1) { 
            flop += DenseFlopCount<value_type>::Trsm(true,  m, m);
          } 
          else if (variant == 2) {
            flop += DenseFlopCount<value_type>::Trsm(true,  m, m);
            flop += DenseFlopCount<value_type>::Trsm(true,  m, n);
          }
          flop += DenseFlopCount<value_type>::Trsm(true,  m, n);
          flop += DenseFlopCount<value_type>::Syrk(m, n);
        }
        break;
      }
      case 2: {
        for (ordinal_type sid=0;sid<_nsupernodes;++sid) {
          auto &s = _h_supernodes(sid);
          const ordinal_type m = s.m, n = s.n - s.m;
          flop += DenseFlopCount<value_type>::LDL(m);
          flop += DenseFlopCount<value_type>::Trsm(true,  m, n);
          flop += DenseFlopCount<value_type>::Syrk(m, n);
        }
        break;
      }
      default: {
        TACHO_TEST_FOR_EXCEPTION(false,
                                 std::logic_error,
                                 "The solution method is not supported");
      }
      }
      const double kilo(1024);
      printf("  FLOPs\n");
      printf("             gflop   for numeric factorization:               %10.3f GFLOP\n", flop/kilo/kilo/kilo);
      printf("             gflop/s for numeric factorization:               %10.3f GFLOP/s\n", flop/stat.t_factor/kilo/kilo/kilo);
      printf("\n");
      printf("  Kernels\n");
      printf("             # of kernels launching:                          %6d\n", stat_level.n_kernel_launching);
      printf("\n");
    }

    inline
    void
    print_stat_solve() override {
      base_type::print_stat_solve();
      printf("  Kernels\n");
      printf("             # of kernels launching:                          %6d\n", stat_level.n_kernel_launching);
      printf("\n");
    }

  public:
    ///
    /// initialization / release
    ///
    inline
    void
    initialize(const ordinal_type device_level_cut,
               const ordinal_type device_factorize_thres,
               const ordinal_type device_solve_thres,
               const ordinal_type verbose = 0) {
      stat_level.n_device_factorize = 0;   stat_level.n_device_solve = 0;
      stat_level.n_team_factorize= 0;      stat_level.n_team_solve = 0;

      Kokkos::Timer timer;

      timer.reset();

      ///
      /// level data structure
      ///

      // # of supernodes
      _nsupernodes = _info.supernodes.extent(0);

      // local host supernodes info 
      _h_supernodes = Kokkos::create_mirror_view_and_copy(host_memory_space(), _info.supernodes);

      // # of levels
      _nlevel = 0;
      {
        for (ordinal_type sid=0;sid<_nsupernodes;++sid) 
          _nlevel = max(_stree_level(sid), _nlevel);
        ++_nlevel;
      }

      // create level ptr
      _h_level_ptr = size_type_array_host("h_level_ptr", _nlevel+1);
      {
        // first count # of supernodes in each level
        for (ordinal_type sid=0;sid<_nsupernodes;++sid) 
          ++_h_level_ptr(_stree_level(sid)+1);

        // scan 
        for (ordinal_type i=0;i<_nlevel;++i) 
          _h_level_ptr(i+1) += _h_level_ptr(i);
      }

      // fill sids
      _h_level_sids = ordinal_type_array_host(do_not_initialize_tag("h_level_sids"), _nsupernodes);
      {
        size_type_array_host tmp_level_ptr(do_not_initialize_tag("tmp_level_ptr"), _h_level_ptr.extent(0));
        Kokkos::deep_copy(tmp_level_ptr, _h_level_ptr);
        for (ordinal_type sid=0;sid<_nsupernodes;++sid) {
          const ordinal_type lvl = _stree_level(sid);
          _h_level_sids(tmp_level_ptr(lvl)++) = sid;
        }
      }
      _level_sids = Kokkos::create_mirror_view_and_copy(exec_memory_space(), _h_level_sids); 
      track_alloc(_level_sids.span()*sizeof(ordinal_type));

      ///
      /// workspace
      ///
      _h_buf_level_ptr = ordinal_type_array_host(do_not_initialize_tag("h_buf_factor_level_ptr"), _nlevel+1);
      {
        _h_buf_level_ptr(0) = 0;
        for (ordinal_type i=0;i<_nlevel;++i) {
          const ordinal_type pbeg = _h_level_ptr(i), pend = _h_level_ptr(i+1);
          _h_buf_level_ptr(i+1) = (pend - pbeg + 1) + _h_buf_level_ptr(i);
        }
      }

      // create workspace for factorization / solve
      _bufsize_factorize = 0;
      _bufsize_solve = 0;
      _h_buf_factor_ptr = size_type_array_host(do_not_initialize_tag("h_buf_factor_ptr"), _h_buf_level_ptr(_nlevel));
      _h_buf_solve_ptr = size_type_array_host(do_not_initialize_tag("h_buf_solve_ptr"), _h_buf_level_ptr(_nlevel)); 
      {
        for (ordinal_type i=0;i<_nlevel;++i) {
          const ordinal_type lbeg = _h_buf_level_ptr(i);
          const ordinal_type pbeg = _h_level_ptr(i), pend = _h_level_ptr(i+1);

          _h_buf_factor_ptr(lbeg) = 0;
          _h_buf_solve_ptr(lbeg) = 0;
          for (ordinal_type p=pbeg,k=(lbeg+1);p<pend;++p,++k) {
            const ordinal_type sid = _h_level_sids(p);
            const auto s = _h_supernodes(sid); 
            const ordinal_type m = s.m, n = s.n, n_m = n-m;
            const ordinal_type schur_work_size = n_m*(n_m+max_factor_team_size);            
            const ordinal_type chol_factor_work_size_variants[3] = { schur_work_size, max(m*m, schur_work_size), m*m+schur_work_size };
            const ordinal_type chol_factor_work_size = chol_factor_work_size_variants[variant];
            const ordinal_type ldl_factor_work_size = chol_factor_work_size_variants[0] + max(32*m, m*n);
            const ordinal_type factor_work_size_variants[2] = { chol_factor_work_size, ldl_factor_work_size };

            const ordinal_type chol_solve_work_size = (variant == 0 ? n_m : n);
            const ordinal_type ldl_solve_work_size = n_m;
            const ordinal_type solve_work_size_variants[2] = { chol_solve_work_size, ldl_solve_work_size };

            const ordinal_type index_work_size = this->getSolutionMethod()-1;
            const ordinal_type factor_work_size = factor_work_size_variants[index_work_size];
            const ordinal_type solve_work_size = solve_work_size_variants[index_work_size];

            _h_buf_factor_ptr(k) = factor_work_size + _h_buf_factor_ptr(k-1);
            _h_buf_solve_ptr(k) = solve_work_size + _h_buf_solve_ptr(k-1);
          }
          const ordinal_type last_idx = lbeg+pend-pbeg;
          _bufsize_factorize = max(_bufsize_factorize, _h_buf_factor_ptr(last_idx));
          _bufsize_solve = max(_bufsize_solve, _h_buf_solve_ptr(last_idx));
        }
      }

      _buf_factor_ptr = Kokkos::create_mirror_view_and_copy(exec_memory_space(), _h_buf_factor_ptr);
      track_alloc(_buf_factor_ptr.span()*sizeof(size_type));

      _buf_solve_ptr = Kokkos::create_mirror_view_and_copy(exec_memory_space(), _h_buf_solve_ptr);
      track_alloc(_buf_solve_ptr.span()*sizeof(size_type));

      _h_buf_solve_nrhs_ptr = size_type_array_host(do_not_initialize_tag("h_buf_solve_nrhs_ptr"), _h_buf_solve_ptr.extent(0));
      _buf_solve_nrhs_ptr = Kokkos::create_mirror_view(exec_memory_space(), _h_buf_solve_nrhs_ptr);
      track_alloc(_buf_solve_nrhs_ptr.span()*sizeof(size_type));

      ///
      /// cuda library initialize
      ///
#if defined(KOKKOS_ENABLE_CUDA)
      if (!_is_cublas_created) {
        _status = cublasCreate(&_handle_blas); checkDeviceBlasStatus("cublasCreate"); _is_cublas_created = true;
      }
      if (!_is_cusolver_dn_created) {
        _status = cusolverDnCreate(&_handle_lapack); checkDeviceLapackStatus("cusolverDnCreate"); _is_cusolver_dn_created = true;
      }
#endif
      stat.t_init = timer.seconds();

      ///
      /// classification of problems
      ///
      timer.reset();

      _device_level_cut = min(device_level_cut, _nlevel);
      _device_factorize_thres = device_factorize_thres;
      _device_solve_thres = device_solve_thres;

      _h_factorize_mode = ordinal_type_array_host(do_not_initialize_tag("h_factorize_mode"), _nsupernodes);
      Kokkos::deep_copy(_h_factorize_mode, -1);

      _h_solve_mode = ordinal_type_array_host(do_not_initialize_tag("h_solve_mode"), _nsupernodes);
      Kokkos::deep_copy(_h_solve_mode, -1);

      if (_device_level_cut > 0) {
        for (ordinal_type lvl=0;lvl<_device_level_cut;++lvl) {
          const ordinal_type 
            pbeg = _h_level_ptr(lvl), 
            pend = _h_level_ptr(lvl+1);
          for (ordinal_type p=pbeg;p<pend;++p) {
            const ordinal_type sid = _h_level_sids(p);
            _h_solve_mode(sid) = 0;
            _h_factorize_mode(sid) = 0;
            ++stat_level.n_device_solve;
            ++stat_level.n_device_factorize;
          }
        }
      }

      _team_serial_level_cut = _nlevel;
      {        
        for (ordinal_type lvl=_device_level_cut;lvl<_team_serial_level_cut;++lvl) {          
          const ordinal_type 
            pbeg = _h_level_ptr(lvl), 
            pend = _h_level_ptr(lvl+1);
          for (ordinal_type p=pbeg;p<pend;++p) {
            const ordinal_type sid = _h_level_sids(p);
            const auto s = _h_supernodes(sid); 
            const ordinal_type m = s.m; //, n_m = s.n-s.m;
            if (m > _device_solve_thres) {// || n > _device_solve_thres) {
              _h_solve_mode(sid) = 0;
              ++stat_level.n_device_solve;
            } else {
              _h_solve_mode(sid) = 1;
              ++stat_level.n_team_solve;
            }
            if (m > _device_factorize_thres) {// || n_m > _device_factorize_thres) {
              _h_factorize_mode(sid) = 0;
              ++stat_level.n_device_factorize;
            } else {
              _h_factorize_mode(sid) = 1;
              ++stat_level.n_team_factorize;
            }
          }
        }
      }

      _factorize_mode = Kokkos::create_mirror_view_and_copy(exec_memory_space(), _h_factorize_mode);
      track_alloc(_factorize_mode.span()*sizeof(ordinal_type));

      _solve_mode = Kokkos::create_mirror_view_and_copy(exec_memory_space(), _h_solve_mode);
      track_alloc(_solve_mode.span()*sizeof(ordinal_type));

      stat.t_mode_classification = timer.seconds();
      if (verbose) {
        switch (this->getSolutionMethod()) {
        case 1: {
          printf("Summary: LevelSetTools-Variant-%d (InitializeCholesky)\n", variant);
          printf("======================================================\n");
          break;
        } 
        case 2: {
          printf("Summary: LevelSetTools (InitializeLDL)\n");
          printf("======================================================\n");          
          break;
        }
        }
        print_stat_init();
      }
    }

    inline
    void
    release(const ordinal_type verbose = 0) override {
      base_type::release(false);
      track_free(_buf_factor_ptr.span()*sizeof(size_type));
      track_free(_buf_solve_ptr.span()*sizeof(size_type));
      track_free(_buf_solve_nrhs_ptr.span()*sizeof(size_type));
      track_free(_buf.span()*sizeof(value_type));
      track_free(_factorize_mode.span()*sizeof(ordinal_type));
      track_free(_solve_mode.span()*sizeof(ordinal_type));
      track_free(_level_sids.span()*sizeof(ordinal_type));
      if (verbose) {
        printf("Summary: LevelSetTools-Variant-%d (Release)\n", variant);
        printf("============================================\n");
        print_stat_memory();
      }
    }

    NumericToolsLevelSet() : base_type() {
      _nlevel = 0;
      _bufsize_factorize = 0;
      _bufsize_solve = 0;
      _nstreams = 0;
      stat_level = stat_level();
    }
    NumericToolsLevelSet(const NumericToolsLevelSet &b) = default;

    NumericToolsLevelSet(const ordinal_type method,
                         // input matrix A
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
                         const ordinal_type_array_host &stree_level,
                         const ordinal_type_array_host &stree_roots)
    : base_type(method, 
                m, ap, aj,
                perm,peri,
                nsupernodes, supernodes,
                gid_ptr, gid_colidx, 
                sid_ptr, sid_colidx, blk_colidx,
                stree_parent, stree_ptr, stree_children,
                stree_level, stree_roots) {
      _nstreams = 0;
#if defined(KOKKOS_ENABLE_CUDA)
      _is_cublas_created = 0;
      _is_cusolver_dn_created = 0;
#endif
    }

    virtual~NumericToolsLevelSet() {
#if defined(KOKKOS_ENABLE_CUDA)
      // destroy previously created streams
      for (ordinal_type i=0;i<_nstreams;++i) {
        _status = cudaStreamDestroy(_cuda_streams[i]); checkDeviceStatus("cudaStreamDestroy");
      }
      _cuda_streams.clear();
      _exec_instances.clear();

      if (_is_cublas_created) {
        _status = cusolverDnDestroy(_handle_lapack); checkDeviceLapackStatus("cusolverDnDestroy");
      }
      if (_is_cusolver_dn_created) {
        _status = cublasDestroy(_handle_blas); checkDeviceBlasStatus("cublasDestroy");
      }
#endif
    }

    inline
    void
    createStream(const ordinal_type nstreams, const ordinal_type verbose = 0) {
#if defined(KOKKOS_ENABLE_CUDA)
      // destroy previously created streams
      for (ordinal_type i=0;i<_nstreams;++i) {
        _status = cudaStreamDestroy(_cuda_streams[i]); checkDeviceStatus("cudaStreamDestroy");
      }
      // new streams
      _nstreams = nstreams;
      //_cuda_streams = cuda_stream_array_host(do_not_initialize_tag("cuda streams"), _nstreams);
      _cuda_streams.clear();
      _cuda_streams.resize(_nstreams);
      for (ordinal_type i=0;i<_nstreams;++i) {
        _status = cudaStreamCreateWithFlags(&_cuda_streams[i], cudaStreamNonBlocking); checkDeviceStatus("cudaStreamCreate");
      }

      _exec_instances.clear();
      _exec_instances.resize(_nstreams);
      for (ordinal_type i=0;i<_nstreams;++i) {
        ExecSpaceFactory<exec_space>::createInstance(_cuda_streams[i], _exec_instances[i]);
      }
      if (verbose) {
        printf("Summary: CreateStream : %3d\n", _nstreams);
        printf("===========================\n");          
      }
#endif
    }

    ///
    /// Device level functions
    ///
    inline
    void
    factorizeCholeskyOnDeviceVar0(const ordinal_type pbeg, 
                                  const ordinal_type pend,
                                  const size_type_array_host &h_buf_factor_ptr,
                                  const value_type_array &work) {
      const value_type one(1), minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA)
      ordinal_type q(0);
#endif
      exec_space exec_instance;
      for (ordinal_type p=pbeg;p<pend;++p) {
        const ordinal_type sid = _h_level_sids(p);
        if (_h_factorize_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA)
          const ordinal_type qid = q%_nstreams;
          const auto mystream = _cuda_streams[qid];
          _status = cublasSetStream(_handle_blas, mystream); checkDeviceBlasStatus("cublasSetStream");
          _status = cusolverDnSetStream(_handle_lapack, mystream); checkDeviceLapackStatus("cusolverDnSetStream");

          exec_instance = _exec_instances[qid];

          const size_type worksize = work.extent(0)/_nstreams;
          value_type_array W(work.data() + worksize*qid, worksize);
          ++q;
#else
          value_type_array W = work;
#endif          
          const auto &s = _h_supernodes(sid);
          {
            const ordinal_type m = s.m, n = s.n, n_m = n-m;
            if (m > 0) {
              value_type *aptr = s.buf;
              UnmanagedViewType<value_type_matrix> ATL(aptr, m, m); aptr += m*m;
              _status = Chol<Uplo::Upper,Algo::OnDevice>
                ::invoke(_handle_lapack, ATL, W); checkDeviceLapackStatus("chol");

              if (n_m > 0) {
                exec_instance.fence();
                UnmanagedViewType<value_type_matrix> ABR(_buf.data()+h_buf_factor_ptr(p-pbeg), n_m, n_m); 
                UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // aptr += m*n_m;
                _status = Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,Algo::OnDevice>
                  ::invoke(_handle_blas, Diag::NonUnit(), one, ATL, ATR); checkDeviceBlasStatus("trsm");
                exec_instance.fence();
                _status = Herk<Uplo::Upper,Trans::ConjTranspose,Algo::OnDevice>
                  ::invoke(_handle_blas, minus_one, ATR, zero, ABR);
                exec_instance.fence();
              }
            }
          }
        }
      }
    }

    inline
    void
    factorizeCholeskyOnDeviceVar1(const ordinal_type pbeg, 
                                  const ordinal_type pend,
                                  const size_type_array_host &h_buf_factor_ptr,
                                  const value_type_array &work) {
      const value_type one(1), minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA)
      ordinal_type q(0);
#endif
      exec_space exec_instance;
      for (ordinal_type p=pbeg;p<pend;++p) {
        const ordinal_type sid = _h_level_sids(p);
        if (_h_factorize_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA)
          const ordinal_type qid = q%_nstreams;
          const auto mystream = _cuda_streams[qid];
          _status = cublasSetStream(_handle_blas, mystream); checkDeviceBlasStatus("cublasSetStream");
          _status = cusolverDnSetStream(_handle_lapack, mystream); checkDeviceLapackStatus("cusolverDnSetStream");

          exec_instance = _exec_instances[qid];

          const size_type worksize = work.extent(0)/_nstreams;
          value_type_array W(work.data() + worksize*qid, worksize);
          ++q;
#else
          value_type_array W = work;
#endif          
          const auto &s = _h_supernodes(sid);
          {
            const ordinal_type m = s.m, n = s.n, n_m = n-m;
            if (m > 0) {
              value_type *aptr = s.buf;
              UnmanagedViewType<value_type_matrix> ATL(aptr, m, m); aptr += m*m;
              _status = Chol<Uplo::Upper,Algo::OnDevice>
                ::invoke(_handle_lapack, ATL, W); checkDeviceLapackStatus("chol");

              value_type *bptr = _buf.data()+h_buf_factor_ptr(p-pbeg);
              UnmanagedViewType<value_type_matrix> T(bptr, m, m);
              _status = SetIdentity<Algo::OnDevice>::invoke(exec_instance, T, one); checkDeviceBlasStatus("SetIdentity");
              exec_instance.fence();
              _status = Trsm<Side::Left,Uplo::Upper,Trans::NoTranspose,Algo::OnDevice>
                ::invoke(_handle_blas, Diag::NonUnit(), one, ATL, T); checkDeviceBlasStatus("trsm");

              if (n_m > 0) {
                exec_instance.fence();
                UnmanagedViewType<value_type_matrix> ABR(bptr, n_m, n_m); 
                UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // aptr += m*n_m;
                _status = Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,Algo::OnDevice>
                  ::invoke(_handle_blas, Diag::NonUnit(), one, ATL, ATR); checkDeviceBlasStatus("trsm");
                exec_instance.fence();
                _status = Copy<Algo::OnDevice>::invoke(exec_instance, ATL, T); checkDeviceBlasStatus("Copy");
                exec_instance.fence();
                _status = Herk<Uplo::Upper,Trans::ConjTranspose,Algo::OnDevice>
                  ::invoke(_handle_blas, minus_one, ATR, zero, ABR);
                exec_instance.fence();
              } else {
                exec_instance.fence();
                _status = Copy<Algo::OnDevice>::invoke(exec_instance, ATL, T); checkDeviceBlasStatus("Copy");
              }
            }
          }
        }
      }
    }

    inline
    void
    factorizeCholeskyOnDeviceVar2(const ordinal_type pbeg, 
                                  const ordinal_type pend,
                                  const size_type_array_host &h_buf_factor_ptr,
                                  const value_type_array &work) {
      const value_type one(1), minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA)
      ordinal_type q(0);
#endif
      exec_space exec_instance;
      for (ordinal_type p=pbeg;p<pend;++p) {
        const ordinal_type sid = _h_level_sids(p);
        if (_h_factorize_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA)
          const ordinal_type qid = q%_nstreams;
          const auto mystream = _cuda_streams[qid];
          _status = cublasSetStream(_handle_blas, mystream); checkDeviceBlasStatus("cublasSetStream");
          _status = cusolverDnSetStream(_handle_lapack, mystream); checkDeviceLapackStatus("cusolverDnSetStream");

          exec_instance = _exec_instances[qid];

          const size_type worksize = work.extent(0)/_nstreams;
          value_type_array W(work.data() + worksize*qid, worksize);
          ++q;
#else
          value_type_array W = work;
#endif          
          const auto &s = _h_supernodes(sid);
          {
            const ordinal_type m = s.m, n = s.n, n_m = n-m;
            if (m > 0) {
              value_type *aptr = s.buf;
              UnmanagedViewType<value_type_matrix> ATL(aptr, m, m); aptr += m*m;
              _status = Chol<Uplo::Upper,Algo::OnDevice>
                ::invoke(_handle_lapack, ATL, W); checkDeviceLapackStatus("chol");

              value_type *bptr = _buf.data()+h_buf_factor_ptr(p-pbeg);
              if (n_m > 0) {
                exec_instance.fence();
                UnmanagedViewType<value_type_matrix> ABR(bptr, n_m, n_m); bptr += ABR.span();
                UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // aptr += m*n_m;
                
                _status = Trsm<Side::Left,Uplo::Upper,Trans::ConjTranspose,Algo::OnDevice>
                  ::invoke(_handle_blas, Diag::NonUnit(), one, ATL, ATR); checkDeviceBlasStatus("trsm");
                exec_instance.fence();
                _status = Herk<Uplo::Upper,Trans::ConjTranspose,Algo::OnDevice>
                  ::invoke(_handle_blas, minus_one, ATR, zero, ABR);
                exec_instance.fence();                

                /// additional things
                UnmanagedViewType<value_type_matrix> T(bptr, m, m); 
                _status = Copy<Algo::OnDevice>::invoke(exec_instance, T, ATL); checkDeviceBlasStatus("Copy");
                exec_instance.fence();
                _status = SetIdentity<Algo::OnDevice>::invoke(exec_instance, ATL, minus_one); checkDeviceBlasStatus("SetIdentity");
                exec_instance.fence();

                UnmanagedViewType<value_type_matrix> AT(ATL.data(), m, n);
                _status = Trsm<Side::Left,Uplo::Upper,Trans::NoTranspose,Algo::OnDevice>
                  ::invoke(_handle_blas, Diag::NonUnit(), minus_one, T, AT); checkDeviceBlasStatus("trsm");
                exec_instance.fence();
              } else {
                exec_instance.fence();
                /// additional things
                UnmanagedViewType<value_type_matrix> T(bptr, m, m); 
                _status = Copy<Algo::OnDevice>::invoke(exec_instance, T, ATL); checkDeviceBlasStatus("Copy");
                exec_instance.fence();
                _status = SetIdentity<Algo::OnDevice>::invoke(exec_instance, ATL, one); checkDeviceBlasStatus("SetIdentity");
                exec_instance.fence();
                _status = Trsm<Side::Left,Uplo::Upper,Trans::NoTranspose,Algo::OnDevice>
                  ::invoke(_handle_blas, Diag::NonUnit(), one, T, ATL); checkDeviceBlasStatus("trsm");
              }
            }
          }
        }
      }
    }

    inline
    void
    factorizeCholeskyOnDevice(const ordinal_type pbeg, 
                              const ordinal_type pend,
                              const size_type_array_host &h_buf_factor_ptr,
                              const value_type_array &work) {
      if      (variant == 0) 
        factorizeCholeskyOnDeviceVar0(pbeg, pend, h_buf_factor_ptr, work); 
      else if (variant == 1)
        factorizeCholeskyOnDeviceVar1(pbeg, pend, h_buf_factor_ptr, work); 
      else if (variant == 2)
        factorizeCholeskyOnDeviceVar2(pbeg, pend, h_buf_factor_ptr, work); 
      else {
        TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, 
                                 "LevelSetTools::factorizeCholeskyOnDevice, algorithm variant is not supported");
      }
    }

    inline
    void
    factorizeLDL_OnDevice(const ordinal_type pbeg, 
                          const ordinal_type pend,
                          const size_type_array_host &h_buf_factor_ptr,
                          const value_type_array &work) {
      const value_type one(1), minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA)
      ordinal_type q(0);
#endif
      exec_space exec_instance;
      for (ordinal_type p=pbeg;p<pend;++p) {
        const ordinal_type sid = _h_level_sids(p);
        if (_h_factorize_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA)
          const ordinal_type qid = q%_nstreams;
          const auto mystream = _cuda_streams[qid];
          _status = cublasSetStream(_handle_blas, mystream); checkDeviceBlasStatus("cublasSetStream");
          _status = cusolverDnSetStream(_handle_lapack, mystream); checkDeviceLapackStatus("cusolverDnSetStream");

          exec_instance = _exec_instances[qid];

          const size_type worksize = work.extent(0)/_nstreams;
          value_type_array W(work.data() + worksize*qid, worksize);
          ++q;
#else
          value_type_array W = work;
#endif          
          const auto &s = _h_supernodes(sid);
          {
            const ordinal_type offs = s.row_begin, m = s.m, n = s.n, n_m = n-m;
            if (m > 0) {
              value_type *aptr = s.buf;
              UnmanagedViewType<value_type_matrix> ATL(aptr, m, m); aptr += m*m;

              _status = Symmetrize<Uplo::Upper,Algo::OnDevice>::invoke(exec_instance, ATL); 
              exec_instance.fence();

              ordinal_type *pivptr = _piv.data() + 4*offs;
              UnmanagedViewType<ordinal_type_array> P(pivptr, 4*m);
              _status = LDL<Uplo::Lower,Algo::OnDevice>::invoke(_handle_lapack, ATL, P, W); checkDeviceLapackStatus("ldl::invoke");
              exec_instance.fence();

              value_type * dptr = _diag.data() + 2*offs;
              UnmanagedViewType<value_type_matrix> D(dptr, m, 2);              
              _status = LDL<Uplo::Lower,Algo::OnDevice>::modify(exec_instance, ATL, P, D); checkDeviceLapackStatus("ldl::modify");
              exec_instance.fence();

              if (n_m > 0) {
                UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // aptr += m*n_m;
                UnmanagedViewType<value_type_matrix> ABR(_buf.data()+h_buf_factor_ptr(p-pbeg), n_m, n_m); 
                UnmanagedViewType<value_type_matrix> STR(ABR.data()+ABR.span(), m, n_m); 

                auto fpiv = ordinal_type_array(P.data()+m, m);
                _status = ApplyPivots<PivotMode::Flame,Side::Left,Direct::Forward,Algo::OnDevice>
                  ::invoke(exec_instance, fpiv, ATR); 
                exec_instance.fence();

                _status = Trsm<Side::Left,Uplo::Lower,Trans::NoTranspose,Algo::OnDevice>
                  ::invoke(_handle_blas, Diag::Unit(), one, ATL, ATR); checkDeviceBlasStatus("trsm");
                exec_instance.fence();

                _status = Copy<Algo::OnDevice>
                  ::invoke(exec_instance, STR, ATR);
                exec_instance.fence();

                _status = Scale2x2_BlockInverseDiagonals<Side::Left,Algo::OnDevice>
                  ::invoke(exec_instance, P, D, ATR);
                exec_instance.fence();

                _status = GemmTriangular<Trans::Transpose,Trans::NoTranspose,Uplo::Upper,Algo::OnDevice>
                  ::invoke(_handle_blas, minus_one, ATR, STR, zero, ABR);
                exec_instance.fence(); checkDeviceBlasStatus("gemm");
              }
            }
          }
        }
      }
    }


    ///
    /// Level set factorize
    ///
    inline 
    void
    factorizeCholesky(const value_type_array &ax,
                      const ordinal_type verbose) {
      constexpr bool is_host = std::is_same<exec_memory_space,Kokkos::HostSpace>::value;
      Kokkos::Timer timer;

      timer.reset();
      value_type_array work;
      {
        _buf = value_type_array(do_not_initialize_tag("buf"), _bufsize_factorize);
        track_alloc(_buf.span()*sizeof(value_type));
        
#if defined (KOKKOS_ENABLE_CUDA)
        value_type_matrix T(NULL, _info.max_supernode_size, _info.max_supernode_size);
        const size_type worksize = Chol<Uplo::Upper,Algo::OnDevice>
          ::invoke(_handle_lapack, T, work); 
        
        work = value_type_array(do_not_initialize_tag("work"), worksize*(_nstreams+1));
        track_alloc(work.span()*sizeof(value_type));
#endif
      }
      stat.t_extra = timer.seconds();

      timer.reset();
      {
        _ax = ax; // matrix values
        _info.copySparseToSuperpanels(_ap, _aj, _ax, _perm, _peri);
      }
      stat.t_copy = timer.seconds();

      stat_level.n_kernel_launching = 0;
      timer.reset();
      { 
        // this should be considered with average problem sizes in levels
        const ordinal_type half_level = _nlevel/2;
        //const ordinal_type team_size_factor[2] = { 64, 16 }, vector_size_factor[2] = { 8, 8};
        //const ordinal_type team_size_factor[2] = { 16, 16 }, vector_size_factor[2] = { 32, 32};
        const ordinal_type team_size_factor[2] = { 64, 64 }, vector_size_factor[2] = { 8, 4};
        const ordinal_type team_size_update[2] = { 16, 8 }, vector_size_update[2] = { 32, 32};
        {
          typedef TeamFunctor_FactorizeChol<supernode_info_type> functor_type;
#if defined(TACHO_TEST_LEVELSET_TOOLS_KERNEL_OVERHEAD)
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::DummyTag> team_policy_factorize;
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::DummyTag> team_policy_update;
#else
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::template FactorizeTag<variant> > team_policy_factor;
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::UpdateTag> team_policy_update;
#endif

          functor_type functor(_info, 
                               _factorize_mode,
                               _level_sids,
                               _buf);

          team_policy_factor policy_factor(1,1,1);
          team_policy_update policy_update(1,1,1);

          {
            for (ordinal_type lvl=(_team_serial_level_cut-1);lvl>=0;--lvl) {
              const ordinal_type 
                pbeg = _h_level_ptr(lvl), 
                pend = _h_level_ptr(lvl+1),
                pcnt = pend - pbeg;

              const range_type range_buf_factor_ptr(_h_buf_level_ptr(lvl), _h_buf_level_ptr(lvl+1));

              const auto buf_factor_ptr = Kokkos::subview(_buf_factor_ptr, range_buf_factor_ptr);
              functor.setRange(pbeg, pend);
              functor.setBufferPtr(buf_factor_ptr);
              if (is_host) {
                policy_factor = team_policy_factor(pcnt, 1, 1);
                policy_update = team_policy_update(pcnt, 1, 1);
              } else {
                const ordinal_type idx = lvl > half_level;
                policy_factor = team_policy_factor(pcnt, team_size_factor[idx], vector_size_factor[idx]);
                policy_update = team_policy_update(pcnt, team_size_update[idx], vector_size_update[idx]);
              }
              if (lvl < _device_level_cut) {
                // do nothing
                //Kokkos::parallel_for("factor lower", policy_factor, functor);
              } else {
                Kokkos::parallel_for("factor", policy_factor, functor);
                ++stat_level.n_kernel_launching;
              }

              const auto h_buf_factor_ptr = Kokkos::subview(_h_buf_factor_ptr, range_buf_factor_ptr);
              factorizeCholeskyOnDevice(pbeg, pend, h_buf_factor_ptr, work); 
              Kokkos::fence();

              Kokkos::parallel_for("update factor", policy_update, functor); 
              ++stat_level.n_kernel_launching;
              exec_space().fence(); //Kokkos::fence();
            }
          }
        } 
      } // end of Cholesky
      stat.t_factor = timer.seconds();

      timer.reset();
      {
#if defined (KOKKOS_ENABLE_CUDA)
        track_free(work.span()*sizeof(value_type));
#endif
        track_free(_buf.span()*sizeof(value_type));
        _buf = value_type_array();
      }
      stat.t_extra += timer.seconds();

      if (verbose) {
        printf("Summary: LevelSetTools-Variant-%d (CholeskyFactorize)\n", variant);
        printf("=====================================================\n");
        print_stat_factor();
      }
    }

    inline
    void
    solveCholeskyLowerOnDeviceVar0(const ordinal_type pbeg, 
                                   const ordinal_type pend,
                                   const size_type_array_host &h_buf_solve_ptr,
                                   const value_type_matrix &t) {
      const ordinal_type nrhs = t.extent(1);
      const value_type minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA)
      ordinal_type q(0);
#endif
      for (ordinal_type p=pbeg;p<pend;++p) {
        const ordinal_type sid = _h_level_sids(p);
        if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA)
          _status = cublasSetStream(_handle_blas, _cuda_streams[q%_nstreams]); checkDeviceStatus("cublasSetStream");
          ++q;
#endif          
          const auto &s = _h_supernodes(sid);
          {
            const ordinal_type m = s.m, n = s.n, n_m = n-m;
            if (m > 0) {
              value_type *aptr = s.buf;
              UnmanagedViewType<value_type_matrix> AL(aptr, m, m); aptr += m*m;

              const ordinal_type offm = s.row_begin;
              auto tT = Kokkos::subview(t, range_type(offm, offm+m), Kokkos::ALL());
              _status = Trsv<Uplo::Upper,Trans::ConjTranspose,Algo::OnDevice>
                ::invoke(_handle_blas, Diag::NonUnit(), AL, tT); checkDeviceBlasStatus("trsv");

              if (n_m > 0) {
                // solve offdiag
                value_type *bptr = _buf.data()+h_buf_solve_ptr(p-pbeg);
                UnmanagedViewType<value_type_matrix> AR(aptr, m, n_m); // aptr += m*n_m;
                UnmanagedViewType<value_type_matrix> bB(bptr, n_m, nrhs);
                _status = Gemv<Trans::ConjTranspose,Algo::OnDevice>
                  ::invoke(_handle_blas, minus_one, AR, tT, zero, bB); checkDeviceBlasStatus("gemv");
              }
            }
          }
        }
      }
    }

    inline
    void
    solveCholeskyLowerOnDeviceVar1(const ordinal_type pbeg, 
                                   const ordinal_type pend,
                                   const size_type_array_host &h_buf_solve_ptr,
                                   const value_type_matrix &t) {
      const ordinal_type nrhs = t.extent(1);
      const value_type one(1), minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA)
      ordinal_type q(0);
#endif
      for (ordinal_type p=pbeg;p<pend;++p) {
        const ordinal_type sid = _h_level_sids(p);
        if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA)
          _status = cublasSetStream(_handle_blas, _cuda_streams[q%_nstreams]); checkDeviceStatus("cublasSetStream");
          ++q;
#endif          
          const auto &s = _h_supernodes(sid);
          {
            const ordinal_type m = s.m, n = s.n, n_m = n-m;
            if (m > 0) {
              value_type *aptr = s.buf;
              UnmanagedViewType<value_type_matrix> AL(aptr, m, m); aptr += m*m;

              value_type *bptr = _buf.data()+h_buf_solve_ptr(p-pbeg);
              UnmanagedViewType<value_type_matrix> bT(bptr, m, nrhs); bptr += m*nrhs; 

              const ordinal_type offm = s.row_begin;
              auto tT = Kokkos::subview(t, range_type(offm, offm+m), Kokkos::ALL());

              _status = Gemv<Trans::ConjTranspose,Algo::OnDevice>
                ::invoke(_handle_blas, one, AL, tT, zero, bT); checkDeviceBlasStatus("gemv");

              if (n_m > 0) {
                // solve offdiag
                UnmanagedViewType<value_type_matrix> AR(aptr, m, n_m); 
                UnmanagedViewType<value_type_matrix> bB(bptr, n_m, nrhs); 

                _status = Gemv<Trans::ConjTranspose,Algo::OnDevice>
                  ::invoke(_handle_blas, minus_one, AR, bT, zero, bB); checkDeviceBlasStatus("gemv");
              }
            }
          }
        }
      }
    }

    inline
    void
    solveCholeskyLowerOnDeviceVar2(const ordinal_type pbeg, 
                                   const ordinal_type pend,
                                   const size_type_array_host &h_buf_solve_ptr,
                                   const value_type_matrix &t) {
      const ordinal_type nrhs = t.extent(1);
      const value_type one(1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA)
      ordinal_type q(0);
#endif
      for (ordinal_type p=pbeg;p<pend;++p) {
        const ordinal_type sid = _h_level_sids(p);
        if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA)
          _status = cublasSetStream(_handle_blas, _cuda_streams[q%_nstreams]); checkDeviceStatus("cublasSetStream");
          ++q;
#endif          
          const auto &s = _h_supernodes(sid);
          {
            const ordinal_type m = s.m, n = s.n;
            if (m > 0) {
              value_type *aptr = s.buf;
              UnmanagedViewType<value_type_matrix> A(aptr, m, n); 

              value_type *bptr = _buf.data()+h_buf_solve_ptr(p-pbeg);
              UnmanagedViewType<value_type_matrix> b(bptr, n, nrhs); 

              const ordinal_type offm = s.row_begin;
              auto tT = Kokkos::subview(t, range_type(offm, offm+m), Kokkos::ALL());

              _status = Gemv<Trans::ConjTranspose,Algo::OnDevice>
                ::invoke(_handle_blas, one, A, tT, zero, b); checkDeviceBlasStatus("gemv");
            }
          }
        }
      }
    }

    inline
    void
    solveCholeskyLowerOnDevice(const ordinal_type pbeg, 
                               const ordinal_type pend,
                               const size_type_array_host &h_buf_solve_ptr,
                               const value_type_matrix &t) {
      if (variant == 0) 
        solveCholeskyLowerOnDeviceVar0(pbeg, pend, h_buf_solve_ptr, t);
      else if (variant == 1)
        solveCholeskyLowerOnDeviceVar1(pbeg, pend, h_buf_solve_ptr, t);
      else if (variant == 2)
        solveCholeskyLowerOnDeviceVar2(pbeg, pend, h_buf_solve_ptr, t);
      else {
        TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, 
                                 "LevelSetTools::solveCholeskyLowerOnDevice, algorithm variant is not supported");
      }
    }

    inline
    void
    solveCholeskyUpperOnDeviceVar0(const ordinal_type pbeg,
                                   const ordinal_type pend,
                                   const size_type_array_host &h_buf_solve_ptr,
                                   const value_type_matrix &t) {
      const ordinal_type nrhs = t.extent(1);
      const value_type minus_one(-1), one(1);
#if defined(KOKKOS_ENABLE_CUDA)
      ordinal_type q(0);
#endif 
      exec_space exec_instance;
      for (ordinal_type p=pbeg;p<pend;++p) {
        const ordinal_type sid = _h_level_sids(p);
        if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA)
          const ordinal_type qid = q%_nstreams;
          const auto mystream = _cuda_streams[qid];
          _status = cublasSetStream(_handle_blas, mystream); checkDeviceStatus("cublasSetStream");
          exec_instance = _exec_instances[qid];
          ++q;
#endif          
          const auto &s = _h_supernodes(sid);
          {
            const ordinal_type m = s.m, n = s.n, n_m = n-m;
            if (m > 0) {
              value_type *aptr = s.buf, *bptr = _buf.data()+h_buf_solve_ptr(p-pbeg);; 
              const UnmanagedViewType<value_type_matrix> AL(aptr, m, m); aptr += m*m;
              const UnmanagedViewType<value_type_matrix> bB(bptr, n_m, nrhs); 

              const ordinal_type offm = s.row_begin;
              const auto tT = Kokkos::subview(t, range_type(offm, offm+m), Kokkos::ALL());

              if (n_m > 0) {
                const UnmanagedViewType<value_type_matrix> AR(aptr, m, n_m); // aptr += m*n;
                _status = Gemv<Trans::NoTranspose,Algo::OnDevice>
                  ::invoke(_handle_blas, minus_one, AR, bB, one, tT); checkDeviceBlasStatus("gemv");
                exec_instance.fence();
              }
              _status = Trsv<Uplo::Upper,Trans::NoTranspose,Algo::OnDevice>
                ::invoke(_handle_blas, Diag::NonUnit(), AL, tT); checkDeviceBlasStatus("trsv");
            }
          }
        }
      }
    }

    inline
    void
    solveCholeskyUpperOnDeviceVar1(const ordinal_type pbeg,
                                   const ordinal_type pend,
                                   const size_type_array_host &h_buf_solve_ptr,
                                   const value_type_matrix &t) {
      const ordinal_type nrhs = t.extent(1);
      const value_type minus_one(-1), one(1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA)
      ordinal_type q(0);
#endif 
      exec_space exec_instance;
      for (ordinal_type p=pbeg;p<pend;++p) {
        const ordinal_type sid = _h_level_sids(p);
        if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA)
          const ordinal_type qid = q%_nstreams;
          const auto mystream = _cuda_streams[qid];
          _status = cublasSetStream(_handle_blas, mystream); checkDeviceStatus("cublasSetStream");
          exec_instance = _exec_instances[qid];
          ++q;
#endif          
          const auto &s = _h_supernodes(sid);
          {
            const ordinal_type m = s.m, n = s.n, n_m = n-m;
            if (m > 0) {
              value_type *aptr = s.buf, *bptr = _buf.data()+h_buf_solve_ptr(p-pbeg);; 
              const UnmanagedViewType<value_type_matrix> AL(aptr, m, m); aptr += m*m;
              const UnmanagedViewType<value_type_matrix> bT(bptr, m, nrhs); bptr += m*nrhs; 

              const ordinal_type offm = s.row_begin;
              const auto tT = Kokkos::subview(t, range_type(offm, offm+m), Kokkos::ALL());

              if (n_m > 0) {
                const UnmanagedViewType<value_type_matrix> AR(aptr, m, n_m); // aptr += m*n;
                const UnmanagedViewType<value_type_matrix> bB(bptr, n_m, nrhs); 
                _status = Gemv<Trans::NoTranspose,Algo::OnDevice>
                  ::invoke(_handle_blas, minus_one, AR, bB, one, tT); checkDeviceBlasStatus("gemv");
                exec_instance.fence();
              }

              _status = Gemv<Trans::NoTranspose,Algo::OnDevice>
                ::invoke(_handle_blas, one, AL, tT, zero, bT); checkDeviceBlasStatus("gemv");

              exec_instance.fence();

              _status = Copy<Algo::OnDevice>::invoke(exec_instance, tT, bT); checkDeviceBlasStatus("Copy");
            }
          }
        }
      }
    }

    inline
    void
    solveCholeskyUpperOnDeviceVar2(const ordinal_type pbeg,
                                   const ordinal_type pend,
                                   const size_type_array_host &h_buf_solve_ptr,
                                   const value_type_matrix &t) {
      const ordinal_type nrhs = t.extent(1);
      const value_type one(1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA)
      ordinal_type q(0);
#endif 
      for (ordinal_type p=pbeg;p<pend;++p) {
        const ordinal_type sid = _h_level_sids(p);
        if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA)
          const ordinal_type qid = q%_nstreams;
          const auto mystream = _cuda_streams[qid];
          _status = cublasSetStream(_handle_blas, mystream); checkDeviceStatus("cublasSetStream");
          ++q;
#endif          
          const auto &s = _h_supernodes(sid);
          {
            const ordinal_type m = s.m, n = s.n;
            if (m > 0 && n > 0) {
              value_type *aptr = s.buf, *bptr = _buf.data()+h_buf_solve_ptr(p-pbeg);; 
              const UnmanagedViewType<value_type_matrix> A(aptr, m, n); 
              const UnmanagedViewType<value_type_matrix> b(bptr, n, nrhs); 
               
              const ordinal_type offm = s.row_begin;
              const auto tT = Kokkos::subview(t, range_type(offm, offm+m), Kokkos::ALL());
               
              _status = Gemv<Trans::NoTranspose,Algo::OnDevice>
                ::invoke(_handle_blas, one, A, b, zero, tT); checkDeviceBlasStatus("gemv");
            }
          }
        }
      }
    }

    inline
    void
    solveCholeskyUpperOnDevice(const ordinal_type pbeg,
                               const ordinal_type pend,
                               const size_type_array_host &h_buf_solve_ptr,
                               const value_type_matrix &t) {
      if (variant == 0) 
        solveCholeskyUpperOnDeviceVar0(pbeg, pend, h_buf_solve_ptr, t);
      else if (variant == 1) 
        solveCholeskyUpperOnDeviceVar1(pbeg, pend, h_buf_solve_ptr, t);
      else if (variant == 2) 
        solveCholeskyUpperOnDeviceVar2(pbeg, pend, h_buf_solve_ptr, t);
      else {
        TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, 
                                 "LevelSetTools::solveCholeskyUpperOnDevice, algorithm variant is not supported");
      }
    }

    inline
    void
    solveLDL_LowerOnDevice(const ordinal_type pbeg, 
                           const ordinal_type pend,
                           const size_type_array_host &h_buf_solve_ptr,
                           const value_type_matrix &t) {
      const ordinal_type nrhs = t.extent(1);
      const value_type minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA)
      ordinal_type q(0);
#endif 
      exec_space exec_instance;
      for (ordinal_type p=pbeg;p<pend;++p) {
        const ordinal_type sid = _h_level_sids(p);
        if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA)
          const ordinal_type qid = q%_nstreams;
          const auto mystream = _cuda_streams[qid];
          _status = cublasSetStream(_handle_blas, mystream); checkDeviceStatus("cublasSetStream");
          exec_instance = _exec_instances[qid];
          ++q;
#endif          
          const auto &s = _h_supernodes(sid);
          {
            const ordinal_type m = s.m, n = s.n, n_m = n-m;
            if (m > 0) {
              value_type *aptr = s.buf;
              UnmanagedViewType<value_type_matrix> AL(aptr, m, m); aptr += m*m;

              const ordinal_type offm = s.row_begin;

              const auto tT = Kokkos::subview(t, range_type(offm, offm+m), Kokkos::ALL());
              const auto fpiv = ordinal_type_array(_piv.data()+4*offm+m, m);
              
              _status = ApplyPivots<PivotMode::Flame,Side::Left,Direct::Forward,Algo::OnDevice> /// row inter-change
                ::invoke(exec_instance, fpiv, tT);
              exec_instance.fence();

              _status = Trsv<Uplo::Lower,Trans::NoTranspose,Algo::OnDevice>
                ::invoke(_handle_blas, Diag::Unit(), AL, tT); checkDeviceBlasStatus("trsv");
              exec_instance.fence();
              if (n_m > 0) {
                value_type *bptr = _buf.data()+h_buf_solve_ptr(p-pbeg);
                UnmanagedViewType<value_type_matrix> AR(aptr, m, n_m); // ptr += m*n_m;
                UnmanagedViewType<value_type_matrix> bB(bptr, n_m, nrhs);
                _status = Gemv<Trans::Transpose,Algo::OnDevice>
                  ::invoke(_handle_blas, minus_one, AR, tT, zero, bB); checkDeviceBlasStatus("gemv");
              }
            }
          }
        }
      }
    }


    inline
    void
    solveLDL_UpperOnDevice(const ordinal_type pbeg,
                           const ordinal_type pend,
                           const size_type_array_host &h_buf_solve_ptr,
                           const value_type_matrix &t) {
      const ordinal_type nrhs = t.extent(1);
      const value_type minus_one(-1), one(1);
#if defined(KOKKOS_ENABLE_CUDA)
      ordinal_type q(0);
#endif 
      exec_space exec_instance;
      for (ordinal_type p=pbeg;p<pend;++p) {
        const ordinal_type sid = _h_level_sids(p);
        if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA)
          const ordinal_type qid = q%_nstreams;
          const auto mystream = _cuda_streams[qid];
          _status = cublasSetStream(_handle_blas, mystream); checkDeviceStatus("cublasSetStream");
          exec_instance = _exec_instances[qid];
          ++q;
#endif          
          const auto &s = _h_supernodes(sid);
          {
            const ordinal_type m = s.m, n = s.n, n_m = n-m;
            if (m > 0) {
              value_type *aptr = s.buf, *bptr = _buf.data()+h_buf_solve_ptr(p-pbeg);; 
              const UnmanagedViewType<value_type_matrix> AL(aptr, m, m); aptr += m*m;
              const UnmanagedViewType<value_type_matrix> bB(bptr, n_m, nrhs); 

              const ordinal_type offm = s.row_begin;
              const auto tT = Kokkos::subview(t, range_type(offm, offm+m), Kokkos::ALL());
              const auto P = ordinal_type_array(_piv.data()+4*offm, 4*m);
              const auto D = value_type_matrix(_diag.data()+2*offm, m, 2);
              _status = Scale2x2_BlockInverseDiagonals<Side::Left,Algo::OnDevice> /// row scaling
                ::invoke(exec_instance, P, D, tT);

              if (n_m > 0) {
                const UnmanagedViewType<value_type_matrix> AR(aptr, m, n_m); // aptr += m*n;
                _status = Gemv<Trans::NoTranspose,Algo::OnDevice>
                  ::invoke(_handle_blas, minus_one, AR, bB, one, tT); checkDeviceBlasStatus("gemv");
                exec_instance.fence();
              }
              _status = Trsv<Uplo::Lower,Trans::Transpose,Algo::OnDevice>
                ::invoke(_handle_blas, Diag::Unit(), AL, tT); checkDeviceBlasStatus("trsv");

              const auto fpiv = ordinal_type_array(P.data()+m, m);
              _status = ApplyPivots<PivotMode::Flame,Side::Left,Direct::Backward,Algo::OnDevice> /// row inter-change
                ::invoke(exec_instance, fpiv, tT);
            }
          }
        }
      }
    }


    inline
    void
    allocateWorkspaceSolve(const ordinal_type nrhs) {
      const size_type buf_extent = _bufsize_solve*nrhs;
      const size_type buf_span = _buf.span();
      
      if (buf_extent != buf_span) {
        _buf = value_type_array(do_not_initialize_tag("buf"), buf_extent);
        track_free(buf_span*sizeof(value_type));
        track_alloc(_buf.span()*sizeof(value_type));
        {
          const Kokkos::RangePolicy<exec_space> policy(0,_buf_solve_ptr.extent(0));
          const auto buf_solve_nrhs_ptr = _buf_solve_nrhs_ptr;
          const auto buf_solve_ptr = _buf_solve_ptr;
          Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const ordinal_type &i) {
              buf_solve_nrhs_ptr(i) = nrhs*buf_solve_ptr(i);
            });
        }
        Kokkos::deep_copy(_h_buf_solve_nrhs_ptr, _buf_solve_nrhs_ptr);
      }
    }


    inline
    void
    solveCholesky(const value_type_matrix &x,   // solution
                  const value_type_matrix &b,   // right hand side
                  const value_type_matrix &t,
                  const ordinal_type verbose) { // temporary workspace (store permuted vectors)
      TACHO_TEST_FOR_EXCEPTION(x.extent(0) != b.extent(0) ||
                               x.extent(1) != b.extent(1) ||
                               x.extent(0) != t.extent(0) ||
                               x.extent(1) != t.extent(1), std::logic_error,
                               "x, b, t, and w dimensions do not match");

      TACHO_TEST_FOR_EXCEPTION(x.data() == b.data() ||
                               x.data() == t.data(), std::logic_error,
                               "x, b, t, and w have the same data pointer");
      constexpr bool is_host = std::is_same<exec_memory_space,Kokkos::HostSpace>::value;

      // solve U^{H} (U x) = b 
      const ordinal_type nrhs = x.extent(1);
      Kokkos::Timer timer;

      stat_level.n_kernel_launching = 0;

      // one-time operation when nrhs is changed
      timer.reset();
      allocateWorkspaceSolve(nrhs);

      // 0. permute and copy b -> t
      applyRowPermutationToDenseMatrix(t, b, _perm);      
      stat.t_extra = timer.seconds();

      timer.reset();
      { 
#if defined(TACHO_ENABLE_SOLVE_CHOLESKY_USE_LIGHT_KERNEL)
        const auto work_item_property = Kokkos::Experimental::WorkItemProperty::HintLightWeight;
#endif
        // this should be considered with average problem sizes in levels
        const ordinal_type half_level = _nlevel/2;
        const ordinal_type team_size_solve[2] = { 64, 16 }, vector_size_solve[2] = { 8, 8};
        const ordinal_type team_size_update[2] = { 128, 32}, vector_size_update[2] = { 1, 1};
        {
          typedef TeamFunctor_SolveLowerChol<supernode_info_type> functor_type;
#if defined(TACHO_TEST_SOLVE_CHOLESKY_KERNEL_OVERHEAD)
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::DummyTag> team_policy_solve;
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::DummyTag> team_policy_update;
#else
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::template SolveTag<variant> > team_policy_solve;
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::template UpdateTag<variant> > team_policy_update;
#endif          
          functor_type functor(_info, 
                               _solve_mode,
                               _level_sids,
                               t,
                               _buf);

          team_policy_solve policy_solve(1,1,1);
          team_policy_update policy_update(1,1,1);

          //  1. U^{H} w = t
          {
            for (ordinal_type lvl=(_team_serial_level_cut-1);lvl>=0;--lvl) {
              const ordinal_type 
                pbeg = _h_level_ptr(lvl), 
                pend = _h_level_ptr(lvl+1),
                pcnt = pend - pbeg;

              const range_type range_solve_buf_ptr(_h_buf_level_ptr(lvl), _h_buf_level_ptr(lvl+1));

              const auto solve_buf_ptr = Kokkos::subview(_buf_solve_nrhs_ptr, range_solve_buf_ptr);
              functor.setRange(pbeg, pend);
              functor.setBufferPtr(solve_buf_ptr);
              if (is_host) {
                policy_solve  = team_policy_solve(pcnt, 1, 1);
                policy_update = team_policy_update(pcnt, 1, 1);
              } else {
                const ordinal_type idx = lvl > half_level;
                policy_solve  = team_policy_solve(pcnt, team_size_solve[idx],  vector_size_solve[idx]);
                policy_update = team_policy_update(pcnt, team_size_update[idx], vector_size_update[idx]);
              }
#if defined(TACHO_ENABLE_SOLVE_CHOLESKY_USE_LIGHT_KERNEL)
              const auto policy_solve_with_work_property = Kokkos::Experimental::require(policy_solve, work_item_property);
              const auto policy_update_with_work_property = Kokkos::Experimental::require(policy_update, work_item_property);
#else
              const auto policy_solve_with_work_property = policy_solve;
              const auto policy_update_with_work_property = policy_update;
#endif
              if (lvl < _device_level_cut) {
                // do nothing
                //Kokkos::parallel_for("solve lower", policy_solve, functor);
              } else {
                Kokkos::parallel_for("solve lower", 
                                     policy_solve_with_work_property, 
                                     functor);
                ++stat_level.n_kernel_launching;
              }
              const auto h_buf_solve_ptr = Kokkos::subview(_h_buf_solve_nrhs_ptr, range_solve_buf_ptr);              
              solveCholeskyLowerOnDevice(pbeg, pend, h_buf_solve_ptr, t); 
              Kokkos::fence();
              
              Kokkos::parallel_for("update lower", 
                                   policy_update_with_work_property, 
                                   functor); 
              ++stat_level.n_kernel_launching;
              exec_space().fence(); //Kokkos::fence();
            }
          }
        } // end of lower tri solve
        
        {
          typedef TeamFunctor_SolveUpperChol<supernode_info_type> functor_type;
#if defined(TACHO_TEST_SOLVE_CHOLESKY_KERNEL_OVERHEAD)
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::DummyTag> team_policy_solve;
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::DummyTag> team_policy_update;
#else
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::template SolveTag<variant> > team_policy_solve;
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::template UpdateTag<variant> > team_policy_update;
#endif
          functor_type functor(_info, 
                               _solve_mode,
                               _level_sids,
                               t,
                               _buf);
          
          team_policy_solve policy_solve(1,1,1);
          team_policy_update policy_update(1,1,1);
          
          //  2. U t = w;
          {
            for (ordinal_type lvl=0;lvl<_team_serial_level_cut;++lvl) {
              const ordinal_type 
                pbeg = _h_level_ptr(lvl), 
                pend = _h_level_ptr(lvl+1),
                pcnt = pend - pbeg;
              
              const range_type range_solve_buf_ptr(_h_buf_level_ptr(lvl), _h_buf_level_ptr(lvl+1));
              const auto solve_buf_ptr = Kokkos::subview(_buf_solve_nrhs_ptr, range_solve_buf_ptr);
              functor.setRange(pbeg, pend);
              functor.setBufferPtr(solve_buf_ptr);
              if (is_host) {
                policy_solve  = team_policy_solve(pcnt, 1, 1);
                policy_update = team_policy_update(pcnt, 1, 1);
              } else {
                const ordinal_type idx = lvl > half_level;
                policy_solve  = team_policy_solve(pcnt, team_size_solve[idx],  vector_size_solve[idx]);
                policy_update = team_policy_update(pcnt, team_size_update[idx], vector_size_update[idx]);
              }
#if defined(TACHO_ENABLE_SOLVE_CHOLESKY_USE_LIGHT_KERNEL)
              const auto policy_solve_with_work_property = Kokkos::Experimental::require(policy_solve, work_item_property);
              const auto policy_update_with_work_property = Kokkos::Experimental::require(policy_update, work_item_property);
#else
              const auto policy_solve_with_work_property = policy_solve;
              const auto policy_update_with_work_property = policy_update;
#endif
              Kokkos::parallel_for("update upper", 
                                   policy_update_with_work_property,
                                   functor);
              ++stat_level.n_kernel_launching;
              exec_space().fence(); //Kokkos::fence();

              if (lvl < _device_level_cut) {
                // do nothing
                //Kokkos::parallel_for("solve upper", policy_solve, functor); 
              } else {
                Kokkos::parallel_for("solve upper", 
                                     policy_solve_with_work_property,
                                     functor);
                ++stat_level.n_kernel_launching;
              }

              const auto h_buf_solve_ptr = Kokkos::subview(_h_buf_solve_nrhs_ptr, range_solve_buf_ptr);
              solveCholeskyUpperOnDevice(pbeg, pend, h_buf_solve_ptr, t);
              Kokkos::fence();
            }
          }
        }/// end of upper tri solve

      } // end of solve
      stat.t_solve = timer.seconds();

      // permute and copy t -> x
      timer.reset();
      applyRowPermutationToDenseMatrix(x, t, _peri);
      stat.t_extra += timer.seconds();

      if (verbose) {
        printf("Summary: LevelSetTools-Variant-%d (Cholesky Solve: %3d)\n", variant, nrhs);
        printf("=======================================================\n");
        print_stat_solve();
      }
    }

    inline
    void
    factorizeLDL(const value_type_array &ax,
                 const ordinal_type verbose) {
      constexpr bool is_host = std::is_same<exec_memory_space,Kokkos::HostSpace>::value;
      Kokkos::Timer timer;
      
      timer.reset();
      value_type_array work;
      {
        _buf = value_type_array(do_not_initialize_tag("buf"), _bufsize_factorize);
        track_alloc(_buf.span()*sizeof(value_type));
        
#if defined (KOKKOS_ENABLE_CUDA)
        value_type_matrix T(NULL, _info.max_supernode_size, _info.max_supernode_size);
        ordinal_type_array P(NULL, _info.max_supernode_size);
        const size_type worksize = LDL<Uplo::Lower,Algo::OnDevice>
          ::invoke(_handle_lapack, T, P, work); 
        work = value_type_array(do_not_initialize_tag("work"), worksize*(_nstreams+1)*max(8, _nstreams));
#else
        const size_type worksize = 32*_info.max_supernode_size;
        work = value_type_array(do_not_initialize_tag("work"), worksize);
#endif
        track_alloc(work.span()*sizeof(value_type));
      }
      stat.t_extra = timer.seconds();
      
      timer.reset();
      {
        _ax = ax; // matrix values
        _info.copySparseToSuperpanels(_ap, _aj, _ax, _perm, _peri);
      }
      stat.t_copy = timer.seconds();

      stat_level.n_kernel_launching = 0;
      timer.reset();
      { 
        // this should be considered with average problem sizes in levels
        const ordinal_type half_level = _nlevel/2;
        //const ordinal_type team_size_factor[2] = { 64, 16 }, vector_size_factor[2] = { 8, 8};
        //const ordinal_type team_size_factor[2] = { 16, 16 }, vector_size_factor[2] = { 32, 32};
#if defined (CUDA_VERSION)
#if (11000 > CUDA_VERSION)
        /// cuda 11.1 below
        const ordinal_type team_size_factor[2] = { 32, 64 }, vector_size_factor[2] = { 8, 4};        
#else 
        /// cuda 11.1 and higher
        const ordinal_type team_size_factor[2] = { 64, 64 }, vector_size_factor[2] = { 8, 4};
#endif
#else
        /// not cuda ... whatever..
        const ordinal_type team_size_factor[2] = { 64, 64 }, vector_size_factor[2] = { 8, 4};
#endif
        const ordinal_type team_size_update[2] = { 16, 8 }, vector_size_update[2] = { 32, 32};
        {
          typedef TeamFunctor_FactorizeLDL<supernode_info_type> functor_type;
#if defined(TACHO_TEST_LEVELSET_TOOLS_KERNEL_OVERHEAD)
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::DummyTag> team_policy_factorize;
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::DummyTag> team_policy_update;
#else
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::FactorizeTag> team_policy_factor;
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::UpdateTag> team_policy_update;
#endif
          functor_type functor(_info, 
                               _factorize_mode,
                               _level_sids,
                               _piv,
                               _diag,
                               _buf);

          team_policy_factor policy_factor(1,1,1);
          team_policy_update policy_update(1,1,1);

          {
            for (ordinal_type lvl=(_team_serial_level_cut-1);lvl>=0;--lvl) {
              const ordinal_type 
                pbeg = _h_level_ptr(lvl), 
                pend = _h_level_ptr(lvl+1),
                pcnt = pend - pbeg;

              const range_type range_buf_factor_ptr(_h_buf_level_ptr(lvl), _h_buf_level_ptr(lvl+1));

              const auto buf_factor_ptr = Kokkos::subview(_buf_factor_ptr, range_buf_factor_ptr);
              functor.setRange(pbeg, pend);
              functor.setBufferPtr(buf_factor_ptr);
              if (is_host) {
                policy_factor = team_policy_factor(pcnt, 1, 1);
                policy_update = team_policy_update(pcnt, 1, 1);
              } else {
                const ordinal_type idx = lvl > half_level;
                policy_factor = team_policy_factor(pcnt, team_size_factor[idx], vector_size_factor[idx]);
                policy_update = team_policy_update(pcnt, team_size_update[idx], vector_size_update[idx]);
              }
              if (lvl < _device_level_cut) {
                // do nothing
                //Kokkos::parallel_for("factor lower", policy_factor, functor);
              } else {
                Kokkos::parallel_for("factor", policy_factor, functor);
                ++stat_level.n_kernel_launching;
              }

              const auto h_buf_factor_ptr = Kokkos::subview(_h_buf_factor_ptr, range_buf_factor_ptr);

              factorizeLDL_OnDevice(pbeg, pend, h_buf_factor_ptr, work); 
              Kokkos::fence();

              Kokkos::parallel_for("update factor", policy_update, functor); 
              ++stat_level.n_kernel_launching;
              exec_space().fence(); //Kokkos::fence();
            }
          }
        }
      } // end of LDL
      stat.t_factor = timer.seconds();

      timer.reset();
      {
#if defined (KOKKOS_ENABLE_CUDA)
        track_free(work.span()*sizeof(value_type));
#endif
        track_free(_buf.span()*sizeof(value_type));
        _buf = value_type_array();
      }
      stat.t_extra += timer.seconds();

      if (verbose) {
        printf("Summary: LevelSetTools (LDL Factorize)\n");
        printf("======================================\n");
        print_stat_factor();
      }
    }

    inline
    void
    solveLDL(const value_type_matrix &x,   // solution
             const value_type_matrix &b,   // right hand side
             const value_type_matrix &t,   // temporary workspace (store permuted vectors)
             const ordinal_type verbose) {
      TACHO_TEST_FOR_EXCEPTION(x.extent(0) != b.extent(0) ||
                               x.extent(1) != b.extent(1) ||
                               x.extent(0) != t.extent(0) ||
                               x.extent(1) != t.extent(1), std::logic_error,
                               "x, b, t, and w dimensions do not match");

      TACHO_TEST_FOR_EXCEPTION(x.data() == b.data() ||
                               x.data() == t.data(), std::logic_error,
                               "x, b, t, and w have the same data pointer");
      constexpr bool is_host = std::is_same<exec_memory_space,Kokkos::HostSpace>::value;

      // solve L D L^{H} x = b 
      const ordinal_type nrhs = x.extent(1);
      Kokkos::Timer timer;

      stat_level.n_kernel_launching = 0;

      // one-time operation when nrhs is changed
      timer.reset();
      allocateWorkspaceSolve(nrhs);

      // 0. permute and copy b -> t
      applyRowPermutationToDenseMatrix(t, b, _perm);      
      stat.t_extra = timer.seconds();

      timer.reset();
      { 
#if defined(TACHO_ENABLE_SOLVE_CHOLESKY_USE_LIGHT_KERNEL)
        const auto work_item_property = Kokkos::Experimental::WorkItemProperty::HintLightWeight;
#endif
        // this should be considered with average problem sizes in levels
        const ordinal_type half_level = _nlevel/2;
#if defined (CUDA_VERSION)
#if (11000 > CUDA_VERSION)
        /// cuda 11.1 below
        const ordinal_type team_size_solve[2] = { 32, 16 }, vector_size_solve[2] = { 8, 8};
#else
        /// cuda 11.1 and higher
        const ordinal_type team_size_solve[2] = { 32, 16 }, vector_size_solve[2] = { 8, 8};
#endif
#else
        /// not cuda whatever...
        const ordinal_type team_size_solve[2] = { 64, 16 }, vector_size_solve[2] = { 8, 8};
#endif
        const ordinal_type team_size_update[2] = { 128, 32}, vector_size_update[2] = { 1, 1};
        {
          typedef TeamFunctor_SolveLowerLDL<supernode_info_type> functor_type;
#if defined(TACHO_TEST_SOLVE_CHOLESKY_KERNEL_OVERHEAD)
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::DummyTag> team_policy_solve;
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::DummyTag> team_policy_update;
#else
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::SolveTag> team_policy_solve;
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::UpdateTag> team_policy_update;
#endif          
          functor_type functor(_info, 
                               _solve_mode,
                               _level_sids,
                               _piv,
                               t,
                               _buf);

          team_policy_solve policy_solve(1,1,1);
          team_policy_update policy_update(1,1,1);

          //  1. L w = t
          {
            for (ordinal_type lvl=(_team_serial_level_cut-1);lvl>=0;--lvl) {
              const ordinal_type 
                pbeg = _h_level_ptr(lvl), 
                pend = _h_level_ptr(lvl+1),
                pcnt = pend - pbeg;

              const range_type range_solve_buf_ptr(_h_buf_level_ptr(lvl), _h_buf_level_ptr(lvl+1));

              const auto solve_buf_ptr = Kokkos::subview(_buf_solve_nrhs_ptr, range_solve_buf_ptr);
              functor.setRange(pbeg, pend);
              functor.setBufferPtr(solve_buf_ptr);
              if (is_host) {
                policy_solve  = team_policy_solve(pcnt, 1, 1);
                policy_update = team_policy_update(pcnt, 1, 1);
              } else {
                const ordinal_type idx = lvl > half_level;
                policy_solve  = team_policy_solve(pcnt, team_size_solve[idx],  vector_size_solve[idx]);
                policy_update = team_policy_update(pcnt, team_size_update[idx], vector_size_update[idx]);
              }
#if defined(TACHO_ENABLE_SOLVE_CHOLESKY_USE_LIGHT_KERNEL)
              const auto policy_solve_with_work_property = Kokkos::Experimental::require(policy_solve, work_item_property);
              const auto policy_update_with_work_property = Kokkos::Experimental::require(policy_update, work_item_property);
#else
              const auto policy_solve_with_work_property = policy_solve;
              const auto policy_update_with_work_property = policy_update;
#endif
              if (lvl < _device_level_cut) {
                // do nothing
                //Kokkos::parallel_for("solve lower", policy_solve, functor);
              } else {
                Kokkos::parallel_for("solve lower", 
                                     policy_solve_with_work_property, 
                                     functor);
                ++stat_level.n_kernel_launching;
              }
              const auto h_buf_solve_ptr = Kokkos::subview(_h_buf_solve_nrhs_ptr, range_solve_buf_ptr);              
              solveLDL_LowerOnDevice(pbeg, pend, h_buf_solve_ptr, t); 
              Kokkos::fence();
              
              Kokkos::parallel_for("update lower", 
                                   policy_update_with_work_property, 
                                   functor); 
              ++stat_level.n_kernel_launching;
              exec_space().fence(); //Kokkos::fence();
            }
          }
        } // end of lower tri solve
        
        {
          typedef TeamFunctor_SolveUpperLDL<supernode_info_type> functor_type;
#if defined(TACHO_TEST_SOLVE_CHOLESKY_KERNEL_OVERHEAD)
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::DummyTag> team_policy_solve;
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::DummyTag> team_policy_update;
#else
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::SolveTag> team_policy_solve;
          typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>,exec_space,
                                     typename functor_type::UpdateTag> team_policy_update;
#endif
          functor_type functor(_info, 
                               _solve_mode,
                               _level_sids,
                               _piv,
                               _diag,
                               t,
                               _buf);
          
          team_policy_solve policy_solve(1,1,1);
          team_policy_update policy_update(1,1,1);
          
          //  2. U t = w;
          {
            for (ordinal_type lvl=0;lvl<_team_serial_level_cut;++lvl) {
              const ordinal_type 
                pbeg = _h_level_ptr(lvl), 
                pend = _h_level_ptr(lvl+1),
                pcnt = pend - pbeg;
              
              const range_type range_solve_buf_ptr(_h_buf_level_ptr(lvl), _h_buf_level_ptr(lvl+1));
              const auto solve_buf_ptr = Kokkos::subview(_buf_solve_nrhs_ptr, range_solve_buf_ptr);
              functor.setRange(pbeg, pend);
              functor.setBufferPtr(solve_buf_ptr);
              if (is_host) {
                policy_solve  = team_policy_solve(pcnt, 1, 1);
                policy_update = team_policy_update(pcnt, 1, 1);
              } else {
                const ordinal_type idx = lvl > half_level;
                policy_solve  = team_policy_solve(pcnt, team_size_solve[idx],  vector_size_solve[idx]);
                policy_update = team_policy_update(pcnt, team_size_update[idx], vector_size_update[idx]);
              }
#if defined(TACHO_ENABLE_SOLVE_CHOLESKY_USE_LIGHT_KERNEL)
              const auto policy_solve_with_work_property = Kokkos::Experimental::require(policy_solve, work_item_property);
              const auto policy_update_with_work_property = Kokkos::Experimental::require(policy_update, work_item_property);
#else
              const auto policy_solve_with_work_property = policy_solve;
              const auto policy_update_with_work_property = policy_update;
#endif
              Kokkos::parallel_for("update upper", 
                                   policy_update_with_work_property,
                                   functor);
              ++stat_level.n_kernel_launching;
              exec_space().fence(); //Kokkos::fence();

              if (lvl < _device_level_cut) {
                // do nothing
                //Kokkos::parallel_for("solve upper", policy_solve, functor); 
              } else {
                Kokkos::parallel_for("solve upper", 
                                     policy_solve_with_work_property,
                                     functor);
                ++stat_level.n_kernel_launching;
              }

              const auto h_buf_solve_ptr = Kokkos::subview(_h_buf_solve_nrhs_ptr, range_solve_buf_ptr);
              solveLDL_UpperOnDevice(pbeg, pend, h_buf_solve_ptr, t);
              Kokkos::fence();
            }
          }
        }/// end of upper tri solve

      } // end of solve
      stat.t_solve = timer.seconds();

      // permute and copy t -> x
      timer.reset();
      applyRowPermutationToDenseMatrix(x, t, _peri);
      stat.t_extra += timer.seconds();

      if (verbose) {
        printf("Summary: LevelSetTools (LDL Solve: %3d)\n", nrhs);
        printf("=======================================\n");
        print_stat_solve();
      }
    }

    inline
    void
    factorize(const value_type_array &ax,
              const ordinal_type verbose = 0) override {
      Kokkos::deep_copy(_superpanel_buf, value_type(0));
      switch (this->getSolutionMethod()) {
      case 1: { /// Cholesky
        factorizeCholesky(ax, verbose);
        break;
      }
      case 2: { /// LDL
        {
          const ordinal_type rlen = 4*_m, plen = _piv.span();
          if (plen < rlen) {
            track_free(_piv.span()*sizeof(ordinal_type));
            _piv = ordinal_type_array("piv", rlen);
            track_alloc(_piv.span()*sizeof(ordinal_type));
          }
        }
        {
          const ordinal_type rlen = 2*_m, dlen = _diag.span();
          if (dlen < rlen) {        
            track_free(_diag.span()*sizeof(value_type));
            _diag = value_type_array("diag", rlen);
            track_alloc(_diag.span()*sizeof(value_type));
          }
        }
        factorizeLDL(ax, verbose);
        break;
      }
      default: {
        TACHO_TEST_FOR_EXCEPTION(false,
                                 std::logic_error,
                                 "The solution method is not supported");
        break;
      }
      }
    }

    inline
    void
    solve(const value_type_matrix &x,   // solution
          const value_type_matrix &b,   // right hand side
          const value_type_matrix &t,   // temporary workspace (store permuted vectors)
          const ordinal_type verbose = 0) override {
      switch (this->getSolutionMethod()) {
      case 1: { /// Cholesky
        solveCholesky(x, b, t, verbose);
        break;
      }
      case 2: { /// LDL
        solveLDL(x, b, t, verbose);
        break;
      }
      default: {
        TACHO_TEST_FOR_EXCEPTION(false,
                                 std::logic_error,
                                 "The solution method is not supported");
        break;
      }
      }
    }
    
  };

}
#endif

