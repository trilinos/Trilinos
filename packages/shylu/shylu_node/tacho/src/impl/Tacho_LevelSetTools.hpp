#ifndef __TACHO_LEVEL_SET_TOOLS_HPP__
#define __TACHO_LEVEL_SET_TOOLS_HPP__

/// \file Tacho_LevelSetTools.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"
#include "Tacho_DenseFlopCount.hpp"

#include "Tacho_Copy.hpp"
#include "Tacho_Copy_OnDevice.hpp"

#include "Tacho_SetIdentity.hpp"
#include "Tacho_SetIdentity_OnDevice.hpp"

#include "Tacho_Chol.hpp"
#include "Tacho_Chol_OnDevice.hpp"

#include "Tacho_Trsm.hpp"
#include "Tacho_Trsm_OnDevice.hpp"

#include "Tacho_Herk.hpp"
#include "Tacho_Herk_OnDevice.hpp"

#include "Tacho_Trsv.hpp"
#include "Tacho_Trsv_OnDevice.hpp"

#include "Tacho_Gemv.hpp"
#include "Tacho_Gemv_OnDevice.hpp"

#include "Tacho_SupernodeInfo.hpp"

#include "Tacho_TeamFunctor_FactorizeChol.hpp"
#include "Tacho_TeamFunctor_SolveLowerChol.hpp"
#include "Tacho_TeamFunctor_SolveUpperChol.hpp"

//#define TACHO_TEST_LEVELSET_TOOLS_KERNEL_OVERHEAD
//#define TACHO_ENABLE_LEVELSET_TOOLS_USE_LIGHT_KERNEL

namespace Tacho {

  ///
  /// Here we do not use a scheduler but all derived types in supernodes 
  /// info are required scheduler  
  ///
  template<typename ValueType, typename SchedulerType, int Var>
  class LevelSetTools {
  public:
    enum { variant = Var,
           max_factor_team_size = 64 };

    typedef ValueType value_type;
    typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;

    ///
    /// device space typedefs
    ///
    typedef SchedulerType scheduler_type;

    typedef typename UseThisDevice<typename scheduler_type::execution_space>::device_type device_type;
    typedef typename device_type::execution_space exec_space;
    typedef typename device_type::memory_space exec_memory_space;

    typedef SupernodeInfo<value_type,scheduler_type> supernode_info_type;

    typedef typename supernode_info_type::ordinal_type_array ordinal_type_array;
    typedef typename supernode_info_type::size_type_array size_type_array;
    typedef typename supernode_info_type::value_type_array value_type_array;

    typedef typename supernode_info_type::value_type_matrix value_type_matrix;
    typedef typename supernode_info_type::supernode_type_array supernode_type_array;
    typedef typename supernode_type_array::HostMirror supernode_type_array_host;

    typedef typename UseThisDevice<Kokkos::DefaultHostExecutionSpace>::device_type host_device_type;
    typedef typename host_device_type::execution_space host_space;
    typedef typename host_device_type::memory_space host_memory_space;

    typedef Kokkos::View<ordinal_type*,host_device_type> ordinal_type_array_host;
    typedef Kokkos::View<size_type*,host_device_type> size_type_array_host;

  private:

    ///
    /// input matrix graph in crs
    ///
    ordinal_type _m;
    size_type_array _ap;
    ordinal_type_array _aj;
    value_type_array _ax;

    // graph ordering input
    ordinal_type_array _perm, _peri;

    ///
    /// supernode info: supernode data structure with "unamanged" view
    /// this is passed into computation algorithm without reference counting
    /// supernode data is managed in the numeric tools
    ///
    supernode_info_type _info;

    // supernode host information for level kernels launching
    ordinal_type _nsupernodes;
    supernode_type_array_host _h_supernodes;
    ConstUnmanagedViewType<ordinal_type_array_host> _h_stree_level;

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
#if defined(KOKKOS_ENABLE_CUDA)
    int _nstreams;
    bool _is_cublas_created, _is_cusolver_dn_created;
    cublasHandle_t _handle_blas;
    cusolverDnHandle_t _handle_lapack;
    using cuda_stream_array_host = std::vector<cudaStream_t>;
    cuda_stream_array_host _cuda_streams;

    using exec_instance_array_host = std::vector<exec_space>;
    exec_instance_array_host _exec_instances;
#else 
    int _nstreams, _handle_blas, _handle_lapack; // dummy handle for convenience
#endif

    ///
    /// statistics
    ///
    struct {
      double t_init, t_mode_classification, t_copy, t_factor, t_solve, t_extra;
      double m_used, m_peak;
      int n_device_factorize, n_team_factorize, n_kernel_launching_factorize;
      int n_device_solve, n_team_solve, n_kernel_launching_solve;
      int n_kernel_launching;
    } stat;

  public:

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

    ///
    /// memory tracking and print stat
    ///

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
      stat.t_solve = 0;
      stat.t_extra = 0;
      stat.m_used = 0;
      stat.m_peak = 0;
    }

    inline
    void
    print_stat_init() {
      printf("  Time\n");
      printf("             time for initialization:                         %10.6f s\n", stat.t_init);
      printf("             time for compute mode classification:            %10.6f s\n", stat.t_mode_classification);
      printf("             total time spent:                                %10.6f s\n", (stat.t_init+stat.t_mode_classification));
      printf("\n");
      printf("  Memory\n");
      printf("             workspace allocated for solve:                   %10.2f MB\n", stat.m_used/1024/1024);
      printf("             peak memory used:                                %10.2f MB\n", stat.m_peak/1024/1024);
      printf("\n");
      printf("  Compute Mode in Factorize with a Threshold(%d)\n", _device_factorize_thres);
      printf("             # of subproblems using device functions:         %6d\n", stat.n_device_factorize);
      printf("             # of subproblems using team functions:           %6d\n", stat.n_team_factorize);
      printf("             total # of subproblems:                          %6d\n", (stat.n_device_factorize+stat.n_team_factorize));
      printf("\n");
      printf("  Compute Mode in Solve with a Threshold(%d)\n", _device_solve_thres);
      printf("             # of subproblems using device functions:         %6d\n", stat.n_device_solve);
      printf("             # of subproblems using team functions:           %6d\n", stat.n_team_solve);
      printf("             total # of subproblems:                          %6d\n", (stat.n_device_solve+stat.n_team_solve));
      printf("\n");
    }

    inline
    void
    print_stat_factorize() {
      double flop = 0;
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

      printf("  Time\n");
      printf("             time for extra work e.g.,workspace allocation:   %10.6f s\n", stat.t_extra);
      printf("             time for copying A into U:                       %10.6f s\n", stat.t_copy);
      printf("             time for factorization:                          %10.6f s\n", stat.t_factor);
      printf("             total time spent:                                %10.6f s\n", (stat.t_extra+stat.t_copy+stat.t_factor));
      printf("\n");
      printf("  Memory\n");
      printf("             memory used in factorization:                    %10.2f MB\n", stat.m_used/1024/1024);
      printf("             peak memory used in factorization:               %10.2f MB\n", stat.m_peak/1024/1024);
      printf("\n");
      printf("  Kernels\n");
      printf("             # of kernels launching:                          %6d\n", stat.n_kernel_launching);
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
      printf("             time for extra work e.g.,workspace and permute:  %10.6f s\n", stat.t_extra);
      printf("             time for solve:                                  %10.6f s\n", stat.t_solve);
      printf("             total time spent:                                %10.6f s\n", (stat.t_solve+stat.t_extra));
      printf("\n");
      printf("  Memory\n");
      printf("             memory used in solve:                            %10.2f MB\n", stat.m_used/1024/1024);
      printf("\n");
      printf("  Kernels\n");
      printf("             # of kernels launching:                          %6d\n", stat.n_kernel_launching);
      printf("\n");
    }

    inline
    void
    print_stat_memory() {
      printf("  Memory\n");
      printf("             memory used now:                                 %10.2f MB\n", stat.m_used/1024/1024);
      printf("             peak memory used:                                %10.2f MB\n", stat.m_peak/1024/1024);
      printf("\n");
    }

    ///
    /// initialization / release
    ///
    inline
    void
    initialize(const ordinal_type device_level_cut,
               const ordinal_type device_factorize_thres,
               const ordinal_type device_solve_thres,
               const ordinal_type verbose = 0) {
      stat.n_device_factorize = 0;   stat.n_device_solve = 0;
      stat.n_team_factorize= 0;      stat.n_team_solve = 0;

      Kokkos::Impl::Timer timer;

      timer.reset();

      ///
      /// level data structure
      ///

      // # of supernodes
      _nsupernodes = _info.supernodes.extent(0);

      // local host supernodes info 
      _h_supernodes = Kokkos::create_mirror_view(host_memory_space(), _info.supernodes);
      Kokkos::deep_copy(_h_supernodes, _info.supernodes);

      // # of levels
      _nlevel = 0;
      {
        for (ordinal_type sid=0;sid<_nsupernodes;++sid) 
          _nlevel = max(_h_stree_level(sid), _nlevel);
        ++_nlevel;
      }

      // create level ptr
      _h_level_ptr = size_type_array_host("h_level_ptr", _nlevel+1);
      {
        // first count # of supernodes in each level
        for (ordinal_type sid=0;sid<_nsupernodes;++sid) 
          ++_h_level_ptr(_h_stree_level(sid)+1);

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
          const ordinal_type lvl = _h_stree_level(sid);
          _h_level_sids(tmp_level_ptr(lvl)++) = sid;
        }
      }
      _level_sids = Kokkos::create_mirror_view(exec_memory_space(), _h_level_sids); 
      Kokkos::deep_copy(_level_sids, _h_level_sids);
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
            const ordinal_type factor_work_size_variants[3] = { schur_work_size, std::max(m*m, schur_work_size), m*m+schur_work_size };
            const ordinal_type factor_work_size = factor_work_size_variants[variant];
            const ordinal_type solve_work_size = variant == 0 ? n_m : n;
            _h_buf_factor_ptr(k) = factor_work_size + _h_buf_factor_ptr(k-1);
            _h_buf_solve_ptr(k) = solve_work_size + _h_buf_solve_ptr(k-1);
          }
          const ordinal_type last_idx = lbeg+pend-pbeg;
          _bufsize_factorize = max(_bufsize_factorize, _h_buf_factor_ptr(last_idx));
          _bufsize_solve = max(_bufsize_solve, _h_buf_solve_ptr(last_idx));
        }
      }

      _buf_factor_ptr = Kokkos::create_mirror_view(exec_memory_space(), _h_buf_factor_ptr);
      Kokkos::deep_copy(_buf_factor_ptr, _h_buf_factor_ptr);
      track_alloc(_buf_factor_ptr.span()*sizeof(size_type));

      _buf_solve_ptr = Kokkos::create_mirror_view(exec_memory_space(), _h_buf_solve_ptr);
      Kokkos::deep_copy(_buf_solve_ptr, _h_buf_solve_ptr);
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
            ++stat.n_device_solve;
            ++stat.n_device_factorize;
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
              ++stat.n_device_solve;
            } else {
              _h_solve_mode(sid) = 1;
              ++stat.n_team_solve;
            }
            if (m > _device_factorize_thres) {// || n_m > _device_factorize_thres) {
              _h_factorize_mode(sid) = 0;
              ++stat.n_device_factorize;
            } else {
              _h_factorize_mode(sid) = 1;
              ++stat.n_team_factorize;
            }
          }
        }
      }

      _factorize_mode = Kokkos::create_mirror_view(exec_memory_space(), _h_factorize_mode);
      Kokkos::deep_copy(_factorize_mode, _h_factorize_mode);
      track_alloc(_factorize_mode.span()*sizeof(ordinal_type));

      _solve_mode = Kokkos::create_mirror_view(exec_memory_space(), _h_solve_mode);
      Kokkos::deep_copy(_solve_mode, _h_solve_mode);
      track_alloc(_solve_mode.span()*sizeof(ordinal_type));

      stat.t_mode_classification = timer.seconds();

      if (verbose) {
        printf("Summary: LevelSetTools-Variant-%d (Initialize)\n", variant);
        printf("===============================================\n");
        print_stat_init();
      }
    }

    inline
    void
    release(const ordinal_type verbose = 0) {
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


    LevelSetTools()
      : _m(0), 
        _nsupernodes(0), 
        _nlevel(0), 
        _bufsize_factorize(0), 
        _bufsize_solve(0), 
        _nstreams(0), 
        stat() {}

    LevelSetTools(const LevelSetTools &b) = default;

    LevelSetTools(const NumericTools<value_type,scheduler_type> &N)
      :
      _m(N.getNumRows()),
      _ap(N.getRowPtr()),
      _aj(N.getCols()),
      _perm(N.getPermutationVector()),
      _peri(N.getInversePermutationVector()),
      _info(N.getSupernodesInfo()),
      _h_stree_level(N.getSupernodesTreeLevel()),
      _nstreams(0)
    {
#if defined(KOKKOS_ENABLE_CUDA)
      _is_cublas_created = 0;
      _is_cusolver_dn_created = 0;
#endif
    }

    virtual~LevelSetTools() {
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
    createStream(const ordinal_type nstreams) {
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
#endif
    }

    ///
    /// Device level functions
    ///
    inline
    void
    factorizeOnDeviceVar0(const ordinal_type pbeg, 
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
    factorizeOnDeviceVar1(const ordinal_type pbeg, 
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
    factorizeOnDeviceVar2(const ordinal_type pbeg, 
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
    factorizeOnDevice(const ordinal_type pbeg, 
                      const ordinal_type pend,
                      const size_type_array_host &h_buf_factor_ptr,
                      const value_type_array &work) {
      if      (variant == 0) 
        factorizeOnDeviceVar0(pbeg, pend, h_buf_factor_ptr, work); 
      else if (variant == 1)
        factorizeOnDeviceVar1(pbeg, pend, h_buf_factor_ptr, work); 
      else if (variant == 2)
        factorizeOnDeviceVar2(pbeg, pend, h_buf_factor_ptr, work); 
      else {
        TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, 
                                 "LevelSetTools::factorizeOnDevice, algorithm variant is not supported");
      }
    }

    ///
    /// Level set factorize
    ///
    inline 
    void
    factorizeCholesky(const value_type_array &ax,
                      const ordinal_type verbose = 0) {
      constexpr bool is_host = std::is_same<exec_memory_space,Kokkos::HostSpace>::value;
      Kokkos::Impl::Timer timer;

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

      stat.n_kernel_launching = 0;
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
                ++stat.n_kernel_launching;
              }

              const auto h_buf_factor_ptr = Kokkos::subview(_h_buf_factor_ptr, range_buf_factor_ptr);
              factorizeOnDevice(pbeg, pend, h_buf_factor_ptr, work); 
              Kokkos::fence();

              Kokkos::parallel_for("update factor", policy_update, functor); 
              ++stat.n_kernel_launching;
              exec_space().fence(); //Kokkos::fence();
            }
          }
        } // end of lower tri solve
      } // end of solve
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
        printf("Summary: LevelSetTools-Variant-%d (Factorize)\n", variant);
        printf("==============================================\n");
        print_stat_factorize();
      }

    }

    inline
    void
    solveLowerOnDeviceVar0(const ordinal_type pbeg, 
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
    solveLowerOnDeviceVar1(const ordinal_type pbeg, 
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
    solveLowerOnDeviceVar2(const ordinal_type pbeg, 
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
    solveLowerOnDevice(const ordinal_type pbeg, 
                       const ordinal_type pend,
                       const size_type_array_host &h_buf_solve_ptr,
                       const value_type_matrix &t) {
      if (variant == 0) 
        solveLowerOnDeviceVar0(pbeg, pend, h_buf_solve_ptr, t);
      else if (variant == 1)
        solveLowerOnDeviceVar1(pbeg, pend, h_buf_solve_ptr, t);
      else if (variant == 2)
        solveLowerOnDeviceVar2(pbeg, pend, h_buf_solve_ptr, t);
      else {
        TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, 
                                 "LevelSetTools::solveLowerOnDevice, algorithm variant is not supported");
      }
    }

    inline
    void
    solveUpperOnDeviceVar0(const ordinal_type pbeg,
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
                Gemv<Trans::NoTranspose,Algo::OnDevice>
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
    solveUpperOnDeviceVar1(const ordinal_type pbeg,
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
                Gemv<Trans::NoTranspose,Algo::OnDevice>
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
    solveUpperOnDeviceVar2(const ordinal_type pbeg,
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
               
              Gemv<Trans::NoTranspose,Algo::OnDevice>
                ::invoke(_handle_blas, one, A, b, zero, tT); checkDeviceBlasStatus("gemv");
            }
          }
        }
      }
    }

    inline
    void
    solveUpperOnDevice(const ordinal_type pbeg,
                       const ordinal_type pend,
                       const size_type_array_host &h_buf_solve_ptr,
                       const value_type_matrix &t) {
      if (variant == 0) 
        solveUpperOnDeviceVar0(pbeg, pend, h_buf_solve_ptr, t);
      else if (variant == 1) 
        solveUpperOnDeviceVar1(pbeg, pend, h_buf_solve_ptr, t);
      else if (variant == 2) 
        solveUpperOnDeviceVar2(pbeg, pend, h_buf_solve_ptr, t);
      else {
        TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, 
                                 "LevelSetTools::solveUpperOnDevice, algorithm variant is not supported");
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
        const Kokkos::RangePolicy<exec_space> policy(0,_buf_solve_ptr.extent(0));
        const auto buf_solve_nrhs_ptr = _buf_solve_nrhs_ptr;
        const auto buf_solve_ptr = _buf_solve_ptr;
        Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const ordinal_type &i) {
            buf_solve_nrhs_ptr(i) = nrhs*buf_solve_ptr(i);
          });
        Kokkos::deep_copy(_h_buf_solve_nrhs_ptr, _buf_solve_nrhs_ptr);
      }
    }


    inline
    void
    solveCholesky(const value_type_matrix &x,   // solution
                  const value_type_matrix &b,   // right hand side
                  const value_type_matrix &t,
                  const ordinal_type verbose = 0) { // temporary workspace (store permuted vectors)
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
      Kokkos::Impl::Timer timer;

      stat.n_kernel_launching = 0;

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
                ++stat.n_kernel_launching;
              }
              const auto h_buf_solve_ptr = Kokkos::subview(_h_buf_solve_nrhs_ptr, range_solve_buf_ptr);              
              solveLowerOnDevice(pbeg, pend, h_buf_solve_ptr, t); 
              Kokkos::fence();
              
              Kokkos::parallel_for("update lower", 
                                   policy_update_with_work_property, 
                                   functor); 
              ++stat.n_kernel_launching;
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
              ++stat.n_kernel_launching;
              exec_space().fence(); //Kokkos::fence();

              if (lvl < _device_level_cut) {
                // do nothing
                //Kokkos::parallel_for("solve upper", policy_solve, functor); 
              } else {
                Kokkos::parallel_for("solve upper", 
                                     policy_solve_with_work_property,
                                     functor);
                ++stat.n_kernel_launching;
              }

              const auto h_buf_solve_ptr = Kokkos::subview(_h_buf_solve_nrhs_ptr, range_solve_buf_ptr);
              solveUpperOnDevice(pbeg, pend, h_buf_solve_ptr, t);
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
        printf("Summary: LevelSetTools-Variant-%d (ParallelSolve: %3d)\n", variant, nrhs);
        printf("=====================================================\n");
        print_stat_solve();
      }
    }

    
  };

}
#endif

