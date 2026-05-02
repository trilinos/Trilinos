// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
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

#include "Tacho_ApplyPermutation.hpp"
#include "Tacho_ApplyPermutation_OnDevice.hpp"

#include "Tacho_Scale2x2_BlockInverseDiagonals.hpp"
#include "Tacho_Scale2x2_BlockInverseDiagonals_OnDevice.hpp"

#include "Tacho_Spmv_OnDevice.hpp"

#include "Tacho_Chol_OnDevice.hpp"
#include "Tacho_NonPivLDL_OnDevice.hpp"
#include "Tacho_LDL_OnDevice.hpp"
#include "Tacho_LU_OnDevice.hpp"

#include "Tacho_GemmTriangular_OnDevice.hpp"
#include "Tacho_Gemm_OnDevice.hpp"
#include "Tacho_Gemv_OnDevice.hpp"
#include "Tacho_Herk_OnDevice.hpp"
#include "Tacho_Trmv_OnDevice.hpp"
#include "Tacho_Trsm_OnDevice.hpp"
#include "Tacho_Trsv_OnDevice.hpp"

#include "Tacho_SupernodeInfo.hpp"
#include "Tacho_TeamFunctor_ExtractCRS.hpp"

#include "Tacho_TeamFunctor_FactorizeChol.hpp"
#include "Tacho_TeamFunctor_SolveLowerChol.hpp"
#include "Tacho_TeamFunctor_SolveUpperChol.hpp"

#include "Tacho_TeamFunctor_FactorizeLDL.hpp"
#include "Tacho_TeamFunctor_SolveLowerLDL.hpp"
#include "Tacho_TeamFunctor_SolveUpperLDL.hpp"

#include "Tacho_TeamFunctor_FactorizeLU.hpp"
#include "Tacho_TeamFunctor_SolveLowerLU.hpp"
#include "Tacho_TeamFunctor_SolveUpperLU.hpp"

//#define TACHO_TEST_LEVELSET_TOOLS_KERNEL_OVERHEAD
//#define TACHO_ENABLE_LEVELSET_TOOLS_USE_LIGHT_KERNEL


namespace Tacho {

template <typename ValueType, typename DeviceType, int Var>
class NumericToolsLevelSet : public NumericToolsBase<ValueType, DeviceType> {
public:
  enum { variant = Var, max_factor_team_size = 64 };

  ///
  /// types
  ///
  using base_type = NumericToolsBase<ValueType, DeviceType>;
  using typename base_type::device_type;
  using typename base_type::exec_memory_space;
  using typename base_type::exec_space;
  using typename base_type::host_memory_space;
  using typename base_type::host_space;
  using typename base_type::ordinal_type_array;
  using typename base_type::ordinal_type_array_host;
  using typename base_type::range_type;
  using typename base_type::size_type_array;
  using typename base_type::size_type_array_host;
  using typename base_type::supernode_info_type;
  using typename base_type::supernode_type_array_host;
  using typename base_type::value_type;
  using typename base_type::mag_type;
  using typename base_type::int_type_array;
  using typename base_type::value_type_array;
  using typename base_type::value_type_matrix;

  using base_type::base_type;

private:
  using base_type::_aj;
  using base_type::_ap;
  using base_type::_ax;
  using base_type::_diag;
  using base_type::_info;
  using base_type::_m;
  using base_type::_nsupernodes;
  using base_type::_peri;
  using base_type::_perm;
  using base_type::_piv;
  using base_type::_stree_level;
  using base_type::_stree_roots;
  using base_type::_superpanel_buf;

  using base_type::print_stat_memory;
  using base_type::reset_stat;
  using base_type::stat;
  using base_type::track_alloc;
  using base_type::track_free;

  using rowptr_view = Kokkos::View<int *, device_type>;
  using colind_view = Kokkos::View<int *, device_type>;
  using nzvals_view = Kokkos::View<value_type *, device_type>;

  // supernode host information for level kernels launching
  supernode_type_array_host _h_supernodes;

  ///
  /// level set infrastructure and tuning parameters
  ///

  // 0: device level function, 1: team policy, 2: team policy recursive
  ordinal_type _device_factorize_thres, _device_solve_thres;
  ordinal_type _device_level_cut, _team_serial_level_cut;

  ordinal_type_array_host _h_factorize_mode, _h_solve_mode;
  ordinal_type_array _factorize_mode, _solve_mode;
  ordinal_type_array_host _h_num_device_calls_factor, _h_num_device_calls_solve;

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
  ordinal_type _nrhs;
  size_type_array_host _h_buf_solve_ptr, _h_buf_solve_nrhs_ptr;
  size_type_array _buf_solve_ptr, _buf_solve_nrhs_ptr;

  // workspace
  size_type _bufsize_factorize, _bufsize_solve;
  size_type _worksize;
  value_type_array _buf;
  value_type_array _work;

  // for using SpMV
  using SpMV_type = SpMV<supernode_info_type>;
  SpMV_type *_spmv;
  bool _keep_zeros;

  // common for host and cuda
  int _status;

  // cuda stream
  int _nstreams;
  bool _team_on_user_stream;

  // workspace for SpMV
  value_type_matrix _w_vec;
#if defined(KOKKOS_ENABLE_CUDA)
  bool _is_cublas_created, _is_cusolver_dn_created;
  cublasHandle_t _handle_blas;
  cusolverDnHandle_t _handle_lapack;
  using blas_handle_type = cublasHandle_t;
  using lapack_handle_type = cusolverDnHandle_t;
  using stream_array_host = std::vector<cudaStream_t>;
  #define getBlasHandle(id)   _handle_blas
  #define getLapackHandle(id) _handle_lapack
#elif defined(KOKKOS_ENABLE_HIP)
  bool _is_rocblas_created;
  rocblas_handle _handle_lapack; // just used for workspace size query
  std::vector<rocblas_handle> _handles;
  using blas_handle_type = rocblas_handle;
  using lapack_handle_type = rocblas_handle;
  using stream_array_host = std::vector<hipStream_t>;
  #define getBlasHandle(id)   _handles[id]
  #define getLapackHandle(id) _handles[id]
#else
  int _handle_blas, _handle_lapack; // dummy handle for convenience
  using blas_handle_type = int;
  using lapack_handle_type = int;
  using stream_array_host = std::vector<int>;
  #define getBlasHandle()   _handle_blas
  #define getLapackHandle() _handle_lapack
#endif
  stream_array_host _streams;

  using exec_instance_array_host = std::vector<exec_space>;
  exec_instance_array_host _exec_instances;

  ///
  /// statistics
  ///
  struct {
    int n_level;
    int n_device_factorize, n_team_factorize, n_kernel_launching_factorize;
    int n_device_solve, n_team_solve, n_kernel_launching_solve;
    int n_kernel_launching;
  } stat_level;

  ///
  /// error check for cuda things
  ///
  inline void checkStatus(const char *func, const char *lib) {
    if (_status != 0) {
      printf("Error: %s, %s returns non-zero status %d\n", lib, func, _status);
      std::runtime_error("checkStatus failed");
    }
  }
  inline void checkDeviceLapackStatus(const char *func) {
#if defined(KOKKOS_ENABLE_CUDA)
    constexpr bool is_host = std::is_same<exec_memory_space, Kokkos::HostSpace>::value;
    checkStatus(func, is_host ? "HostLapack" : "CuSolverDn");
#elif defined(KOKKOS_ENABLE_HIP)
    constexpr bool is_host = std::is_same<exec_memory_space, Kokkos::HostSpace>::value;
    checkStatus(func, is_host ? "HostLapack" : "RocSolver");
#else
    checkStatus(func, "HostLapack");
#endif
  }
  inline void checkDeviceBlasStatus(const char *func) {
#if defined(KOKKOS_ENABLE_CUDA)
    constexpr bool is_host = std::is_same<exec_memory_space, Kokkos::HostSpace>::value;
    checkStatus(func, is_host ? "HostBlas" : "CuBlas");
#elif defined(KOKKOS_ENABLE_HIP)
    constexpr bool is_host = std::is_same<exec_memory_space, Kokkos::HostSpace>::value;
    checkStatus(func, is_host ? "HostBlas" : "RocBlas");
#else
    checkStatus(func, "HostBlas");
#endif
  }
  inline void checkDeviceStatus(const char *func) {
#if defined(KOKKOS_ENABLE_CUDA)
    constexpr bool is_host = std::is_same<exec_memory_space, Kokkos::HostSpace>::value;
    checkStatus(func, is_host ? "Host" : "Cuda");
#elif defined(KOKKOS_ENABLE_HIP)
    constexpr bool is_host = std::is_same<exec_memory_space, Kokkos::HostSpace>::value;
    checkStatus(func, is_host ? "Host" : "HIP");
#else
    checkStatus(func, "Host");
#endif
  }

  inline void print_stat_init() override {
    base_type::print_stat_init();
    const double kilo(1024);
    printf("  Time\n");
    printf("             time for initialization:                         %10.6f s\n", stat.t_init);
    printf("             time for compute mode classification:            %10.6f s\n", stat.t_mode_classification);
    printf("             total time spent:                                %10.6f s\n",
           (stat.t_init + stat.t_mode_classification));
    printf("\n");
    printf("  Memory\n");
    printf("             workspace allocated for solve:                   %10.3f MB\n", stat.m_used / kilo / kilo);
    printf("             peak memory used:                                %10.3f MB\n", stat.m_peak / kilo / kilo);
    printf("\n");
    printf("  Compute Mode in Factorize with a Threshold(%d)\n", _device_factorize_thres);
    printf("             # of levels:                                     %6d\n", _nlevel);
    printf("             # of subproblems using device functions:         %6d\n", stat_level.n_device_factorize);
    printf("             # of subproblems using team functions:           %6d\n", stat_level.n_team_factorize);
    printf("             total # of subproblems:                          %6d\n",
           (stat_level.n_device_factorize + stat_level.n_team_factorize));
    printf("\n");
    printf("  Compute Mode in Solve with a Threshold(%d)\n", _device_solve_thres);
    printf("             # of subproblems using device functions:         %6d\n", stat_level.n_device_solve);
    printf("             # of subproblems using team functions:           %6d\n", stat_level.n_team_solve);
    printf("             total # of subproblems:                          %6d\n",
           (stat_level.n_device_solve + stat_level.n_team_solve));
    printf("\n");
  }

  inline void print_stat_factor() override {
    base_type::print_stat_factor();
    double flop = 0;
    switch (this->getSolutionMethod()) {
    case 0:   /// LDL no-pivot
    case 1: { /// Cholesky
      for (ordinal_type sid = 0; sid < _nsupernodes; ++sid) {
        auto &s = _h_supernodes(sid);
        const ordinal_type m = s.m, n = s.n - s.m;
        flop += DenseFlopCount<value_type>::Chol(m);
        if (variant == 1) {
          flop += DenseFlopCount<value_type>::Trsm(true, m, m);
        } else if (variant == 2 || variant == 3) {
          flop += DenseFlopCount<value_type>::Trsm(true, m, m);
          flop += DenseFlopCount<value_type>::Trsm(true, m, n);
        }
        flop += DenseFlopCount<value_type>::Trsm(true, m, n);
        flop += DenseFlopCount<value_type>::Syrk(m, n);
      }
      break;
    }
    case 2: {
      for (ordinal_type sid = 0; sid < _nsupernodes; ++sid) {
        auto &s = _h_supernodes(sid);
        const ordinal_type m = s.m, n = s.n - s.m;
        flop += DenseFlopCount<value_type>::LDL(m);
        if (variant == 1) {
          flop += DenseFlopCount<value_type>::Trsm(true, m, m);
        } else if (variant == 2 || variant == 3) {
          flop += DenseFlopCount<value_type>::Trsm(true, m, m);
          flop += DenseFlopCount<value_type>::Trsm(true, m, n);
        }
        flop += DenseFlopCount<value_type>::Trsm(true, m, n);
        flop += DenseFlopCount<value_type>::Syrk(m, n);
      }
      break;
    }
    case 3: {
      for (ordinal_type sid = 0; sid < _nsupernodes; ++sid) {
        auto &s = _h_supernodes(sid);
        const ordinal_type m = s.m, n = s.n - s.m;
        flop += DenseFlopCount<value_type>::LU(m, m);
        if (variant == 1) {
          flop += 2 * DenseFlopCount<value_type>::Trsm(true, m, m);
        } else if (variant == 2 || variant == 3) {
          flop += 2 * DenseFlopCount<value_type>::Trsm(true, m, m);
          flop += 2 * DenseFlopCount<value_type>::Trsm(true, m, n);
        }
        flop += 2 * DenseFlopCount<value_type>::Trsm(true, m, n);
        flop += DenseFlopCount<value_type>::Gemm(n, n, m);
      }
      break;
    }
    default: {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, "The solution method is not supported");
    }
    }
    const double kilo(1024);
    printf("  FLOPs\n");
    printf("             gflop   for numeric factorization:               %10.3f GFLOP\n", flop / kilo / kilo / kilo);
    printf("             gflop/s for numeric factorization:               %10.3f GFLOP/s\n",
           flop / stat.t_factor / kilo / kilo / kilo);
    printf("\n");
    printf("  Kernels\n");
    printf("             # of kernels launching:                          %6d\n", stat_level.n_kernel_launching);
    printf("\n");
  }

  inline void print_stat_solve() override {
    base_type::print_stat_solve();
    printf("  Kernels\n");
    printf("             # of kernels launching:                          %6d\n", stat_level.n_kernel_launching);
    printf("\n");
  }

public:

  ///
  /// initialization / release
  ///
  inline void initialize(const ordinal_type device_level_cut, const ordinal_type device_factorize_thres,
                         const ordinal_type device_solve_thres, const int nstreams_in = 1, const bool team_on_user_stream = false,
                         const bool store_transpose = false, const ordinal_type verbose = 0) {
    stat_level.n_level = 0;
    stat_level.n_device_factorize = 0;
    stat_level.n_device_solve = 0;
    stat_level.n_team_factorize = 0;
    stat_level.n_team_solve = 0;

    Kokkos::Timer timer;

    timer.reset();
    // # of streams needs to be at least 1
    const int nstreams = max(1, nstreams_in);

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
      for (ordinal_type sid = 0; sid < _nsupernodes; ++sid)
        _nlevel = max(_stree_level(sid), _nlevel);
      ++_nlevel;
    }
    stat_level.n_level = _nlevel;

    // create level ptr
    _h_level_ptr = size_type_array_host("h_level_ptr", _nlevel + 1);
    {
      // first count # of supernodes in each level
      for (ordinal_type sid = 0; sid < _nsupernodes; ++sid)
        ++_h_level_ptr(_stree_level(sid) + 1);

      // scan
      for (ordinal_type i = 0; i < _nlevel; ++i)
        _h_level_ptr(i + 1) += _h_level_ptr(i);
    }

    // fill sids
    _h_level_sids = ordinal_type_array_host(do_not_initialize_tag("h_level_sids"), _nsupernodes);
    {
      size_type_array_host tmp_level_ptr(do_not_initialize_tag("tmp_level_ptr"), _h_level_ptr.extent(0));
      Kokkos::deep_copy(tmp_level_ptr, _h_level_ptr);
      for (ordinal_type sid = 0; sid < _nsupernodes; ++sid) {
        const ordinal_type lvl = _stree_level(sid);
        _h_level_sids(tmp_level_ptr(lvl)++) = sid;
      }
    }
    _level_sids = Kokkos::create_mirror_view_and_copy(exec_memory_space(), _h_level_sids);
    track_alloc(_level_sids.span() * sizeof(ordinal_type));

    ///
    /// workspace
    ///
    _h_buf_level_ptr = ordinal_type_array_host(do_not_initialize_tag("h_buf_factor_level_ptr"), _nlevel + 1);
    {
      _h_buf_level_ptr(0) = 0;
      for (ordinal_type i = 0; i < _nlevel; ++i) {
        const ordinal_type pbeg = _h_level_ptr(i), pend = _h_level_ptr(i + 1);
        _h_buf_level_ptr(i + 1) = (pend - pbeg + 1) + _h_buf_level_ptr(i);
      }
    }

    // create workspace for factorization / solve
    _bufsize_factorize = 0;
    _bufsize_solve = 0;
    _h_buf_factor_ptr = size_type_array_host(do_not_initialize_tag("h_buf_factor_ptr"), _h_buf_level_ptr(_nlevel));
    _h_buf_solve_ptr = size_type_array_host(do_not_initialize_tag("h_buf_solve_ptr"), _h_buf_level_ptr(_nlevel));
    {
      const ordinal_type method_id = this->getSolutionMethod();
      for (ordinal_type i = 0; i < _nlevel; ++i) {
        const ordinal_type lbeg = _h_buf_level_ptr(i);
        const ordinal_type pbeg = _h_level_ptr(i), pend = _h_level_ptr(i + 1);

        _h_buf_factor_ptr(lbeg) = 0;
        _h_buf_solve_ptr(lbeg) = 0;
        for (ordinal_type p = pbeg, k = (lbeg + 1); p < pend; ++p, ++k) {
          const ordinal_type sid = _h_level_sids(p);
          const auto s = _h_supernodes(sid);
          const ordinal_type m = s.m, n = s.n, n_m = n - m;
          const ordinal_type schur_work_size = n_m * (n_m + max_factor_team_size);
          const ordinal_type chol_factor_work_size_variants[4] = {schur_work_size,
                                                                  max(m * m, schur_work_size),
                                                                  m * m + schur_work_size,
                                                                  m * m + schur_work_size};
          const ordinal_type chol_factor_work_size = chol_factor_work_size_variants[variant] + (method_id == 0 ? m*n_m : 0);
          const ordinal_type ldl_factor_work_size_variant_0 = chol_factor_work_size_variants[0] + max(32 * m, m * n);
          const ordinal_type ldl_factor_work_size_variants[4] = {ldl_factor_work_size_variant_0,
                                                                 max(m * m, ldl_factor_work_size_variant_0 + m * n_m),
                                                                 m * m + ldl_factor_work_size_variant_0 + m * n_m,
                                                                 m * m + ldl_factor_work_size_variant_0 + m * n_m};
          const ordinal_type ldl_factor_work_size = ldl_factor_work_size_variants[variant];
          const ordinal_type lu_factor_work_size_variants[4] = {schur_work_size, max(m * m, schur_work_size),
                                                                m * m + schur_work_size,
                                                                m * m + schur_work_size};
          const ordinal_type lu_factor_work_size = lu_factor_work_size_variants[variant];
          const ordinal_type factor_work_size_variants[4] = {chol_factor_work_size, ldl_factor_work_size,
                                                             lu_factor_work_size,
                                                             lu_factor_work_size};

          const ordinal_type chol_solve_work_size = (variant == 0 ? n_m : n);
          const ordinal_type ldl_solve_work_size = chol_solve_work_size;
          const ordinal_type lu_solve_work_size = chol_solve_work_size;
          const ordinal_type solve_work_size_variants[4] = {chol_solve_work_size, ldl_solve_work_size,
                                                            lu_solve_work_size,
                                                            lu_solve_work_size};

          const ordinal_type index_work_size = (method_id-1 < 0 ? 0 : method_id-1);
          const ordinal_type factor_work_size = factor_work_size_variants[index_work_size];
          const ordinal_type solve_work_size = solve_work_size_variants[index_work_size];

          _h_buf_factor_ptr(k) = factor_work_size + _h_buf_factor_ptr(k - 1);
          _h_buf_solve_ptr(k) = solve_work_size + _h_buf_solve_ptr(k - 1);
        }
        const ordinal_type last_idx = lbeg + pend - pbeg;
        _bufsize_factorize = max(_bufsize_factorize, _h_buf_factor_ptr(last_idx));
        _bufsize_solve = max(_bufsize_solve, _h_buf_solve_ptr(last_idx));
      }
    }

    _buf_factor_ptr = Kokkos::create_mirror_view_and_copy(exec_memory_space(), _h_buf_factor_ptr);
    track_alloc(_buf_factor_ptr.span() * sizeof(size_type));

    _buf_solve_ptr = Kokkos::create_mirror_view_and_copy(exec_memory_space(), _h_buf_solve_ptr);
    track_alloc(_buf_solve_ptr.span() * sizeof(size_type));

    _h_buf_solve_nrhs_ptr =
        size_type_array_host(do_not_initialize_tag("h_buf_solve_nrhs_ptr"), _h_buf_solve_ptr.extent(0));
    Kokkos::deep_copy(_h_buf_solve_nrhs_ptr, _h_buf_solve_ptr);
    _buf_solve_nrhs_ptr = Kokkos::create_mirror_view_and_copy(exec_memory_space(), _h_buf_solve_nrhs_ptr);
    track_alloc(_buf_solve_nrhs_ptr.span() * sizeof(size_type));

    ///
    /// cuda library initialize
    ///
#if defined(KOKKOS_ENABLE_CUDA)
    if (!_is_cublas_created) {
      _status = cublasCreate(&_handle_blas);
      checkDeviceBlasStatus("cublasCreate");
      _is_cublas_created = true;
    }
    if (!_is_cusolver_dn_created) {
      _status = cusolverDnCreate(&_handle_lapack);
      checkDeviceLapackStatus("cusolverDnCreate");
      _is_cusolver_dn_created = true;
    }
#endif
#if defined(KOKKOS_ENABLE_HIP)
    if (!_is_rocblas_created) {
      _status = rocblas_create_handle(&_handle_lapack);
      checkDeviceBlasStatus("rocblasCreate");
      _is_rocblas_created = true;
    }
#endif
    // pre-allocate buf
    _nrhs = 1;
    Kokkos::resize(_buf, max(_bufsize_factorize, _bufsize_solve));
    track_alloc(_buf.span() * sizeof(value_type));
    // pre-allocate work
    _worksize = 0;
    switch (this->getSolutionMethod()) {
    case 0:   /// LDL no-pivot
    case 1: { /// Cholesky
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      value_type_matrix T(NULL, _info.max_supernode_size, _info.max_supernode_size);
      _worksize = Chol<Uplo::Upper, Algo::OnDevice>::invoke(_handle_lapack, T, _work);
#endif
      break;
    }
    case 2: { /// LDL
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      value_type_matrix T(NULL, _info.max_supernode_size, _info.max_supernode_size);
      ordinal_type_array P(NULL, _info.max_supernode_size);
      _worksize = LDL<Uplo::Lower, Algo::OnDevice>::invoke(_handle_lapack, T, P, _work);
#else
      _worksize = 32 * _info.max_supernode_size;
#endif
      break;
    }
    case 3: { /// LU
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      value_type_matrix T(NULL, _info.max_supernode_size, _info.max_num_cols);
      ordinal_type_array P(NULL, std::min(_info.max_supernode_size, _info.max_num_cols));
      _worksize = LU<Algo::OnDevice>::invoke(_handle_lapack, T, P, _work);
#endif
      break;
    }
    }
    if (this->getSolutionMethod() == 0) {
      // worksize for on-device non-pivot LDL
      // the blocksize is used in Tacho_NonPivLDL_OnDevice
      char *nb_env = getenv("TACHO_BLOCK_SIZE");
      const ordinal_type nb = (nb_env == NULL ? 256 : atoi(nb_env));
      _worksize += (nb*_info.max_supernode_size);
    }
    size_type worksize = _worksize * (nstreams + 1);
    Kokkos::resize(_work, worksize);
    track_alloc(_work.span() * sizeof(value_type));
    stat.t_init = timer.seconds();

    ///
    /// classification of problems
    ///
    timer.reset();

    _device_level_cut = min(device_level_cut, _nlevel);
    _device_factorize_thres = device_factorize_thres;
    _device_solve_thres = (variant == 3 ? 0 : device_solve_thres);

    _h_factorize_mode = ordinal_type_array_host(do_not_initialize_tag("h_factorize_mode"), _nsupernodes);
    Kokkos::deep_copy(_h_factorize_mode, -1);

    _h_solve_mode = ordinal_type_array_host(do_not_initialize_tag("h_solve_mode"), _nsupernodes);
    Kokkos::deep_copy(_h_solve_mode, -1);

    _h_num_device_calls_factor = ordinal_type_array_host(do_not_initialize_tag("h_num_device_calls_factor"), _nlevel);
    _h_num_device_calls_solve = ordinal_type_array_host(do_not_initialize_tag("h_num_device_calls_solve"), _nlevel);

    if (_device_level_cut > 0) {
      for (ordinal_type lvl = 0; lvl < _device_level_cut; ++lvl) {
        _h_num_device_calls_solve(lvl) = 0;
        _h_num_device_calls_factor(lvl) = 0;

        const ordinal_type pbeg = _h_level_ptr(lvl), pend = _h_level_ptr(lvl + 1);
        for (ordinal_type p = pbeg; p < pend; ++p) {
          const ordinal_type sid = _h_level_sids(p);
          _h_solve_mode(sid) = 0;
          _h_factorize_mode(sid) = 0;
          ++stat_level.n_device_solve;
          ++stat_level.n_device_factorize;

          const auto s = _h_supernodes(sid);
          const ordinal_type m = s.m;
          if (m > _device_solve_thres) {
            _h_num_device_calls_solve(lvl) ++;
          }
          if (m > _device_factorize_thres) {
            _h_num_device_calls_factor(lvl) ++;
          }
        }
      }
    }

    _team_serial_level_cut = _nlevel;
    {
      for (ordinal_type lvl = _device_level_cut; lvl < _team_serial_level_cut; ++lvl) {
        _h_num_device_calls_solve(lvl) = 0;
        _h_num_device_calls_factor(lvl) = 0;

        const ordinal_type pbeg = _h_level_ptr(lvl), pend = _h_level_ptr(lvl + 1);
        for (ordinal_type p = pbeg; p < pend; ++p) {
          const ordinal_type sid = _h_level_sids(p);
          const auto s = _h_supernodes(sid);
          const ordinal_type m = s.m;    //, n_m = s.n-s.m;
          if (m > _device_solve_thres) { // || n > _device_solve_thres)
            _h_solve_mode(sid) = 0;
            _h_num_device_calls_solve(lvl) ++;
            ++stat_level.n_device_solve;
          } else {
            _h_solve_mode(sid) = 1;
            ++stat_level.n_team_solve;
          }
          if (m > _device_factorize_thres) { // || n_m > _device_factorize_thres)
            _h_factorize_mode(sid) = 0;
            _h_num_device_calls_factor(lvl) ++;
            ++stat_level.n_device_factorize;
          } else {
            _h_factorize_mode(sid) = 1;
            ++stat_level.n_team_factorize;
          }
        }
      }
    }

    _factorize_mode = Kokkos::create_mirror_view_and_copy(exec_memory_space(), _h_factorize_mode);
    track_alloc(_factorize_mode.span() * sizeof(ordinal_type));

    _solve_mode = Kokkos::create_mirror_view_and_copy(exec_memory_space(), _h_solve_mode);
    track_alloc(_solve_mode.span() * sizeof(ordinal_type));

    _team_on_user_stream = team_on_user_stream;
    createStream(nstreams, verbose);
    if (variant == 3 && _keep_zeros) {
      // compress each partitioned inverse at each level into CRS matrix
      setupCRS(store_transpose, verbose);
    }
    stat.t_mode_classification = timer.seconds();

    if (verbose) {
      switch (this->getSolutionMethod()) {
      case 0: {
        printf("Summary: LevelSetTools-Variant-%d (InitializeLDL (no pivot))\n", variant);
        printf("============================================================\n");
        break;
      }
      case 1: {
        printf("Summary: LevelSetTools-Variant-%d (InitializeCholesky)\n", variant);
        printf("======================================================\n");
        break;
      }
      case 2: {
        printf("Summary: LevelSetTools-Variant-%d (InitializeLDL)\n", variant);
        printf("=================================================\n");
        break;
      }
      case 3: {
        printf("Summary: LevelSetTools-Variant-%d (InitializeLU)\n", variant);
        printf("================================================\n");
        break;
      }
      }
      print_stat_init();
      printf("  Execution Mode\n");
      printf("             # of streams:                                      %d\n",nstreams);
      printf("               Team kernels on %s\n",(_team_on_user_stream ? "User Stream-0" : "Default Stream"));
      printf("\n");
      fflush(stdout);
    }
  }

  inline void release(const ordinal_type verbose = 0) override {
    base_type::release(false);
    if (variant == 3) {
      this->releaseCRS(true, verbose);
    }
    track_free(_buf_factor_ptr.span() * sizeof(size_type));
    track_free(_buf_solve_ptr.span() * sizeof(size_type));
    track_free(_buf_solve_nrhs_ptr.span() * sizeof(size_type));
    track_free(_factorize_mode.span() * sizeof(ordinal_type));
    track_free(_solve_mode.span() * sizeof(ordinal_type));
    track_free(_level_sids.span() * sizeof(ordinal_type));

    track_free(_buf.span() * sizeof(value_type));
    track_free(_work.span() * sizeof(value_type));
    Kokkos::resize(_buf, 0);
    Kokkos::resize(_work, 0);

    if (verbose) {
      printf("Summary: LevelSetTools-Variant-%d (Release)\n", variant);
      printf("===========================================\n");
      print_stat_memory();
      fflush(stdout);
    }
  }

  NumericToolsLevelSet() : base_type() {
    _nlevel = 0;
    _bufsize_factorize = 0;
    _bufsize_solve = 0;
    _nstreams = 0;
    _team_on_user_stream = false;
    _spmv = nullptr;
    _keep_zeros = false;
    stat_level = stat_level();
  }
  NumericToolsLevelSet(const NumericToolsLevelSet &b) = default;

  NumericToolsLevelSet(const ordinal_type method,
                       // input matrix A
                       const ordinal_type m, const size_type_array &ap, const ordinal_type_array &aj,
                       // input permutation
                       const ordinal_type_array &perm, const ordinal_type_array &peri,
                       // supernodes
                       const ordinal_type nsupernodes, const ordinal_type_array &supernodes,
                       const size_type_array &gid_ptr, const ordinal_type_array &gid_colidx,
                       const size_type_array &sid_ptr, const ordinal_type_array &sid_colidx,
                       const ordinal_type_array &blk_colidx, const ordinal_type_array &stree_parent,
                       const size_type_array &stree_ptr, const ordinal_type_array &stree_children,
                       const ordinal_type_array_host &stree_level, const ordinal_type_array_host &stree_roots)
      : base_type(method, m, ap, aj, perm, peri, nsupernodes, supernodes, gid_ptr, gid_colidx, sid_ptr, sid_colidx,
                  blk_colidx, stree_parent, stree_ptr, stree_children, stree_level, stree_roots) {
    _keep_zeros = false;
    _nstreams = 0;
    _team_on_user_stream = false;
#if defined(KOKKOS_ENABLE_CUDA)
    _is_cublas_created = 0;
    _is_cusolver_dn_created = 0;
#endif
#if defined(KOKKOS_ENABLE_HIP)
    _is_rocblas_created = 0;
#endif
    if (variant == 3)
      _spmv = new SpMV_type(_keep_zeros);
    else
      _spmv = nullptr;
  }

  virtual ~NumericToolsLevelSet() {
    /// kokkos execution space may fence and it uses the wrapped stream when it is deallocated   
    /// on cuda, deallocting streams first does not cause any errors while hip generates errors.
    /// here, we just follow the consistent destruction process as hip does.
    _exec_instances.clear();

#if defined(KOKKOS_ENABLE_CUDA)
    if (_is_cusolver_dn_created) {
      _status = cusolverDnDestroy(_handle_lapack);
      checkDeviceLapackStatus("cusolverDnDestroy");
    }
    if (_is_cublas_created) {
      _status = cublasDestroy(_handle_blas);
      checkDeviceBlasStatus("cublasDestroy");
    }

    for (ordinal_type i = 0; i < _nstreams; ++i) {
      _status = cudaStreamDestroy(_streams[i]);
      checkDeviceStatus("cudaStreamDestroy");
    }
#endif
#if defined(KOKKOS_ENABLE_HIP)
    if (_is_rocblas_created) {
      _status = rocblas_destroy_handle(_handle_lapack);
      checkDeviceLapackStatus("rocblasDestroy");
    }
    for (ordinal_type i = 0; i < _nstreams; ++i) {
      _status = rocblas_destroy_handle(_handles[i]);
      checkDeviceLapackStatus("rocblasDestroy(handles[i])");
    }
    _handles.clear();

    for (ordinal_type i = 0; i < _nstreams; ++i) {
      _status = hipStreamDestroy(_streams[i]);
      checkDeviceStatus("cudaStreamDestroy");
    }
#endif
    _streams.clear();
    if (_spmv != nullptr) {
      delete _spmv;
      _spmv = nullptr;
    }
  }

  inline void createStream(const ordinal_type nstreams, const ordinal_type verbose = 0) {
    // # of streams needs to be at least 1
    if (nstreams <= 0) return;
    _nstreams = nstreams;

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    if (_streams.size() == size_t(nstreams)) return; // nothing to do

#if defined(KOKKOS_ENABLE_CUDA)
    // destroy previously created streams
    for (size_t i = 0; i < _streams.size(); ++i) {
      _status = cudaStreamDestroy(_streams[i]);
      checkDeviceStatus("cudaStreamDestroy");
    }
    _streams.clear();
    // new streams
    _streams.resize(_nstreams);
    for (ordinal_type i = 0; i < _nstreams; ++i) {
      _status = cudaStreamCreateWithFlags(&_streams[i], cudaStreamNonBlocking);
      checkDeviceStatus("cudaStreamCreate");
    }
#endif
#if defined(KOKKOS_ENABLE_HIP)
    // destroy previously created streams
    for (size_t i = 0; i < _streams.size(); ++i) {
      _status = rocblas_destroy_handle(_handles[i]);
      checkDeviceLapackStatus("rocblasDestroy");
      _status = hipStreamDestroy(_streams[i]);
      checkDeviceStatus("hipStreamDestroy");
    }
    // new streams
    _streams.clear();
    _streams.resize(_nstreams);
    _handles.resize(_nstreams);
    for (ordinal_type i = 0; i < _nstreams; ++i) {
      _status = rocblas_create_handle(&_handles[i]);
      checkDeviceStatus("rocblas_create_handle");
      //_status = hipStreamCreateWithFlags(&_streams[i], hipStreamDefault);
      _status = hipStreamCreateWithFlags(&_streams[i], hipStreamNonBlocking);
      checkDeviceStatus("hipStreamCreate");
      _status = rocblas_set_stream(_handles[i], _streams[i]);
      checkDeviceBlasStatus("rocblasSetStream(handles[qid])");
    }
#endif
    // reinitialize execution instances with the streams
    _exec_instances.clear();
    _exec_instances.resize(_nstreams);
    for (ordinal_type i = 0; i < _nstreams; ++i) {
      ExecSpaceFactory<exec_space>::createInstance(_streams[i], _exec_instances[i]);
    }
#else
    // just one dummy stream
    _streams.clear();
    _streams.resize(1);
    _streams[0] = 0;
    // just one default execution space.
    _exec_instances.clear();
    _exec_instances.resize(1);
    _exec_instances[0] = exec_space();
#endif
    if (verbose) {
      printf("Summary: CreateStream : %3d\n", _nstreams);
      printf("===========================\n");
      fflush(stdout);
    }
  }

  inline void setStreamOnHandle(const ordinal_type qid) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    const auto mystream = _streams[qid];
#if defined(KOKKOS_ENABLE_CUDA)
    _status = cublasSetStream(_handle_blas, mystream);
    checkDeviceBlasStatus("cublasSetStream");

    _status = cusolverDnSetStream(_handle_lapack, mystream);
    checkDeviceLapackStatus("cusolverDnSetStream");
#endif
#if defined(KOKKOS_ENABLE_HIP)
    // > already set in createStream()
    //_status = rocblas_set_stream(_handles[qid], mystream);
    //checkDeviceBlasStatus("rocblasSetStream(handles[qid])");
#endif
#endif
  }

  ///
  /// Device level functions
  ///

  ///
  /// Non-pivot LDL
  ///
  inline int factorizeNoPivotLDLOnDeviceVar0(const ordinal_type pbeg, const ordinal_type pend,
                                             const size_type_array_host &h_buf_factor_ptr,
                                             const value_type_array &work) {
    const value_type one(1), minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    int num_device_calls = 0;
    exec_space exec_instance;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_factorize_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type handle_blas = getBlasHandle(qid);

        setStreamOnHandle(qid);
        exec_instance = _exec_instances[qid];

        const size_type worksize = work.extent(0) / _nstreams;
        value_type_array W(work.data() + worksize * qid, worksize);
        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
        value_type_array W = work;
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            value_type *aptr = s.u_buf;
            UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            aptr += m * m;
            // Calling fall-back default LDL_nopiv<Algo::OnDevice>::invoke with memeber = exec_instance
            _status = LDL_nopiv<Uplo::Upper, Algo::OnDevice>::invoke(handle_blas, exec_instance, ATL, W);
            checkDeviceLapackStatus("chol");

            if (n_m > 0) {
              UnmanagedViewType<value_type_matrix> ABR(_buf.data() + h_buf_factor_ptr(p - pbeg), n_m, n_m);
              UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // aptr += m*n_m;

              // Apply L^{-1} on off-diagonal
              _status = Trsm<Side::Left, Uplo::Upper, Trans::ConjTranspose, Algo::OnDevice>::invoke(
                  handle_blas, Diag::Unit(), one, ATL, ATR);
              checkDeviceBlasStatus("trsm");

              // Save ATR in workspace
              UnmanagedViewType<value_type_matrix> T(_buf.data() + h_buf_factor_ptr(p - pbeg) + (n_m * n_m), m, n_m);
              _status = Copy<Algo::OnDevice>::invoke(exec_instance, T, ATR);

              // Apply D^{-1} on off-diagonal
              _status = Scale_BlockInverseDiagonals<Side::Left, Algo::OnDevice>::invoke(exec_instance, ATL, ATR);

              // ABR = -ATR*W
              _status = GemmTriangular<Trans::Transpose, Trans::NoTranspose, Uplo::Upper, Algo::OnDevice>::invoke(
                  handle_blas, minus_one, ATR, T, zero, ABR);
              checkDeviceBlasStatus("gemm");
            }
            num_device_calls ++;
          }
        }
      }
    }
    return num_device_calls;
  }

  inline int factorizeNoPivotLDLOnDeviceVar1(const ordinal_type pbeg, const ordinal_type pend,
                                             const size_type_array_host &h_buf_factor_ptr,
                                             const value_type_array &work) {
    const value_type one(1), minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    int num_device_calls = 0;
    exec_space exec_instance;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_factorize_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type   handle_blas   = getBlasHandle(qid);

        setStreamOnHandle(qid);
        exec_instance = _exec_instances[qid];

        const size_type worksize = work.extent(0) / _nstreams;
        value_type_array W(work.data() + worksize * qid, worksize);
        ++q;
#else
        blas_handle_type   handle_blas   = getBlasHandle();
        value_type_array W = work;
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            // Factor the diagonal block
            value_type *aptr = s.u_buf;
            UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            aptr += m * m;

            // Calling fall-back default LDL_nopiv<Algo::OnDevice>::invoke with memeber = exec_instance
            _status = LDL_nopiv<Uplo::Upper, Algo::OnDevice>::invoke(handle_blas, exec_instance, ATL, W);
            checkDeviceLapackStatus("chol");

            // Compute inverse of diagonal block (in T)
            value_type *bptr = _buf.data() + h_buf_factor_ptr(p - pbeg);
            UnmanagedViewType<value_type_matrix> T(bptr, m, m); // shared with ABR
            _status = SetIdentity<Algo::OnDevice>::invoke(exec_instance, T, one);
            checkDeviceBlasStatus("SetIdentity");

            _status = Trsm<Side::Left, Uplo::Upper, Trans::NoTranspose, Algo::OnDevice>::invoke(
                handle_blas, Diag::Unit(), one, ATL, T);

            // Copy original diagonal into T
            using policy_type = Kokkos::RangePolicy<exec_space>;
            const auto policy = policy_type(exec_instance, 0, m);
            Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const ordinal_type &i) {
              T(i,i) = ATL(i,i);
            });
            checkDeviceBlasStatus("trsm");
            if (n_m > 0) {
              UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // aptr += m*n_m;
              UnmanagedViewType<value_type_matrix> ABR(bptr, n_m, n_m); // shared with T
              UnmanagedViewType<value_type_matrix> K(bptr+(n_m * n_m), m, n_m);

              _status = Trsm<Side::Left, Uplo::Upper, Trans::ConjTranspose, Algo::OnDevice>::invoke(
                  handle_blas, Diag::Unit(), one, ATL, ATR);
              checkDeviceBlasStatus("trsm");

              // Copy T back to ATL (inverse)
              _status = Copy<Algo::OnDevice>::invoke(exec_instance, ATL, T);
              checkDeviceBlasStatus("Copy");

              // Save ATR in workspace
              _status = Copy<Algo::OnDevice>::invoke(exec_instance, K, ATR);

              // Apply D^{-1} on off-diagonal
              _status = Scale_BlockInverseDiagonals<Side::Left, Algo::OnDevice>::invoke(exec_instance, ATL, ATR);

              // ABR = -ATR*K
              _status = GemmTriangular<Trans::Transpose, Trans::NoTranspose, Uplo::Upper, Algo::OnDevice>::invoke(
                  handle_blas, minus_one, ATR, K, zero, ABR);
              checkDeviceBlasStatus("gemm");
            } else {
              _status = Copy<Algo::OnDevice>::invoke(exec_instance, ATL, T);
              checkDeviceBlasStatus("Copy");
            }
            num_device_calls ++;
          }
        }
      }
    }
    return num_device_calls;
  }

  inline int factorizeNoPivotLDLOnDeviceVar2(const ordinal_type pbeg, const ordinal_type pend,
                                             const size_type_array_host &h_buf_factor_ptr,
                                             const value_type_array &work) {
    const value_type one(1), minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    int num_device_calls = 0;
    exec_space exec_instance;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_factorize_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type   handle_blas   = getBlasHandle(qid);

        setStreamOnHandle(qid);
        exec_instance = _exec_instances[qid];

        const size_type worksize = work.extent(0) / _nstreams;
        value_type_array W(work.data() + worksize * qid, worksize);
        ++q;
#else
        blas_handle_type   handle_blas   = getBlasHandle();
        value_type_array W = work;
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            value_type *aptr = s.u_buf;
            UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            aptr += m * m;

            // Calling fall-back default LDL_nopiv<Algo::OnDevice>::invoke with memeber = exec_instance
            _status = LDL_nopiv<Uplo::Upper, Algo::OnDevice>::invoke(handle_blas, exec_instance, ATL, W);
            checkDeviceLapackStatus("chol");

            value_type *bptr = _buf.data() + h_buf_factor_ptr(p - pbeg);
            if (n_m > 0) {
              UnmanagedViewType<value_type_matrix> ABR(bptr, n_m, n_m);
              bptr += ABR.span();
              UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // aptr += m*n_m;
              {
                _status = Trsm<Side::Left, Uplo::Upper, Trans::ConjTranspose, Algo::OnDevice>::invoke(
                    handle_blas, Diag::Unit(), one, ATL, ATR);
                checkDeviceBlasStatus("trsm");

                // Save ATR in workspace
                UnmanagedViewType<value_type_matrix> T(bptr, m, n_m); // bptr has been shifted
                _status = Copy<Algo::OnDevice>::invoke(exec_instance, T, ATR);

                // Apply D^{-1} on off-diagonal
                _status = Scale_BlockInverseDiagonals<Side::Left, Algo::OnDevice>::invoke(exec_instance, ATL, ATR);

                // ABR = -ATR*W
                _status = GemmTriangular<Trans::Transpose, Trans::NoTranspose, Uplo::Upper, Algo::OnDevice>::invoke(
                    handle_blas, minus_one, ATR, T, zero, ABR);
                checkDeviceBlasStatus("gemm");
              }
              {
                // additional things
                UnmanagedViewType<value_type_matrix> D(bptr, m, m);
                _status = Copy<Algo::OnDevice>::invoke(exec_instance, D, ATL);
                checkDeviceBlasStatus("Copy");

                _status = SetIdentity<Algo::OnDevice>::invoke(exec_instance, ATL, minus_one);
                checkDeviceBlasStatus("SetIdentity");

                UnmanagedViewType<value_type_matrix> AT(ATL.data(), m, n);
                _status = Trsm<Side::Left, Uplo::Upper, Trans::NoTranspose, Algo::OnDevice>::invoke(
                    handle_blas, Diag::Unit(), minus_one, D, AT);
                checkDeviceBlasStatus("trsm");

                // Copy original diagonal into T
                using policy_type = Kokkos::RangePolicy<exec_space>;
                const auto policy = policy_type(exec_instance, 0, m);
                Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const ordinal_type &i) {
                  ATL(i,i) = D(i,i);
                });
              }
            } else {
              /// additional things
              UnmanagedViewType<value_type_matrix> T(bptr, m, m);
              _status = Copy<Algo::OnDevice>::invoke(exec_instance, T, ATL);
              checkDeviceBlasStatus("Copy");

              _status = SetIdentity<Algo::OnDevice>::invoke(exec_instance, ATL, one);
              checkDeviceBlasStatus("SetIdentity");

              _status = Trsm<Side::Left, Uplo::Upper, Trans::NoTranspose, Algo::OnDevice>::invoke(
                  handle_blas, Diag::Unit(), one, T, ATL);

              // Copy original diagonal into T
              using policy_type = Kokkos::RangePolicy<exec_space>;
              const auto policy = policy_type(exec_instance, 0, m);
              Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const ordinal_type &i) {
                ATL(i,i) = T(i,i);
              });
              checkDeviceBlasStatus("trsm");
            }
            num_device_calls ++;
          }
        }
      }
    }
    return num_device_calls;
  }

  inline int factorizeNoPivotLDLOnDevice(const ordinal_type pbeg, const ordinal_type pend,
                                          const size_type_array_host &h_buf_factor_ptr, const value_type_array &work) {
    if (variant == 0)
      return factorizeNoPivotLDLOnDeviceVar0(pbeg, pend, h_buf_factor_ptr, work);
    else if (variant == 1)
      return factorizeNoPivotLDLOnDeviceVar1(pbeg, pend, h_buf_factor_ptr, work);
    else if (variant == 2 || variant == 3)
      return factorizeNoPivotLDLOnDeviceVar2(pbeg, pend, h_buf_factor_ptr, work);
    else {
      std::string msg = "Error: LevelSetTools::factorizeNoPivotLDLOnDevice, algorithm variant ("
                        + std::to_string(variant) + ") is not supported.\n";
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, msg.c_str());
    }
    return 0;
  }

  ///
  /// Cholesky
  ///
  inline int factorizeCholeskyOnDeviceVar0(const ordinal_type pbeg, const ordinal_type pend,
                                           const size_type_array_host &h_buf_factor_ptr,
                                           const value_type_array &work) {
    const value_type one(1), minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    int num_device_calls = 0;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_factorize_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        lapack_handle_type handle_lapack = getLapackHandle(qid);
        setStreamOnHandle(qid);

        const size_type worksize = work.extent(0) / _nstreams;
        value_type_array W(work.data() + worksize * qid, worksize);
        ++q;
#else
        lapack_handle_type handle_lapack = getLapackHandle();
        value_type_array W = work;
#endif
        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m;
          if (m > 0) {
            value_type *aptr = s.u_buf;
            UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            aptr += m * m;

            // On NVIDIA/AMD, calling Chol<ArgUplo, Algo::OnDevice>::invoke with member = handle.
            _status = Chol<Uplo::Upper, Algo::OnDevice>::invoke(handle_lapack, ATL, W);
            checkDeviceLapackStatus("chol");
          }
        }
      }
    }
    #if defined(KOKKOS_ENABLE_HIP)
    Kokkos::fence();
    #endif

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    q = 0;
#endif
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_factorize_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type handle_blas = getBlasHandle(qid);
        setStreamOnHandle(qid);

        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
#endif
        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            value_type *aptr = s.u_buf;
            UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            aptr += m * m;

            if (n_m > 0) {
              UnmanagedViewType<value_type_matrix> ABR(_buf.data() + h_buf_factor_ptr(p - pbeg), n_m, n_m);
              UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // aptr += m*n_m;

              _status = Trsm<Side::Left, Uplo::Upper, Trans::ConjTranspose, Algo::OnDevice>::invoke(
                  handle_blas, Diag::NonUnit(), one, ATL, ATR);
              checkDeviceBlasStatus("trsm");

              _status = Herk<Uplo::Upper, Trans::ConjTranspose, Algo::OnDevice>::invoke(handle_blas, minus_one, ATR,
                                                                                        zero, ABR);
            }
            num_device_calls ++;
          }
        }
      }
    }
    return num_device_calls;
  }

  inline int factorizeCholeskyOnDeviceVar1(const ordinal_type pbeg, const ordinal_type pend,
                                           const size_type_array_host &h_buf_factor_ptr,
                                           const value_type_array &work) {
    const value_type one(1), minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    int num_device_calls = 0;
    exec_space exec_instance;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_factorize_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        lapack_handle_type handle_lapack = getLapackHandle(qid);
        setStreamOnHandle(qid);

        const size_type worksize = work.extent(0) / _nstreams;
        value_type_array W(work.data() + worksize * qid, worksize);
        ++q;
#else
        lapack_handle_type handle_lapack = getLapackHandle();
        value_type_array W = work;
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m;
          if (m > 0) {
            // Factor the diagonal block
            value_type *aptr = s.u_buf;
            UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            aptr += m * m;

            // On NVIDIA/AMD, calling Chol<ArgUplo, Algo::OnDevice>::invoke with member = handle.
            _status = Chol<Uplo::Upper, Algo::OnDevice>::invoke(handle_lapack, ATL, W);
            checkDeviceLapackStatus("chol");
          }
        }
      }
    }
    #if defined(KOKKOS_ENABLE_HIP)
    Kokkos::fence();
    #endif

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    q = 0;
#endif
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_factorize_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type   handle_blas  = getBlasHandle(qid);
        setStreamOnHandle(qid);

        exec_instance = _exec_instances[qid];
        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            value_type *aptr = s.u_buf;
            UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            aptr += m * m;

            // Compute inverse of diagonal block (in T)
            value_type *bptr = _buf.data() + h_buf_factor_ptr(p - pbeg);
            UnmanagedViewType<value_type_matrix> T(bptr, m, m); // shared with ABR
            _status = SetIdentity<Algo::OnDevice>::invoke(exec_instance, T, one);
            checkDeviceBlasStatus("SetIdentity");

            _status = Trsm<Side::Left, Uplo::Upper, Trans::NoTranspose, Algo::OnDevice>::invoke(
                handle_blas, Diag::NonUnit(), one, ATL, T);
            checkDeviceBlasStatus("trsm");

            if (n_m > 0) {
              UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // aptr += m*n_m;
              UnmanagedViewType<value_type_matrix> ABR(bptr, n_m, n_m); // shared with T

              _status = Trsm<Side::Left, Uplo::Upper, Trans::ConjTranspose, Algo::OnDevice>::invoke(
                  handle_blas, Diag::NonUnit(), one, ATL, ATR);
              checkDeviceBlasStatus("trsm");

              // Copy T back to ATL (inverse)
              _status = Copy<Algo::OnDevice>::invoke(exec_instance, ATL, T);
              checkDeviceBlasStatus("Copy");

              // ABR = -ATR'*ATR
              _status = Herk<Uplo::Upper, Trans::ConjTranspose, Algo::OnDevice>::invoke(handle_blas, minus_one, ATR,
                                                                                        zero, ABR);
            } else {
              _status = Copy<Algo::OnDevice>::invoke(exec_instance, ATL, T);
              checkDeviceBlasStatus("Copy");
            }
            num_device_calls ++;
          }
        }
      }
    }
    return num_device_calls;
  }

  inline int factorizeCholeskyOnDeviceVar2(const ordinal_type pbeg, const ordinal_type pend,
                                           const size_type_array_host &h_buf_factor_ptr,
                                           const value_type_array &work) {
    const value_type one(1), minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    int num_device_calls = 0;
    exec_space exec_instance;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_factorize_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        lapack_handle_type handle_lapack = getLapackHandle(qid);
        setStreamOnHandle(qid);

        const size_type worksize = work.extent(0) / _nstreams;
        value_type_array W(work.data() + worksize * qid, worksize);
        ++q;
#else
        lapack_handle_type handle_lapack = getLapackHandle();
        value_type_array W = work;
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m;
          if (m > 0) {
            value_type *aptr = s.u_buf;
            UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            aptr += m * m;

            // On NVIDIA/AMD, calling Chol<ArgUplo, Algo::OnDevice>::invoke with member = handle.
            _status = Chol<Uplo::Upper, Algo::OnDevice>::invoke(handle_lapack, ATL, W);
            checkDeviceLapackStatus("chol");
          }
        }
      }
    }
    #if defined(KOKKOS_ENABLE_HIP)
    Kokkos::fence();
    #endif

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    q = 0;
#endif
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_factorize_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type   handle_blas = getBlasHandle(qid);
        setStreamOnHandle(qid);

        exec_instance = _exec_instances[qid];
        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            value_type *aptr = s.u_buf;
            UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            aptr += m * m;

            value_type *bptr = _buf.data() + h_buf_factor_ptr(p - pbeg);
            if (n_m > 0) {
              UnmanagedViewType<value_type_matrix> ABR(bptr, n_m, n_m);
              bptr += ABR.span();
              UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // aptr += m*n_m;

              _status = Trsm<Side::Left, Uplo::Upper, Trans::ConjTranspose, Algo::OnDevice>::invoke(
                  handle_blas, Diag::NonUnit(), one, ATL, ATR);
              checkDeviceBlasStatus("trsm");

              _status = Herk<Uplo::Upper, Trans::ConjTranspose, Algo::OnDevice>::invoke(handle_blas, minus_one, ATR,
                                                                                        zero, ABR);

              /// additional things
              UnmanagedViewType<value_type_matrix> T(bptr, m, m);
              _status = Copy<Algo::OnDevice>::invoke(exec_instance, T, ATL);
              checkDeviceBlasStatus("Copy");

              _status = SetIdentity<Algo::OnDevice>::invoke(exec_instance, ATL, minus_one);
              checkDeviceBlasStatus("SetIdentity");

              UnmanagedViewType<value_type_matrix> AT(ATL.data(), m, n);
              _status = Trsm<Side::Left, Uplo::Upper, Trans::NoTranspose, Algo::OnDevice>::invoke(
                  handle_blas, Diag::NonUnit(), minus_one, T, AT);
              checkDeviceBlasStatus("trsm");
            } else {
              /// additional things
              UnmanagedViewType<value_type_matrix> T(bptr, m, m);
              _status = Copy<Algo::OnDevice>::invoke(exec_instance, T, ATL);
              checkDeviceBlasStatus("Copy");

              _status = SetIdentity<Algo::OnDevice>::invoke(exec_instance, ATL, one);
              checkDeviceBlasStatus("SetIdentity");

              _status = Trsm<Side::Left, Uplo::Upper, Trans::NoTranspose, Algo::OnDevice>::invoke(
                  handle_blas, Diag::NonUnit(), one, T, ATL);
              checkDeviceBlasStatus("trsm");
            }
            num_device_calls ++;
          }
        }
      }
    }
    return num_device_calls;
  }

  inline int factorizeCholeskyOnDevice(const ordinal_type pbeg, const ordinal_type pend,
                                       const size_type_array_host &h_buf_factor_ptr, const value_type_array &work) {
    if (variant == 0)
      return factorizeCholeskyOnDeviceVar0(pbeg, pend, h_buf_factor_ptr, work);
    else if (variant == 1)
      return factorizeCholeskyOnDeviceVar1(pbeg, pend, h_buf_factor_ptr, work);
    else if (variant == 2 || variant == 3)
      return factorizeCholeskyOnDeviceVar2(pbeg, pend, h_buf_factor_ptr, work);
    else {
      std::string msg = "Error: LevelSetTools::factorizeCholeskyOnDevice, algorithm variant ("
                        + std::to_string(variant) + ") is not supported.\n";
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, msg.c_str());
    }
    return 0;
  }

  /// 
  /// LDL
  ///
  inline void factorizeLDL_OnDeviceVar0(const ordinal_type pbeg, const ordinal_type pend,
                                        const size_type_array_host &h_buf_factor_ptr, const value_type_array &work) {
    const value_type one(1), minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    exec_space exec_instance;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_factorize_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type   handle_blas   = getBlasHandle(qid);
        lapack_handle_type handle_lapack = getLapackHandle(qid);

        setStreamOnHandle(qid);
        exec_instance = _exec_instances[qid];

        const size_type worksize = work.extent(0) / _nstreams;
        value_type_array W(work.data() + worksize * qid, worksize);
        ++q;
#else
        blas_handle_type   handle_blas   = getBlasHandle();
        lapack_handle_type handle_lapack = getLapackHandle();
        value_type_array W = work;
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type offs = s.row_begin, m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            value_type *aptr = s.u_buf;
            UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            aptr += m * m;

            _status = Symmetrize<Uplo::Upper, Algo::OnDevice>::invoke(exec_instance, ATL);

            ordinal_type *pivptr = _piv.data() + 4 * offs;
            UnmanagedViewType<ordinal_type_array> P(pivptr, 4 * m);
            _status = LDL<Uplo::Lower, Algo::OnDevice>::invoke(handle_lapack, ATL, P, W);
            checkDeviceLapackStatus("ldl::invoke");

            value_type *dptr = _diag.data() + 2 * offs;
            UnmanagedViewType<value_type_matrix> D(dptr, m, 2);
            _status = LDL<Uplo::Lower, Algo::OnDevice>::modify(exec_instance, ATL, P, D);
            checkDeviceLapackStatus("ldl::modify");

            if (n_m > 0) {
              UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // aptr += m*n_m;
              UnmanagedViewType<value_type_matrix> ABR(_buf.data() + h_buf_factor_ptr(p - pbeg), n_m, n_m);
              UnmanagedViewType<value_type_matrix> STR(ABR.data() + ABR.span(), m, n_m);

              auto fpiv = ordinal_type_array(P.data() + m, m);
              _status = ApplyPivots<PivotMode::Flame, Side::Left, Direct::Forward, Algo::OnDevice>::invoke(
                  exec_instance, fpiv, ATR);

              _status = Trsm<Side::Left, Uplo::Lower, Trans::NoTranspose, Algo::OnDevice>::invoke(
                  handle_blas, Diag::Unit(), one, ATL, ATR);
              checkDeviceBlasStatus("trsm");

              _status = Copy<Algo::OnDevice>::invoke(exec_instance, STR, ATR);

              _status = Scale2x2_BlockInverseDiagonals<Side::Left, Algo::OnDevice>::invoke(exec_instance, P, D, ATR);

              _status = GemmTriangular<Trans::Transpose, Trans::NoTranspose, Uplo::Upper, Algo::OnDevice>::invoke(
                  handle_blas, minus_one, ATR, STR, zero, ABR);
              checkDeviceBlasStatus("gemm");
            }
          }
        }
      }
    }
  }

  inline void factorizeLDL_OnDeviceVar1(const ordinal_type pbeg, const ordinal_type pend,
                                        const size_type_array_host &h_buf_factor_ptr, const value_type_array &work) {
    const value_type one(1), minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    exec_space exec_instance;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_factorize_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type   handle_blas   = getBlasHandle(qid);
        lapack_handle_type handle_lapack = getLapackHandle(qid);

        setStreamOnHandle(qid);
        exec_instance = _exec_instances[qid];

        const size_type worksize = work.extent(0) / _nstreams;
        value_type_array W(work.data() + worksize * qid, worksize);
        ++q;
#else
        blas_handle_type   handle_blas   = getBlasHandle();
        lapack_handle_type handle_lapack = getLapackHandle();
        value_type_array W = work;
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type offs = s.row_begin, m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            value_type *aptr = s.u_buf;
            UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            aptr += m * m;

            _status = Symmetrize<Uplo::Upper, Algo::OnDevice>::invoke(exec_instance, ATL);

            ordinal_type *pivptr = _piv.data() + 4 * offs;
            UnmanagedViewType<ordinal_type_array> P(pivptr, 4 * m);
            _status = LDL<Uplo::Lower, Algo::OnDevice>::invoke(handle_lapack, ATL, P, W);
            checkDeviceLapackStatus("ldl::invoke");

            value_type *dptr = _diag.data() + 2 * offs;
            UnmanagedViewType<value_type_matrix> D(dptr, m, 2);
            _status = LDL<Uplo::Lower, Algo::OnDevice>::modify(exec_instance, ATL, P, D);
            checkDeviceLapackStatus("ldl::modify");

            value_type *bptr = _buf.data() + h_buf_factor_ptr(p - pbeg);
            UnmanagedViewType<value_type_matrix> T(bptr, m, m);

            _status = SetIdentity<Algo::OnDevice>::invoke(exec_instance, T, one);
            checkDeviceBlasStatus("SetIdentity");

            _status = Trsm<Side::Left, Uplo::Lower, Trans::NoTranspose, Algo::OnDevice>::invoke(
                handle_blas, Diag::Unit(), one, ATL, T);
            checkDeviceBlasStatus("trsm");

            if (n_m > 0) {
              UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // aptr += m*n_m;
              UnmanagedViewType<value_type_matrix> ABR(bptr, n_m, n_m);

              const ordinal_type used_span = max(ABR.span(), T.span());
              UnmanagedViewType<value_type_matrix> STR(ABR.data() + used_span, m, n_m);

              ConstUnmanagedViewType<ordinal_type_array> perm(P.data() + 2 * m, m);
              _status = ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(exec_instance, ATR,
                                                                                                 perm, STR);

              _status = Trsm<Side::Left, Uplo::Lower, Trans::NoTranspose, Algo::OnDevice>::invoke(
                  handle_blas, Diag::Unit(), one, ATL, STR);
              checkDeviceBlasStatus("trsm");

              _status = Copy<Algo::OnDevice>::invoke(exec_instance, ATL, T);
              _status = Copy<Algo::OnDevice>::invoke(exec_instance, ATR, STR);

              _status = Scale2x2_BlockInverseDiagonals<Side::Left, Algo::OnDevice>::invoke(exec_instance, P, D, ATR);

              _status = GemmTriangular<Trans::Transpose, Trans::NoTranspose, Uplo::Upper, Algo::OnDevice>::invoke(
                  handle_blas, minus_one, ATR, STR, zero, ABR);
              checkDeviceBlasStatus("gemm");
            } else {
              _status = Copy<Algo::OnDevice>::invoke(exec_instance, ATL, T);
            }
          }
        }
      }
    }
  }

  inline void factorizeLDL_OnDeviceVar2(const ordinal_type pbeg, const ordinal_type pend,
                                        const size_type_array_host &h_buf_factor_ptr, const value_type_array &work) {
    const value_type one(1), minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    exec_space exec_instance;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_factorize_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type   handle_blas   = getBlasHandle(qid);
        lapack_handle_type handle_lapack = getLapackHandle(qid);

        setStreamOnHandle(qid);
        exec_instance = _exec_instances[qid];

        const size_type worksize = work.extent(0) / _nstreams;
        value_type_array W(work.data() + worksize * qid, worksize);
        ++q;
#else
        blas_handle_type   handle_blas   = getBlasHandle();
        lapack_handle_type handle_lapack = getLapackHandle();
        value_type_array W = work;
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type offs = s.row_begin, m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            value_type *aptr = s.u_buf;
            UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            aptr += m * m;

            _status = Symmetrize<Uplo::Upper, Algo::OnDevice>::invoke(exec_instance, ATL);

            ordinal_type *pivptr = _piv.data() + 4 * offs;
            UnmanagedViewType<ordinal_type_array> P(pivptr, 4 * m);
            _status = LDL<Uplo::Lower, Algo::OnDevice>::invoke(handle_lapack, ATL, P, W);
            checkDeviceLapackStatus("ldl::invoke");

            value_type *dptr = _diag.data() + 2 * offs;
            UnmanagedViewType<value_type_matrix> D(dptr, m, 2);
            _status = LDL<Uplo::Lower, Algo::OnDevice>::modify(exec_instance, ATL, P, D);
            checkDeviceLapackStatus("ldl::modify");

            value_type *bptr = _buf.data() + h_buf_factor_ptr(p - pbeg);

            if (n_m > 0) {
              UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // aptr += m*n_m;
              UnmanagedViewType<value_type_matrix> ABR(bptr, n_m, n_m);
              UnmanagedViewType<value_type_matrix> T(bptr + ABR.span(), m, m);

              const ordinal_type used_span = ABR.span() + T.span();
              UnmanagedViewType<value_type_matrix> STR(bptr + used_span, m, n_m);

              // STR := L^{-1} * ATR
              ConstUnmanagedViewType<ordinal_type_array> perm(P.data() + 2 * m, m);
              _status = ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(exec_instance, ATR,
                                                                                                 perm, STR);

              _status = Trsm<Side::Left, Uplo::Lower, Trans::NoTranspose, Algo::OnDevice>::invoke(
                  handle_blas, Diag::Unit(), one, ATL, STR);
              checkDeviceBlasStatus("trsm");

              // ATR := (LT)^{-1} * ATR
              _status = Copy<Algo::OnDevice>::invoke(exec_instance, ATR, STR);
              _status = Scale2x2_BlockInverseDiagonals<Side::Left, Algo::OnDevice>::invoke(exec_instance, P, D, ATR);

              // ABR := ATR^T * STR = ((LT)^{-1}*ATR)^T * (L^{-1}*ATR)
              _status = GemmTriangular<Trans::Transpose, Trans::NoTranspose, Uplo::Upper, Algo::OnDevice>::invoke(
                  handle_blas, minus_one, ATR, STR, zero, ABR);
              checkDeviceBlasStatus("gemm");

              // AT = ATL^{-1} [I, ATR] (= L^{-1} where A = LTL^T and A^{-1} = L^{-T} T^{-1} L^{-1} = (TL^{-T})^{-1) * L^{-1})
              //                                                             = (solveLDL_Upper_varian2 with Scale2x2_BlockInverseDiagonals) * (solveLDL_Lower_variant2)
              _status = Copy<Algo::OnDevice>::invoke(exec_instance, T, ATL);
              _status = Symmetrize<Uplo::Lower, Algo::OnDevice>::invoke(exec_instance, T);
              _status = SetIdentity<Algo::OnDevice>::invoke(exec_instance, ATL, minus_one);

              UnmanagedViewType<value_type_matrix> AT(ATL.data(), m, n);
              _status = Trsm<Side::Left, Uplo::Upper, Trans::NoTranspose, Algo::OnDevice>::invoke(
                  handle_blas, Diag::Unit(), minus_one, T, AT);
            } else {
              UnmanagedViewType<value_type_matrix> T(bptr, m, m);
              _status = Copy<Algo::OnDevice>::invoke(exec_instance, T, ATL);

              _status = SetIdentity<Algo::OnDevice>::invoke(exec_instance, ATL, one);

              _status = Trsm<Side::Left, Uplo::Lower, Trans::NoTranspose, Algo::OnDevice>::invoke(
                  handle_blas, Diag::Unit(), one, T, ATL);
            }
          }
        }
      }
    }
  }

  inline void factorizeLDL_OnDevice(const ordinal_type pbeg, const ordinal_type pend,
                                    const size_type_array_host &h_buf_factor_ptr, const value_type_array &work) {
    if (variant == 0)
      factorizeLDL_OnDeviceVar0(pbeg, pend, h_buf_factor_ptr, work);
    else if (variant == 1)
      factorizeLDL_OnDeviceVar1(pbeg, pend, h_buf_factor_ptr, work);
    else if (variant == 2 || variant == 3)
      factorizeLDL_OnDeviceVar2(pbeg, pend, h_buf_factor_ptr, work);
    else {
      std::string msg = "Error: LevelSetTools::factorizeLDL_OnDevice, algorithm variant ("
                        + std::to_string(variant) + ") is not supported.\n";
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, msg.c_str());
    }
  }

  /// 
  /// LU
  ///
  inline int factorizeLU_OnDeviceVar0(const ordinal_type pbeg, const ordinal_type pend,
                                      const size_type_array_host &h_buf_factor_ptr, const value_type_array &work) {
    const value_type one(1), minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    int num_device_calls = 0;
    exec_space exec_instance;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_factorize_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type   handle_blas   = getBlasHandle(qid);
        lapack_handle_type handle_lapack = getLapackHandle(qid);
        setStreamOnHandle(qid);
        exec_instance = _exec_instances[qid];

        const size_type worksize = work.extent(0) / _nstreams;
        value_type_array W(work.data() + worksize * qid, worksize);
        ++q;
#else
        blas_handle_type   handle_blas   = getBlasHandle();
        lapack_handle_type handle_lapack = getLapackHandle();
        value_type_array W = work;
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type offs = s.row_begin, m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            value_type *uptr = s.u_buf;
            UnmanagedViewType<value_type_matrix> AT(uptr, m, n);

            ordinal_type *pivptr = _piv.data() + 4 * offs;
            UnmanagedViewType<ordinal_type_array> P(pivptr, 4 * m);
            _status = LU<Algo::OnDevice>::invoke(handle_lapack, AT, P, W);
            checkDeviceLapackStatus("lu::invoke");

            _status = LU<Algo::OnDevice>::modify(exec_instance, m, P);
            checkDeviceLapackStatus("lu::modify");

            if (n_m > 0) {
              UnmanagedViewType<value_type_matrix> ATL(uptr, m, m);
              uptr += m * m;
              UnmanagedViewType<value_type_matrix> ATR(uptr, m, n_m);

              UnmanagedViewType<value_type_matrix> AL(s.l_buf, n, m);
              const auto ABL = Kokkos::subview(AL, range_type(m, n), Kokkos::ALL());
              UnmanagedViewType<value_type_matrix> ABR(_buf.data() + h_buf_factor_ptr(p - pbeg), n_m, n_m);

              _status = Trsm<Side::Right, Uplo::Upper, Trans::NoTranspose, Algo::OnDevice>::invoke(
                  handle_blas, Diag::NonUnit(), one, ATL, ABL);
              checkDeviceBlasStatus("trsm");

              _status = Gemm<Trans::NoTranspose, Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, minus_one,
                                                                                             ABL, ATR, zero, ABR);
              checkDeviceBlasStatus("gemm");
            }
            num_device_calls ++;
          }
        }
      }
    }
    return num_device_calls;
  }

  inline int factorizeLU_OnDeviceVar1(const ordinal_type pbeg, const ordinal_type pend,
                                      const size_type_array_host &h_buf_factor_ptr, const value_type_array &work) {
    const value_type one(1), minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    int num_device_calls = 0;
    exec_space exec_instance;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_factorize_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type   handle_blas   = getBlasHandle(qid);
        lapack_handle_type handle_lapack = getLapackHandle(qid);

        setStreamOnHandle(qid);
        exec_instance = _exec_instances[qid];

        const size_type worksize = work.extent(0) / _nstreams;
        value_type_array W(work.data() + worksize * qid, worksize);
        ++q;
#else
        blas_handle_type   handle_blas   = getBlasHandle();
        lapack_handle_type handle_lapack = getLapackHandle();
        value_type_array W = work;
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type offs = s.row_begin, m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            UnmanagedViewType<value_type_matrix> AT(s.u_buf, m, n);

            ordinal_type *pivptr = _piv.data() + 4 * offs;
            UnmanagedViewType<ordinal_type_array> P(pivptr, 4 * m);
            _status = LU<Algo::OnDevice>::invoke(handle_lapack, AT, P, W);
            checkDeviceLapackStatus("lu::invoke");

            _status = LU<Algo::OnDevice>::modify(exec_instance, m, P);
            checkDeviceLapackStatus("lu::modify");

            value_type *bptr = _buf.data() + h_buf_factor_ptr(p - pbeg);
            UnmanagedViewType<value_type_matrix> T(bptr, m, m);

            if (n_m > 0) {
              UnmanagedViewType<value_type_matrix> ATL(s.u_buf, m, m);
              UnmanagedViewType<value_type_matrix> ATR(s.u_buf + ATL.span(), m, n_m);

              UnmanagedViewType<value_type_matrix> AL(s.l_buf, n, m);
              const auto ATL2 = Kokkos::subview(AL, range_type(0, m), Kokkos::ALL());
              const auto ABL = Kokkos::subview(AL, range_type(m, n), Kokkos::ALL());

              UnmanagedViewType<value_type_matrix> ABR(bptr, n_m, n_m);

              _status = Copy<Algo::OnDevice>::invoke(exec_instance, T, ATL);
              _status = Trsm<Side::Right, Uplo::Upper, Trans::NoTranspose, Algo::OnDevice>::invoke(
                  handle_blas, Diag::NonUnit(), one, ATL, ABL);
              checkDeviceBlasStatus("trsm");

              _status = SetIdentity<Algo::OnDevice>::invoke(exec_instance, ATL, one);
              _status = SetIdentity<Algo::OnDevice>::invoke(exec_instance, ATL2, one);

              _status = Trsm<Side::Left, Uplo::Upper, Trans::NoTranspose, Algo::OnDevice>::invoke(
                  handle_blas, Diag::NonUnit(), one, T, ATL);
              checkDeviceBlasStatus("trsm");
              _status = Trsm<Side::Left, Uplo::Lower, Trans::NoTranspose, Algo::OnDevice>::invoke(
                  handle_blas, Diag::Unit(), one, T, ATL2);
              checkDeviceBlasStatus("trsm");

              _status = Gemm<Trans::NoTranspose, Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, minus_one,
                                                                                             ABL, ATR, zero, ABR);
              checkDeviceBlasStatus("gemm");
            } else {
              UnmanagedViewType<value_type_matrix> ATL(s.u_buf, m, m);
              UnmanagedViewType<value_type_matrix> ATL2(s.l_buf, m, m);

              _status = Copy<Algo::OnDevice>::invoke(exec_instance, T, ATL);

              _status = SetIdentity<Algo::OnDevice>::invoke(exec_instance, ATL, one);
              _status = SetIdentity<Algo::OnDevice>::invoke(exec_instance, ATL2, one);

              _status = Trsm<Side::Left, Uplo::Upper, Trans::NoTranspose, Algo::OnDevice>::invoke(
                  handle_blas, Diag::NonUnit(), one, T, ATL);
              checkDeviceBlasStatus("trsm");
              _status = Trsm<Side::Left, Uplo::Lower, Trans::NoTranspose, Algo::OnDevice>::invoke(
                  handle_blas, Diag::Unit(), one, T, ATL2);
              checkDeviceBlasStatus("trsm");
            }
            num_device_calls ++;
          }
        }
      }
    }
    return num_device_calls;
  }

  inline int factorizeLU_OnDeviceVar2(const ordinal_type pbeg, const ordinal_type pend,
                                      const size_type_array_host &h_buf_factor_ptr, const value_type_array &work) {
    const value_type one(1), minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    int num_device_calls = 0;
    exec_space exec_instance;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_factorize_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type   handle_blas   = getBlasHandle(qid);
        lapack_handle_type handle_lapack = getLapackHandle(qid);
        setStreamOnHandle(qid);
        exec_instance = _exec_instances[qid];

        const size_type worksize = work.extent(0) / _nstreams;
        value_type_array W(work.data() + worksize * qid, worksize);
        ++q;
#else
        blas_handle_type   handle_blas   = getBlasHandle();
        lapack_handle_type handle_lapack = getLapackHandle();
        value_type_array W = work;
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type offs = s.row_begin, m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            UnmanagedViewType<value_type_matrix> AT(s.u_buf, m, n);

            ordinal_type *pivptr = _piv.data() + 4 * offs;
            UnmanagedViewType<ordinal_type_array> P(pivptr, 4 * m);
            _status = LU<Algo::OnDevice>::invoke(handle_lapack, AT, P, W);
            checkDeviceLapackStatus("lu::invoke");

            _status = LU<Algo::OnDevice>::modify(exec_instance, m, P);
            checkDeviceLapackStatus("lu::modify");

            value_type *bptr = _buf.data() + h_buf_factor_ptr(p - pbeg);

            if (n_m > 0) {
              UnmanagedViewType<value_type_matrix> ATL(s.u_buf, m, m);
              UnmanagedViewType<value_type_matrix> ATR(s.u_buf + ATL.span(), m, n_m);

              UnmanagedViewType<value_type_matrix> AL(s.l_buf, n, m);
              const auto ATL2 = Kokkos::subview(AL, range_type(0, m), Kokkos::ALL());
              const auto ABL = Kokkos::subview(AL, range_type(m, n), Kokkos::ALL());

              UnmanagedViewType<value_type_matrix> ABR(bptr, n_m, n_m);
              UnmanagedViewType<value_type_matrix> T(bptr + ABR.span(), m, m);

              _status = Copy<Algo::OnDevice>::invoke(exec_instance, T, ATL);
              _status = Trsm<Side::Right, Uplo::Upper, Trans::NoTranspose, Algo::OnDevice>::invoke(
                  handle_blas, Diag::NonUnit(), one, ATL, ABL);
              checkDeviceBlasStatus("trsm");

              _status = SetIdentity<Algo::OnDevice>::invoke(exec_instance, ATL, minus_one);
              _status = SetIdentity<Algo::OnDevice>::invoke(exec_instance, ATL2, minus_one);

              _status = Gemm<Trans::NoTranspose, Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, minus_one,
                                                                                             ABL, ATR, zero, ABR);
              checkDeviceBlasStatus("gemm");

              _status = Trsm<Side::Left, Uplo::Upper, Trans::NoTranspose, Algo::OnDevice>::invoke(
                  handle_blas, Diag::NonUnit(), minus_one, T, AT);
              checkDeviceBlasStatus("trsm");
              _status = Trsm<Side::Right, Uplo::Lower, Trans::NoTranspose, Algo::OnDevice>::invoke(
                  handle_blas, Diag::Unit(), minus_one, T, AL);
              checkDeviceBlasStatus("trsm");

            } else {
              UnmanagedViewType<value_type_matrix> ATL(s.u_buf, m, m);
              UnmanagedViewType<value_type_matrix> ATL2(s.l_buf, m, m);
              UnmanagedViewType<value_type_matrix> T(bptr, m, m);

              _status = Copy<Algo::OnDevice>::invoke(exec_instance, T, ATL);

              _status = SetIdentity<Algo::OnDevice>::invoke(exec_instance, ATL, one);
              _status = SetIdentity<Algo::OnDevice>::invoke(exec_instance, ATL2, one);

              _status = Trsm<Side::Left, Uplo::Upper, Trans::NoTranspose, Algo::OnDevice>::invoke(
                  handle_blas, Diag::NonUnit(), one, T, ATL);
              checkDeviceBlasStatus("trsm");
              _status = Trsm<Side::Left, Uplo::Lower, Trans::NoTranspose, Algo::OnDevice>::invoke(
                  handle_blas, Diag::Unit(), one, T, ATL2);
              checkDeviceBlasStatus("trsm");
            }
            num_device_calls ++;
          }
        }
      }
    }
    return num_device_calls;
  }

  inline int factorizeLU_OnDevice(const ordinal_type pbeg, const ordinal_type pend,
                                  const size_type_array_host &h_buf_factor_ptr, const value_type_array &work) {
    if (variant == 0)
      return factorizeLU_OnDeviceVar0(pbeg, pend, h_buf_factor_ptr, work);
    else if (variant == 1)
      return factorizeLU_OnDeviceVar1(pbeg, pend, h_buf_factor_ptr, work);
    else if (variant == 2 || variant == 3)
      return factorizeLU_OnDeviceVar2(pbeg, pend, h_buf_factor_ptr, work);
    else {
      std::string msg = "Error: LevelSetTools::factorizeLU_OnDevice, algorithm variant ("
                        + std::to_string(variant) + ") is not supported.\n";
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, msg.c_str());
    }
    return 0;
  }

  ///
  /// Extract CRS for SpMV
  ///
  inline void setupCRS(bool store_transpose, bool verbose) {
    const int method = this->getSolutionMethod();
    _spmv->Setup(store_transpose, verbose, method, _m, _team_serial_level_cut,
                _h_level_ptr, _level_sids, _h_level_sids, _h_solve_mode, _info, _h_supernodes, _solve_mode, _piv,
                _streams[0],
                _w_vec);
  }

  inline void loadCRS(bool store_transpose, bool verbose) {
    const int method = this->getSolutionMethod();
    _spmv->Load(store_transpose, verbose, method, _m,
                _h_level_ptr, _level_sids, _h_level_sids, _h_solve_mode, _info, _h_supernodes, _solve_mode, _piv,
                _streams[0], _w_vec);
  }

  inline void extractCRS(bool store_transpose, bool verbose) {
    if (verbose) {
      printf("LevelSetTools:extractCRS\n");
      printf("========================\n");
      if (store_transpose) printf( "Store Transpose\n" );
    }
    if (!_keep_zeros) {
      setupCRS(store_transpose, verbose);
    }
    loadCRS(store_transpose, verbose);
    if (verbose) {
      printf("\n");
    }
  }

  /// 
  /// Release CRS extracted for SpMV
  ///
  inline void releaseCRS(bool release_all, bool verbose) {
    _spmv->Release(release_all, verbose,
                  _h_level_ptr, _h_level_sids, _h_supernodes);
  }


  /// 
  /// Functions for Solve
  ///
  inline void allocateWorkspaceSolve(const ordinal_type nrhs) {
    if (variant == 3) {
    } else {
      const size_type buf_extent = _bufsize_solve * nrhs;
      const size_type buf_span = _buf.span();

      if (buf_extent > buf_span) {
        if (_buf.span() > 0) track_free(buf_span * sizeof(value_type));
        Kokkos::resize(_buf, buf_extent);
        track_alloc(_buf.span() * sizeof(value_type));
      }
      if (nrhs > _nrhs) {
        // update **pointer** to solver-workspace with differet nrhs
        const Kokkos::RangePolicy<exec_space> policy(0, _buf_solve_ptr.extent(0));
        const auto buf_solve_nrhs_ptr = _buf_solve_nrhs_ptr;
        const auto buf_solve_ptr = _buf_solve_ptr;
        Kokkos::parallel_for(
            policy, KOKKOS_LAMBDA(const ordinal_type &i) { buf_solve_nrhs_ptr(i) = nrhs * buf_solve_ptr(i); });
        Kokkos::deep_copy(_h_buf_solve_nrhs_ptr, _buf_solve_nrhs_ptr);
        _nrhs = nrhs;
      }
    }
  }

  /// 
  /// Non-pivot LDL Lower
  ///
  inline int solveNoPivotLDLLowerOnDeviceVar0(const ordinal_type pbeg, const ordinal_type pend,
                                              const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const ordinal_type nrhs = t.extent(1);
    const value_type minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    int num_device_calls = 0;
    exec_space exec_instance;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type handle_blas = getBlasHandle(qid);
        setStreamOnHandle(qid);
        exec_instance = _exec_instances[qid];
        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            value_type *aptr = s.u_buf;
            UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            aptr += ATL.span();

            const ordinal_type offm = s.row_begin;
            auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());

            _status =
                Trsv<Uplo::Upper, Trans::ConjTranspose, Algo::OnDevice>::invoke(handle_blas, Diag::Unit(), ATL, tT);
            checkDeviceBlasStatus("trsv");

            if (n_m > 0) {
              // solve offdiag
              value_type *bptr = _buf.data() + h_buf_solve_ptr(p - pbeg);
              UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // aptr += m*n_m;
              UnmanagedViewType<value_type_matrix> bB(bptr, n_m, nrhs);

              _status = Gemv<Trans::ConjTranspose, Algo::OnDevice>::invoke(handle_blas, minus_one, ATR, tT, zero, bB);
              checkDeviceBlasStatus("gemv");
            }

            // Apply D^{-1} on off-diagonal
            _status = Scale_BlockInverseDiagonals<Side::Left, Algo::OnDevice>::invoke(exec_instance, ATL, tT);
            checkDeviceBlasStatus("scale");

            // increment num device calls
            num_device_calls ++;
          }
        }
      }
    }
    return num_device_calls;
  }

  inline int solveNoPivotLDLLowerOnDeviceVar1(const ordinal_type pbeg, const ordinal_type pend,
                                              const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const ordinal_type nrhs = t.extent(1);
    const value_type minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    int num_device_calls = 0;
    exec_space exec_instance;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type handle_blas = getBlasHandle(qid);
        setStreamOnHandle(qid);
        exec_instance = _exec_instances[qid];
        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            value_type *aptr = s.u_buf;
            UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            aptr += m * m;

            value_type *bptr = _buf.data() + h_buf_solve_ptr(p - pbeg);
            UnmanagedViewType<value_type_matrix> bT(bptr, m, nrhs);
            bptr += m * nrhs;

            const ordinal_type offm = s.row_begin;
            const auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());

             // Copy current block (TODO: should be in Trmv?)
            _status = Copy<Algo::OnDevice>::invoke(exec_instance, bT, tT);
            checkDeviceBlasStatus("copy");

            // Solve diag (ATL is square)
            _status = Trmv<Uplo::Upper, Trans::ConjTranspose, Algo::OnDevice>::invoke(handle_blas, Diag::Unit(), ATL, tT, bT);
            checkDeviceBlasStatus("gemv");

            if (n_m > 0) {
              // Update offdiag
              UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m);
              UnmanagedViewType<value_type_matrix> bB(bptr, n_m, nrhs);

              _status = Gemv<Trans::ConjTranspose, Algo::OnDevice>::invoke(handle_blas, minus_one, ATR, bT, zero, bB);
              checkDeviceBlasStatus("gemv");
            }

            // Apply D^{-1} on off-diagonal
            _status = Scale_BlockInverseDiagonals<Side::Left, Algo::OnDevice>::invoke(exec_instance, ATL, bT);
            checkDeviceBlasStatus("scale");

            // increment num device calls
            num_device_calls ++;
          }
        }
      }
    }
    return num_device_calls;
  }

  inline int solveNoPivotLDLLowerOnDeviceVar2(const ordinal_type pbeg, const ordinal_type pend,
                                              const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const ordinal_type nrhs = t.extent(1);

    exec_space exec_instance;
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    int num_device_calls = 0;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type handle_blas = getBlasHandle(qid);
        setStreamOnHandle(qid);
        exec_instance = _exec_instances[qid];
        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n;
          if (m > 0) {
            value_type *aptr = s.u_buf;
            UnmanagedViewType<value_type_matrix> AT(aptr, m, n);

            value_type *bptr = _buf.data() + h_buf_solve_ptr(p - pbeg);
            UnmanagedViewType<value_type_matrix> b(bptr, n, nrhs);

            const ordinal_type offm = s.row_begin;
            auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());

            // Copy current block from t to buffer (NOTE: should be in Trmv?)
            const auto dx = Kokkos::subview(b,  range_type(0, m), Kokkos::ALL());
            _status = Copy<Algo::OnDevice>::invoke(exec_instance, dx, tT);
            checkDeviceBlasStatus("copy");

            // AT is short-wide: b := AT'\t : compute current block, and then use it to update off-diagonal
            _status = Trmv<Uplo::Upper, Trans::ConjTranspose, Algo::OnDevice>::invoke(handle_blas, Diag::Unit(), AT, tT, b);
            checkDeviceBlasStatus("trmv");

            // Apply D^{-1} on diagonal (already updated using this block, also n >= m so diagonal block is stored first in aptr)
            UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            _status = Scale_BlockInverseDiagonals<Side::Left, Algo::OnDevice>::invoke(exec_instance, ATL, dx);
            checkDeviceBlasStatus("scale");

            // increment num device calls
            num_device_calls ++;
          }
        }
      }
    }
    return num_device_calls;
  }

  inline int solveNoPivotLDLLowerOnDevice(const ordinal_type lvl, const ordinal_type nlvls,
                                          const ordinal_type pbeg, const ordinal_type pend,
                                          const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    if (variant == 0)
      return solveNoPivotLDLLowerOnDeviceVar0(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 1)
      return solveNoPivotLDLLowerOnDeviceVar1(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 2)
      return solveNoPivotLDLLowerOnDeviceVar2(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 3) {
      // apply L^{-1}
      auto &s0 = _h_supernodes(_h_level_sids(pbeg));
      _spmv->ApplyL_OnDevice(lvl, s0, t, _w_vec);
      {
        // apply D^{-1}
        const ordinal_type nrhs = t.extent(1);
        auto matY = (lvl == 0 ? t : ((nlvls-1-lvl)%2 == 0 ? _w_vec : t));
        auto &s0 = _h_supernodes(_h_level_sids(pbeg));
        const UnmanagedViewType<nzvals_view> matD(s0.nzvalsD, _m);

        using policy_type = Kokkos::RangePolicy<exec_space>;
        const auto policy = policy_type(_exec_instances[0], 0, _m);
        Kokkos::parallel_for(
            policy, KOKKOS_LAMBDA(const ordinal_type &i) {
              for (ordinal_type j = 0; j < nrhs; j++) matY(i,j) = matY(i, j) / matD(i); 
            });
      }
    } else {
      std::string msg = "Error: LevelSetTools::solveNoPivotLDLLowerOnDevice, algorithm variant ("
                        + std::to_string(variant) + "is not supported.\n";
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, msg.c_str());
    }
    return 0;
  }

  /// 
  /// Non-pivot LDL Upper
  ///
  inline int solveNoPivotLDLUpperOnDeviceVar0(const ordinal_type pbeg, const ordinal_type pend,
                                              const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const ordinal_type nrhs = t.extent(1);
    const value_type minus_one(-1), one(1);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    int num_device_calls = 0;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type handle_blas = getBlasHandle(qid);
        setStreamOnHandle(qid);
        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            value_type *aptr = s.u_buf, *bptr = _buf.data() + h_buf_solve_ptr(p - pbeg);
            ;
            const UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            aptr += m * m;
            const UnmanagedViewType<value_type_matrix> bB(bptr, n_m, nrhs);

            const ordinal_type offm = s.row_begin;
            const auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());

            if (n_m > 0) {
              const UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // aptr += m*n;
              _status = Gemv<Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, minus_one, ATR, bB, one, tT);
              checkDeviceBlasStatus("gemv");
            }
            _status =
                Trsv<Uplo::Upper, Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, Diag::Unit(), ATL, tT);
            checkDeviceBlasStatus("trsv");

            // increment num device calls
            num_device_calls ++;
          }
        }
      }
    }
    return num_device_calls;
  }

  inline int solveNoPivotLDLUpperOnDeviceVar1(const ordinal_type pbeg, const ordinal_type pend,
                                              const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const ordinal_type nrhs = t.extent(1);
    const value_type minus_one(-1), one(1);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    int num_device_calls = 0;
    exec_space exec_instance;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type handle_blas = getBlasHandle(qid);
        setStreamOnHandle(qid);
        exec_instance = _exec_instances[qid];
        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
#endif
        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            value_type *aptr = s.u_buf, *bptr = _buf.data() + h_buf_solve_ptr(p - pbeg);

            const UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            aptr += m * m;
            const UnmanagedViewType<value_type_matrix> bT(bptr, m, nrhs);
            bptr += m * nrhs;

            const ordinal_type offm = s.row_begin;
            const auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());

            // Update with off-diagonal
            if (n_m > 0) {
              const UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // aptr += m*n;
              const UnmanagedViewType<value_type_matrix> bB(bptr, n_m, nrhs);
              _status = Gemv<Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, minus_one, ATR, bB, one, tT);
              checkDeviceBlasStatus("gemv");
            }

            // Copy current block (NOTE: should be in Trmv?)
            _status = Copy<Algo::OnDevice>::invoke(exec_instance, bT, tT);
            checkDeviceBlasStatus("Copy");

            // Solve with diagonal block
            // invoke Gemv_OnDevice with handle_blas as member (which will call blas_/cublas_/rocblas_invoke)
            _status = Trmv<Uplo::Upper, Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, Diag::Unit(), ATL, tT, bT);
            checkDeviceBlasStatus("gemv");

            _status = Copy<Algo::OnDevice>::invoke(exec_instance, tT, bT);
            checkDeviceBlasStatus("Copy");

            // increment num device calls
            num_device_calls ++;
          }
        }
      }
    }
    return num_device_calls;
  }

  inline int solveNoPivotLDLUpperOnDeviceVar2(const ordinal_type pbeg, const ordinal_type pend,
                                              const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const ordinal_type nrhs = t.extent(1);
    exec_space exec_instance;
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    int num_device_calls = 0;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type handle_blas = getBlasHandle(qid);
        setStreamOnHandle(qid);
        exec_instance = _exec_instances[qid];
        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n;
          if (m > 0 && n > 0) {
            value_type *aptr = s.u_buf, *bptr = _buf.data() + h_buf_solve_ptr(p - pbeg);
 
            const UnmanagedViewType<value_type_matrix> AT(aptr, m, n);
            const UnmanagedViewType<value_type_matrix> b(bptr, n, nrhs);

            const ordinal_type offm = s.row_begin;
            const auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());

            // Copy current block (NOTE: should be in Trmv?)
            const auto dx = Kokkos::subview(b, range_type(0, m), Kokkos::ALL());
            _status = Copy<Algo::OnDevice>::invoke(exec_instance, tT, dx);
            checkDeviceBlasStatus("copy");

            _status = Trmv<Uplo::Upper, Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, Diag::Unit(), AT, b, tT);
            checkDeviceBlasStatus("trmv");

            // increment num device calls
            num_device_calls ++;
          }
        }
      }
    }
    return num_device_calls;
  }

  inline int solveNoPivotLDLUpperOnDevice(const ordinal_type lvl, const ordinal_type nlvls,
                                          const ordinal_type pbeg, const ordinal_type pend,
                                          const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    if (variant == 0)
      return solveNoPivotLDLUpperOnDeviceVar0(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 1)
      return solveNoPivotLDLUpperOnDeviceVar1(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 2)
      return solveNoPivotLDLUpperOnDeviceVar2(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 3) {
      auto &s0 = _h_supernodes(_h_level_sids(pbeg));
      _spmv->ApplyU_OnDevice(lvl, s0, t, _w_vec);
    } else {
      std::string msg = "Error: LevelSetTools::solveNoPivotLDLUpperOnDevice, algorithm variant ("
                        + std::to_string(variant) + ") is not supported.\n";
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, msg.c_str());
    }
    return 0;
  }

  /// 
  /// Cholesky Lower
  ///
  inline int solveCholeskyLowerOnDeviceVar0(const ordinal_type pbeg, const ordinal_type pend,
                                            const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const ordinal_type nrhs = t.extent(1);
    const value_type minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    int num_device_calls = 0;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type handle_blas = getBlasHandle(qid);
        setStreamOnHandle(qid);
        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            value_type *aptr = s.u_buf;
            UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            aptr += ATL.span();

            const ordinal_type offm = s.row_begin;
            auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());
            _status =
                Trsv<Uplo::Upper, Trans::ConjTranspose, Algo::OnDevice>::invoke(handle_blas, Diag::NonUnit(), ATL, tT);
            checkDeviceBlasStatus("trsv");

            if (n_m > 0) {
              // solve offdiag
              value_type *bptr = _buf.data() + h_buf_solve_ptr(p - pbeg);
              UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // aptr += m*n_m;
              UnmanagedViewType<value_type_matrix> bB(bptr, n_m, nrhs);

              _status = Gemv<Trans::ConjTranspose, Algo::OnDevice>::invoke(handle_blas, minus_one, ATR, tT, zero, bB);
              checkDeviceBlasStatus("gemv");
            }
            // increment num device calls
            num_device_calls ++;
          }
        }
      }
    }
    return num_device_calls;
  }

  inline int solveCholeskyLowerOnDeviceVar1(const ordinal_type pbeg, const ordinal_type pend,
                                            const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const ordinal_type nrhs = t.extent(1);
    const value_type one(1), minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    int num_device_calls = 0;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type handle_blas = getBlasHandle(qid);
        setStreamOnHandle(qid);
        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            value_type *aptr = s.u_buf;
            UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            aptr += m * m;

            value_type *bptr = _buf.data() + h_buf_solve_ptr(p - pbeg);
            UnmanagedViewType<value_type_matrix> bT(bptr, m, nrhs);
            bptr += m * nrhs;

            const ordinal_type offm = s.row_begin;
            const auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());

            // Solve diag
            _status = Gemv<Trans::ConjTranspose, Algo::OnDevice>::invoke(handle_blas, one, ATL, tT, zero, bT);
            checkDeviceBlasStatus("gemv");

            if (n_m > 0) {
              // Update offdiag
              UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m);
              UnmanagedViewType<value_type_matrix> bB(bptr, n_m, nrhs);

              _status = Gemv<Trans::ConjTranspose, Algo::OnDevice>::invoke(handle_blas, minus_one, ATR, bT, zero, bB);
              checkDeviceBlasStatus("gemv");
            }
            // increment num device calls
            num_device_calls ++;
          }
        }
      }
    }
    return num_device_calls;
  }

  inline int solveCholeskyLowerOnDeviceVar2(const ordinal_type pbeg, const ordinal_type pend,
                                            const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const value_type one(1), zero(0);
    const ordinal_type nrhs = t.extent(1);

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    int num_device_calls = 0;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        //const ordinal_type qid = 0; //q % _nstreams;
        const ordinal_type qid = q % _nstreams;
        blas_handle_type handle_blas = getBlasHandle(qid);
        setStreamOnHandle(qid);
        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n;
          if (m > 0) {
            value_type *aptr = s.u_buf;
            UnmanagedViewType<value_type_matrix> AT(aptr, m, n);

            value_type *bptr = _buf.data() + h_buf_solve_ptr(p - pbeg);
            UnmanagedViewType<value_type_matrix> b(bptr, n, nrhs);

            const ordinal_type offm = s.row_begin;
            auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());

            _status = Gemv<Trans::ConjTranspose, Algo::OnDevice>::invoke(handle_blas, one, AT, tT, zero, b);
            checkDeviceBlasStatus("gemv");

            // increment num device calls
            num_device_calls ++;
          }
        }
      }
    }
    return num_device_calls;
  }

  inline int solveCholeskyLowerOnDevice(const ordinal_type lvl, const ordinal_type nlvls,
                                        const ordinal_type pbeg, const ordinal_type pend,
                                        const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    if (variant == 0)
      return solveCholeskyLowerOnDeviceVar0(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 1)
      return solveCholeskyLowerOnDeviceVar1(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 2)
      return solveCholeskyLowerOnDeviceVar2(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 3) {
      auto &s0 = _h_supernodes(_h_level_sids(pbeg));
      _spmv->ApplyL_OnDevice(lvl, s0, t, _w_vec);
    } else {
      std::string msg = "Error: LevelSetTools::solveCholeskyLowerOnDevice, algorithm variant ("
                        + std::to_string(variant) + ") is not supported.\n";
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, msg.c_str());
    }
    return 0;
  }

  /// 
  /// Cholesky Upper
  ///
  inline int solveCholeskyUpperOnDeviceVar0(const ordinal_type pbeg, const ordinal_type pend,
                                            const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const ordinal_type nrhs = t.extent(1);
    const value_type minus_one(-1), one(1);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    int num_device_calls = 0;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type handle_blas = getBlasHandle(qid);
        setStreamOnHandle(qid);
        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            value_type *aptr = s.u_buf, *bptr = _buf.data() + h_buf_solve_ptr(p - pbeg);

            const UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            aptr += m * m;
            const UnmanagedViewType<value_type_matrix> bB(bptr, n_m, nrhs);

            const ordinal_type offm = s.row_begin;
            const auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());

            if (n_m > 0) {
              const UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // aptr += m*n;
              _status = Gemv<Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, minus_one, ATR, bB, one, tT);
              checkDeviceBlasStatus("gemv");
            }
            _status =
                Trsv<Uplo::Upper, Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, Diag::NonUnit(), ATL, tT);
            checkDeviceBlasStatus("trsv");

            // increment num device calls
            num_device_calls ++;
          }
        }
      }
    }
    return num_device_calls;
  }

  inline int solveCholeskyUpperOnDeviceVar1(const ordinal_type pbeg, const ordinal_type pend,
                                            const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const ordinal_type nrhs = t.extent(1);
    const value_type minus_one(-1), one(1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    int num_device_calls = 0;
    exec_space exec_instance;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type handle_blas = getBlasHandle(qid);
        setStreamOnHandle(qid);
        exec_instance = _exec_instances[qid];
        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
#endif
        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            value_type *aptr = s.u_buf, *bptr = _buf.data() + h_buf_solve_ptr(p - pbeg);

            const UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            aptr += m * m;
            const UnmanagedViewType<value_type_matrix> bT(bptr, m, nrhs);
            bptr += m * nrhs;

            const ordinal_type offm = s.row_begin;
            const auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());

            // Update with off-diagonal
            if (n_m > 0) {
              const UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // aptr += m*n;
              const UnmanagedViewType<value_type_matrix> bB(bptr, n_m, nrhs);
              _status = Gemv<Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, minus_one, ATR, bB, one, tT);
              checkDeviceBlasStatus("gemv");
            }

            // Solve with diagonal block
            // invoke Gemv_OnDevice with handle_blas as member (which will call blas_/cublas_/rocblas_invoke)
            _status = Gemv<Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, one, ATL, tT, zero, bT);
            checkDeviceBlasStatus("gemv");

            _status = Copy<Algo::OnDevice>::invoke(exec_instance, tT, bT);
            checkDeviceBlasStatus("Copy");

            // increment num device calls
            num_device_calls ++;
          }
        }
      }
    }
    return num_device_calls;
  }

  inline int solveCholeskyUpperOnDeviceVar2(const ordinal_type pbeg, const ordinal_type pend,
                                            const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const value_type one(1), zero(0);
    const ordinal_type nrhs = t.extent(1);

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    int num_device_calls = 0;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        //const ordinal_type qid = 0; //q % _nstreams;
        blas_handle_type handle_blas = getBlasHandle(qid);
        setStreamOnHandle(qid);
        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n;
          if (m > 0 && n > 0) {
            value_type *aptr = s.u_buf, *bptr = _buf.data() + h_buf_solve_ptr(p - pbeg);
 
            const UnmanagedViewType<value_type_matrix> AT(aptr, m, n);
            const UnmanagedViewType<value_type_matrix> b(bptr, n, nrhs);

            const ordinal_type offm = s.row_begin;
            const auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());

            _status = Gemv<Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, one, AT, b, zero, tT);
            checkDeviceBlasStatus("gemv");

            // increment num device calls
            num_device_calls ++;
          }
        }
      }
    }
    return num_device_calls;
  }

  inline int solveCholeskyUpperOnDevice(const ordinal_type lvl, const ordinal_type nlvls,
                                        const ordinal_type pbeg, const ordinal_type pend,
                                        const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    if (variant == 0)
      return solveCholeskyUpperOnDeviceVar0(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 1)
      return solveCholeskyUpperOnDeviceVar1(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 2)
      return solveCholeskyUpperOnDeviceVar2(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 3) {
      auto &s0 = _h_supernodes(_h_level_sids(pbeg));
      _spmv->ApplyU_OnDevice(lvl, s0, t, _w_vec);
    } else {
      std::string msg = "Error: LevelSetTools::solveCholeskyUpperOnDevice, algorithm variant ("
                        + std::to_string(variant) + ") is not supported.\n";
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, msg.c_str());
    }
    return 0;
  }

  /// 
  /// LDL Lower
  ///
  inline void solveLDL_LowerOnDeviceVar0(const ordinal_type pbeg, const ordinal_type pend,
                                         const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const ordinal_type nrhs = t.extent(1);
    const value_type minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    exec_space exec_instance;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type handle_blas = getBlasHandle(qid);
        setStreamOnHandle(qid);
        exec_instance = _exec_instances[qid];
        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            value_type *aptr = s.u_buf;
            UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            aptr += m * m;

            const ordinal_type offm = s.row_begin;

            const auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());
            if (!s.do_not_apply_pivots) {
              const auto fpiv = ordinal_type_array(_piv.data() + 4 * offm + m, m);
              _status = ApplyPivots<PivotMode::Flame, Side::Left, Direct::Forward, Algo::OnDevice> /// row inter-change
                  ::invoke(exec_instance, fpiv, tT);
            }

            _status =
                Trsv<Uplo::Lower, Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, Diag::Unit(), ATL, tT);
            checkDeviceBlasStatus("trsv");
            if (n_m > 0) {
              value_type *bptr = _buf.data() + h_buf_solve_ptr(p - pbeg);
              UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // ptr += m*n_m;
              UnmanagedViewType<value_type_matrix> bB(bptr, n_m, nrhs);
              _status = Gemv<Trans::Transpose, Algo::OnDevice>::invoke(handle_blas, minus_one, ATR, tT, zero, bB);
              checkDeviceBlasStatus("gemv");
            }
          }
        }
      }
    }
  }

  inline void solveLDL_LowerOnDeviceVar1(const ordinal_type pbeg, const ordinal_type pend,
                                         const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const ordinal_type nrhs = t.extent(1);
    const value_type one(1), minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    exec_space exec_instance;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type handle_blas = getBlasHandle(qid);
        setStreamOnHandle(qid);
        exec_instance = _exec_instances[qid];
        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            value_type *aptr = s.u_buf;
            UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            aptr += m * m;

            value_type *bptr = _buf.data() + h_buf_solve_ptr(p - pbeg);
            UnmanagedViewType<value_type_matrix> bT(bptr, m, nrhs);
            bptr += m * nrhs;

            const ordinal_type offm = s.row_begin;

            const auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());

            if (s.do_not_apply_pivots) {
              _status = Copy<Algo::OnDevice>::invoke(exec_instance, bT, tT);
            } else {
              ConstUnmanagedViewType<ordinal_type_array> perm(_piv.data() + 4 * offm + 2 * m, m);
              _status =
                  ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(exec_instance, tT, perm, bT);
            }

            _status = Gemv<Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, one, ATL, bT, zero, tT);
            checkDeviceBlasStatus("gemv");

            if (n_m > 0) {
              UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m);
              UnmanagedViewType<value_type_matrix> bB(bptr, n_m, nrhs);

              _status = Gemv<Trans::Transpose, Algo::OnDevice>::invoke(handle_blas, minus_one, ATR, tT, zero, bB);
              checkDeviceBlasStatus("gemv");
            }
          }
        }
      }
    }
  }

  inline void solveLDL_LowerOnDeviceVar2(const ordinal_type pbeg, const ordinal_type pend,
                                         const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const ordinal_type nrhs = t.extent(1);
    const value_type one(1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    exec_space exec_instance;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type handle_blas = getBlasHandle(qid);
        setStreamOnHandle(qid);
        exec_instance = _exec_instances[qid];
        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n;
          if (m > 0) {
            value_type *aptr = s.u_buf, *bptr = _buf.data() + h_buf_solve_ptr(p - pbeg);
            UnmanagedViewType<value_type_matrix> AT(aptr, m, n);
            UnmanagedViewType<value_type_matrix> b(bptr, n, nrhs);

            const ordinal_type offm = s.row_begin;
            const auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());

            if (!s.do_not_apply_pivots) {
              UnmanagedViewType<value_type_matrix> bT(bptr, m, nrhs);
              ConstUnmanagedViewType<ordinal_type_array> perm(_piv.data() + 4 * offm + 2 * m, m);
              _status = Copy<Algo::OnDevice>::invoke(exec_instance, bT, tT);
              _status =
                  ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(exec_instance, bT, perm, tT);
            }

            _status = Gemv<Trans::Transpose, Algo::OnDevice>::invoke(handle_blas, one, AT, tT, zero, b);
            checkDeviceBlasStatus("gemv");
          }
        }
      }
    }
  }

  inline void solveLDL_LowerOnDevice(const ordinal_type lvl, const ordinal_type nlvls,
                                     const ordinal_type pbeg, const ordinal_type pend,
                                     const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    if (variant == 0)
      solveLDL_LowerOnDeviceVar0(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 1)
      solveLDL_LowerOnDeviceVar1(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 2)
      solveLDL_LowerOnDeviceVar2(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 3) {
      auto &s0 = _h_supernodes(_h_level_sids(pbeg));
      _spmv->ApplyL_OnDevice(lvl, s0, t, _w_vec);
    } else {
      std::string msg = "Error: LevelSetTools::solveLDL_LowerOnDevice, algorithm variant ("
                        + std::to_string(variant) + ") is not supported.\n";
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, msg.c_str());
    }
  }

  /// 
  /// LDL Upper
  ///
  inline void solveLDL_UpperOnDeviceVar0(const ordinal_type pbeg, const ordinal_type pend,
                                         const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const ordinal_type nrhs = t.extent(1);
    const value_type minus_one(-1), one(1);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    exec_space exec_instance;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type handle_blas = getBlasHandle(qid);
        setStreamOnHandle(qid);
        exec_instance = _exec_instances[qid];
        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            value_type *aptr = s.u_buf, *bptr = _buf.data() + h_buf_solve_ptr(p - pbeg);
            ;
            const UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            aptr += m * m;
            const UnmanagedViewType<value_type_matrix> bB(bptr, n_m, nrhs);

            const ordinal_type offm = s.row_begin;
            const auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());
            const auto P = ordinal_type_array(_piv.data() + 4 * offm, 4 * m);
            const auto D = value_type_matrix(_diag.data() + 2 * offm, m, 2);
            _status = Scale2x2_BlockInverseDiagonals<Side::Left, Algo::OnDevice> /// row scaling
                ::invoke(exec_instance, P, D, tT);

            if (n_m > 0) {
              const UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // aptr += m*n;
              _status = Gemv<Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, minus_one, ATR, bB, one, tT);
              checkDeviceBlasStatus("gemv");
            }
            _status = Trsv<Uplo::Lower, Trans::Transpose, Algo::OnDevice>::invoke(handle_blas, Diag::Unit(), ATL, tT);
            checkDeviceBlasStatus("trsv");

            if (!s.do_not_apply_pivots) {
              const auto fpiv = ordinal_type_array(P.data() + m, m);
              _status = ApplyPivots<PivotMode::Flame, Side::Left, Direct::Backward, Algo::OnDevice> /// row inter-change
                  ::invoke(exec_instance, fpiv, tT);
            }
          }
        }
      }
    }
  }

  inline void solveLDL_UpperOnDeviceVar1(const ordinal_type pbeg, const ordinal_type pend,
                                         const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const ordinal_type nrhs = t.extent(1);
    const value_type minus_one(-1), one(1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    exec_space exec_instance;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type handle_blas = getBlasHandle(qid);
        setStreamOnHandle(qid);
        exec_instance = _exec_instances[qid];
        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            value_type *aptr = s.u_buf, *bptr = _buf.data() + h_buf_solve_ptr(p - pbeg);

            const UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            aptr += ATL.span();
            const UnmanagedViewType<value_type_matrix> bT(bptr, m, nrhs);
            bptr += bT.span();

            const ordinal_type offm = s.row_begin;
            const auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());
            const auto P = ordinal_type_array(_piv.data() + 4 * offm, 4 * m);
            const auto D = value_type_matrix(_diag.data() + 2 * offm, m, 2);
            _status = Scale2x2_BlockInverseDiagonals<Side::Left, Algo::OnDevice> /// row scaling
                ::invoke(exec_instance, P, D, tT);

            if (n_m > 0) {
              const UnmanagedViewType<value_type_matrix> bB(bptr, n_m, nrhs);
              const UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // aptr += m*n;
              _status = Gemv<Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, minus_one, ATR, bB, one, tT);
              checkDeviceBlasStatus("gemv");
            }

            _status = Gemv<Trans::Transpose, Algo::OnDevice>::invoke(handle_blas, one, ATL, tT, zero, bT);
            checkDeviceBlasStatus("gemv");

            if (s.do_not_apply_pivots) {
              _status = Copy<Algo::OnDevice>::invoke(exec_instance, tT, bT);
            } else {
              ConstUnmanagedViewType<ordinal_type_array> peri(P.data() + 3 * m, m);
              _status =
                  ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(exec_instance, bT, peri, tT);
            }
          }
        }
      }
    }
  }

  inline void solveLDL_UpperOnDeviceVar2(const ordinal_type pbeg, const ordinal_type pend,
                                         const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const ordinal_type nrhs = t.extent(1);
    const value_type one(1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    exec_space exec_instance;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type handle_blas = getBlasHandle(qid);
        setStreamOnHandle(qid);
        exec_instance = _exec_instances[qid];
        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n;
          if (m > 0 && n > 0) {
            value_type *aptr = s.u_buf, *bptr = _buf.data() + h_buf_solve_ptr(p - pbeg);

            const UnmanagedViewType<value_type_matrix> AT(aptr, m, n);
            const UnmanagedViewType<value_type_matrix> b(bptr, n, nrhs);

            const ordinal_type offm = s.row_begin;
            const auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());
            const auto bT = Kokkos::subview(b, range_type(0, m), Kokkos::ALL());

            ConstUnmanagedViewType<ordinal_type_array> P(_piv.data() + offm * 4, m * 4);
            ConstUnmanagedViewType<value_type_matrix> D(_diag.data() + offm * 2, m, 2);

            _status = Scale2x2_BlockInverseDiagonals<Side::Left, Algo::OnDevice> /// row scaling
                ::invoke(exec_instance, P, D, bT);

            _status = Gemv<Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, one, AT, b, zero, tT);

            if (!s.do_not_apply_pivots) {
              _status = Copy<Algo::OnDevice>::invoke(exec_instance, bT, tT);

              ConstUnmanagedViewType<ordinal_type_array> peri(P.data() + 3 * m, m);
              _status =
                  ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(exec_instance, bT, peri, tT);
            }
          }
        }
      }
    }
  }

  inline void solveLDL_DiagOnDevice(const ordinal_type pbeg, const ordinal_type pend,
                                    const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    exec_space exec_instance;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      //if (_h_solve_mode(sid) == 0)
      {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        setStreamOnHandle(qid);
        exec_instance = _exec_instances[qid];
        ++q;
#endif

        const auto &s = _h_supernodes(sid);
        {
          if (s.m > 0) {
            const ordinal_type m = s.m;
            const ordinal_type offm = s.row_begin;
            const auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());
            const auto P = ordinal_type_array(_piv.data() + 4 * offm, 4 * m);
            const auto D = value_type_matrix(_diag.data() + 2 * offm, m, 2);
            _status = Scale2x2_BlockInverseDiagonals<Side::Left, Algo::OnDevice> /// row scaling
                ::invoke(exec_instance, P, D, tT);
          }
        }
      }
    }
  }

  inline void solveLDL_UpperOnDevice(const ordinal_type lvl, const ordinal_type nlvls,
                                     const ordinal_type pbeg, const ordinal_type pend,
                                     const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    if (variant == 0) {
      solveLDL_UpperOnDeviceVar0(pbeg, pend, h_buf_solve_ptr, t);
    } else if (variant == 1) {
      solveLDL_UpperOnDeviceVar1(pbeg, pend, h_buf_solve_ptr, t);
    } else if (variant == 2) {
      solveLDL_UpperOnDeviceVar2(pbeg, pend, h_buf_solve_ptr, t);
    } else if (variant == 3) {
      // NOTE: merge D into U, or convert D into Crs?
      solveLDL_DiagOnDevice(pbeg, pend, h_buf_solve_ptr, (lvl%2 == 0 ? t : _w_vec));
      Kokkos::fence();
      auto &s0 = _h_supernodes(_h_level_sids(pbeg));
      _spmv->ApplyU_OnDevice(lvl, s0, t, _w_vec);
    } else {
      std::string msg = "Error: LevelSetTools::solveLDL_UpperOnDevice, algorithm variant ("
                        + std::to_string(variant) + ") is not supported.\n";
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, msg.c_str());
    }
  }

  /// 
  /// LU Lower
  ///
  inline int solveLU_LowerOnDeviceVar0(const ordinal_type pbeg, const ordinal_type pend,
                                       const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const ordinal_type nrhs = t.extent(1);
    const value_type minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    int num_device_calls = 0;
    exec_space exec_instance;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type handle_blas = getBlasHandle(qid);
        exec_instance = _exec_instances[qid];
        setStreamOnHandle(qid);
        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            UnmanagedViewType<value_type_matrix> ATL(s.u_buf, m, m);

            const ordinal_type offm = s.row_begin;
            const auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());
            const auto fpiv = ordinal_type_array(_piv.data() + 4 * offm + m, m);

            _status = ApplyPivots<PivotMode::Flame, Side::Left, Direct::Forward, Algo::OnDevice> /// row inter-change
                ::invoke(exec_instance, fpiv, tT);

            _status =
                Trsv<Uplo::Lower, Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, Diag::Unit(), ATL, tT);
            checkDeviceBlasStatus("trsv");

            if (n_m > 0) {
              value_type *bptr = _buf.data() + h_buf_solve_ptr(p - pbeg);
              UnmanagedViewType<value_type_matrix> AL(s.l_buf, n, m);
              const auto ABL = Kokkos::subview(AL, range_type(m, n), Kokkos::ALL());
              UnmanagedViewType<value_type_matrix> bB(bptr, n_m, nrhs);
              _status = Gemv<Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, minus_one, ABL, tT, zero, bB);
              checkDeviceBlasStatus("gemv");
            }
            num_device_calls ++;
          }
        }
      }
    }
    return num_device_calls;
  }

  inline int solveLU_LowerOnDeviceVar1(const ordinal_type pbeg, const ordinal_type pend,
                                       const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const ordinal_type nrhs = t.extent(1);
    const value_type one(1), minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    int num_device_calls = 0;
    exec_space exec_instance;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type handle_blas = getBlasHandle(qid);
        setStreamOnHandle(qid);
        exec_instance = _exec_instances[qid];
        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            UnmanagedViewType<value_type_matrix> AL(s.l_buf, n, m);
            const auto ATL = Kokkos::subview(AL, range_type(0, m), Kokkos::ALL());

            value_type *bptr = _buf.data() + h_buf_solve_ptr(p - pbeg);
            UnmanagedViewType<value_type_matrix> bT(bptr, m, nrhs);

            const ordinal_type offm = s.row_begin;

            const auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());
            ConstUnmanagedViewType<ordinal_type_array> perm(_piv.data() + 4 * offm + 2 * m, m);

            if (s.do_not_apply_pivots) {
              _status = Copy<Algo::OnDevice>::invoke(exec_instance, bT, tT);
            } else {
              _status =
                  ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(exec_instance, tT, perm, bT);
            }

            _status = Gemv<Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, one, ATL, bT, zero, tT);
            checkDeviceBlasStatus("gemv");

            if (n_m > 0) {
              const auto ABL = Kokkos::subview(AL, range_type(m, n), Kokkos::ALL());
              UnmanagedViewType<value_type_matrix> bB(bptr + bT.span(), n_m, nrhs);
              _status = Gemv<Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, minus_one, ABL, tT, zero, bB);
              checkDeviceBlasStatus("gemv");
            }
            num_device_calls ++;
          }
        }
      }
    }
    return num_device_calls;
  }

  inline int solveLU_LowerOnDeviceVar2(const ordinal_type pbeg, const ordinal_type pend,
                                       const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const ordinal_type nrhs = t.extent(1);
    const value_type one(1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    int num_device_calls = 0;
    exec_space exec_instance;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type handle_blas = getBlasHandle(qid);
        setStreamOnHandle(qid);
        exec_instance = _exec_instances[qid];
        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n;
          if (m > 0) {
            UnmanagedViewType<value_type_matrix> AL(s.l_buf, n, m);

            value_type *bptr = _buf.data() + h_buf_solve_ptr(p - pbeg);
            UnmanagedViewType<value_type_matrix> b(bptr, n, nrhs);

            const ordinal_type offm = s.row_begin;
            const auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());

            if (!s.do_not_apply_pivots) {
              UnmanagedViewType<value_type_matrix> bT(bptr, m, nrhs);
              ConstUnmanagedViewType<ordinal_type_array> perm(_piv.data() + 4 * offm + 2 * m, m);
              _status = Copy<Algo::OnDevice>::invoke(exec_instance, bT, tT);
              _status =
                  ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(exec_instance, bT, perm, tT);
            }
            _status = Gemv<Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, one, AL, tT, zero, b);
            checkDeviceBlasStatus("gemv");
            num_device_calls ++;
          }
        }
      }
    }
    return num_device_calls;
  }

  inline int solveLU_LowerOnDevice(const ordinal_type lvl, const ordinal_type nlvls,
                                   const ordinal_type pbeg, const ordinal_type pend,
                                   const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    if (variant == 0)
      return solveLU_LowerOnDeviceVar0(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 1)
      return solveLU_LowerOnDeviceVar1(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 2)
      return solveLU_LowerOnDeviceVar2(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 3) {
      // L (stored by cols) incorporate partial-pivoting
      auto &s0 = _h_supernodes(_h_level_sids(pbeg));
      _spmv->ApplyL_OnDevice(lvl, s0, t, _w_vec);
    } else {
      std::string msg = "Error: LevelSetTools::solveLU_LowerOnDevice, algorithm variant ("
                        + std::to_string(variant) + ") is not supported.\n";
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, msg.c_str());
    }
    return 0;
  }

  /// 
  /// LU Upper
  ///
  inline int solveLU_UpperOnDeviceVar0(const ordinal_type pbeg, const ordinal_type pend,
                                       const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const ordinal_type nrhs = t.extent(1);
    const value_type minus_one(-1), one(1);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    int num_device_calls = 0;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type handle_blas = getBlasHandle(qid);
        setStreamOnHandle(qid);
        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            value_type *uptr = s.u_buf, *bptr = _buf.data() + h_buf_solve_ptr(p - pbeg);
            ;
            const UnmanagedViewType<value_type_matrix> ATL(uptr, m, m);
            uptr += m * m;
            const UnmanagedViewType<value_type_matrix> bB(bptr, n_m, nrhs);

            const ordinal_type offm = s.row_begin;
            const auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());
            if (n_m > 0) {
              const UnmanagedViewType<value_type_matrix> ATR(uptr, m, n_m); // uptr += m*n;
              _status = Gemv<Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, minus_one, ATR, bB, one, tT);
              checkDeviceBlasStatus("gemv");
            }
            _status =
                Trsv<Uplo::Upper, Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, Diag::NonUnit(), ATL, tT);
            checkDeviceBlasStatus("trsv");
          }
          num_device_calls ++;
        }
      }
    }
    return num_device_calls;
  }

  inline int solveLU_UpperOnDeviceVar1(const ordinal_type pbeg, const ordinal_type pend,
                                       const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const ordinal_type nrhs = t.extent(1);
    const value_type minus_one(-1), one(1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    int num_device_calls = 0;
    exec_space exec_instance;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type handle_blas = getBlasHandle(qid);
        setStreamOnHandle(qid);
        exec_instance = _exec_instances[qid];
        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            value_type *bptr = _buf.data() + h_buf_solve_ptr(p - pbeg);

            const UnmanagedViewType<value_type_matrix> ATL(s.u_buf, m, m);
            const UnmanagedViewType<value_type_matrix> bT(bptr, m, nrhs);

            const ordinal_type offm = s.row_begin;
            const auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());

            _status = Copy<Algo::OnDevice>::invoke(exec_instance, bT, tT);

            if (n_m > 0) {
              const UnmanagedViewType<value_type_matrix> ATR(s.u_buf + ATL.span(), m, n_m);
              const UnmanagedViewType<value_type_matrix> bB(bptr + bT.span(), n_m, nrhs);
              _status = Gemv<Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, minus_one, ATR, bB, one, bT);
              checkDeviceBlasStatus("gemv");
            }

            _status = Gemv<Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, one, ATL, bT, zero, tT);
            checkDeviceBlasStatus("gemv");
            num_device_calls ++;
          }
        }
      }
    }
    return num_device_calls;
  }

  inline int solveLU_UpperOnDeviceVar2(const ordinal_type pbeg, const ordinal_type pend,
                                       const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const ordinal_type nrhs = t.extent(1);
    const value_type one(1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
    int num_device_calls = 0;
    for (ordinal_type p = pbeg; p < pend; ++p) {
      const ordinal_type sid = _h_level_sids(p);
      if (_h_solve_mode(sid) == 0) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        const ordinal_type qid = q % _nstreams;
        blas_handle_type handle_blas = getBlasHandle(qid);
        setStreamOnHandle(qid);
        ++q;
#else
        blas_handle_type handle_blas = getBlasHandle();
#endif

        const auto &s = _h_supernodes(sid);
        {
          const ordinal_type m = s.m, n = s.n;
          if (m > 0) {
            value_type *bptr = _buf.data() + h_buf_solve_ptr(p - pbeg);

            const UnmanagedViewType<value_type_matrix> AT(s.u_buf, m, n);
            const UnmanagedViewType<value_type_matrix> b(bptr, n, nrhs);

            const ordinal_type offm = s.row_begin;
            const auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());
            _status = Gemv<Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, one, AT, b, zero, tT);
            checkDeviceBlasStatus("gemv");
            num_device_calls ++;
          }
        }
      }
    }
    return num_device_calls;
  }

  inline int solveLU_UpperOnDevice(const ordinal_type lvl, const ordinal_type nlvls,
                                   const ordinal_type pbeg, const ordinal_type pend,
                                   const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    if (variant == 0)
      return solveLU_UpperOnDeviceVar0(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 1)
      return solveLU_UpperOnDeviceVar1(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 2)
      return solveLU_UpperOnDeviceVar2(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 3) {
      auto &s0 = _h_supernodes(_h_level_sids(pbeg));
      _spmv->ApplyU_OnDevice(lvl, s0, t, _w_vec);
    } else {
      std::string msg = "Error: LevelSetTools::solveLU_UpperOnDevice, algorithm variant ("
                        + std::to_string(variant) + ") is not supported.\n";
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, msg.c_str());
    }
    return 0;
  }

  ///
  /// Level set factorize & solve
  ///

  /// 
  /// Cholesky
  ///
  inline void factorizeCholesky(const value_type_array &ax, const bool store_transpose,
                                const mag_type shift, const mag_type pivot_tol,
                                const ordinal_type verbose) {

    constexpr bool is_host = std::is_same<exec_memory_space, Kokkos::HostSpace>::value;
    Kokkos::Timer timer;
    Kokkos::Timer tick;
    double time_parallel = 0.0;
    double time_device = 0.0;
    double time_update = 0.0;

    timer.reset();
    if (_buf.span() < size_t(_bufsize_factorize)) {
      if (_buf.span() > 0) track_free(_buf.span() * sizeof(value_type));
      Kokkos::resize(_buf, _bufsize_factorize);
      track_alloc(_buf.span() * sizeof(value_type));
    }

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    {
      size_type worksize = _worksize * (_nstreams + 1);
      if (size_t(worksize) > _work.span()) {
        if (_work.span() > 0) track_free(_work.span() * sizeof(value_type));
        Kokkos::resize(_work, worksize);
        track_alloc(_work.span() * sizeof(value_type));
      }
    }
#endif
    stat.t_extra = timer.seconds();

    timer.reset();
    {
      _ax = ax; // matrix values
      constexpr bool copy_to_l_buf(false);
      _info.copySparseToSuperpanels(copy_to_l_buf, _ap, _aj, _ax, shift, _perm, _peri);
      exec_space().fence(); // wait for copy
    }
    stat.t_copy = timer.seconds();

    stat_level.n_kernel_launching = 0;
    timer.reset();
    {
      // this should be considered with average problem sizes in levels
      const ordinal_type half_level = _nlevel / 2;
      // const ordinal_type team_size_factor[2] = { 64, 16 }, vector_size_factor[2] = { 8, 8};
      // const ordinal_type team_size_factor[2] = { 16, 16 }, vector_size_factor[2] = { 32, 32};
      const ordinal_type team_size_factor[2] = {64, 64}, vector_size_factor[2] = {8, 4};
      const ordinal_type team_size_update[2] = {16, 8}, vector_size_update[2] = {32, 32};
      // returned value from team Chol
      colind_view d_rval("rval",1);
      auto h_rval = Kokkos::create_mirror_view(host_memory_space(), d_rval);
      {
        typedef TeamFunctor_FactorizeChol<supernode_info_type> functor_type;
        functor_type functor(_info, _factorize_mode, _level_sids, _buf, d_rval.data());
        if (pivot_tol > 0.0) {
          functor.setDiagPertubationTol(pivot_tol);
        }
        if (this->getSolutionMethod() == 0) {
          functor.setIndefiniteFactorization(true);
        }

#if defined(TACHO_TEST_LEVELSET_TOOLS_KERNEL_OVERHEAD)
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space, typename functor_type::DummyTag>
            team_policy_factorize;
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space, typename functor_type::DummyTag>
            team_policy_update;
#else
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space,
                                   typename functor_type::template FactorizeTag<variant>>
            team_policy_factor;
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space, typename functor_type::UpdateTag>
            team_policy_update;
#endif
        team_policy_factor policy_factor(1, 1, 1);
        team_policy_update policy_update(1, 1, 1);

        // get max vector size
        const ordinal_type vmax = policy_factor.vector_length_max();
        {
          for (ordinal_type lvl = (_team_serial_level_cut - 1); lvl >= 0; --lvl) {
            const ordinal_type pbeg = _h_level_ptr(lvl), pend = _h_level_ptr(lvl + 1), pcnt = pend - pbeg;

            const range_type range_buf_factor_ptr(_h_buf_level_ptr(lvl), _h_buf_level_ptr(lvl + 1));

            const auto buf_factor_ptr = Kokkos::subview(_buf_factor_ptr, range_buf_factor_ptr);
            functor.setRange(pbeg, pend);
            functor.setBufferPtr(buf_factor_ptr);
            if (is_host) {
              policy_factor = team_policy_factor(pcnt, 1, 1);
              policy_update = team_policy_update(pcnt, 1, 1);
            } else {
              const ordinal_type idx = lvl > half_level;
              // pick vector sizes
              ordinal_type vsize_factor = std::min(vector_size_factor[idx],vmax);
              ordinal_type vsize_update = std::min(vector_size_update[idx],vmax);
              // pick teamm sizes
              policy_factor = team_policy_factor(pcnt, 1, vsize_factor);
              policy_update = team_policy_update(pcnt, 1, vsize_update);
              const ordinal_type factor_tmax = policy_factor.team_size_max(functor, Kokkos::ParallelForTag());
              const ordinal_type update_tmax = policy_update.team_size_max(functor, Kokkos::ParallelForTag());;

              // create policies
              policy_factor = team_policy_factor(pcnt, std::min(team_size_factor[idx],factor_tmax), vsize_factor);
              policy_update = team_policy_update(pcnt, std::min(team_size_update[idx],update_tmax), vsize_update);
            }
            if (lvl < _device_level_cut) {
              // do nothing
              // Kokkos::parallel_for("factor lower", policy_factor, functor);
            } else {
              if (verbose) {
                Kokkos::fence(); tick.reset();
              }
              Kokkos::parallel_for("factor", policy_factor, functor);
              if (verbose) {
                Kokkos::fence(); time_parallel += tick.seconds();
              }
              ++stat_level.n_kernel_launching;
            }

            if (verbose) {
              Kokkos::fence(); tick.reset();
            }
            const auto h_buf_factor_ptr = Kokkos::subview(_h_buf_factor_ptr, range_buf_factor_ptr);
            if (this->getSolutionMethod() == 0) {
              factorizeNoPivotLDLOnDevice(pbeg, pend, h_buf_factor_ptr, _work);
            } else {
              factorizeCholeskyOnDevice(pbeg, pend, h_buf_factor_ptr, _work);
            }
            if (verbose) {
              Kokkos::fence(); time_device += tick.seconds();
              tick.reset();
            }
            Kokkos::deep_copy(h_rval, d_rval);
            int rval = h_rval(0);
            if (rval != 0) {
              TACHO_TEST_FOR_EXCEPTION(true, std::runtime_error, "POTRF (team) returns non-zero error code.");
            }
            if (_h_num_device_calls_factor(lvl) > 0)
              Kokkos::fence(); // sync device-calls before calling update

            Kokkos::parallel_for("update factor", policy_update, functor);
            if (lvl == 0 || _h_num_device_calls_factor(lvl-1) > 0)
              exec_space().fence(); // sync default, before next device factor calls
            if (verbose) {
              Kokkos::fence(); time_update += tick.seconds();
            }
            ++stat_level.n_kernel_launching;
          }
        }
      }
    } // end of Cholesky
    stat.t_factor = timer.seconds();
    timer.reset();
    if (variant == 3) {
      // compress each partitioned inverse at each level into CRS matrix
      extractCRS(store_transpose, verbose);
    }
    stat.t_extra += timer.seconds();

    if (verbose) {
      printf("Summary: LevelSetTools-Variant-%d (CholeskyFactorize)\n", variant);
      printf("=====================================================\n");
      printf( "\n  ** Team = %f s, Device = %f s, Update = %f s **\n",time_parallel,time_device,time_update );
      printf( "\n" );
      print_stat_factor();
      fflush(stdout);
    }
  }

  inline void solveCholesky(const value_type_matrix &x, // solution
                            const value_type_matrix &b, // right hand side
                            const value_type_matrix &t,
                            const ordinal_type verbose) { // temporary workspace (store permuted vectors)

    TACHO_TEST_FOR_EXCEPTION(x.extent(0) != b.extent(0) || x.extent(1) != b.extent(1) || x.extent(0) != t.extent(0) ||
                                 x.extent(1) != t.extent(1),
                             std::logic_error, "x, b, t, and w dimensions do not match");

    TACHO_TEST_FOR_EXCEPTION(x.data() == b.data() || x.data() == t.data(), std::logic_error,
                             "x, b, t, and w have the same data pointer");
    constexpr bool is_host = std::is_same<exec_memory_space, Kokkos::HostSpace>::value;

    // solve U^{H} (U x) = b
    const ordinal_type nrhs = x.extent(1);
    Kokkos::Timer timer;

    stat_level.n_kernel_launching = 0;

    // one-time operation when nrhs is changed
    timer.reset();
    allocateWorkspaceSolve(nrhs);

    // execution spaces for non-device calls
    const bool need_fence = (!_team_on_user_stream || _nstreams > 1);
    const auto perm_exec_instance = _exec_instances[0];
    const auto team_exec_instance = (_team_on_user_stream ? _exec_instances[0] : exec_space());

    // 0. permute (from METIS) and copy b -> t
    ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(perm_exec_instance, b, _perm, t);
    if (variant != 3 && need_fence) perm_exec_instance.fence();
    stat.t_extra = timer.seconds();

    timer.reset();
    {
#if defined(TACHO_ENABLE_SOLVE_CHOLESKY_USE_LIGHT_KERNEL)
      const auto work_item_property = Kokkos::Experimental::WorkItemProperty::HintLightWeight;
#endif
      // this should be considered with average problem sizes in levels
      const ordinal_type half_level = _nlevel / 2;
      const ordinal_type team_size_solve[2] = {16, 16}, vector_size_solve[2] = {8, 8};
      const ordinal_type team_size_update[2] = {128, 32}, vector_size_update[2] = {1, 1};
      {
        typedef TeamFunctor_SolveLowerChol<supernode_info_type> functor_type;
#if defined(TACHO_TEST_SOLVE_CHOLESKY_KERNEL_OVERHEAD)
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space, typename functor_type::DummyTag>
            team_policy_solve;
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space, typename functor_type::DummyTag>
            team_policy_update;
#else
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space,
                                   typename functor_type::template SolveTag<variant>>
            team_policy_solve;
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space,
                                   typename functor_type::template UpdateTag<variant>>
            team_policy_update;
#endif
        functor_type functor(_info, _solve_mode, _level_sids, t, _buf);

        if (this->getSolutionMethod() == 0) {
          functor.setIndefiniteFactorization(true);
        }

        team_policy_solve policy_solve(1, 1, 1);
        team_policy_update policy_update(1, 1, 1);

        //  1. U^{H} w = t
        {
          for (ordinal_type lvl = (_team_serial_level_cut - 1); lvl >= 0; --lvl) {
            const ordinal_type pbeg = _h_level_ptr(lvl), pend = _h_level_ptr(lvl + 1), pcnt = pend - pbeg;

            const range_type range_solve_buf_ptr(_h_buf_level_ptr(lvl), _h_buf_level_ptr(lvl + 1));

            const auto solve_buf_ptr = Kokkos::subview(_buf_solve_nrhs_ptr, range_solve_buf_ptr);
            functor.setRange(pbeg, pend);
            functor.setBufferPtr(solve_buf_ptr);
            if (is_host) {
              policy_solve = team_policy_solve(pcnt, 1, 1);
              policy_update = team_policy_update(pcnt, 1, 1);
            } else {
              const ordinal_type idx = lvl > half_level;
              policy_solve = team_policy_solve(team_exec_instance, pcnt, team_size_solve[idx], vector_size_solve[idx]);
              policy_update = team_policy_update(team_exec_instance, pcnt, team_size_update[idx], vector_size_update[idx]);
            }
#if defined(TACHO_ENABLE_SOLVE_CHOLESKY_USE_LIGHT_KERNEL)
            const auto policy_solve_with_work_property =
                Kokkos::Experimental::require(policy_solve, work_item_property);
            const auto policy_update_with_work_property =
                Kokkos::Experimental::require(policy_update, work_item_property);
#else
            const auto policy_solve_with_work_property = policy_solve;
            const auto policy_update_with_work_property = policy_update;
#endif
            if (variant != 3) {
              if (lvl < _device_level_cut) {
                // do nothing
                // Kokkos::parallel_for("solve lower", policy_solve, functor);
              } else {
                Kokkos::parallel_for("solve lower", policy_solve_with_work_property, functor);
                ++stat_level.n_kernel_launching;
              }
            }
            const auto h_buf_solve_ptr = Kokkos::subview(_h_buf_solve_nrhs_ptr, range_solve_buf_ptr);
            if (this->getSolutionMethod() == 0) {
              solveNoPivotLDLLowerOnDevice(lvl, _team_serial_level_cut, pbeg, pend, h_buf_solve_ptr, t);
            } else {
              solveCholeskyLowerOnDevice(lvl, _team_serial_level_cut, pbeg, pend, h_buf_solve_ptr, t);
            }
            if (variant != 3) {
              if (need_fence && _h_num_device_calls_solve(lvl) > 0)
                Kokkos::fence(); // synch device calls before batched update

              // copy from buffer to t
              Kokkos::parallel_for("update lower", policy_update_with_work_property, functor);
              if (need_fence && (lvl == 0 || _h_num_device_calls_solve(lvl-1) > 0))
                team_exec_instance.fence(); // sync default, for next device solve calls
              ++stat_level.n_kernel_launching;
            }
          }
        }
      } // end of lower tri solve
      {
        typedef TeamFunctor_SolveUpperChol<supernode_info_type> functor_type;
#if defined(TACHO_TEST_SOLVE_CHOLESKY_KERNEL_OVERHEAD)
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space, typename functor_type::DummyTag>
            team_policy_solve;
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space, typename functor_type::DummyTag>
            team_policy_update;
#else
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space,
                                   typename functor_type::template SolveTag<variant>>
            team_policy_solve;
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space,
                                   typename functor_type::template UpdateTag<variant>>
            team_policy_update;
#endif
        functor_type functor(_info, _solve_mode, _level_sids, t, _buf);
        if (this->getSolutionMethod() == 0) {
          functor.setIndefiniteFactorization(true);
        }
        team_policy_solve policy_solve(1, 1, 1);
        team_policy_update policy_update(1, 1, 1);

        //  2. U t = w;
        {
          for (ordinal_type lvl = 0; lvl < _team_serial_level_cut; ++lvl) {
            const ordinal_type pbeg = _h_level_ptr(lvl), pend = _h_level_ptr(lvl + 1), pcnt = pend - pbeg;

            const range_type range_solve_buf_ptr(_h_buf_level_ptr(lvl), _h_buf_level_ptr(lvl + 1));
            const auto solve_buf_ptr = Kokkos::subview(_buf_solve_nrhs_ptr, range_solve_buf_ptr);
            functor.setRange(pbeg, pend);
            functor.setBufferPtr(solve_buf_ptr);
            if (is_host) {
              policy_solve = team_policy_solve(pcnt, 1, 1);
              policy_update = team_policy_update(pcnt, 1, 1);
            } else {
              const ordinal_type idx = lvl > half_level;
              policy_solve = team_policy_solve(team_exec_instance, pcnt, team_size_solve[idx], vector_size_solve[idx]);
              policy_update = team_policy_update(team_exec_instance, pcnt, team_size_update[idx], vector_size_update[idx]);
            }
#if defined(TACHO_ENABLE_SOLVE_CHOLESKY_USE_LIGHT_KERNEL)
            const auto policy_solve_with_work_property =
                Kokkos::Experimental::require(policy_solve, work_item_property);
            const auto policy_update_with_work_property =
                Kokkos::Experimental::require(policy_update, work_item_property);
#else
            const auto policy_solve_with_work_property = policy_solve;
            const auto policy_update_with_work_property = policy_update;
#endif
            if (variant != 3) {
              Kokkos::parallel_for("update upper", policy_update_with_work_property, functor);
              if (need_fence && _h_num_device_calls_solve(lvl) > 0)
                team_exec_instance.fence(); // sync default, before next device solve-upper calls
              ++stat_level.n_kernel_launching;

              if (lvl < _device_level_cut) {
                // do nothing
                // Kokkos::parallel_for("solve upper", policy_solve, functor);
              } else {
                Kokkos::parallel_for("solve upper", policy_solve_with_work_property, functor);
                ++stat_level.n_kernel_launching; // no need to synch because synched after next update if needed
              }
            }
            const auto h_buf_solve_ptr = Kokkos::subview(_h_buf_solve_nrhs_ptr, range_solve_buf_ptr);
            if (this->getSolutionMethod() == 0) {
              solveNoPivotLDLUpperOnDevice(lvl, _team_serial_level_cut, pbeg, pend, h_buf_solve_ptr, t);
            } else {
              solveCholeskyUpperOnDevice(lvl, _team_serial_level_cut, pbeg, pend, h_buf_solve_ptr, t);
            }
            if (need_fence && _h_num_device_calls_solve(lvl) > 0)
              Kokkos::fence(); // synch device calls before next update
          }
        }
      } /// end of upper tri solve

    } // end of solve
    stat.t_solve = timer.seconds();

    // permute (from METIS) and copy t -> x
    if (variant != 3 && need_fence) Kokkos::fence(); // synch user or default streams
    timer.reset();
    ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(perm_exec_instance, t, _peri, x);
    perm_exec_instance.fence();
    stat.t_extra += timer.seconds();

    if (verbose) {
      printf("Summary: LevelSetTools-Variant-%d (Cholesky Solve: %3d)\n", variant, nrhs);
      printf("=======================================================\n");
      print_stat_solve();
      fflush(stdout);
    }
  }

  /// 
  /// LDL
  ///
  inline void factorizeLDL(const value_type_array &ax, const bool store_transpose, const mag_type shift,
                           const ordinal_type verbose) {
    constexpr bool is_host = std::is_same<exec_memory_space, Kokkos::HostSpace>::value;
    Kokkos::Timer timer;
    Kokkos::Timer tick;
    double time_parallel = 0.0;
    double time_device = 0.0;
    double time_update = 0.0;

    timer.reset();
    if (_buf.span() < size_t(_bufsize_factorize)) {
      if (_buf.span() > 0) track_free(_buf.span() * sizeof(value_type));
      Kokkos::resize(_buf, _bufsize_factorize);
      track_alloc(_buf.span() * sizeof(value_type));
    }
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    {
      size_type worksize = _worksize * (_nstreams + 1);
      if (size_t(worksize) > _work.span()) {
        if (_work.span() > 0) track_free(_work.span() * sizeof(value_type));
        Kokkos::resize(_work, worksize);
        track_alloc(_work.span() * sizeof(value_type));
      }
    }
#endif
    stat.t_extra = timer.seconds();

    timer.reset();
    {
      _ax = ax; // matrix values
      constexpr bool copy_to_l_buf(false);
      _info.copySparseToSuperpanels(copy_to_l_buf, _ap, _aj, _ax, shift, _perm, _peri);
      exec_space().fence(); // wait for copy
    }
    stat.t_copy = timer.seconds();

    stat_level.n_kernel_launching = 0;
    timer.reset();
    {
      // this should be considered with average problem sizes in levels
      const ordinal_type half_level = _nlevel / 2;
#if defined(CUDA_VERSION)
#if (11000 > CUDA_VERSION)
      /// cuda 11.1 below
      const ordinal_type team_size_factor[2] = {32, 64}, vector_size_factor[2] = {8, 4};
#else
      /// cuda 11.1 and higher
      const ordinal_type team_size_factor[2] = {64, 64}, vector_size_factor[2] = {8, 4};
#endif
#else
      /// not cuda
      const ordinal_type team_size_factor[2] = {64, 64}, vector_size_factor[2] = {8, 4};
#endif
      const ordinal_type team_size_update[2] = {16, 8},  vector_size_update[2] = {32, 32};
      // returned value from team LDL
      colind_view d_rval("rval",1);
      auto h_rval = Kokkos::create_mirror_view(host_memory_space(), d_rval);
      {
        typedef TeamFunctor_FactorizeLDL<supernode_info_type> functor_type;
        functor_type functor(_info, _factorize_mode, _level_sids, _piv, _diag, _buf, d_rval.data());

#if defined(TACHO_TEST_LEVELSET_TOOLS_KERNEL_OVERHEAD)
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space, typename functor_type::DummyTag>
            team_policy_factorize;
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space, typename functor_type::DummyTag>
            team_policy_update;
#else
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space,
                                   typename functor_type::template FactorizeTag<variant>>
            team_policy_factor;
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space, typename functor_type::UpdateTag>
            team_policy_update;
#endif

        // get max vector length
        team_policy_factor policy_factor(1, 1, 1);
        team_policy_update policy_update(1, 1, 1);
        const ordinal_type vmax = policy_factor.vector_length_max();
        {
          for (ordinal_type lvl = (_team_serial_level_cut - 1); lvl >= 0; --lvl) {
            const ordinal_type pbeg = _h_level_ptr(lvl), pend = _h_level_ptr(lvl + 1), pcnt = pend - pbeg;

            const range_type range_buf_factor_ptr(_h_buf_level_ptr(lvl), _h_buf_level_ptr(lvl + 1));

            const auto buf_factor_ptr = Kokkos::subview(_buf_factor_ptr, range_buf_factor_ptr);
            functor.setRange(pbeg, pend);
            functor.setBufferPtr(buf_factor_ptr);
            if (is_host) {
              policy_factor = team_policy_factor(pcnt, 1, 1);
              policy_update = team_policy_update(pcnt, 1, 1);
            } else {
              const ordinal_type idx = lvl > half_level;
              // pick vector sizes
              ordinal_type vsize_factor = std::min(vector_size_factor[idx],vmax);
              ordinal_type vsize_update = std::min(vector_size_update[idx],vmax);
              // pick teamm sizes
              policy_factor = team_policy_factor(pcnt, 1, vsize_factor);
              policy_update = team_policy_update(pcnt, 1, vsize_update);
              const ordinal_type factor_tmax = policy_factor.team_size_max(functor, Kokkos::ParallelForTag());
              const ordinal_type update_tmax = policy_update.team_size_max(functor, Kokkos::ParallelForTag());

              policy_factor = team_policy_factor(pcnt, std::min(team_size_factor[idx],factor_tmax), vsize_factor);
              policy_update = team_policy_update(pcnt, std::min(team_size_update[idx],update_tmax), vsize_update);
            }
            if (lvl < _device_level_cut) {
              // do nothing
              // Kokkos::parallel_for("factor lower", policy_factor, functor);
            } else {
              if (verbose) {
                Kokkos::fence(); tick.reset();
              }
              Kokkos::parallel_for("factor", policy_factor, functor);
              if (verbose) {
                Kokkos::fence(); time_parallel += tick.seconds();
              }
              ++stat_level.n_kernel_launching;
            }

            const auto h_buf_factor_ptr = Kokkos::subview(_h_buf_factor_ptr, range_buf_factor_ptr);

            if (verbose) {
              Kokkos::fence(); tick.reset();
            }
            factorizeLDL_OnDevice(pbeg, pend, h_buf_factor_ptr, _work);
            if (verbose) {
              Kokkos::fence(); time_device += tick.seconds();
              tick.reset();
            }
            Kokkos::deep_copy(h_rval, d_rval);
            int rval = h_rval(0);
            if (rval != 0) {
              TACHO_TEST_FOR_EXCEPTION(rval, std::runtime_error, "SYTRF (team) returns non-zero error code.");
            }
            Kokkos::fence(); // sync device calls

            Kokkos::parallel_for("update factor", policy_update, functor);
            if (verbose) {
              Kokkos::fence(); time_update += tick.seconds();
            }
            exec_space().fence();
            ++stat_level.n_kernel_launching;
          }
          // NOTE: device info not needed on host?
          //const auto exec_instance = exec_space();
          //Kokkos::deep_copy(exec_instance, _h_supernodes, _info.supernodes);
        }
      }
    } // end of LDL
    stat.t_factor = timer.seconds();
    timer.reset();
    if (variant == 3) {
      // compress each partitioned inverse at each level into CRS matrix
      extractCRS(store_transpose, verbose);
    }
    stat.t_extra += timer.seconds();

    if (verbose) {
      printf("Summary: LevelSetTools-Variant-%d (LDL Factorize)\n", variant);
      printf("=================================================\n");
      printf( "\n  ** Team = %f s, Device = %f s, Update = %f s **\n\n",time_parallel,time_device,time_update );
      print_stat_factor();
      fflush(stdout);
    }
  }

  inline void solveLDL(const value_type_matrix &x, // solution
                       const value_type_matrix &b, // right hand side
                       const value_type_matrix &t, // temporary workspace (store permuted vectors)
                       const ordinal_type verbose) {
    TACHO_TEST_FOR_EXCEPTION(x.extent(0) != b.extent(0) || x.extent(1) != b.extent(1) || x.extent(0) != t.extent(0) ||
                                 x.extent(1) != t.extent(1),
                             std::logic_error, "x, b, t, and w dimensions do not match");

    TACHO_TEST_FOR_EXCEPTION(x.data() == b.data() || x.data() == t.data(), std::logic_error,
                             "x, b, t, and w have the same data pointer");
    constexpr bool is_host = std::is_same<exec_memory_space, Kokkos::HostSpace>::value;

    // solve L D L^{H} x = b
    const ordinal_type nrhs = x.extent(1);
    Kokkos::Timer timer;

    stat_level.n_kernel_launching = 0;

    // one-time operation when nrhs is changed
    timer.reset();
    allocateWorkspaceSolve(nrhs);

    // execution spaces for non-device calls
    const bool need_fence = (!_team_on_user_stream || _nstreams > 1);
    const auto perm_exec_instance = _exec_instances[0];
    const auto team_exec_instance = (_team_on_user_stream ? _exec_instances[0] : exec_space());

    // 0. permute (from METIS) and copy b -> t
    ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(perm_exec_instance, b, _perm, t);
    if(variant != 3 && need_fence) perm_exec_instance.fence();
    stat.t_extra = timer.seconds();

    timer.reset();
    {
#if defined(TACHO_ENABLE_SOLVE_CHOLESKY_USE_LIGHT_KERNEL)
      const auto work_item_property = Kokkos::Experimental::WorkItemProperty::HintLightWeight;
#endif
      // this should be considered with average problem sizes in levels
      const ordinal_type half_level = _nlevel / 2;
#if defined(CUDA_VERSION)
      /// cuda
      const ordinal_type team_size_solve[2] = {16, 16}, vector_size_solve[2] = {8, 8};
#else
      /// not cuda
      const ordinal_type team_size_solve[2] = {64, 16}, vector_size_solve[2] = {8, 8};
#endif
      const ordinal_type team_size_update[2] = {128, 32}, vector_size_update[2] = {1, 1};
      {
        typedef TeamFunctor_SolveLowerLDL<supernode_info_type> functor_type;
#if defined(TACHO_TEST_SOLVE_CHOLESKY_KERNEL_OVERHEAD)
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space, typename functor_type::DummyTag>
            team_policy_solve;
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space, typename functor_type::DummyTag>
            team_policy_update;
#else
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space,
                                   typename functor_type::template SolveTag<variant>>
            team_policy_solve;
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space,
                                   typename functor_type::template UpdateTag<variant>>
            team_policy_update;
#endif
        functor_type functor(_info, _solve_mode, _level_sids, _piv, t, _buf);

        team_policy_solve policy_solve(1, 1, 1);
        team_policy_update policy_update(1, 1, 1);

        //  1. L^{-1}t = t
        {
          for (ordinal_type lvl = (_team_serial_level_cut - 1); lvl >= 0; --lvl) {
            const ordinal_type pbeg = _h_level_ptr(lvl), pend = _h_level_ptr(lvl + 1), pcnt = pend - pbeg;

            const range_type range_solve_buf_ptr(_h_buf_level_ptr(lvl), _h_buf_level_ptr(lvl + 1));

            const auto solve_buf_ptr = Kokkos::subview(_buf_solve_nrhs_ptr, range_solve_buf_ptr);
            functor.setRange(pbeg, pend);
            functor.setBufferPtr(solve_buf_ptr);
            if (is_host) {
              policy_solve = team_policy_solve(pcnt, 1, 1);
              policy_update = team_policy_update(pcnt, 1, 1);
            } else {
              const ordinal_type idx = lvl > half_level;
              policy_solve = team_policy_solve(team_exec_instance, pcnt, team_size_solve[idx], vector_size_solve[idx]);
              policy_update = team_policy_update(team_exec_instance, pcnt, team_size_update[idx], vector_size_update[idx]);
            }
#if defined(TACHO_ENABLE_SOLVE_CHOLESKY_USE_LIGHT_KERNEL)
            const auto policy_solve_with_work_property =
                Kokkos::Experimental::require(policy_solve, work_item_property);
            const auto policy_update_with_work_property =
                Kokkos::Experimental::require(policy_update, work_item_property);
#else
            const auto policy_solve_with_work_property = policy_solve;
            const auto policy_update_with_work_property = policy_update;
#endif
            if (variant != 3) {
              if (lvl < _device_level_cut) {
                // do nothing
                // Kokkos::parallel_for("solve lower", policy_solve, functor);
              } else {
                Kokkos::parallel_for("solve lower", policy_solve_with_work_property, functor);
                ++stat_level.n_kernel_launching;
              }
            }
            const auto h_buf_solve_ptr = Kokkos::subview(_h_buf_solve_nrhs_ptr, range_solve_buf_ptr);
            solveLDL_LowerOnDevice(lvl, _team_serial_level_cut, pbeg, pend, h_buf_solve_ptr, t);
            if (need_fence && _h_num_device_calls_solve(lvl) > 0)
              Kokkos::fence(); // fence solve on device before updating on default

            if (variant != 3) {
              Kokkos::parallel_for("update lower", policy_update_with_work_property, functor);
              ++stat_level.n_kernel_launching;
              if (need_fence && (lvl == 0 || _h_num_device_calls_solve(lvl-1) > 0))
                team_exec_instance.fence(); // synch update on default before next solve on device
            }
          }
        }
      } // end of lower tri solve

      {
        typedef TeamFunctor_SolveUpperLDL<supernode_info_type> functor_type;
#if defined(TACHO_TEST_SOLVE_CHOLESKY_KERNEL_OVERHEAD)
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space, typename functor_type::DummyTag>
            team_policy_solve;
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space, typename functor_type::DummyTag>
            team_policy_update;
#else
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space,
                                   typename functor_type::template SolveTag<variant>>
            team_policy_solve;
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space,
                                   typename functor_type::template UpdateTag<variant>>
            team_policy_update;
#endif
        functor_type functor(_info, _solve_mode, _level_sids, _piv, _diag, t, _buf);

        team_policy_solve policy_solve(1, 1, 1);
        team_policy_update policy_update(1, 1, 1);

        //  2. U^{-1} t = t;
        {
          for (ordinal_type lvl = 0; lvl < _team_serial_level_cut; ++lvl) {
            const ordinal_type pbeg = _h_level_ptr(lvl), pend = _h_level_ptr(lvl + 1), pcnt = pend - pbeg;

            const range_type range_solve_buf_ptr(_h_buf_level_ptr(lvl), _h_buf_level_ptr(lvl + 1));
            const auto solve_buf_ptr = Kokkos::subview(_buf_solve_nrhs_ptr, range_solve_buf_ptr);
            functor.setRange(pbeg, pend);
            functor.setBufferPtr(solve_buf_ptr);
            if (is_host) {
              policy_solve = team_policy_solve(pcnt, 1, 1);
              policy_update = team_policy_update(pcnt, 1, 1);
            } else {
              const ordinal_type idx = lvl > half_level;
              policy_solve = team_policy_solve(team_exec_instance, pcnt, team_size_solve[idx], vector_size_solve[idx]);
              policy_update = team_policy_update(team_exec_instance, pcnt, team_size_update[idx], vector_size_update[idx]);
            }
#if defined(TACHO_ENABLE_SOLVE_CHOLESKY_USE_LIGHT_KERNEL)
            const auto policy_solve_with_work_property =
                Kokkos::Experimental::require(policy_solve, work_item_property);
            const auto policy_update_with_work_property =
                Kokkos::Experimental::require(policy_update, work_item_property);
#else
            const auto policy_solve_with_work_property = policy_solve;
            const auto policy_update_with_work_property = policy_update;
#endif
            if (variant != 3) {
              Kokkos::parallel_for("update upper", policy_update_with_work_property, functor);
              ++stat_level.n_kernel_launching;
              if (need_fence && _h_num_device_calls_solve(lvl) > 0)
                team_exec_instance.fence(); // synch update befor solve on device

              if (lvl < _device_level_cut) {
                // do nothing
                // Kokkos::parallel_for("solve upper", policy_solve, functor);
              } else {
                Kokkos::parallel_for("solve upper", policy_solve_with_work_property, functor);
                ++stat_level.n_kernel_launching;
              }
            }

            const auto h_buf_solve_ptr = Kokkos::subview(_h_buf_solve_nrhs_ptr, range_solve_buf_ptr);
            solveLDL_UpperOnDevice(lvl, _team_serial_level_cut, pbeg, pend, h_buf_solve_ptr, t);
            if (need_fence && _h_num_device_calls_solve(lvl) > 0)
              Kokkos::fence(); // synch solve on device before next update
          }
        }
      } /// end of upper tri solve
    } // end of solve
    stat.t_solve = timer.seconds();

    // permute (from METIS) and copy t -> x
    if (variant != 3 && need_fence) Kokkos::fence(); // synch user or default streams
    timer.reset();
    ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(perm_exec_instance, t, _peri, x);
    perm_exec_instance.fence();
    stat.t_extra += timer.seconds();

    if (verbose) {
      printf("Summary: LevelSetTools-Variant-%d (LDL Solve: %3d)\n", variant, nrhs);
      printf("==================================================\n");
      print_stat_solve();
      fflush(stdout);
    }
  }

  /// 
  /// LU
  ///
  inline void factorizeLU(const value_type_array &ax, const bool store_transpose,
                          const mag_type shift, const mag_type pivot_tol, const ordinal_type verbose) {
    constexpr bool is_host = std::is_same<exec_memory_space, Kokkos::HostSpace>::value;
    Kokkos::Timer timer;
    Kokkos::Timer tick;
    double time_parallel = 0.0;
    double time_device = 0.0;
    double time_update = 0.0;

    timer.reset();
    if (_buf.span() < size_t(_bufsize_factorize)) {
      if (_buf.span() > 0) track_free(_buf.span() * sizeof(value_type));
      Kokkos::resize(_buf, _bufsize_factorize);
      track_alloc(_buf.span() * sizeof(value_type));
    }

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    {
      size_type worksize = _worksize * (_nstreams + 1);
      if (size_t(worksize) > _work.span()) {
        if (_work.span() > 0) track_free(_work.span() * sizeof(value_type));
        Kokkos::resize(_work, worksize);
        track_alloc(_work.span() * sizeof(value_type));
      }
    }
#endif
    stat.t_extra = timer.seconds();

    timer.reset();
    {
      _ax = ax; // matrix values
      constexpr bool copy_to_l_buf(true);
      _info.copySparseToSuperpanels(copy_to_l_buf, _ap, _aj, _ax, shift, _perm, _peri);
      exec_space().fence(); // wait for copy
    }
    stat.t_copy = timer.seconds();

    stat_level.n_kernel_launching = 0;
    timer.reset();
    {
      // this should be considered with average problem sizes in levels
      const ordinal_type half_level = _nlevel / 2;
#if defined(CUDA_VERSION)
#if (11000 > CUDA_VERSION)
      /// cuda 11.1 below
      const ordinal_type team_size_factor[2] = {32, 64}, vector_size_factor[2] = {8, 4};
#else
      /// cuda 11.1 and higher
      const ordinal_type team_size_factor[2] = {64, 64}, vector_size_factor[2] = {8, 4};
#endif
#else
      /// not cuda
      const ordinal_type team_size_factor[2] = {64, 64}, vector_size_factor[2] = {8, 4};
#endif
      const ordinal_type team_size_update[2] = {16, 8},  vector_size_update[2] = {32, 32};

      // returned value from team LU
      colind_view d_rval("rval",1);
      auto h_rval = Kokkos::create_mirror_view(host_memory_space(), d_rval);
      {
        typedef TeamFunctor_FactorizeLU<supernode_info_type> functor_type;
#if defined(TACHO_TEST_LEVELSET_TOOLS_KERNEL_OVERHEAD)
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space, typename functor_type::DummyTag>
            team_policy_factorize;
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space, typename functor_type::DummyTag>
            team_policy_update;
#else
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space,
                                   typename functor_type::template FactorizeTag<variant>>
            team_policy_factor;
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space, typename functor_type::UpdateTag>
            team_policy_update;
#endif
        team_policy_factor policy_factor(1, 1, 1);
        team_policy_update policy_update(1, 1, 1);
        functor_type functor(_info, _factorize_mode, _level_sids, _piv, _buf, d_rval.data());
        if (pivot_tol > 0.0) {
          functor.setDiagPertubationTol(pivot_tol);
        }
        // get max vector length
        const ordinal_type vmax = policy_factor.vector_length_max();
        {
          for (ordinal_type lvl = (_team_serial_level_cut - 1); lvl >= 0; --lvl) {
            const ordinal_type pbeg = _h_level_ptr(lvl), pend = _h_level_ptr(lvl + 1), pcnt = pend - pbeg;

            const range_type range_buf_factor_ptr(_h_buf_level_ptr(lvl), _h_buf_level_ptr(lvl + 1));
            const auto buf_factor_ptr = Kokkos::subview(_buf_factor_ptr, range_buf_factor_ptr);
            functor.setRange(pbeg, pend);
            functor.setBufferPtr(buf_factor_ptr);
            if (is_host) {
              policy_factor = team_policy_factor(pcnt, 1, 1);
              policy_update = team_policy_update(pcnt, 1, 1);
            } else {
              const ordinal_type idx = lvl > half_level;
              // pick vector sizes
              ordinal_type vsize_factor = std::min(vector_size_factor[idx],vmax);
              ordinal_type vsize_update = std::min(vector_size_update[idx],vmax);
              // pick teamm sizes
              policy_factor = team_policy_factor(pcnt, 1, vsize_factor);
              policy_update = team_policy_update(pcnt, 1, vsize_update);
              const ordinal_type factor_tmax = policy_factor.team_size_max(functor, Kokkos::ParallelForTag());
              const ordinal_type update_tmax = policy_update.team_size_max(functor, Kokkos::ParallelForTag());

              // create policies
              policy_factor = team_policy_factor(pcnt, std::min(team_size_factor[idx],factor_tmax), vsize_factor);
              policy_update = team_policy_update(pcnt, std::min(team_size_update[idx],update_tmax), vsize_update);
            }
            if (lvl < _device_level_cut) {
              // do nothing
              // Kokkos::parallel_for("factor lower", policy_factor, functor);
            } else {
              if (verbose) {
                Kokkos::fence(); tick.reset();
              }
              Kokkos::parallel_for("factor", policy_factor, functor);
              if (verbose) {
                Kokkos::fence(); time_parallel += tick.seconds();
              }
              ++stat_level.n_kernel_launching;
            }

            if (verbose) {
              Kokkos::fence(); tick.reset();
            }
            const auto h_buf_factor_ptr = Kokkos::subview(_h_buf_factor_ptr, range_buf_factor_ptr);
            factorizeLU_OnDevice(pbeg, pend, h_buf_factor_ptr, _work);
            if (verbose) {
              Kokkos::fence(); time_device += tick.seconds();
              tick.reset();
            }
            Kokkos::deep_copy(h_rval, d_rval);
            int rval = h_rval(0);
            if (rval != 0) {
              TACHO_TEST_FOR_EXCEPTION(rval, std::runtime_error, "GETRF (team) returns non-zero error code.");
            }
            if (_h_num_device_calls_factor(lvl) > 0)
              Kokkos::fence(); // sync device calls before update

            Kokkos::parallel_for("update factor", policy_update, functor);
            if (verbose) {
              Kokkos::fence(); time_update += tick.seconds();
            }
            if (lvl == 0 || _h_num_device_calls_factor(lvl-1) > 0)
              exec_space().fence(); // synch default before the next device factor call
            ++stat_level.n_kernel_launching;
          }
          // NOTE: device info not needed on host?
          //const auto exec_instance = exec_space();
          //Kokkos::deep_copy(exec_instance, _h_supernodes, _info.supernodes);
        }
      }
    } // end of LU
    stat.t_factor = timer.seconds();
    timer.reset();
    if (variant == 3) {
      // compress each partitioned inverse at each level into CRS matrix
      extractCRS(store_transpose, verbose);
    }
    stat.t_extra += timer.seconds();

    if (verbose) {
      printf("Summary: LevelSetTools-Variant-%d (LU Factorize)\n", variant);
      printf("================================================\n");
      printf( "\n  ** Team = %f s, Device = %f s, Update = %f s (%d streams) **\n",time_parallel,time_device,time_update,_nstreams );
      printf( "\n" );
      print_stat_factor();
      fflush(stdout);
    }
  }

  inline void solveLU(const value_type_matrix &x, // solution
                      const value_type_matrix &b, // right hand side
                      const value_type_matrix &t, // temporary workspace (store permuted vectors)
                      const ordinal_type verbose) {
    TACHO_TEST_FOR_EXCEPTION(x.extent(0) != b.extent(0) || x.extent(1) != b.extent(1) || x.extent(0) != t.extent(0) ||
                                 x.extent(1) != t.extent(1),
                             std::logic_error, "x, b, t, and w dimensions do not match");

    TACHO_TEST_FOR_EXCEPTION(x.data() == b.data() || x.data() == t.data(), std::logic_error,
                             "x, b, t, and w have the same data pointer");
    constexpr bool is_host = std::is_same<exec_memory_space, Kokkos::HostSpace>::value;

    // solve LU x = b
    const ordinal_type nrhs = x.extent(1);
    Kokkos::Timer timer;

    stat_level.n_kernel_launching = 0;

    // one-time operation when nrhs is changed
    timer.reset();
    allocateWorkspaceSolve(nrhs);

    // execution spaces for non-device calls
    const bool need_fence = (!_team_on_user_stream || _nstreams > 1);
    const auto perm_exec_instance = _exec_instances[0];
    const auto team_exec_instance = (_team_on_user_stream ? _exec_instances[0] : exec_space());

    // 0. permute (from METIS) and copy b -> t
    ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(perm_exec_instance, b, _perm, t);
    if (variant != 3 && need_fence) perm_exec_instance.fence();
    stat.t_extra = timer.seconds();

    timer.reset();
    {
#if defined(TACHO_ENABLE_SOLVE_CHOLESKY_USE_LIGHT_KERNEL)
      const auto work_item_property = Kokkos::Experimental::WorkItemProperty::HintLightWeight;
#endif
      // this should be considered with average problem sizes in levels
      const ordinal_type half_level = _nlevel / 2;
#if defined(CUDA_VERSION)
      /// cuda
      const ordinal_type team_size_solve[2] = {16, 16}, vector_size_solve[2] = {8, 8};
#else
      /// not cuda
      const ordinal_type team_size_solve[2] = {64, 16}, vector_size_solve[2] = {8, 8};
#endif
      const ordinal_type team_size_update[2] = {128, 32}, vector_size_update[2] = {1, 1};
      {
        typedef TeamFunctor_SolveLowerLU<supernode_info_type> functor_type;
#if defined(TACHO_TEST_SOLVE_CHOLESKY_KERNEL_OVERHEAD)
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space, typename functor_type::DummyTag>
            team_policy_solve;
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space, typename functor_type::DummyTag>
            team_policy_update;
#else
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space,
                                   typename functor_type::template SolveTag<variant>>
            team_policy_solve;
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space,
                                   typename functor_type::template UpdateTag<variant>>
            team_policy_update;
#endif
        functor_type functor(_info, _solve_mode, _level_sids, _piv, t, _buf);

        team_policy_solve policy_solve(1, 1, 1);
        team_policy_update policy_update(1, 1, 1);

        //  1. L w = t
        {
          for (ordinal_type lvl = (_team_serial_level_cut - 1); lvl >= 0; --lvl) {
            const ordinal_type pbeg = _h_level_ptr(lvl), pend = _h_level_ptr(lvl + 1), pcnt = pend - pbeg;
            const range_type range_solve_buf_ptr(_h_buf_level_ptr(lvl), _h_buf_level_ptr(lvl + 1));

            const auto solve_buf_ptr = Kokkos::subview(_buf_solve_nrhs_ptr, range_solve_buf_ptr);
            functor.setRange(pbeg, pend);
            functor.setBufferPtr(solve_buf_ptr);
            if (is_host) {
              policy_solve = team_policy_solve(pcnt, 1, 1);
              policy_update = team_policy_update(pcnt, 1, 1);
            } else {
              const ordinal_type idx = lvl > half_level;
              policy_solve = team_policy_solve(team_exec_instance, pcnt, team_size_solve[idx], vector_size_solve[idx]);
              policy_update = team_policy_update(team_exec_instance, pcnt, team_size_update[idx], vector_size_update[idx]);
            }
#if defined(TACHO_ENABLE_SOLVE_CHOLESKY_USE_LIGHT_KERNEL)
            const auto policy_solve_with_work_property =
                Kokkos::Experimental::require(policy_solve, work_item_property);
            const auto policy_update_with_work_property =
                Kokkos::Experimental::require(policy_update, work_item_property);
#else
            const auto policy_solve_with_work_property = policy_solve;
            const auto policy_update_with_work_property = policy_update;
#endif
            if (variant != 3) {
              if (lvl < _device_level_cut) {
                // do nothing
                // Kokkos::parallel_for("solve lower", policy_solve, functor);
              } else {
                Kokkos::parallel_for("solve lower", policy_solve_with_work_property, functor);
                ++stat_level.n_kernel_launching;
              }
            }
            const auto h_buf_solve_ptr = Kokkos::subview(_h_buf_solve_nrhs_ptr, range_solve_buf_ptr);
            solveLU_LowerOnDevice(lvl, _team_serial_level_cut, pbeg, pend, h_buf_solve_ptr, t);
            if (variant != 3) {
              if (need_fence && _h_num_device_calls_solve(lvl) > 0)
                Kokkos::fence(); // synch solve device before update

              Kokkos::parallel_for("update lower", policy_update_with_work_property, functor);
              ++stat_level.n_kernel_launching;
              if (need_fence && (lvl == 0 || _h_num_device_calls_solve(lvl-1) > 0))
                team_exec_instance.fence(); // synch update on default before next solve on device
            }
          }
        }
      } // end of lower tri solve

      {
        typedef TeamFunctor_SolveUpperLU<supernode_info_type> functor_type;
#if defined(TACHO_TEST_SOLVE_CHOLESKY_KERNEL_OVERHEAD)
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space, typename functor_type::DummyTag>
            team_policy_solve;
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space, typename functor_type::DummyTag>
            team_policy_update;
#else
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space,
                                   typename functor_type::template SolveTag<variant>>
            team_policy_solve;
        typedef Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space,
                                   typename functor_type::template UpdateTag<variant>>
            team_policy_update;
#endif
        functor_type functor(_info, _solve_mode, _level_sids, t, _buf);

        team_policy_solve policy_solve(1, 1, 1);
        team_policy_update policy_update(1, 1, 1);

        //  2. U t = w;
        {
          for (ordinal_type lvl = 0; lvl < _team_serial_level_cut; ++lvl) {
            const ordinal_type pbeg = _h_level_ptr(lvl), pend = _h_level_ptr(lvl + 1), pcnt = pend - pbeg;

            const range_type range_solve_buf_ptr(_h_buf_level_ptr(lvl), _h_buf_level_ptr(lvl + 1));
            const auto solve_buf_ptr = Kokkos::subview(_buf_solve_nrhs_ptr, range_solve_buf_ptr);
            functor.setRange(pbeg, pend);
            functor.setBufferPtr(solve_buf_ptr);
            if (is_host) {
              policy_solve = team_policy_solve(pcnt, 1, 1);
              policy_update = team_policy_update(pcnt, 1, 1);
            } else {
              const ordinal_type idx = lvl > half_level;
              policy_solve = team_policy_solve(team_exec_instance, pcnt, team_size_solve[idx], vector_size_solve[idx]);
              policy_update = team_policy_update(team_exec_instance, pcnt, team_size_update[idx], vector_size_update[idx]);
            }
#if defined(TACHO_ENABLE_SOLVE_CHOLESKY_USE_LIGHT_KERNEL)
            const auto policy_solve_with_work_property =
                Kokkos::Experimental::require(policy_solve, work_item_property);
            const auto policy_update_with_work_property =
                Kokkos::Experimental::require(policy_update, work_item_property);
#else
            const auto policy_solve_with_work_property = policy_solve;
            const auto policy_update_with_work_property = policy_update;
#endif
            if (variant != 3) {
              Kokkos::parallel_for("update upper", policy_update_with_work_property, functor);
              ++stat_level.n_kernel_launching;
              if (need_fence && _h_num_device_calls_solve(lvl) > 0)
                team_exec_instance.fence(); // synch update on default before solve on device

              if (lvl < _device_level_cut) {
                // do nothing
                // Kokkos::parallel_for("solve upper", policy_solve, functor);
              } else {
                Kokkos::parallel_for("solve upper", policy_solve_with_work_property, functor);
                ++stat_level.n_kernel_launching;
              }
            }
            const auto h_buf_solve_ptr = Kokkos::subview(_h_buf_solve_nrhs_ptr, range_solve_buf_ptr);
            solveLU_UpperOnDevice(lvl, _team_serial_level_cut, pbeg, pend, h_buf_solve_ptr, t);
            if (variant != 3) {
              if (need_fence && _h_num_device_calls_solve(lvl) > 0)
                Kokkos::fence(); // synch solve on device before the next update on default
              //else if (lvl == _team_serial_level_cut-1 || _h_num_device_calls_solve(lvl+1) > 0)
              //  exec_space().fence(); // synch bached solve-upper calls
            }
          }
        }
      } /// end of upper tri solve

    } // end of solve
    stat.t_solve = timer.seconds();

    // permute (from METIS) and copy t -> x
    if (variant != 3 && need_fence) Kokkos::fence();
    timer.reset();
    ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(perm_exec_instance, t, _peri, x);
    perm_exec_instance.fence();
    stat.t_extra += timer.seconds();

    if (verbose) {
      printf("Summary: LevelSetTools-Variant-%d (LU Solve: %3d)\n", variant, nrhs);
      printf("=================================================\n");
      print_stat_solve();
    }
  }

  /// 
  /// Main Entry Functions
  ///
  inline void factorize(const value_type_array &ax, const bool store_transpose,
                        const mag_type shift, const mag_type pivot_tol = 0.0,
                        const ordinal_type verbose = 0) override {
    Kokkos::deep_copy(_superpanel_buf, value_type(0));
    switch (this->getSolutionMethod()) {
    case 0:   /// LDL no-pivot
    case 1: { /// Cholesky
      factorizeCholesky(ax, store_transpose, shift, pivot_tol, verbose);
      break;
    }
    case 2: { /// LDL
      {
        const ordinal_type rlen = 4 * _m, plen = _piv.span();
        if (plen < rlen) {
          track_free(_piv.span() * sizeof(ordinal_type));
          _piv = ordinal_type_array("piv", rlen);
          track_alloc(_piv.span() * sizeof(ordinal_type));
        }
      }
      {
        const ordinal_type rlen = 2 * _m, dlen = _diag.span();
        if (dlen < rlen) {
          track_free(_diag.span() * sizeof(value_type));
          _diag = value_type_array("diag", rlen);
          track_alloc(_diag.span() * sizeof(value_type));
        }
      }
      factorizeLDL(ax, store_transpose, shift, verbose);
      break;
    }
    case 3: { /// LU
      {
        const ordinal_type rlen = 4 * _m, plen = _piv.span();
        if (plen < rlen) {
          track_free(_piv.span() * sizeof(ordinal_type));
          _piv = ordinal_type_array("piv", rlen);
          track_alloc(_piv.span() * sizeof(ordinal_type));
        }
      }
      factorizeLU(ax, store_transpose, shift, pivot_tol, verbose);
      break;
    }
    default: {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, "The solution method is not supported");
      break;
    }
    }
  }

  inline void solve(const value_type_matrix &x, // solution
                    const value_type_matrix &b, // right hand side
                    const value_type_matrix &t, // temporary workspace (store permuted vectors)
                    const ordinal_type verbose = 0) override {
    switch (this->getSolutionMethod()) {
    case 0:   /// LDL no-pivot
    case 1: { /// Cholesky
      solveCholesky(x, b, t, verbose);
      break;
    }
    case 2: { /// LDL
      solveLDL(x, b, t, verbose);
      break;
    }
    case 3: { /// LU
      solveLU(x, b, t, verbose);
      break;
    }
    default: {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error, "The solution method is not supported");
      break;
    }
    }
  }
};

} // namespace Tacho
#endif
