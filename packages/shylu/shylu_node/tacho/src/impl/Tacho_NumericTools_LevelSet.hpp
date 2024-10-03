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

#include "Tacho_Chol_OnDevice.hpp"
#include "Tacho_GemmTriangular_OnDevice.hpp"
#include "Tacho_Gemm_OnDevice.hpp"
#include "Tacho_Gemv_OnDevice.hpp"
#include "Tacho_Herk_OnDevice.hpp"
#include "Tacho_LDL_OnDevice.hpp"
#include "Tacho_LU_OnDevice.hpp"
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

#if (defined(KOKKOS_ENABLE_CUDA) && defined(TACHO_HAVE_CUSPARSE))
  // SpMV flag
  #if (CUSPARSE_VERSION >= 11400)
    #define TACHO_CUSPARSE_SPMV_ALG CUSPARSE_SPMV_ALG_DEFAULT
  #else
    #define TACHO_CUSPARSE_SPMV_ALG CUSPARSE_MV_ALG_DEFAULT
  #endif
  // SpMM flag
  #if (CUSPARSE_VERSION >= 11000)
    #define TACHO_CUSPARSE_SPMM_ALG CUSPARSE_SPMM_ALG_DEFAULT
  #else
    #define TACHO_CUSPARSE_SPMM_ALG CUSPARSE_MM_ALG_DEFAULT
  #endif
#elif defined(KOKKOS_ENABLE_HIP)
  #if (ROCM_VERSION >= 60000)
    #define tacho_rocsparse_spmv rocsparse_spmv
  #elif (ROCM_VERSION >= 50400)
    #define tacho_rocsparse_spmv rocsparse_spmv_ex
  #else
    #define tacho_rocsparse_spmv rocsparse_spmv
  #endif
#endif

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

  // for using SpMV
  rowptr_view rowptrU;
  colind_view colindU;
  nzvals_view nzvalsU;

  rowptr_view rowptrL;
  colind_view colindL;
  nzvals_view nzvalsL;

  // common for host and cuda
  int _status;

  // cuda stream
  int _nstreams;

  // workspace for SpMV
  bool _is_spmv_extracted;
  value_type_matrix _w_vec;
  value_type_array  buffer_U;
  value_type_array  buffer_L;
#if defined(KOKKOS_ENABLE_CUDA)
  bool _is_cublas_created, _is_cusolver_dn_created;
  cublasHandle_t _handle_blas;
  cusolverDnHandle_t _handle_lapack;
  #if defined(TACHO_HAVE_CUSPARSE)
  // workspace for SpMV
  // (separte for U and L, so that we can "destroy" without waiting for the other)
  cusparseDnMatDescr_t matL, matU, matW;
  cusparseDnVecDescr_t vecL, vecU, vecW;
  cusparseHandle_t cusparseHandle;
  #endif

  using blas_handle_type = cublasHandle_t;
  using lapack_handle_type = cusolverDnHandle_t;
  using stream_array_host = std::vector<cudaStream_t>;
  #define getBlasHandle(id)   _handle_blas
  #define getLapackHandle(id) _handle_lapack
#elif defined(KOKKOS_ENABLE_HIP)
  bool _is_rocblas_created;
  rocblas_handle _handle_blas;
  rocblas_handle _handle_lapack;
  std::vector<rocblas_handle> _handles;
  // workspace for SpMV
  rocsparse_dnmat_descr matL, matU, matW;
  rocsparse_dnvec_descr vecL, vecU, vecW;
  rocsparse_handle rocsparseHandle;

  using blas_handle_type = rocblas_handle;
  using lapack_handle_type = rocblas_handle;
  using stream_array_host = std::vector<hipStream_t>;
  #define getBlasHandle(id)   _handles[id]
  #define getLapackHandle(id) _handles[id]
#else
  int _handle_blas, _handle_lapack; // dummy handle for convenience
  using blas_handle_type = int;
  using lapack_handle_type = int;
  #define getBlasHandle()   _handle_blas
  #define getLapackHandle() _handle_lapack
#endif

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
  stream_array_host _streams;
  using exec_instance_array_host = std::vector<exec_space>;
  exec_instance_array_host _exec_instances;
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
    case 1: {
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
      TACHO_TEST_FOR_EXCEPTION(false, std::logic_error, "The solution method is not supported");
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
                         const ordinal_type device_solve_thres, const ordinal_type verbose = 0) {
    stat_level.n_device_factorize = 0;
    stat_level.n_device_solve = 0;
    stat_level.n_team_factorize = 0;
    stat_level.n_team_solve = 0;

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
      for (ordinal_type sid = 0; sid < _nsupernodes; ++sid)
        _nlevel = max(_stree_level(sid), _nlevel);
      ++_nlevel;
    }

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
          const ordinal_type chol_factor_work_size_variants[4] = {schur_work_size, max(m * m, schur_work_size),
                                                                  m * m + schur_work_size,
                                                                  m * m + schur_work_size};
          const ordinal_type chol_factor_work_size = chol_factor_work_size_variants[variant];
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

          const ordinal_type index_work_size = this->getSolutionMethod() - 1;
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
    _buf_solve_nrhs_ptr = Kokkos::create_mirror_view(exec_memory_space(), _h_buf_solve_nrhs_ptr);
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
      _status = rocblas_create_handle(&_handle_blas);
      checkDeviceBlasStatus("rocblasCreate");
      _is_rocblas_created = true;
      _handle_lapack = _handle_blas;
    }
#endif
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

    if (_device_level_cut > 0) {
      for (ordinal_type lvl = 0; lvl < _device_level_cut; ++lvl) {
        const ordinal_type pbeg = _h_level_ptr(lvl), pend = _h_level_ptr(lvl + 1);
        for (ordinal_type p = pbeg; p < pend; ++p) {
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
      for (ordinal_type lvl = _device_level_cut; lvl < _team_serial_level_cut; ++lvl) {
        const ordinal_type pbeg = _h_level_ptr(lvl), pend = _h_level_ptr(lvl + 1);
        for (ordinal_type p = pbeg; p < pend; ++p) {
          const ordinal_type sid = _h_level_sids(p);
          const auto s = _h_supernodes(sid);
          const ordinal_type m = s.m;    //, n_m = s.n-s.m;
          if (m > _device_solve_thres) { // || n > _device_solve_thres)
            _h_solve_mode(sid) = 0;
            ++stat_level.n_device_solve;
          } else {
            _h_solve_mode(sid) = 1;
            ++stat_level.n_team_solve;
          }
          if (m > _device_factorize_thres) { // || n_m > _device_factorize_thres)
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
    track_alloc(_factorize_mode.span() * sizeof(ordinal_type));

    _solve_mode = Kokkos::create_mirror_view_and_copy(exec_memory_space(), _h_solve_mode);
    track_alloc(_solve_mode.span() * sizeof(ordinal_type));

    stat.t_mode_classification = timer.seconds();
    if (verbose) {
      switch (this->getSolutionMethod()) {
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
      fflush(stdout);
    }
  }

  inline void release(const ordinal_type verbose = 0) override {
    base_type::release(false);
    if (variant == 3) {
      this->releaseCRS(true);
    }
    track_free(_buf_factor_ptr.span() * sizeof(size_type));
    track_free(_buf_solve_ptr.span() * sizeof(size_type));
    track_free(_buf_solve_nrhs_ptr.span() * sizeof(size_type));
    track_free(_buf.span() * sizeof(value_type));
    track_free(_factorize_mode.span() * sizeof(ordinal_type));
    track_free(_solve_mode.span() * sizeof(ordinal_type));
    track_free(_level_sids.span() * sizeof(ordinal_type));
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
    _nstreams = 0;
    _is_spmv_extracted = 0;
#if defined(KOKKOS_ENABLE_CUDA)
    _is_cublas_created = 0;
    _is_cusolver_dn_created = 0;
#endif
#if defined(KOKKOS_ENABLE_HIP)
    _is_rocblas_created = 0;
#endif
  }

  virtual ~NumericToolsLevelSet() {
#if defined(KOKKOS_ENABLE_CUDA)
    /// kokkos execution space may fence and it uses the wrapped stream when it is deallocated   
    /// on cuda, deallocting streams first does not cause any errors while hip generates errors.
    /// here, we just follow the consistent destruction process as hip does.
    _exec_instances.clear();

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
    _streams.clear();
#endif
#if defined(KOKKOS_ENABLE_HIP)
    /// kokkos execution space may fence and it uses the wrapped stream when it is deallocated   
    _exec_instances.clear();

    if (_is_rocblas_created) {
      _status = rocblas_destroy_handle(_handle_blas);
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
    _streams.clear();
#endif
  }

  inline void createStream(const ordinal_type nstreams, const ordinal_type verbose = 0) {
#if defined(KOKKOS_ENABLE_CUDA)
    // destroy previously created streams
    for (ordinal_type i = 0; i < _nstreams; ++i) {
      _status = cudaStreamDestroy(_streams[i]);
      checkDeviceStatus("cudaStreamDestroy");
    }
    // new streams
    _nstreams = nstreams;
    _streams.clear();
    _streams.resize(_nstreams);
    for (ordinal_type i = 0; i < _nstreams; ++i) {
      _status = cudaStreamCreateWithFlags(&_streams[i], cudaStreamNonBlocking);
      checkDeviceStatus("cudaStreamCreate");
    }
#endif
#if defined(KOKKOS_ENABLE_HIP)
    // destroy previously created streams
    for (ordinal_type i = 0; i < _nstreams; ++i) {
      _status = rocblas_destroy_handle(_handles[i]);
      checkDeviceLapackStatus("rocblasDestroy");

      _status = hipStreamDestroy(_streams[i]);
      checkDeviceStatus("hipStreamDestroy");
    }
    // new streams
    _nstreams = nstreams;
    _streams.clear();
    _streams.resize(_nstreams);
    _handles.resize(_nstreams);
    for (ordinal_type i = 0; i < _nstreams; ++i) {
      _status = rocblas_create_handle(&_handles[i]);
      checkDeviceStatus("rocblas_create_handle");
      _status = hipStreamCreateWithFlags(&_streams[i], hipStreamNonBlocking);
      checkDeviceStatus("hipStreamCreate");
    }
#endif
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    _exec_instances.clear();
    _exec_instances.resize(_nstreams);
    for (ordinal_type i = 0; i < _nstreams; ++i) {
      ExecSpaceFactory<exec_space>::createInstance(_streams[i], _exec_instances[i]);
    }

    if (verbose) {
      printf("Summary: CreateStream : %3d\n", _nstreams);
      printf("===========================\n");
      fflush(stdout);
    }
#endif
  }

  inline void setStreamOnHandle(const ordinal_type qid) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    // const ordinal_type qid = q % _nstreams;
    const auto mystream = _streams[qid];
#if defined(KOKKOS_ENABLE_CUDA)
    _status = cublasSetStream(_handle_blas, mystream);
    checkDeviceBlasStatus("cublasSetStream");

    _status = cusolverDnSetStream(_handle_lapack, mystream);
    checkDeviceLapackStatus("cusolverDnSetStream");
#endif
#if defined(KOKKOS_ENABLE_HIP)
    _status = rocblas_set_stream(_handle_blas, mystream);
    checkDeviceBlasStatus("rocblasSetStream(handle_blas)");
    _status = rocblas_set_stream(_handles[qid], mystream);
    checkDeviceBlasStatus("rocblasSetStream(handles[qid])");
#endif
#endif
  }

  ///
  /// Device level functions
  ///
  inline void factorizeCholeskyOnDeviceVar0(const ordinal_type pbeg, const ordinal_type pend,
                                            const size_type_array_host &h_buf_factor_ptr,
                                            const value_type_array &work) {
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
          const ordinal_type m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            value_type *aptr = s.u_buf;
            UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            aptr += m * m;
            _status = Chol<Uplo::Upper, Algo::OnDevice>::invoke(handle_lapack, ATL, W);
            checkDeviceLapackStatus("chol");

            if (n_m > 0) {
              UnmanagedViewType<value_type_matrix> ABR(_buf.data() + h_buf_factor_ptr(p - pbeg), n_m, n_m);
              UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // aptr += m*n_m;
              _status = Trsm<Side::Left, Uplo::Upper, Trans::ConjTranspose, Algo::OnDevice>::invoke(
                  handle_blas, Diag::NonUnit(), one, ATL, ATR);
              checkDeviceBlasStatus("trsm");

              _status = Herk<Uplo::Upper, Trans::ConjTranspose, Algo::OnDevice>::invoke(handle_blas, minus_one, ATR,
                                                                                        zero, ABR);
            }
          }
        }
      }
    }
  }

  inline void factorizeCholeskyOnDeviceVar1(const ordinal_type pbeg, const ordinal_type pend,
                                            const size_type_array_host &h_buf_factor_ptr,
                                            const value_type_array &work) {
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
          const ordinal_type m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            value_type *aptr = s.u_buf;
            UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            aptr += m * m;
            _status = Chol<Uplo::Upper, Algo::OnDevice>::invoke(handle_lapack, ATL, W);
            checkDeviceLapackStatus("chol");

            value_type *bptr = _buf.data() + h_buf_factor_ptr(p - pbeg);
            UnmanagedViewType<value_type_matrix> T(bptr, m, m);
            _status = SetIdentity<Algo::OnDevice>::invoke(exec_instance, T, one);
            checkDeviceBlasStatus("SetIdentity");

            _status = Trsm<Side::Left, Uplo::Upper, Trans::NoTranspose, Algo::OnDevice>::invoke(
                handle_blas, Diag::NonUnit(), one, ATL, T);
            checkDeviceBlasStatus("trsm");

            if (n_m > 0) {
              UnmanagedViewType<value_type_matrix> ABR(bptr, n_m, n_m);
              UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // aptr += m*n_m;
              _status = Trsm<Side::Left, Uplo::Upper, Trans::ConjTranspose, Algo::OnDevice>::invoke(
                  handle_blas, Diag::NonUnit(), one, ATL, ATR);
              checkDeviceBlasStatus("trsm");

              _status = Copy<Algo::OnDevice>::invoke(exec_instance, ATL, T);
              checkDeviceBlasStatus("Copy");

              _status = Herk<Uplo::Upper, Trans::ConjTranspose, Algo::OnDevice>::invoke(handle_blas, minus_one, ATR,
                                                                                        zero, ABR);
            } else {
              _status = Copy<Algo::OnDevice>::invoke(exec_instance, ATL, T);
              checkDeviceBlasStatus("Copy");
            }
          }
        }
      }
    }
  }

  inline void factorizeCholeskyOnDeviceVar2(const ordinal_type pbeg, const ordinal_type pend,
                                            const size_type_array_host &h_buf_factor_ptr,
                                            const value_type_array &work) {
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
          const ordinal_type m = s.m, n = s.n, n_m = n - m;
          if (m > 0) {
            value_type *aptr = s.u_buf;
            UnmanagedViewType<value_type_matrix> ATL(aptr, m, m);
            aptr += m * m;
            _status = Chol<Uplo::Upper, Algo::OnDevice>::invoke(handle_lapack, ATL, W);
            checkDeviceLapackStatus("chol");

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
          }
        }
      }
    }
  }

  inline void factorizeCholeskyOnDevice(const ordinal_type pbeg, const ordinal_type pend,
                                        const size_type_array_host &h_buf_factor_ptr, const value_type_array &work) {
    if (variant == 0)
      factorizeCholeskyOnDeviceVar0(pbeg, pend, h_buf_factor_ptr, work);
    else if (variant == 1)
      factorizeCholeskyOnDeviceVar1(pbeg, pend, h_buf_factor_ptr, work);
    else if (variant == 2 || variant == 3)
      factorizeCholeskyOnDeviceVar2(pbeg, pend, h_buf_factor_ptr, work);
    else {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "LevelSetTools::factorizeCholeskyOnDevice, algorithm variant is not supported");
    }
  }

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

              ConstUnmanagedViewType<ordinal_type_array> perm(P.data() + 2 * m, m);
              _status = ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(exec_instance, ATR,
                                                                                                 perm, STR);

              _status = Trsm<Side::Left, Uplo::Lower, Trans::NoTranspose, Algo::OnDevice>::invoke(
                  handle_blas, Diag::Unit(), one, ATL, STR);
              checkDeviceBlasStatus("trsm");

              _status = Copy<Algo::OnDevice>::invoke(exec_instance, T, ATL);
              _status = Copy<Algo::OnDevice>::invoke(exec_instance, ATR, STR);

              _status = Symmetrize<Uplo::Lower, Algo::OnDevice>::invoke(exec_instance, T);
              _status = SetIdentity<Algo::OnDevice>::invoke(exec_instance, ATL, minus_one);
              _status = Scale2x2_BlockInverseDiagonals<Side::Left, Algo::OnDevice>::invoke(exec_instance, P, D, ATR);

              _status = GemmTriangular<Trans::Transpose, Trans::NoTranspose, Uplo::Upper, Algo::OnDevice>::invoke(
                  handle_blas, minus_one, ATR, STR, zero, ABR);
              checkDeviceBlasStatus("gemm");

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
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "LevelSetTools::factorizeLDL_OnDevice, algorithm variant is not supported");
    }
  }

  inline void factorizeLU_OnDeviceVar0(const ordinal_type pbeg, const ordinal_type pend,
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
          }
        }
      }
    }
  }

  inline void factorizeLU_OnDeviceVar1(const ordinal_type pbeg, const ordinal_type pend,
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
          }
        }
      }
    }
  }

  inline void factorizeLU_OnDeviceVar2(const ordinal_type pbeg, const ordinal_type pend,
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
          }
        }
      }
    }
  }

  inline void factorizeLU_OnDevice(const ordinal_type pbeg, const ordinal_type pend,
                                   const size_type_array_host &h_buf_factor_ptr, const value_type_array &work) {
    if (variant == 0)
      factorizeLU_OnDeviceVar0(pbeg, pend, h_buf_factor_ptr, work);
    else if (variant == 1)
      factorizeLU_OnDeviceVar1(pbeg, pend, h_buf_factor_ptr, work);
    else if (variant == 2 || variant == 3)
      factorizeLU_OnDeviceVar2(pbeg, pend, h_buf_factor_ptr, work);
    else {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "LevelSetTools::factorizeLU_OnDevice, algorithm variant is not supported");
    }
  }

  inline void extractCRS(bool lu) {
#if (defined(KOKKOS_ENABLE_CUDA) && defined(TACHO_HAVE_CUSPARSE)) || \
     defined(KOKKOS_ENABLE_HIP)

    const ordinal_type nrhs = 1;
    const ordinal_type m = _m;
    const value_type one(1);
    const value_type zero(0);

    // ========================
    // free CRS, 
    // if it has been extracted
#if defined(KOKKOS_ENABLE_HIP)
    this->releaseCRS(!lu);
#else
    this->releaseCRS(true);
#endif

    // ========================
    // workspace
    Kokkos::resize(_w_vec, m, nrhs);
    int ldw = _w_vec.stride(1);
#if defined(KOKKOS_ENABLE_CUDA)
    cudaDataType computeType;
    if (std::is_same<value_type, double>::value) {
      computeType = CUDA_R_64F;
    } else if (std::is_same<value_type, float>::value) {
      computeType = CUDA_R_32F;
    } else {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "LevelSetTools::solveCholeskyLowerOnDevice: ComputeSPMV only supported double or float");
    }
    // create cusparse handle
    cusparseCreate(&cusparseHandle);
    // attach to Cusparse data struct
    cusparseCreateDnMat(&matW, m, nrhs, ldw, (void*)(_w_vec.data()), computeType, CUSPARSE_ORDER_COL);
    cusparseCreateDnVec(&vecW, m, (void*)(_w_vec.data()), computeType);
    // also to T, to be destroyed before each SpMV call
    cusparseCreateDnMat(&matL, m, nrhs, ldw, (void*)(_w_vec.data()), computeType, CUSPARSE_ORDER_COL);
    cusparseCreateDnVec(&vecL, m, (void*)(_w_vec.data()), computeType);
    cusparseCreateDnMat(&matU, m, nrhs, ldw, (void*)(_w_vec.data()), computeType, CUSPARSE_ORDER_COL);
    cusparseCreateDnVec(&vecU, m, (void*)(_w_vec.data()), computeType);
    // vectors used for preprocessing
#ifdef USE_SPMM_FOR_WORKSPACE_SIZE
    cusparseDnMatDescr_t vecX, vecY;
    const ordinal_type ldx = _w_vec.stride(1);
    cusparseCreateDnMat(&vecX, m, nrhs, ldx, _w_vec.data(), computeType, CUSPARSE_ORDER_COL);
    cusparseCreateDnMat(&vecY, m, nrhs, ldx, _w_vec.data(), computeType, CUSPARSE_ORDER_COL);
#else
    cusparseDnVecDescr_t vecX, vecY;
    cusparseCreateDnVec(&vecX, m, _w_vec.data(), computeType);
    cusparseCreateDnVec(&vecY, m, _w_vec.data(), computeType);
#endif
#elif defined(KOKKOS_ENABLE_HIP)
    rocsparse_datatype rocsparse_compute_type = rocsparse_datatype_f64_r;
    if (std::is_same<value_type, float>::value) {
      rocsparse_compute_type = rocsparse_datatype_f32_r;
    }
    // create rocsparse handle
    rocsparse_create_handle(&rocsparseHandle);
    // attach to Rocsparse data struct
    rocsparse_create_dnmat_descr(&matW, m, nrhs, ldw, (void*)(_w_vec.data()), rocsparse_compute_type, rocsparse_order_column);
    rocsparse_create_dnvec_descr(&vecW, m, (void*)(_w_vec.data()), rocsparse_compute_type);
    // also to T, to be destroyed before each SpMV call
    rocsparse_create_dnmat_descr(&matL, m, nrhs, ldw, (void*)(_w_vec.data()), rocsparse_compute_type, rocsparse_order_column);
    rocsparse_create_dnvec_descr(&vecL, m, (void*)(_w_vec.data()), rocsparse_compute_type);
    rocsparse_create_dnmat_descr(&matU, m, nrhs, ldw, (void*)(_w_vec.data()), rocsparse_compute_type, rocsparse_order_column);
    rocsparse_create_dnvec_descr(&vecU, m, (void*)(_w_vec.data()), rocsparse_compute_type);
    // vectors used for preprocessing
    rocsparse_dnvec_descr vecX, vecY;
    rocsparse_create_dnvec_descr(&vecX, m, (void*)_w_vec.data(), rocsparse_compute_type);
    rocsparse_create_dnvec_descr(&vecY, m, (void*)_w_vec.data(), rocsparse_compute_type);
#endif

    // allocate rowptrs
    Kokkos::resize(rowptrU, _team_serial_level_cut*(1+m));
    Kokkos::resize(rowptrL, _team_serial_level_cut*(1+m));
    Kokkos::deep_copy(rowptrL, 0);
    // counting nnz, first, so that we can allocate in NumericalTool
    size_t ptr = 0;
    size_t nnzU = 0;
    size_t nnzL = 0;
    typedef TeamFunctor_ExtractCrs<supernode_info_type> functor_type;
    for (ordinal_type lvl = 0; lvl < _team_serial_level_cut; ++lvl) {
      const ordinal_type pbeg = _h_level_ptr(lvl), pend = _h_level_ptr(lvl + 1);

      // the first supernode in this lvl (where the CRS matrix is stored)
      auto &s0 = _h_supernodes(_h_level_sids(pbeg));
#if defined(KOKKOS_ENABLE_HIP)
      s0.spmv_explicit_transpose = true;
#else
      s0.spmv_explicit_transpose = false; // okay for SpMV, though may not for SpMM
#endif

      #define TACHO_INSERT_DIAGONALS
      // NOTE: this needs extra vector-entry copy for the non-active rows at each level for solve (copy t to w, and w back to t)
      //       but it seems faster on AMD 250X GPU, and not much performance impact on V100
      #define TACHO_INSERT_DIAGONALS
      // ========================
      // count nnz / row
      auto d_rowptrU = Kokkos::subview(rowptrU, range_type(ptr, ptr+m+1));
      s0.rowptrU = d_rowptrU.data();

      functor_type extractor_crs(_info, _solve_mode, _level_sids);
      extractor_crs.setGlobalSize(m);
      extractor_crs.setRange(pbeg, pend);
      extractor_crs.setRowPtr(s0.rowptrU);
      {
        using team_policy_type = Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space,
                                                    typename functor_type::ExtractPtrTag>;
        team_policy_type team_policy((pend-pbeg)+1, Kokkos::AUTO());

        Kokkos::parallel_for("extract rowptr", team_policy, extractor_crs);
        exec_space().fence();
      }

      // ========================
      // shift to generate rowptr
      using range_type = Kokkos::pair<int, int>;
      {
        using range_policy_type = Kokkos::RangePolicy<exec_space>;
        Kokkos::parallel_scan("shiftRowptr", range_policy_type(0, m+1), rowptr_sum(s0.rowptrU));
        exec_space().fence();
        // get nnz
        auto d_nnz = Kokkos::subview(d_rowptrU, range_type(m, m+1));
        auto h_nnz = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_nnz);
        s0.nnzU = h_nnz(0);
        nnzU += s0.nnzU;
      }

      if (lu) {
        // get nnz per row (L is stored by column)
        auto d_rowptrL = Kokkos::subview(rowptrL, range_type(ptr, ptr+m+1));
        s0.rowptrL = d_rowptrL.data();
        {
          using team_policy_type = Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space,
                                                      typename functor_type::ExtractPtrColTag>;
          team_policy_type team_policy((pend-pbeg)+1, Kokkos::AUTO());

          extractor_crs.setRowPtr(s0.rowptrL);
          Kokkos::parallel_for("extract rowptr L", team_policy, extractor_crs);
          exec_space().fence();
        }
        {
          // convert to offset
          using range_policy_type = Kokkos::RangePolicy<exec_space>;
          Kokkos::parallel_scan("shiftRowptr L", range_policy_type(0, m+1), rowptr_sum(s0.rowptrL));
          exec_space().fence();
          // get nnz (on CPU for now)
          auto d_nnz = Kokkos::subview(d_rowptrL, range_type(m, m+1));
          auto h_nnz = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_nnz);
          s0.nnzL = h_nnz(0);
          nnzL += s0.nnzL;
        }
        s0.spmv_explicit_transpose = true;
      } else if (s0.spmv_explicit_transpose) {
        // ========================
        // explicitly form transpose
        s0.nnzL = s0.nnzU;
        auto d_rowptrL = Kokkos::subview(rowptrL, range_type(ptr, ptr+m+1));
        s0.rowptrL = d_rowptrL.data();
        nnzL += s0.nnzL;
      }
      ptr += (1+m);
    }

    // allocate (TODO: move to symbolic)
    if (nnzU) {
      Kokkos::resize(colindU, nnzU);
      Kokkos::resize(nzvalsU, nnzU);
    } 
    if (nnzL) {
      Kokkos::resize(colindL, nnzL);
      Kokkos::resize(nzvalsL, nnzL);
    }

    // load nonzero val/ind
    ptr = 0;
    nnzU = 0;
    nnzL = 0;
    for (ordinal_type lvl = 0; lvl < _team_serial_level_cut; ++lvl) {
      const ordinal_type pbeg = _h_level_ptr(lvl), pend = _h_level_ptr(lvl + 1);

      // the first supernode in this lvl (where the CRS matrix is stored)
      auto &s0 = _h_supernodes(_h_level_sids(pbeg));

      // ========================
      // assign memory
      auto d_rowptrU = Kokkos::subview(rowptrU, range_type(ptr, ptr+m+1));
      auto d_colindU = Kokkos::subview(colindU, range_type(nnzU, nnzU+s0.nnzU));
      auto d_nzvalsU = Kokkos::subview(nzvalsU, range_type(nnzU, nnzU+s0.nnzU));
      s0.colindU = d_colindU.data();
      s0.nzvalsU = d_nzvalsU.data();
      nnzU += s0.nnzU;

      // ========================
      // extract nonzero element
      functor_type extractor_crs(_info, _solve_mode, _level_sids);
      extractor_crs.setGlobalSize(m);
      extractor_crs.setRange(pbeg, pend);
      extractor_crs.setRowPtr(s0.rowptrU);
      extractor_crs.setCrsView(s0.colindU, s0.nzvalsU);
      {
        using team_policy_type = Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space,
                                                    typename functor_type::ExtractValTag>;
        team_policy_type team_policy((pend-pbeg)+1, Kokkos::AUTO());

        // >> launch functor to extract nonzero entries
        Kokkos::parallel_for("extract nzvals", team_policy, extractor_crs);
        exec_space().fence();
      }

      // ========================
      // shift back (TODO: shift first to avoid this)
      {
        //  copy to CPU, for now
        auto h_rowptr = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_rowptrU);
        for (ordinal_type i = m; i > 0 ; i--) h_rowptr(i) = h_rowptr(i-1);
        h_rowptr(0) = 0;
        Kokkos::deep_copy(d_rowptrU, h_rowptr);
      }

      if (lu) {
        auto d_rowptrL = Kokkos::subview(rowptrL, range_type(ptr, ptr+m+1));
        auto d_colindL = Kokkos::subview(colindL, range_type(nnzL, nnzL+s0.nnzL));
        auto d_nzvalsL = Kokkos::subview(nzvalsL, range_type(nnzL, nnzL+s0.nnzL));
        s0.colindL = d_colindL.data();
        s0.nzvalsL = d_nzvalsL.data();
        nnzL += s0.nnzL;

        // ========================
        // insert nonzeros
        extractor_crs.setRowPtr(s0.rowptrL);
        extractor_crs.setCrsView(s0.colindL, s0.nzvalsL);
        extractor_crs.setPivPtr(_piv);
        {
          using team_policy_type = Kokkos::TeamPolicy<Kokkos::Schedule<Kokkos::Static>, exec_space,
                                                      typename functor_type::ExtractValColTag>;
          team_policy_type team_policy((pend-pbeg)+1, Kokkos::AUTO());

          // >> launch functor to extract nonzero entries
          Kokkos::parallel_for("extract nzvals L", team_policy, extractor_crs);
          exec_space().fence();
        }
        // ========================
        // shift back
        // (TODO: shift first to avoid this)
        {
          //  copy to CPU, for now
          auto h_rowptr = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_rowptrL);
          for (ordinal_type i = m; i > 0 ; i--) h_rowptr(i) = h_rowptr(i-1);
          h_rowptr(0) = 0;
          Kokkos::deep_copy(d_rowptrL, h_rowptr);
        }
      } else if (s0.spmv_explicit_transpose) {
        // ========================
        // transpose
        // >> generate rowptr
        extractor_crs.setRowPtrT(s0.rowptrL);
        {
          // >> count nnz / row (transpose)
          using team_policy_type = Kokkos::RangePolicy<typename functor_type::TransPtrTag, exec_space>;
          team_policy_type team_policy(0, m);
          Kokkos::parallel_for("transpose pointer", team_policy, extractor_crs);
        }
        {
          // >> accumulate to generate rowptr (transpose)
          using range_policy_type = Kokkos::RangePolicy<exec_space>;
          Kokkos::parallel_scan("shiftRowptrT", range_policy_type(0, m+1), rowptr_sum(s0.rowptrL));
          exec_space().fence();
        }

        s0.nnzL = s0.nnzU;
        auto d_colindL = Kokkos::subview(colindL, range_type(nnzL, nnzL+s0.nnzL));
        auto d_nzvalsL = Kokkos::subview(nzvalsL, range_type(nnzL, nnzL+s0.nnzL));
        s0.colindL = d_colindL.data();
        s0.nzvalsL = d_nzvalsL.data();
        nnzL += s0.nnzL;
 
        // ========================
        // >> copy into transpose-matrix
        extractor_crs.setRowPtrT(s0.rowptrL);
        extractor_crs.setCrsViewT(s0.colindL, s0.nzvalsL);
        {
          using team_policy_type = Kokkos::RangePolicy<typename functor_type::TransMatTag, exec_space>;
          team_policy_type team_policy(0, m);
          Kokkos::parallel_for("transpose pointer", team_policy, extractor_crs);
        }
        // ========================
        // shift back
        // (TODO: shift first to avoid this)
        {
          // copy to CPU, for now
          auto d_rowptrL = Kokkos::subview(rowptrL, range_type(ptr, ptr+m+1));
          auto h_rowptr = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_rowptrL);
          for (ordinal_type i = m; i > 0 ; i--) h_rowptr(i) = h_rowptr(i-1);
          h_rowptr(0) = 0;
          Kokkos::deep_copy(d_rowptrL, h_rowptr);
        }
      }
      ptr += (1+m);

      // ========================
      // create NVIDIA/AMD data structures for SpMV
      size_t buffer_size_L = 0;
      size_t buffer_size_U = 0;
      value_type alpha = one;
      #ifdef TACHO_INSERT_DIAGONALS
      value_type beta = zero;
      #else
      value_type beta = one;
      #endif
#if defined(KOKKOS_ENABLE_CUDA)
      // create matrix
      cusparseCreateCsr(&s0.U_cusparse, m, m, s0.nnzU,
                        s0.rowptrU, s0.colindU, s0.nzvalsU,
                        CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                        CUSPARSE_INDEX_BASE_ZERO, computeType);

#ifdef USE_SPMM_FOR_WORKSPACE_SIZE
      cusparseSpMM_bufferSize(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE,
                              &alpha, s0.U_cusparse, vecX, &beta, vecY,
                              computeType, TACHO_CUSPARSE_SPMM_ALG, &buffer_size_U);
#else
      cusparseSpMV_bufferSize(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, s0.U_cusparse, vecX, &beta, vecY,
                              computeType, TACHO_CUSPARSE_SPMV_ALG, &buffer_size_U);
#endif
      if (s0.spmv_explicit_transpose) {
        // create matrix (transpose(U) or L)
        cusparseCreateCsr(&s0.L_cusparse, m, m, s0.nnzL,
                          s0.rowptrL, s0.colindL, s0.nzvalsL,
                          CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                          CUSPARSE_INDEX_BASE_ZERO, computeType);
        // workspace size
#ifdef USE_SPMM_FOR_WORKSPACE_SIZE
        cusparseSpMM_bufferSize(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                &alpha, s0.L_cusparse, vecX, &beta, vecY,
                                computeType, TACHO_CUSPARSE_SPMM_ALG, &buffer_size_L);
#else
        cusparseSpMV_bufferSize(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, s0.L_cusparse, vecX, &beta, vecY,
                                computeType, TACHO_CUSPARSE_SPMV_ALG, &buffer_size_L);
#endif
      } else {
        // create matrix (L_cusparse stores the same ptrs as descrU, but optimized for trans)
        s0.nnzL = s0.nnzU;
        cusparseCreateCsr(&s0.L_cusparse, m, m, s0.nnzL,
                          s0.rowptrU, s0.colindU, s0.nzvalsU,
                          CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                          CUSPARSE_INDEX_BASE_ZERO, computeType);
        // workspace size for transpose SpMV
#ifdef USE_SPMM_FOR_WORKSPACE_SIZE
        cusparseSpMM_bufferSize(cusparseHandle, CUSPARSE_OPERATION_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                &alpha, s0.L_cusparse, vecX, &beta, vecY,
                                computeType, TACHO_CUSPARSE_SPMM_ALG, &buffer_size_L);
#else
        cusparseSpMV_bufferSize(cusparseHandle, CUSPARSE_OPERATION_TRANSPOSE, &alpha, s0.L_cusparse, vecX, &beta, vecY,
                                computeType, TACHO_CUSPARSE_SPMV_ALG, &buffer_size_L);
#endif
      }
      // allocate workspace
      if (buffer_size_U > buffer_U.extent(0)) {
        Kokkos::resize(buffer_U, buffer_size_U);
      }
      if (buffer_size_L > buffer_L.extent(0)) {
        Kokkos::resize(buffer_L, buffer_size_L);
      }
#elif defined(KOKKOS_ENABLE_HIP)
      // create matrix
      rocsparse_create_csr_descr(&(s0.descrU), m, m, s0.nnzU,
                                 s0.rowptrU, s0.colindU, s0.nzvalsU,
                                 rocsparse_indextype_i32, rocsparse_indextype_i32, rocsparse_index_base_zero, rocsparse_compute_type);
      // workspace
      tacho_rocsparse_spmv
          (rocsparseHandle, rocsparse_operation_none,
           &alpha, s0.descrU, vecX, &beta, vecY,
           rocsparse_compute_type, rocsparse_spmv_alg_default,
           #if ROCM_VERSION >= 50400
           rocsparse_spmv_stage_buffer_size,
           #endif
           &buffer_size_U, nullptr);
      // allocate workspace
      if (buffer_size_U > buffer_U.extent(0)) {
        Kokkos::resize(buffer_U, buffer_size_U);
      }
      #if ROCM_VERSION >= 50400
      // preprocess
      buffer_size_U = buffer_U.extent(0);
      tacho_rocsparse_spmv
          (rocsparseHandle, rocsparse_operation_none,
           &alpha, s0.descrU, vecX, &beta, vecY,
           rocsparse_compute_type, rocsparse_spmv_alg_default,
           rocsparse_spmv_stage_preprocess,
           &buffer_size_U, (void*)buffer_U.data());
      #endif
      if (s0.spmv_explicit_transpose) {
        // create matrix (transpose)
        rocsparse_create_csr_descr(&(s0.descrL), m, m, s0.nnzL,
                                   s0.rowptrL, s0.colindL, s0.nzvalsL,
                                   rocsparse_indextype_i32, rocsparse_indextype_i32, rocsparse_index_base_zero, rocsparse_compute_type);
        // workspace
        tacho_rocsparse_spmv
          (rocsparseHandle, rocsparse_operation_none,
           &alpha, s0.descrL, vecX, &beta, vecY,
           rocsparse_compute_type, rocsparse_spmv_alg_default,
           #if ROCM_VERSION >= 50400
           rocsparse_spmv_stage_buffer_size,
           #endif
           &buffer_size_L, nullptr);
        // allocate workspace
        if (buffer_size_L > buffer_L.extent(0)) {
          Kokkos::resize(buffer_L, buffer_size_L);
        }
        #if ROCM_VERSION >= 50400
        // preprocess
        buffer_size_L = buffer_L.extent(0);
        tacho_rocsparse_spmv
          (rocsparseHandle, rocsparse_operation_none,
           &alpha, s0.descrL, vecX, &beta, vecY,
           rocsparse_compute_type, rocsparse_spmv_alg_default,
            rocsparse_spmv_stage_preprocess,
           &buffer_size_L, (void*)buffer_L.data());
        #endif
      } else {
        // create matrix, transpose (L_cusparse stores the same ptrs as descrU, but optimized for trans)
        rocsparse_create_csr_descr(&(s0.descrL), m, m, s0.nnzL,
                                   s0.rowptrU, s0.colindU, s0.nzvalsU,
                                   rocsparse_indextype_i32, rocsparse_indextype_i32, rocsparse_index_base_zero, rocsparse_compute_type);
        // workspace (transpose)
        tacho_rocsparse_spmv
          (rocsparseHandle, rocsparse_operation_transpose,
           &alpha, s0.descrL, vecX, &beta, vecY,
           rocsparse_compute_type, rocsparse_spmv_alg_default,
           #if ROCM_VERSION >= 50400
           rocsparse_spmv_stage_buffer_size,
           #endif
           &buffer_size_L, nullptr);
        // allcate workspace
        if (buffer_size_L > buffer_L.extent(0)) {
          Kokkos::resize(buffer_L, buffer_size_L);
        }
        #if ROCM_VERSION >= 50400
        // preprocess
        buffer_size_L = buffer_L.extent(0);
        tacho_rocsparse_spmv
          (rocsparseHandle, rocsparse_operation_transpose,
           &alpha, s0.descrL, vecX, &beta, vecY,
           rocsparse_compute_type, rocsparse_spmv_alg_default,
            rocsparse_spmv_stage_preprocess,
           &buffer_size_L, (void*)buffer_L.data());
        #endif
      } 
#endif
    }

#if defined(KOKKOS_ENABLE_CUDA)
#ifdef USE_SPMM_FOR_WORKSPACE_SIZE
    cusparseDestroyDnMat(vecX);
    cusparseDestroyDnMat(vecY);
#else
    cusparseDestroyDnVec(vecX);
    cusparseDestroyDnVec(vecY);
#endif
#elif defined(KOKKOS_ENABLE_HIP)
    rocsparse_destroy_dnvec_descr(vecX);
    rocsparse_destroy_dnvec_descr(vecY);
#endif
    _is_spmv_extracted = 1;
#endif
  }

  inline void releaseCRS(bool release_all) {
    if(_is_spmv_extracted) {
      Kokkos::fence();
      if (release_all) {
        for (ordinal_type lvl = 0; lvl < _team_serial_level_cut; ++lvl) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
          const ordinal_type pbeg = _h_level_ptr(lvl);
#endif
          // the first supernode in this lvl (where the CRS matrix is stored)
#if defined(KOKKOS_ENABLE_CUDA)
          auto &s0 = _h_supernodes(_h_level_sids(pbeg));
          cusparseDestroySpMat(s0.U_cusparse);
          cusparseDestroySpMat(s0.L_cusparse);
#elif defined(KOKKOS_ENABLE_HIP)
          auto &s0 = _h_supernodes(_h_level_sids(pbeg));
          rocsparse_destroy_spmat_descr(s0.descrU);
          rocsparse_destroy_spmat_descr(s0.descrL);
#endif
        }
      }
#if defined(TACHO_HAVE_CUSPARSE) && defined(KOKKOS_ENABLE_CUDA)
      cusparseDestroy(cusparseHandle);
      cusparseDestroyDnMat(matL);
      cusparseDestroyDnVec(vecL);
      cusparseDestroyDnMat(matU);
      cusparseDestroyDnVec(vecU);
      cusparseDestroyDnMat(matW);
      cusparseDestroyDnVec(vecW); 
#elif defined(KOKKOS_ENABLE_HIP)
      rocsparse_destroy_handle(rocsparseHandle);
      rocsparse_destroy_dnmat_descr(matL);
      rocsparse_destroy_dnvec_descr(vecL);
      rocsparse_destroy_dnmat_descr(matU);
      rocsparse_destroy_dnvec_descr(vecU);
      rocsparse_destroy_dnmat_descr(matW);
      rocsparse_destroy_dnvec_descr(vecW); 
#endif
      _is_spmv_extracted = 0;
    }
  }

  ///
  /// Level set factorize
  ///
  inline void factorizeCholesky(const value_type_array &ax, const ordinal_type verbose) {
    constexpr bool is_host = std::is_same<exec_memory_space, Kokkos::HostSpace>::value;
    Kokkos::Timer timer;
    Kokkos::Timer tick;
    double time_parallel = 0.0;
    double time_device = 0.0;
    double time_update = 0.0;

    timer.reset();
    value_type_array work;
    {
      _buf = value_type_array(do_not_initialize_tag("buf"), _bufsize_factorize);
      track_alloc(_buf.span() * sizeof(value_type));

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      value_type_matrix T(NULL, _info.max_supernode_size, _info.max_supernode_size);
      const size_type worksize = Chol<Uplo::Upper, Algo::OnDevice>::invoke(_handle_lapack, T, work);

      work = value_type_array(do_not_initialize_tag("work"), worksize * (_nstreams + 1));
      track_alloc(work.span() * sizeof(value_type));
#endif
    }
    stat.t_extra = timer.seconds();

    timer.reset();
    {
      _ax = ax; // matrix values
      constexpr bool copy_to_l_buf(false);
      _info.copySparseToSuperpanels(copy_to_l_buf, _ap, _aj, _ax, _perm, _peri);
    }
    if (_nstreams > 1) {
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
      {
        typedef TeamFunctor_FactorizeChol<supernode_info_type> functor_type;
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

        int rval = 0;
        team_policy_factor policy_factor(1, 1, 1);
        team_policy_update policy_update(1, 1, 1);
        functor_type functor(_info, _factorize_mode, _level_sids, _buf, &rval);

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
              // get max teamm size
              policy_factor = team_policy_factor(pcnt, 1, std::min(vector_size_factor[idx],vmax));
              policy_update = team_policy_update(pcnt, 1, std::min(vector_size_update[idx],vmax));
              const ordinal_type factor_tmax = policy_factor.team_size_max(functor, Kokkos::ParallelForTag());
              const ordinal_type update_tmax = policy_update.team_size_max(functor, Kokkos::ParallelForTag());;

              // create policies
              policy_factor = team_policy_factor(pcnt, std::min(team_size_factor[idx],factor_tmax), std::min(vector_size_factor[idx],vmax));
              policy_update = team_policy_update(pcnt, std::min(team_size_update[idx],update_tmax), std::min(vector_size_update[idx],vmax));
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
            factorizeCholeskyOnDevice(pbeg, pend, h_buf_factor_ptr, work);
            if (verbose) {
              Kokkos::fence(); time_device += tick.seconds();
              tick.reset();
            }
            Kokkos::fence();
            if (rval != 0) {
              TACHO_TEST_FOR_EXCEPTION(rval, std::runtime_error, "POTRF (team) returns non-zero error code.");
            }

            Kokkos::parallel_for("update factor", policy_update, functor);
            if (verbose) {
              Kokkos::fence(); time_update += tick.seconds();
            }
            ++stat_level.n_kernel_launching;
            exec_space().fence(); // Kokkos::fence();
          }
        }
      }
    } // end of Cholesky
    stat.t_factor = timer.seconds();
    timer.reset();
    if (variant == 3) {
      // compress each partitioned inverse at each level into CRS matrix
      bool lu = false;
      extractCRS(lu);
    }
    {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      track_free(work.span() * sizeof(value_type));
#endif
      track_free(_buf.span() * sizeof(value_type));
      _buf = value_type_array();
    }
    stat.t_extra += timer.seconds();

    if (verbose) {
      printf("Summary: LevelSetTools-Variant-%d (CholeskyFactorize)\n", variant);
      printf("=====================================================\n");
      printf( "\n  ** Team = %f s, Device = %f s, Update = %f s **\n\n",time_parallel,time_device,time_update );
      print_stat_factor();
      fflush(stdout);
    }
  }

  inline void solveCholeskyLowerOnDeviceVar0(const ordinal_type pbeg, const ordinal_type pend,
                                             const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const ordinal_type nrhs = t.extent(1);
    const value_type minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
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
          }
        }
      }
    }
  }

  inline void solveCholeskyLowerOnDeviceVar1(const ordinal_type pbeg, const ordinal_type pend,
                                             const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const ordinal_type nrhs = t.extent(1);
    const value_type one(1), minus_one(-1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
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

            _status = Gemv<Trans::ConjTranspose, Algo::OnDevice>::invoke(handle_blas, one, ATL, tT, zero, bT);
            checkDeviceBlasStatus("gemv");

            if (n_m > 0) {
              // solve offdiag
              UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m);
              UnmanagedViewType<value_type_matrix> bB(bptr, n_m, nrhs);

              _status = Gemv<Trans::ConjTranspose, Algo::OnDevice>::invoke(handle_blas, minus_one, ATR, bT, zero, bB);
              checkDeviceBlasStatus("gemv");
            }
          }
        }
      }
    }
  }

  inline void solveGenericLowerOnDeviceVar2_SpMV(const ordinal_type lvl, const ordinal_type nlvls,
                                                 const ordinal_type pbeg, const ordinal_type pend,
                                                 const value_type_matrix &t) {
#if (defined(KOKKOS_ENABLE_CUDA) && defined(TACHO_HAVE_CUSPARSE)) || \
     defined(KOKKOS_ENABLE_HIP)
    const ordinal_type m = t.extent(0);
    const ordinal_type nrhs = t.extent(1);
#if defined(KOKKOS_ENABLE_CUDA)
    cudaDataType computeType = CUDA_R_64F;
    if (std::is_same<value_type, float>::value) {
      computeType = CUDA_R_32F;
    } else if (!std::is_same<value_type, double>::value) {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "LevelSetTools::solveCholeskyLowerOnDevice: ComputeSPMV only supported double or float");
    }
#elif defined(KOKKOS_ENABLE_HIP)
    rocsparse_datatype rocsparse_compute_type = rocsparse_datatype_f64_r;
    if (std::is_same<value_type, float>::value) {
      rocsparse_compute_type = rocsparse_datatype_f32_r;
    } else if (!std::is_same<value_type, double>::value) {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "LevelSetTools::solveCholeskyLowerOnDevice: ComputeSPMV only supported double or float");
    }
#endif
    #ifdef TACHO_INSERT_DIAGONALS
    // compute t = L^{-1}*w
    const value_type alpha (1);
    const value_type beta  (0);
    if (_w_vec.extent(1) != nrhs) {
      // expand workspace
      Kokkos::resize(_w_vec, m, nrhs);
      // attach to Cusparse/Rocsparse data struct
      int ldw = _w_vec.stride(1);
#if defined(KOKKOS_ENABLE_CUDA)
      // destroy previous
      cusparseDestroyDnMat(matW);
      cusparseDestroyDnVec(vecW);
      // create new
      cusparseCreateDnMat(&matW, m, nrhs, ldw, (void*)(_w_vec.data()), computeType, CUSPARSE_ORDER_COL);
      cusparseCreateDnVec(&vecW, m, (void*)(_w_vec.data()), computeType);
#elif defined(KOKKOS_ENABLE_HIP)
      // destroy previous
      rocsparse_destroy_dnmat_descr(matW);
      rocsparse_destroy_dnvec_descr(vecW);
      // create new
      rocsparse_create_dnmat_descr(&matW, m, nrhs, ldw, (void*)(_w_vec.data()), rocsparse_compute_type, rocsparse_order_column);
      rocsparse_create_dnvec_descr(&vecW, m, (void*)(_w_vec.data()), rocsparse_compute_type);
#endif
    }
    const ordinal_type ldt = t.stride(1);
    auto &s0 = _h_supernodes(_h_level_sids(pbeg));
    #else
    exit(0);
    #endif
#if defined(KOKKOS_ENABLE_CUDA)
    // Desctory old CSR
    cusparseDestroySpMat(s0.L_cusparse);
    // Re-create CuSparse CSR
    if (s0.spmv_explicit_transpose) {
      cusparseCreateCsr(&s0.L_cusparse, m, m, s0.nnzL,
                        s0.rowptrL, s0.colindL, s0.nzvalsL,
                        CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                        CUSPARSE_INDEX_BASE_ZERO, computeType);
    } else {
      cusparseCreateCsr(&s0.L_cusparse, m, m, s0.nnzU,
                        s0.rowptrU, s0.colindU, s0.nzvalsU,
                        CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                        CUSPARSE_INDEX_BASE_ZERO, computeType);
    }
    // Call SpMV/SPMM
    cusparseStatus_t status;
    cusparseOperation_t opL = (s0.spmv_explicit_transpose ? CUSPARSE_OPERATION_NON_TRANSPOSE : CUSPARSE_OPERATION_TRANSPOSE);
    if (nrhs > 1) {
      if (lvl == nlvls-1) {
        // start : destroy previous
        cusparseDestroyDnMat(matL);
        // start : create DnMat for T
        cusparseCreateDnMat(&matL, m, nrhs, ldt, (void*)(t.data()), computeType, CUSPARSE_ORDER_COL);
      }
      // create vectors
      auto matX = ((nlvls-1-lvl)%2 == 0 ? matL : matW);
      auto matY = ((nlvls-1-lvl)%2 == 0 ? matW : matL);
      // SpMM
      status = cusparseSpMM(cusparseHandle, opL, CUSPARSE_OPERATION_NON_TRANSPOSE,
                            &alpha, s0.L_cusparse, 
                                    matX,
                            &beta,  matY,
                            computeType, TACHO_CUSPARSE_SPMM_ALG, (void*)buffer_L.data());
    } else {
      if (lvl == nlvls-1) {
        // start : destroy previous
        cusparseDestroyDnVec(vecL);
        // start : create DnMat for T
        cusparseCreateDnVec(&vecL, m, (void*)(t.data()), computeType);
      }
      // create vectors
      auto vecX = ((nlvls-1-lvl)%2 == 0 ? vecL : vecW);
      auto vecY = ((nlvls-1-lvl)%2 == 0 ? vecW : vecL);
      // SpMV
      status = cusparseSpMV(cusparseHandle, opL,
                            &alpha, s0.L_cusparse, 
                                    vecX,
                            &beta,  vecY,
                            computeType, TACHO_CUSPARSE_SPMV_ALG, (void*)buffer_L.data());
    }
    if (CUSPARSE_STATUS_SUCCESS != status) {
      printf( " Failed cusparseSpMV for SpMV (lower)\n" );
    }
#elif defined(KOKKOS_ENABLE_HIP)
    rocsparse_status status;
    if (nrhs > 1) {
      if (lvl == nlvls-1) {
        // start : destroy previous
        rocsparse_destroy_dnmat_descr(matL);
        // start : create DnMat for T
        rocsparse_create_dnmat_descr(&matL, m, nrhs, ldt, (void*)(t.data()), rocsparse_compute_type, rocsparse_order_column);
      }
      // create vectors
      auto vecX = ((nlvls-1-lvl)%2 == 0 ? matL : matW);
      auto vecY = ((nlvls-1-lvl)%2 == 0 ? matW : matL);
      if (s0.spmv_explicit_transpose) {
        size_t buffer_size_L = buffer_L.extent(0);
        status = rocsparse_spmm(rocsparseHandle, rocsparse_operation_none, rocsparse_operation_none,
                                &alpha, s0.descrL, vecX, &beta, vecY,
                                rocsparse_compute_type, rocsparse_spmm_alg_default,
                                rocsparse_spmm_stage_compute,
                                &buffer_size_L, (void*)buffer_L.data());
      } else {
        size_t buffer_size_L = buffer_L.extent(0);
        status = rocsparse_spmm(rocsparseHandle, rocsparse_operation_transpose, rocsparse_operation_none,
                                &alpha, s0.descrL, vecX, &beta, vecY, // dscrL stores the same ptrs as descrU, but optimized for trans
                                rocsparse_compute_type, rocsparse_spmm_alg_default,
                                rocsparse_spmm_stage_compute,
                                &buffer_size_L, (void*)buffer_L.data());
      }
    } else {
      if (lvl == nlvls-1) {
        // start : destroy previous
        rocsparse_destroy_dnvec_descr(vecL);
        // start : create DnVec for T
        rocsparse_create_dnvec_descr(&vecL, m, (void*)(t.data()), rocsparse_compute_type);
      }
      size_t buffer_size_L = buffer_L.extent(0);
      auto vecX = ((nlvls-1-lvl)%2 == 0 ? vecL : vecW);
      auto vecY = ((nlvls-1-lvl)%2 == 0 ? vecW : vecL);
      if (s0.spmv_explicit_transpose) {
        status = tacho_rocsparse_spmv
          (rocsparseHandle, rocsparse_operation_none,
           &alpha, s0.descrL, vecX, &beta, vecY,
           rocsparse_compute_type, rocsparse_spmv_alg_default,
           #if ROCM_VERSION >= 50400
           rocsparse_spmv_stage_compute,
           #endif
           &buffer_size_L, (void*)buffer_L.data());
      } else {
        status = tacho_rocsparse_spmv
          (rocsparseHandle, rocsparse_operation_transpose,
           &alpha, s0.descrL, vecX, &beta, vecY, // dscrL stores the same ptrs as descrU, but optimized for trans
           rocsparse_compute_type, rocsparse_spmv_alg_default,
           #if ROCM_VERSION >= 50400
           rocsparse_spmv_stage_compute,
           #endif
           &buffer_size_L, (void*)buffer_L.data());
      }
    }
    if (rocsparse_status_success != status) {
      printf( " Failed rocsparse_spmv for L\n" );
    }
#else
    const value_type zero(0);
    auto h_w = Kokkos::create_mirror_view_and_copy(host_memory_space(), ((nlvls-1-lvl)%2 == 0 ? t : _w_vec));
    auto h_t = Kokkos::create_mirror_view(host_memory_space(), ((nlvls-1-lvl)%2 == 0 ? _w_vec : t));
    Kokkos::deep_copy(h_t, zero);

    if (s0.spmv_explicit_transpose) {
      UnmanagedViewType<int_type_array>    d_rowptrL(s0.rowptrL, m+1);
      UnmanagedViewType<int_type_array>    d_colindL(s0.colindL, s0.nnzL);
      UnmanagedViewType<value_type_array>  d_nzvalsL(s0.nzvalsL, s0.nnzL);
      auto h_rowptr = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_rowptrL);
      auto h_colind = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_colindL);
      auto h_nzvals = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_nzvalsL);
      for (ordinal_type i = 0; i < m ; i++) {
        for (int k = h_rowptr(i); k < h_rowptr(i+1); k++) {
          for (int j = 0; j < nrhs; j++) {
            h_t(i, j) += h_nzvals(k) * h_w(h_colind(k), j);
          }
        }
      }
    } else {
      UnmanagedViewType<int_type_array>    d_rowptrU(s0.rowptrU, m+1);
      UnmanagedViewType<int_type_array>    d_colindU(s0.colindU, s0.nnzU);
      UnmanagedViewType<value_type_array>  d_nzvalsU(s0.nzvalsU, s0.nnzU);
      auto h_rowptr = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_rowptrU);
      auto h_colind = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_colindU);
      auto h_nzvals = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_nzvalsU);
      for (ordinal_type i = 0; i < m ; i++) {
        for (int k = h_rowptr(i); k < h_rowptr(i+1); k++) {
          for (int j = 0; j < nrhs; j++) {
            h_t(h_colind(k), j) += h_nzvals(k) * h_w(i, j);
          }
        }
      }
    }
#endif
    if (lvl == 0) {
      // end : copy to output
      if ((nlvls-1)%2 == 0) {
        Kokkos::deep_copy(t, _w_vec);
      }
    }
#endif
  }

  inline void solveCholeskyLowerOnDeviceVar2(const ordinal_type pbeg, const ordinal_type pend,
                                             const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const ordinal_type nrhs = t.extent(1);
    const value_type one(1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
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
            value_type *aptr = s.u_buf;
            UnmanagedViewType<value_type_matrix> AT(aptr, m, n);

            value_type *bptr = _buf.data() + h_buf_solve_ptr(p - pbeg);
            UnmanagedViewType<value_type_matrix> b(bptr, n, nrhs);

            const ordinal_type offm = s.row_begin;
            auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());

            _status = Gemv<Trans::ConjTranspose, Algo::OnDevice>::invoke(handle_blas, one, AT, tT, zero, b);
            checkDeviceBlasStatus("gemv");
          }
        }
      }
    }
  }

  inline void solveCholeskyLowerOnDevice(const ordinal_type lvl, const ordinal_type nlvls,
                                         const ordinal_type pbeg, const ordinal_type pend,
                                         const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    if (variant == 0)
      solveCholeskyLowerOnDeviceVar0(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 1)
      solveCholeskyLowerOnDeviceVar1(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 2)
      solveCholeskyLowerOnDeviceVar2(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 3) {
      solveGenericLowerOnDeviceVar2_SpMV(lvl, nlvls, pbeg, pend, t);
    } else {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "LevelSetTools::solveCholeskyLowerOnDevice, algorithm variant is not supported");
    }
  }

  inline void solveCholeskyUpperOnDeviceVar0(const ordinal_type pbeg, const ordinal_type pend,
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
                Trsv<Uplo::Upper, Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, Diag::NonUnit(), ATL, tT);
            checkDeviceBlasStatus("trsv");
          }
        }
      }
    }
  }

  inline void solveCholeskyUpperOnDeviceVar1(const ordinal_type pbeg, const ordinal_type pend,
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
            aptr += m * m;
            const UnmanagedViewType<value_type_matrix> bT(bptr, m, nrhs);
            bptr += m * nrhs;

            const ordinal_type offm = s.row_begin;
            const auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());

            if (n_m > 0) {
              const UnmanagedViewType<value_type_matrix> ATR(aptr, m, n_m); // aptr += m*n;
              const UnmanagedViewType<value_type_matrix> bB(bptr, n_m, nrhs);
              _status = Gemv<Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, minus_one, ATR, bB, one, tT);
              checkDeviceBlasStatus("gemv");
            }

            _status = Gemv<Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, one, ATL, tT, zero, bT);
            checkDeviceBlasStatus("gemv");

            _status = Copy<Algo::OnDevice>::invoke(exec_instance, tT, bT);
            checkDeviceBlasStatus("Copy");
          }
        }
      }
    }
  }

  inline void solveGenericUpperOnDeviceVar2_SpMV(const ordinal_type lvl, const ordinal_type nlvls,
                                                 const ordinal_type pbeg, const ordinal_type pend,
                                                 const value_type_matrix &t) {
#if (defined(KOKKOS_ENABLE_CUDA) && defined(TACHO_HAVE_CUSPARSE)) || \
     defined(KOKKOS_ENABLE_HIP)
    const ordinal_type m = t.extent(0);
    const ordinal_type nrhs = t.extent(1);

    auto &s0 = _h_supernodes(_h_level_sids(pbeg));

    #ifdef TACHO_INSERT_DIAGONALS
    // x = t & y = w (lvl = 0,2,4)
    // compute t = L^{-1}*w
    const value_type alpha (1);
    const value_type beta  (0);
    const ordinal_type ldt = t.stride(1);
    #else
    exit(0);
    #endif
#if defined(KOKKOS_ENABLE_CUDA)
    cudaDataType computeType = CUDA_R_64F;
    if (std::is_same<value_type, float>::value) {
      computeType = CUDA_R_32F;
    } else if (!std::is_same<value_type, double>::value) {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "LevelSetTools::solveCholeskyLowerOnDevice: ComputeSPMV only supported double or float");
    }

    cusparseStatus_t status;
    // Desctory old CSR
    cusparseDestroySpMat(s0.U_cusparse);
    // Re-create CuSparse CSR
    cusparseCreateCsr(&s0.U_cusparse, m, m, s0.nnzU,
                      s0.rowptrU, s0.colindU, s0.nzvalsU,
                      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                      CUSPARSE_INDEX_BASE_ZERO, computeType);

    // Call SpMV/SPMM
    if (nrhs > 1) {
      if (lvl == 0) {
        // start : destroy previous
        cusparseDestroyDnMat(matU);
        // start : create DnMat for T
        cusparseCreateDnMat(&matU, m, nrhs, ldt, (void*)(t.data()), computeType, CUSPARSE_ORDER_COL);
      }
      auto vecX = (lvl%2 == 0 ? matU : matW);
      auto vecY = (lvl%2 == 0 ? matW : matU);
      // SpMM
      status = cusparseSpMM(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE,
                            &alpha, s0.U_cusparse, 
                                    vecX,
                            &beta,  vecY,
                            computeType, TACHO_CUSPARSE_SPMM_ALG, (void*)buffer_U.data());
    } else {
      if (lvl == 0) {
        // start : destroy previous
        cusparseDestroyDnVec(vecU);
        // start : create DnMat for T
        cusparseCreateDnVec(&vecU, m, (void*)(t.data()), computeType);
      }
      auto vecX = (lvl%2 == 0 ? vecU : vecW);
      auto vecY = (lvl%2 == 0 ? vecW : vecU);
      // SpMV
      status = cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                            &alpha, s0.U_cusparse, 
                                    vecX,
                            &beta,  vecY,
                            computeType, TACHO_CUSPARSE_SPMV_ALG, (void*)buffer_U.data());
    }
    if (CUSPARSE_STATUS_SUCCESS != status) {
       printf( " Failed cusparseSpMV for SpMV (upper)\n" );
    }
#elif defined(KOKKOS_ENABLE_HIP)
    rocsparse_datatype rocsparse_compute_type = rocsparse_datatype_f64_r;
    if (std::is_same<value_type, float>::value) {
      rocsparse_compute_type = rocsparse_datatype_f32_r;
    } else if (!std::is_same<value_type, double>::value) {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "LevelSetTools::solveCholeskyLowerOnDevice: ComputeSPMV only supported double or float");
    }
    size_t buffer_size_U = buffer_U.extent(0);
    rocsparse_status status;
    if (nrhs > 1) {
      if (lvl == 0) {
        // start : create DnMat for T
        rocsparse_destroy_dnmat_descr(matU);
        rocsparse_create_dnmat_descr(&matU, m, nrhs, ldt, (void*)(t.data()), rocsparse_compute_type, rocsparse_order_column);
      }
      auto vecX = (lvl%2 == 0 ? matU : matW);
      auto vecY = (lvl%2 == 0 ? matW : matU);
      status = rocsparse_spmm(rocsparseHandle, rocsparse_operation_none, rocsparse_operation_none,
                              &alpha, s0.descrU, vecX, &beta, vecY,
                              rocsparse_compute_type, rocsparse_spmm_alg_default,
                              rocsparse_spmm_stage_compute,
                              &buffer_size_U, (void*)buffer_U.data());
    } else {
      if (lvl == 0) {
        // start : create DnVec for T
        rocsparse_destroy_dnvec_descr(vecU);
        rocsparse_create_dnvec_descr(&vecU, m, (void*)(t.data()), rocsparse_compute_type);
      }
      auto vecX = (lvl%2 == 0 ? vecU : vecW);
      auto vecY = (lvl%2 == 0 ? vecW : vecU);
      status = tacho_rocsparse_spmv
          (rocsparseHandle, rocsparse_operation_none,
           &alpha, s0.descrU, vecX, &beta, vecY,
           rocsparse_compute_type, rocsparse_spmv_alg_default,
           #if ROCM_VERSION >= 50400
           rocsparse_spmv_stage_compute,
           #endif
           &buffer_size_U, (void*)buffer_U.data());
    }
    if (rocsparse_status_success != status) {
      printf( " Failed rocsparse_spmv for U\n" );
    }
#else
    const value_type zero(0);
    auto h_w = Kokkos::create_mirror_view_and_copy(host_memory_space(), (lvl%2 == 0 ? t : _w_vec));
    auto h_t = Kokkos::create_mirror_view(host_memory_space(), (lvl%2 == 0 ? _w_vec : t));
    Kokkos::deep_copy(h_t, zero);

    UnmanagedViewType<int_type_array>    d_rowptrU(s0.rowptrU, m+1);
    UnmanagedViewType<int_type_array>    d_colindU(s0.colindU, s0.nnzU);
    UnmanagedViewType<value_type_array>  d_nzvalsU(s0.nzvalsU, s0.nnzU);
    auto h_rowptr = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_rowptrU);
    auto h_colind = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_colindU);
    auto h_nzvals = Kokkos::create_mirror_view_and_copy(host_memory_space(), d_nzvalsU);

    for (ordinal_type i = 0; i < m ; i++) {
      for (int k = h_rowptr(i); k < h_rowptr(i+1); k++) {
        for (int j = 0; j < nrhs; j++) {
          h_t(i, j) += h_nzvals(k) * h_w(h_colind(k), j);
        }
      }
    }
    Kokkos::deep_copy(t, h_t);
#endif
    if (lvl == nlvls-1) {
      // end : copy to output
      if (lvl%2 == 0) {
        Kokkos::deep_copy(t, _w_vec);
      }
    }
#endif
  }

  inline void solveCholeskyUpperOnDeviceVar2(const ordinal_type pbeg, const ordinal_type pend,
                                             const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    const ordinal_type nrhs = t.extent(1);
    const value_type one(1), zero(0);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    ordinal_type q(0);
#endif
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
          if (m > 0 && n > 0) {
            value_type *aptr = s.u_buf, *bptr = _buf.data() + h_buf_solve_ptr(p - pbeg);
            ;
            const UnmanagedViewType<value_type_matrix> AT(aptr, m, n);
            const UnmanagedViewType<value_type_matrix> b(bptr, n, nrhs);

            const ordinal_type offm = s.row_begin;
            const auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());

            _status = Gemv<Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, one, AT, b, zero, tT);
            checkDeviceBlasStatus("gemv");
          }
        }
      }
    }
  }

  inline void solveCholeskyUpperOnDevice(const ordinal_type lvl, const ordinal_type nlvls,
                                         const ordinal_type pbeg, const ordinal_type pend,
                                         const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    if (variant == 0)
      solveCholeskyUpperOnDeviceVar0(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 1)
      solveCholeskyUpperOnDeviceVar1(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 2)
      solveCholeskyUpperOnDeviceVar2(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 3) {
      solveGenericUpperOnDeviceVar2_SpMV(lvl, nlvls, pbeg, pend, t);
    } else {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "LevelSetTools::solveCholeskyUpperOnDevice, algorithm variant is not supported");
    }
  }

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

  inline void solveLDL_LowerOnDevice(const ordinal_type pbeg, const ordinal_type pend,
                                     const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    if (variant == 0)
      solveLDL_LowerOnDeviceVar0(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 1)
      solveLDL_LowerOnDeviceVar1(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 2)
      solveLDL_LowerOnDeviceVar2(pbeg, pend, h_buf_solve_ptr, t);
    else {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "LevelSetTools::solveLDL_LowerOnDevice, algorithm variant is not supported");
    }
  }

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
            const UnmanagedViewType<value_type_matrix> bT(bptr, m, nrhs);

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

  inline void solveLDL_UpperOnDevice(const ordinal_type pbeg, const ordinal_type pend,
                                     const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    if (variant == 0)
      solveLDL_UpperOnDeviceVar0(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 1)
      solveLDL_UpperOnDeviceVar1(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 2)
      solveLDL_UpperOnDeviceVar2(pbeg, pend, h_buf_solve_ptr, t);
    else {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "LevelSetTools::solveLDL_UpperOnDevice, algorithm variant is not supported");
    }
  }

  inline void solveLU_LowerOnDeviceVar0(const ordinal_type pbeg, const ordinal_type pend,
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
          }
        }
      }
    }
  }

  inline void solveLU_LowerOnDeviceVar1(const ordinal_type pbeg, const ordinal_type pend,
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
          }
        }
      }
    }
  }

  inline void solveLU_LowerOnDeviceVar2(const ordinal_type pbeg, const ordinal_type pend,
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
          }
        }
      }
    }
  }

  inline void solveLU_LowerOnDevice(const ordinal_type lvl, const ordinal_type nlvls,
                                    const ordinal_type pbeg, const ordinal_type pend,
                                    const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    if (variant == 0)
      solveLU_LowerOnDeviceVar0(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 1)
      solveLU_LowerOnDeviceVar1(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 2)
      solveLU_LowerOnDeviceVar2(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 3) {
      solveGenericLowerOnDeviceVar2_SpMV(lvl, nlvls, pbeg, pend, t);
    } else {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "LevelSetTools::solveLU_LowerOnDevice, algorithm variant is not supported");
    }
  }

  inline void solveLU_UpperOnDeviceVar0(const ordinal_type pbeg, const ordinal_type pend,
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
        }
      }
    }
  }
  inline void solveLU_UpperOnDeviceVar1(const ordinal_type pbeg, const ordinal_type pend,
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
          }
        }
      }
    }
  }

  inline void solveLU_UpperOnDeviceVar2(const ordinal_type pbeg, const ordinal_type pend,
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
            value_type *bptr = _buf.data() + h_buf_solve_ptr(p - pbeg);

            const UnmanagedViewType<value_type_matrix> AT(s.u_buf, m, n);
            const UnmanagedViewType<value_type_matrix> b(bptr, n, nrhs);

            const ordinal_type offm = s.row_begin;
            const auto tT = Kokkos::subview(t, range_type(offm, offm + m), Kokkos::ALL());
            _status = Gemv<Trans::NoTranspose, Algo::OnDevice>::invoke(handle_blas, one, AT, b, zero, tT);
            checkDeviceBlasStatus("gemv");
          }
        }
      }
    }
  }

  inline void solveLU_UpperOnDevice(const ordinal_type lvl, const ordinal_type nlvls,
                                    const ordinal_type pbeg, const ordinal_type pend,
                                    const size_type_array_host &h_buf_solve_ptr, const value_type_matrix &t) {
    if (variant == 0)
      solveLU_UpperOnDeviceVar0(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 1)
      solveLU_UpperOnDeviceVar1(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 2)
      solveLU_UpperOnDeviceVar2(pbeg, pend, h_buf_solve_ptr, t);
    else if (variant == 3) {
      solveGenericUpperOnDeviceVar2_SpMV(lvl, nlvls, pbeg, pend, t);
    } else {
      TACHO_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "LevelSetTools::solveLU_UpperOnDevice, algorithm variant is not supported");
    }
  }

  inline void allocateWorkspaceSolve(const ordinal_type nrhs) {
    if (variant == 3) {
    } else {
      const size_type buf_extent = _bufsize_solve * nrhs;
      const size_type buf_span = _buf.span();

      if (buf_extent != buf_span) {
        _buf = value_type_array(do_not_initialize_tag("buf"), buf_extent);
        track_free(buf_span * sizeof(value_type));
        track_alloc(_buf.span() * sizeof(value_type));
        {
          const Kokkos::RangePolicy<exec_space> policy(0, _buf_solve_ptr.extent(0));
          const auto buf_solve_nrhs_ptr = _buf_solve_nrhs_ptr;
          const auto buf_solve_ptr = _buf_solve_ptr;
          Kokkos::parallel_for(
              policy, KOKKOS_LAMBDA(const ordinal_type &i) { buf_solve_nrhs_ptr(i) = nrhs * buf_solve_ptr(i); });
        }
        Kokkos::deep_copy(_h_buf_solve_nrhs_ptr, _buf_solve_nrhs_ptr);
      }
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

    // 0. permute and copy b -> t
    const auto exec_instance = exec_space();
    ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(exec_instance, b, _perm, t);
    stat.t_extra = timer.seconds();

    timer.reset();
    {
#if defined(TACHO_ENABLE_SOLVE_CHOLESKY_USE_LIGHT_KERNEL)
      const auto work_item_property = Kokkos::Experimental::WorkItemProperty::HintLightWeight;
#endif
      // this should be considered with average problem sizes in levels
      const ordinal_type half_level = _nlevel / 2;
      const ordinal_type team_size_solve[2] = {64, 16}, vector_size_solve[2] = {8, 8};
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
              policy_solve = team_policy_solve(pcnt, team_size_solve[idx], vector_size_solve[idx]);
              policy_update = team_policy_update(pcnt, team_size_update[idx], vector_size_update[idx]);
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
            if (lvl < _device_level_cut) {
              // do nothing
              // Kokkos::parallel_for("solve lower", policy_solve, functor);
            } else {
              Kokkos::parallel_for("solve lower", policy_solve_with_work_property, functor);
              ++stat_level.n_kernel_launching;
            }
            const auto h_buf_solve_ptr = Kokkos::subview(_h_buf_solve_nrhs_ptr, range_solve_buf_ptr);
            solveCholeskyLowerOnDevice(lvl, _team_serial_level_cut, pbeg, pend, h_buf_solve_ptr, t);

            if (variant != 3) {
              Kokkos::fence();
              Kokkos::parallel_for("update lower", policy_update_with_work_property, functor);
              ++stat_level.n_kernel_launching;
              exec_space().fence();
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
              policy_solve = team_policy_solve(pcnt, team_size_solve[idx], vector_size_solve[idx]);
              policy_update = team_policy_update(pcnt, team_size_update[idx], vector_size_update[idx]);
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
              exec_space().fence();

              if (lvl < _device_level_cut) {
                // do nothing
                // Kokkos::parallel_for("solve upper", policy_solve, functor);
              } else {
                Kokkos::parallel_for("solve upper", policy_solve_with_work_property, functor);
                ++stat_level.n_kernel_launching;
              }
            }
            const auto h_buf_solve_ptr = Kokkos::subview(_h_buf_solve_nrhs_ptr, range_solve_buf_ptr);
            solveCholeskyUpperOnDevice(lvl, _team_serial_level_cut, pbeg, pend, h_buf_solve_ptr, t);
            if (variant != 3) {
              Kokkos::fence();
            }
          }
        }
      } /// end of upper tri solve

    } // end of solve
    stat.t_solve = timer.seconds();

    // permute and copy t -> x
    timer.reset();
    ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(exec_instance, t, _peri, x);
    stat.t_extra += timer.seconds();

    if (verbose) {
      printf("Summary: LevelSetTools-Variant-%d (Cholesky Solve: %3d)\n", variant, nrhs);
      printf("=======================================================\n");
      print_stat_solve();
      fflush(stdout);
    }
  }

  inline void factorizeLDL(const value_type_array &ax, const ordinal_type verbose) {
    constexpr bool is_host = std::is_same<exec_memory_space, Kokkos::HostSpace>::value;
    Kokkos::Timer timer;
    Kokkos::Timer tick;
    double time_parallel = 0.0;
    double time_device = 0.0;
    double time_update = 0.0;

    timer.reset();
    value_type_array work;
    {
      _buf = value_type_array(do_not_initialize_tag("buf"), _bufsize_factorize);
      track_alloc(_buf.span() * sizeof(value_type));

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      value_type_matrix T(NULL, _info.max_supernode_size, _info.max_supernode_size);
      ordinal_type_array P(NULL, _info.max_supernode_size);
      const size_type worksize = LDL<Uplo::Lower, Algo::OnDevice>::invoke(_handle_lapack, T, P, work);

      work = value_type_array(do_not_initialize_tag("work"), worksize * (_nstreams + 1) * max(8, _nstreams));
#else
      const size_type worksize = 32 * _info.max_supernode_size;
      work = value_type_array(do_not_initialize_tag("work"), worksize);
#endif
      track_alloc(work.span() * sizeof(value_type));
    }
    stat.t_extra = timer.seconds();

    timer.reset();
    {
      _ax = ax; // matrix values
      constexpr bool copy_to_l_buf(false);
      _info.copySparseToSuperpanels(copy_to_l_buf, _ap, _aj, _ax, _perm, _peri);
    }
    if (_nstreams > 1) {
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
      {
        typedef TeamFunctor_FactorizeLDL<supernode_info_type> functor_type;
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
        int rval = 0;
        team_policy_factor policy_factor(1, 1, 1);
        team_policy_update policy_update(1, 1, 1);
        functor_type functor(_info, _factorize_mode, _level_sids, _piv, _diag, _buf, &rval);

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
              // get max teamm sizes
              policy_factor = team_policy_factor(pcnt, 1, std::min(vector_size_factor[idx],vmax));
              policy_update = team_policy_update(pcnt, 1, std::min(vector_size_update[idx],vmax));
              const ordinal_type factor_tmax = policy_factor.team_size_max(functor, Kokkos::ParallelForTag());
              const ordinal_type update_tmax = policy_update.team_size_max(functor, Kokkos::ParallelForTag());

              policy_factor = team_policy_factor(pcnt, std::min(team_size_factor[idx],factor_tmax), std::min(vector_size_factor[idx],vmax));
              policy_update = team_policy_update(pcnt, std::min(team_size_update[idx],update_tmax), std::min(vector_size_update[idx],vmax));
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
            factorizeLDL_OnDevice(pbeg, pend, h_buf_factor_ptr, work);
            if (verbose) {
              Kokkos::fence(); time_device += tick.seconds();
              tick.reset();
            }
            Kokkos::fence();
            if (rval != 0) {
              TACHO_TEST_FOR_EXCEPTION(rval, std::runtime_error, "SYTRF (team) returns non-zero error code.");
            }

            Kokkos::parallel_for("update factor", policy_update, functor);
            if (verbose) {
              Kokkos::fence(); time_update += tick.seconds();
            }
            ++stat_level.n_kernel_launching;
            exec_space().fence();
          }
          const auto exec_instance = exec_space();
          Kokkos::deep_copy(exec_instance, _h_supernodes, _info.supernodes);
        }
      }
    } // end of LDL
    stat.t_factor = timer.seconds();
    timer.reset();
    {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      track_free(work.span() * sizeof(value_type));
#endif
      track_free(_buf.span() * sizeof(value_type));
      _buf = value_type_array();
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

    // 0. permute and copy b -> t
    const auto exec_instance = exec_space();
    ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(exec_instance, b, _perm, t);
    stat.t_extra = timer.seconds();

    timer.reset();
    {
#if defined(TACHO_ENABLE_SOLVE_CHOLESKY_USE_LIGHT_KERNEL)
      const auto work_item_property = Kokkos::Experimental::WorkItemProperty::HintLightWeight;
#endif
      // this should be considered with average problem sizes in levels
      const ordinal_type half_level = _nlevel / 2;
#if defined(CUDA_VERSION)
#if (11000 > CUDA_VERSION)
      /// cuda 11.1 below
      const ordinal_type team_size_solve[2] = {32, 16}, vector_size_solve[2] = {8, 8};
#else
      /// cuda 11.1 and higher
      const ordinal_type team_size_solve[2] = {32, 16}, vector_size_solve[2] = {8, 8};
#endif
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
              policy_solve = team_policy_solve(pcnt, team_size_solve[idx], vector_size_solve[idx]);
              policy_update = team_policy_update(pcnt, team_size_update[idx], vector_size_update[idx]);
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
            if (lvl < _device_level_cut) {
              // do nothing
              // Kokkos::parallel_for("solve lower", policy_solve, functor);
            } else {
              Kokkos::parallel_for("solve lower", policy_solve_with_work_property, functor);
              ++stat_level.n_kernel_launching;
            }
            const auto h_buf_solve_ptr = Kokkos::subview(_h_buf_solve_nrhs_ptr, range_solve_buf_ptr);
            solveLDL_LowerOnDevice(pbeg, pend, h_buf_solve_ptr, t);
            Kokkos::fence();

            Kokkos::parallel_for("update lower", policy_update_with_work_property, functor);
            ++stat_level.n_kernel_launching;
            exec_space().fence();
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
              policy_solve = team_policy_solve(pcnt, team_size_solve[idx], vector_size_solve[idx]);
              policy_update = team_policy_update(pcnt, team_size_update[idx], vector_size_update[idx]);
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
            Kokkos::parallel_for("update upper", policy_update_with_work_property, functor);
            ++stat_level.n_kernel_launching;
            exec_space().fence();

            if (lvl < _device_level_cut) {
              // do nothing
              // Kokkos::parallel_for("solve upper", policy_solve, functor);
            } else {
              Kokkos::parallel_for("solve upper", policy_solve_with_work_property, functor);
              ++stat_level.n_kernel_launching;
            }

            const auto h_buf_solve_ptr = Kokkos::subview(_h_buf_solve_nrhs_ptr, range_solve_buf_ptr);
            solveLDL_UpperOnDevice(pbeg, pend, h_buf_solve_ptr, t);
            Kokkos::fence();
          }
        }
      } /// end of upper tri solve

    } // end of solve
    stat.t_solve = timer.seconds();

    // permute and copy t -> x
    timer.reset();
    ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(exec_instance, t, _peri, x);
    stat.t_extra += timer.seconds();

    if (verbose) {
      printf("Summary: LevelSetTools-Variant-%d (LDL Solve: %3d)\n", variant, nrhs);
      printf("==================================================\n");
      print_stat_solve();
      fflush(stdout);
    }
  }

  inline void factorizeLU(const value_type_array &ax, const ordinal_type verbose) {
    constexpr bool is_host = std::is_same<exec_memory_space, Kokkos::HostSpace>::value;
    Kokkos::Timer timer;
    Kokkos::Timer tick;
    double time_parallel = 0.0;
    double time_device = 0.0;
    double time_update = 0.0;

    timer.reset();
    value_type_array work;
    {
      _buf = value_type_array(do_not_initialize_tag("buf"), _bufsize_factorize);
      track_alloc(_buf.span() * sizeof(value_type));

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      // NOTE : move this to symbolic with the actual max worksize?
      value_type_matrix T(NULL, _info.max_supernode_size, _info.max_num_cols);
      ordinal_type_array P(NULL, std::min(_info.max_supernode_size, _info.max_num_cols));
      const size_type worksize = LU<Algo::OnDevice>::invoke(_handle_lapack, T, P, work);

      work = value_type_array(do_not_initialize_tag("work"), worksize * (_nstreams + 1));
      // work = value_type_array(do_not_initialize_tag("work"), worksize*_nstreams);
#endif
      track_alloc(work.span() * sizeof(value_type));
    }
    stat.t_extra = timer.seconds();

    timer.reset();
    {
      _ax = ax; // matrix values
      constexpr bool copy_to_l_buf(true);
      _info.copySparseToSuperpanels(copy_to_l_buf, _ap, _aj, _ax, _perm, _peri);
    }
    if (_nstreams > 1) {
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
        int rval = 0;
        team_policy_factor policy_factor(1, 1, 1);
        team_policy_update policy_update(1, 1, 1);
        functor_type functor(_info, _factorize_mode, _level_sids, _piv, _buf, &rval);

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
              // get max teamm sizes
              policy_factor = team_policy_factor(pcnt, 1, std::min(vector_size_factor[idx],vmax));
              policy_update = team_policy_update(pcnt, 1, std::min(vector_size_update[idx],vmax));
              const ordinal_type factor_tmax = policy_factor.team_size_max(functor, Kokkos::ParallelForTag());
              const ordinal_type update_tmax = policy_update.team_size_max(functor, Kokkos::ParallelForTag());

              // create policies
              policy_factor = team_policy_factor(pcnt, std::min(team_size_factor[idx],factor_tmax), std::min(vector_size_factor[idx],vmax));
              policy_update = team_policy_update(pcnt, std::min(team_size_update[idx],update_tmax), std::min(vector_size_update[idx],vmax));
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
            factorizeLU_OnDevice(pbeg, pend, h_buf_factor_ptr, work);
            if (verbose) {
              Kokkos::fence(); time_device += tick.seconds();
              tick.reset();
            }
            Kokkos::fence();
            if (rval != 0) {
              TACHO_TEST_FOR_EXCEPTION(rval, std::runtime_error, "GETRF (team) returns non-zero error code.");
            }

            Kokkos::parallel_for("update factor", policy_update, functor);
            if (verbose) {
              Kokkos::fence(); time_update += tick.seconds();
            }
            ++stat_level.n_kernel_launching;
            exec_space().fence();
          }
          const auto exec_instance = exec_space();
          Kokkos::deep_copy(exec_instance, _h_supernodes, _info.supernodes);
        }
      }
    } // end of LU
    stat.t_factor = timer.seconds();
    timer.reset();
    if (variant == 3) {
      // compress each partitioned inverse at each level into CRS matrix
      bool lu = true;
      extractCRS(lu);
    }
    {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      track_free(work.span() * sizeof(value_type));
#endif
      track_free(_buf.span() * sizeof(value_type));
      _buf = value_type_array();
    }
    stat.t_extra += timer.seconds();

    if (verbose) {
      printf("Summary: LevelSetTools-Variant-%d (LU Factorize)\n", variant);
      printf("================================================\n");
      printf( "\n  ** Team = %f s, Device = %f s, Update = %f s **\n\n",time_parallel,time_device,time_update );
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

    // 0. permute and copy b -> t
    const auto exec_instance = exec_space();
    ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(exec_instance, b, _perm, t);
    exec_instance.fence();
    stat.t_extra = timer.seconds();

    timer.reset();
    {
#if defined(TACHO_ENABLE_SOLVE_CHOLESKY_USE_LIGHT_KERNEL)
      const auto work_item_property = Kokkos::Experimental::WorkItemProperty::HintLightWeight;
#endif
      // this should be considered with average problem sizes in levels
      const ordinal_type half_level = _nlevel / 2;
#if defined(CUDA_VERSION)
#if (11000 > CUDA_VERSION)
      /// cuda 11.1 below
      const ordinal_type team_size_solve[2] = {32, 16}, vector_size_solve[2] = {8, 8};
#else
      /// cuda 11.1 and higher
      const ordinal_type team_size_solve[2] = {32, 16}, vector_size_solve[2] = {8, 8};
#endif
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
              policy_solve = team_policy_solve(pcnt, team_size_solve[idx], vector_size_solve[idx]);
              policy_update = team_policy_update(pcnt, team_size_update[idx], vector_size_update[idx]);
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
              Kokkos::fence();

              Kokkos::parallel_for("update lower", policy_update_with_work_property, functor);
              ++stat_level.n_kernel_launching;
              exec_space().fence();
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
              policy_solve = team_policy_solve(pcnt, team_size_solve[idx], vector_size_solve[idx]);
              policy_update = team_policy_update(pcnt, team_size_update[idx], vector_size_update[idx]);
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
              exec_space().fence();

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
              Kokkos::fence();
            }
          }
        }
      } /// end of upper tri solve

    } // end of solve
    stat.t_solve = timer.seconds();

    // permute and copy t -> x
    timer.reset();
    ApplyPermutation<Side::Left, Trans::NoTranspose, Algo::OnDevice>::invoke(exec_instance, t, _peri, x);
    exec_instance.fence();
    stat.t_extra += timer.seconds();

    if (verbose) {
      printf("Summary: LevelSetTools-Variant-%d (LU Solve: %3d)\n", variant, nrhs);
      printf("=================================================\n");
      print_stat_solve();
    }
  }

  inline void factorize(const value_type_array &ax, const ordinal_type verbose = 0) override {
    Kokkos::deep_copy(_superpanel_buf, value_type(0));
    switch (this->getSolutionMethod()) {
    case 1: { /// Cholesky
      factorizeCholesky(ax, verbose);
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
      factorizeLDL(ax, verbose);
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
      factorizeLU(ax, verbose);
      break;
    }
    default: {
      TACHO_TEST_FOR_EXCEPTION(false, std::logic_error, "The solution method is not supported");
      break;
    }
    }
  }

  inline void solve(const value_type_matrix &x, // solution
                    const value_type_matrix &b, // right hand side
                    const value_type_matrix &t, // temporary workspace (store permuted vectors)
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
    case 3: { /// LU
      solveLU(x, b, t, verbose);
      break;
    }
    default: {
      TACHO_TEST_FOR_EXCEPTION(false, std::logic_error, "The solution method is not supported");
      break;
    }
    }
  }
};

} // namespace Tacho
#endif
