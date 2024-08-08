//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
#ifndef KOKKOSBLAS3_GEMM_PERF_TEST_H_
#define KOKKOSBLAS3_GEMM_PERF_TEST_H_

// #include <complex.h>

#include "Kokkos_MathematicalFunctions.hpp"

#include "KokkosBlas3_common.hpp"
#include <Kokkos_Random.hpp>

#include <KokkosBlas3_gemm.hpp>

#include "KokkosBatched_HostLevel_Gemm.hpp"
#include "KokkosBatched_Gemm_Decl.hpp"
#include "KokkosBatched_Util.hpp"
#include "gtest/gtest.h"  // EXPECT_NEAR
#include "KokkosKernels_TestUtils.hpp"

#include <chrono>

#if defined(KOKKOSKERNELS_ENABLE_TPL_ARMPL)
#include "armpl.h"
#else
using armpl_int_t = int;
#endif  // KOKKOSKERNELS_ENABLE_TPL_ARMPL

// #define GEMM_PERF_TEST_DEBUG

// Forward declarations
void do_gemm_serial_blas(options_t options);
void do_gemm_serial_batched(options_t options);
void do_gemm_serial_batched_blocked(options_t options);
// void do_gemm_experiment(options_t options);

// void do_gemm_serial_blas_parallel(options_t options);
// Not valid! The KokkosBlas::gemm function may take the entire device per
// invocation!
void do_gemm_heuristic_batched_parallel(options_t options);
void do_gemm_serial_batched_parallel(options_t options);
void do_gemm_serial_batched_blocked_parallel(options_t options);
void do_gemm_serial_simd_batched_parallel(options_t options);
void do_gemm_serial_simd_batched_blocked_parallel(options_t options);
void do_gemm_serial_batched_compact_mkl_parallel(options_t options);
void do_gemm_team_batched_parallel(options_t options);
void do_gemm_team_batched_blocked_parallel(options_t options);
void do_gemm_team_vector_batched_parallel(options_t options);
void do_gemm_team_vector_batched_blocked_parallel(options_t options);
void do_gemm_team_simd_batched_parallel(options_t options);
void do_gemm_team_simd_batched_blocked_parallel(options_t options);
void do_gemm_experiment_parallel(options_t options);

struct SerialTag {};
struct SerialBatchDim3Tag {};
struct SerialSimdTag {};
struct SerialSimdBatchDim3Tag {};
struct TeamTag {};
struct TeamBatchDim3Tag {};
struct TeamVectorTag {};
struct TeamVectorBatchDim3Tag {};
struct TeamSimdTag {};
struct TeamSimdBatchDim4Tag {};
struct LayoutLeftTag {};
struct LayoutRightTag {};
struct SimdCpuTag {};

// gemm invoke table
void (*do_gemm_invoke[LOOP_N][TEST_N])(options_t) = {
    {
        do_gemm_serial_blas,                                     // BLAS
        do_gemm_serial_batched, do_gemm_serial_batched_blocked,  // Serial
        NULL, NULL,                                              // Team
        NULL, NULL,                                              // TeamVector
        NULL, NULL,                                              // TeamSimd
        NULL,                                                    // Serial Experiment
    },
    {
        NULL,  // BLAS
        do_gemm_heuristic_batched_parallel,
        do_gemm_serial_batched_parallel,  // Serial
        do_gemm_serial_batched_blocked_parallel, do_gemm_serial_simd_batched_parallel,
        do_gemm_serial_simd_batched_blocked_parallel, do_gemm_serial_batched_compact_mkl_parallel,
        do_gemm_team_batched_parallel,
        do_gemm_team_batched_blocked_parallel,       // Team
        do_gemm_team_vector_batched_parallel, NULL,  // TeamVector
        do_gemm_team_simd_batched_parallel,
        do_gemm_team_simd_batched_blocked_parallel,  // TeamSimd
        do_gemm_experiment_parallel                  // Parallel Experiment
    }};

/*************************** Test types and defaults **************************/
#define DEFAULT_GEMM_ARGS "NN"
#define DEFAULT_GEMM_ALPHA 1.0
#define DEFAULT_GEMM_BETA 1.0

using view_type_3d = Kokkos::View<default_scalar ***, default_layout, default_device>;
using view_type_4d = Kokkos::View<default_scalar ****, default_layout, default_device>;
using view_type_5d = Kokkos::View<default_scalar *****, default_layout, default_device>;

// Construct the vector type
using memory_space             = typename default_device::execution_space::memory_space;
constexpr int simd_vector_size = KokkosBatched::DefaultVectorLength<default_scalar, memory_space>::value;
constexpr int simd_internal_vector_size =
    KokkosBatched::DefaultInternalVectorLength<default_scalar, memory_space>::value;
using vector_type          = KokkosBatched::Vector<KokkosBatched::SIMD<default_scalar>, simd_vector_size>;
using internal_vector_type = KokkosBatched::Vector<KokkosBatched::SIMD<default_scalar>, simd_internal_vector_size>;
using vector_view_type_3d  = Kokkos::View<vector_type ***, default_layout, default_device>;
using internal_vector_view_type_4d = Kokkos::View<internal_vector_type ****, default_layout, default_device>;

struct batched_params {
  int team_size;
  int vector_len;
};
typedef struct batched_params batched_params_t;

/**
 * @brief struct gemm_simd_args encapsulates the data types required
 * for allocating and passing a single matrix to the KokkosBatched gemm
 * kernels. To invoke gemm on a batch of matrices, three instances of this
 * struct are required, one for each matrix, A, B, and C.
 *
 * @var  vec_3d: 3-rank view type used for allocating the underlying data.
 *               A reference must be kept to this object to ensure the
 *               data is not free'd by the C++ runtime.
 * @var  mat_4d: 4-rank view type used for populating the simd view with
                 random values.
 * @var ivec_4d: 4-rank view type used for passing to math kernels. This
 *               view type is used for leveraging simd instructions on
 *               both the host and device.
 */
struct gemm_simd_args {
  vector_view_type_3d vec_3d;
  view_type_4d mat_4d;
  internal_vector_view_type_4d ivec_4d;
};
typedef struct gemm_simd_args gemm_simd_args_t;

/**
 * @brief struct gemm_armpl_args encapsulates the data required for allocating
 * and passing a single matrix to the armpl interleave gemm api. To invoke gemm,
 * three instances of this struct are needed, one for each matrix: A, B, and C.
 *
 * @var mat: A flat copy of the 3-rank matrix in interleaved format.
 * @var jstrd: See
 * https://developer.arm.com/documentation/101004/2100/Interleave-batch-functions/armpl-dgemm-interleave-batch.
 * @var istrd: See
 * https://developer.arm.com/documentation/101004/2100/Interleave-batch-functions/armpl-dgemm-interleave-batch.
 * @var bstrd: See
 * https://developer.arm.com/documentation/101004/2100/Interleave-batch-functions/armpl-dgemm-interleave-batch.
 */
struct gemm_armpl_args {
  default_scalar *mat;
  armpl_int_t jstrd, istrd, bstrd;
};
typedef struct gemm_armpl_args gemm_armpl_args_t;

/**
 * @brief struct gemm_args are common arguments passed to
 * both gemm implementations in the KokkosBlas and KokkosBatched
 * namespaces throughout these performance tests.
 *
 * @var transA: transpose type for A matrix.
 *              supported types:   'n' - no transpose, 't' - transpose.
 *              unsupported types: 'c' - conjugate transpose.
 * @var transB: transpose type for B matrix.
 *              supported types:   'n' - no transpose, 't' - transpose.
 *              unsupported types: 'c' - conjugate transpose.
 * @var alpha: scalar applied to A matrix.
 * @var beta:  scalar applied to B matrix.
 * @var A:     3-rank view type used in all non-simd tests.
 * @var B:     3-rank view type used in all non-simd tests.
 * @var C:     3-rank view type used in all non-simd tests.
 * @var bp:    team_size and vector_length for tests that use
 * Kokkos::TeamPolicy.
 * @var Av:    3-rank and 4-rank vector view types for simd tests.
 * @var Bv:    3-rank and 4-rank vector view types for simd tests.
 * @var Cv:    3-rank and 4-rank vector view types for simd tests.
 * @var A_pl:  flat array with strides for armpl interleave API.
 * @var B_pl:  flat array with strides for armpl interleave API.
 * @var C_pl:  flat array with strides for armpl interleave API.
 * @var ninter: number of interleaved matrices (sub-batch size) for armpl
 * interleave API.
 * @var nbatch: number of batches of interleaves matrices for armpl interleave
 * API.
 */
struct gemm_args {
  char transA, transB;
  default_scalar alpha;
  default_scalar beta;
  view_type_3d A, B, C;
  batched_params_t bp;
  // Below are matrices for simd tests
  gemm_simd_args_t Av, Bv, Cv;
  // Below are matrices and args for armpl interleave batch tests
  gemm_armpl_args_t A_pl, B_pl, C_pl;
  armpl_int_t ninter;
  armpl_int_t nbatch;
  matrix_dims_t dims;
};
typedef struct gemm_args gemm_args_t;

static std::string gemm_csv_header_str =
    "algorithm,vector_type,transAtransB,alpha,beta,team_size,vector_len,loop_"
    "type,A_dims,B_"
    "dims,C_dims,warm_up_n,"
    "iter,total_time(s),average_time(s),FLOPS,GFLOP/average_time(s)";

/*************************** Internal helper fns **************************/
// Flop count formula from lapack working note 41:
// http://www.icl.utk.edu/~mgates3/docs/lawn41.pdf
static inline double __gemm_flop_count(double a_m, double a_n, double b_n) {
  // TODO: if not Kokkos::complex.
  if (std::is_same<double, default_scalar>::value || std::is_same<float, default_scalar>::value ||
      std::is_same<Kokkos::Experimental::half_t, default_scalar>::value ||
      std::is_same<Kokkos::Experimental::bhalf_t, default_scalar>::value)
    return 2 * a_m * b_n * a_n;
  else
    // For complex, we need to count 2 flops for each add and 6 flops for each
    // multiply.
    return (2 + 6) * a_m * b_n * a_n;
}

static inline std::string __gemm_output_dim_string(options_t options, matrix_dim_t dim) {
  std::string x   = "x";
  std::string ret = std::to_string(dim.m) + x + std::to_string(dim.n);

  if (options.blas_args.batch_size_last_dim)
    return ret + x + std::to_string(dim.k);
  else
    return std::to_string(dim.k) + x + ret;
}

static void __gemm_output_csv_row(options_t options, gemm_args_t gemm_args, double time_in_seconds,
                                  const char *experiment_name = nullptr, const char *team_size = nullptr,
                                  const char *vec_len = nullptr, const char *vec_type = nullptr) {
  std::string algo_name = !experiment_name ? test_e_str[options.test] : std::string(experiment_name);
  std::string ts        = !team_size ? std::to_string(gemm_args.bp.team_size) : std::string(team_size);
  std::string vlen      = !vec_len ? std::to_string(gemm_args.bp.vector_len) : std::string(vec_len);
  std::string vtype     = !vec_type ? internal_vector_type::label() : std::string(vec_type);
  if (options.blas_args.use_auto) ts = vlen = "Kokkos::AUTO";

  double flops;
  double gflops;
  double average_time = time_in_seconds / options.n;

  if (options.verify) return;

  flops = gemm_args.dims.a.k * __gemm_flop_count(gemm_args.dims.a.m, gemm_args.dims.a.n, gemm_args.dims.b.n);

  gflops = flops / 1e9;

  options.out[0] << algo_name << "," << vtype << "," << options.blas_args.gemm.gemm_args << ","
                 << static_cast<double>(options.blas_args.gemm.alpha) << ","
                 << static_cast<double>(options.blas_args.gemm.beta) << "," << ts << "," << vlen << ","
                 << loop_e_str[options.loop] << "," << __gemm_output_dim_string(options, gemm_args.dims.a) << ","
                 << __gemm_output_dim_string(options, gemm_args.dims.b) << ","
                 << __gemm_output_dim_string(options, gemm_args.dims.c) << "," << options.warm_up_n << "," << options.n
                 << "," << time_in_seconds << "," << time_in_seconds / options.n << "," << flops << ","
                 << gflops / average_time << std::endl;
}

#ifdef PERF_TEST_DEBUG
static void __print_gemm_perf_test_options(options_t options) {
  printf("options.test      = %s\n", test_e_str[options.test].c_str());
  printf("options.loop      = %s\n", loop_e_str[options.loop].c_str());
  printf("options.start     = %dx%d,%dx%d\n", options.start.a.m, options.start.a.n, options.start.b.m,
         options.start.b.n);
  printf("options.stop      = %dx%d,%dx%d\n", options.stop.a.m, options.stop.a.n, options.stop.b.m, options.stop.b.n);
  printf("options.step      = %d\n", options.step);
  printf("options.warm_up_n = %d\n", options.warm_up_n);
  printf("options.n         = %d\n", options.n);
  printf("options.blas_args.gemm.gemm_args = %s\n", options.blas_args.gemm.gemm_args.c_str());
  printf("options.out_file  = %s\n", options.out_file.c_str());
  if (std::is_same<double, default_scalar>::value)
    printf("options.alpha     = %lf\n", options.blas_args.gemm.alpha);
  else if (std::is_same<float, default_scalar>::value)
    printf("options.alpha     = %f\n", options.blas_args.gemm.alpha);
  return;
}
#else
static void __print_gemm_perf_test_options(options_t /*options*/) { return; }
#endif  // PERF_TEST_DEBUG

/*************************** Internal templated fns **************************/
#if !defined(KOKKOS_ENABLE_CUDA)
template <class scalar_type, class vta, class vtb, class device_type>
void __do_gemm_serial_blas(options_t options, gemm_args_t gemm_args) {
  // Need to take subviews on the device
  Kokkos::Timer timer;

  STATUS;

  auto __do_loop = [](uint32_t n, gemm_args_t _gemm_args, bool batch_size_last_dim) {
    for (uint32_t i = 0; i < n; ++i) {
      for (int j = 0; j < _gemm_args.dims.c.k; j++) {
        auto A = Kokkos::subview(_gemm_args.A, j, Kokkos::ALL(), Kokkos::ALL());
        auto B = Kokkos::subview(_gemm_args.B, j, Kokkos::ALL(), Kokkos::ALL());
        auto C = Kokkos::subview(_gemm_args.C, j, Kokkos::ALL(), Kokkos::ALL());
        if (batch_size_last_dim) {
          A = Kokkos::subview(_gemm_args.A, Kokkos::ALL(), Kokkos::ALL(), j);
          B = Kokkos::subview(_gemm_args.B, Kokkos::ALL(), Kokkos::ALL(), j);
          C = Kokkos::subview(_gemm_args.C, Kokkos::ALL(), Kokkos::ALL(), j);
        }

        KokkosBlas::gemm(&_gemm_args.transA, &_gemm_args.transB, _gemm_args.alpha, A, B, _gemm_args.beta, C);
      }
    }
  };
  __do_loop(options.warm_up_n, gemm_args, options.blas_args.batch_size_last_dim);
  Kokkos::fence();

  timer.reset();
  __do_loop(options.n, gemm_args, options.blas_args.batch_size_last_dim);
  Kokkos::fence();

  __gemm_output_csv_row(options, gemm_args, timer.seconds());
  return;
}
#else
template <class scalar_type, class vta, class vtb, class device_type>
void __do_gemm_serial_blas(options_t /*options*/, gemm_args_t /*gemm_args*/) {
  std::cerr << std::string(__func__) << " disabled since KOKKOS_ENABLE_CUDA is defined." << std::endl;
  return;
}
#endif  // !KOKKOS_ENABLE_CUDA

#if !defined(KOKKOS_ENABLE_CUDA)
template <class TransAType, class TransBType, class AlgoType>
void __do_gemm_serial_batched_template(options_t options, gemm_args_t gemm_args) {
  // Need to take subviews on the device
  Kokkos::Timer timer;

  auto __do_loop = [](uint32_t n, gemm_args_t _gemm_args, bool batch_size_last_dim) {
    for (uint32_t i = 0; i < n; ++i) {
      for (int j = 0; j < _gemm_args.dims.c.k; j++) {
        auto A = Kokkos::subview(_gemm_args.A, j, Kokkos::ALL(), Kokkos::ALL());
        auto B = Kokkos::subview(_gemm_args.B, j, Kokkos::ALL(), Kokkos::ALL());
        auto C = Kokkos::subview(_gemm_args.C, j, Kokkos::ALL(), Kokkos::ALL());
        if (batch_size_last_dim) {
          A = Kokkos::subview(_gemm_args.A, Kokkos::ALL(), Kokkos::ALL(), j);
          B = Kokkos::subview(_gemm_args.B, Kokkos::ALL(), Kokkos::ALL(), j);
          C = Kokkos::subview(_gemm_args.C, Kokkos::ALL(), Kokkos::ALL(), j);
        }

        KokkosBatched::SerialGemm<TransAType, TransBType, AlgoType>::invoke(_gemm_args.alpha, A, B, _gemm_args.beta, C);
      }
    }
  };

  __do_loop(options.warm_up_n, gemm_args, options.blas_args.batch_size_last_dim);
  Kokkos::fence();

  timer.reset();
  __do_loop(options.n, gemm_args, options.blas_args.batch_size_last_dim);
  Kokkos::fence();
  __gemm_output_csv_row(options, gemm_args, timer.seconds());
}
#else
template <class TransAType, class TransBType, class AlgoType>
void __do_gemm_serial_batched_template(options_t /*options*/, gemm_args_t /*gemm_args*/) {
  std::cerr << std::string(__func__) << " disabled since KOKKOS_ENABLE_CUDA is defined." << std::endl;
}
#endif  // !KOKKOS_ENABLE_CUDA

template <class scalar_type, class vta, class vtb, class vtc, class device_type, class algo_type>
void __do_gemm_serial_batched(options_t options, gemm_args_t gemm_args) {
  char a  = toupper(gemm_args.transA);
  char b  = toupper(gemm_args.transB);
  using N = KokkosBatched::Trans::NoTranspose;
  using T = KokkosBatched::Trans::Transpose;
  // using C = KokkosBatched::Trans::ConjTranspose;

  STATUS;

  if (a == 'N' && b == 'N') {
    __do_gemm_serial_batched_template<N, N, algo_type>(options, gemm_args);
  } else if (a == 'N' && b == 'T') {
    __do_gemm_serial_batched_template<N, T, algo_type>(options, gemm_args);
    //} else if (a == 'N' && b == 'C') {
    //  __do_gemm_serial_batched_template<N, C, algo_type>(options, gemm_args);
  } else if (a == 'T' && b == 'N') {
    __do_gemm_serial_batched_template<T, N, algo_type>(options, gemm_args);
  } else if (a == 'T' && b == 'T') {
    __do_gemm_serial_batched_template<T, T, algo_type>(options, gemm_args);
    //} else if (a == 'T' && b == 'C') {
    //  __do_gemm_serial_batched_template<T, C, algo_type>(options, gemm_args);
    //} else if (a == 'C' && b == 'N') {
    //  __do_gemm_serial_batched_template<C, N, algo_type>(options, gemm_args);
    //} else if (a == 'C' && b == 'T') {
    //  __do_gemm_serial_batched_template<C, T, algo_type>(options, gemm_args);
    //} else if (a == 'C' && b == 'C') {
    //  __do_gemm_serial_batched_template<C, C, algo_type>(options, gemm_args);
  } else {
    FATAL_ERROR("Bad gemm_args TransA or TransB value");
  }
  return;
}

template <class algo_tag, class blocking_type, class device_type, class algo_mode = void>
void __do_gemm_parallel_batched_heuristic_template(options_t options, gemm_args_t gemm_args) {
  KokkosBatched::BatchedGemmHandle batchedGemmHandle(KokkosBatched::BaseHeuristicAlgos::SQUARE);
  char a  = toupper(gemm_args.transA);
  char b  = toupper(gemm_args.transB);
  using N = KokkosBatched::Trans::NoTranspose;
  using T = KokkosBatched::Trans::Transpose;
  // using C = KokkosBatched::Trans::ConjTranspose;
  using KokkosBatched::BatchLayout;

  STATUS;
  if (a == 'N' && b == 'N') {
    if constexpr (std::is_same_v<typename vector_view_type_3d::array_layout, Kokkos::LayoutLeft>) {
      if (options.use_simd) {
        KokkosBatched::BatchedGemm<N, N, BatchLayout::Right>(&batchedGemmHandle, gemm_args.alpha, gemm_args.Av.vec_3d,
                                                             gemm_args.Bv.vec_3d, gemm_args.beta, gemm_args.Cv.vec_3d);
      } else if constexpr (std::is_same_v<typename view_type_3d::array_layout, Kokkos::LayoutLeft>) {
        KokkosBatched::BatchedGemm<N, N, BatchLayout::Right>(&batchedGemmHandle, gemm_args.alpha, gemm_args.A,
                                                             gemm_args.B, gemm_args.beta, gemm_args.C);
      }
    } else if (options.use_simd) {
      KokkosBatched::BatchedGemm<N, N, BatchLayout::Left>(&batchedGemmHandle, gemm_args.alpha, gemm_args.Av.vec_3d,
                                                          gemm_args.Bv.vec_3d, gemm_args.beta, gemm_args.Cv.vec_3d);
    } else if constexpr (std::is_same_v<typename view_type_3d::array_layout, Kokkos::LayoutRight>) {
      KokkosBatched::BatchedGemm<N, N, BatchLayout::Left>(&batchedGemmHandle, gemm_args.alpha, gemm_args.A, gemm_args.B,
                                                          gemm_args.beta, gemm_args.C);
    }
  } else if (a == 'N' && b == 'T') {
    if constexpr (std::is_same_v<typename vector_view_type_3d::array_layout, Kokkos::LayoutLeft>) {
      if (options.use_simd) {
        KokkosBatched::BatchedGemm<N, T, BatchLayout::Right>(&batchedGemmHandle, gemm_args.alpha, gemm_args.Av.vec_3d,
                                                             gemm_args.Bv.vec_3d, gemm_args.beta, gemm_args.Cv.vec_3d);
      } else if constexpr (std::is_same_v<typename view_type_3d::array_layout, Kokkos::LayoutLeft>) {
        KokkosBatched::BatchedGemm<N, T, BatchLayout::Right>(&batchedGemmHandle, gemm_args.alpha, gemm_args.A,
                                                             gemm_args.B, gemm_args.beta, gemm_args.C);
      }
    } else if (options.use_simd) {
      KokkosBatched::BatchedGemm<N, T, BatchLayout::Left>(&batchedGemmHandle, gemm_args.alpha, gemm_args.Av.vec_3d,
                                                          gemm_args.Bv.vec_3d, gemm_args.beta, gemm_args.Cv.vec_3d);
    } else if constexpr (std::is_same_v<typename view_type_3d::array_layout, Kokkos::LayoutRight>) {
      KokkosBatched::BatchedGemm<N, T, BatchLayout::Left>(&batchedGemmHandle, gemm_args.alpha, gemm_args.A, gemm_args.B,
                                                          gemm_args.beta, gemm_args.C);
    }
    //} else if (a == 'N' && b == 'C') {
    //  __do_gemm_serial_batched_template<N, C, algo_type>(options, gemm_args);
  } else if (a == 'T' && b == 'N') {
    if constexpr (std::is_same_v<typename vector_view_type_3d::array_layout, Kokkos::LayoutLeft>) {
      if (options.use_simd) {
        KokkosBatched::BatchedGemm<T, N, BatchLayout::Right>(&batchedGemmHandle, gemm_args.alpha, gemm_args.Av.vec_3d,
                                                             gemm_args.Bv.vec_3d, gemm_args.beta, gemm_args.Cv.vec_3d);
      } else if constexpr (std::is_same_v<typename view_type_3d::array_layout, Kokkos::LayoutLeft>) {
        KokkosBatched::BatchedGemm<T, N, BatchLayout::Right>(&batchedGemmHandle, gemm_args.alpha, gemm_args.A,
                                                             gemm_args.B, gemm_args.beta, gemm_args.C);
      }
    } else if (options.use_simd) {
      KokkosBatched::BatchedGemm<T, N, BatchLayout::Left>(&batchedGemmHandle, gemm_args.alpha, gemm_args.Av.vec_3d,
                                                          gemm_args.Bv.vec_3d, gemm_args.beta, gemm_args.Cv.vec_3d);
    } else if constexpr (std::is_same_v<typename view_type_3d::array_layout, Kokkos::LayoutRight>) {
      KokkosBatched::BatchedGemm<T, N, BatchLayout::Left>(&batchedGemmHandle, gemm_args.alpha, gemm_args.A, gemm_args.B,
                                                          gemm_args.beta, gemm_args.C);
    }
  } else if (a == 'T' && b == 'T') {
    if constexpr (std::is_same_v<typename vector_view_type_3d::array_layout, Kokkos::LayoutLeft>) {
      if (options.use_simd) {
        KokkosBatched::BatchedGemm<T, T, BatchLayout::Right>(&batchedGemmHandle, gemm_args.alpha, gemm_args.Av.vec_3d,
                                                             gemm_args.Bv.vec_3d, gemm_args.beta, gemm_args.Cv.vec_3d);
      } else if constexpr (std::is_same_v<typename view_type_3d::array_layout, Kokkos::LayoutLeft>) {
        KokkosBatched::BatchedGemm<T, T, BatchLayout::Right>(&batchedGemmHandle, gemm_args.alpha, gemm_args.A,
                                                             gemm_args.B, gemm_args.beta, gemm_args.C);
      }
    } else if (options.use_simd) {
      KokkosBatched::BatchedGemm<T, T, BatchLayout::Left>(&batchedGemmHandle, gemm_args.alpha, gemm_args.Av.vec_3d,
                                                          gemm_args.Bv.vec_3d, gemm_args.beta, gemm_args.Cv.vec_3d);
    } else if constexpr (std::is_same_v<typename view_type_3d::array_layout, Kokkos::LayoutRight>) {
      KokkosBatched::BatchedGemm<T, T, BatchLayout::Left>(&batchedGemmHandle, gemm_args.alpha, gemm_args.A, gemm_args.B,
                                                          gemm_args.beta, gemm_args.C);
    }
    //} else if (a == 'T' && b == 'C') {
    //  __do_gemm_serial_batched_template<T, C, algo_type>(options, gemm_args);
    //} else if (a == 'C' && b == 'N') {
    //  __do_gemm_serial_batched_template<C, N, algo_type>(options, gemm_args);
    //} else if (a == 'C' && b == 'T') {
    //  __do_gemm_serial_batched_template<C, T, algo_type>(options, gemm_args);
    //} else if (a == 'C' && b == 'C') {
    //  __do_gemm_serial_batched_template<C, C, algo_type>(options, gemm_args);
  } else {
    FATAL_ERROR("Bad gemm_args TransA or TransB value");
  }
}

template <class algo_tag, class blocking_type, class device_type, class algo_mode = void>
void __do_gemm_parallel_batched_heuristic(options_t options, gemm_args_t gemm_args) {
  Kokkos::Timer timer;

  for (uint32_t i = 0; i < options.warm_up_n; ++i)
    __do_gemm_parallel_batched_heuristic_template<algo_tag, blocking_type, device_type, algo_mode>(options, gemm_args);
  Kokkos::fence();

  timer.reset();
  for (uint32_t i = 0; i < options.n; ++i)
    __do_gemm_parallel_batched_heuristic_template<algo_tag, blocking_type, device_type, algo_mode>(options, gemm_args);
  Kokkos::fence();

  __gemm_output_csv_row(options, gemm_args, timer.seconds(), nullptr, "-", "-", "-");
}

template <class TransAType, class TransBType, class BlockingType>
struct parallel_batched_gemm_range_policy {
  gemm_args_t gemm_args_;

  parallel_batched_gemm_range_policy(gemm_args_t gemm_args) : gemm_args_(gemm_args) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const SerialTag &, const int &i) const {
    auto svA = Kokkos::subview(gemm_args_.A, i, Kokkos::ALL(), Kokkos::ALL());
    auto svB = Kokkos::subview(gemm_args_.B, i, Kokkos::ALL(), Kokkos::ALL());
    auto svC = Kokkos::subview(gemm_args_.C, i, Kokkos::ALL(), Kokkos::ALL());

    KokkosBatched::SerialGemm<TransAType, TransBType, BlockingType>::invoke(gemm_args_.alpha, svA, svB, gemm_args_.beta,
                                                                            svC);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const SerialBatchDim3Tag &, const int &i) const {
    auto svA = Kokkos::subview(gemm_args_.A, Kokkos::ALL(), Kokkos::ALL(), i);
    auto svB = Kokkos::subview(gemm_args_.B, Kokkos::ALL(), Kokkos::ALL(), i);
    auto svC = Kokkos::subview(gemm_args_.C, Kokkos::ALL(), Kokkos::ALL(), i);

    KokkosBatched::SerialGemm<TransAType, TransBType, BlockingType>::invoke(gemm_args_.alpha, svA, svB, gemm_args_.beta,
                                                                            svC);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const SerialSimdTag &, const int &i) const {
    auto svA = Kokkos::subview(gemm_args_.Av.vec_3d, i, Kokkos::ALL(), Kokkos::ALL());
    auto svB = Kokkos::subview(gemm_args_.Bv.vec_3d, i, Kokkos::ALL(), Kokkos::ALL());
    auto svC = Kokkos::subview(gemm_args_.Cv.vec_3d, i, Kokkos::ALL(), Kokkos::ALL());

    KokkosBatched::SerialGemm<TransAType, TransBType, BlockingType>::invoke(gemm_args_.alpha, svA, svB, gemm_args_.beta,
                                                                            svC);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const SerialSimdBatchDim3Tag &, const int &i) const {
    auto svA = Kokkos::subview(gemm_args_.Av.vec_3d, Kokkos::ALL(), Kokkos::ALL(), i);
    auto svB = Kokkos::subview(gemm_args_.Bv.vec_3d, Kokkos::ALL(), Kokkos::ALL(), i);
    auto svC = Kokkos::subview(gemm_args_.Cv.vec_3d, Kokkos::ALL(), Kokkos::ALL(), i);

    KokkosBatched::SerialGemm<TransAType, TransBType, BlockingType>::invoke(gemm_args_.alpha, svA, svB, gemm_args_.beta,
                                                                            svC);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TeamTag &, const int & /*i*/) const {
    Kokkos::abort("TeamTag not supported using RangePolicy.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TeamBatchDim3Tag &, const int & /*i*/) const {
    Kokkos::abort("TeamBatchDim3Tag not supported using RangePolicy.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TeamVectorTag &, const int & /*i*/) const {
    Kokkos::abort("TeamVectorTag not supported using RangePolicy.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TeamVectorBatchDim3Tag &, const int & /*i*/) const {
    Kokkos::abort("TeamVectorBatchDim3Tag not supported using RangePolicy.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TeamSimdTag &, const int & /*i*/) const {
    Kokkos::abort("TeamSimdTag not supported using RangePolicy.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TeamSimdBatchDim4Tag &, const int & /*i*/) const {
    Kokkos::abort("TeamSimdBatchDim4Tag not supported using RangePolicy.");
  }
};

template <class MemberType, class TransAType, class TransBType, class BlockingType, class AlgoMode = void>
struct parallel_batched_gemm {
  gemm_args_t gemm_args_;

  parallel_batched_gemm(gemm_args_t gemm_args) : gemm_args_(gemm_args) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const SerialTag &, const MemberType &member) const {
    auto i   = member.league_rank();
    auto svA = Kokkos::subview(gemm_args_.A, i, Kokkos::ALL(), Kokkos::ALL());
    auto svB = Kokkos::subview(gemm_args_.B, i, Kokkos::ALL(), Kokkos::ALL());
    auto svC = Kokkos::subview(gemm_args_.C, i, Kokkos::ALL(), Kokkos::ALL());

    KokkosBatched::SerialGemm<TransAType, TransBType, BlockingType>::invoke(gemm_args_.alpha, svA, svB, gemm_args_.beta,
                                                                            svC);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const SerialBatchDim3Tag &, const MemberType &member) const {
    auto i   = member.league_rank();
    auto svA = Kokkos::subview(gemm_args_.A, Kokkos::ALL(), Kokkos::ALL(), i);
    auto svB = Kokkos::subview(gemm_args_.B, Kokkos::ALL(), Kokkos::ALL(), i);
    auto svC = Kokkos::subview(gemm_args_.C, Kokkos::ALL(), Kokkos::ALL(), i);

    KokkosBatched::SerialGemm<TransAType, TransBType, BlockingType>::invoke(gemm_args_.alpha, svA, svB, gemm_args_.beta,
                                                                            svC);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TeamTag &, const MemberType &member) const {
    auto i   = member.league_rank();
    auto svA = Kokkos::subview(gemm_args_.A, i, Kokkos::ALL(), Kokkos::ALL());
    auto svB = Kokkos::subview(gemm_args_.B, i, Kokkos::ALL(), Kokkos::ALL());
    auto svC = Kokkos::subview(gemm_args_.C, i, Kokkos::ALL(), Kokkos::ALL());

    KokkosBatched::TeamGemm<MemberType, TransAType, TransBType, BlockingType>::invoke(member, gemm_args_.alpha, svA,
                                                                                      svB, gemm_args_.beta, svC);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TeamBatchDim3Tag &, const MemberType &member) const {
    auto i   = member.league_rank();
    auto svA = Kokkos::subview(gemm_args_.A, Kokkos::ALL(), Kokkos::ALL(), i);
    auto svB = Kokkos::subview(gemm_args_.B, Kokkos::ALL(), Kokkos::ALL(), i);
    auto svC = Kokkos::subview(gemm_args_.C, Kokkos::ALL(), Kokkos::ALL(), i);

    KokkosBatched::TeamGemm<MemberType, TransAType, TransBType, BlockingType>::invoke(member, gemm_args_.alpha, svA,
                                                                                      svB, gemm_args_.beta, svC);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TeamVectorTag &, const MemberType &member) const {
    auto team_idx = member.league_rank();
    auto svA      = Kokkos::subview(gemm_args_.A, team_idx, Kokkos::ALL(), Kokkos::ALL());
    auto svB      = Kokkos::subview(gemm_args_.B, team_idx, Kokkos::ALL(), Kokkos::ALL());
    auto svC      = Kokkos::subview(gemm_args_.C, team_idx, Kokkos::ALL(), Kokkos::ALL());

    KokkosBatched::TeamVectorGemm<MemberType, TransAType, TransBType, BlockingType>::invoke(
        member, gemm_args_.alpha, svA, svB, gemm_args_.beta, svC);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TeamVectorBatchDim3Tag &, const MemberType &member) const {
    auto team_idx = member.league_rank();
    auto svA      = Kokkos::subview(gemm_args_.A, Kokkos::ALL(), Kokkos::ALL(), team_idx);
    auto svB      = Kokkos::subview(gemm_args_.B, Kokkos::ALL(), Kokkos::ALL(), team_idx);
    auto svC      = Kokkos::subview(gemm_args_.C, Kokkos::ALL(), Kokkos::ALL(), team_idx);

    KokkosBatched::TeamVectorGemm<MemberType, TransAType, TransBType, BlockingType>::invoke(
        member, gemm_args_.alpha, svA, svB, gemm_args_.beta, svC);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TeamSimdTag &, const MemberType &member) const {
    auto i = member.league_rank();
    Kokkos::parallel_for(
        Kokkos::ThreadVectorRange(member, gemm_args_.Cv.ivec_4d.extent(3)), [&](const int &vector_lane) {
          auto svA = Kokkos::subview(gemm_args_.Av.ivec_4d, i, Kokkos::ALL(), Kokkos::ALL(), vector_lane);
          auto svB = Kokkos::subview(gemm_args_.Bv.ivec_4d, i, Kokkos::ALL(), Kokkos::ALL(), vector_lane);
          auto svC = Kokkos::subview(gemm_args_.Cv.ivec_4d, i, Kokkos::ALL(), Kokkos::ALL(), vector_lane);

          KokkosBatched::Gemm<MemberType, TransAType, TransBType, AlgoMode, BlockingType>::invoke(
              member, gemm_args_.alpha, svA, svB, gemm_args_.beta, svC);
        });
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TeamSimdBatchDim4Tag &, const MemberType &member) const {
    auto i = member.league_rank();
    Kokkos::parallel_for(
        Kokkos::ThreadVectorRange(member, gemm_args_.Cv.ivec_4d.extent(0)), [&](const int &vector_lane) {
          auto svA = Kokkos::subview(gemm_args_.Av.ivec_4d, vector_lane, Kokkos::ALL(), Kokkos::ALL(), i);
          auto svB = Kokkos::subview(gemm_args_.Bv.ivec_4d, vector_lane, Kokkos::ALL(), Kokkos::ALL(), i);
          auto svC = Kokkos::subview(gemm_args_.Cv.ivec_4d, vector_lane, Kokkos::ALL(), Kokkos::ALL(), i);

          KokkosBatched::Gemm<MemberType, TransAType, TransBType, AlgoMode, BlockingType>::invoke(
              member, gemm_args_.alpha, svA, svB, gemm_args_.beta, svC);
        });
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const SerialSimdTag &, const MemberType & /*member*/) const {
    Kokkos::abort("SerialSimdTag not supported using RangePolicy.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const SerialSimdBatchDim3Tag &, const MemberType & /*member*/) const {
    Kokkos::abort("SerialSimdBatchDim3Tag not supported using RangePolicy.");
  }
};

template <class TransAType, class TransBType, class BlockingType, class AlgoTag, class device_type>
void __do_gemm_parallel_batched_template_range_policy(options_t options, gemm_args_t gemm_args) {
  using execution_space = typename device_type::execution_space;
  using policy_type     = Kokkos::RangePolicy<AlgoTag, execution_space>;
  using functor_type    = parallel_batched_gemm_range_policy<TransAType, TransBType, BlockingType>;

  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  auto batch_size    = options.start.c.k;
  Kokkos::Timer timer;

  STATUS;

  functor_type parallel_batched_gemm_functor(gemm_args);

  if (std::is_same<AlgoTag, SerialSimdTag>::value || std::is_same<AlgoTag, SerialSimdBatchDim3Tag>::value) {
    batch_size = options.blas_args.batch_size_last_dim ? gemm_args.Cv.vec_3d.extent(2) : gemm_args.Cv.vec_3d.extent(0);
  }

  for (uint32_t i = 0; i < warm_up_n; i++) {
    Kokkos::parallel_for("parallelBatchedWarmUpLoopGemm", policy_type(0, batch_size), parallel_batched_gemm_functor);
    Kokkos::fence();
  }

  timer.reset();
  for (uint32_t i = 0; i < n; i++) {
    Kokkos::parallel_for("parallelBatchedTimedLoopGemm", policy_type(0, batch_size), parallel_batched_gemm_functor);
    Kokkos::fence();
  }

  __gemm_output_csv_row(options, gemm_args, timer.seconds());

  return;
}

template <class TransAType, class TransBType, class BlockingType, class AlgoTag, class device_type,
          class algo_mode = void>
void __do_gemm_parallel_batched_template(options_t options, gemm_args_t gemm_args) {
  using execution_space = typename device_type::execution_space;
  using policy_type     = Kokkos::TeamPolicy<AlgoTag, execution_space>;
  using member_type     = typename policy_type::member_type;
  using functor_type    = parallel_batched_gemm<member_type, TransAType, TransBType, BlockingType, algo_mode>;

  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  auto league_size   = options.start.c.k;
  auto team_size     = gemm_args.bp.team_size;
  auto vector_len    = gemm_args.bp.vector_len;
  Kokkos::Timer timer;

  if (std::is_same<AlgoTag, SerialTag>::value || std::is_same<AlgoTag, SerialBatchDim3Tag>::value ||
      std::is_same<AlgoTag, SerialSimdTag>::value || std::is_same<AlgoTag, SerialSimdBatchDim3Tag>::value) {
    return __do_gemm_parallel_batched_template_range_policy<TransAType, TransBType, BlockingType, AlgoTag, device_type>(
        options, gemm_args);
  }

  if (std::is_same<AlgoTag, TeamSimdTag>::value || std::is_same<AlgoTag, TeamSimdBatchDim4Tag>::value) {
    league_size =
        options.blas_args.batch_size_last_dim ? gemm_args.Cv.ivec_4d.extent(3) : gemm_args.Cv.ivec_4d.extent(0);
    vector_len = simd_vector_size / simd_internal_vector_size;  // TODO: use bp.vector_len?
  }

  STATUS;

  functor_type parallel_batched_gemm_functor(gemm_args);

  if (options.blas_args.use_auto) {
    for (uint32_t i = 0; i < warm_up_n; i++) {
      Kokkos::parallel_for("parallelBatchedWarmUpLoopGemm", policy_type(league_size, Kokkos::AUTO, Kokkos::AUTO),
                           parallel_batched_gemm_functor);
      Kokkos::fence();
    }

    timer.reset();
    for (uint32_t i = 0; i < n; i++) {
      Kokkos::parallel_for("parallelBatchedTimedLoopGemm", policy_type(league_size, Kokkos::AUTO, Kokkos::AUTO),
                           parallel_batched_gemm_functor);
      Kokkos::fence();
    }
  } else {
    for (uint32_t i = 0; i < warm_up_n; i++) {
      Kokkos::parallel_for("parallelBatchedWarmUpLoopGemm", policy_type(league_size, team_size, vector_len),
                           parallel_batched_gemm_functor);
      Kokkos::fence();
    }

    timer.reset();
    for (uint32_t i = 0; i < n; i++) {
      Kokkos::parallel_for("parallelBatchedTimedLoopGemm", policy_type(league_size, team_size, vector_len),
                           parallel_batched_gemm_functor);
      Kokkos::fence();
    }
  }

  __gemm_output_csv_row(options, gemm_args, timer.seconds());

  return;
}

template <class algo_tag, class blocking_type, class device_type, class algo_mode = void>
void __do_gemm_parallel_batched(options_t options, gemm_args_t gemm_args) {
  char a  = gemm_args.transA;
  char b  = gemm_args.transB;
  using N = KokkosBatched::Trans::NoTranspose;
  using T = KokkosBatched::Trans::Transpose;
  // using C = KokkosBatched::Trans::ConjTranspose;

  STATUS;

  if (a == 'N' && b == 'N') {
    __do_gemm_parallel_batched_template<N, N, blocking_type, algo_tag, device_type, algo_mode>(options, gemm_args);
  } else if (a == 'N' && b == 'T') {
    __do_gemm_parallel_batched_template<N, T, blocking_type, algo_tag, device_type, algo_mode>(options, gemm_args);
    //} else if (a == 'N' && b == 'C') {
    //  __do_gemm_parallel_batched_template<N, C, blocking_type, algo_tag,
    //  device_type>(options, gemm_args);
  } else if (a == 'T' && b == 'N') {
    __do_gemm_parallel_batched_template<T, N, blocking_type, algo_tag, device_type, algo_mode>(options, gemm_args);
  } else if (a == 'T' && b == 'T') {
    __do_gemm_parallel_batched_template<T, T, blocking_type, algo_tag, device_type, algo_mode>(options, gemm_args);
    //} else if (a == 'T' && b == 'C') {
    //  __do_gemm_parallel_batched_template<T, C, blocking_type, algo_tag,
    //  device_type>(options, gemm_args);
    //} else if (a == 'C' && b == 'N') {
    //  __do_gemm_parallel_batched_template<C, N, blocking_type, algo_tag,
    //  device_type>(options, gemm_args);
    //} else if (a == 'C' && b == 'T') {
    //  __do_gemm_parallel_batched_template<C, T, blocking_type, algo_tag,
    //  device_type>(options, gemm_args);
    //} else if (a == 'C' && b == 'C') {
    //  __do_gemm_parallel_batched_template<C, C, blocking_type, algo_tag,
    //  device_type>(options, gemm_args);
  } else {
    FATAL_ERROR("Bad gemm_args TransA or TransB value");
  }

  return;
}

template <class TransAType, class TransBType, class BlockingType>
struct parallel_batched_gemm_experiment1 {
  gemm_args_t gemm_args_;

  parallel_batched_gemm_experiment1(gemm_args_t gemm_args) : gemm_args_(gemm_args) {}

  KOKKOS_INLINE_FUNCTION

  void operator()(const SerialTag &, const int &i) const {
    auto svA = Kokkos::subview(gemm_args_.A, i, Kokkos::ALL(), Kokkos::ALL());
    auto svB = Kokkos::subview(gemm_args_.B, i, Kokkos::ALL(), Kokkos::ALL());
    auto svC = Kokkos::subview(gemm_args_.C, i, Kokkos::ALL(), Kokkos::ALL());

    // Uses two serial for-loops internally
    KokkosBatched::SerialGemm<TransAType, TransBType, BlockingType>::invoke(gemm_args_.alpha, svA, svB, gemm_args_.beta,
                                                                            svC);
  }
};

/**
 * 1. parallel_for(rangePolicy<Kokkos::DefaultExecutionSpace>(N)): serialGemm
 *
 */
template <class TransAType, class TransBType, class BlockingType, class device_type>
void __do_gemm_parallel_experiment1(options_t options, gemm_args_t gemm_args) {
  using execution_space = typename device_type::execution_space;
  using policy_type     = Kokkos::RangePolicy<SerialTag, execution_space>;
  using functor_type    = parallel_batched_gemm_experiment1<TransAType, TransBType, BlockingType>;

  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  auto k             = options.start.c.k;
  Kokkos::Timer timer;
  STATUS;

  functor_type experiment1_functor(gemm_args);

  for (uint32_t i = 0; i < warm_up_n; ++i) {
    Kokkos::parallel_for("parallelBatchedUntimedExperiment1Gemm", policy_type(0, k), experiment1_functor);
  }
  Kokkos::fence();

  timer.reset();
  for (uint32_t i = 0; i < n; ++i) {
    Kokkos::parallel_for("parallelBatchedTimedExperiment1Gemm", policy_type(0, k), experiment1_functor);
  }
  Kokkos::fence();

  __gemm_output_csv_row(options, gemm_args, timer.seconds(), "experiment1");
  return;
}

template <class TransAType, class TransBType, class BlockingType, class MemberType>
struct parallel_batched_gemm_experiment2_3_4 {
  gemm_args_t gemm_args_;

  parallel_batched_gemm_experiment2_3_4(gemm_args_t gemm_args) : gemm_args_(gemm_args) {}

  // Experiment 2
  KOKKOS_INLINE_FUNCTION
  void operator()(const TeamVectorTag &, const MemberType &member) const {
    auto i   = member.league_rank();
    auto svA = Kokkos::subview(gemm_args_.A, i, Kokkos::ALL(), Kokkos::ALL());
    auto svB = Kokkos::subview(gemm_args_.B, i, Kokkos::ALL(), Kokkos::ALL());
    auto svC = Kokkos::subview(gemm_args_.C, i, Kokkos::ALL(), Kokkos::ALL());

    // Uses TeamThreadRange over C-rows
    //        ThreadVectorRange over C-cols
    KokkosBatched::TeamVectorGemm<MemberType, TransAType, TransBType, BlockingType>::invoke(
        member, gemm_args_.alpha, svA, svB, gemm_args_.beta, svC);
  }

  // Experiment 3
  KOKKOS_INLINE_FUNCTION
  void operator()(const LayoutLeftTag &, const MemberType &member) const {
    auto team_idx = member.league_rank();
    auto svA      = Kokkos::subview(gemm_args_.A, team_idx, Kokkos::ALL(), Kokkos::ALL());
    auto svB      = Kokkos::subview(gemm_args_.B, team_idx, Kokkos::ALL(), Kokkos::ALL());
    auto svC      = Kokkos::subview(gemm_args_.C, team_idx, Kokkos::ALL(), Kokkos::ALL());

    // TeamThreadRange:   splits the index range over the threads of the team
    // ThreadVectorRange: splits the index range over the vector lanes of the
    // calling thread

    auto svC_cols = svC.extent(1);
    // In a given team, for each vector lane, compute zero or more output
    // columns of C depending on the index range
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, svC_cols), [&](const int &lane_idx) {
      auto svB_col = Kokkos::subview(svB, Kokkos::ALL(), lane_idx);
      auto svC_col = Kokkos::subview(svC, Kokkos::ALL(), lane_idx);
      // TeamGemm Calls TeamThreadRange over M*N meaning the flat M*N array
      // is split over all threads of the team
      KokkosBatched::TeamGemm<MemberType, TransAType, TransBType, BlockingType>::invoke(
          member, gemm_args_.alpha, svA, svB_col, gemm_args_.beta, svC_col);
    });
  }

  // TODO: Why is this faster than the LayoutLeftTag operator above for both
  // LayoutLeft and LayoutRight? Experiment 4
  KOKKOS_INLINE_FUNCTION
  void operator()(const LayoutRightTag &, const MemberType &member) const {
    auto team_idx = member.league_rank();
    auto svA      = Kokkos::subview(gemm_args_.A, team_idx, Kokkos::ALL(), Kokkos::ALL());
    auto svB      = Kokkos::subview(gemm_args_.B, team_idx, Kokkos::ALL(), Kokkos::ALL());
    auto svC      = Kokkos::subview(gemm_args_.C, team_idx, Kokkos::ALL(), Kokkos::ALL());

    // TeamThreadRange:   splits the index range over the threads of the team
    // ThreadVectorRange: splits the index range over the vector lanes of the
    // calling thread

    auto svC_rows = svC.extent(0);
    // In a given team, for each vector lane, compute zero or more output rows
    // of C depending on the index range
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, svC_rows), [&](const int &lane_idx) {
      auto svA_row = Kokkos::subview(svA, lane_idx, Kokkos::ALL());
      auto svC_row = Kokkos::subview(svC, lane_idx, Kokkos::ALL());
      // TeamGemm Calls TeamThreadRange over M*N meaning the flat M*N array
      // is split over all threads of the team
      KokkosBatched::TeamGemm<MemberType, TransAType, TransBType, BlockingType>::invoke(
          member, gemm_args_.alpha, svA_row, svB, gemm_args_.beta, svC_row);
    });
  }
};

/**
 * 2. case a)
 * parallel_for(teamPolicy): TeamVectorGemm
 *
 */
template <class TransAType, class TransBType, class BlockingType, class device_type>
void __do_gemm_parallel_experiment2(options_t options, gemm_args_t gemm_args) {
  using execution_space = typename device_type::execution_space;
  using policy_type     = Kokkos::TeamPolicy<TeamVectorTag, execution_space>;
  using member_type     = typename policy_type::member_type;
  using functor_type    = parallel_batched_gemm_experiment2_3_4<TransAType, TransBType, BlockingType, member_type>;

  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  auto league_size   = options.start.c.k;
  Kokkos::Timer timer;
  STATUS;

  functor_type experiment2_functor(gemm_args);

  auto team_size  = gemm_args.bp.team_size;
  auto vector_len = gemm_args.bp.vector_len;

  for (uint32_t i = 0; i < warm_up_n; ++i) {
    Kokkos::parallel_for("parallelBatchedUntimedExperiment2Gemm", policy_type(league_size, team_size, vector_len),
                         experiment2_functor);
  }
  Kokkos::fence();

  timer.reset();

  for (uint32_t i = 0; i < n; ++i) {
    Kokkos::parallel_for("parallelBatchedTimedExperiment2Gemm", policy_type(league_size, team_size, vector_len),
                         experiment2_functor);
  }
  Kokkos::fence();

  __gemm_output_csv_row(options, gemm_args, timer.seconds(), "experiment2");
  return;
}

/**
 * 3. case b)
 *    parallel_for(teamPolicy):
 *      parallel_for(TeamThreadRange):
 *         VectorGemm
 *
 * VectorGemm has not been implemented!
 * I think this experiment can be removed. TeamGemm calls TeamThreadRange
 * internally! TeamVectorGemm calls both TeamThreadRange and ThreadVectorRange
 * internally!
 */
template <class TransAType, class TransBType, class BlockingType, class device_type>
void __do_gemm_parallel_experiment3(options_t options, gemm_args_t gemm_args) {
  using execution_space = typename device_type::execution_space;
  // using layout_tag = std::conditional<std::is_same<default_layout,
  // Kokkos::LayoutLeft>::value, LayoutLeftTag, LayoutRightTag>::type;
  using policy_type  = Kokkos::TeamPolicy<LayoutLeftTag, execution_space>;
  using member_type  = typename policy_type::member_type;
  using functor_type = parallel_batched_gemm_experiment2_3_4<TransAType, TransBType, BlockingType, member_type>;

  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  auto league_size   = options.start.c.k;
  Kokkos::Timer timer;
  STATUS;

  functor_type experiment3_functor(gemm_args);

  auto team_size  = gemm_args.bp.team_size;
  auto vector_len = gemm_args.bp.vector_len;

  for (uint32_t i = 0; i < warm_up_n; ++i) {
    Kokkos::parallel_for("parallelBatchedUntimedExperiment3Gemm", policy_type(league_size, team_size, vector_len),
                         experiment3_functor);
  }
  Kokkos::fence();

  timer.reset();

  for (uint32_t i = 0; i < n; ++i) {
    Kokkos::parallel_for("parallelBatchedTimedExperiment3Gemm", policy_type(league_size, team_size, vector_len),
                         experiment3_functor);
  }
  Kokkos::fence();

  __gemm_output_csv_row(options, gemm_args, timer.seconds(), "experiment3");
  return;
}

/**
 * 4. case c)
 * parallel_for(teamPolicy):
 *      parallel_for(ThreadVectorRange)
 *        TeamGemm
 */
template <class TransAType, class TransBType, class BlockingType, class device_type>
void __do_gemm_parallel_experiment4(options_t options, gemm_args_t gemm_args) {
  using execution_space = typename device_type::execution_space;
  // using layout_tag = std::conditional<std::is_same<default_layout,
  // Kokkos::LayoutLeft>::value, LayoutLeftTag, LayoutRightTag>::type;
  using policy_type  = Kokkos::TeamPolicy<LayoutRightTag, execution_space>;
  using member_type  = typename policy_type::member_type;
  using functor_type = parallel_batched_gemm_experiment2_3_4<TransAType, TransBType, BlockingType, member_type>;

  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  auto league_size   = options.start.c.k;
  Kokkos::Timer timer;
  STATUS;

  functor_type experiment4_functor(gemm_args);

  auto team_size  = gemm_args.bp.team_size;
  auto vector_len = gemm_args.bp.vector_len;

  for (uint32_t i = 0; i < warm_up_n; ++i) {
    Kokkos::parallel_for("parallelBatchedUntimedExperiment4Gemm", policy_type(league_size, team_size, vector_len),
                         experiment4_functor);
  }
  Kokkos::fence();

  timer.reset();

  for (uint32_t i = 0; i < n; ++i) {
    Kokkos::parallel_for("parallelBatchedTimedExperiment4Gemm", policy_type(league_size, team_size, vector_len),
                         experiment4_functor);
  }
  Kokkos::fence();

  __gemm_output_csv_row(options, gemm_args, timer.seconds(), "experiment4");
  return;
}

template <class SimdViewType, class TransAType, class TransBType, class BlockingType>
class parallel_batched_gemm_experiment5 {
 private:
  SimdViewType &A, &B, &C;
  gemm_args_t gemm_args;

 public:
  parallel_batched_gemm_experiment5(SimdViewType &_A, SimdViewType &_B, SimdViewType &_C, gemm_args_t _gemm_args)
      : A(_A), B(_B), C(_C), gemm_args(_gemm_args) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const SimdCpuTag &, const int &i) const {
    auto svA = Kokkos::subview(A, i, Kokkos::ALL(), Kokkos::ALL());
    auto svB = Kokkos::subview(B, i, Kokkos::ALL(), Kokkos::ALL());
    auto svC = Kokkos::subview(C, i, Kokkos::ALL(), Kokkos::ALL());

    // Uses two serial for-loops internally
    KokkosBatched::SerialGemm<TransAType, TransBType, BlockingType>::invoke(gemm_args.alpha, svA, svB, gemm_args.beta,
                                                                            svC);
  }
};

/**
 * 5.
 * parallel_for(RangePolicy<Kokkos:DefaultHostExecutionSpace>(N/vl+(N%vl>0)>):
 * serialGemm
 *
 * Not portable to GPU
 */
#if !defined(KOKKOS_ENABLE_CUDA) && !defined(KOKKOS_ENABLE_HIP)
template <class TransAType, class TransBType, class BlockingType, class device_type>
void __do_gemm_parallel_experiment5(options_t options, gemm_args_t gemm_args) {
  using execution_space = typename device_type::execution_space;
  using policy_type     = Kokkos::RangePolicy<SimdCpuTag, execution_space>;

  // Construct the SimdType
  using scalar_type    = typename view_type_3d::value_type;
  constexpr int vl     = KokkosBatched::DefaultVectorLength<scalar_type, execution_space>::value;
  using simd_type      = KokkosBatched::Vector<KokkosBatched::SIMD<scalar_type>, simd_vector_size>;
  using simd_view_type = Kokkos::View<simd_type ***, default_layout, default_device>;
  using functor_type   = parallel_batched_gemm_experiment5<simd_view_type, TransAType, TransBType, BlockingType>;

  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  auto k             = options.start.c.k;
  Kokkos::Timer timer;
  auto simd_batch_size = k / vl + (k % vl > 0);
  STATUS;

  // Increases each array size by sizeof(scalar_type) * (vl-1) bytes!
  simd_view_type A("A", simd_batch_size, gemm_args.A.extent(0), gemm_args.A.extent(1));
  simd_view_type B("B", simd_batch_size, gemm_args.B.extent(0), gemm_args.B.extent(1));
  simd_view_type C("C", simd_batch_size, gemm_args.C.extent(0), gemm_args.C.extent(1));

  // uint64_t seed =
  //     std::chrono::high_resolution_clock::now().time_since_epoch().count();
  // Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(seed);
  // Kokkos::fill_random(A, rand_pool,
  // Kokkos::rand<Kokkos::Random_XorShift64<execution_space>,
  // simd_type>::max()); Kokkos::fill_random(B, rand_pool,
  // Kokkos::rand<Kokkos::Random_XorShift64<execution_space>,
  // simd_type>::max()); Kokkos::fill_random(C, rand_pool,
  // Kokkos::rand<Kokkos::Random_XorShift64<execution_space>,
  // simd_type>::max()); execution_space::fence();

  functor_type experiment5_functor(A, B, C, gemm_args);

  for (uint32_t i = 0; i < warm_up_n; ++i) {
    Kokkos::parallel_for("parallelBatchedUntimedExperiment5Gemm", policy_type(0, simd_batch_size), experiment5_functor);
  }
  Kokkos::fence();

  timer.reset();

  for (uint32_t i = 0; i < n; ++i) {
    Kokkos::parallel_for("parallelBatchedTimedExperiment5Gemm", policy_type(0, simd_batch_size), experiment5_functor);
  }
  Kokkos::fence();

  __gemm_output_csv_row(options, gemm_args, timer.seconds(), "experiment5");
  return;
}
#else
template <class TransAType, class TransBType, class BlockingType, class device_type>
void __do_gemm_parallel_experiment5(options_t /*options*/, gemm_args_t /*gemm_args*/) {
  std::cerr << std::string(__func__) << " disabled since KOKKOS_ENABLE_CUDA or KOKKOS_ENABLE_HIP is defined."
            << std::endl;
  return;
}
#endif  // !KOKKOS_ENABLE_CUDA || !KOKKOS_ENABLE_HIP

template <class MemberType, class SimdViewType, class TransAType, class TransBType, class BlockingType>
class parallel_batched_gemm_experiment6 {
 private:
  SimdViewType &A, &B, &C;
  gemm_args_t gemm_args;

 public:
  parallel_batched_gemm_experiment6(SimdViewType &_A, SimdViewType &_B, SimdViewType &_C, gemm_args_t _gemm_args)
      : A(_A), B(_B), C(_C), gemm_args(_gemm_args) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const MemberType &member) const {
    auto i   = member.league_rank();
    auto svA = Kokkos::subview(A, i, Kokkos::ALL(), Kokkos::ALL());
    auto svB = Kokkos::subview(B, i, Kokkos::ALL(), Kokkos::ALL());
    auto svC = Kokkos::subview(C, i, Kokkos::ALL(), Kokkos::ALL());

    // Uses two serial for-loops internally
    KokkosBatched::TeamVectorGemm<MemberType, TransAType, TransBType, BlockingType>::invoke(
        member, gemm_args.alpha, svA, svB, gemm_args.beta, svC);
  }
};

#if 0
template <class TransAType, class TransBType, class BlockingType,
          class device_type>
void __do_gemm_parallel_experiment6(options_t options, gemm_args_t gemm_args) {
  using execution_space = typename device_type::execution_space;
  using policy_type     = Kokkos::TeamPolicy<execution_space>;
  using member_type     = typename policy_type::member_type;

  // Construct the vector type
  using scalar_type = typename view_type_3d::value_type;
  constexpr int vl =
      KokkosBatched::DefaultVectorLength<scalar_type, execution_space>::value;
  constexpr int il =
      KokkosBatched::DefaultInternalVectorLength<scalar_type, execution_space>::value;
  using view_type = Kokkos::View<scalar_type***[vl], default_layout, default_device>;
  using vector_view_type = Kokkos::View<vector_type***, default_layout, default_device>;
  using internal_vector_view_type = Kokkos::View<internal_vector_type***, default_layout, default_device>;
  using functor_type =
      parallel_batched_gemm_experiment6<member_type, internal_vector_view_type,
                                        TransAType, TransBType, BlockingType>;

  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  auto k             = options.start.c.k;
  Kokkos::Timer timer;
  auto simd_batch_size = k / vl + (k % vl > 0);
  STATUS;

  // Construct matrices
  vector_view_type A_vector("A_vector", simd_batch_size, gemm_args.A.extent(0), gemm_args.A.extent(1));
  view_type A((scalar_type *)A_vector.data(), simd_batch_size, gemm_args.A.extent(0), gemm_args.A.extent(1));
  internal_vector_view_type A_vector_internal(A_vector.data(), simd_batch_size, gemm_args.A.extent(0), gemm_args.A.extent(1));

  vector_view_type B_vector("B_vector", simd_batch_size, gemm_args.B.extent(0), gemm_args.B.extent(1));
  view_type B((scalar_type *)B_vector.data(), simd_batch_size, gemm_args.B.extent(0), gemm_args.B.extent(1));
  internal_vector_view_type B_vector_internal(B_vector.data(), simd_batch_size, gemm_args.B.extent(0), gemm_args.B.extent(1));

  vector_view_type C_vector("C_vector", simd_batch_size, gemm_args.C.extent(0), gemm_args.C.extent(1));
  view_type C((scalar_type *)C_vector.data(), simd_batch_size, gemm_args.C.extent(0), gemm_args.C.extent(1));
  internal_vector_view_type C_vector_internal(C_vector.data(), simd_batch_size, gemm_args.C.extent(0), gemm_args.C.extent(1));

  uint64_t seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(seed);
  Kokkos::fill_random(A, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, scalar_type>::max());
  Kokkos::fill_random(B, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, scalar_type>::max());
  Kokkos::fill_random(C, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, scalar_type>::max());
  Kokkos::fence();

  functor_type experiment6_functor(A_vector_internal, B_vector_internal, C_vector_internal, gemm_args);

  for (uint32_t i = 0; i < warm_up_n; ++i) {
    Kokkos::parallel_for("parallelBatchedUntimedExperiment6Gemm",
                         policy_type(simd_batch_size, Kokkos::AUTO, vl/il), experiment6_functor);
    Kokkos::fence();
  }

  timer.reset();
  for (uint32_t i = 0; i < n; ++i) {
    Kokkos::parallel_for("parallelBatchedTimedExperiment6Gemm",
                         policy_type(simd_batch_size, Kokkos::AUTO, vl/il), experiment6_functor);
    Kokkos::fence();
  }

  __gemm_output_csv_row(options, gemm_args, timer.seconds(), "experiment6");
  return;
}
#else
template <class TransAType, class TransBType, class BlockingType, class device_type>
void __do_gemm_parallel_experiment6(options_t /*options*/, gemm_args_t /*gemm_args*/) {
  return;
}
#endif

/**
 * examples/armpl_dgemm_interleave_batch_c_example.c was used as a reference
 * when writing this.
 **/
template <class TransAType, class TransBType, class BlockingType, class device_type>
void __do_gemm_armpl(options_t options, gemm_args_t gemm_args) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ARMPL) && ARMPL_BUILD >= 1058
  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  Kokkos::Timer timer;
  char transa = std::is_same<TransAType, KokkosBatched::Trans::NoTranspose>::value ? 'N' : 'T';
  char transb = std::is_same<TransBType, KokkosBatched::Trans::NoTranspose>::value ? 'N' : 'T';

  if (!std::is_same<default_scalar, double>::value) FATAL_ERROR("only double scalars are supported!");

  STATUS;

  for (uint32_t i = 0; i < warm_up_n; i++) {
    armpl_dgemm_interleave_batch(gemm_args.ninter, gemm_args.nbatch, transa, transb, gemm_args.dims.c.m,
                                 gemm_args.dims.c.n, gemm_args.dims.a.n, gemm_args.alpha, gemm_args.A_pl.mat,
                                 gemm_args.A_pl.bstrd, gemm_args.A_pl.istrd, gemm_args.A_pl.jstrd, gemm_args.B_pl.mat,
                                 gemm_args.B_pl.bstrd, gemm_args.B_pl.istrd, gemm_args.B_pl.jstrd, gemm_args.beta,
                                 gemm_args.C_pl.mat, gemm_args.C_pl.bstrd, gemm_args.C_pl.istrd, gemm_args.C_pl.jstrd);
  }

  timer.reset();
  for (uint32_t i = 0; i < n; i++) {
    armpl_dgemm_interleave_batch(gemm_args.ninter, gemm_args.nbatch, transa, transb, gemm_args.dims.c.m,
                                 gemm_args.dims.c.n, gemm_args.dims.a.n, gemm_args.alpha, gemm_args.A_pl.mat,
                                 gemm_args.A_pl.bstrd, gemm_args.A_pl.istrd, gemm_args.A_pl.jstrd, gemm_args.B_pl.mat,
                                 gemm_args.B_pl.bstrd, gemm_args.B_pl.istrd, gemm_args.B_pl.jstrd, gemm_args.beta,
                                 gemm_args.C_pl.mat, gemm_args.C_pl.bstrd, gemm_args.C_pl.istrd, gemm_args.C_pl.jstrd);
  }

  __gemm_output_csv_row(options, gemm_args, timer.seconds());
#else
  // Cast to void to supress unused param warnings
  (void)options;
  (void)gemm_args;
  std::cerr << std::string(__func__) << " disabled since KOKKOSKERNELS_ENABLE_TPL_ARMPL is undefined." << std::endl;
#endif  // KOKKOSKERNELS_ENABLE_TPL_ARMPL
  return;
}

/**
 * Check difference of scalars expected and actual at indexes i,j,k
 * @var expected: The expected result.
 * @var actual:   The actual result.
 * @var epsilon:  The tolerance to use when comparing.
 * @return true if the comparison fails and false if the comparison succeeds.
 */
template <class ViewType>
static inline bool __gemm_print_compare_failure(ViewType h_expected, ViewType h_actual, int i, int j, int k,
                                                double epsilon) {
  STATUS;
  auto diff = std::fabs(static_cast<double>(h_expected(i, j, k) - h_actual(i, j, k)));

  if (diff > epsilon) {
    printf("fabs(expected(%d,%d,%d):%g - actual(%d,%d,%d):%g):%g > epsilon:%g\n", i, j, k,
           static_cast<double>(h_expected(i, j, k)), i, j, k, static_cast<double>(h_actual(i, j, k)), diff, epsilon);
    FATAL_ERROR("Comparison failure!");
    return true;
  }
  return false;
}

/**
 * Compare all values of expected with all values of actual.
 * @var expected: the expected results
 * @var actual:   the actual results
 * @return false if expected matches actual within epsilon, otherwise true.
 */
template <class ScalarType, class LayoutType>
static inline bool __gemm_do_compare(view_type_3d expected, view_type_3d actual) {
  double epsilon = Test::epsilon<ScalarType>::value * 1e3;
  STATUS;

  typename view_type_3d::HostMirror h_expected = Kokkos::create_mirror_view(expected);
  typename view_type_3d::HostMirror h_actual   = Kokkos::create_mirror_view(actual);

  // Copy to host for comparision
  Kokkos::deep_copy(h_expected, expected);
  Kokkos::deep_copy(h_actual, actual);
  Kokkos::fence();

  if (std::is_same<LayoutType, Kokkos::LayoutRight>::value) {
    for (size_t i = 0; i < h_expected.extent(0); i++) {
      for (size_t j = 0; j < h_expected.extent(1); j++) {
        for (size_t k = 0; k < h_expected.extent(2); k++) {
          if (__gemm_print_compare_failure<decltype(h_expected)>(h_expected, h_actual, i, j, k, epsilon)) return true;
        }
      }
    }
  }

  if (std::is_same<LayoutType, Kokkos::LayoutLeft>::value) {
    for (size_t k = 0; k < h_expected.extent(2); k++) {
      for (size_t j = 0; j < h_expected.extent(1); j++) {
        for (size_t i = 0; i < h_expected.extent(0); i++) {
          if (__gemm_print_compare_failure<decltype(h_expected)>(h_expected, h_actual, i, j, k, epsilon)) return true;
        }
      }
    }
  }

  return false;
}

template <class dstViewType>
static inline void __gemm_copy_simd_view_to_3d_view(gemm_simd_args_t src, dstViewType dst, options_t options) {
  // clang-format off
  // Related issue: https://github.com/kokkos/kokkos-kernels/issues/998
  //   CUDA VERSION 10.2.2 generates a compiler error:
  //     KokkosBlas3_gemm_perf_test.hpp: error: ‘h_subview_type_2d’ was not declared in this scope
  // clang-format on
#if (CUDA_VERSION != 10020)
  using dst_scalar_type = typename dstViewType::value_type;
  using src_scalar_type = typename view_type_5d::value_type;
  size_t remainder, vector_batch_size, simd_batch_size, last_batch;
  bool data_layout_same_as_3d_view        = false;
  typename dstViewType::HostMirror h_dst  = Kokkos::create_mirror_view(dst);
  typename view_type_4d::HostMirror h_src = Kokkos::create_mirror_view(src.mat_4d);
  Kokkos::deep_copy(h_src, src.mat_4d);
  Kokkos::fence();

  if (options.blas_args.batch_size_last_dim) {
    remainder         = dst.extent(2) % simd_internal_vector_size;
    vector_batch_size = src.ivec_4d.extent(0);
    simd_batch_size   = src.ivec_4d.extent(3);
    last_batch        = dst.extent(2);
    if (std::is_same<default_layout, Kokkos::LayoutRight>::value && remainder == 0) data_layout_same_as_3d_view = true;

  } else {
    remainder         = dst.extent(0) % simd_internal_vector_size;
    vector_batch_size = src.ivec_4d.extent(3);
    simd_batch_size   = src.ivec_4d.extent(0);
    last_batch        = dst.extent(0);
    if (std::is_same<default_layout, Kokkos::LayoutLeft>::value && remainder == 0) data_layout_same_as_3d_view = true;
  }

  // When the batch_size is a multiple of the simd_vector_size and the
  // batch_size dimension is nearest to the simd_vector_size dimension, each
  // 2-rank matrix lies in the correct location and the data can simply be cast
  // to the 3d view.
  if (data_layout_same_as_3d_view) {
    // We can just re-cast the data to the 3d view but we'll copy it for
    // verification
    memcpy(h_dst.data(), h_src.data(), sizeof(dst_scalar_type) * dst.extent(0) * dst.extent(1) * dst.extent(2));
    Kokkos::deep_copy(dst, h_dst);
    Kokkos::fence();
    return;
  }

  // If the remainder is 0, we have simd_vector_size sub-batches to copy out...
  // this is a bad data access pattern but for these perf_tests we will support
  // it. If the remainder is non-zero, we have simd_vector_size sub-batches +
  // remainder to copy out.
  remainder += simd_internal_vector_size;

  // Views needed for slow manual copy
  using h_view_type_5d    = Kokkos::View<src_scalar_type *****, default_layout, Kokkos::HostSpace>;
  using h_subview_type_2d = Kokkos::View<src_scalar_type **, Kokkos::LayoutStride, Kokkos::HostSpace>;
  using h_subview_type_3d = Kokkos::View<src_scalar_type ***, Kokkos::LayoutStride, Kokkos::HostSpace>;
  using h_subview_type_4d = Kokkos::View<src_scalar_type ****, Kokkos::LayoutStride, Kokkos::HostSpace>;
  h_view_type_5d h_src_raw;
  h_subview_type_4d h_sv0;
  h_subview_type_3d h_sv1;
  h_subview_type_2d h_sv2;

  // TODO: Clean everything below this point up...
  if (std::is_same<default_layout, Kokkos::LayoutRight>::value)
    h_src_raw = h_view_type_5d((src_scalar_type *)h_src.data(), src.ivec_4d.extent(0), src.ivec_4d.extent(1),
                               src.ivec_4d.extent(2), src.ivec_4d.extent(3), simd_internal_vector_size);
  else
    h_src_raw = h_view_type_5d((src_scalar_type *)h_src.data(), simd_internal_vector_size, src.ivec_4d.extent(0),
                               src.ivec_4d.extent(1), src.ivec_4d.extent(2), src.ivec_4d.extent(3));

  // The below loops copies each corresponding 2-rank matrix within the simd
  // view back to the 3-rank view.
  for (size_t simd_internal_vec_idx = 0; simd_internal_vec_idx < remainder; simd_internal_vec_idx++) {
    if (std::is_same<default_layout, Kokkos::LayoutRight>::value)
      h_sv0 =
          Kokkos::subview(h_src_raw, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), simd_internal_vec_idx);
    else
      h_sv0 =
          Kokkos::subview(h_src_raw, simd_internal_vec_idx, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());

    for (size_t vector_batch_idx = 0; vector_batch_idx < vector_batch_size; vector_batch_idx++) {
      if (options.blas_args.batch_size_last_dim)
        h_sv1 = Kokkos::subview(h_sv0, vector_batch_idx, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
      else
        h_sv1 = Kokkos::subview(h_sv0, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), vector_batch_idx);
      for (size_t simd_batch_size_idx = 0; simd_batch_size_idx < simd_batch_size; simd_batch_size_idx++) {
        if (options.blas_args.batch_size_last_dim)
          h_sv2 = Kokkos::subview(h_sv1, Kokkos::ALL(), Kokkos::ALL(), simd_batch_size_idx);
        else
          h_sv2 = Kokkos::subview(h_sv1, simd_batch_size_idx, Kokkos::ALL(), Kokkos::ALL());
        for (size_t m = 0; m < src.ivec_4d.extent(1); m++) {
          for (size_t n = 0; n < src.ivec_4d.extent(2); n++) {
            if (options.blas_args.batch_size_last_dim)
              h_dst(m, n, simd_internal_vec_idx + simd_batch_size_idx + vector_batch_idx) = h_sv2(m, n);
            else
              h_dst(simd_internal_vec_idx + simd_batch_size_idx + vector_batch_idx, m, n) = h_sv2(m, n);
          }
        }
        if (simd_internal_vec_idx + simd_batch_size_idx + vector_batch_idx == last_batch - 1) goto out;
      }
    }
  }
out:
  Kokkos::deep_copy(dst, h_dst);
  Kokkos::fence();
#else
  // Avoid unused parameter warnings:
  (void)src;
  (void)dst;
  (void)options;

  Kokkos::abort("Cannot perform simd verification with cuda/10.2.2, rerun with -v 0");
#endif  // #if (CUDA_VERSION != 10020)
}

/**
 * Compare all values of expected with all values of actual.
 * @var expected: the expected results
 * @var actual:   the actual results
 * @return false if expected matches actual within epsilon, otherwise true.
 */
template <class ScalarType, class LayoutType>
static inline bool __gemm_do_compare(view_type_3d expected, gemm_simd_args_t actual, options_t options) {
  decltype(expected) actual_data("actual_data", expected.extent(0), expected.extent(1), expected.extent(2));

  STATUS;

  // Copy the simd view to a 3d view for comparision.
  // NOTE: The raw results are different when batch_size % simd_vector_size !=
  // 0. Also note that when batch_size % simd_vector_size != 0, the simd
  // operation calculates results that we do not require. So, we end up running
  // an extra batch_size % simd_vector_size GEMMs!
  __gemm_copy_simd_view_to_3d_view(actual, actual_data, options);
  return __gemm_do_compare<ScalarType, LayoutType>(expected, actual_data);
}

template <class ScalarType, class LayoutType, class DeviceType>
static inline void __gemm_do_verify(options_t options, gemm_args_t gemm_args, void (*fn)(options_t, gemm_args_t)) {
  using execution_space = typename DeviceType::execution_space;
  // Just create "expected" types using non-simd types.
  decltype(gemm_args.C) C_expected;
  decltype(gemm_args.A) A_expected;
  decltype(gemm_args.B) B_expected;
  STATUS;

  if (options.blas_args.batch_size_last_dim) {
    C_expected = decltype(C_expected)("C_expected", gemm_args.dims.c.m, gemm_args.dims.c.n, gemm_args.dims.c.k);
    A_expected = decltype(A_expected)("A_expected", gemm_args.dims.a.m, gemm_args.dims.a.n, gemm_args.dims.a.k);
    B_expected = decltype(B_expected)("B_expected", gemm_args.dims.b.m, gemm_args.dims.b.n, gemm_args.dims.b.k);
  } else {
    C_expected = decltype(C_expected)("C_expected", gemm_args.dims.c.k, gemm_args.dims.c.m, gemm_args.dims.c.n);
    A_expected = decltype(A_expected)("A_expected", gemm_args.dims.a.k, gemm_args.dims.a.m, gemm_args.dims.a.n);
    B_expected = decltype(B_expected)("B_expected", gemm_args.dims.b.k, gemm_args.dims.b.m, gemm_args.dims.b.n);
  }

  // Initialize "expected" matrices.
  if (gemm_args.C.data() != nullptr) {
    Kokkos::deep_copy(C_expected, gemm_args.C);
    Kokkos::deep_copy(A_expected, gemm_args.A);
    Kokkos::deep_copy(B_expected, gemm_args.B);

    Kokkos::fence();  // Ensure that deep_copy has completed

    // Check that initial values match
    if (__gemm_do_compare<ScalarType, LayoutType>(C_expected, gemm_args.C)) FATAL_ERROR("Inital values mismatch!");
  } else if (gemm_args.Cv.vec_3d.data() != nullptr) {
    __gemm_copy_simd_view_to_3d_view<decltype(C_expected)>(gemm_args.Cv, C_expected, options);
    __gemm_copy_simd_view_to_3d_view<decltype(A_expected)>(gemm_args.Av, A_expected, options);
    __gemm_copy_simd_view_to_3d_view<decltype(B_expected)>(gemm_args.Bv, B_expected, options);

    // Check that initial values match
    if (__gemm_do_compare<ScalarType, LayoutType>(C_expected, gemm_args.Cv, options))
      FATAL_ERROR("Inital values mismatch!");
  } else {
    FATAL_ERROR("Input arguments are empty!");
  }

  // Populate "expected" matrices via VanillaGemm
  Test::Functor_BatchedVanillaGEMM<decltype(A_expected), decltype(B_expected), decltype(C_expected), execution_space>
      vgemm;
  vgemm.A_t = toupper(gemm_args.transA) == 'T';
  vgemm.B_t = toupper(gemm_args.transB) == 'T';
  vgemm.A_c = vgemm.B_c     = false;
  vgemm.batch_size_last_dim = options.blas_args.batch_size_last_dim;
  vgemm.A                   = A_expected;
  vgemm.B                   = B_expected;
  vgemm.C                   = C_expected;
  vgemm.alpha               = gemm_args.alpha;
  vgemm.beta                = gemm_args.beta;
  vgemm.run();  // Compute C_expected

  // Run routine with warm_up_n = 1 and n = 0.
  auto warm_up_n_bak = options.warm_up_n;
  options.warm_up_n  = 1;
  auto n_bak         = options.n;
  options.n          = 0;
  fn(options, gemm_args);

  Kokkos::fence();  // Redundant fence.

  // Check the result
  if (gemm_args.C.data() != nullptr) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_ARMPL) && ARMPL_BUILD >= 1058
    if (options.test == EXPERIMENT) {
      using view_type_2d = Kokkos::View<default_scalar **, Kokkos::LayoutStride, default_device>;
      view_type_2d C;
      for (int ib = 0; ib < gemm_args.nbatch; ++ib) {
        for (int i = 0; i < gemm_args.ninter; ++i) {
          if (options.blas_args.batch_size_last_dim) {
            C = Kokkos::subview(gemm_args.C, Kokkos::ALL(), Kokkos::ALL(), ib * gemm_args.ninter + i);
          } else {
            C = Kokkos::subview(gemm_args.C, ib * gemm_args.ninter + i, Kokkos::ALL(), Kokkos::ALL());
          }
          auto info = armpl_dge_deinterleave(gemm_args.ninter, i, gemm_args.dims.c.m, gemm_args.dims.c.n, C.data(),
                                             C.stride(0), C.stride(1), &gemm_args.C_pl.mat[gemm_args.C_pl.bstrd * ib],
                                             gemm_args.C_pl.istrd, gemm_args.C_pl.jstrd);
          if (info != ARMPL_STATUS_SUCCESS) {
            FATAL_ERROR("armpl_dge_deinterleave (C)\n");
          }
        }
      }
    }
#endif  // KOKKOSKERNELS_ENABLE_TPL_ARMPL && ARMPL_BUILD >= 1058
    if (__gemm_do_compare<ScalarType, LayoutType>(C_expected, gemm_args.C)) FATAL_ERROR("Result value mismatch!");
  }

  if (gemm_args.Cv.vec_3d.data() != nullptr) {
    if (__gemm_do_compare<ScalarType, LayoutType>(C_expected, gemm_args.Cv, options))
      FATAL_ERROR("Result value mismatch!");
  }

  // Run actual timed test.
  options.verify    = false;  // Set verify to false for csv output.
  options.warm_up_n = warm_up_n_bak;
  options.n         = n_bak;
  fn(options, gemm_args);

  // Reset verify for next matrix size.
  options.verify = true;
}

/*************************** Internal setup fns **************************/
template <class scalar_type, class vta, class vtb, class vtc, class device_type>
gemm_args_t __do_setup(options_t options, matrix_dims_t dims) {
  using execution_space = typename device_type::execution_space;

  gemm_args_t gemm_args;
  uint64_t seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(seed);
  STATUS;

  gemm_args.A_pl.mat = nullptr;
  gemm_args.B_pl.mat = nullptr;
  gemm_args.C_pl.mat = nullptr;
  gemm_args.dims     = dims;
  gemm_args.transA   = options.blas_args.gemm.gemm_args.c_str()[0];
  gemm_args.transB   = options.blas_args.gemm.gemm_args.c_str()[1];
  if (options.use_simd) {
    // Calculate the batch size for simd views
    auto a_simd_batch_size = dims.a.k / simd_vector_size + (dims.a.k % simd_vector_size > 0);
    auto b_simd_batch_size = dims.b.k / simd_vector_size + (dims.b.k % simd_vector_size > 0);
    auto c_simd_batch_size = dims.c.k / simd_vector_size + (dims.c.k % simd_vector_size > 0);

    // Reference gemm simd arguments for allocating A, B, and C matrices
    gemm_simd_args_t &A = gemm_args.Av, &B = gemm_args.Bv, &C = gemm_args.Cv;

    if (options.blas_args.batch_size_last_dim) {
      // Construct simd matrices with batch_size in the last dimension (better
      // for LayoutLeft views)
      A.vec_3d  = vector_view_type_3d("A_vector", dims.a.m, dims.a.n, a_simd_batch_size);
      A.mat_4d  = view_type_4d((scalar_type *)A.vec_3d.data(), simd_vector_size, dims.a.m, dims.a.n, a_simd_batch_size);
      A.ivec_4d = internal_vector_view_type_4d((internal_vector_type *)A.mat_4d.data(),
                                               simd_vector_size / simd_internal_vector_size, dims.a.m, dims.a.n,
                                               a_simd_batch_size);

      B.vec_3d  = vector_view_type_3d("B_vector", dims.b.m, dims.b.n, b_simd_batch_size);
      B.mat_4d  = view_type_4d((scalar_type *)B.vec_3d.data(), simd_vector_size, dims.b.m, dims.b.n, b_simd_batch_size);
      B.ivec_4d = internal_vector_view_type_4d((internal_vector_type *)B.mat_4d.data(),
                                               simd_vector_size / simd_internal_vector_size, dims.b.m, dims.b.n,
                                               b_simd_batch_size);

      C.vec_3d  = vector_view_type_3d("C_vector", dims.c.m, dims.c.n, c_simd_batch_size);
      C.mat_4d  = view_type_4d((scalar_type *)C.vec_3d.data(), simd_vector_size, dims.c.m, dims.c.n, c_simd_batch_size);
      C.ivec_4d = internal_vector_view_type_4d((internal_vector_type *)C.mat_4d.data(),
                                               simd_vector_size / simd_internal_vector_size, dims.c.m, dims.c.n,
                                               c_simd_batch_size);

    } else {
      // Construct simd matrices with batch_size in the first dimension (better
      // for LayoutRight views)
      A.vec_3d  = vector_view_type_3d("A_vector", a_simd_batch_size, dims.a.m, dims.a.n);
      A.mat_4d  = view_type_4d((scalar_type *)A.vec_3d.data(), a_simd_batch_size, dims.a.m, dims.a.n, simd_vector_size);
      A.ivec_4d = internal_vector_view_type_4d((internal_vector_type *)A.mat_4d.data(), a_simd_batch_size, dims.a.m,
                                               dims.a.n, simd_vector_size / simd_internal_vector_size);

      B.vec_3d  = vector_view_type_3d("B_vector", b_simd_batch_size, dims.b.m, dims.b.n);
      B.mat_4d  = view_type_4d((scalar_type *)B.vec_3d.data(), b_simd_batch_size, dims.b.m, dims.b.n, simd_vector_size);
      B.ivec_4d = internal_vector_view_type_4d((internal_vector_type *)B.mat_4d.data(), b_simd_batch_size, dims.b.m,
                                               dims.b.n, simd_vector_size / simd_internal_vector_size);

      C.vec_3d  = vector_view_type_3d("C_vector", c_simd_batch_size, dims.c.m, dims.c.n);
      C.mat_4d  = view_type_4d((scalar_type *)C.vec_3d.data(), c_simd_batch_size, dims.c.m, dims.c.n, simd_vector_size);
      C.ivec_4d = internal_vector_view_type_4d((internal_vector_type *)C.mat_4d.data(), c_simd_batch_size, dims.c.m,
                                               dims.c.n, simd_vector_size / simd_internal_vector_size);
    }

    // Use the non-simd 4-rank view type to randomly populate the gemm simd
    // arguments
    using tmp_view_type_4d = Kokkos::View<double ****, default_layout, default_device>;
    tmp_view_type_4d tmpA("tmpA", gemm_args.Av.mat_4d.extent(0), gemm_args.Av.mat_4d.extent(1),
                          gemm_args.Av.mat_4d.extent(2), gemm_args.Av.mat_4d.extent(3));
    Kokkos::fill_random(tmpA, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, double>::max());
    tmp_view_type_4d tmpB("tmpB", gemm_args.Bv.mat_4d.extent(0), gemm_args.Bv.mat_4d.extent(1),
                          gemm_args.Bv.mat_4d.extent(2), gemm_args.Bv.mat_4d.extent(3));
    Kokkos::fill_random(tmpB, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, double>::max());
    tmp_view_type_4d tmpC("tmpC", gemm_args.Cv.mat_4d.extent(0), gemm_args.Cv.mat_4d.extent(1),
                          gemm_args.Cv.mat_4d.extent(2), gemm_args.Cv.mat_4d.extent(3));
    Kokkos::fill_random(tmpC, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, double>::max());
    Kokkos::fence();
    Kokkos::deep_copy(gemm_args.Av.mat_4d, tmpA);
    Kokkos::deep_copy(gemm_args.Bv.mat_4d, tmpB);
    Kokkos::deep_copy(gemm_args.Cv.mat_4d, tmpC);
    Kokkos::fence();
  } else {
    if (options.blas_args.batch_size_last_dim) {
      gemm_args.A = vta("gemm_args.A", dims.a.m, dims.a.n, dims.a.k);
      gemm_args.B = vtb("gemm_args.B", dims.b.m, dims.b.n, dims.b.k);
      gemm_args.C = vtc("gemm_args.C", dims.c.m, dims.c.n, dims.c.k);
    } else {
      gemm_args.A = vta("gemm_args.A", dims.a.k, dims.a.m, dims.a.n);
      gemm_args.B = vtb("gemm_args.B", dims.b.k, dims.b.m, dims.b.n);
      gemm_args.C = vtc("gemm_args.C", dims.c.k, dims.c.m, dims.c.n);
    }

    using tmp_view_type_3d = Kokkos::View<double ***, default_layout, default_device>;
    tmp_view_type_3d tmpA("tmpA", gemm_args.A.extent(0), gemm_args.A.extent(1), gemm_args.A.extent(2));
    Kokkos::fill_random(tmpA, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, double>::max());
    tmp_view_type_3d tmpB("tmpB", gemm_args.B.extent(0), gemm_args.B.extent(1), gemm_args.B.extent(2));
    Kokkos::fill_random(tmpB, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, double>::max());
    tmp_view_type_3d tmpC("tmpC", gemm_args.C.extent(0), gemm_args.C.extent(1), gemm_args.C.extent(2));
    Kokkos::fill_random(tmpC, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, double>::max());

    Kokkos::fence();
    Kokkos::deep_copy(gemm_args.A, tmpA);
    Kokkos::deep_copy(gemm_args.B, tmpB);
    Kokkos::deep_copy(gemm_args.C, tmpC);
    Kokkos::fence();
  }

#if defined(KOKKOSKERNELS_ENABLE_TPL_ARMPL) && ARMPL_BUILD >= 1058
  if (options.test == EXPERIMENT) {
    armpl_int_t bstrd_A, istrd_A, jstrd_A, bstrd_B, istrd_B, jstrd_B, bstrd_C, istrd_C, jstrd_C;

    armpl_int_t ninter = (armpl_int_t)options.ninter, nbatch = gemm_args.dims.c.k / ninter;

    if (gemm_args.dims.c.k % ninter) FATAL_ERROR("batch size must be evenly divisible by ninter!");

    jstrd_A = ninter;
    istrd_A = jstrd_A * gemm_args.dims.a.n;
    bstrd_A = istrd_A * gemm_args.dims.a.m;

    jstrd_B = ninter;
    istrd_B = jstrd_B * gemm_args.dims.b.n;
    bstrd_B = istrd_B * gemm_args.dims.b.m;

    jstrd_C = ninter;
    istrd_C = jstrd_C * gemm_args.dims.c.n;
    bstrd_C = istrd_C * gemm_args.dims.c.m;

    gemm_args.ninter     = ninter;
    gemm_args.nbatch     = nbatch;
    gemm_args.A_pl.jstrd = jstrd_A;
    gemm_args.B_pl.jstrd = jstrd_B;
    gemm_args.C_pl.jstrd = jstrd_C;
    gemm_args.A_pl.istrd = istrd_A;
    gemm_args.B_pl.istrd = istrd_B;
    gemm_args.C_pl.istrd = istrd_C;
    gemm_args.A_pl.bstrd = bstrd_A;
    gemm_args.B_pl.bstrd = bstrd_B;
    gemm_args.C_pl.bstrd = bstrd_C;

    default_scalar *A_p = (default_scalar *)malloc(sizeof(default_scalar) * bstrd_A * nbatch);
    default_scalar *B_p = (default_scalar *)malloc(sizeof(default_scalar) * bstrd_B * nbatch);
    default_scalar *C_p = (default_scalar *)malloc(sizeof(default_scalar) * bstrd_C * nbatch);

    using view_type_2d = Kokkos::View<default_scalar **, Kokkos::LayoutStride, default_device>;
    view_type_2d A, B, C;

    // Populate interleave-batch matrices
    for (int ib = 0; ib < nbatch; ++ib) {
      for (int i = 0; i < ninter; ++i) {
        if (options.blas_args.batch_size_last_dim) {
          A = Kokkos::subview(gemm_args.A, Kokkos::ALL(), Kokkos::ALL(), ib * ninter + i);
          B = Kokkos::subview(gemm_args.B, Kokkos::ALL(), Kokkos::ALL(), ib * ninter + i);
          C = Kokkos::subview(gemm_args.C, Kokkos::ALL(), Kokkos::ALL(), ib * ninter + i);
        } else {
          A = Kokkos::subview(gemm_args.A, ib * ninter + i, Kokkos::ALL(), Kokkos::ALL());
          B = Kokkos::subview(gemm_args.B, ib * ninter + i, Kokkos::ALL(), Kokkos::ALL());
          C = Kokkos::subview(gemm_args.C, ib * ninter + i, Kokkos::ALL(), Kokkos::ALL());
        }

        auto info = armpl_dge_interleave(ninter, i, gemm_args.dims.a.m, gemm_args.dims.a.n, A.data(), A.stride(0),
                                         A.stride(1), &A_p[bstrd_A * ib], istrd_A, jstrd_A);
        if (info != ARMPL_STATUS_SUCCESS) {
          FATAL_ERROR("armpl_dge_interleave (A)\n");
        }
        info = armpl_dge_interleave(ninter, i, gemm_args.dims.b.m, gemm_args.dims.b.n, B.data(), B.stride(0),
                                    B.stride(1), &B_p[bstrd_B * ib], istrd_B, jstrd_B);
        if (info != ARMPL_STATUS_SUCCESS) {
          FATAL_ERROR("armpl_dge_interleave (B)\n");
        }
        info = armpl_dge_interleave(ninter, i, gemm_args.dims.c.m, gemm_args.dims.c.n, C.data(), C.stride(0),
                                    C.stride(1), &C_p[bstrd_C * ib], istrd_C, jstrd_C);
        if (info != ARMPL_STATUS_SUCCESS) {
          FATAL_ERROR("armpl_dge_interleave (C)\n");
        }
      }
    }

    gemm_args.A_pl.mat = A_p;
    gemm_args.B_pl.mat = B_p;
    gemm_args.C_pl.mat = C_p;
  }
#endif  // KOKKOSKERNELS_ENABLE_TPL_ARMPL && ARMPL_BUILD >= 1058

  gemm_args.alpha         = options.blas_args.gemm.alpha;
  gemm_args.beta          = options.blas_args.gemm.beta;
  gemm_args.bp.team_size  = options.blas_args.team_size;
  gemm_args.bp.vector_len = options.blas_args.vector_len;

  Kokkos::fence();  // Ensure that fill_random has completed.

  return gemm_args;
}

/*************************** Interal run helper fns **************************/
void __do_loop_and_invoke(options_t options, void (*fn)(options_t, gemm_args_t)) {
  matrix_dims_t cur_dims;
  gemm_args_t gemm_args;
  STATUS;

  __print_gemm_perf_test_options(options);
  std::cout << "SCALAR:" << typeid(default_scalar).name() << ", LAYOUT:" << typeid(default_layout).name()
            << ", DEVICE:" << typeid(default_device).name() << ", SPACE:" << typeid(memory_space).name() << std::endl;

  options.out[0] << gemm_csv_header_str << std::endl;

  for (cur_dims = options.start;
       cur_dims.a.m <= options.stop.a.m && cur_dims.a.n <= options.stop.a.n && cur_dims.b.m <= options.stop.b.m &&
       cur_dims.b.n <= options.stop.b.n && cur_dims.c.m <= options.stop.c.m && cur_dims.c.n <= options.stop.c.n;
       cur_dims.a.m += options.step, cur_dims.a.n += options.step, cur_dims.b.m += options.step,
      cur_dims.b.n += options.step, cur_dims.c.m += options.step, cur_dims.c.n += options.step) {
    gemm_args = __do_setup<default_scalar, view_type_3d, view_type_3d, view_type_3d, default_device>(options, cur_dims);

    if (options.verify) {
      __gemm_do_verify<default_scalar, default_layout, default_device>(options, gemm_args, fn);
    } else {
      fn(options, gemm_args);
    }

    if (gemm_args.A_pl.mat != nullptr) {
      free(gemm_args.A_pl.mat);
      free(gemm_args.B_pl.mat);
      free(gemm_args.C_pl.mat);
    }
  }
  return;
}

/*************************** External fns **************************/
void do_gemm_serial_blas(options_t options) {
  STATUS;
  __do_loop_and_invoke(options, __do_gemm_serial_blas<default_scalar, view_type_3d, view_type_3d, default_device>);
  return;
}

void do_gemm_serial_batched(options_t options) {
  STATUS;
  __do_loop_and_invoke(options, __do_gemm_serial_batched<default_scalar, view_type_3d, view_type_3d, view_type_3d,
                                                         default_device, KokkosBatched::Algo::Gemm::Unblocked>);
  return;
}

void do_gemm_serial_batched_blocked(options_t options) {
  STATUS;
  __do_loop_and_invoke(options, __do_gemm_serial_batched<default_scalar, view_type_3d, view_type_3d, view_type_3d,
                                                         default_device, KokkosBatched::Algo::Gemm::Blocked>);
  return;
}

void do_gemm_heuristic_batched_parallel(options_t options) {
  STATUS;
  if (options.blas_args.use_auto) {
    fprintf(stderr, "ERROR: --test=%s does not support --use_auto=%d\n", test_e_str[options.test].c_str(),
            (int)options.blas_args.use_auto);
    exit(-EINVAL);
  }

  __do_loop_and_invoke(options, __do_gemm_parallel_batched_heuristic<void, void, default_device>);
  return;
}

void do_gemm_serial_batched_parallel(options_t options) {
  STATUS;
  if (options.blas_args.batch_size_last_dim)
    __do_loop_and_invoke(
        options, __do_gemm_parallel_batched<SerialBatchDim3Tag, KokkosBatched::Algo::Gemm::Unblocked, default_device>);
  else
    __do_loop_and_invoke(options,
                         __do_gemm_parallel_batched<SerialTag, KokkosBatched::Algo::Gemm::Unblocked, default_device>);
  return;
}

void do_gemm_serial_batched_blocked_parallel(options_t options) {
  STATUS;
  if (options.blas_args.batch_size_last_dim)
    __do_loop_and_invoke(
        options, __do_gemm_parallel_batched<SerialBatchDim3Tag, KokkosBatched::Algo::Gemm::Blocked, default_device>);
  else
    __do_loop_and_invoke(options,
                         __do_gemm_parallel_batched<SerialTag, KokkosBatched::Algo::Gemm::Blocked, default_device>);
  return;
}

void do_gemm_serial_simd_batched_parallel(options_t options) {
  STATUS;
  // SerialBatchDim3Tag
  // SerialSimdTag
  options.use_simd = true;
  if (options.blas_args.batch_size_last_dim)
    __do_loop_and_invoke(options, __do_gemm_parallel_batched<TeamSimdBatchDim4Tag, KokkosBatched::Algo::Gemm::Unblocked,
                                                             default_device, KokkosBatched::Mode::Serial>);
  else
    __do_loop_and_invoke(options, __do_gemm_parallel_batched<TeamSimdTag, KokkosBatched::Algo::Gemm::Unblocked,
                                                             default_device, KokkosBatched::Mode::Serial>);
  return;
}

void do_gemm_serial_simd_batched_blocked_parallel(options_t options) {
  STATUS;
  // SerialBatchDim3Tag
  // SerialSimdTag
  options.use_simd = true;
  if (options.blas_args.batch_size_last_dim)
    __do_loop_and_invoke(options, __do_gemm_parallel_batched<TeamSimdBatchDim4Tag, KokkosBatched::Algo::Gemm::Blocked,
                                                             default_device, KokkosBatched::Mode::Serial>);
  else
    __do_loop_and_invoke(options, __do_gemm_parallel_batched<TeamSimdTag, KokkosBatched::Algo::Gemm::Blocked,
                                                             default_device, KokkosBatched::Mode::Serial>);
  return;
}

#if defined(__KOKKOSBATCHED_ENABLE_INTEL_MKL__) && defined(__KOKKOSBATCHED_ENABLE_INTEL_MKL_BATCHED__) && \
    defined(__KOKKOSBATCHED_ENABLE_INTEL_MKL_COMPACT_BATCHED__)
void do_gemm_serial_batched_compact_mkl_parallel(options_t options) {
  STATUS;
  if (options.blas_args.batch_size_last_dim)
    __do_loop_and_invoke(
        options,
        __do_gemm_parallel_batched<SerialSimdBatchDim3Tag, KokkosBatched::Algo::Gemm::CompactMKL, default_device>);
  else
    __do_loop_and_invoke(
        options, __do_gemm_parallel_batched<SerialSimdTag, KokkosBatched::Algo::Gemm::CompactMKL, default_device>);
  return;
}
#else
void do_gemm_serial_batched_compact_mkl_parallel(options_t) {
  STATUS;
#if !defined(__KOKKOSBATCHED_ENABLE_INTEL_MKL__)
  std::cerr << std::string(__func__) << " disabled since __KOKKOSBATCHED_ENABLE_INTEL_MKL__ is undefined." << std::endl;
#elif !defined(__KOKKOSBATCHED_ENABLE_INTEL_MKL_BATCHED__)
  std::cerr << std::string(__func__)
            << " disabled since __KOKKOSBATCHED_ENABLE_INTEL_MKL_BATCHED__ is "
               "undefined."
            << std::endl;
#elif !defined(__KOKKOSBATCHED_ENABLE_INTEL_MKL_COMPACT_BATCHED__)
  std::cerr << std::string(__func__)
            << " disabled since __KOKKOSBATCHED_ENABLE_INTEL_MKL_COMPACT_BATCHED__ "
               "is undefined."
            << std::endl;
#endif
  return;
}
#endif

void do_gemm_team_batched_parallel(options_t options) {
  STATUS;
  if (options.blas_args.batch_size_last_dim)
    __do_loop_and_invoke(
        options, __do_gemm_parallel_batched<TeamBatchDim3Tag, KokkosBatched::Algo::Gemm::Unblocked, default_device>);
  else
    __do_loop_and_invoke(options,
                         __do_gemm_parallel_batched<TeamTag, KokkosBatched::Algo::Gemm::Unblocked, default_device>);
  return;
}

void do_gemm_team_batched_blocked_parallel(options_t options) {
  STATUS;
  if (options.blas_args.batch_size_last_dim)
    __do_loop_and_invoke(
        options, __do_gemm_parallel_batched<TeamBatchDim3Tag, KokkosBatched::Algo::Gemm::Blocked, default_device>);
  else
    __do_loop_and_invoke(options,
                         __do_gemm_parallel_batched<TeamTag, KokkosBatched::Algo::Gemm::Blocked, default_device>);
  return;
}

void do_gemm_team_vector_batched_parallel(options_t options) {
  STATUS;
  if (options.blas_args.batch_size_last_dim)
    __do_loop_and_invoke(
        options,
        __do_gemm_parallel_batched<TeamVectorBatchDim3Tag, KokkosBatched::Algo::Gemm::Unblocked, default_device>);
  else
    __do_loop_and_invoke(
        options, __do_gemm_parallel_batched<TeamVectorTag, KokkosBatched::Algo::Gemm::Unblocked, default_device>);
  return;
}

void do_gemm_team_simd_batched_parallel(options_t options) {
  STATUS;
  options.use_simd = true;
  if (options.blas_args.batch_size_last_dim)
    __do_loop_and_invoke(options, __do_gemm_parallel_batched<TeamSimdBatchDim4Tag, KokkosBatched::Algo::Gemm::Unblocked,
                                                             default_device, KokkosBatched::Mode::Team>);
  else
    __do_loop_and_invoke(options, __do_gemm_parallel_batched<TeamSimdTag, KokkosBatched::Algo::Gemm::Unblocked,
                                                             default_device, KokkosBatched::Mode::Team>);
  return;
}

void do_gemm_team_simd_batched_blocked_parallel(options_t options) {
  STATUS;
  options.use_simd = true;
  if (options.blas_args.batch_size_last_dim)
    __do_loop_and_invoke(options, __do_gemm_parallel_batched<TeamSimdBatchDim4Tag, KokkosBatched::Algo::Gemm::Blocked,
                                                             default_device, KokkosBatched::Mode::Team>);
  else
    __do_loop_and_invoke(options, __do_gemm_parallel_batched<TeamSimdTag, KokkosBatched::Algo::Gemm::Blocked,
                                                             default_device, KokkosBatched::Mode::Team>);
  return;
}

// Blocked algo not yet implemented for TeamVectorGemm.
/* void do_gemm_team_vector_batched_blocked_parallel(options_t options) {
  STATUS;
  __do_loop_and_invoke(
      options, __do_gemm_parallel_batched<TeamVectorTag,
KokkosBatched::Algo::Gemm::Blocked, default_device>); return;
} */

void do_gemm_experiment_parallel(options_t options) {
  STATUS;
  using TransAType   = KokkosBatched::Trans::NoTranspose;
  using TransBType   = KokkosBatched::Trans::NoTranspose;
  using BlockingType = KokkosBatched::Algo::Gemm::Unblocked;

  // __do_loop_and_invoke(
  //     options, __do_gemm_parallel_experiment1<TransAType, TransBType,
  //                                             BlockingType, default_device>);
  // __do_loop_and_invoke(
  //     options, __do_gemm_parallel_experiment2<TransAType, TransBType,
  //                                             BlockingType, default_device>);
  // __do_loop_and_invoke(
  //     options, __do_gemm_parallel_experiment3<TransAType, TransBType,
  //                                             BlockingType, default_device>);
  // __do_loop_and_invoke(
  //     options, __do_gemm_parallel_experiment4<TransAType, TransBType,
  //                                             BlockingType, default_device>);
  // __do_loop_and_invoke(
  //     options, __do_gemm_parallel_experiment5<TransAType, TransBType,
  //                                             BlockingType, default_device>);
  // __do_loop_and_invoke(
  //     options, __do_gemm_parallel_experiment6<TransAType, TransBType,
  //                                             BlockingType, default_device>);
  __do_loop_and_invoke(options, __do_gemm_armpl<TransAType, TransBType, BlockingType, default_device>);
}

#endif  // KOKKOSBLAS3_GEMM_PERF_TEST_H_
