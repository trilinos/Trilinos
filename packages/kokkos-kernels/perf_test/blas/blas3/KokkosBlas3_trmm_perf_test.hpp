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
#ifndef KOKKOSBLAS3_TRMM_PERF_TEST_H_
#define KOKKOSBLAS3_TRMM_PERF_TEST_H_

// #include <complex.h>
#include "KokkosBlas3_common.hpp"

#include <Kokkos_Random.hpp>

#include <KokkosBlas3_trmm.hpp>

#include "KokkosBatched_Trmm_Decl.hpp"
#include "KokkosBatched_Trmm_Serial_Impl.hpp"
#include "KokkosBatched_Util.hpp"

#include <chrono>

// #define PERF_TEST_DEBUG

// Forward declarations
void do_trmm_serial_blas(options_t options);
void do_trmm_serial_batched(options_t options);
void do_trmm_parallel_blas(options_t options);
void do_trmm_parallel_batched(options_t options);

// trmm invoke table
void (*do_trmm_invoke[LOOP_N][TEST_N])(options_t) = {{do_trmm_serial_blas, do_trmm_serial_batched},
                                                     {do_trmm_parallel_blas, do_trmm_parallel_batched}};

/*************************** Test types and defaults **************************/
#define DEFAULT_TRMM_ARGS "LUNU"
#define DEFAULT_TRMM_ALPHA 1.0

/**
 * The KokkosBatched::SerialTrmm implementation performs dot products on
 * non-zero elements of the triangular matrices. The flop calculation below
 * assumes KokkosBatched::SerialTrmm is being used. Since the dot products
 * do a multiply and add we can calculate the flops for any element in the last
 * column of the LHS to be 2*columns_LHS, any element in the last-1 column of
 * the LHS to be 2*(columns_LHS-1), and so on. We do this for every row of the
 * LHS giving us this flop count: flops = columns_LHS * (columns_LHS + 1) flops
 * = (flops / 2) * 2 flops = flops * rows_LHS
 */
static inline int __trmm_impl_flop_count(char side, int b_m, int b_n, int /*a_m*/, int /*a_n*/) {
  int flops;

  if (side == 'L' || side == 'l') {
    flops = (b_m * (b_m + 1)) * b_n;
  } else {
    flops = (b_n * (b_n + 1)) * b_m;
  }

  if (std::is_same<double, default_scalar>::value || std::is_same<float, default_scalar>::value ||
      std::is_same<Kokkos::Experimental::half_t, default_scalar>::value)
    return flops;

  // Account for 6 additional flops when complex numbers are used.
  // Above we have counted 1 flop for each add and 1 flop for each multiply.
  // For complex, we need to count 2 flops for each add and 6 flops for each
  // multiply.
  return flops * 4;
}

// Flop count formula from lapack working note 41:
// http://www.icl.utk.edu/~mgates3/docs/lawn41.pdf
static inline double __trmm_flop_count(char side, double b_m, double b_n, double /*a_m*/, double /*a_n*/) {
  double flops;

  if (side == 'L' || side == 'l') {
    flops = b_m * b_m * b_n;
  } else {
    flops = b_n * b_n * b_m;
  }

  if (std::is_same<double, default_scalar>::value || std::is_same<float, default_scalar>::value ||
      std::is_same<Kokkos::Experimental::half_t, default_scalar>::value)
    return flops;

  // Account for 6 additional flops when complex numbers are used.
  // Above we have counted 1 flop for each add and 1 flop for each multiply.
  // For complex, we need to count 2 flops for each add and 6 flops for each
  // multiply.
  return flops * 4;
}

using view_type_3d = Kokkos::View<default_scalar***, default_layout, default_device>;
struct trmm_args {
  char side, uplo, trans, diag;
  default_scalar alpha;
  view_type_3d A, B;
};
typedef struct trmm_args trmm_args_t;

static std::string trmm_csv_header_str =
    "algorithm,side-uplo-trans-diag,alpha,loop_type,A_dims,B_dims,warm_up_n,"
    "iter,total_time(s),average_time(s),FLOPS,GFLOP/"
    "average_time(s),min_achieved_bandwidth(GB/s),max_achieved_bandwidth(GB/s)";

/*************************** Internal helper fns **************************/
static void __trmm_output_csv_row(options_t options, trmm_args_t trmm_args, double time_in_seconds) {
  double flops = trmm_args.A.extent(0) * __trmm_flop_count(trmm_args.side, trmm_args.B.extent(1), trmm_args.B.extent(2),
                                                           trmm_args.A.extent(1), trmm_args.A.extent(2));
  double gflops       = flops / 1e9;
  double average_time = time_in_seconds / options.n;
  double gbytes_in_matrix =
      (trmm_args.B.extent(0) * trmm_args.B.extent(1) * trmm_args.B.extent(2) * sizeof(default_scalar)) / 1e9;
  double min_memory_transactions, max_memory_transactions;

  // Assuming infinite cache size
  // We have to read A and B into the cache once and then write
  // B back out to main memory once.
  min_memory_transactions = 3;

  // Assuming no register or real caching
  // We have to go out to memory for every element we read from A and B as well
  // as every element we write to B. We use the trmm flops from lapack note 41
  // and multiple by 3/2 to account for the write to B since this flop count is
  // for one multiply and one add.
  if (trmm_args.side == 'l' || trmm_args.side == 'L')
    max_memory_transactions = trmm_args.B.extent(1) * trmm_args.B.extent(1) * trmm_args.B.extent(2) * (3. / 2.);
  else
    max_memory_transactions = trmm_args.B.extent(2) * trmm_args.B.extent(2) * trmm_args.B.extent(1) * (3. / 2.);

  options.out[0] << test_e_str[options.test] << "," << options.blas_args.trmm.trmm_args << ","
                 << static_cast<double>(options.blas_args.trmm.alpha) << "," << loop_e_str[options.loop] << ","
                 << trmm_args.A.extent(0) << "x" << trmm_args.A.extent(1) << "x" << trmm_args.A.extent(2) << ","
                 << trmm_args.B.extent(0) << "x" << trmm_args.B.extent(1) << "x" << trmm_args.B.extent(2) << ","
                 << options.warm_up_n << "," << options.n << "," << time_in_seconds << "," << average_time << ","
                 << flops << "," << gflops / average_time << ","
                 << (gbytes_in_matrix * min_memory_transactions) / average_time << ","
                 << (gbytes_in_matrix * max_memory_transactions) / average_time << std::endl;
}

#ifdef PERF_TEST_DEBUG
static void __print_trmm_perf_test_options(options_t options) {
  printf("options.test      = %s\n", test_e_str[options.test].c_str());
  printf("options.loop      = %s\n", loop_e_str[options.loop].c_str());
  printf("options.start     = %dx%d,%dx%d\n", options.start.a.m, options.start.a.n, options.start.b.m,
         options.start.b.n);
  printf("options.stop      = %dx%d,%dx%d\n", options.stop.a.m, options.stop.a.n, options.stop.b.m, options.stop.b.n);
  printf("options.step      = %d\n", options.step);
  printf("options.warm_up_n = %d\n", options.warm_up_n);
  printf("options.n         = %d\n", options.n);
  printf("options.blas_args.trmm.trmm_args = %s\n", options.blas_args.trmm.trmm_args.c_str());
  printf("options.out_file  = %s\n", options.out_file.c_str());
  if (std::is_same<double, default_scalar>::value)
    printf("options.alpha     = %lf\n", options.blas_args.trmm.alpha);
  else if (std::is_same<float, default_scalar>::value)
    printf("options.alpha     = %f\n", options.blas_args.trmm.alpha);
  return;
}
#else
static void __print_trmm_perf_test_options(options_t /*options*/) { return; }
#endif  // PERF_TEST_DEBUG

/*************************** Internal templated fns **************************/
// Need to take subviews on the device
#if !defined(KOKKOS_ENABLE_CUDA) && !defined(KOKKOS_ENABLE_HIP) && !defined(KOKKOS_ENABLE_OPENMPTARGET)
template <class scalar_type, class vta, class vtb, class device_type>
void __do_trmm_serial_blas(options_t options, trmm_args_t trmm_args) {
  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  Kokkos::Timer timer;

  STATUS;

  for (uint32_t j = 0; j < warm_up_n; ++j) {
    for (int i = 0; i < options.start.a.k; ++i) {
      auto A = Kokkos::subview(trmm_args.A, i, Kokkos::ALL(), Kokkos::ALL());
      auto B = Kokkos::subview(trmm_args.B, i, Kokkos::ALL(), Kokkos::ALL());

      KokkosBlas::trmm(&trmm_args.side, &trmm_args.uplo, &trmm_args.trans, &trmm_args.diag, trmm_args.alpha, A, B);
    }
    // Fence after submitting each batch operation
    Kokkos::fence();
  }

  timer.reset();
  for (uint32_t j = 0; j < n; ++j) {
    for (int i = 0; i < options.start.a.k; ++i) {
      auto A = Kokkos::subview(trmm_args.A, i, Kokkos::ALL(), Kokkos::ALL());
      auto B = Kokkos::subview(trmm_args.B, i, Kokkos::ALL(), Kokkos::ALL());

      KokkosBlas::trmm(&trmm_args.side, &trmm_args.uplo, &trmm_args.trans, &trmm_args.diag, trmm_args.alpha, A, B);
    }
    // Fence after submitting each batch operation
    Kokkos::fence();
  }
  __trmm_output_csv_row(options, trmm_args, timer.seconds());
  return;
}
#else
template <class scalar_type, class vta, class vtb, class device_type>
void __do_trmm_serial_blas(options_t /*options*/, trmm_args_t /*trmm_args*/) {
  std::cerr << std::string(__func__)
            << " disabled since KOKKOS_ENABLE_CUDA or "
               "KOKKOS_ENABLE_OPENMPTARGET is defined."
            << std::endl;
  return;
}
#endif  // !KOKKOS_ENABLE_CUDA && !KOKKOS_ENABLE_OPENMPTARGET

// Need to take subviews on the device
#if !defined(KOKKOS_ENABLE_CUDA) && !defined(KOKKOS_ENABLE_HIP) && !defined(KOKKOS_ENABLE_OPENMPTARGET)
template <class side, class uplo, class trans, class diag>
void __do_trmm_serial_batched_template(options_t options, trmm_args_t trmm_args) {
  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  Kokkos::Timer timer;
  using tag = KokkosBatched::Algo::Trmm::Unblocked;

  for (uint32_t j = 0; j < warm_up_n; ++j) {
    for (int i = 0; i < options.start.a.k; ++i) {
      auto A = Kokkos::subview(trmm_args.A, i, Kokkos::ALL(), Kokkos::ALL());
      auto B = Kokkos::subview(trmm_args.B, i, Kokkos::ALL(), Kokkos::ALL());

      KokkosBatched::SerialTrmm<side, uplo, trans, diag, tag>::invoke(trmm_args.alpha, A, B);
    }
    // Fence after submitting each batch operation
    Kokkos::fence();
  }

  timer.reset();
  for (uint32_t j = 0; j < n; ++j) {
    for (int i = 0; i < options.start.a.k; ++i) {
      auto A = Kokkos::subview(trmm_args.A, i, Kokkos::ALL(), Kokkos::ALL());
      auto B = Kokkos::subview(trmm_args.B, i, Kokkos::ALL(), Kokkos::ALL());

      KokkosBatched::SerialTrmm<side, uplo, trans, diag, tag>::invoke(trmm_args.alpha, A, B);
    }
    // Fence after submitting each batch operation
    Kokkos::fence();
  }
  __trmm_output_csv_row(options, trmm_args, timer.seconds());
}
#else
template <class side, class uplo, class trans, class diag>
void __do_trmm_serial_batched_template(options_t /*options*/, trmm_args_t /*trmm_args*/) {
  std::cerr << std::string(__func__)
            << " disabled since KOKKOS_ENABLE_CUDA or "
               "KOKKOS_ENABLE_OPENMPTARGET is defined."
            << std::endl;
}
#endif  // !KOKKOS_ENABLE_CUDA && !KOKKOS_ENABLE_OPENMPTARGET

template <class scalar_type, class vta, class vtb, class device_type>
void __do_trmm_serial_batched(options_t options, trmm_args_t trmm_args) {
  char __side = tolower(trmm_args.side), __uplo = tolower(trmm_args.uplo), __trans = tolower(trmm_args.trans);
  //__diag = tolower(diag[0]);

  using KokkosBatched::Diag;
  using KokkosBatched::Side;
  using KokkosBatched::Trans;
  using KokkosBatched::Uplo;

  STATUS;

  //// Lower non-transpose ////
  if (__side == 'l' && __uplo == 'l' && __trans == 'n') {
    STATUS;
    __do_trmm_serial_batched_template<Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::Unit>(options, trmm_args);
  }
  if (__side == 'r' && __uplo == 'l' && __trans == 'n') {
    STATUS;
    __do_trmm_serial_batched_template<Side::Right, Uplo::Lower, Trans::NoTranspose, Diag::Unit>(options, trmm_args);
  }
  //// Lower transpose /////
  // Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (__side == 'l' && __uplo == 'l' && __trans == 't') {
    STATUS;
    __do_trmm_serial_batched_template<Side::Left, Uplo::Lower, Trans::Transpose, Diag::Unit>(options, trmm_args);
  }
  if (__side == 'r' && __uplo == 'l' && __trans == 't') {
    STATUS;
    __do_trmm_serial_batched_template<Side::Right, Uplo::Lower, Trans::Transpose, Diag::Unit>(options, trmm_args);
  }
  //// Lower conjugate-transpose ////
  // Conjugate-Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (__side == 'l' && __uplo == 'l' && __trans == 'c') {
    STATUS;
    __do_trmm_serial_batched_template<Side::Left, Uplo::Lower, Trans::ConjTranspose, Diag::Unit>(options, trmm_args);
  }
  if (__side == 'r' && __uplo == 'l' && __trans == 'c') {
    STATUS;
    __do_trmm_serial_batched_template<Side::Right, Uplo::Lower, Trans::ConjTranspose, Diag::Unit>(options, trmm_args);
  }
  //// Upper non-transpose ////
  if (__side == 'l' && __uplo == 'u' && __trans == 'n') {
    STATUS;
    __do_trmm_serial_batched_template<Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::Unit>(options, trmm_args);
  }
  if (__side == 'r' && __uplo == 'u' && __trans == 'n') {
    STATUS;
    __do_trmm_serial_batched_template<Side::Right, Uplo::Upper, Trans::NoTranspose, Diag::Unit>(options, trmm_args);
  }
  //// Upper transpose
  // Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (__side == 'l' && __uplo == 'u' && __trans == 't') {
    STATUS;
    __do_trmm_serial_batched_template<Side::Left, Uplo::Upper, Trans::Transpose, Diag::Unit>(options, trmm_args);
  }
  if (__side == 'r' && __uplo == 'u' && __trans == 't') {
    STATUS;
    __do_trmm_serial_batched_template<Side::Right, Uplo::Upper, Trans::Transpose, Diag::Unit>(options, trmm_args);
  }

  //// Upper conjugate-transpose ////
  // Conjugate-Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (__side == 'l' && __uplo == 'u' && __trans == 'c') {
    STATUS;
    __do_trmm_serial_batched_template<Side::Left, Uplo::Upper, Trans::ConjTranspose, Diag::Unit>(options, trmm_args);
  }
  if (__side == 'r' && __uplo == 'u' && __trans == 'c') {
    STATUS;
    __do_trmm_serial_batched_template<Side::Right, Uplo::Upper, Trans::ConjTranspose, Diag::Unit>(options, trmm_args);
  }

  return;
}

#if !defined(KOKKOS_ENABLE_CUDA) && !defined(KOKKOS_ENABLE_HIP) && !defined(KOKKOS_ENABLE_OPENMPTARGET)
template <class ExecutionSpace>
struct parallel_blas_trmm {
  trmm_args_t trmm_args_;

  parallel_blas_trmm(trmm_args_t trmm_args) : trmm_args_(trmm_args) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i) const {
    auto svA = Kokkos::subview(trmm_args_.A, i, Kokkos::ALL(), Kokkos::ALL());
    auto svB = Kokkos::subview(trmm_args_.B, i, Kokkos::ALL(), Kokkos::ALL());

    KokkosBlas::trmm(&trmm_args_.side, &trmm_args_.uplo, &trmm_args_.trans, &trmm_args_.diag, trmm_args_.alpha, svA,
                     svB);
  }
};
#endif  // !KOKKOSKERNELS_ENABLE_DEVICE

template <class scalar_type, class vta, class vtb, class device_type>
void __do_trmm_parallel_blas(options_t options, trmm_args_t trmm_args) {
// TODO: Note why this is disabled on CUDA, OPENMPTARGET and HIP
#if !defined(KOKKOS_ENABLE_CUDA) && !defined(KOKKOS_ENABLE_HIP) && !defined(KOKKOS_ENABLE_OPENMPTARGET)
  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  Kokkos::Timer timer;
  using execution_space = typename device_type::execution_space;
  using functor_type    = parallel_blas_trmm<execution_space>;
  functor_type parallel_blas_trmm_functor(trmm_args);

  STATUS;

  for (uint32_t j = 0; j < warm_up_n; ++j) {
    Kokkos::parallel_for("parallelBlasWarmUpLoopTrmm", Kokkos::RangePolicy<execution_space>(0, options.start.a.k),
                         parallel_blas_trmm_functor);
    // Fence after each batch operation
    Kokkos::fence();
  }

  timer.reset();
  for (uint32_t j = 0; j < n; ++j) {
    Kokkos::parallel_for("parallelBlasTimedLoopTrmm", Kokkos::RangePolicy<execution_space>(0, options.start.a.k),
                         parallel_blas_trmm_functor);
    // Fence after each batch operation
    Kokkos::fence();
  }
  __trmm_output_csv_row(options, trmm_args, timer.seconds());
#else
  std::cerr << std::string(__func__)
            << " disabled since KOKKOS_ENABLE_CUDA, KOKKOS_ENABLE_HIP "
               "or KOKKOS_ENABLE_OPENMPTARGET is defined."
            << std::endl;
  __trmm_output_csv_row(options, trmm_args, -1);
#endif  // !KOKKOS_ENABLE_DEVICE
  return;
}

template <class side, class uplo, class trans, class diag, class tag, class ExecutionSpace>
struct parallel_batched_trmm {
  trmm_args_t trmm_args_;

  parallel_batched_trmm(trmm_args_t trmm_args) : trmm_args_(trmm_args) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i) const {
    auto svA = Kokkos::subview(trmm_args_.A, i, Kokkos::ALL(), Kokkos::ALL());
    auto svB = Kokkos::subview(trmm_args_.B, i, Kokkos::ALL(), Kokkos::ALL());

    KokkosBatched::SerialTrmm<side, uplo, trans, diag, tag>::invoke(trmm_args_.alpha, svA, svB);
  }
};

template <class side, class uplo, class trans, class diag, class device_type>
void __do_trmm_parallel_batched_template(options_t options, trmm_args_t trmm_args) {
  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  Kokkos::Timer timer;
  using tag             = KokkosBatched::Algo::Trmm::Unblocked;
  using execution_space = typename device_type::execution_space;
  using functor_type    = parallel_batched_trmm<side, uplo, trans, diag, tag, execution_space>;
  functor_type parallel_batched_trmm_functor(trmm_args);

  STATUS;

  for (uint32_t j = 0; j < warm_up_n; ++j) {
    Kokkos::parallel_for("parallelBatchedWarmUpLoopTrmm", Kokkos::RangePolicy<execution_space>(0, options.start.a.k),
                         parallel_batched_trmm_functor);
    // Fence after each batch operation
    Kokkos::fence();
  }

  timer.reset();
  for (uint32_t j = 0; j < n; ++j) {
    Kokkos::parallel_for("parallelBatchedTimedLoopTrmm", Kokkos::RangePolicy<execution_space>(0, options.start.a.k),
                         parallel_batched_trmm_functor);
    // Fence after each batch operation
    Kokkos::fence();
  }
  __trmm_output_csv_row(options, trmm_args, timer.seconds());

  return;
}

template <class scalar_type, class vta, class vtb, class device_type>
void __do_trmm_parallel_batched(options_t options, trmm_args_t trmm_args) {
  char __side = tolower(trmm_args.side), __uplo = tolower(trmm_args.uplo), __trans = tolower(trmm_args.trans);
  //__diag = tolower(diag[0]);

  using KokkosBatched::Diag;
  using KokkosBatched::Side;
  using KokkosBatched::Trans;
  using KokkosBatched::Uplo;

  STATUS;

  //// Lower non-transpose ////
  if (__side == 'l' && __uplo == 'l' && __trans == 'n') {
    STATUS;
    __do_trmm_parallel_batched_template<Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::Unit, device_type>(
        options, trmm_args);
  }
  if (__side == 'r' && __uplo == 'l' && __trans == 'n') {
    STATUS;
    __do_trmm_parallel_batched_template<Side::Right, Uplo::Lower, Trans::NoTranspose, Diag::Unit, device_type>(
        options, trmm_args);
  }
  //// Lower transpose /////
  // Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (__side == 'l' && __uplo == 'l' && __trans == 't') {
    STATUS;
    __do_trmm_parallel_batched_template<Side::Left, Uplo::Lower, Trans::Transpose, Diag::Unit, device_type>(options,
                                                                                                            trmm_args);
  }
  if (__side == 'r' && __uplo == 'l' && __trans == 't') {
    STATUS;
    __do_trmm_parallel_batched_template<Side::Right, Uplo::Lower, Trans::Transpose, Diag::Unit, device_type>(options,
                                                                                                             trmm_args);
  }
  //// Lower conjugate-transpose ////
  // Conjugate-Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (__side == 'l' && __uplo == 'l' && __trans == 'c') {
    STATUS;
    __do_trmm_parallel_batched_template<Side::Left, Uplo::Lower, Trans::ConjTranspose, Diag::Unit, device_type>(
        options, trmm_args);
  }
  if (__side == 'r' && __uplo == 'l' && __trans == 'c') {
    STATUS;
    __do_trmm_parallel_batched_template<Side::Right, Uplo::Lower, Trans::ConjTranspose, Diag::Unit, device_type>(
        options, trmm_args);
  }
  //// Upper non-transpose ////
  if (__side == 'l' && __uplo == 'u' && __trans == 'n') {
    STATUS;
    __do_trmm_parallel_batched_template<Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::Unit, device_type>(
        options, trmm_args);
  }
  if (__side == 'r' && __uplo == 'u' && __trans == 'n') {
    STATUS;
    __do_trmm_parallel_batched_template<Side::Right, Uplo::Upper, Trans::NoTranspose, Diag::Unit, device_type>(
        options, trmm_args);
  }
  //// Upper transpose
  // Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (__side == 'l' && __uplo == 'u' && __trans == 't') {
    STATUS;
    __do_trmm_parallel_batched_template<Side::Left, Uplo::Upper, Trans::Transpose, Diag::Unit, device_type>(options,
                                                                                                            trmm_args);
  }
  if (__side == 'r' && __uplo == 'u' && __trans == 't') {
    STATUS;
    __do_trmm_parallel_batched_template<Side::Right, Uplo::Upper, Trans::Transpose, Diag::Unit, device_type>(options,
                                                                                                             trmm_args);
  }

  //// Upper conjugate-transpose ////
  // Conjugate-Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (__side == 'l' && __uplo == 'u' && __trans == 'c') {
    STATUS;
    __do_trmm_parallel_batched_template<Side::Left, Uplo::Upper, Trans::ConjTranspose, Diag::Unit, device_type>(
        options, trmm_args);
  }
  if (__side == 'r' && __uplo == 'u' && __trans == 'c') {
    STATUS;
    __do_trmm_parallel_batched_template<Side::Right, Uplo::Upper, Trans::ConjTranspose, Diag::Unit, device_type>(
        options, trmm_args);
  }

  return;
}

/*************************** Internal setup fns **************************/
template <class scalar_type, class vta, class vtb, class device_type>
trmm_args_t __do_setup(options_t options, matrix_dims_t dim) {
  using execution_space = typename device_type::execution_space;

  trmm_args_t trmm_args;
  uint64_t seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(seed);
  decltype(dim.a.m) min_dim = dim.a.m < dim.a.n ? dim.a.m : dim.a.n;
  typename vta::HostMirror host_A;
  STATUS;

  trmm_args.side  = options.blas_args.trmm.trmm_args.c_str()[0];
  trmm_args.uplo  = options.blas_args.trmm.trmm_args.c_str()[1];
  trmm_args.trans = options.blas_args.trmm.trmm_args.c_str()[2];
  trmm_args.diag  = options.blas_args.trmm.trmm_args.c_str()[3];
  trmm_args.A     = vta("trmm_args.A", dim.a.k, dim.a.m, dim.a.n);
  trmm_args.B     = vtb("trmm_args.B", dim.b.k, dim.b.m, dim.b.n);
  trmm_args.alpha = options.blas_args.trmm.alpha;
  host_A          = Kokkos::create_mirror_view(trmm_args.A);

  {
    Kokkos::View<double***, default_layout, default_device> tmp("tmp", trmm_args.A.extent(0), trmm_args.A.extent(1),
                                                                trmm_args.A.extent(2));
    Kokkos::fill_random(tmp, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, double>::max());
    Kokkos::deep_copy(host_A, tmp);
  }

  if (trmm_args.uplo == 'U' || trmm_args.uplo == 'u') {
    // Make A upper triangular
    for (int k = 0; k < dim.a.k; ++k) {
      auto A = Kokkos::subview(host_A, k, Kokkos::ALL(), Kokkos::ALL());
      for (int i = 1; i < dim.a.m; i++) {
        for (int j = 0; j < i; j++) {
          A(i, j) = scalar_type(0);
        }
      }
    }
  } else {
    // Make A lower triangular
    // Kokkos::parallel_for("toLowerLoop", options.n, KOKKOS_LAMBDA (const int&
    // i) {
    for (int k = 0; k < dim.a.k; ++k) {
      auto A = Kokkos::subview(host_A, k, Kokkos::ALL(), Kokkos::ALL());
      for (int i = 0; i < dim.a.m - 1; i++) {
        for (int j = i + 1; j < dim.a.n; j++) {
          A(i, j) = scalar_type(0);
        }
      }
    }
  }

  if (trmm_args.diag == 'U' || trmm_args.diag == 'u') {
    for (int k = 0; k < dim.a.k; ++k) {
      auto A = Kokkos::subview(host_A, k, Kokkos::ALL(), Kokkos::ALL());
      for (int i = 0; i < min_dim; i++) {
        A(i, i) = scalar_type(1);
      }
    }
  }
  Kokkos::deep_copy(trmm_args.A, host_A);

  {
    Kokkos::View<double***, default_layout, default_device> tmp("tmp", trmm_args.B.extent(0), trmm_args.B.extent(1),
                                                                trmm_args.B.extent(2));
    Kokkos::fill_random(tmp, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, double>::max());
    Kokkos::deep_copy(trmm_args.B, tmp);
  }

  return trmm_args;
}

/*************************** Interal run helper fns **************************/
void __do_loop_and_invoke(options_t options, void (*fn)(options_t, trmm_args_t)) {
  matrix_dims_t cur_dims;
  trmm_args_t trmm_args;
  STATUS;

  __print_trmm_perf_test_options(options);
  std::cout << "SCALAR:" << typeid(default_scalar).name() << ", LAYOUT:" << typeid(default_layout).name()
            << ", DEVICE:" << typeid(default_device).name() << std::endl;

  options.out[0] << trmm_csv_header_str << std::endl;

  for (cur_dims = options.start; cur_dims.a.m <= options.stop.a.m && cur_dims.a.n <= options.stop.a.n &&
                                 cur_dims.b.m <= options.stop.b.m && cur_dims.b.n <= options.stop.b.n;
       cur_dims.a.m += options.step, cur_dims.a.n += options.step, cur_dims.b.m += options.step,
      cur_dims.b.n += options.step) {
    trmm_args = __do_setup<default_scalar, view_type_3d, view_type_3d, default_device>(options, cur_dims);
    fn(options, trmm_args);
  }
  return;
}

/*************************** External fns **************************/
void do_trmm_serial_blas(options_t options) {
  STATUS;
  __do_loop_and_invoke(options, __do_trmm_serial_blas<default_scalar, view_type_3d, view_type_3d, default_device>);
  return;
}

void do_trmm_serial_batched(options_t options) {
  STATUS;
  __do_loop_and_invoke(options, __do_trmm_serial_batched<default_scalar, view_type_3d, view_type_3d, default_device>);
  return;
}

void do_trmm_parallel_blas(options_t options) {
  STATUS;
  __do_loop_and_invoke(options, __do_trmm_parallel_blas<default_scalar, view_type_3d, view_type_3d, default_device>);
  return;
}

void do_trmm_parallel_batched(options_t options) {
  STATUS;
  __do_loop_and_invoke(options, __do_trmm_parallel_batched<default_scalar, view_type_3d, view_type_3d, default_device>);
  return;
}

#endif  // KOKKOSBLAS3_TRMM_PERF_TEST_H_
