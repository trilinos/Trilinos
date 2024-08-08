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
#ifndef KOKKOSBLAS_TRTRI_PERF_TEST_H_
#define KOKKOSBLAS_TRTRI_PERF_TEST_H_

// #include <complex.h>
#include "KokkosBlas_common.hpp"

#include <Kokkos_Random.hpp>

#include <KokkosLapack_trtri.hpp>

#include "KokkosBatched_Trtri_Decl.hpp"
#include "KokkosBatched_Trtri_Serial_Impl.hpp"
#include "KokkosBatched_Util.hpp"

#include <chrono>

// #define TRTRI_PERF_TEST_DEBUG

// Forward declarations
void do_trtri_serial_blas(options_t options);
void do_trtri_serial_batched(options_t options);
void do_trtri_parallel_blas(options_t options);
void do_trtri_parallel_batched(options_t options);

// trtri invoke table
void (*do_trtri_invoke[LOOP_N][TEST_N])(options_t) = {{do_trtri_serial_blas, do_trtri_serial_batched},
                                                      {do_trtri_parallel_blas, do_trtri_parallel_batched}};

/*************************** Print macros **************************/
#ifdef TRTRI_PERF_TEST_DEBUG
#define STATUS printf("STATUS: %s:%d.\n", __func__, __LINE__);
#else
#define STATUS
#endif  // TRTRI_PERF_TEST_DEBUG

/*************************** Test types and defaults **************************/
#define DEFAULT_TRTRI_ARGS "UU"

/**
 * The KokkosBatched::SerialTrtri implementation performs trmm and scal on
 * subblocks of the A matrix. a_m subblocks are selected.
 */
static inline double __trtri_impl_flop_count(double a_m, double /*a_n*/) {
  double flop_count = 0;
  double flops_per_div, flops_per_mul, flops_per_add;

  if (std::is_same<double, default_scalar>::value || std::is_same<float, default_scalar>::value ||
      std::is_same<Kokkos::Experimental::half_t, default_scalar>::value) {
    flops_per_div = 1;
    flops_per_mul = 1;
    flops_per_add = 1;
  } else {
    // For complex, we need to count 2 flops for each add and 6 flops for each
    // multiply or divide.
    flops_per_div = 6;
    flops_per_mul = 6;
    flops_per_add = 2;
  }

  for (int i = 0; i < a_m; i++) {
    flop_count += flops_per_div;                                          // 1 / A[i,j]
    flop_count += ((i * (i + 1)) / 2) * (flops_per_mul + flops_per_add);  // TRMM FLOPS
    flop_count += i * flops_per_mul;                                      // SCAL FLOPS
  }

  return flop_count;
}

// Flop count formula from lapack working note 41:
// http://www.icl.utk.edu/~mgates3/docs/lawn41.pdf
static inline double __trtri_flop_count(double a_m, double a_n) {
  double flops;
  double flops_per_mul;
  double flops_per_add;

  if (a_m != a_n) {
    fprintf(stderr, "%s:%d:ERROR: a_m != a_n.\n", __FILE__, __LINE__);
    exit(255);
  }

  if (std::is_same<double, default_scalar>::value || std::is_same<float, default_scalar>::value ||
      std::is_same<Kokkos::Experimental::half_t, default_scalar>::value) {
    flops_per_mul = 1;
    flops_per_add = 1;
  } else {
    // For complex, we need to count 2 flops for each add and 6 flops for each
    // multiply.
    flops_per_mul = 6;
    flops_per_add = 2;
  }

  flops = (1. / 6. * a_n * a_n * a_n + 1. / 2. * a_n * a_n + 1. / 3. * a_n) * flops_per_mul +
          (1. / 6. * a_n * a_n * a_n - 1. / 2. * a_n * a_n + 1. / 3. * a_n) * flops_per_add;

  return flops;
}

using view_type_3d = Kokkos::View<default_scalar***, default_layout, default_device>;
struct trtri_args {
  char uplo, diag;
  view_type_3d A;
};
typedef struct trtri_args trtri_args_t;

static std::string trtri_csv_header_str =
    "algorithm,side-uplo-trans-diag,loop_type,A_dims,warm_up_n,iter,"
    "total_time(s),average_time(s),FLOPS,GFLOP/average_time(s)";

/*************************** Internal helper fns **************************/
static void __trtri_output_csv_row(options_t options, trtri_args_t trtri_args, double time_in_seconds) {
  double flops        = trtri_args.A.extent(0) * __trtri_flop_count(trtri_args.A.extent(1), trtri_args.A.extent(2));
  double gflops       = flops / 1e9;
  double average_time = time_in_seconds / options.n;

  options.out[0] << test_e_str[options.test] << "," << options.blas_args.trtri.trtri_args << ","
                 << loop_e_str[options.loop] << "," << trtri_args.A.extent(0) << "x" << trtri_args.A.extent(1) << "x"
                 << trtri_args.A.extent(2) << "," << options.warm_up_n << "," << options.n << "," << time_in_seconds
                 << "," << average_time << "," << flops << "," << gflops / average_time << std::endl;
}

#ifdef TRTRI_PERF_TEST_DEBUG
static void __print_trtri_perf_test_options(options_t options) {
  printf("options.test      = %s\n", test_e_str[options.test].c_str());
  printf("options.loop      = %s\n", loop_e_str[options.loop].c_str());
  printf("options.start     = %dx%d,%dx%d\n", options.start.a.m, options.start.a.n, options.start.b.m,
         options.start.b.n);
  printf("options.stop      = %dx%d,%d,%d\n", options.stop.a.m, options.stop.a.n, options.stop.b.m, options.stop.b.n);
  printf("options.step      = %d\n", options.step);
  printf("options.warm_up_n = %d\n", options.warm_up_n);
  printf("options.n         = %d\n", options.n);
  printf("options.blas_args.trtri.trtri_args = %s\n", options.blas_args.trtri.trtri_args.c_str());
  printf("options.out_file  = %s\n", options.out_file.c_str());
  std::cout << "SCALAR:" << typeid(default_scalar).name() << ", LAYOUT:" << typeid(default_layout).name()
            << ", DEVICE:." << typeid(default_device).name() << std::endl;
#else
static void __print_trtri_perf_test_options(options_t) {
#endif  // TRTRI_PERF_TEST_DEBUG
  return;
}

/*************************** Internal templated fns **************************/
#if !defined(KOKKOS_ENABLE_CUDA) && !defined(KOKKOS_ENABLE_HIP) && !defined(KOKKOS_ENABLE_OPENMPTARGET)
template <class scalar_type, class vta, class device_type>
void __do_trtri_serial_blas(options_t options, trtri_args_t trtri_args) {
  // Need to take subviews on the device
  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  Kokkos::Timer timer;

  STATUS;

  for (uint32_t j = 0; j < warm_up_n; ++j) {
    for (int i = 0; i < options.start.a.k; ++i) {
      auto A = Kokkos::subview(trtri_args.A, i, Kokkos::ALL(), Kokkos::ALL());

      KokkosLapack::trtri(&trtri_args.uplo, &trtri_args.diag, A);
    }
    // Fence after each batch operation
    Kokkos::fence();
  }

  timer.reset();
  for (uint32_t j = 0; j < n; ++j) {
    for (int i = 0; i < options.start.a.k; ++i) {
      auto A = Kokkos::subview(trtri_args.A, i, Kokkos::ALL(), Kokkos::ALL());

      KokkosLapack::trtri(&trtri_args.uplo, &trtri_args.diag, A);
    }
    // Fence after each batch operation
    Kokkos::fence();
  }
  __trtri_output_csv_row(options, trtri_args, timer.seconds());
  return;
}
#else
template <class scalar_type, class vta, class device_type>
void __do_trtri_serial_blas(options_t /*options*/, trtri_args_t /*trtri_args*/) {
  std::cerr << std::string(__func__) << " disabled since KOKKOS_ENABLE_DEVICE is defined." << std::endl;
  return;
}
#endif  // !KOKKOS_ENABLE_CUDA

// Need to take subviews on the device
#if !defined(KOKKOS_ENABLE_CUDA) && !defined(KOKKOS_ENABLE_HIP) && !defined(KOKKOS_ENABLE_OPENMPTARGET)
template <class uplo, class diag>
void __do_trtri_serial_batched_template(options_t options, trtri_args_t trtri_args) {
  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  Kokkos::Timer timer;
  using tag = KokkosBatched::Algo::Trtri::Unblocked;

  for (uint32_t j = 0; j < warm_up_n; ++j) {
    for (int i = 0; i < options.start.a.k; ++i) {
      auto A = Kokkos::subview(trtri_args.A, i, Kokkos::ALL(), Kokkos::ALL());

      KokkosBatched::SerialTrtri<uplo, diag, tag>::invoke(A);
    }
    // Fence after each batch operation
    Kokkos::fence();
  }

  timer.reset();
  for (uint32_t j = 0; j < n; ++j) {
    for (int i = 0; i < options.start.a.k; ++i) {
      auto A = Kokkos::subview(trtri_args.A, i, Kokkos::ALL(), Kokkos::ALL());

      KokkosBatched::SerialTrtri<uplo, diag, tag>::invoke(A);
    }
    // Fence after each batch operation
    Kokkos::fence();
  }
  __trtri_output_csv_row(options, trtri_args, timer.seconds());
}
#else
template <class uplo, class diag>
void __do_trtri_serial_batched_template(options_t /*options*/, trtri_args_t /*trtri_args*/) {
  std::cerr << std::string(__func__) << " disabled since KOKKOS_ENABLE_DEVICE is defined." << std::endl;
}
#endif  // !KOKKOS_ENABLE_CUDA

template <class scalar_type, class vta, class device_type>
void __do_trtri_serial_batched(options_t options, trtri_args_t trtri_args) {
  using KokkosBatched::Diag;
  using KokkosBatched::Uplo;

  char __uplo = tolower(trtri_args.uplo), __diag = tolower(trtri_args.diag);

  STATUS;

  //// Lower ////
  if (__uplo == 'l') {
    if (__diag == 'u') {
      __do_trtri_serial_batched_template<Uplo::Lower, Diag::Unit>(options, trtri_args);
    } else {
      __do_trtri_serial_batched_template<Uplo::Lower, Diag::NonUnit>(options, trtri_args);
    }
  } else {
    //// Upper ////
    if (__diag == 'u') {
      __do_trtri_serial_batched_template<Uplo::Upper, Diag::Unit>(options, trtri_args);
    } else {
      __do_trtri_serial_batched_template<Uplo::Upper, Diag::NonUnit>(options, trtri_args);
    }
  }

  return;
}

#if !defined(KOKKOS_ENABLE_CUDA) && !defined(KOKKOS_ENABLE_HIP) && !defined(KOKKOS_ENABLE_OPENMPTARGET)
template <class ExecutionSpace>
struct parallel_blas_trtri {
  trtri_args_t trtri_args_;

  parallel_blas_trtri(trtri_args_t trtri_args) : trtri_args_(trtri_args) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i) const {
    auto svA = Kokkos::subview(trtri_args_.A, i, Kokkos::ALL(), Kokkos::ALL());

    KokkosLapack::trtri(&trtri_args_.uplo, &trtri_args_.diag, svA);
  }
};
#endif  // !KOKKOS_ENABLE_CUDA && !KOKKOS_ENABLE_HIP &&
        // !KOKKOS_ENABLE_OPENMPTARGET

template <class scalar_type, class vta, class device_type>
void __do_trtri_parallel_blas(options_t options, trtri_args_t trtri_args) {
#if !defined(KOKKOS_ENABLE_CUDA) && !defined(KOKKOS_ENABLE_HIP) && !defined(KOKKOS_ENABLE_OPENMPTARGET)
  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  Kokkos::Timer timer;
  using execution_space = typename device_type::execution_space;
  using functor_type    = parallel_blas_trtri<execution_space>;
  functor_type parallel_blas_trtri_functor(trtri_args);

  STATUS;

  for (uint32_t i = 0; i < warm_up_n; ++i) {
    Kokkos::parallel_for("parallelBlasWarmUpLoopTrtri", Kokkos::RangePolicy<execution_space>(0, options.start.a.k),
                         parallel_blas_trtri_functor);
    // Fence after each batch operation
    Kokkos::fence();
  }

  timer.reset();
  for (uint32_t i = 0; i < n; ++i) {
    Kokkos::parallel_for("parallelBlasTimedLoopTrtri", Kokkos::RangePolicy<execution_space>(0, options.start.a.k),
                         parallel_blas_trtri_functor);
    // Fence after each batch operation
    Kokkos::fence();
  }
  __trtri_output_csv_row(options, trtri_args, timer.seconds());
#else
  std::cerr << std::string(__func__)
            << " disabled since KOKKOS_ENABLE_CUDA, KOKKOS_ENABLE_HIP or "
               "KOKKOS_ENABLE_OPENMPTARGET is defined."
            << std::endl;
  __trtri_output_csv_row(options, trtri_args, -1);
#endif  // !KOKKOS_ENABLE_CUDA && !KOKKOS_ENABLE_HIP &&
        // !defined(KOKKOS_ENABLE_OPENMPTARGET)
  return;
}

template <class uplo, class diag, class tag, class ExecutionSpace>
struct parallel_batched_trtri {
  trtri_args_t trtri_args_;

  parallel_batched_trtri(trtri_args_t trtri_args) : trtri_args_(trtri_args) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i) const {
    auto svA = Kokkos::subview(trtri_args_.A, i, Kokkos::ALL(), Kokkos::ALL());

    KokkosBatched::SerialTrtri<uplo, diag, tag>::invoke(svA);
  }
};

template <class uplo, class diag, class device_type>
void __do_trtri_parallel_batched_template(options_t options, trtri_args_t trtri_args) {
  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  Kokkos::Timer timer;
  using tag             = KokkosBatched::Algo::Trtri::Unblocked;
  using execution_space = typename device_type::execution_space;
  using functor_type    = parallel_batched_trtri<uplo, diag, tag, execution_space>;
  functor_type parallel_batched_trtri_functor(trtri_args);

  STATUS;

  for (uint32_t i = 0; i < warm_up_n; ++i) {
    Kokkos::parallel_for("parallelBatchedWarmUpLoopTrtri", Kokkos::RangePolicy<execution_space>(0, options.start.a.k),
                         parallel_batched_trtri_functor);
    // Fence after each batch operation
    Kokkos::fence();
  }

  timer.reset();

  for (uint32_t i = 0; i < n; ++i) {
    Kokkos::parallel_for("parallelBatchedTimedLoopTrtri", Kokkos::RangePolicy<execution_space>(0, options.start.a.k),
                         parallel_batched_trtri_functor);
    // Fence after each batch operation
    Kokkos::fence();
  }
  __trtri_output_csv_row(options, trtri_args, timer.seconds());

  return;
}

template <class scalar_type, class vta, class device_type>
void __do_trtri_parallel_batched(options_t options, trtri_args_t trtri_args) {
  using KokkosBatched::Diag;
  using KokkosBatched::Uplo;

  char __uplo = tolower(trtri_args.uplo), __diag = tolower(trtri_args.diag);

  STATUS;

  //// Lower ////
  if (__uplo == 'l') {
    if (__diag == 'u') {
      __do_trtri_parallel_batched_template<Uplo::Lower, Diag::Unit, device_type>(options, trtri_args);
    } else {
      __do_trtri_parallel_batched_template<Uplo::Lower, Diag::NonUnit, device_type>(options, trtri_args);
    }
  } else {
    //// Upper ////
    if (__diag == 'u') {
      __do_trtri_parallel_batched_template<Uplo::Upper, Diag::Unit, device_type>(options, trtri_args);
    } else {
      __do_trtri_parallel_batched_template<Uplo::Upper, Diag::NonUnit, device_type>(options, trtri_args);
    }
  }

  return;
}

/*************************** Internal setup fns **************************/
template <class scalar_type, class vta, class device_type>
trtri_args_t __do_setup(options_t options, matrix_dims_t dim) {
  using execution_space = typename device_type::execution_space;

  trtri_args_t trtri_args;
  uint64_t seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(seed);
  decltype(dim.a.m) min_dim = dim.a.m < dim.a.n ? dim.a.m : dim.a.n;
  typename vta::HostMirror host_A;
  STATUS;

  trtri_args.uplo = options.blas_args.trtri.trtri_args.c_str()[0];
  trtri_args.diag = options.blas_args.trtri.trtri_args.c_str()[1];
  trtri_args.A    = vta("trtri_args.A", dim.a.k, dim.a.m, dim.a.n);
  host_A          = Kokkos::create_mirror_view(trtri_args.A);

  Kokkos::fill_random(trtri_args.A, rand_pool,
                      Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, scalar_type>::max());
  Kokkos::deep_copy(host_A, trtri_args.A);

  if (trtri_args.uplo == 'U' || trtri_args.uplo == 'u') {
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

  if (trtri_args.diag == 'U' || trtri_args.diag == 'u') {
    for (int k = 0; k < dim.a.k; ++k) {
      auto A = Kokkos::subview(host_A, k, Kokkos::ALL(), Kokkos::ALL());
      for (int i = 0; i < min_dim; i++) {
        A(i, i) = scalar_type(1);
      }
    }
  }

  Kokkos::deep_copy(trtri_args.A, host_A);

  return trtri_args;
}

/*************************** Interal run helper fns **************************/
void __do_loop_and_invoke(options_t options, void (*fn)(options_t, trtri_args_t)) {
  matrix_dims_t cur_dims;
  trtri_args_t trtri_args;
  STATUS;

  __print_trtri_perf_test_options(options);
  std::cout << "SCALAR:" << typeid(default_scalar).name() << ", LAYOUT:" << typeid(default_layout).name()
            << ", DEVICE:." << typeid(default_device).name() << std::endl;

  options.out[0] << trtri_csv_header_str << std::endl;

  for (cur_dims = options.start; cur_dims.a.m <= options.stop.a.m && cur_dims.a.n <= options.stop.a.n &&
                                 cur_dims.b.m <= options.stop.b.m && cur_dims.b.n <= options.stop.b.n;
       cur_dims.a.m += options.step, cur_dims.a.n += options.step, cur_dims.b.m += options.step,
      cur_dims.b.n += options.step) {
    trtri_args = __do_setup<default_scalar, view_type_3d, default_device>(options, cur_dims);
    fn(options, trtri_args);
  }
  return;
}

/*************************** External fns **************************/
void do_trtri_serial_blas(options_t options) {
  STATUS;
  __do_loop_and_invoke(options, __do_trtri_serial_blas<default_scalar, view_type_3d, default_device>);
  return;
}

void do_trtri_serial_batched(options_t options) {
  STATUS;
  __do_loop_and_invoke(options, __do_trtri_serial_batched<default_scalar, view_type_3d, default_device>);
  return;
}

void do_trtri_parallel_blas(options_t options) {
  STATUS;
  __do_loop_and_invoke(options, __do_trtri_parallel_blas<default_scalar, view_type_3d, default_device>);
  return;
}

void do_trtri_parallel_batched(options_t options) {
  STATUS;
  __do_loop_and_invoke(options, __do_trtri_parallel_batched<default_scalar, view_type_3d, default_device>);
  return;
}

#endif  // KOKKOSBLAS_TRTRI_PERF_TEST_H_
