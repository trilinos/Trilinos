/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#ifndef KOKKOSBLAS_TRTRI_PERF_TEST_H_
#define KOKKOSBLAS_TRTRI_PERF_TEST_H_

//#include <complex.h>
#include "KokkosBlas_common.hpp"

#include <Kokkos_Random.hpp>

#include <KokkosBlas_trtri.hpp>

#include "KokkosBatched_Trtri_Decl.hpp"
#include "KokkosBatched_Trtri_Serial_Impl.hpp"
#include "KokkosBatched_Util.hpp"

//#define TRTRI_PERF_TEST_DEBUG

// Forward declarations
void do_trtri_serial_blas(options_t options);
void do_trtri_serial_batched(options_t options);
void do_trtri_parallel_blas(options_t options);
void do_trtri_parallel_batched(options_t options);

// trtri invoke table
void (*do_trtri_invoke[LOOP_N][TEST_N])(options_t) = {
    {do_trtri_serial_blas, do_trtri_serial_batched},
    {do_trtri_parallel_blas, do_trtri_parallel_batched}};

/*************************** Print macros **************************/
#ifdef TRTRI_PERF_TEST_DEBUG
#define STATUS printf("STATUS: %s:%d.\n", __func__, __LINE__);
#else
#define STATUS
#endif  // TRTRI_PERF_TEST_DEBUG

/*************************** Test types and defaults **************************/
#define DEFAULT_TRTRI_ARGS "UU"

using view_type_3d =
    Kokkos::View<default_scalar***, default_layout, default_device>;
struct trtri_args {
  char uplo, diag;
  view_type_3d A;
};
typedef struct trtri_args trtri_args_t;

static std::string trtri_csv_header_str =
    "algorithm,side-uplo-trans-diag,alpha,loop_type,A_dims,warm_up_n,iter,"
    "total_time(s),average_time(s)";

/*************************** Internal helper fns **************************/
static void __trtri_output_csv_row(options_t options, trtri_args_t trtri_args,
                                   double time_in_seconds) {
  options.out[0] << test_e_str[options.test] << ","
                 << options.blas_args.trtri.trtri_args << ","
                 << loop_e_str[options.loop] << "," << trtri_args.A.extent(1)
                 << "x" << trtri_args.A.extent(2) << "," << options.warm_up_n
                 << "," << options.n << "," << time_in_seconds << ","
                 << time_in_seconds / options.n << std::endl;
}

static void __print_trtri_perf_test_options(options_t options) {
#ifdef TRTRI_PERF_TEST_DEBUG
  printf("options.test      = %s\n", test_e_str[options.test].c_str());
  printf("options.loop      = %s\n", loop_e_str[options.loop].c_str());
  printf("options.start     = %dx%d,%dx%d\n", options.start.a.m,
         options.start.a.n, options.start.b.m, options.start.b.n);
  printf("options.stop      = %dx%d,%d,%d\n", options.stop.a.m,
         options.stop.a.n, options.stop.b.m, options.stop.b.n);
  printf("options.step      = %d\n", options.step);
  printf("options.warm_up_n = %d\n", options.warm_up_n);
  printf("options.n         = %d\n", options.n);
  printf("options.blas_args.trtri.trtri_args = %s\n",
         options.blas_args.trtri.trtri_args.c_str());
  printf("options.out_file  = %s\n", options.out_file.c_str());
  std::cout << "SCALAR:" << typeid(default_scalar).name()
            << ", LAYOUT:" << typeid(default_layout).name() << ", DEVICE:."
            << typeid(default_device).name() << std::endl;
#endif  // TRTRI_PERF_TEST_DEBUG
  return;
}

/*************************** Internal templated fns **************************/
template <class scalar_type, class vta, class device_type>
void __do_trtri_serial_blas(options_t options, trtri_args_t trtri_args) {
// Need to take subviews on the device
#if !defined(KOKKOS_ENABLE_CUDA)
  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  Kokkos::Timer timer;

  STATUS;

  for (uint32_t i = 0; i < warm_up_n; ++i) {
    auto A = Kokkos::subview(trtri_args.A, i, Kokkos::ALL(), Kokkos::ALL());

    KokkosBlas::trtri(&trtri_args.uplo, &trtri_args.diag, A);
  }

  timer.reset();
  for (uint32_t i = 0; i < n; ++i) {
    auto A = Kokkos::subview(trtri_args.A, i, Kokkos::ALL(), Kokkos::ALL());

    KokkosBlas::trtri(&trtri_args.uplo, &trtri_args.diag, A);
  }
  Kokkos::fence();
  __trtri_output_csv_row(options, trtri_args, timer.seconds());
#else
  std::cerr << std::string(__func__)
            << " disabled since KOKKOS_ENABLE_CUDA is defined." << std::endl;
#endif  // !KOKKOS_ENABLE_CUDA
  return;
}

template <class uplo, class diag>
void __do_trtri_serial_batched_template(options_t options,
                                        trtri_args_t trtri_args) {
// Need to take subviews on the device
#if !defined(KOKKOS_ENABLE_CUDA)
  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  Kokkos::Timer timer;
  using tag = Algo::Trtri::Unblocked;

  for (uint32_t i = 0; i < warm_up_n; ++i) {
    auto A = Kokkos::subview(trtri_args.A, i, Kokkos::ALL(), Kokkos::ALL());

    SerialTrtri<uplo, diag, tag>::invoke(A);
  }

  timer.reset();
  for (uint32_t i = 0; i < n; ++i) {
    auto A = Kokkos::subview(trtri_args.A, i, Kokkos::ALL(), Kokkos::ALL());

    SerialTrtri<uplo, diag, tag>::invoke(A);
  }
  Kokkos::fence();
  __trtri_output_csv_row(options, trtri_args, timer.seconds());
#else
  std::cerr << std::string(__func__)
            << " disabled since KOKKOS_ENABLE_CUDA is defined." << std::endl;
#endif  // !KOKKOS_ENABLE_CUDA
}

template <class scalar_type, class vta, class device_type>
void __do_trtri_serial_batched(options_t options, trtri_args_t trtri_args) {
  char __uplo = tolower(trtri_args.uplo), __diag = tolower(trtri_args.diag);

  STATUS;

  //// Lower ////
  if (__uplo == 'l') {
    if (__diag == 'u') {
      __do_trtri_serial_batched_template<Uplo::Lower, Diag::Unit>(options,
                                                                  trtri_args);
    } else {
      __do_trtri_serial_batched_template<Uplo::Lower, Diag::NonUnit>(
          options, trtri_args);
    }
  } else {
    //// Upper ////
    if (__diag == 'u') {
      __do_trtri_serial_batched_template<Uplo::Upper, Diag::Unit>(options,
                                                                  trtri_args);
    } else {
      __do_trtri_serial_batched_template<Uplo::Upper, Diag::NonUnit>(
          options, trtri_args);
    }
  }

  return;
}

#if !defined(KOKKOS_ENABLE_CUDA)
template <class ExecutionSpace>
struct parallel_blas_trtri {
  trtri_args_t trtri_args_;

  parallel_blas_trtri(trtri_args_t trtri_args) : trtri_args_(trtri_args) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i) const {
    auto svA = Kokkos::subview(trtri_args_.A, i, Kokkos::ALL(), Kokkos::ALL());

    KokkosBlas::trtri(&trtri_args_.uplo, &trtri_args_.diag, svA);
  }
};
#endif  // !KOKKOS_ENABLE_CUDA

template <class scalar_type, class vta, class device_type>
void __do_trtri_parallel_blas(options_t options, trtri_args_t trtri_args) {
#if !defined(KOKKOS_ENABLE_CUDA)
  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  Kokkos::Timer timer;
  using execution_space = typename device_type::execution_space;
  using functor_type    = parallel_blas_trtri<execution_space>;
  functor_type parallel_blas_trtri_functor(trtri_args);

  STATUS;

  Kokkos::parallel_for("parallelBlasWarmUpLoopTrtri",
                       Kokkos::RangePolicy<execution_space>(0, warm_up_n),
                       parallel_blas_trtri_functor);
  Kokkos::fence();

  timer.reset();
  Kokkos::parallel_for("parallelBlasTimedLoopTrtri",
                       Kokkos::RangePolicy<execution_space>(0, n),
                       parallel_blas_trtri_functor);
  Kokkos::fence();
  __trtri_output_csv_row(options, trtri_args, timer.seconds());
#else
  std::cerr << std::string(__func__)
            << " disabled since KOKKOS_ENABLE_CUDA is defined." << std::endl;
  __trtri_output_csv_row(options, trtri_args, -1);
#endif  // !KOKKOS_ENABLE_CUDA
  return;
}

template <class uplo, class diag, class tag, class ExecutionSpace>
struct parallel_batched_trtri {
  trtri_args_t trtri_args_;

  parallel_batched_trtri(trtri_args_t trtri_args) : trtri_args_(trtri_args) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i) const {
    auto svA = Kokkos::subview(trtri_args_.A, i, Kokkos::ALL(), Kokkos::ALL());

    SerialTrtri<uplo, diag, tag>::invoke(svA);
  }
};

template <class uplo, class diag, class device_type>
void __do_trtri_parallel_batched_template(options_t options,
                                          trtri_args_t trtri_args) {
  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  Kokkos::Timer timer;
  using tag             = Algo::Trtri::Unblocked;
  using execution_space = typename device_type::execution_space;
  using functor_type = parallel_batched_trtri<uplo, diag, tag, execution_space>;
  functor_type parallel_batched_trtri_functor(trtri_args);

  STATUS;

  Kokkos::parallel_for("parallelBatchedWarmUpLoopTrtri",
                       Kokkos::RangePolicy<execution_space>(0, warm_up_n),
                       parallel_batched_trtri_functor);
  Kokkos::fence();

  timer.reset();
  Kokkos::parallel_for("parallelBatchedTimedLoopTrtri",
                       Kokkos::RangePolicy<execution_space>(0, n),
                       parallel_batched_trtri_functor);
  Kokkos::fence();
  __trtri_output_csv_row(options, trtri_args, timer.seconds());

  return;
}

template <class scalar_type, class vta, class device_type>
void __do_trtri_parallel_batched(options_t options, trtri_args_t trtri_args) {
  char __uplo = tolower(trtri_args.uplo), __diag = tolower(trtri_args.diag);

  STATUS;

  //// Lower ////
  if (__uplo == 'l') {
    if (__diag == 'u') {
      __do_trtri_parallel_batched_template<Uplo::Lower, Diag::Unit,
                                           device_type>(options, trtri_args);
    } else {
      __do_trtri_parallel_batched_template<Uplo::Lower, Diag::NonUnit,
                                           device_type>(options, trtri_args);
    }
  } else {
    //// Upper ////
    if (__diag == 'u') {
      __do_trtri_parallel_batched_template<Uplo::Upper, Diag::Unit,
                                           device_type>(options, trtri_args);
    } else {
      __do_trtri_parallel_batched_template<Uplo::Upper, Diag::NonUnit,
                                           device_type>(options, trtri_args);
    }
  }

  return;
}

/*************************** Internal setup fns **************************/
template <class scalar_type, class vta, class device_type>
trtri_args_t __do_setup(options_t options, matrix_dims_t dim) {
  using execution_space = typename device_type::execution_space;

  trtri_args_t trtri_args;
  uint64_t seed = Kokkos::Impl::clock_tic();
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(seed);
  decltype(dim.a.m) min_dim = dim.a.m < dim.a.n ? dim.a.m : dim.a.n;
  typename vta::HostMirror host_A;
  STATUS;

  trtri_args.uplo = options.blas_args.trtri.trtri_args.c_str()[0];
  trtri_args.diag = options.blas_args.trtri.trtri_args.c_str()[1];
  trtri_args.A    = vta("trtri_args.A", options.n, dim.a.m, dim.a.n);
  host_A          = Kokkos::create_mirror_view(trtri_args.A);

  Kokkos::fill_random(trtri_args.A, rand_pool,
                      Kokkos::rand<Kokkos::Random_XorShift64<execution_space>,
                                   scalar_type>::max());
  Kokkos::deep_copy(host_A, trtri_args.A);

  if (trtri_args.uplo == 'U' || trtri_args.uplo == 'u') {
    // Make A upper triangular
    for (uint32_t k = 0; k < options.n; ++k) {
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
    for (uint32_t k = 0; k < options.n; ++k) {
      auto A = Kokkos::subview(host_A, k, Kokkos::ALL(), Kokkos::ALL());
      for (int i = 0; i < dim.a.m - 1; i++) {
        for (int j = i + 1; j < dim.a.n; j++) {
          A(i, j) = scalar_type(0);
        }
      }
    }
  }

  if (trtri_args.diag == 'U' || trtri_args.diag == 'u') {
    for (uint32_t k = 0; k < options.n; ++k) {
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
void __do_loop_and_invoke(options_t options,
                          void (*fn)(options_t, trtri_args_t)) {
  matrix_dims_t cur_dims;
  trtri_args_t trtri_args;
  STATUS;

  __print_trtri_perf_test_options(options);
  std::cout << "SCALAR:" << typeid(default_scalar).name()
            << ", LAYOUT:" << typeid(default_layout).name() << ", DEVICE:."
            << typeid(default_device).name() << std::endl;

  options.out[0] << trtri_csv_header_str << std::endl;

  for (cur_dims = options.start;
       cur_dims.a.m <= options.stop.a.m && cur_dims.a.n <= options.stop.a.n &&
       cur_dims.b.m <= options.stop.b.m && cur_dims.b.n <= options.stop.b.n;
       cur_dims.a.m *= options.step, cur_dims.a.n *= options.step,
      cur_dims.b.m *= options.step, cur_dims.b.n *= options.step) {
    trtri_args = __do_setup<default_scalar, view_type_3d, default_device>(
        options, cur_dims);
    fn(options, trtri_args);
  }
  return;
}

/*************************** External fns **************************/
void do_trtri_serial_blas(options_t options) {
  STATUS;
  __do_loop_and_invoke(
      options,
      __do_trtri_serial_blas<default_scalar, view_type_3d, default_device>);
  return;
}

void do_trtri_serial_batched(options_t options) {
  STATUS;
  __do_loop_and_invoke(
      options,
      __do_trtri_serial_batched<default_scalar, view_type_3d, default_device>);
  return;
}

void do_trtri_parallel_blas(options_t options) {
  STATUS;
  __do_loop_and_invoke(
      options,
      __do_trtri_parallel_blas<default_scalar, view_type_3d, default_device>);
  return;
}

void do_trtri_parallel_batched(options_t options) {
  STATUS;
  __do_loop_and_invoke(options,
                       __do_trtri_parallel_batched<default_scalar, view_type_3d,
                                                   default_device>);
  return;
}

#endif  // KOKKOSBLAS_TRTRI_PERF_TEST_H_
