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
#ifndef KOKKOSBLAS3_TRMM_PERF_TEST_H_
#define KOKKOSBLAS3_TRMM_PERF_TEST_H_

//#include <complex.h>
#include "KokkosBlas3_common.hpp"

#include <Kokkos_Random.hpp>

#include <KokkosBlas3_trmm.hpp>

#include "KokkosBatched_Trmm_Decl.hpp"
#include "KokkosBatched_Trmm_Serial_Impl.hpp"
#include "KokkosBatched_Util.hpp"

//#define TRMM_PERF_TEST_DEBUG

// Forward declarations
void do_trmm_serial_blas(options_t options);
void do_trmm_serial_batched(options_t options);
void do_trmm_parallel_blas(options_t options);
void do_trmm_parallel_batched(options_t options);

// trmm invoke table
void (*do_trmm_invoke[LOOP_N][TEST_N])(options_t) = {
    {do_trmm_serial_blas, do_trmm_serial_batched},
    {do_trmm_parallel_blas, do_trmm_parallel_batched}};

/*************************** Print macros **************************/
#ifdef TRMM_PERF_TEST_DEBUG
#define STATUS printf("STATUS: %s:%d.\n", __func__, __LINE__);
#else
#define STATUS
#endif  // TRMM_PERF_TEST_DEBUG

/*************************** Test types and defaults **************************/
#define DEFAULT_TRMM_ARGS "LUNU"
#define DEFAULT_TRMM_ALPHA 1.0

using view_type_3d =
    Kokkos::View<default_scalar***, default_layout, default_device>;
struct trmm_args {
  char side, uplo, trans, diag;
  default_scalar alpha;
  view_type_3d A, B;
};
typedef struct trmm_args trmm_args_t;

static std::string trmm_csv_header_str =
    "algorithm,side-uplo-trans-diag,alpha,loop_type,A_dims,B_dims,warm_up_n,"
    "iter,total_time(s),average_time(s)";

/*************************** Internal helper fns **************************/
static void __trmm_output_csv_row(options_t options, trmm_args_t trmm_args,
                                  double time_in_seconds) {
  options.out[0] << test_e_str[options.test] << ","
                 << options.blas_args.trmm.trmm_args << ","
                 << options.blas_args.trmm.alpha << ","
                 << loop_e_str[options.loop] << "," << trmm_args.A.extent(1)
                 << "x" << trmm_args.A.extent(2) << "," << trmm_args.B.extent(1)
                 << "x" << trmm_args.B.extent(2) << "," << options.warm_up_n
                 << "," << options.n << "," << time_in_seconds << ","
                 << time_in_seconds / options.n << std::endl;
}

static void __print_trmm_perf_test_options(options_t options) {
#ifdef TRMM_PERF_TEST_DEBUG
  printf("options.test      = %s\n", test_e_str[options.test].c_str());
  printf("options.loop      = %s\n", loop_e_str[options.loop].c_str());
  printf("options.start     = %dx%d,%dx%d\n", options.start.a.m,
         options.start.a.n, options.start.b.m, options.start.b.n);
  printf("options.stop      = %dx%d,%dx%d\n", options.stop.a.m,
         options.stop.a.n, options.stop.b.m, options.stop.b.n);
  printf("options.step      = %d\n", options.step);
  printf("options.warm_up_n = %d\n", options.warm_up_n);
  printf("options.n         = %d\n", options.n);
  printf("options.blas_args.trmm.trmm_args = %s\n",
         options.blas_args.trmm.trmm_args.c_str());
  printf("options.out_file  = %s\n", options.out_file.c_str());
  if (std::is_same<double, default_scalar>::value)
    printf("options.alpha     = %lf\n", options.blas_args.trmm.alpha);
  else if (std::is_same<float, default_scalar>::value)
    printf("options.alpha     = %f\n", options.blas_args.trmm.alpha);
#endif  // TRMM_PERF_TEST_DEBUG
  return;
}

/*************************** Internal templated fns **************************/
template <class scalar_type, class vta, class vtb, class device_type>
void __do_trmm_serial_blas(options_t options, trmm_args_t trmm_args) {
// Need to take subviews on the device
#if !defined(KOKKOS_ENABLE_CUDA)
  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  Kokkos::Timer timer;

  STATUS;

  for (uint32_t i = 0; i < warm_up_n; ++i) {
    auto A = Kokkos::subview(trmm_args.A, i, Kokkos::ALL(), Kokkos::ALL());
    auto B = Kokkos::subview(trmm_args.B, i, Kokkos::ALL(), Kokkos::ALL());

    KokkosBlas::trmm(&trmm_args.side, &trmm_args.uplo, &trmm_args.trans,
                     &trmm_args.diag, trmm_args.alpha, A, B);
  }

  timer.reset();
  for (uint32_t i = 0; i < n; ++i) {
    auto A = Kokkos::subview(trmm_args.A, i, Kokkos::ALL(), Kokkos::ALL());
    auto B = Kokkos::subview(trmm_args.B, i, Kokkos::ALL(), Kokkos::ALL());

    KokkosBlas::trmm(&trmm_args.side, &trmm_args.uplo, &trmm_args.trans,
                     &trmm_args.diag, trmm_args.alpha, A, B);
  }
  Kokkos::fence();
  __trmm_output_csv_row(options, trmm_args, timer.seconds());
#else
  std::cerr << std::string(__func__)
            << " disabled since KOKKOS_ENABLE_CUDA is defined." << std::endl;
#endif  // !KOKKOS_ENABLE_CUDA
  return;
}

template <class side, class uplo, class trans, class diag>
void __do_trmm_serial_batched_template(options_t options,
                                       trmm_args_t trmm_args) {
// Need to take subviews on the device
#if !defined(KOKKOS_ENABLE_CUDA)
  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  Kokkos::Timer timer;
  using tag = Algo::Trmm::Unblocked;

  for (uint32_t i = 0; i < warm_up_n; ++i) {
    auto A = Kokkos::subview(trmm_args.A, i, Kokkos::ALL(), Kokkos::ALL());
    auto B = Kokkos::subview(trmm_args.B, i, Kokkos::ALL(), Kokkos::ALL());

    SerialTrmm<side, uplo, trans, diag, tag>::invoke(trmm_args.alpha, A, B);
  }

  timer.reset();
  for (uint32_t i = 0; i < n; ++i) {
    auto A = Kokkos::subview(trmm_args.A, i, Kokkos::ALL(), Kokkos::ALL());
    auto B = Kokkos::subview(trmm_args.B, i, Kokkos::ALL(), Kokkos::ALL());

    SerialTrmm<side, uplo, trans, diag, tag>::invoke(trmm_args.alpha, A, B);
  }
  Kokkos::fence();
  __trmm_output_csv_row(options, trmm_args, timer.seconds());
#else
  std::cerr << std::string(__func__)
            << " disabled since KOKKOS_ENABLE_CUDA is defined." << std::endl;
#endif  // !KOKKOS_ENABLE_CUDA
}

template <class scalar_type, class vta, class vtb, class device_type>
void __do_trmm_serial_batched(options_t options, trmm_args_t trmm_args) {
  char __side = tolower(trmm_args.side), __uplo = tolower(trmm_args.uplo),
       __trans = tolower(trmm_args.trans);
  //__diag = tolower(diag[0]);

  STATUS;

  //// Lower non-transpose ////
  if (__side == 'l' && __uplo == 'l' && __trans == 'n') {
    STATUS;
    __do_trmm_serial_batched_template<Side::Left, Uplo::Lower,
                                      Trans::NoTranspose, Diag::Unit>(
        options, trmm_args);
  }
  if (__side == 'r' && __uplo == 'l' && __trans == 'n') {
    STATUS;
    __do_trmm_serial_batched_template<Side::Right, Uplo::Lower,
                                      Trans::NoTranspose, Diag::Unit>(
        options, trmm_args);
  }
  //// Lower transpose /////
  // Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (__side == 'l' && __uplo == 'l' && __trans == 't') {
    STATUS;
    __do_trmm_serial_batched_template<Side::Left, Uplo::Lower, Trans::Transpose,
                                      Diag::Unit>(options, trmm_args);
  }
  if (__side == 'r' && __uplo == 'l' && __trans == 't') {
    STATUS;
    __do_trmm_serial_batched_template<Side::Right, Uplo::Lower,
                                      Trans::Transpose, Diag::Unit>(options,
                                                                    trmm_args);
  }
  //// Lower conjugate-transpose ////
  // Conjugate-Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (__side == 'l' && __uplo == 'l' && __trans == 'c') {
    STATUS;
    __do_trmm_serial_batched_template<Side::Left, Uplo::Lower,
                                      Trans::ConjTranspose, Diag::Unit>(
        options, trmm_args);
  }
  if (__side == 'r' && __uplo == 'l' && __trans == 'c') {
    STATUS;
    __do_trmm_serial_batched_template<Side::Right, Uplo::Lower,
                                      Trans::ConjTranspose, Diag::Unit>(
        options, trmm_args);
  }
  //// Upper non-transpose ////
  if (__side == 'l' && __uplo == 'u' && __trans == 'n') {
    STATUS;
    __do_trmm_serial_batched_template<Side::Left, Uplo::Upper,
                                      Trans::NoTranspose, Diag::Unit>(
        options, trmm_args);
  }
  if (__side == 'r' && __uplo == 'u' && __trans == 'n') {
    STATUS;
    __do_trmm_serial_batched_template<Side::Right, Uplo::Upper,
                                      Trans::NoTranspose, Diag::Unit>(
        options, trmm_args);
  }
  //// Upper transpose
  // Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (__side == 'l' && __uplo == 'u' && __trans == 't') {
    STATUS;
    __do_trmm_serial_batched_template<Side::Left, Uplo::Upper, Trans::Transpose,
                                      Diag::Unit>(options, trmm_args);
  }
  if (__side == 'r' && __uplo == 'u' && __trans == 't') {
    STATUS;
    __do_trmm_serial_batched_template<Side::Right, Uplo::Upper,
                                      Trans::Transpose, Diag::Unit>(options,
                                                                    trmm_args);
  }

  //// Upper conjugate-transpose ////
  // Conjugate-Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (__side == 'l' && __uplo == 'u' && __trans == 'c') {
    STATUS;
    __do_trmm_serial_batched_template<Side::Left, Uplo::Upper,
                                      Trans::ConjTranspose, Diag::Unit>(
        options, trmm_args);
  }
  if (__side == 'r' && __uplo == 'u' && __trans == 'c') {
    STATUS;
    __do_trmm_serial_batched_template<Side::Right, Uplo::Upper,
                                      Trans::ConjTranspose, Diag::Unit>(
        options, trmm_args);
  }

  return;
}

#if !defined(KOKKOS_ENABLE_CUDA)
template <class ExecutionSpace>
struct parallel_blas_trmm {
  trmm_args_t trmm_args_;

  parallel_blas_trmm(trmm_args_t trmm_args) : trmm_args_(trmm_args) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i) const {
    auto svA = Kokkos::subview(trmm_args_.A, i, Kokkos::ALL(), Kokkos::ALL());
    auto svB = Kokkos::subview(trmm_args_.B, i, Kokkos::ALL(), Kokkos::ALL());

    KokkosBlas::trmm(&trmm_args_.side, &trmm_args_.uplo, &trmm_args_.trans,
                     &trmm_args_.diag, trmm_args_.alpha, svA, svB);
  }
};
#endif  // !KOKKOS_ENABLE_CUDA

template <class scalar_type, class vta, class vtb, class device_type>
void __do_trmm_parallel_blas(options_t options, trmm_args_t trmm_args) {
#if !defined(KOKKOS_ENABLE_CUDA)
  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  Kokkos::Timer timer;
  using execution_space = typename device_type::execution_space;
  using functor_type    = parallel_blas_trmm<execution_space>;
  functor_type parallel_blas_trmm_functor(trmm_args);

  STATUS;

  Kokkos::parallel_for("parallelBlasWarmUpLoopTrmm",
                       Kokkos::RangePolicy<execution_space>(0, warm_up_n),
                       parallel_blas_trmm_functor);
  Kokkos::fence();

  timer.reset();
  Kokkos::parallel_for("parallelBlasTimedLoopTrmm",
                       Kokkos::RangePolicy<execution_space>(0, n),
                       parallel_blas_trmm_functor);
  Kokkos::fence();
  __trmm_output_csv_row(options, trmm_args, timer.seconds());
#else
  std::cerr << std::string(__func__)
            << " disabled since KOKKOS_ENABLE_CUDA is defined." << std::endl;
  __trmm_output_csv_row(options, trmm_args, -1);
#endif  // !KOKKOS_ENABLE_CUDA
  return;
}

template <class side, class uplo, class trans, class diag, class tag,
          class ExecutionSpace>
struct parallel_batched_trmm {
  trmm_args_t trmm_args_;

  parallel_batched_trmm(trmm_args_t trmm_args) : trmm_args_(trmm_args) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i) const {
    auto svA = Kokkos::subview(trmm_args_.A, i, Kokkos::ALL(), Kokkos::ALL());
    auto svB = Kokkos::subview(trmm_args_.B, i, Kokkos::ALL(), Kokkos::ALL());

    SerialTrmm<side, uplo, trans, diag, tag>::invoke(trmm_args_.alpha, svA,
                                                     svB);
  }
};

template <class side, class uplo, class trans, class diag, class device_type>
void __do_trmm_parallel_batched_template(options_t options,
                                         trmm_args_t trmm_args) {
  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  Kokkos::Timer timer;
  using tag             = Algo::Trmm::Unblocked;
  using execution_space = typename device_type::execution_space;
  using functor_type =
      parallel_batched_trmm<side, uplo, trans, diag, tag, execution_space>;
  functor_type parallel_batched_trmm_functor(trmm_args);

  STATUS;

  Kokkos::parallel_for("parallelBatchedWarmUpLoopTrmm",
                       Kokkos::RangePolicy<execution_space>(0, warm_up_n),
                       parallel_batched_trmm_functor);
  Kokkos::fence();

  timer.reset();
  Kokkos::parallel_for("parallelBatchedTimedLoopTrmm",
                       Kokkos::RangePolicy<execution_space>(0, n),
                       parallel_batched_trmm_functor);
  Kokkos::fence();
  __trmm_output_csv_row(options, trmm_args, timer.seconds());

  return;
}

template <class scalar_type, class vta, class vtb, class device_type>
void __do_trmm_parallel_batched(options_t options, trmm_args_t trmm_args) {
  char __side = tolower(trmm_args.side), __uplo = tolower(trmm_args.uplo),
       __trans = tolower(trmm_args.trans);
  //__diag = tolower(diag[0]);

  STATUS;

  //// Lower non-transpose ////
  if (__side == 'l' && __uplo == 'l' && __trans == 'n') {
    STATUS;
    __do_trmm_parallel_batched_template<
        Side::Left, Uplo::Lower, Trans::NoTranspose, Diag::Unit, device_type>(
        options, trmm_args);
  }
  if (__side == 'r' && __uplo == 'l' && __trans == 'n') {
    STATUS;
    __do_trmm_parallel_batched_template<
        Side::Right, Uplo::Lower, Trans::NoTranspose, Diag::Unit, device_type>(
        options, trmm_args);
  }
  //// Lower transpose /////
  // Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (__side == 'l' && __uplo == 'l' && __trans == 't') {
    STATUS;
    __do_trmm_parallel_batched_template<
        Side::Left, Uplo::Lower, Trans::Transpose, Diag::Unit, device_type>(
        options, trmm_args);
  }
  if (__side == 'r' && __uplo == 'l' && __trans == 't') {
    STATUS;
    __do_trmm_parallel_batched_template<
        Side::Right, Uplo::Lower, Trans::Transpose, Diag::Unit, device_type>(
        options, trmm_args);
  }
  //// Lower conjugate-transpose ////
  // Conjugate-Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (__side == 'l' && __uplo == 'l' && __trans == 'c') {
    STATUS;
    __do_trmm_parallel_batched_template<
        Side::Left, Uplo::Lower, Trans::ConjTranspose, Diag::Unit, device_type>(
        options, trmm_args);
  }
  if (__side == 'r' && __uplo == 'l' && __trans == 'c') {
    STATUS;
    __do_trmm_parallel_batched_template<Side::Right, Uplo::Lower,
                                        Trans::ConjTranspose, Diag::Unit,
                                        device_type>(options, trmm_args);
  }
  //// Upper non-transpose ////
  if (__side == 'l' && __uplo == 'u' && __trans == 'n') {
    STATUS;
    __do_trmm_parallel_batched_template<
        Side::Left, Uplo::Upper, Trans::NoTranspose, Diag::Unit, device_type>(
        options, trmm_args);
  }
  if (__side == 'r' && __uplo == 'u' && __trans == 'n') {
    STATUS;
    __do_trmm_parallel_batched_template<
        Side::Right, Uplo::Upper, Trans::NoTranspose, Diag::Unit, device_type>(
        options, trmm_args);
  }
  //// Upper transpose
  // Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (__side == 'l' && __uplo == 'u' && __trans == 't') {
    STATUS;
    __do_trmm_parallel_batched_template<
        Side::Left, Uplo::Upper, Trans::Transpose, Diag::Unit, device_type>(
        options, trmm_args);
  }
  if (__side == 'r' && __uplo == 'u' && __trans == 't') {
    STATUS;
    __do_trmm_parallel_batched_template<
        Side::Right, Uplo::Upper, Trans::Transpose, Diag::Unit, device_type>(
        options, trmm_args);
  }

  //// Upper conjugate-transpose ////
  // Conjugate-Transpose A by simply swapping the dimensions (extent) and stride
  // parameters
  if (__side == 'l' && __uplo == 'u' && __trans == 'c') {
    STATUS;
    __do_trmm_parallel_batched_template<
        Side::Left, Uplo::Upper, Trans::ConjTranspose, Diag::Unit, device_type>(
        options, trmm_args);
  }
  if (__side == 'r' && __uplo == 'u' && __trans == 'c') {
    STATUS;
    __do_trmm_parallel_batched_template<Side::Right, Uplo::Upper,
                                        Trans::ConjTranspose, Diag::Unit,
                                        device_type>(options, trmm_args);
  }

  return;
}

/*************************** Internal setup fns **************************/
template <class scalar_type, class vta, class vtb, class device_type>
trmm_args_t __do_setup(options_t options, matrix_dims_t dim) {
  using execution_space = typename device_type::execution_space;

  trmm_args_t trmm_args;
  uint64_t seed = Kokkos::Impl::clock_tic();
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(seed);
  decltype(dim.a.m) min_dim = dim.a.m < dim.a.n ? dim.a.m : dim.a.n;
  typename vta::HostMirror host_A;
  STATUS;

  trmm_args.side  = options.blas_args.trmm.trmm_args.c_str()[0];
  trmm_args.uplo  = options.blas_args.trmm.trmm_args.c_str()[1];
  trmm_args.trans = options.blas_args.trmm.trmm_args.c_str()[2];
  trmm_args.diag  = options.blas_args.trmm.trmm_args.c_str()[3];
  trmm_args.A     = vta("trmm_args.A", options.n, dim.a.m, dim.a.n);
  trmm_args.B     = vtb("trmm_args.B", options.n, dim.b.m, dim.b.n);
  trmm_args.alpha = options.blas_args.trmm.alpha;
  host_A          = Kokkos::create_mirror_view(trmm_args.A);

  Kokkos::fill_random(trmm_args.A, rand_pool,
                      Kokkos::rand<Kokkos::Random_XorShift64<execution_space>,
                                   scalar_type>::max());
  Kokkos::deep_copy(host_A, trmm_args.A);

  if (trmm_args.uplo == 'U' || trmm_args.uplo == 'u') {
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

  if (trmm_args.diag == 'U' || trmm_args.diag == 'u') {
    for (uint32_t k = 0; k < options.n; ++k) {
      auto A = Kokkos::subview(host_A, k, Kokkos::ALL(), Kokkos::ALL());
      for (int i = 0; i < min_dim; i++) {
        A(i, i) = scalar_type(1);
      }
    }
  }
  Kokkos::deep_copy(trmm_args.A, host_A);

  Kokkos::fill_random(trmm_args.B, rand_pool,
                      Kokkos::rand<Kokkos::Random_XorShift64<execution_space>,
                                   scalar_type>::max());

  return trmm_args;
}

/*************************** Interal run helper fns **************************/
void __do_loop_and_invoke(options_t options,
                          void (*fn)(options_t, trmm_args_t)) {
  matrix_dims_t cur_dims;
  trmm_args_t trmm_args;
  STATUS;

  __print_trmm_perf_test_options(options);
  std::cout << "SCALAR:" << typeid(default_scalar).name()
            << ", LAYOUT:" << typeid(default_layout).name()
            << ", DEVICE:" << typeid(default_device).name() << std::endl;

  options.out[0] << trmm_csv_header_str << std::endl;

  for (cur_dims = options.start;
       cur_dims.a.m <= options.stop.a.m && cur_dims.a.n <= options.stop.a.n &&
       cur_dims.b.m <= options.stop.b.m && cur_dims.b.n <= options.stop.b.n;
       cur_dims.a.m *= options.step, cur_dims.a.n *= options.step,
      cur_dims.b.m *= options.step, cur_dims.b.n *= options.step) {
    trmm_args =
        __do_setup<default_scalar, view_type_3d, view_type_3d, default_device>(
            options, cur_dims);
    fn(options, trmm_args);
  }
  return;
}

/*************************** External fns **************************/
void do_trmm_serial_blas(options_t options) {
  STATUS;
  __do_loop_and_invoke(
      options, __do_trmm_serial_blas<default_scalar, view_type_3d, view_type_3d,
                                     default_device>);
  return;
}

void do_trmm_serial_batched(options_t options) {
  STATUS;
  __do_loop_and_invoke(options,
                       __do_trmm_serial_batched<default_scalar, view_type_3d,
                                                view_type_3d, default_device>);
  return;
}

void do_trmm_parallel_blas(options_t options) {
  STATUS;
  __do_loop_and_invoke(options,
                       __do_trmm_parallel_blas<default_scalar, view_type_3d,
                                               view_type_3d, default_device>);
  return;
}

void do_trmm_parallel_batched(options_t options) {
  STATUS;
  __do_loop_and_invoke(
      options, __do_trmm_parallel_batched<default_scalar, view_type_3d,
                                          view_type_3d, default_device>);
  return;
}

#endif  // KOKKOSBLAS3_TRMM_PERF_TEST_H_
