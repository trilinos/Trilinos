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
#ifndef KOKKOSBLAS3_GEMM_PERF_TEST_H_
#define KOKKOSBLAS3_GEMM_PERF_TEST_H_

//#include <complex.h>
#include "KokkosBlas3_common.hpp"

#include <Kokkos_Random.hpp>

#include <KokkosBlas3_gemm.hpp>

#include "KokkosBatched_Gemm_Decl.hpp"
#include "KokkosBatched_Gemm_Serial_Impl.hpp"
//#include "KokkosBatched_Gemm_Team_Impl.hpp"
//#include "KokkosBatched_Gemm_TeamVector_Impl.hpp"
#include "KokkosBatched_Util.hpp"

//#define GEMM_PERF_TEST_DEBUG

// Forward declarations
void do_gemm_serial_blas(options_t options);
void do_gemm_serial_batched(options_t options);
void do_gemm_serial_batched_blocked(options_t options);
// void do_gemm_experiment(options_t options);

// void do_gemm_serial_blas_parallel(options_t options);
// Not valid! The KokkosBlas::gemm function may take the entire device per
// invocation!
void do_gemm_serial_batched_parallel(options_t options);
void do_gemm_serial_batched_blocked_parallel(options_t options);
void do_gemm_team_batched_parallel(options_t options);
void do_gemm_team_batched_blocked_parallel(options_t options);
void do_gemm_team_vector_batched_parallel(options_t options);
void do_gemm_team_vector_batched_blocked_parallel(options_t options);
void do_gemm_experiment_parallel(options_t options);

struct SerialTag {};
struct TeamTag {};
struct TeamVectorTag {};
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
        NULL  // Serial Experiment
    },
    {
        NULL,  // BLAS
        do_gemm_serial_batched_parallel,
        do_gemm_serial_batched_blocked_parallel,  // Serial
        do_gemm_team_batched_parallel,
        do_gemm_team_batched_blocked_parallel,       // Team
        do_gemm_team_vector_batched_parallel, NULL,  // TeamVector
        do_gemm_experiment_parallel                  // Parallel Experiment
    }};

/*************************** Test types and defaults **************************/
#define DEFAULT_GEMM_ARGS "NN"
#define DEFAULT_GEMM_ALPHA 1.0

using view_type_3d =
    Kokkos::View<default_scalar ***, default_layout, default_device>;

struct batched_params {
  int team_size;
  int vector_len;
};
typedef struct batched_params batched_params_t;

struct gemm_args {
  char transA, transB;
  default_scalar alpha;
  default_scalar beta;
  view_type_3d A, B, C;
  batched_params_t bp;
};
typedef struct gemm_args gemm_args_t;

static std::string gemm_csv_header_str =
    "algorithm,transAtransB,alpha,beta,team_size,vector_len,loop_type,A_dims,B_"
    "dims,C_dims,warm_up_n,"
    "iter,total_time(s),average_time(s)";

/*************************** Internal helper fns **************************/
static void __gemm_output_csv_row(options_t options, gemm_args_t gemm_args,
                                  double time_in_seconds,
                                  const char *experiment_name = nullptr) {
  std::string algo_name = test_e_str[options.test];
  if (experiment_name) algo_name = std::string(experiment_name);

  options.out[0] << algo_name << "," << options.blas_args.gemm.gemm_args << ","
                 << options.blas_args.gemm.alpha << ","
                 << options.blas_args.gemm.beta << "," << gemm_args.bp.team_size
                 << "," << gemm_args.bp.vector_len << ","
                 << loop_e_str[options.loop] << "," << gemm_args.A.extent(0)
                 << "x" << gemm_args.A.extent(1) << "x" << gemm_args.A.extent(2)
                 << "," << gemm_args.B.extent(0) << "x" << gemm_args.B.extent(1)
                 << "x" << gemm_args.B.extent(2) << "," << gemm_args.C.extent(0)
                 << "x" << gemm_args.C.extent(1) << "x" << gemm_args.C.extent(2)
                 << "," << options.warm_up_n << "," << options.n << ","
                 << time_in_seconds << "," << time_in_seconds / options.n
                 << std::endl;
}

static void __print_gemm_perf_test_options(options_t options) {
#ifdef PERF_TEST_DEBUG
  printf("options.test      = %s\n", test_e_str[options.test].c_str());
  printf("options.loop      = %s\n", loop_e_str[options.loop].c_str());
  printf("options.start     = %dx%d,%dx%d\n", options.start.a.m,
         options.start.a.n, options.start.b.m, options.start.b.n);
  printf("options.stop      = %dx%d,%dx%d\n", options.stop.a.m,
         options.stop.a.n, options.stop.b.m, options.stop.b.n);
  printf("options.step      = %d\n", options.step);
  printf("options.warm_up_n = %d\n", options.warm_up_n);
  printf("options.n         = %d\n", options.n);
  printf("options.blas_args.gemm.gemm_args = %s\n",
         options.blas_args.gemm.gemm_args.c_str());
  printf("options.out_file  = %s\n", options.out_file.c_str());
  if (std::is_same<double, default_scalar>::value)
    printf("options.alpha     = %lf\n", options.blas_args.gemm.alpha);
  else if (std::is_same<float, default_scalar>::value)
    printf("options.alpha     = %f\n", options.blas_args.gemm.alpha);
#endif  // PERF_TEST_DEBUG
  return;
}

/*************************** Internal templated fns **************************/
template <class scalar_type, class vta, class vtb, class device_type>
void __do_gemm_serial_blas(options_t options, gemm_args_t gemm_args) {
// Need to take subviews on the device
#if !defined(KOKKOS_ENABLE_CUDA)
  Kokkos::Timer timer;

  STATUS;

  auto __do_loop = [](uint32_t n, gemm_args_t _gemm_args) {
    for (uint32_t i = 0; i < n; ++i) {
      auto A = Kokkos::subview(_gemm_args.A, i, Kokkos::ALL(), Kokkos::ALL());
      auto B = Kokkos::subview(_gemm_args.B, i, Kokkos::ALL(), Kokkos::ALL());
      auto C = Kokkos::subview(_gemm_args.C, i, Kokkos::ALL(), Kokkos::ALL());

      KokkosBlas::gemm(&_gemm_args.transA, &_gemm_args.transB, _gemm_args.alpha,
                       A, B, _gemm_args.beta, C);
    }
  };
  __do_loop(options.warm_up_n, gemm_args);
  Kokkos::fence();

  timer.reset();
  __do_loop(options.n, gemm_args);
  Kokkos::fence();

  __gemm_output_csv_row(options, gemm_args, timer.seconds());
#else
  std::cerr << std::string(__func__)
            << " disabled since KOKKOS_ENABLE_CUDA is defined." << std::endl;
#endif  // !KOKKOS_ENABLE_CUDA
  return;
}

template <class TransAType, class TransBType, class AlgoType>
void __do_gemm_serial_batched_template(options_t options,
                                       gemm_args_t gemm_args) {
// Need to take subviews on the device
#if !defined(KOKKOS_ENABLE_CUDA)
  Kokkos::Timer timer;

  auto __do_loop = [](uint32_t n, gemm_args_t _gemm_args) {
    for (uint32_t i = 0; i < n; ++i) {
      auto A = Kokkos::subview(_gemm_args.A, i, Kokkos::ALL(), Kokkos::ALL());
      auto B = Kokkos::subview(_gemm_args.B, i, Kokkos::ALL(), Kokkos::ALL());
      auto C = Kokkos::subview(_gemm_args.C, i, Kokkos::ALL(), Kokkos::ALL());

      SerialGemm<TransAType, TransBType, AlgoType>::invoke(
          _gemm_args.alpha, A, B, _gemm_args.beta, C);
    }
  };

  __do_loop(options.warm_up_n, gemm_args);
  Kokkos::fence();

  timer.reset();
  __do_loop(options.n, gemm_args);
  Kokkos::fence();
  __gemm_output_csv_row(options, gemm_args, timer.seconds());
#else
  std::cerr << std::string(__func__)
            << " disabled since KOKKOS_ENABLE_CUDA is defined." << std::endl;
#endif  // !KOKKOS_ENABLE_CUDA
}

template <class scalar_type, class vta, class vtb, class vtc, class device_type,
          class algo_type>
void __do_gemm_serial_batched(options_t options, gemm_args_t gemm_args) {
  char a  = gemm_args.transA;
  char b  = gemm_args.transB;
  using N = Trans::NoTranspose;
  using T = Trans::Transpose;
  // using C = Trans::ConjTranspose;

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

#if !defined(KOKKOS_ENABLE_CUDA)
template <class ExecutionSpace>
struct parallel_blas_gemm {
  gemm_args_t gemm_args_;

  parallel_blas_gemm(gemm_args_t gemm_args) : gemm_args_(gemm_args) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int &i) const {
    auto svA = Kokkos::subview(gemm_args_.A, i, Kokkos::ALL(), Kokkos::ALL());
    auto svB = Kokkos::subview(gemm_args_.B, i, Kokkos::ALL(), Kokkos::ALL());
    auto svC = Kokkos::subview(gemm_args_.C, i, Kokkos::ALL(), Kokkos::ALL());

    KokkosBlas::gemm(&gemm_args_.transA, &gemm_args_.transB, gemm_args_.alpha,
                     svA, svB, gemm_args_.beta, svC);
  }
};
#endif  // !KOKKOS_ENABLE_CUDA

template <class scalar_type, class vta, class vtb, class device_type>
void __do_gemm_parallel_blas(options_t options, gemm_args_t gemm_args) {
#if !defined(KOKKOS_ENABLE_CUDA) && !defined(KOKKOS_ENABLE_HIP)
  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  Kokkos::Timer timer;
  using execution_space = typename device_type::execution_space;
  using functor_type    = parallel_blas_gemm<execution_space>;
  functor_type parallel_blas_gemm_functor(gemm_args);

  STATUS;

  Kokkos::parallel_for("parallelBlasWarmUpLoopGemm",
                       Kokkos::RangePolicy<execution_space>(0, warm_up_n),
                       parallel_blas_gemm_functor);
  Kokkos::fence();

  timer.reset();
  Kokkos::parallel_for("parallelBlasTimedLoopGemm",
                       Kokkos::RangePolicy<execution_space>(0, n),
                       parallel_blas_gemm_functor);
  Kokkos::fence();
  __gemm_output_csv_row(options, gemm_args, timer.seconds());
#else
  std::cerr << std::string(__func__)
            << " disabled since KOKKOS_ENABLE_CUDA is defined." << std::endl;
  __gemm_output_csv_row(options, gemm_args, -1);
#endif  // !KOKKOS_ENABLE_CUDA
  return;
}

template <class MemberType, class TransAType, class TransBType,
          class BlockingType>
struct parallel_batched_gemm {
  gemm_args_t gemm_args_;

  parallel_batched_gemm(gemm_args_t gemm_args) : gemm_args_(gemm_args) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const SerialTag &, const MemberType &member) const {
    auto i   = member.league_rank();
    auto svA = Kokkos::subview(gemm_args_.A, i, Kokkos::ALL(), Kokkos::ALL());
    auto svB = Kokkos::subview(gemm_args_.B, i, Kokkos::ALL(), Kokkos::ALL());
    auto svC = Kokkos::subview(gemm_args_.C, i, Kokkos::ALL(), Kokkos::ALL());

    KokkosBatched::SerialGemm<TransAType, TransBType, BlockingType>::invoke(
        gemm_args_.alpha, svA, svB, gemm_args_.beta, svC);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TeamTag &, const MemberType &member) const {
    auto i   = member.league_rank();
    auto svA = Kokkos::subview(gemm_args_.A, i, Kokkos::ALL(), Kokkos::ALL());
    auto svB = Kokkos::subview(gemm_args_.B, i, Kokkos::ALL(), Kokkos::ALL());
    auto svC = Kokkos::subview(gemm_args_.C, i, Kokkos::ALL(), Kokkos::ALL());

    KokkosBatched::TeamGemm<MemberType, TransAType, TransBType,
                            BlockingType>::invoke(member, gemm_args_.alpha, svA,
                                                  svB, gemm_args_.beta, svC);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TeamVectorTag &, const MemberType &member) const {
    auto team_idx = member.league_rank();
    auto svA =
        Kokkos::subview(gemm_args_.A, team_idx, Kokkos::ALL(), Kokkos::ALL());
    auto svB =
        Kokkos::subview(gemm_args_.B, team_idx, Kokkos::ALL(), Kokkos::ALL());
    auto svC =
        Kokkos::subview(gemm_args_.C, team_idx, Kokkos::ALL(), Kokkos::ALL());

    KokkosBatched::TeamVectorGemm<MemberType, TransAType, TransBType,
                                  BlockingType>::invoke(member,
                                                        gemm_args_.alpha, svA,
                                                        svB, gemm_args_.beta,
                                                        svC);
  }
};

template <class TransAType, class TransBType, class BlockingType, class AlgoTag,
          class device_type>
void __do_gemm_parallel_batched_template(options_t options,
                                         gemm_args_t gemm_args) {
  using execution_space = typename device_type::execution_space;
  using policy_type     = Kokkos::TeamPolicy<AlgoTag, execution_space>;
  using member_type     = typename policy_type::member_type;
  using functor_type =
      parallel_batched_gemm<member_type, TransAType, TransBType, BlockingType>;

  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  auto league_size   = options.start.c.k;
  Kokkos::Timer timer;

  STATUS;

  functor_type parallel_batched_gemm_functor(gemm_args);
  auto team_size  = gemm_args.bp.team_size;
  auto vector_len = gemm_args.bp.vector_len;

  for (uint32_t i = 0; i < warm_up_n; i++) {
    Kokkos::parallel_for("parallelBatchedWarmUpLoopGemm",
                         policy_type(league_size, team_size, vector_len),
                         parallel_batched_gemm_functor);
  }
  Kokkos::fence();

  timer.reset();
  for (uint32_t i = 0; i < n; i++) {
    Kokkos::parallel_for("parallelBatchedTimedLoopGemm",
                         policy_type(league_size, team_size, vector_len),
                         parallel_batched_gemm_functor);
  }
  Kokkos::fence();

  __gemm_output_csv_row(options, gemm_args, timer.seconds());

  return;
}

template <class algo_tag, class blocking_type, class device_type>
void __do_gemm_parallel_batched(options_t options, gemm_args_t gemm_args) {
  char a  = gemm_args.transA;
  char b  = gemm_args.transB;
  using N = Trans::NoTranspose;
  using T = Trans::Transpose;
  // using C = Trans::ConjTranspose;

  STATUS;

  if (a == 'N' && b == 'N') {
    __do_gemm_parallel_batched_template<N, N, blocking_type, algo_tag,
                                        device_type>(options, gemm_args);
  } else if (a == 'N' && b == 'T') {
    __do_gemm_parallel_batched_template<N, T, blocking_type, algo_tag,
                                        device_type>(options, gemm_args);
    //} else if (a == 'N' && b == 'C') {
    //  __do_gemm_parallel_batched_template<N, C, blocking_type, algo_tag,
    //  device_type>(options, gemm_args);
  } else if (a == 'T' && b == 'N') {
    __do_gemm_parallel_batched_template<T, N, blocking_type, algo_tag,
                                        device_type>(options, gemm_args);
  } else if (a == 'T' && b == 'T') {
    __do_gemm_parallel_batched_template<T, T, blocking_type, algo_tag,
                                        device_type>(options, gemm_args);
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

  parallel_batched_gemm_experiment1(gemm_args_t gemm_args)
      : gemm_args_(gemm_args) {}

  KOKKOS_INLINE_FUNCTION

  void operator()(const SerialTag &, const int &i) const {
    auto svA = Kokkos::subview(gemm_args_.A, i, Kokkos::ALL(), Kokkos::ALL());
    auto svB = Kokkos::subview(gemm_args_.B, i, Kokkos::ALL(), Kokkos::ALL());
    auto svC = Kokkos::subview(gemm_args_.C, i, Kokkos::ALL(), Kokkos::ALL());

    // Uses two serial for-loops internally
    KokkosBatched::SerialGemm<TransAType, TransBType, BlockingType>::invoke(
        gemm_args_.alpha, svA, svB, gemm_args_.beta, svC);
  }
};

/**
 * 1. parallel_for(rangePolicy<Kokkos::DefaultExecutionSpace>(N)): serialGemm
 *
 */
template <class TransAType, class TransBType, class BlockingType,
          class device_type>
void __do_gemm_parallel_experiment1(options_t options, gemm_args_t gemm_args) {
  using execution_space = typename device_type::execution_space;
  using policy_type     = Kokkos::RangePolicy<SerialTag, execution_space>;
  using functor_type =
      parallel_batched_gemm_experiment1<TransAType, TransBType, BlockingType>;

  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  auto k             = options.start.c.k;
  Kokkos::Timer timer;
  STATUS;

  functor_type experiment1_functor(gemm_args);

  for (uint32_t i = 0; i < warm_up_n; ++i) {
    Kokkos::parallel_for("parallelBatchedUntimedExperiment1Gemm",
                         policy_type(0, k), experiment1_functor);
  }
  Kokkos::fence();

  timer.reset();
  for (uint32_t i = 0; i < n; ++i) {
    Kokkos::parallel_for("parallelBatchedTimedExperiment1Gemm",
                         policy_type(0, k), experiment1_functor);
  }
  Kokkos::fence();

  __gemm_output_csv_row(options, gemm_args, timer.seconds(), "experiment1");
  return;
}

template <class TransAType, class TransBType, class BlockingType,
          class MemberType>
struct parallel_batched_gemm_experiment2_3_4 {
  gemm_args_t gemm_args_;

  parallel_batched_gemm_experiment2_3_4(gemm_args_t gemm_args)
      : gemm_args_(gemm_args) {}

  // Experiment 2
  KOKKOS_INLINE_FUNCTION
  void operator()(const TeamVectorTag &, const MemberType &member) const {
    auto i   = member.league_rank();
    auto svA = Kokkos::subview(gemm_args_.A, i, Kokkos::ALL(), Kokkos::ALL());
    auto svB = Kokkos::subview(gemm_args_.B, i, Kokkos::ALL(), Kokkos::ALL());
    auto svC = Kokkos::subview(gemm_args_.C, i, Kokkos::ALL(), Kokkos::ALL());

    // Uses TeamThreadRange over C-rows
    //        ThreadVectorRange over C-cols
    KokkosBatched::TeamVectorGemm<MemberType, TransAType, TransBType,
                                  BlockingType>::invoke(member,
                                                        gemm_args_.alpha, svA,
                                                        svB, gemm_args_.beta,
                                                        svC);
  }

  // Experiment 3
  KOKKOS_INLINE_FUNCTION
  void operator()(const LayoutLeftTag &, const MemberType &member) const {
    auto team_idx = member.league_rank();
    auto svA =
        Kokkos::subview(gemm_args_.A, team_idx, Kokkos::ALL(), Kokkos::ALL());
    auto svB =
        Kokkos::subview(gemm_args_.B, team_idx, Kokkos::ALL(), Kokkos::ALL());
    auto svC =
        Kokkos::subview(gemm_args_.C, team_idx, Kokkos::ALL(), Kokkos::ALL());

    // TeamThreadRange:   splits the index range over the threads of the team
    // ThreadVectorRange: splits the index range over the vector lanes of the
    // calling thread

    auto svC_cols = svC.extent(1);
    // In a given team, for each vector lane, compute zero or more output
    // columns of C depending on the index range
    Kokkos::parallel_for(
        Kokkos::ThreadVectorRange(member, svC_cols), [&](const int &lane_idx) {
          auto svB_col = Kokkos::subview(svB, Kokkos::ALL(), lane_idx);
          auto svC_col = Kokkos::subview(svC, Kokkos::ALL(), lane_idx);
          // TeamGemm Calls TeamThreadRange over M*N meaning the flat M*N array
          // is split over all threads of the team
          KokkosBatched::TeamGemm<MemberType, TransAType, TransBType,
                                  BlockingType>::invoke(member,
                                                        gemm_args_.alpha, svA,
                                                        svB_col,
                                                        gemm_args_.beta,
                                                        svC_col);
        });
  }

  // TODO: Why is this faster than the LayoutLeftTag operator above for both
  // LayoutLeft and LayoutRight? Experiment 4
  KOKKOS_INLINE_FUNCTION
  void operator()(const LayoutRightTag &, const MemberType &member) const {
    auto team_idx = member.league_rank();
    auto svA =
        Kokkos::subview(gemm_args_.A, team_idx, Kokkos::ALL(), Kokkos::ALL());
    auto svB =
        Kokkos::subview(gemm_args_.B, team_idx, Kokkos::ALL(), Kokkos::ALL());
    auto svC =
        Kokkos::subview(gemm_args_.C, team_idx, Kokkos::ALL(), Kokkos::ALL());

    // TeamThreadRange:   splits the index range over the threads of the team
    // ThreadVectorRange: splits the index range over the vector lanes of the
    // calling thread

    auto svC_rows = svC.extent(0);
    // In a given team, for each vector lane, compute zero or more output rows
    // of C depending on the index range
    Kokkos::parallel_for(
        Kokkos::ThreadVectorRange(member, svC_rows), [&](const int &lane_idx) {
          auto svA_row = Kokkos::subview(svA, lane_idx, Kokkos::ALL());
          auto svC_row = Kokkos::subview(svC, lane_idx, Kokkos::ALL());
          // TeamGemm Calls TeamThreadRange over M*N meaning the flat M*N array
          // is split over all threads of the team
          KokkosBatched::TeamGemm<MemberType, TransAType, TransBType,
                                  BlockingType>::invoke(member,
                                                        gemm_args_.alpha,
                                                        svA_row, svB,
                                                        gemm_args_.beta,
                                                        svC_row);
        });
  }
};

/**
 * 2. case a)
 * parallel_for(teamPolicy): TeamVectorGemm
 *
 */
template <class TransAType, class TransBType, class BlockingType,
          class device_type>
void __do_gemm_parallel_experiment2(options_t options, gemm_args_t gemm_args) {
  using execution_space = typename device_type::execution_space;
  using policy_type     = Kokkos::TeamPolicy<TeamVectorTag, execution_space>;
  using member_type     = typename policy_type::member_type;
  using functor_type =
      parallel_batched_gemm_experiment2_3_4<TransAType, TransBType,
                                            BlockingType, member_type>;

  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  auto league_size   = options.start.c.k;
  Kokkos::Timer timer;
  STATUS;

  functor_type experiment2_functor(gemm_args);

  auto team_size  = gemm_args.bp.team_size;
  auto vector_len = gemm_args.bp.vector_len;

  for (uint32_t i = 0; i < warm_up_n; ++i) {
    Kokkos::parallel_for("parallelBatchedUntimedExperiment2Gemm",
                         policy_type(league_size, team_size, vector_len),
                         experiment2_functor);
  }
  Kokkos::fence();

  timer.reset();

  for (uint32_t i = 0; i < n; ++i) {
    Kokkos::parallel_for("parallelBatchedTimedExperiment2Gemm",
                         policy_type(league_size, team_size, vector_len),
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
template <class TransAType, class TransBType, class BlockingType,
          class device_type>
void __do_gemm_parallel_experiment3(options_t options, gemm_args_t gemm_args) {
  using execution_space = typename device_type::execution_space;
  // using layout_tag = std::conditional<std::is_same<default_layout,
  // Kokkos::LayoutLeft>::value, LayoutLeftTag, LayoutRightTag>::type;
  using policy_type = Kokkos::TeamPolicy<LayoutLeftTag, execution_space>;
  using member_type = typename policy_type::member_type;
  using functor_type =
      parallel_batched_gemm_experiment2_3_4<TransAType, TransBType,
                                            BlockingType, member_type>;

  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  auto league_size   = options.start.c.k;
  Kokkos::Timer timer;
  STATUS;

  functor_type experiment3_functor(gemm_args);

  auto team_size  = gemm_args.bp.team_size;
  auto vector_len = gemm_args.bp.vector_len;

  for (uint32_t i = 0; i < warm_up_n; ++i) {
    Kokkos::parallel_for("parallelBatchedUntimedExperiment3Gemm",
                         policy_type(league_size, team_size, vector_len),
                         experiment3_functor);
  }
  Kokkos::fence();

  timer.reset();

  for (uint32_t i = 0; i < n; ++i) {
    Kokkos::parallel_for("parallelBatchedTimedExperiment3Gemm",
                         policy_type(league_size, team_size, vector_len),
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
template <class TransAType, class TransBType, class BlockingType,
          class device_type>
void __do_gemm_parallel_experiment4(options_t options, gemm_args_t gemm_args) {
  using execution_space = typename device_type::execution_space;
  // using layout_tag = std::conditional<std::is_same<default_layout,
  // Kokkos::LayoutLeft>::value, LayoutLeftTag, LayoutRightTag>::type;
  using policy_type = Kokkos::TeamPolicy<LayoutRightTag, execution_space>;
  using member_type = typename policy_type::member_type;
  using functor_type =
      parallel_batched_gemm_experiment2_3_4<TransAType, TransBType,
                                            BlockingType, member_type>;

  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  auto league_size   = options.start.c.k;
  Kokkos::Timer timer;
  STATUS;

  functor_type experiment4_functor(gemm_args);

  auto team_size  = gemm_args.bp.team_size;
  auto vector_len = gemm_args.bp.vector_len;

  for (uint32_t i = 0; i < warm_up_n; ++i) {
    Kokkos::parallel_for("parallelBatchedUntimedExperiment4Gemm",
                         policy_type(league_size, team_size, vector_len),
                         experiment4_functor);
  }
  Kokkos::fence();

  timer.reset();

  for (uint32_t i = 0; i < n; ++i) {
    Kokkos::parallel_for("parallelBatchedTimedExperiment4Gemm",
                         policy_type(league_size, team_size, vector_len),
                         experiment4_functor);
  }
  Kokkos::fence();

  __gemm_output_csv_row(options, gemm_args, timer.seconds(), "experiment4");
  return;
}

template <class SimdViewType, class TransAType, class TransBType,
          class BlockingType>
class parallel_batched_gemm_experiment5 {
 private:
  SimdViewType &A, &B, &C;
  gemm_args_t gemm_args;

 public:
  parallel_batched_gemm_experiment5(SimdViewType &_A, SimdViewType &_B,
                                    SimdViewType &_C, gemm_args_t _gemm_args)
      : A(_A), B(_B), C(_C), gemm_args(_gemm_args) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const SimdCpuTag &, const int &i) const {
    auto svA = Kokkos::subview(A, i, Kokkos::ALL(), Kokkos::ALL());
    auto svB = Kokkos::subview(B, i, Kokkos::ALL(), Kokkos::ALL());
    auto svC = Kokkos::subview(C, i, Kokkos::ALL(), Kokkos::ALL());

    // Uses two serial for-loops internally
    KokkosBatched::SerialGemm<TransAType, TransBType, BlockingType>::invoke(
        gemm_args.alpha, svA, svB, gemm_args.beta, svC);
  }
};

/**
 * 5.
 * parallel_for(RangePolicy<Kokkos:DefaultHostExecutionSpace>(N/vl+(N%vl>0)>):
 * serialGemm
 *
 * Not portable to GPU
 */
template <class TransAType, class TransBType, class BlockingType,
          class device_type>
void __do_gemm_parallel_experiment5(options_t options, gemm_args_t gemm_args) {
#if !defined(KOKKOS_ENABLE_CUDA) && !defined(KOKKOS_ENABLE_HIP)
  using execution_space = typename device_type::execution_space;
  using policy_type     = Kokkos::RangePolicy<SimdCpuTag, execution_space>;

  // Construct the SimdType
  using scalar_type = typename view_type_3d::value_type;
  constexpr int vl =
      KokkosBatched::DefaultVectorLength<scalar_type, execution_space>::value;
  using simd_type = KokkosBatched::Vector<KokkosBatched::SIMD<scalar_type>, vl>;
  using simd_view_type =
      Kokkos::View<simd_type ***, default_layout, default_device>;
  using functor_type =
      parallel_batched_gemm_experiment5<simd_view_type, TransAType, TransBType,
                                        BlockingType>;

  uint32_t warm_up_n = options.warm_up_n;
  uint32_t n         = options.n;
  auto k             = options.start.c.k;
  Kokkos::Timer timer;
  auto simd_batch_size = k / vl + (k % vl > 0);
  STATUS;

  // Increases each array size by sizeof(scalar_type) * (vl-1) bytes!
  simd_view_type A("A", simd_batch_size, gemm_args.A.extent(0),
                   gemm_args.A.extent(1));
  simd_view_type B("B", simd_batch_size, gemm_args.B.extent(0),
                   gemm_args.B.extent(1));
  simd_view_type C("C", simd_batch_size, gemm_args.C.extent(0),
                   gemm_args.C.extent(1));

  // uint64_t seed = Kokkos::Impl::clock_tic();
  // Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(seed);
  // Kokkos::fill_random(A, rand_pool,
  // Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, simd_type>::max());
  // Kokkos::fill_random(B, rand_pool,
  // Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, simd_type>::max());
  // Kokkos::fill_random(C, rand_pool,
  // Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, simd_type>::max());
  // execution_space::fence();

  functor_type experiment5_functor(A, B, C, gemm_args);

  for (uint32_t i = 0; i < warm_up_n; ++i) {
    Kokkos::parallel_for("parallelBatchedUntimedExperiment5Gemm",
                         policy_type(0, simd_batch_size), experiment5_functor);
  }
  Kokkos::fence();

  timer.reset();

  for (uint32_t i = 0; i < n; ++i) {
    Kokkos::parallel_for("parallelBatchedTimedExperiment5Gemm",
                         policy_type(0, simd_batch_size), experiment5_functor);
  }
  Kokkos::fence();

  __gemm_output_csv_row(options, gemm_args, timer.seconds(), "experiment5");
#else
  std::cerr
      << std::string(__func__)
      << " disabled since KOKKOS_ENABLE_CUDA or KOKKOS_ENABLE_HIP is defined."
      << std::endl;
#endif  // !KOKKOS_ENABLE_CUDA || !KOKKOS_ENABLE_HIP
  return;
}

/*************************** Internal setup fns **************************/
template <class scalar_type, class vta, class vtb, class vtc, class device_type>
gemm_args_t __do_setup(options_t options, matrix_dims_t dim) {
  using execution_space = typename device_type::execution_space;

  gemm_args_t gemm_args;
  uint64_t seed = Kokkos::Impl::clock_tic();
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(seed);
  STATUS;

  gemm_args.transA        = options.blas_args.gemm.gemm_args.c_str()[0];
  gemm_args.transB        = options.blas_args.gemm.gemm_args.c_str()[1];
  gemm_args.A             = vta("gemm_args.A", dim.a.k, dim.a.m, dim.a.n);
  gemm_args.B             = vtb("gemm_args.B", dim.b.k, dim.b.m, dim.b.n);
  gemm_args.C             = vtc("gemm_args.C", dim.c.k, dim.c.m, dim.c.n);
  gemm_args.alpha         = options.blas_args.gemm.alpha;
  gemm_args.alpha         = options.blas_args.gemm.beta;
  gemm_args.bp.team_size  = options.blas_args.team_size;
  gemm_args.bp.vector_len = options.blas_args.vector_len;

  Kokkos::fill_random(gemm_args.A, rand_pool,
                      Kokkos::rand<Kokkos::Random_XorShift64<execution_space>,
                                   scalar_type>::max());
  Kokkos::fill_random(gemm_args.B, rand_pool,
                      Kokkos::rand<Kokkos::Random_XorShift64<execution_space>,
                                   scalar_type>::max());
  Kokkos::fill_random(gemm_args.C, rand_pool,
                      Kokkos::rand<Kokkos::Random_XorShift64<execution_space>,
                                   scalar_type>::max());

  return gemm_args;
}

/*************************** Interal run helper fns **************************/
void __do_loop_and_invoke(options_t options,
                          void (*fn)(options_t, gemm_args_t)) {
  matrix_dims_t cur_dims;
  gemm_args_t gemm_args;
  STATUS;

  __print_gemm_perf_test_options(options);
  std::cout << "SCALAR:" << typeid(default_scalar).name()
            << ", LAYOUT:" << typeid(default_layout).name()
            << ", DEVICE:" << typeid(default_device).name() << std::endl;

  options.out[0] << gemm_csv_header_str << std::endl;

  for (cur_dims = options.start;
       cur_dims.a.m <= options.stop.a.m && cur_dims.a.n <= options.stop.a.n &&
       cur_dims.b.m <= options.stop.b.m && cur_dims.b.n <= options.stop.b.n &&
       cur_dims.c.m <= options.stop.c.m && cur_dims.c.n <= options.stop.c.n;
       cur_dims.a.m *= options.step, cur_dims.a.n *= options.step,
      cur_dims.b.m *= options.step, cur_dims.b.n *= options.step,
      cur_dims.c.m *= options.step, cur_dims.c.n *= options.step) {
    gemm_args = __do_setup<default_scalar, view_type_3d, view_type_3d,
                           view_type_3d, default_device>(options, cur_dims);
    fn(options, gemm_args);
  }
  return;
}

/*************************** External fns **************************/
void do_gemm_serial_blas(options_t options) {
  STATUS;
  __do_loop_and_invoke(
      options, __do_gemm_serial_blas<default_scalar, view_type_3d, view_type_3d,
                                     default_device>);
  return;
}

void do_gemm_serial_batched(options_t options) {
  STATUS;
  __do_loop_and_invoke(
      options, __do_gemm_serial_batched<default_scalar, view_type_3d,
                                        view_type_3d, view_type_3d,
                                        default_device, Algo::Gemm::Unblocked>);
  return;
}

void do_gemm_serial_batched_blocked(options_t options) {
  STATUS;
  __do_loop_and_invoke(
      options, __do_gemm_serial_batched<default_scalar, view_type_3d,
                                        view_type_3d, view_type_3d,
                                        default_device, Algo::Gemm::Blocked>);
  return;
}

void do_gemm_serial_batched_parallel(options_t options) {
  STATUS;
  __do_loop_and_invoke(
      options, __do_gemm_parallel_batched<SerialTag, Algo::Gemm::Unblocked,
                                          default_device>);
  return;
}

void do_gemm_serial_batched_blocked_parallel(options_t options) {
  STATUS;
  __do_loop_and_invoke(
      options, __do_gemm_parallel_batched<SerialTag, Algo::Gemm::Blocked,
                                          default_device>);
  return;
}

void do_gemm_team_batched_parallel(options_t options) {
  STATUS;
  __do_loop_and_invoke(
      options, __do_gemm_parallel_batched<TeamTag, Algo::Gemm::Unblocked,
                                          default_device>);
  return;
}

void do_gemm_team_batched_blocked_parallel(options_t options) {
  STATUS;
  __do_loop_and_invoke(
      options,
      __do_gemm_parallel_batched<TeamTag, Algo::Gemm::Blocked, default_device>);
  return;
}

void do_gemm_team_vector_batched_parallel(options_t options) {
  STATUS;
  __do_loop_and_invoke(
      options, __do_gemm_parallel_batched<TeamVectorTag, Algo::Gemm::Unblocked,
                                          default_device>);
  return;
}

/* void do_gemm_team_vector_batched_blocked_parallel(options_t options) {
  STATUS;
  __do_loop_and_invoke(
      options, __do_gemm_parallel_batched<TeamVectorTag, Algo::Gemm::Blocked,
default_device>); return;
} */

void do_gemm_experiment_parallel(options_t options) {
  STATUS;
  using TransAType   = Trans::NoTranspose;
  using TransBType   = Trans::NoTranspose;
  using BlockingType = Algo::Gemm::Unblocked;

  __do_loop_and_invoke(
      options, __do_gemm_parallel_experiment1<TransAType, TransBType,
                                              BlockingType, default_device>);
  __do_loop_and_invoke(
      options, __do_gemm_parallel_experiment2<TransAType, TransBType,
                                              BlockingType, default_device>);
  __do_loop_and_invoke(
      options, __do_gemm_parallel_experiment3<TransAType, TransBType,
                                              BlockingType, default_device>);
  __do_loop_and_invoke(
      options, __do_gemm_parallel_experiment4<TransAType, TransBType,
                                              BlockingType, default_device>);
  __do_loop_and_invoke(
      options, __do_gemm_parallel_experiment5<TransAType, TransBType,
                                              BlockingType, default_device>);
}

#endif  // KOKKOSBLAS3_GEMM_PERF_TEST_H_
