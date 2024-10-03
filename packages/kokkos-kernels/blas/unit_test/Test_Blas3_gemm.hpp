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
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosBlas3_gemm.hpp>
#include <KokkosKernels_TestUtils.hpp>

#include <chrono>

namespace Test {

template <class ViewTypeA, class ViewTypeB, class ViewTypeC, class ExecutionSpace>
struct gemm_VanillaGEMM {
  bool A_t, B_t, A_c, B_c;
  int N, K;
  ViewTypeA A;
  ViewTypeB B;
  ViewTypeC C;

  typedef typename ViewTypeA::value_type ScalarA;
  typedef typename ViewTypeB::value_type ScalarB;
  typedef typename ViewTypeC::value_type ScalarC;
  typedef Kokkos::ArithTraits<ScalarC> APT;
  typedef typename APT::mag_type mag_type;
  ScalarA alpha;
  ScalarC beta;

  KOKKOS_INLINE_FUNCTION
  void operator()(const typename Kokkos::TeamPolicy<ExecutionSpace>::member_type& team) const {
// GNU COMPILER BUG WORKAROUND
#if defined(KOKKOS_COMPILER_GNU) && !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
    int i = team.league_rank();
#else
    const int i = team.league_rank();
#endif
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [&](const int& j) {
      ScalarC C_ij = 0.0;

      // GNU 5.3, 5.4 and 6.1 (and maybe more) crash with another nested lambda
      // here

#if defined(KOKKOS_COMPILER_GNU) && !defined(KOKKOS_COMPILER_NVCC)
      for (int k = 0; k < K; k++) {
        ScalarA A_ik = A_t ? (A_c ? APT::conj(A(k, i)) : A(k, i)) : A(i, k);
        ScalarB B_kj = B_t ? (B_c ? APT::conj(B(j, k)) : B(j, k)) : B(k, j);
        C_ij += A_ik * B_kj;
      }
#else
        Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(team,K), [&] (const int& k, ScalarC& lsum) {
           ScalarA A_ik = A_t?(A_c?APT::conj(A(k,i)):A(k,i)):A(i,k);
           ScalarB B_kj = B_t?(B_c?APT::conj(B(j,k)):B(j,k)):B(k,j);
           lsum += A_ik*B_kj;
        },C_ij);
#endif

      C(i, j) = beta * C(i, j) + alpha * C_ij;
    });
  }
};

template <class ViewTypeA, class ViewTypeB, class ViewTypeC>
void build_matrices(const int M, const int N, const int K, const typename ViewTypeA::value_type alpha, ViewTypeA& A,
                    ViewTypeB& B, const typename ViewTypeA::value_type beta, ViewTypeC& C, ViewTypeC& Cref) {
  using execution_space = typename TestDevice::execution_space;
  using ScalarA         = typename ViewTypeA::non_const_value_type;
  using ScalarB         = typename ViewTypeB::non_const_value_type;
  using ScalarC         = typename ViewTypeC::non_const_value_type;

  A    = ViewTypeA("A", M, K);
  B    = ViewTypeB("B", K, N);
  C    = ViewTypeC("C", M, N);
  Cref = ViewTypeC("Cref", M, N);

  // (SA 11 Dec 2019) Max (previously: 10) increased to detect the bug in
  // Trilinos issue #6418
  const uint64_t seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(seed);
  Kokkos::fill_random(
      A, rand_pool,
      Kokkos::rand<typename Kokkos::Random_XorShift64_Pool<execution_space>::generator_type, ScalarA>::max());
  Kokkos::fill_random(
      B, rand_pool,
      Kokkos::rand<typename Kokkos::Random_XorShift64_Pool<execution_space>::generator_type, ScalarB>::max());
  Kokkos::fill_random(
      C, rand_pool,
      Kokkos::rand<typename Kokkos::Random_XorShift64_Pool<execution_space>::generator_type, ScalarC>::max());

  Kokkos::deep_copy(Cref, C);
  Kokkos::fence();

  struct Test::gemm_VanillaGEMM<ViewTypeA, ViewTypeB, ViewTypeC, execution_space> vgemm;
  vgemm.A_t   = false;
  vgemm.B_t   = false;
  vgemm.A_c   = false;
  vgemm.B_c   = false;
  vgemm.N     = N;
  vgemm.K     = K;
  vgemm.A     = A;
  vgemm.B     = B;
  vgemm.C     = Cref;
  vgemm.alpha = alpha;
  vgemm.beta  = beta;

  Kokkos::parallel_for("KokkosBlas::Test::gemm_VanillaGEMM",
                       Kokkos::TeamPolicy<execution_space>(
                           M, Kokkos::AUTO, KokkosKernels::Impl::kk_get_max_vector_size<execution_space>()),
                       vgemm);
  Kokkos::fence();
}

template <class ViewTypeC, class ExecutionSpace>
struct DiffGEMM {
  int N;
  ViewTypeC C, C2;

  typedef typename ViewTypeC::value_type ScalarC;
  typedef Kokkos::ArithTraits<ScalarC> APT;
  typedef typename APT::mag_type mag_type;

  KOKKOS_INLINE_FUNCTION
  void operator()(const typename Kokkos::TeamPolicy<ExecutionSpace>::member_type& team, mag_type& diff) const {
    const int i       = team.league_rank();
    mag_type diff_row = 0;
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, N),
        [&](const int& j, mag_type& diff_ij) {
          // printf("A (%i %i) (%i %i) (%i
          // %i)\n",C.extent(0),C.extent(1),C2.extent(0),C2.extent(1),i,j);
          diff_ij += APT::abs(C(i, j) - C2(i, j));
          // printf("B (%i %i) (%i %i) (%i
          // %i)\n",C.extent(0),C.extent(1),C2.extent(0),C2.extent(1),i,j);
        },
        diff_row);
    Kokkos::single(Kokkos::PerTeam(team), [&]() { diff += diff_row; });
  }
};

template <class ViewTypeA, class ViewTypeB, class ViewTypeC, class Device>
void impl_test_gemm(const char* TA, const char* TB, int M, int N, int K, typename ViewTypeA::value_type alpha,
                    typename ViewTypeC::value_type beta) {
  bool A_t = (TA[0] != 'N') && (TA[0] != 'n');
  bool B_t = (TB[0] != 'N') && (TB[0] != 'n');
  bool A_c = (TA[0] == 'C') || (TA[0] == 'c');
  bool B_c = (TB[0] == 'C') || (TB[0] == 'c');
  typedef typename ViewTypeA::device_type::execution_space execution_space;
  typedef typename ViewTypeA::value_type ScalarA;
  typedef typename ViewTypeB::value_type ScalarB;
  typedef typename ViewTypeC::value_type ScalarC;
  typedef Kokkos::ArithTraits<ScalarC> APT;
  typedef typename APT::mag_type mag_type;

  double machine_eps = APT::epsilon();

  ViewTypeA A("A", A_t ? K : M, A_t ? M : K);
  ViewTypeB B("B", B_t ? N : K, B_t ? K : N);
  ViewTypeC C("C", M, N);
  ViewTypeC C2("C", M, N);

  const uint64_t seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(seed);

  // (SA 11 Dec 2019) Max (previously: 10) increased to detect the bug in
  // Trilinos issue #6418
  Kokkos::fill_random(
      A, rand_pool,
      Kokkos::rand<typename Kokkos::Random_XorShift64_Pool<execution_space>::generator_type, ScalarA>::max());
  Kokkos::fill_random(
      B, rand_pool,
      Kokkos::rand<typename Kokkos::Random_XorShift64_Pool<execution_space>::generator_type, ScalarB>::max());
  Kokkos::fill_random(
      C, rand_pool,
      Kokkos::rand<typename Kokkos::Random_XorShift64_Pool<execution_space>::generator_type, ScalarC>::max());

  Kokkos::deep_copy(C2, C);
  Kokkos::fence();

  struct gemm_VanillaGEMM<ViewTypeA, ViewTypeB, ViewTypeC, execution_space> vgemm;
  vgemm.A_t   = A_t;
  vgemm.B_t   = B_t;
  vgemm.A_c   = A_c;
  vgemm.B_c   = B_c;
  vgemm.N     = N;
  vgemm.K     = K;
  vgemm.A     = A;
  vgemm.B     = B;
  vgemm.C     = C2;
  vgemm.alpha = alpha;
  vgemm.beta  = beta;

  Kokkos::parallel_for("KokkosBlas::Test::gemm_VanillaGEMM",
                       Kokkos::TeamPolicy<execution_space>(
                           M, Kokkos::AUTO, KokkosKernels::Impl::kk_get_max_vector_size<execution_space>()),
                       vgemm);

  KokkosBlas::gemm(TA, TB, alpha, A, B, beta, C);

  mag_type diff_C = 0;
  struct DiffGEMM<ViewTypeC, execution_space> diffgemm;
  diffgemm.N  = N;
  diffgemm.C  = C;
  diffgemm.C2 = C2;

  Kokkos::parallel_reduce("KokkosBlas::Test::DiffGEMM", Kokkos::TeamPolicy<execution_space>(M, Kokkos::AUTO), diffgemm,
                          diff_C);

  if (N != 0 && M != 0) {
    int K_eff             = (K == 0) ? 1 : K;
    double diff_C_average = diff_C / (N * M);
    // Expected Result: Random Walk in the least significant bit (i.e. ~
    // sqrt(K)*eps eps scales with the total sum and has a factor in it for the
    // accuracy of the operations -> eps = K * 75 * machine_eps * 7
    double diff_C_expected = 1.0 * sqrt(K_eff) * K_eff * 75 * machine_eps * 7;

    if ((diff_C_average >= 1.05 * diff_C_expected)) {
      printf("Result: %e %e\n", diff_C_average, diff_C_expected);
    }
    EXPECT_TRUE((diff_C_average < 1.05 * diff_C_expected));
  }
}

template <typename Scalar, typename Layout, typename Device>
void impl_test_stream_gemm_psge2(const int M, const int N, const int K, const Scalar alpha, const Scalar beta) {
  using execution_space = typename Device::execution_space;
  using ViewTypeA       = Kokkos::View<Scalar**, Layout, Device>;
  using ViewTypeB       = Kokkos::View<Scalar**, Layout, Device>;
  using ViewTypeC       = Kokkos::View<Scalar**, Layout, Device>;
  using ScalarC         = typename ViewTypeC::value_type;
  using APT             = Kokkos::ArithTraits<ScalarC>;
  using mag_type        = typename APT::mag_type;

  const char tA[]          = {"N"};
  const char tB[]          = {"N"};
  const double machine_eps = APT::epsilon();

  ViewTypeA A1, A2;
  ViewTypeB B1, B2;
  ViewTypeC C1, C1ref, C2, C2ref;

  Test::build_matrices(M, N, K, alpha, A1, B1, beta, C1, C1ref);
  Test::build_matrices(N, M, K, alpha, A2, B2, beta, C2, C2ref);

  auto instances = Kokkos::Experimental::partition_space(execution_space(), 1, 1);
  KokkosBlas::gemm(instances[0], tA, tB, alpha, A1, B1, beta, C1);
  KokkosBlas::gemm(instances[1], tA, tB, alpha, A2, B2, beta, C2);
  Kokkos::fence();

  mag_type diff_C1 = 0;
  struct Test::DiffGEMM<ViewTypeC, execution_space> diffgemm1;
  diffgemm1.N  = N;
  diffgemm1.C  = C1;
  diffgemm1.C2 = C1ref;

  Kokkos::parallel_reduce("KokkosBlas::Test::DiffGEMM1",
                          Kokkos::TeamPolicy<execution_space>(
                              M, Kokkos::AUTO, KokkosKernels::Impl::kk_get_max_vector_size<execution_space>()),
                          diffgemm1, diff_C1);

  mag_type diff_C2 = 0;
  struct Test::DiffGEMM<ViewTypeC, execution_space> diffgemm2;
  diffgemm2.N  = M;
  diffgemm2.C  = C2;
  diffgemm2.C2 = C2ref;

  Kokkos::parallel_reduce("KokkosBlas::Test::DiffGEMM2",
                          Kokkos::TeamPolicy<execution_space>(
                              N, Kokkos::AUTO, KokkosKernels::Impl::kk_get_max_vector_size<execution_space>()),
                          diffgemm2, diff_C2);
  Kokkos::fence();

  if (N != 0 && M != 0) {
    int K_eff = (K == 0) ? 1 : K;
    // Expected Result: Random Walk in the least significant bit (i.e. ~
    // sqrt(K)*eps eps scales with the total sum and has a factor in it for the
    // accuracy of the operations -> eps = K * 75 * machine_eps * 7
    const double diff_C_expected = 1.0 * sqrt(K_eff) * K_eff * 75 * machine_eps * 7;

    const double diff_C1_average = diff_C1 / (N * M);
    if ((diff_C1_average >= 1.05 * diff_C_expected)) {
      printf("Result: %e %e\n", diff_C1_average, diff_C_expected);
    }
    EXPECT_TRUE((diff_C1_average < 1.05 * diff_C_expected));

    const double diff_C2_average = diff_C2 / (N * M);
    if ((diff_C2_average >= 1.05 * diff_C_expected)) {
      printf("Result: %e %e\n", diff_C2_average, diff_C_expected);
    }
    EXPECT_TRUE((diff_C2_average < 1.05 * diff_C_expected));
  }
}
}  // namespace Test

template <typename Scalar, typename Layout>
void test_gemm() {
  typedef typename TestDevice::execution_space execution_space;
  typedef Kokkos::View<Scalar**, Layout, TestDevice> view_type_a;
  typedef Kokkos::View<Scalar**, Layout, TestDevice> view_type_b;
  typedef Kokkos::View<Scalar**, Layout, TestDevice> view_type_c;
  std::vector<const char*> modes = {"N", "T"};
  if (std::is_same<Scalar, Kokkos::complex<float>>::value || std::is_same<Scalar, Kokkos::complex<double>>::value)
    modes.push_back("C");
  Scalar alpha              = 4.5;
  std::vector<Scalar> betas = {0.0, 3.0};
  for (Scalar beta : betas) {
    for (auto amode : modes) {
      for (auto bmode : modes) {
        Test::impl_test_gemm<view_type_a, view_type_b, view_type_c, TestDevice>(amode, bmode, 0, 0, 0, alpha, beta);
        // BMK: N = 1 exercises the special GEMV code path in GEMM (currently,
        // only for modes N/N)
        Test::impl_test_gemm<view_type_a, view_type_b, view_type_c, TestDevice>(amode, bmode, 50, 1, 40, alpha, beta);
        // LBV: K = 0 exercise the quick return code path in GEMM
        Test::impl_test_gemm<view_type_a, view_type_b, view_type_c, TestDevice>(amode, bmode, 20, 14, 0, alpha, beta);
        Test::impl_test_gemm<view_type_a, view_type_b, view_type_c, TestDevice>(amode, bmode, 13, 15, 17, alpha, beta);
        Test::impl_test_gemm<view_type_a, view_type_b, view_type_c, TestDevice>(amode, bmode, 179, 15, 211, alpha,
                                                                                beta);
        Test::impl_test_gemm<view_type_a, view_type_b, view_type_c, TestDevice>(amode, bmode, 12, 3071, 517, alpha,
                                                                                beta);
      }
    }
  }
  auto pool_size = execution_space().concurrency();
  if (pool_size >= 2) {
    Test::impl_test_stream_gemm_psge2<Scalar, Layout, TestDevice>(53, 42, 17, 4.5,
                                                                  3.0);                  // General code path
    Test::impl_test_stream_gemm_psge2<Scalar, Layout, TestDevice>(13, 1, 17, 4.5, 3.0);  // gemv based gemm code path
    Test::impl_test_stream_gemm_psge2<Scalar, Layout, TestDevice>(7, 13, 17, 4.5,
                                                                  3.0);  // dot based gemm code path
  }
}

template <typename Scalar>
void test_gemm_enabled_layouts() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  test_gemm<Scalar, Kokkos::LayoutLeft>();
#endif
#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  test_gemm<Scalar, Kokkos::LayoutRight>();
#endif
}

template <class scalar1, class scalar2>
void test_gemm_mixed_scalars() {
  using Matrix1 = Kokkos::View<scalar1**, TestDevice>;
  using Matrix2 = Kokkos::View<scalar2**, TestDevice>;

  const int dim1 = 400, dim2 = 1000;

  Matrix1 A("A", dim1, dim1);
  Matrix1 B("B", dim2, dim2);
  Matrix1 C("C", dim2, dim1);
  Matrix2 D("D", dim2, dim1);

  Kokkos::deep_copy(A, Kokkos::ArithTraits<scalar1>::one());
  Kokkos::deep_copy(B, Kokkos::ArithTraits<scalar1>::one());
  Kokkos::deep_copy(C, Kokkos::ArithTraits<scalar1>::one());

  KokkosBlas::gemm(TestDevice(), "N", "N", 1.0, D, A, 0.0, C);
  KokkosBlas::gemm(TestDevice(), "N", "T", 1.0, C, D, 0.0, B);
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, gemm_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemm_float");
  test_gemm_enabled_layouts<float>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, gemm_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemm_double");
  test_gemm_enabled_layouts<double>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, gemm_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemm_complex_double");
  test_gemm_enabled_layouts<Kokkos::complex<double>>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, gemm_complex_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemm_complex_float");
  test_gemm_enabled_layouts<Kokkos::complex<float>>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) && !defined(KOKKOSKERNELS_ETI_ONLY)
TEST_F(TestCategory, gemm_mixed_scalars_complex_double_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemm_mixed_complex_double_double");
  test_gemm_mixed_scalars<Kokkos::complex<double>, double>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT) && !defined(KOKKOSKERNELS_ETI_ONLY)
TEST_F(TestCategory, gemm_mixed_scalar_complex_float_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemm_mixed_complex_float_float");
  test_gemm_mixed_scalars<Kokkos::complex<float>, float>();
  Kokkos::Profiling::popRegion();
}
#endif
