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
#include <KokkosLapack_trtri.hpp>
#include <KokkosKernels_TestUtils.hpp>

#include <chrono>

namespace Test {

template <class ViewTypeA, class ExecutionSpace>
struct UnitDiagTRTRI {
  ViewTypeA A_;
  using ScalarA = typename ViewTypeA::value_type;

  UnitDiagTRTRI(const ViewTypeA& A) : A_(A) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i) const { A_(i, i) = ScalarA(1); }
};
template <class ViewTypeA, class ExecutionSpace>
struct NonUnitDiagTRTRI {
  ViewTypeA A_;
  using ScalarA = typename ViewTypeA::value_type;

  NonUnitDiagTRTRI(const ViewTypeA& A) : A_(A) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i) const { A_(i, i) = A_(i, i) + 10; }
};
template <class ViewTypeA, class ViewTypeB, class ViewTypeC, class ExecutionSpace>
struct VanillaGEMM {
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

template <class ViewTypeA, class Device>
int impl_test_trtri(int bad_diag_idx, const char* uplo, const char* diag, const int M, const int N) {
  using execution_space = typename ViewTypeA::device_type::execution_space;
  using ScalarA         = typename ViewTypeA::value_type;
  using APT             = Kokkos::ArithTraits<ScalarA>;
  using mag_type        = typename APT::mag_type;

  double machine_eps         = APT::epsilon();
  const mag_type eps         = 1.0e8 * machine_eps;  //~1e-13 for double
  bool is_A_lower_triangular = (uplo[0] == 'L') || (uplo[0] == 'l');
  int ret;
  ViewTypeA A("A", M, N);
  ViewTypeA A_original("A_original", M, N);
  ViewTypeA A_I("A_I", M, N);  // is I taken...?
  uint64_t seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  ScalarA beta  = ScalarA(0);
  ScalarA cur_check_val;  // Either 1 or 0, to check A_I

  // const int As0 = A.stride(0), As1 = A.stride(1);
  // const int Ae0 = A.extent(0), Ae1 = A.extent(1);
  // printf("KokkosLapack::trtri test for %c %c, M %d, N %d, eps %g, ViewType:
  // %s, A.stride(0): %d, A.stride(1): %d, A.extent(0): %d, A.extent(1): %d
  // START\n", uplo[0],diag[0],M,N,eps,typeid(ViewTypeA).name(), As0, As1, Ae0,
  // Ae1); fflush(stdout);

  typename ViewTypeA::HostMirror host_A = Kokkos::create_mirror_view(A);
  typename ViewTypeA::HostMirror host_I = Kokkos::create_mirror_view(A);

  if (M != N || bad_diag_idx > 0) {
    if (bad_diag_idx > 0) {
      for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
          if (i == j)
            host_A(i, j) = ScalarA(1);
          else
            host_A(i, j) = ScalarA(0);
        }
      }
      // Set just 1 value in the diagonal to 0.
      if (M > 0 && N > 0) host_A(bad_diag_idx - 1, bad_diag_idx - 1) = ScalarA(0);
      Kokkos::deep_copy(A, host_A);
    }
    return KokkosLapack::trtri(uplo, diag, A);
  }

  // If M is greater than 100 and A is an unit triangluar matrix, make A the
  // identity matrix due to large rounding errors in unit matrices
  bool M_gt_100 = (M > 100) && ((diag[0] == 'U') || (diag[0] == 'u'));

  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(seed);

  // Initialize A with deterministic random numbers
  Kokkos::fill_random(A, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, ScalarA>::max());
  if ((diag[0] == 'U') || (diag[0] == 'u')) {
    using functor_type = UnitDiagTRTRI<ViewTypeA, execution_space>;
    functor_type udtrtri(A);
    // Initialize As diag with 1s
    Kokkos::parallel_for("KokkosLapack::Test::UnitDiagTRTRI", Kokkos::RangePolicy<execution_space>(0, M), udtrtri);
  } else {  //(diag[0]=='N')||(diag[0]=='n')
    using functor_type = NonUnitDiagTRTRI<ViewTypeA, execution_space>;
    functor_type nudtrtri(A);
    // Initialize As diag with A(i,i)+10
    Kokkos::parallel_for("KokkosLapack::Test::NonUnitDiagTRTRI", Kokkos::RangePolicy<execution_space>(0, M), nudtrtri);
  }
  Kokkos::fence();
  Kokkos::deep_copy(host_A, A);

  // Make host_A a lower triangle
  if (is_A_lower_triangular || M_gt_100) {
    for (int i = 0; i < M - 1; i++)
      for (int j = i + 1; j < N; j++) host_A(i, j) = ScalarA(0);
  }
  if (!is_A_lower_triangular || M_gt_100) {
    // Make host_A a upper triangle
    for (int i = 1; i < M; i++)
      for (int j = 0; j < i; j++) host_A(i, j) = ScalarA(0);
  }
  Kokkos::deep_copy(A, host_A);
  Kokkos::deep_copy(A_original, A);

#if 0
    Kokkos::deep_copy(host_A, A);
    printf("host_A:\n");
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
          printf("%*.13lf ", 20, host_A(i,j));
        }
        printf("\n");
    }
#endif

  // A = A^-1
  ret = KokkosLapack::trtri(uplo, diag, A);
  Kokkos::fence();

  if (ret) {
    printf("KokkosLapack::trtri(%c, %c, %s) returned %d\n", uplo[0], diag[0], typeid(ViewTypeA).name(), ret);
    return ret;
  }

#if 0
    Kokkos::deep_copy(host_A, A);
    printf("host_A:\n");
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
          printf("%*.13lf ", 20, host_A(i,j));
        }
        printf("\n");
    }
#endif

  // A_I = A * A_original
  struct VanillaGEMM<ViewTypeA, ViewTypeA, ViewTypeA, execution_space> vgemm;
  vgemm.A_t   = false;
  vgemm.B_t   = false;
  vgemm.A_c   = false;
  vgemm.B_c   = false;
  vgemm.N     = N;
  vgemm.K     = M;
  vgemm.A     = A;
  vgemm.B     = A_original;
  vgemm.C     = A_I;  // out
  vgemm.alpha = ScalarA(1);
  vgemm.beta  = beta;
  Kokkos::parallel_for("KokkosLapack::Test::VanillaGEMM",
                       Kokkos::TeamPolicy<execution_space>(
                           M, Kokkos::AUTO, KokkosKernels::Impl::kk_get_max_vector_size<execution_space>()),
                       vgemm);
  Kokkos::fence();
  Kokkos::deep_copy(host_I, A_I);

#if 0
    printf("host_I:\n");
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
          printf("%*.13lf ", 20, host_I(i,j));
        }
        printf("\n");
    }
#endif

  bool test_flag = true;
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      // Set check value
      cur_check_val = (i == j) ? ScalarA(1) : ScalarA(0);  // APT::abs(host_A(i,j));

      // Check how close |A_I - cur_check_val| is to 0.
      if (APT::abs(APT::abs(host_I(i, j)) - cur_check_val) > eps) {
        test_flag = false;
        // printf("   Error: eps ( %g ), host_I ( %.15f ) != cur_check_val (
        // %.15f ) (abs result-cur_check_val %g) at (i %d, j %d)\n", eps,
        // host_I(i,j), cur_check_val, APT::abs(host_I(i,j) - cur_check_val), i,
        // j);
        break;
      }
    }
    if (!test_flag) break;
  }
  EXPECT_EQ(test_flag, true);
  return ret;
}
}  // namespace Test

template <class ScalarA, class Device>
int test_trtri(const char* mode) {
  int ret;
  int bad_diag_idx = -1;
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  using view_type_a_layout_left = Kokkos::View<ScalarA**, Kokkos::LayoutLeft, Device>;

  ret = Test::impl_test_trtri<view_type_a_layout_left, Device>(bad_diag_idx, &mode[0], &mode[1], 0, 0);
  EXPECT_EQ(ret, 0);

  ret = Test::impl_test_trtri<view_type_a_layout_left, Device>(bad_diag_idx, &mode[0], &mode[1], 1, 1);
  EXPECT_EQ(ret, 0);

  ret = Test::impl_test_trtri<view_type_a_layout_left, Device>(bad_diag_idx, &mode[0], &mode[1], 15, 15);
  EXPECT_EQ(ret, 0);

  ret = Test::impl_test_trtri<view_type_a_layout_left, Device>(bad_diag_idx, &mode[0], &mode[1], 100, 100);
  EXPECT_EQ(ret, 0);

  // Rounding errors with randomly generated matrices begin here where M>100, so
  // we pass in A=I
  ret = Test::impl_test_trtri<view_type_a_layout_left, Device>(bad_diag_idx, &mode[0], &mode[1], 273, 273);
  EXPECT_EQ(ret, 0);

  // Only non-unit matrices could be singular.
  if (mode[1] == 'N' || mode[1] == 'n') {
    bad_diag_idx = 2;  // 1-index based
    ret          = Test::impl_test_trtri<view_type_a_layout_left, Device>(bad_diag_idx, &mode[0], &mode[1], 2, 2);
    EXPECT_EQ(ret, bad_diag_idx);
    bad_diag_idx = -1;
  }

  // One time check, disabled due to runtime throw instead of return here
  // ret = Test::impl_test_trtri<view_type_a_layout_left,
  // Device>(&mode[0],&mode[1],1031,731); EXPECT_NE(ret, 0);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  using view_type_a_layout_right = Kokkos::View<ScalarA**, Kokkos::LayoutRight, Device>;

  ret = Test::impl_test_trtri<view_type_a_layout_right, Device>(bad_diag_idx, &mode[0], &mode[1], 0, 0);
  EXPECT_EQ(ret, 0);

  ret = Test::impl_test_trtri<view_type_a_layout_right, Device>(bad_diag_idx, &mode[0], &mode[1], 1, 1);
  EXPECT_EQ(ret, 0);

  ret = Test::impl_test_trtri<view_type_a_layout_right, Device>(bad_diag_idx, &mode[0], &mode[1], 15, 15);
  EXPECT_EQ(ret, 0);

  ret = Test::impl_test_trtri<view_type_a_layout_right, Device>(bad_diag_idx, &mode[0], &mode[1], 100, 100);
  EXPECT_EQ(ret, 0);

  // Rounding errors with randomly generated matrices begin here where M>100, so
  // we pass in A=I
  ret = Test::impl_test_trtri<view_type_a_layout_right, Device>(bad_diag_idx, &mode[0], &mode[1], 273, 273);
  EXPECT_EQ(ret, 0);

  // Only non-unit matrices could be singular.
  if (mode[1] == 'N' || mode[1] == 'n') {
    bad_diag_idx = 2;  // 1-index based
    ret          = Test::impl_test_trtri<view_type_a_layout_right, Device>(bad_diag_idx, &mode[0], &mode[1], 2, 2);
    EXPECT_EQ(ret, bad_diag_idx);
    bad_diag_idx = -1;
  }
#endif

  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, trtri_float) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::trtri_float");
  test_trtri<float, TestDevice>("UN");
  test_trtri<float, TestDevice>("UU");
  test_trtri<float, TestDevice>("LN");
  test_trtri<float, TestDevice>("LU");
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, trtri_double) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::trtri_double");
  test_trtri<double, TestDevice>("UN");
  test_trtri<double, TestDevice>("UU");
  test_trtri<double, TestDevice>("LN");
  test_trtri<double, TestDevice>("LU");
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, trtri_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::trtri_complex_double");
  test_trtri<Kokkos::complex<double>, TestDevice>("UN");
  test_trtri<Kokkos::complex<double>, TestDevice>("UU");
  test_trtri<Kokkos::complex<double>, TestDevice>("LN");
  test_trtri<Kokkos::complex<double>, TestDevice>("LU");
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, trtri_complex_float) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::trtri_complex_float");
  test_trtri<Kokkos::complex<float>, TestDevice>("UN");
  test_trtri<Kokkos::complex<float>, TestDevice>("UU");
  test_trtri<Kokkos::complex<float>, TestDevice>("LN");
  test_trtri<Kokkos::complex<float>, TestDevice>("LU");
  Kokkos::Profiling::popRegion();
}
#endif
