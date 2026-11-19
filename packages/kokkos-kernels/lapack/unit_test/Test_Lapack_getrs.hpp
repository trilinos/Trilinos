// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef TEST_LAPACK_GETRS_HPP
#define TEST_LAPACK_GETRS_HPP

#if (defined(TEST_CUDA_LAPACK_CPP) && defined(KOKKOSKERNELS_ENABLE_TPL_CUSOLVER)) ||               \
    (defined(TEST_HIP_LAPACK_CPP) && defined(KOKKOSKERNELS_ENABLE_TPL_ROCSOLVER)) ||               \
    ((defined(KOKKOSKERNELS_ENABLE_TPL_LAPACK) || defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)) && \
     (defined(TEST_OPENMP_LAPACK_CPP) || defined(TEST_SERIAL_LAPACK_CPP) || defined(TEST_THREADS_LAPACK_CPP)))

#include <gtest/gtest.h>
#include <KokkosLapack_getrf.hpp>
#include <KokkosLapack_getrs.hpp>
#include <KokkosKernels_ArithTraits.hpp>
#include <KokkosKernels_TestUtils.hpp>

namespace Test {

template <class MatrixType>
void impl_test_getrs_sym() {
  using Device         = typename MatrixType::device_type;
  using IpivType       = Kokkos::View<int *, Device>;
  using InfoType       = Kokkos::View<int *, Device>;
  using scalar_type    = typename MatrixType::non_const_value_type;
  using ExecutionSpace = typename Device::execution_space;

  constexpr int m = 4, n = 4;

  const scalar_type zero  = KokkosKernels::ArithTraits<scalar_type>::zero();
  const scalar_type one   = KokkosKernels::ArithTraits<scalar_type>::one();
  const scalar_type two   = one + one;
  const scalar_type three = two + one;
  const scalar_type four  = two + two;

  MatrixType A("matrix A", m, n);
  auto A_h  = Kokkos::create_mirror_view(A);
  A_h(0, 0) = two;
  A_h(0, 1) = -one;
  A_h(0, 2) = zero;
  A_h(0, 3) = zero;
  A_h(1, 0) = -one;
  A_h(1, 1) = two;
  A_h(1, 2) = -one;
  A_h(1, 3) = zero;
  A_h(2, 0) = zero;
  A_h(2, 1) = -one;
  A_h(2, 2) = two;
  A_h(2, 3) = -one;
  A_h(3, 0) = zero;
  A_h(3, 1) = zero;
  A_h(3, 2) = -one;
  A_h(3, 3) = two;
  Kokkos::deep_copy(A, A_h);

  const int min_mn = Kokkos::min(A.extent(0), A.extent(1));
  IpivType ipiv("LU pivots", min_mn);
  InfoType info("LU info", 1);

  KokkosLapack::getrf(ExecutionSpace(), A, ipiv, info);
  Kokkos::fence();

  // Single rhs
  MatrixType b1("b1", n, 1);
  auto b1_h = Kokkos::create_mirror_view(b1);
  auto tol  = min_mn * m * n * KokkosKernels::ArithTraits<scalar_type>::eps();

  for (auto mode : {"N", "T", "C"}) {
    b1_h(0, 0) = one;
    b1_h(1, 0) = zero;
    b1_h(2, 0) = zero;
    b1_h(3, 0) = one;
    Kokkos::deep_copy(b1, b1_h);

    KokkosLapack::getrs(ExecutionSpace(), mode, A, ipiv, b1, info);
    Kokkos::fence();

    Kokkos::deep_copy(b1_h, b1);
    for (int idx = 0; idx < b1.extent_int(0); ++idx) {
      EXPECT_NEAR_KK_REL(b1_h(idx, 0), one, tol);
    }
  }

  // Three rhs
  MatrixType b3("b3", n, 3);
  auto b3_h = Kokkos::create_mirror_view(b3);

  for (auto mode : {"N", "T", "C"}) {
    b3_h(0, 0) = one;
    b3_h(0, 1) = zero;
    b3_h(0, 2) = -two;
    b3_h(1, 0) = zero;
    b3_h(1, 1) = zero;
    b3_h(1, 2) = -two;
    b3_h(2, 0) = zero;
    b3_h(2, 1) = zero;
    b3_h(2, 2) = -two;
    b3_h(3, 0) = one;
    b3_h(3, 1) = 5 * one;
    b3_h(3, 2) = 23 * one;
    Kokkos::deep_copy(b3, b3_h);

    KokkosLapack::getrs(ExecutionSpace(), mode, A, ipiv, b3, info);
    Kokkos::fence();

    Kokkos::deep_copy(b3_h, b3);
    EXPECT_NEAR_KK_REL(b3_h(0, 0), one, tol);
    EXPECT_NEAR_KK_REL(b3_h(1, 0), one, tol);
    EXPECT_NEAR_KK_REL(b3_h(2, 0), one, tol);
    EXPECT_NEAR_KK_REL(b3_h(3, 0), one, tol);

    EXPECT_NEAR_KK_REL(b3_h(0, 1), one, tol);
    EXPECT_NEAR_KK_REL(b3_h(1, 1), two, tol);
    EXPECT_NEAR_KK_REL(b3_h(2, 1), three, tol);
    EXPECT_NEAR_KK_REL(b3_h(3, 1), four, tol);

    EXPECT_NEAR_KK_REL(b3_h(0, 2), one * one, tol);
    EXPECT_NEAR_KK_REL(b3_h(1, 2), two * two, tol);
    EXPECT_NEAR_KK_REL(b3_h(2, 2), three * three, tol);
    EXPECT_NEAR_KK_REL(b3_h(3, 2), four * four, tol);
  }
}

template <class MatrixType>
void impl_test_getrs_unsym() {
  using Device         = typename MatrixType::device_type;
  using IpivType       = Kokkos::View<int *, Device>;
  using InfoType       = Kokkos::View<int *, Device>;
  using scalar_type    = typename MatrixType::non_const_value_type;
  using ExecutionSpace = typename Device::execution_space;

  constexpr int n = 4;

  const scalar_type zero  = KokkosKernels::ArithTraits<scalar_type>::zero();
  const scalar_type one   = KokkosKernels::ArithTraits<scalar_type>::one();
  const scalar_type two   = one + one;
  const scalar_type three = two + one;

  MatrixType A("matrix A", n, n);
  auto A_h  = Kokkos::create_mirror_view(A);
  A_h(0, 0) = two;
  A_h(0, 1) = three;
  A_h(0, 2) = one;
  A_h(0, 3) = 5 * one;
  A_h(1, 0) = 6 * one;
  A_h(1, 1) = 13 * one;
  A_h(1, 2) = 5 * one;
  A_h(1, 3) = 19 * one;
  A_h(2, 0) = 2 * one;
  A_h(2, 1) = 19 * one;
  A_h(2, 2) = 10 * one;
  A_h(2, 3) = 23 * one;
  A_h(3, 0) = 4 * one;
  A_h(3, 1) = 10 * one;
  A_h(3, 2) = 11 * one;
  A_h(3, 3) = 31 * one;
  Kokkos::deep_copy(A, A_h);

  IpivType ipiv("LU pivots", n);
  InfoType info("LU info", 1);

  KokkosLapack::getrf(ExecutionSpace(), A, ipiv, info);
  Kokkos::fence();

  auto tol = 31 * n * n * n * KokkosKernels::ArithTraits<scalar_type>::eps();

  // Two rhs
  MatrixType b2("b3", n, 2);
  auto b2_h = Kokkos::create_mirror_view(b2);

  {
    b2_h(0, 0) = -6 * one;
    b2_h(0, 1) = zero;
    b2_h(1, 0) = -14 * one;
    b2_h(1, 1) = -10 * one;
    b2_h(2, 0) = 6 * one;
    b2_h(2, 1) = -41 * one;
    b2_h(3, 0) = -45 * one;
    b2_h(3, 1) = -14 * one;
    Kokkos::deep_copy(b2, b2_h);

    KokkosLapack::getrs(ExecutionSpace(), "N", A, ipiv, b2, info);
    Kokkos::fence();

    Kokkos::deep_copy(b2_h, b2);
    EXPECT_NEAR_KK_REL(b2_h(0, 0), -one, tol);
    EXPECT_NEAR_KK_REL(b2_h(1, 0), three, tol);
    EXPECT_NEAR_KK_REL(b2_h(2, 0), two, tol);
    EXPECT_NEAR_KK_REL(b2_h(3, 0), -three, tol);

    EXPECT_NEAR_KK_REL(b2_h(0, 1), two, tol);
    EXPECT_NEAR_KK_REL(b2_h(1, 1), -two, tol);
    EXPECT_NEAR_KK_REL(b2_h(2, 1), -three, tol);
    EXPECT_NEAR_KK_REL(b2_h(3, 1), one, tol);
  }

  // Test transpose modes
  for (auto mode : {"T", "C"}) {
    b2_h(0, 0) = 8 * one;
    b2_h(0, 1) = -10 * one;
    b2_h(1, 0) = 44 * one;
    b2_h(1, 1) = -67 * one;
    b2_h(2, 0) = 1 * one;
    b2_h(2, 1) = -27 * one;
    b2_h(3, 0) = 5 * one;
    b2_h(3, 1) = -66 * one;
    Kokkos::deep_copy(b2, b2_h);

    KokkosLapack::getrs(ExecutionSpace(), mode, A, ipiv, b2, info);
    Kokkos::fence();

    Kokkos::deep_copy(b2_h, b2);
    EXPECT_NEAR_KK_REL(b2_h(0, 0), -one, tol);
    EXPECT_NEAR_KK_REL(b2_h(1, 0), three, tol);
    EXPECT_NEAR_KK_REL(b2_h(2, 0), two, tol);
    EXPECT_NEAR_KK_REL(b2_h(3, 0), -three, tol);

    EXPECT_NEAR_KK_REL(b2_h(0, 1), two, tol);
    EXPECT_NEAR_KK_REL(b2_h(1, 1), -two, tol);
    EXPECT_NEAR_KK_REL(b2_h(2, 1), -three, tol);
    EXPECT_NEAR_KK_REL(b2_h(3, 1), one, tol);
  }
}

template <class AMatrixType>
void impl_test_getrs(const int n) {
  using Device         = typename AMatrixType::device_type;
  using IpivType       = Kokkos::View<int *, Device>;
  using InfoType       = Kokkos::View<int *, Device>;
  using scalar_type    = typename AMatrixType::non_const_value_type;
  using ExecutionSpace = typename Device::execution_space;

  ExecutionSpace space{};

  const auto tol = 100 * n * n * n * KokkosKernels::ArithTraits<scalar_type>::eps();

  AMatrixType A("matrix A", n, n);

  IpivType ipiv("LU pivots", n);
  InfoType info("LU info", 1);

  Kokkos::Random_XorShift64_Pool<ExecutionSpace> rand_pool(13718);
  Kokkos::fill_random(A, rand_pool, 100);

  for (int nrhs = 1; nrhs < 5; ++nrhs) {
    AMatrixType x("x", n, nrhs), b("b", n, nrhs);
    Kokkos::fill_random(x, rand_pool, 100);

    // Compute b = 0.0 * b + 1.0 * A * x
    KokkosBlas::gemm(space, "N", "N", 1.0, A, x, 0.0, b);

    // Decompose A into L and U factors
    KokkosLapack::getrf(space, A, ipiv, info);

    // Solve A * x = b for x but overwriting the results in b
    KokkosLapack::getrs(space, "N", A, ipiv, b, info);
    Kokkos::fence();

    auto b_h = Kokkos::create_mirror_view(b);
    Kokkos::deep_copy(b_h, b);
    auto x_h = Kokkos::create_mirror_view(x);
    Kokkos::deep_copy(x_h, x);

    for (int rowIdx = 0; rowIdx < b.extent_int(0); ++rowIdx) {
      for (int colIdx = 0; colIdx < b.extent_int(1); ++colIdx) {
        EXPECT_NEAR_KK_REL(b_h(rowIdx, colIdx), x_h(rowIdx, colIdx), tol);
      }
    }
  }
}

}  // namespace Test

template <class Scalar, class Device>
void test_getrs() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  using view_type_a = Kokkos::View<Scalar **, Kokkos::LayoutLeft, Device>;

  Test::impl_test_getrs_sym<view_type_a>();
  Test::impl_test_getrs_unsym<view_type_a>();

  Test::impl_test_getrs<view_type_a>(0);
  Test::impl_test_getrs<view_type_a>(1);
  Test::impl_test_getrs<view_type_a>(2);
  Test::impl_test_getrs<view_type_a>(4);
  Test::impl_test_getrs<view_type_a>(100);
#endif
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, getrs_float) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::getrs_float");
  test_getrs<float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, getrs_double) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::getrs_double");
  test_getrs<double, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, getrs_complex_float) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::getrs_complex_float");
  test_getrs<Kokkos::complex<float>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, getrs_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::getrs_complex_double");
  test_getrs<Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#endif  // CUDA+CUSOLVER or HIP+ROCSOLVER or LAPACK+HOST
#endif  // TEST_LAPACK_GETRS_HPP
