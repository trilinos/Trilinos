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
/// \author Brian Kelley (bmkelle@sandia.gov)

#include "KokkosBatched_SVD_Decl.hpp"             //For testing overall kernel
#include "KokkosBatched_SVD_Serial_Internal.hpp"  //For unit testing individual components
#include "KokkosBatched_SetIdentity_Decl.hpp"

namespace Test {
template <typename Scalar>
Scalar svdEpsilon() {
  throw std::runtime_error("Unsupported scalar type");
}

template <>
double svdEpsilon() {
  return 1e-13;
}

template <>
float svdEpsilon() {
  return 2e-6f;
}
}  // namespace Test

// NOTE: simpleDot and simpleNorm2 currently support only real scalars (OK since
// SVD does as well)
template <typename V1, typename V2>
typename V1::non_const_value_type simpleDot(const V1& v1, const V2& v2) {
  using Scalar = typename V1::non_const_value_type;
  Scalar d;
  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<typename V1::execution_space>(0, v1.extent(0)),
      KOKKOS_LAMBDA(int i, Scalar& ld) { ld += v1(i) * v2(i); }, d);
  return d;
}
template <typename V>
typename V::non_const_value_type simpleNorm2(const V& v) {
  return Kokkos::sqrt(simpleDot(v, v));
}

// Check that all columns of X are unit length and pairwise orthogonal
template <typename Mat>
void verifyOrthogonal(const Mat& X) {
  using Scalar = typename Mat::non_const_value_type;
  int k        = X.extent(1);
  for (int i = 0; i < k; i++) {
    auto col1  = Kokkos::subview(X, Kokkos::ALL(), i);
    double len = simpleNorm2(col1);
    Test::EXPECT_NEAR_KK(len, 1.0, Test::svdEpsilon<Scalar>());
    for (int j = 0; j < i; j++) {
      auto col2 = Kokkos::subview(X, Kokkos::ALL(), j);
      double d  = Kokkos::ArithTraits<Scalar>::abs(simpleDot(col1, col2));
      Test::EXPECT_NEAR_KK(d, 0.0, Test::svdEpsilon<Scalar>());
    }
  }
}

template <typename AView, typename UView, typename VtView, typename SigmaView>
void verifySVD(const AView& A, const UView& U, const VtView& Vt, const SigmaView& sigma) {
  using Scalar = typename AView::non_const_value_type;
  using KAT    = Kokkos::ArithTraits<Scalar>;
  // Check that U/V columns are unit length and orthogonal, and that U *
  // diag(sigma) * V^T == A
  int m       = A.extent(0);
  int n       = A.extent(1);
  int maxrank = std::min(m, n);
  verifyOrthogonal(U);
  // NOTE: V^T being square and orthonormal implies that V is, so we don't have
  // to transpose it here.
  verifyOrthogonal(Vt);
  Kokkos::View<Scalar**, typename AView::device_type> usvt("USV^T", m, n);
  for (int i = 0; i < maxrank; i++) {
    auto Ucol  = Kokkos::subview(U, Kokkos::ALL(), Kokkos::make_pair<int>(i, i + 1));
    auto Vtrow = Kokkos::subview(Vt, Kokkos::make_pair<int>(i, i + 1), Kokkos::ALL());
    Test::vanillaGEMM(sigma(i), Ucol, Vtrow, 1.0, usvt);
  }
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      Test::EXPECT_NEAR_KK(usvt(i, j), A(i, j), Test::svdEpsilon<Scalar>());
    }
  }
  // Make sure all singular values are positive
  for (int i = 0; i < maxrank; i++) {
    EXPECT_GE(sigma(i), KAT::zero());
  }
  // Make sure singular values are in descending order
  for (int i = 0; i < maxrank - 1; i++) {
    EXPECT_GE(sigma(i), sigma(i + 1));
  }
}

template <typename Matrix, typename Vector>
struct SerialSVDFunctor_Full {
  SerialSVDFunctor_Full(const Matrix& A_, const Matrix& U_, const Matrix& Vt_, const Vector& sigma_,
                        const Vector& work_)
      : A(A_), U(U_), Vt(Vt_), sigma(sigma_), work(work_) {}

  // NOTE: this functor is only meant to be launched with a single element range
  // policy
  KOKKOS_INLINE_FUNCTION void operator()(int) const {
    KokkosBatched::SerialSVD::invoke(KokkosBatched::SVD_USV_Tag(), A, U, sigma, Vt, work);
  }

  Matrix A;
  Matrix U;
  Matrix Vt;
  Vector sigma;
  Vector work;
};

template <typename Matrix, typename Vector>
struct SerialSVDFunctor_SingularValuesOnly {
  SerialSVDFunctor_SingularValuesOnly(const Matrix& A_, const Vector& sigma_, const Vector& work_)
      : A(A_), sigma(sigma_), work(work_) {}

  // NOTE: this functor is only meant to be launched with a single element range
  // policy
  KOKKOS_INLINE_FUNCTION void operator()(int) const {
    KokkosBatched::SerialSVD::invoke(KokkosBatched::SVD_S_Tag(), A, sigma, work);
  }

  Matrix A;
  Vector sigma;
  Vector work;
};

template <typename Matrix>
Matrix randomMatrixWithRank(int m, int n, int rank) {
  Matrix A("A", m, n);
  if (rank == Kokkos::min(m, n)) {
    // A is full-rank so as a shortcut, fill it with random values directly.
    Kokkos::Random_XorShift64_Pool<typename Matrix::device_type> rand_pool(13318);
    Kokkos::fill_random(A, rand_pool, -1.0, 1.0);
  } else {
    // A is rank-deficient, so compute it as a product of two random matrices
    using MatrixHost = typename Matrix::HostMirror;
    auto Ahost       = Kokkos::create_mirror_view(A);
    Kokkos::Random_XorShift64_Pool<Kokkos::DefaultHostExecutionSpace> rand_pool(13318);
    MatrixHost U("U", m, rank);
    MatrixHost Vt("Vt", rank, n);
    Kokkos::fill_random(U, rand_pool, -1.0, 1.0);
    Kokkos::fill_random(Vt, rand_pool, -1.0, 1.0);
    Test::vanillaGEMM(1.0, U, Vt, 0.0, Ahost);
    Kokkos::deep_copy(A, Ahost);
  }
  return A;
}

template <typename Matrix>
Matrix randomMatrixWithRank(int m, int n) {
  return randomMatrixWithRank<Matrix>(m, n, Kokkos::min(m, n));
}

template <typename Scalar, typename Layout, typename Device>
void testSerialSVD(int m, int n, int rank) {
  using Matrix    = Kokkos::View<Scalar**, Layout, Device>;
  using Vector    = Kokkos::View<Scalar*, Device>;
  using ExecSpace = typename Device::execution_space;
  Matrix A        = randomMatrixWithRank<Matrix>(m, n, rank);
  // Fill U, Vt, sigma with nonzeros as well to make sure they are properly
  // overwritten
  Matrix U("U", m, m);
  Matrix Vt("Vt", n, n);
  int maxrank = std::min(m, n);
  Vector sigma("sigma", maxrank);
  Vector work("work", std::max(m, n));
  // Fill these views with an arbitrary value, to make sure SVD
  // doesn't rely on them being zero initialized.
  Kokkos::deep_copy(U, -5.0);
  Kokkos::deep_copy(Vt, -5.0);
  Kokkos::deep_copy(sigma, -5.0);
  Kokkos::deep_copy(work, -5.0);
  // Make a copy of A (before SVD) for verification, since the original will be
  // overwritten
  typename Matrix::HostMirror Acopy("Acopy", m, n);
  Kokkos::deep_copy(Acopy, A);
  // Run the SVD
  Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0, 1),
                       SerialSVDFunctor_Full<Matrix, Vector>(A, U, Vt, sigma, work));
  // Get the results back
  auto Uhost     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), U);
  auto Vthost    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Vt);
  auto sigmaHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), sigma);
  // Verify the SVD is correct
  verifySVD(Acopy, Uhost, Vthost, sigmaHost);
}

template <typename Scalar, typename Layout, typename Device>
void testSerialSVD(int m, int n) {
  testSerialSVD<Scalar, Layout, Device>(m, n, Kokkos::min(m, n));
}

template <typename Scalar, typename Layout, typename Device>
void testSerialSVDSingularValuesOnly(int m, int n) {
  using Matrix    = Kokkos::View<Scalar**, Layout, Device>;
  using Vector    = Kokkos::View<Scalar*, Device>;
  using ExecSpace = typename Device::execution_space;
  Matrix A        = randomMatrixWithRank<Matrix>(m, n);
  // Fill U, Vt, sigma with nonzeros as well to make sure they are properly
  // overwritten
  Matrix U("U", m, m);
  Matrix Vt("Vt", n, n);
  int maxrank = std::min(m, n);
  Vector sigma1("sigma", maxrank);
  Vector sigma2("sigma", maxrank);
  Vector work("work", std::max(m, n));
  Kokkos::deep_copy(U, -5.0);
  Kokkos::deep_copy(Vt, -5.0);
  Kokkos::deep_copy(sigma1, -5.0);
  Kokkos::deep_copy(sigma2, -7.0);
  Kokkos::deep_copy(work, -5.0);
  // Make a copy of A (before SVD) for verification, since the original will be
  // overwritten
  typename Matrix::HostMirror Acopy("Acopy", m, n);
  Kokkos::deep_copy(Acopy, A);
  // Run the SVD (full mode)
  Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0, 1),
                       SerialSVDFunctor_Full<Matrix, Vector>(A, U, Vt, sigma1, work));
  Kokkos::deep_copy(A, Acopy);
  // Run the same SVD (singular values only mode)
  Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0, 1),
                       SerialSVDFunctor_SingularValuesOnly<Matrix, Vector>(A, sigma2, work));
  auto sigma1Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), sigma1);
  auto sigma2Host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), sigma2);
  // Make sure they match
  for (int i = 0; i < maxrank; i++) {
    Test::EXPECT_NEAR_KK(sigma1Host(i), sigma2Host(i), Test::svdEpsilon<Scalar>());
  }
}

// Test the bidiagonal n*n SVD step where the last diagonal entry is 0
template <typename Scalar, typename Layout>
void testSerialSVDZeroLastRow(int n) {
  // Generate a bidiagonal matrix
  using Matrix = Kokkos::View<Scalar**, Layout, Kokkos::HostSpace>;
  using KAT    = Kokkos::ArithTraits<Scalar>;
  Matrix B     = randomMatrixWithRank<Matrix>(n, n);
  // Zero out entries to make B bidiagonal
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i != j && i + 1 != j) {
        B(i, j) = KAT::zero();
      }
    }
  }
  // Also zero out the final diagonal to test this routine
  B(n - 1, n - 1) = KAT::zero();
  Matrix Vt("Vt", n, n);
  KokkosBatched::SerialSetIdentity::invoke<Matrix>(Vt);
  // Compute the initial product to make sure it's maintained by the routine
  Matrix BVt("UBVt", n, n);
  Test::vanillaGEMM(1.0, B, Vt, 0.0, BVt);
  // Run the routine (just on host)
  KokkosBatched::SerialSVDInternal::svdZeroLastColumn<Scalar>(B.data(), n, B.stride(0), B.stride(1), n, Vt.data(),
                                                              Vt.stride(0), Vt.stride(1));
  // Check that B is still bidiagonal (to a tight tolerance, but not exactly
  // zero)
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i != j && i + 1 != j) {
        Test::EXPECT_NEAR_KK(B(i, j), KAT::zero(), Test::svdEpsilon<Scalar>());
      }
    }
  }
  // Check that the last superdiagonal is now zero
  Test::EXPECT_NEAR_KK(B(n - 2, n - 1), KAT::zero(), Test::svdEpsilon<Scalar>());
  // Check that the product is still maintained
  Matrix BVt2("UBVt", n, n);
  Test::vanillaGEMM(1.0, B, Vt, 0.0, BVt2);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      Test::EXPECT_NEAR_KK(BVt(i, j), BVt2(i, j), Test::svdEpsilon<Scalar>());
    }
  }
  // Check that Vt is still orthogonal
  verifyOrthogonal(Vt);
}

// Test bidiagonal n*n SVD step where some diagonal i (not the last) is 0.
template <typename Scalar, typename Layout>
void testSerialSVDZeroDiagonal(int n, int row) {
  // Generate a bidiagonal matrix
  using Matrix = Kokkos::View<Scalar**, Layout, Kokkos::HostSpace>;
  using KAT    = Kokkos::ArithTraits<Scalar>;
  int m        = n + 2;  // Make U somewhat bigger to make sure the Givens transforms
                         // are applied correctly
  Matrix B = randomMatrixWithRank<Matrix>(m, n);
  // Zero out entries to make B bidiagonal
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      if (i != j && i + 1 != j) {
        B(i, j) = KAT::zero();
      }
    }
  }
  // Also zero out a diagonal to test this routine
  B(row, row) = KAT::zero();
  Matrix U("U", m, m);
  KokkosBatched::SerialSetIdentity::invoke<Matrix>(U);
  // Compute the initial product to make sure it's maintained by the routine
  Matrix UB("UB", m, n);
  Test::vanillaGEMM(1.0, U, B, 0.0, UB);
  // Run the routine (just on host)
  KokkosBatched::SerialSVDInternal::svdZeroRow<Scalar>(row, B.data(), n, B.stride(0), B.stride(1), U.data(), m,
                                                       U.stride(0), U.stride(1));
  // Check that B is still bidiagonal (to a tight tolerance, but not exactly
  // zero)
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      if (i != j && i + 1 != j) {
        Test::EXPECT_NEAR_KK(B(i, j), KAT::zero(), Test::svdEpsilon<Scalar>());
      }
    }
  }
  // Check that row's diagonal is now zero
  Test::EXPECT_NEAR_KK(B(row, row), KAT::zero(), Test::svdEpsilon<Scalar>());
  // Check that the product is still maintained
  Matrix UB2("UB", m, n);
  Test::vanillaGEMM(1.0, U, B, 0.0, UB2);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      Test::EXPECT_NEAR_KK(UB(i, j), UB2(i, j), Test::svdEpsilon<Scalar>());
    }
  }
  // Check that U is still orthogonal
  verifyOrthogonal(U);
}

template <typename Scalar, typename Layout, typename Device>
void testSVD() {
  testSerialSVD<Scalar, Layout, Device>(0, 0);
  testSerialSVD<Scalar, Layout, Device>(1, 0);
  testSerialSVD<Scalar, Layout, Device>(0, 1);
  testSerialSVD<Scalar, Layout, Device>(2, 2);
  testSerialSVD<Scalar, Layout, Device>(2, 2, 1);
  testSerialSVD<Scalar, Layout, Device>(10, 8);
  testSerialSVD<Scalar, Layout, Device>(8, 10);
  testSerialSVD<Scalar, Layout, Device>(10, 1);
  testSerialSVD<Scalar, Layout, Device>(1, 10, 0);
  testSerialSVD<Scalar, Layout, Device>(10, 8, 3);
  testSerialSVD<Scalar, Layout, Device>(8, 10, 4);
  testSerialSVD<Scalar, Layout, Device>(8, 10, 7);
  // Test some important internal routines which are not called often
  testSerialSVDZeroLastRow<Scalar, Layout>(10);
  testSerialSVDZeroDiagonal<Scalar, Layout>(10, 3);
  // Test the mode that just computes singular values
  testSerialSVDSingularValuesOnly<Scalar, Layout, Device>(10, 8);
}

template <typename ViewT>
KOKKOS_INLINE_FUNCTION constexpr auto Determinant(ViewT F)
    -> std::enable_if_t<Kokkos::is_view<ViewT>::value && ViewT::rank == 2, double> {
  return (F(0, 0) * F(1, 1) * F(2, 2) + F(0, 1) * F(1, 2) * F(2, 0) + F(0, 2) * F(1, 0) * F(2, 1) -
          (F(0, 2) * F(1, 1) * F(2, 0) + F(0, 1) * F(1, 0) * F(2, 2) + F(0, 0) * F(1, 2) * F(2, 1)));
}

template <typename ExeSpace, typename ViewT>
void GenerateTestData(ViewT data) {
  using memory_space = typename ExeSpace::memory_space;
  // finite difference should return dPK2dU. So, we can analyze two cases.
  Kokkos::Random_XorShift64_Pool<memory_space> random(13718);
  Kokkos::fill_random(data, random, 1.0);
  Kokkos::parallel_for(
      Kokkos::RangePolicy<ExeSpace>(0, data.extent(0)), KOKKOS_LAMBDA(int i) {
        auto data_i = Kokkos::subview(data, i, Kokkos::ALL(), Kokkos::ALL());
        while (Determinant(data_i) < 0.5) {
          data_i(0, 0) += 1.0;
          data_i(1, 1) += 1.0;
          data_i(2, 2) += 1.0;
        }
      });
}

template <typename Scalar, typename Layout, typename Device, int N = 3>
void testIssue1786() {
  using execution_space   = typename Device::execution_space;
  using memory_space      = typename Device::memory_space;
  constexpr int num_tests = 4;
  Kokkos::View<Scalar* [3][3], Layout, memory_space> matrices("data", num_tests);
  GenerateTestData<execution_space>(matrices);
  Kokkos::View<Scalar* [N][N], Layout, memory_space> Us("Us", matrices.extent(0));
  Kokkos::View<Scalar* [N], Layout, memory_space> Ss("Ss", matrices.extent(0));
  Kokkos::View<Scalar* [N][N], Layout, memory_space> Vts("Vts", matrices.extent(0));
  // Make sure the 2nd dimension of works is contiguous
  Kokkos::View<Scalar* [N], Kokkos::LayoutRight, memory_space> works("works", matrices.extent(0));
  Kokkos::View<Scalar* [N][N], Layout, memory_space> matrices_copy("matrices_copy", matrices.extent(0));
  // make a copy of the input data to avoid overwriting it
  Kokkos::deep_copy(matrices_copy, matrices);
  auto policy = Kokkos::RangePolicy<execution_space>(0, matrices.extent(0));
  Kokkos::parallel_for(
      "polar decomposition", policy, KOKKOS_LAMBDA(int i) {
        auto matrix_copy = Kokkos::subview(matrices_copy, i, Kokkos::ALL(), Kokkos::ALL());
        auto U           = Kokkos::subview(Us, i, Kokkos::ALL(), Kokkos::ALL());
        auto S           = Kokkos::subview(Ss, i, Kokkos::ALL());
        auto Vt          = Kokkos::subview(Vts, i, Kokkos::ALL(), Kokkos::ALL());
        auto work        = Kokkos::subview(works, i, Kokkos::ALL());
        KokkosBatched::SerialSVD::invoke(KokkosBatched::SVD_USV_Tag{}, matrix_copy, U, S, Vt, work);
      });

  auto Us_h       = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, Us);
  auto Ss_h       = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, Ss);
  auto Vts_h      = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, Vts);
  auto matrices_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, matrices);
  for (int i = 0; i < num_tests; i++) {
    auto A  = Kokkos::subview(matrices_h, i, Kokkos::ALL(), Kokkos::ALL());
    auto U  = Kokkos::subview(Us_h, i, Kokkos::ALL(), Kokkos::ALL());
    auto S  = Kokkos::subview(Ss_h, i, Kokkos::ALL());
    auto Vt = Kokkos::subview(Vts_h, i, Kokkos::ALL(), Kokkos::ALL());
    verifySVD(A, U, Vt, S);
  }
}

// Generate specific test cases
template <typename Scalar, typename Layout, typename Device>
Kokkos::View<Scalar**, Layout, Device> getTestCase(int testCase) {
  using MatrixHost = Kokkos::View<Scalar**, Layout, Kokkos::HostSpace>;
  MatrixHost Ahost;
  int m, n;
  switch (testCase) {
    case 0:
      // Issue #2344 case 1
      m           = 3;
      n           = 3;
      Ahost       = MatrixHost("A0", m, n);
      Ahost(1, 0) = 3.58442287931538747e-02;
      Ahost(1, 1) = 3.81743062695684907e-02;
      Ahost(2, 2) = -5.55555555555555733e-02;
      break;
    case 1:
      // Test a matrix that is strictly lower triangular (so the diagonal
      // is zero)
      m     = 8;
      n     = 8;
      Ahost = MatrixHost("A1", m, n);
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < i; j++) {
          Ahost(i, j) = 1;
        }
      }
      break;
    case 2:
      // Test a matrix that's already diagonal, except for one superdiagonal in the middle
      m     = 10;
      n     = 5;
      Ahost = MatrixHost("A2", m, n);
      for (int i = 0; i < n; i++) Ahost(i, i) = 1.0;
      Ahost(2, 3) = 2.2;
      break;
    case 3:
      // Test a matrix that is already bidiagonal, and has a zero diagonal in the middle
      m     = 10;
      n     = 7;
      Ahost = MatrixHost("A3", m, n);
      for (int i = 0; i < n; i++) Ahost(i, i) = 1.0;
      for (int i = 0; i < n - 1; i++) Ahost(i, i + 1) = 0.7;
      Ahost(4, 4) = 0;
      break;
    case 4: {
      // Issue #2344 case 2
      m           = 3;
      n           = 4;
      Ahost       = MatrixHost("A4", m, n);
      Ahost(0, 0) = -2.0305040121856084e-02;
      Ahost(1, 0) = 0.0000000000000000e+00;
      Ahost(2, 0) = 0.0000000000000000e+00;
      Ahost(0, 1) = -0.0000000000000000e+00;
      Ahost(1, 1) = -0.0000000000000000e+00;
      Ahost(2, 1) = 1.9506119814028472e-02;
      Ahost(0, 2) = -2.0305040121856091e-02;
      Ahost(1, 2) = 0.0000000000000000e+00;
      Ahost(2, 2) = 0.0000000000000000e+00;
      Ahost(0, 3) = -0.0000000000000000e+00;
      Ahost(1, 3) = -0.0000000000000000e+00;
      Ahost(2, 3) = 1.9506119814028472e-02;
      break;
    }
    case 5: {
      // Test with all-zero matrix
      m     = 17;
      n     = 19;
      Ahost = MatrixHost("A5", m, n);
      break;
    }
    default: throw std::runtime_error("Test case out of bounds.");
  }
  Kokkos::View<Scalar**, Layout, Device> A(Ahost.label(), m, n);
  Kokkos::deep_copy(A, Ahost);
  return A;
}

template <typename Scalar, typename Layout, typename Device>
void testSpecialCases() {
  using Matrix    = Kokkos::View<Scalar**, Layout, Device>;
  using Vector    = Kokkos::View<Scalar*, Device>;
  using ExecSpace = typename Device::execution_space;
  for (int i = 0; i < 6; i++) {
    Matrix A = getTestCase<Scalar, Layout, Device>(i);
    int m    = A.extent(0);
    int n    = A.extent(1);
    Matrix U("U", m, m);
    Matrix Vt("Vt", n, n);
    int maxrank = std::min(m, n);
    Vector sigma("sigma", maxrank);
    Vector work("work", std::max(m, n));
    Kokkos::deep_copy(U, -5.0);
    Kokkos::deep_copy(Vt, -5.0);
    Kokkos::deep_copy(sigma, -5.0);
    Kokkos::deep_copy(work, -5.0);
    // Make a copy of A (before SVD) for verification, since the original will be
    // overwritten
    typename Matrix::HostMirror Acopy("Acopy", m, n);
    Kokkos::deep_copy(Acopy, A);
    // Run the SVD
    Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0, 1),
                         SerialSVDFunctor_Full<Matrix, Vector>(A, U, Vt, sigma, work));
    // Get the results back
    auto Uhost     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), U);
    auto Vthost    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Vt);
    auto sigmaHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), sigma);

    // Verify the SVD is correct
    verifySVD(Acopy, Uhost, Vthost, sigmaHost);
  }
}

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, batched_scalar_serial_svd_double) {
  // Test general SVD on a few different input sizes (full rank randomized)
  testSVD<double, Kokkos::LayoutLeft, TestDevice>();
  testSVD<double, Kokkos::LayoutRight, TestDevice>();
  testIssue1786<double, Kokkos::LayoutLeft, TestDevice>();
  testIssue1786<double, Kokkos::LayoutRight, TestDevice>();
  testSpecialCases<double, Kokkos::LayoutLeft, TestDevice>();
  testSpecialCases<double, Kokkos::LayoutRight, TestDevice>();
}
#endif

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, batched_scalar_serial_svd_float) {
  // Test general SVD on a few different input sizes (full rank randomized)
  testSVD<float, Kokkos::LayoutLeft, TestDevice>();
  testSVD<float, Kokkos::LayoutRight, TestDevice>();
  testIssue1786<float, Kokkos::LayoutLeft, TestDevice>();
  testIssue1786<float, Kokkos::LayoutRight, TestDevice>();
  testSpecialCases<float, Kokkos::LayoutLeft, TestDevice>();
  testSpecialCases<float, Kokkos::LayoutRight, TestDevice>();
}
#endif
