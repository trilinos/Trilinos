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

template <typename Vector>
double simpleNorm2(const Vector& v) {
  using Scalar = typename Vector::non_const_value_type;
  using KAT    = Kokkos::ArithTraits<Scalar>;
  auto vhost   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);
  double d     = 0;
  for (size_t i = 0; i < v.extent(0); i++) {
    double m = KAT::abs(vhost(i));
    d += m * m;
  }
  return Kokkos::Experimental::sqrt(d);
}

template <typename V1, typename V2>
typename V1::non_const_value_type simpleDot(const V1& v1, const V2& v2) {
  using Scalar = typename V1::non_const_value_type;
  using KAT    = Kokkos::ArithTraits<Scalar>;
  auto v1host  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v1);
  auto v2host  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v2);
  typename V1::non_const_value_type val = KAT::zero();
  for (size_t i = 0; i < v1.extent(0); i++) {
    val += v1host(i) * v2host(i);
  }
  return val;
}

// Check that all columns of X are unit length and pairwise orthogonal
template <typename Mat>
void verifyOrthogonal(const Mat& X) {
  using value_type = typename Mat::non_const_value_type;
  int k            = X.extent(1);
  for (int i = 0; i < k; i++) {
    auto col1  = Kokkos::subview(X, Kokkos::ALL(), i);
    double len = simpleNorm2(col1);
    Test::EXPECT_NEAR_KK(len, 1.0, Test::svdEpsilon<value_type>());
    for (int j = 0; j < i; j++) {
      auto col2 = Kokkos::subview(X, Kokkos::ALL(), j);
      double d  = Kokkos::ArithTraits<value_type>::abs(simpleDot(col1, col2));
      Test::EXPECT_NEAR_KK(d, 0.0, Test::svdEpsilon<value_type>());
    }
  }
}

template <typename AView, typename UView, typename VtView, typename SigmaView>
void verifySVD(const AView& A, const UView& U, const VtView& Vt,
               const SigmaView& sigma) {
  using value_type = typename AView::non_const_value_type;
  using KAT        = Kokkos::ArithTraits<value_type>;
  // Check that U/V columns are unit length and orthogonal, and that U *
  // diag(sigma) * V^T == A
  int m       = A.extent(0);
  int n       = A.extent(1);
  int maxrank = std::min(m, n);
  verifyOrthogonal(U);
  // NOTE: V^T being square and orthonormal implies that V is, so we don't have
  // to transpose it here.
  verifyOrthogonal(Vt);
  AView usvt("USV^T", m, n);
  for (int i = 0; i < maxrank; i++) {
    auto Ucol =
        Kokkos::subview(U, Kokkos::ALL(), Kokkos::make_pair<int>(i, i + 1));
    auto Vtrow =
        Kokkos::subview(Vt, Kokkos::make_pair<int>(i, i + 1), Kokkos::ALL());
    Test::vanillaGEMM(sigma(i), Ucol, Vtrow, 1.0, usvt);
  }
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      Test::EXPECT_NEAR_KK(usvt(i, j), A(i, j), Test::svdEpsilon<value_type>());
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

template <typename Matrix>
Matrix createRandomMatrix(int m, int n, int deficiency, double maxval = 1.0) {
  using Scalar = typename Matrix::non_const_value_type;
  Matrix mat("A", m, n);
  auto mhost = Kokkos::create_mirror_view(mat);
  // Fill mat with random values first
  if (maxval != 0.0) {
    Kokkos::Random_XorShift64_Pool<Kokkos::DefaultHostExecutionSpace> rand_pool(
        13718);
    Scalar minrand, maxrand;
    Test::getRandomBounds<Scalar>(maxval, minrand, maxrand);
    Kokkos::fill_random(mhost, rand_pool, minrand, maxrand);
  }
  // Apply the rank deficiency.
  // If m < n, make some rows a multiple of the first row.
  // Otherwise, make some columns a multiple of the first column.
  if (m < n) {
    for (int i = 0; i < deficiency; i++) {
      // make row i + 1 a multiple of row 0
      for (int j = 0; j < n; j++) {
        mhost(i + 1, j) = (double)(i + 2) * mhost(0, j);
      }
    }
  } else {
    for (int i = 0; i < deficiency; i++) {
      // make col i + 1 a multiple of col 0
      for (int j = 0; j < m; j++) {
        mhost(j, i + 1) = (double)(i + 2) * mhost(j, 0);
      }
    }
  }
  Kokkos::deep_copy(mat, mhost);
  return mat;
}

template <typename Matrix, typename Vector>
struct SerialSVDFunctor_Full {
  SerialSVDFunctor_Full(const Matrix& A_, const Matrix& U_, const Matrix& Vt_,
                        const Vector& sigma_, const Vector& work_)
      : A(A_), U(U_), Vt(Vt_), sigma(sigma_), work(work_) {}

  // NOTE: this functor is only meant to be launched with a single element range
  // policy
  KOKKOS_INLINE_FUNCTION void operator()(int) const {
    KokkosBatched::SerialSVD::invoke(KokkosBatched::SVD_USV_Tag(), A, U, sigma,
                                     Vt, work);
  }

  Matrix A;
  Matrix U;
  Matrix Vt;
  Vector sigma;
  Vector work;
};

template <typename Matrix, typename Vector>
struct SerialSVDFunctor_SingularValuesOnly {
  SerialSVDFunctor_SingularValuesOnly(const Matrix& A_, const Vector& sigma_,
                                      const Vector& work_)
      : A(A_), sigma(sigma_), work(work_) {}

  // NOTE: this functor is only meant to be launched with a single element range
  // policy
  KOKKOS_INLINE_FUNCTION void operator()(int) const {
    KokkosBatched::SerialSVD::invoke(KokkosBatched::SVD_S_Tag(), A, sigma,
                                     work);
  }

  Matrix A;
  Vector sigma;
  Vector work;
};

template <typename Scalar, typename Layout, typename Device>
void testSerialSVD(int m, int n, int deficiency, double maxval = 1.0) {
  using Matrix    = Kokkos::View<Scalar**, Layout, Device>;
  using Vector    = Kokkos::View<Scalar*, Device>;
  using ExecSpace = typename Device::execution_space;
  Matrix A        = createRandomMatrix<Matrix>(m, n, deficiency, maxval);
  // Fill U, Vt, sigma with nonzeros as well to make sure they are properly
  // overwritten
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
  Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecSpace>(0, 1),
      SerialSVDFunctor_Full<Matrix, Vector>(A, U, Vt, sigma, work));
  // Get the results back
  auto Uhost  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), U);
  auto Vthost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Vt);
  auto sigmaHost =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), sigma);
  // Verify the SVD is correct
  verifySVD(Acopy, Uhost, Vthost, sigmaHost);
}

template <typename Scalar, typename Layout, typename Device>
void testSerialSVDSingularValuesOnly(int m, int n) {
  using Matrix    = Kokkos::View<Scalar**, Layout, Device>;
  using Vector    = Kokkos::View<Scalar*, Device>;
  using ExecSpace = typename Device::execution_space;
  Matrix A        = createRandomMatrix<Matrix>(m, n, 0);
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
  Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecSpace>(0, 1),
      SerialSVDFunctor_Full<Matrix, Vector>(A, U, Vt, sigma1, work));
  Kokkos::deep_copy(A, Acopy);
  // Run the same SVD (singular values only mode)
  Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecSpace>(0, 1),
      SerialSVDFunctor_SingularValuesOnly<Matrix, Vector>(A, sigma2, work));
  auto sigma1Host =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), sigma1);
  auto sigma2Host =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), sigma2);
  // Make sure they match
  for (int i = 0; i < maxrank; i++) {
    Test::EXPECT_NEAR_KK(sigma1Host(i), sigma2Host(i),
                         Test::svdEpsilon<Scalar>());
  }
}

// Test the bidiagonal n*n SVD step where the last diagonal entry is 0
template <typename Scalar, typename Layout>
void testSerialSVDZeroLastRow(int n) {
  // Generate a bidiagonal matrix
  using Matrix = Kokkos::View<Scalar**, Layout, Kokkos::HostSpace>;
  using KAT    = Kokkos::ArithTraits<Scalar>;
  Matrix B     = createRandomMatrix<Matrix>(n, n, 0, 1.0);
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
  KokkosBatched::SerialSVDInternal::svdZeroLastColumn<Scalar>(
      B.data(), n, B.stride(0), B.stride(1), Vt.data(), Vt.stride(0),
      Vt.stride(1));
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
  Test::EXPECT_NEAR_KK(B(n - 2, n - 1), KAT::zero(),
                       Test::svdEpsilon<Scalar>());
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
  int m = n + 2;  // Make U somewhat bigger to make sure the Givens transforms
                  // are applied correctly
  Matrix B = createRandomMatrix<Matrix>(m, n, 0, 1.0);
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
  KokkosBatched::SerialSVDInternal::svdZeroRow<Scalar>(
      row, B.data(), n, B.stride(0), B.stride(1), U.data(), m, U.stride(0),
      U.stride(1));
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
  testSerialSVD<Scalar, Layout, Device>(0, 0, 0);
  testSerialSVD<Scalar, Layout, Device>(1, 0, 0);
  testSerialSVD<Scalar, Layout, Device>(0, 1, 0);
  testSerialSVD<Scalar, Layout, Device>(2, 2, 0);
  testSerialSVD<Scalar, Layout, Device>(2, 2, 1);
  testSerialSVD<Scalar, Layout, Device>(10, 8, 0);
  testSerialSVD<Scalar, Layout, Device>(8, 10, 0);
  testSerialSVD<Scalar, Layout, Device>(10, 1, 0);
  testSerialSVD<Scalar, Layout, Device>(1, 10, 0);
  testSerialSVD<Scalar, Layout, Device>(10, 8, 3);
  testSerialSVD<Scalar, Layout, Device>(8, 10, 4);
  // Test with all-zero matrix
  testSerialSVD<Scalar, Layout, Device>(8, 10, 0, 0.0);
  // Test some important internal routines which are not called often
  testSerialSVDZeroLastRow<Scalar, Layout>(10);
  testSerialSVDZeroDiagonal<Scalar, Layout>(10, 3);
  // Test the mode that just computes singular values
  testSerialSVDSingularValuesOnly<Scalar, Layout, Device>(10, 8);
}

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, batched_scalar_serial_svd_double) {
  // Test general SVD on a few different input sizes (full rank randomized)
  testSVD<double, Kokkos::LayoutLeft, TestExecSpace>();
  testSVD<double, Kokkos::LayoutRight, TestExecSpace>();
}
#endif

#if defined(KOKKOSKERNELS_INST_FLOAT)
TEST_F(TestCategory, batched_scalar_serial_svd_float) {
  // Test general SVD on a few different input sizes (full rank randomized)
  testSVD<float, Kokkos::LayoutLeft, TestExecSpace>();
  testSVD<float, Kokkos::LayoutRight, TestExecSpace>();
}
#endif
