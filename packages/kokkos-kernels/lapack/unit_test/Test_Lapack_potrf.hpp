// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

// Only enable this test where KokkosLapack supports potrf:
// CUDA+CUSOLVER, HIP+ROCSOLVER and HOST+LAPACK/ACCELERATE
#if (defined(TEST_CUDA_LAPACK_CPP) && defined(KOKKOSKERNELS_ENABLE_TPL_CUSOLVER)) ||               \
    (defined(TEST_HIP_LAPACK_CPP) && defined(KOKKOSKERNELS_ENABLE_TPL_ROCSOLVER)) ||               \
    ((defined(KOKKOSKERNELS_ENABLE_TPL_LAPACK) || defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)) && \
     (defined(TEST_OPENMP_LAPACK_CPP) || defined(TEST_SERIAL_LAPACK_CPP) || defined(TEST_THREADS_LAPACK_CPP)))

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

#include <KokkosBlas3_gemm.hpp>
#include <KokkosLapack_potrf.hpp>
#include <KokkosKernels_TestUtils.hpp>
#include <KokkosKernels_TestMatrixUtils.hpp>

namespace Test {

/// \brief Test Cholesky factorization (potrf) for a given matrix size and uplo.
///
/// For N==3: uses a small SPD matrix with known Cholesky factors.
/// For other N: constructs A = B^H*B + N*I (strictly SPD), runs potrf,
///   reconstructs and verifies it matches the original A.
///
/// uplo='L': A = L * L^H  (lower triangular)
/// uplo='U': A = U^H * U  (upper triangular)
template <class AViewType, class Device>
void impl_test_potrf(int N, const char uplo) {
  using ScalarType      = typename AViewType::value_type;
  using ats             = KokkosKernels::ArithTraits<ScalarType>;
  using MagnitudeType   = typename ats::mag_type;
  using execution_space = typename Device::execution_space;
  using ALayout_t       = typename AViewType::array_layout;

  // potrf TPL specializations require column-major storage (LayoutLeft)
  if constexpr (!std::is_same_v<ALayout_t, Kokkos::LayoutLeft>) return;

  const char uplo_arr[] = {uplo, '\0'};

  // ********************************************************************
  // Absolute tolerance for the N==3 known-value check
  // ********************************************************************
  MagnitudeType absTol = ats::epsilon() * MagnitudeType(1000);
  if constexpr (std::is_same_v<MagnitudeType, float>) {
    absTol = MagnitudeType(1e-4);
  }

  // ====================================================================
  // Small test with known values (N == 3)
  //
  //   A = | 4  2  2 |
  //       | 2  5  3 |
  //       | 2  3  6 |
  //
  //   L (uplo='L'):  L = | 2  0  0 |    A = L * L^H
  //                      | 1  2  0 |
  //                      | 1  1  2 |
  //
  //   U (uplo='U'):  U = | 2  1  1 |    A = U^H * U
  //                      | 0  2  1 |
  //                      | 0  0  2 |
  // ====================================================================
  if (N == 3) {
    AViewType A("A", 3, 3);
    typename AViewType::host_mirror_type h_A = Kokkos::create_mirror_view(A);

    // clang-format off
    std::vector<std::vector<ScalarType>> A_data = {
        {4., 2., 2.},
        {2., 5., 3.},
        {2., 3., 6.}};
    // clang-format on
    Test::fill_view_from_fixture(A, A_data);

    KokkosLapack::potrf(uplo_arr, A);
    Kokkos::fence();
    Kokkos::deep_copy(h_A, A);

    bool test_flag = true;
    if (uplo == 'L' || uplo == 'l') {
      // Expected lower Cholesky factor; check lower triangle only
      // clang-format off
      const std::vector<std::vector<ScalarType>> refL = {
          {2., 0., 0.},
          {1., 2., 0.},
          {1., 1., 2.}};
      // clang-format on
      for (int i = 0; (i < 3) && test_flag; ++i) {
        for (int j = 0; (j <= i) && test_flag; ++j) {
          if (ats::abs(h_A(i, j) - refL[i][j]) > absTol) {
            std::cout << "potrf('L') N=3 lower-triangle check FAILED"
                      << " i=" << i << " j=" << j << " h_A(i,j)=" << h_A(i, j) << " expected=" << refL[i][j]
                      << " |diff|=" << ats::abs(h_A(i, j) - refL[i][j]) << " absTol=" << absTol << std::endl;
            test_flag = false;
          }
        }
      }
    } else {
      // Expected upper Cholesky factor; check upper triangle only
      // clang-format off
      const std::vector<std::vector<ScalarType>> refU = {
          {2., 1., 1.},
          {0., 2., 1.},
          {0., 0., 2.}};
      // clang-format on
      for (int i = 0; (i < 3) && test_flag; ++i) {
        for (int j = i; (j < 3) && test_flag; ++j) {
          if (ats::abs(h_A(i, j) - refU[i][j]) > absTol) {
            std::cout << "potrf('U') N=3 upper-triangle check FAILED"
                      << " i=" << i << " j=" << j << " h_A(i,j)=" << h_A(i, j) << " expected=" << refU[i][j]
                      << " |diff|=" << ats::abs(h_A(i, j) - refU[i][j]) << " absTol=" << absTol << std::endl;
            test_flag = false;
          }
        }
      }
    }
    ASSERT_EQ(test_flag, true);
    return;
  }

  // ====================================================================
  // General random test
  //
  //  1. Fill B with random values.
  //  2. Compute A = B^H * B + N * I  (strictly Hermitian positive definite).
  //  3. Save Aorig = A.
  //  4. Factorize: potrf(uplo, A)  ->  triangular factor stored in A.
  //  5. Zero the unused triangle to isolate the factor.
  //  6. Reconstruct: L*L^H (lower) or U^H*U (upper) and verify vs Aorig.
  // ====================================================================

  // Relative tolerance: accumulated error scales with N
  MagnitudeType relTol = ats::epsilon() * MagnitudeType(N) * MagnitudeType(100);
  if constexpr (std::is_same_v<MagnitudeType, float>) {
    relTol = MagnitudeType(N) * MagnitudeType(1e-4);
  }

  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);

  AViewType B("B", N, N);
  AViewType A("A", N, N);
  AViewType Aorig("Aorig", N, N);

  // ----------------------------------------------------------------
  // Construct A = B^H * B + N * I on the host, then copy to device
  // ----------------------------------------------------------------
  typename AViewType::host_mirror_type h_B     = Kokkos::create_mirror_view(B);
  typename AViewType::host_mirror_type h_A     = Kokkos::create_mirror_view(A);
  typename AViewType::host_mirror_type h_Aorig = Kokkos::create_mirror_view(Aorig);

  Kokkos::fill_random(B, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, ScalarType>::max());
  Kokkos::deep_copy(h_B, B);

  // h_A = B^H * B
  KokkosBlas::gemm("C", "N", ScalarType(1), h_B, h_B, ScalarType(0), h_A);
  // h_A += N * I  (ensures strict positive definiteness)
  for (int i = 0; i < N; ++i) {
    h_A(i, i) += ScalarType(N);
  }

  Kokkos::deep_copy(A, h_A);
  Kokkos::deep_copy(h_Aorig, h_A);  // save original before factorization
  Kokkos::deep_copy(Aorig, h_Aorig);

  // ----------------------------------------------------------------
  // Factorize
  // ----------------------------------------------------------------
  KokkosLapack::potrf(uplo_arr, A);
  Kokkos::fence();
  Kokkos::deep_copy(h_A, A);

  // ----------------------------------------------------------------
  // Isolate the triangular factor and reconstruct
  // ----------------------------------------------------------------
  typename AViewType::host_mirror_type h_Areconst("Areconst", N, N);
  if (uplo == 'L' || uplo == 'l') {
    // Zero strict upper triangle to isolate L
    for (int i = 0; i < N; ++i)
      for (int j = i + 1; j < N; ++j) h_A(i, j) = ats::zero();
    // Reconstruct: L * L^H
    KokkosBlas::gemm("N", "C", ScalarType(1), h_A, h_A, ScalarType(0), h_Areconst);
  } else {
    // Zero strict lower triangle to isolate U
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < i; ++j) h_A(i, j) = ats::zero();
    // Reconstruct: U^H * U
    KokkosBlas::gemm("C", "N", ScalarType(1), h_A, h_A, ScalarType(0), h_Areconst);
  }

  // ----------------------------------------------------------------
  // Verify: |Areconst(i,j) - Aorig(i,j)| / max(|Aorig(i,j)|, 1) < relTol
  // ----------------------------------------------------------------
  bool test_flag = true;
  for (int i = 0; (i < N) && test_flag; ++i) {
    for (int j = 0; (j < N) && test_flag; ++j) {
      MagnitudeType diff  = ats::abs(h_Areconst(i, j) - h_Aorig(i, j));
      MagnitudeType scale = std::max(ats::abs(h_Aorig(i, j)), MagnitudeType(1));
      if (diff > relTol * scale) {
        std::cout << "potrf('" << uplo << "') reconstruction check FAILED"
                  << " N=" << N << " i=" << i << " j=" << j << " Areconst=" << h_Areconst(i, j)
                  << " Aorig=" << h_Aorig(i, j) << " |diff|=" << diff << " relTol*scale=" << relTol * scale
                  << std::endl;
        test_flag = false;
      }
    }
  }
  ASSERT_EQ(test_flag, true);
}

}  // namespace Test

template <class Scalar, class Device>
void test_potrf() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  using AViewType = Kokkos::View<Scalar**, Kokkos::LayoutLeft, Device>;

  Test::impl_test_potrf<AViewType, Device>(3, 'L');
  Test::impl_test_potrf<AViewType, Device>(10, 'L');
  Test::impl_test_potrf<AViewType, Device>(100, 'L');
  Test::impl_test_potrf<AViewType, Device>(3, 'U');
  Test::impl_test_potrf<AViewType, Device>(10, 'U');
  Test::impl_test_potrf<AViewType, Device>(100, 'U');
#endif
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, potrf_float) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::potrf_float");
  test_potrf<float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, potrf_double) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::potrf_double");
  test_potrf<double, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, potrf_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::potrf_complex_double");
  test_potrf<Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, potrf_complex_float) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::potrf_complex_float");
  test_potrf<Kokkos::complex<float>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#endif  // CUDA+CUSOLVER or HIP+ROCSOLVER or LAPACK/ACCELERATE+HOST
