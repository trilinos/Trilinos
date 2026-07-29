// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

// Only enable this test where KokkosLapack supports potrf/potrs:
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
#include <KokkosLapack_potrs.hpp>
#include <KokkosKernels_TestUtils.hpp>
#include <KokkosKernels_TestMatrixUtils.hpp>

namespace Test {

/// \brief Test potrs (Cholesky solve) for a given matrix size, uplo, and nrhs.
///
/// Strategy:
///  1. Construct a strictly SPD matrix A.
///  2. Choose an exact solution X_exact (known or random).
///  3. Compute B = A * X_exact.
///  4. Factor A with potrf(uplo, A).
///  5. Solve with potrs(uplo, A, B)  -> B becomes X_solve.
///  6. Verify X_solve ≈ X_exact.
///
/// For N==3, a fixed SPD matrix with known values is used so that rounding
/// errors can be checked against a tight absolute tolerance.
/// For N>3, a random B is used to build A = B^H*B + N*I.
template <class AViewType, class BViewType, class Device>
void impl_test_potrs(int N, int nrhs, const char uplo) {
  using ScalarA         = typename AViewType::non_const_value_type;
  using ScalarB         = typename BViewType::value_type;
  using ats             = KokkosKernels::ArithTraits<ScalarA>;
  using MagnitudeType   = typename ats::mag_type;
  using execution_space = typename Device::execution_space;
  using ALayout_t       = typename AViewType::array_layout;

  static_assert(std::is_same_v<ScalarA, ScalarB>, "A and B must have the same scalar type");

  // potrs TPL specializations require column-major storage
  if constexpr (!std::is_same_v<ALayout_t, Kokkos::LayoutLeft>) return;

  const char uplo_arr[] = {uplo, '\0'};

  MagnitudeType absTol = ats::epsilon() * MagnitudeType(1000);
  if constexpr (std::is_same_v<MagnitudeType, float>) absTol = MagnitudeType(1e-3);

  // ====================================================================
  // Small test with known values (N == 3)
  //
  //   A = | 4  2  2 |
  //       | 2  5  3 |
  //       | 2  3  6 |
  //
  //   X_exact (2 RHS):  col-0 = [1, 1, 1]^T,  col-1 = [1, 0, -1]^T
  //
  //   B = A * X_exact:
  //     col-0: [4+2+2, 2+5+3, 2+3+6] = [8, 10, 11]
  //     col-1: [4+0-2, 2+0-3, 2+0-6] = [2, -1, -4]
  // ====================================================================
  if (N == 3) {
    typename AViewType::non_const_type A("A", 3, 3);
    BViewType B("B", 3, 2);

    // clang-format off
    std::vector<std::vector<ScalarA>> A_data = {
        {4., 2., 2.},
        {2., 5., 3.},
        {2., 3., 6.}};
    // B = A * X_exact (column-major), X_exact = [[1,1],[1,0],[1,-1]]
    std::vector<std::vector<ScalarB>> B_data = {
        { 8.,  2.},
        {10., -1.},
        {11., -4.}};
    // clang-format on
    Test::fill_view_from_fixture(A, A_data);
    Test::fill_view_from_fixture(B, B_data);

    // Factor (A is non-const, suitable for potrf)
    KokkosLapack::potrf(uplo_arr, A);

    // Solve: B -> X_solve  (Aconst is a const view of the same factored data)
    AViewType Aconst(A);
    KokkosLapack::potrs(uplo_arr, Aconst, B);

    Kokkos::fence();
    typename BViewType::host_mirror_type h_B = Kokkos::create_mirror_view(B);
    Kokkos::deep_copy(h_B, B);

    // Expected: X_exact
    // clang-format off
    const std::vector<std::vector<ScalarB>> ref = {
        {1.,  1.},
        {1.,  0.},
        {1., -1.}};
    // clang-format on
    bool test_flag = true;
    for (int i = 0; (i < 3) && test_flag; ++i) {
      for (int j = 0; (j < 2) && test_flag; ++j) {
        if (ats::abs(h_B(i, j) - ref[i][j]) > absTol) {
          std::cout << "potrs('" << uplo << "') N=3 solve check FAILED"
                    << " i=" << i << " j=" << j << " got=" << h_B(i, j) << " expected=" << ref[i][j]
                    << " |diff|=" << ats::abs(h_B(i, j) - ref[i][j]) << " absTol=" << absTol << std::endl;
          test_flag = false;
        }
      }
    }
    ASSERT_EQ(test_flag, true);
    return;
  }

  // ====================================================================
  // General random test
  //
  //  1. Fill Btmp with random values.
  //  2. Compute A = Btmp^H * Btmp + N * I  (strictly SPD).
  //  3. Choose random X_exact  (N x nrhs).
  //  4. Compute RHS = A * X_exact.
  //  5. Factor A with potrf.
  //  6. Solve with potrs; result overwrites RHS.
  //  7. Verify RHS ≈ X_exact element-wise.
  // ====================================================================

  MagnitudeType relTol = ats::epsilon() * MagnitudeType(N) * MagnitudeType(100);
  if constexpr (std::is_same_v<MagnitudeType, float>) relTol = MagnitudeType(N) * MagnitudeType(1e-4);

  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);

  // Build a non-const A view for factorization
  typename AViewType::non_const_type Atmp("Atmp", N, N);
  typename AViewType::non_const_type Btmp("Btmp", N, N);
  BViewType X_exact("X_exact", N, nrhs);
  BViewType RHS("RHS", N, nrhs);

  // Host mirrors
  auto h_Btmp    = Kokkos::create_mirror_view(Btmp);
  auto h_Atmp    = Kokkos::create_mirror_view(Atmp);
  auto h_X_exact = Kokkos::create_mirror_view(X_exact);

  // Random Btmp for constructing A
  Kokkos::fill_random(Btmp, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, ScalarA>::max());
  Kokkos::deep_copy(h_Btmp, Btmp);

  // h_Atmp = Btmp^H * Btmp
  KokkosBlas::gemm("C", "N", ScalarA(1), h_Btmp, h_Btmp, ScalarA(0), h_Atmp);
  // Add N*I for strict positive definiteness
  for (int i = 0; i < N; ++i) h_Atmp(i, i) += ScalarA(N);

  Kokkos::deep_copy(Atmp, h_Atmp);

  // Random X_exact
  Kokkos::fill_random(X_exact, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, ScalarB>::max());
  Kokkos::deep_copy(h_X_exact, X_exact);

  // RHS = A * X_exact  (on host, then copy to device)
  auto h_RHS = Kokkos::create_mirror_view(RHS);
  KokkosBlas::gemm("N", "N", ScalarA(1), h_Atmp, h_X_exact, ScalarA(0), h_RHS);
  Kokkos::deep_copy(RHS, h_RHS);

  // Factor A
  KokkosLapack::potrf(uplo_arr, Atmp);

  // Solve: RHS -> X_solve
  AViewType Aconst(Atmp);
  KokkosLapack::potrs(uplo_arr, Aconst, RHS);

  Kokkos::fence();
  Kokkos::deep_copy(h_RHS, RHS);

  // Verify X_solve ≈ X_exact
  bool test_flag = true;
  for (int i = 0; (i < N) && test_flag; ++i) {
    for (int j = 0; (j < nrhs) && test_flag; ++j) {
      MagnitudeType diff  = ats::abs(h_RHS(i, j) - h_X_exact(i, j));
      MagnitudeType scale = std::max(ats::abs(h_X_exact(i, j)), MagnitudeType(1));
      if (diff > relTol * scale) {
        std::cout << "potrs('" << uplo << "') solve check FAILED"
                  << " N=" << N << " nrhs=" << nrhs << " i=" << i << " j=" << j << " X_solve=" << h_RHS(i, j)
                  << " X_exact=" << h_X_exact(i, j) << " |diff|=" << diff << " relTol*scale=" << relTol * scale
                  << std::endl;
        test_flag = false;
      }
    }
  }
  ASSERT_EQ(test_flag, true);
}

}  // namespace Test

template <class Scalar, class Device>
void test_potrs() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  using AViewType = Kokkos::View<const Scalar**, Kokkos::LayoutLeft, Device>;
  using BViewType = Kokkos::View<Scalar**, Kokkos::LayoutLeft, Device>;

  // Small known-value test
  Test::impl_test_potrs<AViewType, BViewType, Device>(3, 2, 'L');
  Test::impl_test_potrs<AViewType, BViewType, Device>(3, 2, 'U');

  // Random tests: single and multiple RHS, lower and upper
  Test::impl_test_potrs<AViewType, BViewType, Device>(10, 1, 'L');
  Test::impl_test_potrs<AViewType, BViewType, Device>(10, 3, 'L');
  Test::impl_test_potrs<AViewType, BViewType, Device>(10, 1, 'U');
  Test::impl_test_potrs<AViewType, BViewType, Device>(10, 3, 'U');
  Test::impl_test_potrs<AViewType, BViewType, Device>(100, 5, 'L');
  Test::impl_test_potrs<AViewType, BViewType, Device>(100, 5, 'U');
#endif
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, potrs_float) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::potrs_float");
  test_potrs<float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, potrs_double) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::potrs_double");
  test_potrs<double, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, potrs_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::potrs_complex_double");
  test_potrs<Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, potrs_complex_float) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::potrs_complex_float");
  test_potrs<Kokkos::complex<float>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#endif  // CUDA+CUSOLVER or HIP+ROCSOLVER or LAPACK/ACCELERATE+HOST
