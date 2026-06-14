// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

// Only enable this test where KokkosLapack supports geqrf:
// CUDA+CUSOLVER, HIP+ROCSOLVER and HOST+LAPACK
#if (defined(TEST_CUDA_LAPACK_CPP) && defined(KOKKOSKERNELS_ENABLE_TPL_CUSOLVER)) ||               \
    (defined(TEST_HIP_LAPACK_CPP) && defined(KOKKOSKERNELS_ENABLE_TPL_ROCSOLVER)) ||               \
    ((defined(KOKKOSKERNELS_ENABLE_TPL_LAPACK) || defined(KOKKOSKERNELS_ENABLE_TPL_ACCELERATE)) && \
     (defined(TEST_OPENMP_LAPACK_CPP) || defined(TEST_SERIAL_LAPACK_CPP) || defined(TEST_THREADS_LAPACK_CPP)))

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

#include <KokkosBlas2_ger.hpp>
#include <KokkosBlas3_gemm.hpp>
#include <KokkosLapack_geqrf.hpp>
#include <KokkosLapack_gemqr.hpp>
#include <KokkosLapack_gegqr.hpp>
#include <KokkosKernels_TestUtils.hpp>

namespace Test {

template <class ViewTypeA, class ViewTypeTau, class Device>
void impl_test_gegqr(int m, int n) {
  using ALayout_t       = typename ViewTypeA::array_layout;
  using ViewTypeInfo    = Kokkos::View<int*, ALayout_t, Device>;
  using execution_space = typename Device::execution_space;
  using ScalarA         = typename ViewTypeA::non_const_value_type;
  using ats             = KokkosKernels::ArithTraits<ScalarA>;

  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);

  const int minMN = std::min(m, n);

  // ********************************************************************
  // Create device views
  // ********************************************************************
  ViewTypeA A("A", m, n);
  ViewTypeA Aorig("Aorig", m, n);
  ViewTypeTau Tau("Tau", minMN);
  ViewTypeInfo Info("Info", 1);

  // ********************************************************************
  // Create host mirrors of device views
  // ********************************************************************
  typename ViewTypeA::host_mirror_type h_A       = Kokkos::create_mirror_view(A);
  typename ViewTypeA::host_mirror_type h_Aorig   = Kokkos::create_mirror_view(Aorig);
  typename ViewTypeTau::host_mirror_type h_tau   = Kokkos::create_mirror_view(Tau);
  typename ViewTypeInfo::host_mirror_type h_info = Kokkos::create_mirror_view(Info);

  // ********************************************************************
  // Initialize data
  // ********************************************************************
  Kokkos::fill_random(A, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, ScalarA>::max());
  Kokkos::deep_copy(h_A, A);

  Kokkos::deep_copy(Aorig, A);
  Kokkos::deep_copy(h_Aorig, h_A);
  Kokkos::fence();

  // ********************************************************************
  // Perform the QR factorization
  // ********************************************************************
  execution_space space{};
  try {
    KokkosLapack::geqrf(space, A, Tau, Info);
  } catch (const std::runtime_error& e) {
    std::cout << "KokkosLapack::geqrf(): caught exception '" << e.what() << "'" << std::endl;
    FAIL();
    return;
  }
  Kokkos::fence();

  Kokkos::deep_copy(h_info, Info);
  EXPECT_EQ(h_info[0], 0) << "Failed geqrf() test: Info[0] = " << h_info[0];

  // ********************************************************************
  // Get the results
  // ********************************************************************
  Kokkos::deep_copy(h_A, A);
  Kokkos::deep_copy(h_tau, Tau);

  typename ats::mag_type absTol;
  if constexpr (std::is_same_v<typename ats::mag_type, float>) {
    absTol = 5.e-5;
  } else {
    absTol = 1.e-8;
  }

  // ********************************************************************
  // Compute Q, R, and QR
  // ********************************************************************
  ViewTypeA Qref("Qref", m, m);
  ViewTypeA QR("QR", m, n);

  typename ViewTypeA::host_mirror_type h_Qref = Kokkos::create_mirror_view(Qref);
  typename ViewTypeA::host_mirror_type h_QR   = Kokkos::create_mirror_view(QR);

  // Load identity matrix in Qref
  for (int idx = 0; idx < m; ++idx) {
    h_Qref(idx, idx) = 1.0;
  }
  Kokkos::deep_copy(Qref, h_Qref);

  // Apply Q stored in A to our Qref that is currently set as the identity
  KokkosLapack::gemqr(space, "L", "N", A, Tau, Qref, Info);
  Kokkos::fence();

  Kokkos::deep_copy(h_Qref, Qref);

  // Finally, compute Q using gegqr
  KokkosLapack::gegqr(space, minMN, A, Tau, Info);
  Kokkos::deep_copy(h_A, A);

  // ********************************************************************
  // Check that A = QR
  // ********************************************************************
  {
    bool test_flag_QR = true;
    for (int i = 0; (i < m) && test_flag_QR; ++i) {
      for (int j = 0; (j < n) && test_flag_QR; ++j) {
        if ((ats::abs(h_Qref(i, j) - h_A(i, j)) > absTol)) {
          std::cout << "Q checking"
                    << ", m = " << m << ", n = " << n << ", i = " << i << ", j = " << j
                    << ", h_A(i,j) = " << std::setprecision(16) << h_A(i, j)
                    << ", h_Qref(i,j) = " << std::setprecision(16) << h_Qref(i, j)
                    << ", |diff| = " << std::setprecision(16) << ats::abs(h_Qref(i, j) - h_A(i, j))
                    << ", absTol = " << std::setprecision(16) << absTol << std::endl;
          test_flag_QR = false;
        }
      }
    }
    ASSERT_EQ(test_flag_QR, true);
  }
}

template <class ViewTypeA, class ViewTypeTau, class Device>
void computeQ_analytic() {
  using ALayout_t       = typename ViewTypeA::array_layout;
  using ViewTypeInfo    = Kokkos::View<int*, ALayout_t, Device>;
  using execution_space = typename Device::execution_space;
  using Scalar          = typename ViewTypeA::value_type;

  ViewTypeA A("A", 3, 3);
  ViewTypeTau Tau("tau", 3);
  ViewTypeInfo Info("Info", 1);
  ViewTypeA Qref("Q ref", 3, 3);

  typename ViewTypeA::host_mirror_type h_A    = Kokkos::create_mirror_view(A);
  typename ViewTypeA::host_mirror_type h_Qref = Kokkos::create_mirror_view(Qref);

  h_A(0, 0) = 12;
  h_A(0, 1) = -51;
  h_A(0, 2) = 4;
  h_A(1, 0) = 6;
  h_A(1, 1) = 167;
  h_A(1, 2) = -68;
  h_A(2, 0) = -4;
  h_A(2, 1) = 24;
  h_A(2, 2) = -41;
  Kokkos::deep_copy(A, h_A);

  h_Qref(0, 0) = -6. / 7.;
  h_Qref(0, 1) = 69. / 175.;
  h_Qref(0, 2) = 58. / 175.;
  h_Qref(1, 0) = -3. / 7.;
  h_Qref(1, 1) = -158. / 175.;
  h_Qref(1, 2) = -6. / 175.;
  h_Qref(2, 0) = 2. / 7.;
  h_Qref(2, 1) = -6. / 35.;
  h_Qref(2, 2) = 33. / 35.;
  Kokkos::deep_copy(Qref, h_Qref);

  try {
    execution_space space{};
    KokkosLapack::geqrf(space, A, Tau, Info);
    KokkosLapack::gegqr(space, 3, A, Tau, Info);
  } catch (const std::runtime_error& e) {
    std::cout << "KokkosLapack::gegqr(): caught exception '" << e.what() << "'" << std::endl;
    FAIL();
    return;
  }
  Kokkos::fence();

  Kokkos::deep_copy(h_A, A);
  for (int rowIdx = 0; rowIdx < 3; ++rowIdx) {
    for (int colIdx = 0; colIdx < 3; ++colIdx) {
      Test::EXPECT_NEAR_KK_REL(h_Qref(rowIdx, colIdx), h_A(rowIdx, colIdx),
                               10 * KokkosKernels::ArithTraits<Scalar>::eps());
    }
  }
}

}  // namespace Test

template <class Scalar, class Device>
void test_gegqr() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  using view_type_a_ll   = Kokkos::View<Scalar**, Kokkos::LayoutLeft, Device>;
  using view_type_tau_ll = Kokkos::View<Scalar*, Kokkos::LayoutLeft, Device>;

  Test::computeQ_analytic<view_type_a_ll, view_type_tau_ll, Device>();

  Test::impl_test_gegqr<view_type_a_ll, view_type_tau_ll, Device>(1, 1);
  Test::impl_test_gegqr<view_type_a_ll, view_type_tau_ll, Device>(2, 1);
  Test::impl_test_gegqr<view_type_a_ll, view_type_tau_ll, Device>(2, 2);
  Test::impl_test_gegqr<view_type_a_ll, view_type_tau_ll, Device>(3, 1);
  Test::impl_test_gegqr<view_type_a_ll, view_type_tau_ll, Device>(3, 2);
  Test::impl_test_gegqr<view_type_a_ll, view_type_tau_ll, Device>(3, 3);

  // Note that {un,or}gqr require M >= N >= 0
  // see https://www.netlib.org/lapack/double/dorgqr.f
  // so not testing 70x100 with this rountine!
  Test::impl_test_gegqr<view_type_a_ll, view_type_tau_ll, Device>(100, 70);
  Test::impl_test_gegqr<view_type_a_ll, view_type_tau_ll, Device>(100, 100);
#endif
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, gegqr_float) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::gegqr_float");
  test_gegqr<float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, gegqr_double) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::gegqr_double");
  test_gegqr<double, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, gegqr_complex_float) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::gegqr_complex_float");
  test_gegqr<Kokkos::complex<float>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, gegqr_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::gegqr_complex_double");
  test_gegqr<Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#endif  // CUDA+CUSOLVER or HIP+ROCSOLVER or LAPACK+HOST
