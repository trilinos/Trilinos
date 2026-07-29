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
#include <KokkosKernels_TestUtils.hpp>

namespace Test {

template <class ViewTypeA, class ViewTypeTau, class Device>
void impl_test_gemqr(int m, int n) {
  using ALayout_t       = typename ViewTypeA::array_layout;
  using ViewTypeInfo    = Kokkos::View<int*, ALayout_t, Device>;
  using execution_space = typename Device::execution_space;
  using ScalarA         = typename ViewTypeA::value_type;
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
  if ((m == 3) && (n == 3)) {
    h_A(0, 0) = ScalarA(12.);
    h_A(0, 1) = ScalarA(-51.);
    h_A(0, 2) = ScalarA(4.);

    h_A(1, 0) = ScalarA(6.);
    h_A(1, 1) = ScalarA(167.);
    h_A(1, 2) = ScalarA(-68.);

    h_A(2, 0) = ScalarA(-4.);
    h_A(2, 1) = ScalarA(24.);
    h_A(2, 2) = ScalarA(-41.);

    Kokkos::deep_copy(A, h_A);
  } else {
    Kokkos::fill_random(A, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, ScalarA>::max());
    Kokkos::deep_copy(h_A, A);
  }

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

  typename KokkosKernels::ArithTraits<typename ViewTypeA::non_const_value_type>::mag_type absTol(1.e-8);
  if constexpr (std::is_same_v<typename KokkosKernels::ArithTraits<typename ViewTypeA::non_const_value_type>::mag_type,
                               float>) {
    absTol = 5.e-5;
  }

  // ********************************************************************
  // Check outputs h_A and h_tau
  // ********************************************************************
  if ((m == 3) && (n == 3)) {
    Kokkos::View<ScalarA**, Kokkos::HostSpace> refMatrix("ref matrix", m, n);
    Kokkos::View<ScalarA*, Kokkos::HostSpace> refTau("ref tau", m);

    refMatrix(0, 0) = ScalarA(-14.);
    refMatrix(0, 1) = ScalarA(-21.);
    refMatrix(0, 2) = ScalarA(14.);

    refMatrix(1, 0) = ScalarA(0.2307692307692308);
    refMatrix(1, 1) = ScalarA(-175.);
    refMatrix(1, 2) = ScalarA(70.);

    refMatrix(2, 0) = ScalarA(-0.1538461538461539);
    refMatrix(2, 1) = ScalarA(1. / 18.);
    refMatrix(2, 2) = ScalarA(-35.);

    refTau(0) = ScalarA(1.857142857142857);
    refTau(1) = ScalarA(1.993846153846154);
    refTau(2) = ScalarA(0.);

    {
      bool test_flag_A = true;
      for (int i = 0; (i < m) && test_flag_A; ++i) {
        for (int j = 0; (j < n) && test_flag_A; ++j) {
          if (ats::abs(h_A(i, j) - refMatrix(i, j)) > absTol) {
            std::cout << "h_Aoutput checking"
                      << ", m = " << m << ", n = " << n << ", i = " << i << ", j = " << j
                      << ", h_Aoutput(i,j) = " << std::setprecision(16) << h_A(i, j)
                      << ", refMatrix(i,j) = " << std::setprecision(16) << refMatrix(i, j)
                      << ", |diff| = " << std::setprecision(16) << ats::abs(h_A(i, j) - refMatrix(i, j))
                      << ", absTol = " << std::setprecision(16) << absTol << std::endl;
            test_flag_A = false;
          }
        }
      }
      ASSERT_EQ(test_flag_A, true);
    }

    {
      bool test_flag_tau = true;
      for (int i = 0; (i < m) && test_flag_tau; ++i) {
        if (ats::abs(h_tau(i) - refTau(i)) > absTol) {
          std::cout << "tau checking"
                    << ", m = " << m << ", n = " << n << ", i = " << i << ", h_tau(i,j) = " << std::setprecision(16)
                    << h_tau(i) << ", refTau(i,j) = " << std::setprecision(16) << refTau(i)
                    << ", |diff| = " << std::setprecision(16) << ats::abs(h_tau(i) - refTau(i))
                    << ", absTol = " << std::setprecision(16) << absTol << std::endl;
          test_flag_tau = false;
        }
      }
      ASSERT_EQ(test_flag_tau, true);
    }
  }

  // ********************************************************************
  // Compute Q, R, and QR
  // ********************************************************************
  ViewTypeA Q("Q", m, m);
  ViewTypeA R("R", m, n);
  ViewTypeA QR("QR", m, n);

  typename ViewTypeA::host_mirror_type h_Q  = Kokkos::create_mirror_view(Q);
  typename ViewTypeA::host_mirror_type h_R  = Kokkos::create_mirror_view(R);
  typename ViewTypeA::host_mirror_type h_QR = Kokkos::create_mirror_view(QR);

  // Load identity matrix in Q
  for (int idx = 0; idx < m; ++idx) {
    h_Q(idx, idx) = 1.0;
  }
  Kokkos::deep_copy(Q, h_Q);

  // Load R from A
  for (int rowIdx = 0; rowIdx < minMN; ++rowIdx) {
    for (int colIdx = 0; colIdx < n; ++colIdx) {
      if (rowIdx <= colIdx) {
        h_R(rowIdx, colIdx) = h_A(rowIdx, colIdx);
      }
    }
  }
  Kokkos::deep_copy(R, h_R);

  // Apply Q stored in A to our Q that is currently set as the identity
  KokkosLapack::gemqr(space, "L", "N", A, Tau, Q, Info);
  Kokkos::deep_copy(h_Q, Q);

  // Recompute A from Q and R factors
  KokkosBlas::gemm("N", "N", 1., Q, R, 0., QR);
  Kokkos::deep_copy(h_QR, QR);

  // ********************************************************************
  // Check Q, R, and QR
  // ********************************************************************
  if ((m == 3) && (n == 3)) {
    Kokkos::View<ScalarA**, typename Kokkos::DefaultHostExecutionSpace::memory_space> refQ("ref Q", m, n);
    Kokkos::View<ScalarA**, typename Kokkos::DefaultHostExecutionSpace::memory_space> refR("ref R", m, n);

    refQ(0, 0) = ScalarA(-6. / 7.);
    refQ(0, 1) = ScalarA(69. / 175.);
    refQ(0, 2) = ScalarA(58. / 175.);

    refQ(1, 0) = ScalarA(-3. / 7.);
    refQ(1, 1) = ScalarA(-158. / 175.);
    refQ(1, 2) = ScalarA(-6. / 175.);

    refQ(2, 0) = ScalarA(2. / 7.);
    refQ(2, 1) = ScalarA(-6. / 35.);
    refQ(2, 2) = ScalarA(33. / 35.);

    refR(0, 0) = ScalarA(-14.);
    refR(0, 1) = ScalarA(-21.);
    refR(0, 2) = ScalarA(14.);

    refR(1, 1) = ScalarA(-175.);
    refR(1, 2) = ScalarA(70.);

    refR(2, 2) = ScalarA(-35.);

    {
      bool test_flag_Q = true;
      for (int i = 0; (i < m) && test_flag_Q; ++i) {
        for (int j = 0; (j < m) && test_flag_Q; ++j) {
          if (ats::abs(h_Q(i, j) - refQ(i, j)) > absTol) {
            std::cout << "Q checking"
                      << ", m = " << m << ", n = " << n << ", i = " << i << ", j = " << j
                      << ", h_Q(i,j) = " << std::setprecision(16) << h_Q(i, j)
                      << ", refQ(i,j) = " << std::setprecision(16) << refQ(i, j)
                      << ", |diff| = " << std::setprecision(16) << ats::abs(h_Q(i, j) - refQ(i, j))
                      << ", absTol = " << std::setprecision(16) << absTol << std::endl;
            test_flag_Q = false;
          }
        }
      }
      ASSERT_EQ(test_flag_Q, true);
    }

    {
      bool test_flag_R = true;
      for (int i = 0; (i < m) && test_flag_R; ++i) {
        for (int j = 0; (j < n) && test_flag_R; ++j) {
          if (ats::abs(h_R(i, j) - refR(i, j)) > absTol) {
            std::cout << "R checking"
                      << ", m = " << m << ", n = " << n << ", i = " << i << ", j = " << j
                      << ", h_R(i,j) = " << std::setprecision(16) << h_R(i, j)
                      << ", refR(i,j) = " << std::setprecision(16) << refR(i, j)
                      << ", |diff| = " << std::setprecision(16) << ats::abs(h_R(i, j) - refR(i, j))
                      << ", absTol = " << std::setprecision(16) << absTol << std::endl;
            test_flag_R = false;
          }
        }
      }
      ASSERT_EQ(test_flag_R, true);
    }
  }

  // ********************************************************************
  // Check that A = QR
  // ********************************************************************
  {
    bool test_flag_QR = true;
    for (int i = 0; (i < m) && test_flag_QR; ++i) {
      for (int j = 0; (j < n) && test_flag_QR; ++j) {
        if ((ats::abs(h_QR(i, j) - h_Aorig(i, j)) > absTol)) {
          std::cout << "QR checking"
                    << ", m = " << m << ", n = " << n << ", i = " << i << ", j = " << j
                    << ", h_Aorig(i,j) = " << std::setprecision(16) << h_Aorig(i, j)
                    << ", h_QR(i,j) = " << std::setprecision(16) << h_QR(i, j) << ", |diff| = " << std::setprecision(16)
                    << ats::abs(h_QR(i, j) - h_Aorig(i, j)) << ", absTol = " << std::setprecision(16) << absTol
                    << std::endl;
          test_flag_QR = false;
        }
      }
    }
    ASSERT_EQ(test_flag_QR, true);
  }
}

template <class ViewTypeA, class ViewTypeTau, class Device>
void applyQ_analytic() {
  using ALayout_t       = typename ViewTypeA::array_layout;
  using ViewTypeInfo    = Kokkos::View<int*, ALayout_t, Device>;
  using execution_space = typename Device::execution_space;
  using Scalar          = typename ViewTypeA::value_type;

  ViewTypeA A("A", 3, 3);
  ViewTypeTau Tau("tau", 3);
  ViewTypeInfo Info("Info", 1);
  ViewTypeA Q("Q", 3, 3);
  ViewTypeA Qref("Q ref", 3, 3);

  typename ViewTypeA::host_mirror_type h_A    = Kokkos::create_mirror_view(A);
  typename ViewTypeA::host_mirror_type h_Q    = Kokkos::create_mirror_view(Q);
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

  // Store the identity so once Q is applied to
  // this matrix we will recover the entries of Q
  h_Q(0, 0) = 1;
  h_Q(0, 1) = 0;
  h_Q(0, 2) = 0;
  h_Q(1, 0) = 0;
  h_Q(1, 1) = 1;
  h_Q(1, 2) = 0;
  h_Q(2, 0) = 0;
  h_Q(2, 1) = 0;
  h_Q(2, 2) = 1;
  Kokkos::deep_copy(Q, h_Q);

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
    KokkosLapack::gemqr(space, "L", "N", A, Tau, Q, Info);
  } catch (const std::runtime_error& e) {
    std::cout << "KokkosLapack::gemqr(): caught exception '" << e.what() << "'" << std::endl;
    FAIL();
    return;
  }
  Kokkos::fence();

  Kokkos::deep_copy(h_Q, Q);
  for (int rowIdx = 0; rowIdx < 3; ++rowIdx) {
    for (int colIdx = 0; colIdx < 3; ++colIdx) {
      Test::EXPECT_NEAR_KK_REL(h_Qref(rowIdx, colIdx), h_Q(rowIdx, colIdx),
                               10 * KokkosKernels::ArithTraits<Scalar>::eps());
    }
  }
}

}  // namespace Test

template <class Scalar, class Device>
void test_gemqr() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  using view_type_a_ll   = Kokkos::View<Scalar**, Kokkos::LayoutLeft, Device>;
  using view_type_tau_ll = Kokkos::View<Scalar*, Kokkos::LayoutLeft, Device>;

  Test::applyQ_analytic<view_type_a_ll, view_type_tau_ll, Device>();

  Test::impl_test_gemqr<view_type_a_ll, view_type_tau_ll, Device>(1, 1);
  Test::impl_test_gemqr<view_type_a_ll, view_type_tau_ll, Device>(2, 1);
  Test::impl_test_gemqr<view_type_a_ll, view_type_tau_ll, Device>(2, 2);
  Test::impl_test_gemqr<view_type_a_ll, view_type_tau_ll, Device>(3, 1);
  Test::impl_test_gemqr<view_type_a_ll, view_type_tau_ll, Device>(3, 2);
  Test::impl_test_gemqr<view_type_a_ll, view_type_tau_ll, Device>(3, 3);

  Test::impl_test_gemqr<view_type_a_ll, view_type_tau_ll, Device>(100, 70);
  Test::impl_test_gemqr<view_type_a_ll, view_type_tau_ll, Device>(70, 100);
  Test::impl_test_gemqr<view_type_a_ll, view_type_tau_ll, Device>(100, 100);
#endif
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, gemqr_float) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::gemqr_float");
  test_gemqr<float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, gemqr_double) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::gemqr_double");
  test_gemqr<double, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, gemqr_complex_float) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::gemqr_complex_float");
  test_gemqr<Kokkos::complex<float>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, gemqr_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::gemqr_complex_double");
  test_gemqr<Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#endif  // CUDA+CUSOLVER or HIP+ROCSOLVER or LAPACK+HOST
