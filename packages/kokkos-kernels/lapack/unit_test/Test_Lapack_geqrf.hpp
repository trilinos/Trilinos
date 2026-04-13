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
#include <KokkosKernels_TestUtils.hpp>

namespace Test {

template <class ViewTypeA, class ViewTypeTau>
void getQR(int const m, int const n, typename ViewTypeA::host_mirror_type const& h_A,
           typename ViewTypeTau::host_mirror_type const& h_tau, typename ViewTypeA::host_mirror_type& h_Q,
           typename ViewTypeA::host_mirror_type& h_R, typename ViewTypeA::host_mirror_type& h_QR) {
  using ScalarA = typename ViewTypeA::value_type;

  // ********************************************************************
  // Populate h_R
  // ********************************************************************
  for (int i = 0; i < m; ++i) {
    for (int j(0); j < n; ++j) {
      if (i <= j) {
        h_R(i, j) = h_A(i, j);
      } else {
        h_R(i, j) = KokkosKernels::ArithTraits<ScalarA>::zero();
      }
    }
  }

  // ********************************************************************
  // Instantiate the m x m identity matrix h_I
  // ********************************************************************
  ViewTypeA I("I", m, m);
  typename ViewTypeA::host_mirror_type h_I = Kokkos::create_mirror_view(I);
  for (int i = 0; i < m; ++i) {
    h_I(i, i) = ScalarA(1.);
  }

  // ********************************************************************
  // Compute h_Q
  // ********************************************************************
  int minMN(std::min(m, n));
  ViewTypeTau v("v", m);
  typename ViewTypeTau::host_mirror_type h_v = Kokkos::create_mirror_view(v);

  ViewTypeA Qk("Qk", m, m);
  typename ViewTypeA::host_mirror_type h_Qk = Kokkos::create_mirror_view(Qk);

  ViewTypeA auxM("auxM", m, m);
  typename ViewTypeA::host_mirror_type h_auxM = Kokkos::create_mirror_view(auxM);

  // Q = H(0) H(1) . . . H(min(M,N)-1), where for k=0,1,...,min(m,n)-1:
  //   H(k) = I - Tau(k) * v * v**H, and
  //   v is a vector of size m with:
  //     v(0:k-1) = 0,
  //     v(k)     = 1,
  //     v(k+1:m-1) = A(k+1:m-1,k).
  for (int k = 0; k < minMN; ++k) {
    Kokkos::deep_copy(h_v, KokkosKernels::ArithTraits<ScalarA>::zero());
    h_v[k] = 1.;
    for (int index(k + 1); index < m; ++index) {
      h_v[index] = h_A(index, k);
    }

    // Rank-1 update of a general matrix: A = A + alpha * x * y^{T,H}.
    // void ger( const char                                   trans[]
    //         , const typename AViewType::const_value_type & alpha
    //         , const XViewType                            & x
    //         , const YViewType                            & y
    //         , const AViewType                            & A
    //         );
    Kokkos::deep_copy(h_Qk, h_I);
    KokkosBlas::ger("H", -h_tau[k], h_v, h_v, h_Qk);

    // Dense matrix-matrix multiply: C = beta*C + alpha*op(A)*op(B).
    // void gemm( const char                             transA[]
    //          , const char                             transB[]
    //          , typename AViewType::const_value_type & alpha
    //          , const AViewType                      & A
    //          , const BViewType                      & B
    //          , typename CViewType::const_value_type & beta
    //          , const CViewType                      & C
    //          );
    if (k == 0) {
      Kokkos::deep_copy(h_Q, h_Qk);
    } else {
      Kokkos::deep_copy(h_auxM, KokkosKernels::ArithTraits<ScalarA>::zero());
      KokkosBlas::gemm("N", "N", 1., h_Q, h_Qk, 0., h_auxM);
      Kokkos::deep_copy(h_Q, h_auxM);
    }
  }  // for k

  // ********************************************************************
  // Check that Q^H Q = I
  // ********************************************************************
  {
    Kokkos::deep_copy(h_auxM, KokkosKernels::ArithTraits<ScalarA>::zero());
    KokkosBlas::gemm("C", "N", 1., h_Q, h_Q, 0., h_auxM);

    typename KokkosKernels::ArithTraits<typename ViewTypeA::non_const_value_type>::mag_type absTol(1.e-8);
    if constexpr (std::is_same_v<
                      typename KokkosKernels::ArithTraits<typename ViewTypeA::non_const_value_type>::mag_type, float>) {
      absTol = 5.e-5;
    }

    using ats          = KokkosKernels::ArithTraits<ScalarA>;
    bool test_flag_QHQ = true;
    for (int i = 0; (i < m) && test_flag_QHQ; ++i) {
      for (int j = 0; (j < m) && test_flag_QHQ; ++j) {
        if (ats::abs(h_auxM(i, j) - h_I(i, j)) > absTol) {
          std::cout << "QHQ checking"
                    << ", m = " << m << ", n = " << n << ", i = " << i << ", j = " << j
                    << ", h_auxM(i,j) = " << std::setprecision(16) << h_auxM(i, j)
                    << ", h_I(i,j) = " << std::setprecision(16) << h_I(i, j) << ", |diff| = " << std::setprecision(16)
                    << ats::abs(h_auxM(i, j) - h_I(i, j)) << ", absTol = " << std::setprecision(16) << absTol
                    << std::endl;
          test_flag_QHQ = false;
        }
      }
    }
    ASSERT_EQ(test_flag_QHQ, true);
  }

  // ********************************************************************
  // Compute h_QR
  // ********************************************************************
  Kokkos::deep_copy(h_QR, KokkosKernels::ArithTraits<ScalarA>::zero());
  KokkosBlas::gemm("N", "N", 1., h_Q, h_R, 0., h_QR);
}

template <class ViewTypeA, class ViewTypeTau, class Device>
void impl_test_geqrf(int m, int n) {
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
  try {
    execution_space space{};
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

  getQR<ViewTypeA, ViewTypeTau>(m, n, h_A, h_tau, h_Q, h_R, h_QR);

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
        for (int j = 0; (j < n) && test_flag_Q; ++j) {
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
        if (ats::abs(h_QR(i, j) - h_Aorig(i, j)) > absTol) {
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

}  // namespace Test

template <class Scalar, class Device>
void test_geqrf() {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  using view_type_a_ll   = Kokkos::View<Scalar**, Kokkos::LayoutLeft, Device>;
  using view_type_tau_ll = Kokkos::View<Scalar*, Kokkos::LayoutLeft, Device>;

  Test::impl_test_geqrf<view_type_a_ll, view_type_tau_ll, Device>(1, 1);
  Test::impl_test_geqrf<view_type_a_ll, view_type_tau_ll, Device>(2, 1);
  Test::impl_test_geqrf<view_type_a_ll, view_type_tau_ll, Device>(2, 2);
  Test::impl_test_geqrf<view_type_a_ll, view_type_tau_ll, Device>(3, 1);
  Test::impl_test_geqrf<view_type_a_ll, view_type_tau_ll, Device>(3, 2);
  Test::impl_test_geqrf<view_type_a_ll, view_type_tau_ll, Device>(3, 3);

  Test::impl_test_geqrf<view_type_a_ll, view_type_tau_ll, Device>(100, 70);
  Test::impl_test_geqrf<view_type_a_ll, view_type_tau_ll, Device>(70, 100);
  Test::impl_test_geqrf<view_type_a_ll, view_type_tau_ll, Device>(100, 100);
#endif
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, geqrf_float) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::geqrf_float");
  test_geqrf<float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, geqrf_double) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::geqrf_double");
  test_geqrf<double, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, geqrf_complex_float) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::geqrf_complex_float");
  test_geqrf<Kokkos::complex<float>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, geqrf_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::geqrf_complex_double");
  test_geqrf<Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#endif  // CUDA+CUSOLVER or HIP+ROCSOLVER or LAPACK+HOST
