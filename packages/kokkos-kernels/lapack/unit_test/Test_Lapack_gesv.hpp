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

// only enable this test where KokkosLapack supports gesv:
// CUDA+(MAGMA or CUSOLVER), HIP+ROCSOLVER and HOST+LAPACK
#if (defined(TEST_CUDA_LAPACK_CPP) &&                                                            \
     (defined(KOKKOSKERNELS_ENABLE_TPL_MAGMA) || defined(KOKKOSKERNELS_ENABLE_TPL_CUSOLVER))) || \
    (defined(TEST_HIP_LAPACK_CPP) && defined(KOKKOSKERNELS_ENABLE_TPL_ROCSOLVER)) ||             \
    (defined(KOKKOSKERNELS_ENABLE_TPL_LAPACK) &&                                                 \
     (defined(TEST_OPENMP_LAPACK_CPP) || defined(TEST_SERIAL_LAPACK_CPP) || defined(TEST_THREADS_LAPACK_CPP)))

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

#include <KokkosLapack_gesv.hpp>
#include <KokkosBlas2_gemv.hpp>
#include <KokkosBlas3_gemm.hpp>
#include <KokkosKernels_TestUtils.hpp>

namespace Test {

template <class ViewTypeA, class ViewTypeB, class Device, bool MAGMA>
void impl_test_gesv(const char* mode, const char* padding, int N) {
  using execution_space = typename Device::execution_space;
  using ScalarA         = typename ViewTypeA::value_type;
  using ats             = Kokkos::ArithTraits<ScalarA>;

  execution_space space{};

  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);

  int ldda, lddb;

  if (padding[0] == 'Y') {  // rounded up to multiple of 32
    ldda = ((N + 32 - 1) / 32) * 32;
    lddb = ldda;
  } else {
    ldda = N;
    lddb = N;
  }

  // Create device views
  ViewTypeA A("A", ldda, N);
  ViewTypeB X0("X0", N);
  ViewTypeB B("B", lddb);

  // Create host mirrors of device views.
  typename ViewTypeB::HostMirror h_X0 = Kokkos::create_mirror_view(X0);
  typename ViewTypeB::HostMirror h_B  = Kokkos::create_mirror(B);

  // Initialize data.
  Kokkos::fill_random(A, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, ScalarA>::max());
  Kokkos::fill_random(X0, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, ScalarA>::max());

  // Generate RHS B = A*X0.
  ScalarA alpha = 1.0;
  ScalarA beta  = 0.0;

  KokkosBlas::gemv("N", alpha, A, X0, beta, B);
  Kokkos::fence();

  // Deep copy device view to host view.
  Kokkos::deep_copy(h_X0, X0);

  // Allocate IPIV view on host
  using ViewTypeP = typename std::conditional<MAGMA, Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::HostSpace>,
                                              Kokkos::View<int*, Kokkos::LayoutLeft, execution_space>>::type;
  ViewTypeP ipiv;
  int Nt = 0;
  if (mode[0] == 'Y') {
    Nt   = N;
    ipiv = ViewTypeP("IPIV", Nt);
  }

  // Solve.
  try {
    KokkosLapack::gesv(space, A, B, ipiv);
  } catch (const std::runtime_error& error) {
    // Check for expected runtime errors due to:
    // no-pivoting case (note: only MAGMA supports no-pivoting interface)
    // and no-tpl case
    bool nopivot_runtime_err = false;
    bool notpl_runtime_err   = false;
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA   // have MAGMA TPL
#ifdef KOKKOSKERNELS_ENABLE_TPL_LAPACK  // and have LAPACK TPL
    nopivot_runtime_err = (!std::is_same<typename Device::memory_space, Kokkos::CudaSpace>::value) &&
                          (ipiv.extent(0) == 0) && (ipiv.data() == nullptr);
    notpl_runtime_err = false;
#else
    notpl_runtime_err = true;
#endif
#else                                   // not have MAGMA TPL
#ifdef KOKKOSKERNELS_ENABLE_TPL_LAPACK  // but have LAPACK TPL
    nopivot_runtime_err = (ipiv.extent(0) == 0) && (ipiv.data() == nullptr);
    notpl_runtime_err   = false;
#else
    notpl_runtime_err = true;
#endif
#endif
    if (!nopivot_runtime_err && !notpl_runtime_err) FAIL();
    return;
  }
  Kokkos::fence();

  // Get the solution vector.
  Kokkos::deep_copy(h_B, B);

  // Checking vs ref on CPU, this eps is about 10^-9
  typedef typename ats::mag_type mag_type;
  const mag_type eps = 3.0e7 * ats::epsilon();
  bool test_flag     = true;
  for (int i = 0; i < N; i++) {
    if (ats::abs(h_B(i) - h_X0(i)) > eps) {
      test_flag = false;
      printf(
          "    Error %d, pivot %c, padding %c: result( %.15lf ) !="
          "solution( %.15lf ) at (%d), error=%.15e, eps=%.15e\n",
          N, mode[0], padding[0], ats::abs(h_B(i)), ats::abs(h_X0(i)), int(i), ats::abs(h_B(i) - h_X0(i)), eps);
      break;
    }
  }
  ASSERT_EQ(test_flag, true);
}

template <class ViewTypeA, class ViewTypeB, class Device, bool MAGMA>
void impl_test_gesv_mrhs(const char* mode, const char* padding, int N, int nrhs) {
  using execution_space = typename Device::execution_space;
  using ScalarA         = typename ViewTypeA::value_type;
  using ats             = Kokkos::ArithTraits<ScalarA>;

  execution_space space{};

  Kokkos::Random_XorShift64_Pool<execution_space> rand_pool(13718);

  int ldda, lddb;

  if (padding[0] == 'Y') {  // rounded up to multiple of 32
    ldda = ((N + 32 - 1) / 32) * 32;
    lddb = ldda;
  } else {
    ldda = N;
    lddb = N;
  }

  // Create device views
  ViewTypeA A("A", ldda, N);
  ViewTypeB X0("X0", N, nrhs);
  ViewTypeB B("B", lddb, nrhs);

  // Create host mirrors of device views.
  typename ViewTypeB::HostMirror h_X0 = Kokkos::create_mirror_view(X0);
  typename ViewTypeB::HostMirror h_B  = Kokkos::create_mirror(B);

  // Initialize data.
  Kokkos::fill_random(A, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, ScalarA>::max());
  Kokkos::fill_random(X0, rand_pool, Kokkos::rand<Kokkos::Random_XorShift64<execution_space>, ScalarA>::max());

  // Generate RHS B = A*X0.
  ScalarA alpha = 1.0;
  ScalarA beta  = 0.0;

  KokkosBlas::gemm("N", "N", alpha, A, X0, beta, B);
  Kokkos::fence();

  // Deep copy device view to host view.
  Kokkos::deep_copy(h_X0, X0);

  // Allocate IPIV view on host
  using ViewTypeP = typename std::conditional<MAGMA, Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::HostSpace>,
                                              Kokkos::View<int*, Kokkos::LayoutLeft, execution_space>>::type;
  ViewTypeP ipiv;
  int Nt = 0;
  if (mode[0] == 'Y') {
    Nt   = N;
    ipiv = ViewTypeP("IPIV", Nt);
  }

  // Solve.
  try {
    KokkosLapack::gesv(space, A, B, ipiv);
  } catch (const std::runtime_error& error) {
    // Check for expected runtime errors due to:
    // no-pivoting case (note: only MAGMA supports no-pivoting interface)
    // and no-tpl case
    bool nopivot_runtime_err = false;
    bool notpl_runtime_err   = false;
#ifdef KOKKOSKERNELS_ENABLE_TPL_MAGMA   // have MAGMA TPL
#ifdef KOKKOSKERNELS_ENABLE_TPL_LAPACK  // and have LAPACK TPL
    nopivot_runtime_err = (!std::is_same<typename Device::memory_space, Kokkos::CudaSpace>::value) &&
                          (ipiv.extent(0) == 0) && (ipiv.data() == nullptr);
    notpl_runtime_err = false;
#else
    notpl_runtime_err = true;
#endif
#else                                   // not have MAGMA TPL
#ifdef KOKKOSKERNELS_ENABLE_TPL_LAPACK  // but have LAPACK TPL
    nopivot_runtime_err = (ipiv.extent(0) == 0) && (ipiv.data() == nullptr);
    notpl_runtime_err   = false;
#else
    notpl_runtime_err = true;
#endif
#endif
    if (!nopivot_runtime_err && !notpl_runtime_err) FAIL();
    return;
  }
  Kokkos::fence();

  // Get the solution vector.
  Kokkos::deep_copy(h_B, B);

  // Checking vs ref on CPU, this eps is about 10^-9
  typedef typename ats::mag_type mag_type;
  const mag_type eps = 1.0e7 * ats::epsilon();
  bool test_flag     = true;
  for (int j = 0; j < nrhs; j++) {
    for (int i = 0; i < N; i++) {
      if (ats::abs(h_B(i, j) - h_X0(i, j)) > eps) {
        test_flag = false;
        // printf( "    Error %d, pivot %c, padding %c: result( %.15lf ) !=
        // solution( %.15lf ) at (%ld) at rhs %d\n", N, mode[0], padding[0],
        // ats::abs(h_B(i,j)), ats::abs(h_X0(i,j)), i, j );
        break;
      }
    }
    if (test_flag == false) break;
  }
  ASSERT_EQ(test_flag, true);
}

}  // namespace Test

template <class Scalar, class Device>
int test_gesv(const char* mode) {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  using view_type_a_ll = Kokkos::View<Scalar**, Kokkos::LayoutLeft, Device>;
  using view_type_b_ll = Kokkos::View<Scalar*, Kokkos::LayoutLeft, Device>;

#if (defined(TEST_CUDA_LAPACK_CPP) && defined(KOKKOSKERNELS_ENABLE_TPL_CUSOLVER)) || \
    (defined(TEST_HIP_LAPACK_CPP) && defined(KOKKOSKERNELS_ENABLE_TPL_ROCSOLVER)) || \
    (defined(KOKKOSKERNELS_ENABLE_TPL_LAPACK) &&                                     \
     (defined(TEST_OPENMP_LAPACK_CPP) || defined(TEST_SERIAL_LAPACK_CPP) || defined(TEST_THREADS_LAPACK_CPP)))
  Test::impl_test_gesv<view_type_a_ll, view_type_b_ll, Device, false>(&mode[0], "N", 2);     // no padding
  Test::impl_test_gesv<view_type_a_ll, view_type_b_ll, Device, false>(&mode[0], "N", 13);    // no padding
  Test::impl_test_gesv<view_type_a_ll, view_type_b_ll, Device, false>(&mode[0], "N", 179);   // no padding
  Test::impl_test_gesv<view_type_a_ll, view_type_b_ll, Device, false>(&mode[0], "N", 64);    // no padding
  Test::impl_test_gesv<view_type_a_ll, view_type_b_ll, Device, false>(&mode[0], "N", 1024);  // no padding

#elif defined(KOKKOSKERNELS_ENABLE_TPL_MAGMA) && defined(KOKKOS_ENABLE_CUDA)
  if constexpr (std::is_same_v<Kokkos::Cuda, typename Device::execution_space>) {
    Test::impl_test_gesv<view_type_a_ll, view_type_b_ll, Device, true>(&mode[0], "N", 2);     // no padding
    Test::impl_test_gesv<view_type_a_ll, view_type_b_ll, Device, true>(&mode[0], "N", 13);    // no padding
    Test::impl_test_gesv<view_type_a_ll, view_type_b_ll, Device, true>(&mode[0], "N", 179);   // no padding
    Test::impl_test_gesv<view_type_a_ll, view_type_b_ll, Device, true>(&mode[0], "N", 64);    // no padding
    Test::impl_test_gesv<view_type_a_ll, view_type_b_ll, Device, true>(&mode[0], "N", 1024);  // no padding

    Test::impl_test_gesv<view_type_a_ll, view_type_b_ll, Device, true>(&mode[0], "Y",
                                                                       13);  // padding
    Test::impl_test_gesv<view_type_a_ll, view_type_b_ll, Device, true>(&mode[0], "Y",
                                                                       179);  // padding
  }
#endif
#endif

  // Supress unused parameters on CUDA10
  (void)mode;
  return 1;
}

template <class Scalar, class Device>
int test_gesv_mrhs(const char* mode) {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  using view_type_a_ll = Kokkos::View<Scalar**, Kokkos::LayoutLeft, Device>;
  using view_type_b_ll = Kokkos::View<Scalar**, Kokkos::LayoutLeft, Device>;

#if (defined(TEST_CUDA_LAPACK_CPP) && defined(KOKKOSKERNELS_ENABLE_TPL_CUSOLVER)) || \
    (defined(TEST_HIP_LAPACK_CPP) && defined(KOKKOSKERNELS_ENABLE_TPL_ROCSOLVER)) || \
    (defined(KOKKOSKERNELS_ENABLE_TPL_LAPACK) &&                                     \
     (defined(TEST_OPENMP_LAPACK_CPP) || defined(TEST_SERIAL_LAPACK_CPP) || defined(TEST_THREADS_LAPACK_CPP)))
  Test::impl_test_gesv_mrhs<view_type_a_ll, view_type_b_ll, Device, false>(&mode[0], "N", 2, 5);     // no padding
  Test::impl_test_gesv_mrhs<view_type_a_ll, view_type_b_ll, Device, false>(&mode[0], "N", 13, 5);    // no padding
  Test::impl_test_gesv_mrhs<view_type_a_ll, view_type_b_ll, Device, false>(&mode[0], "N", 179, 5);   // no padding
  Test::impl_test_gesv_mrhs<view_type_a_ll, view_type_b_ll, Device, false>(&mode[0], "N", 64, 5);    // no padding
  Test::impl_test_gesv_mrhs<view_type_a_ll, view_type_b_ll, Device, false>(&mode[0], "N", 1024, 5);  // no padding

// When appropriate run MAGMA specific tests
#elif defined(KOKKOSKERNELS_ENABLE_TPL_MAGMA) && defined(KOKKOS_ENABLE_CUDA)
  if constexpr (std::is_same_v<Kokkos::Cuda, typename Device::execution_space>) {
    Test::impl_test_gesv_mrhs<view_type_a_ll, view_type_b_ll, Device, true>(&mode[0], "N", 2, 5);     // no padding
    Test::impl_test_gesv_mrhs<view_type_a_ll, view_type_b_ll, Device, true>(&mode[0], "N", 13, 5);    // no padding
    Test::impl_test_gesv_mrhs<view_type_a_ll, view_type_b_ll, Device, true>(&mode[0], "N", 179, 5);   // no padding
    Test::impl_test_gesv_mrhs<view_type_a_ll, view_type_b_ll, Device, true>(&mode[0], "N", 64, 5);    // no padding
    Test::impl_test_gesv_mrhs<view_type_a_ll, view_type_b_ll, Device, true>(&mode[0], "N", 1024, 5);  // no padding

    Test::impl_test_gesv_mrhs<view_type_a_ll, view_type_b_ll, Device, true>(&mode[0], "Y", 13, 5);   // padding
    Test::impl_test_gesv_mrhs<view_type_a_ll, view_type_b_ll, Device, true>(&mode[0], "Y", 179, 5);  // padding
  }
#endif
#endif

  // Supress unused parameters on CUDA10
  (void)mode;
  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, gesv_float) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::gesv_float");
  test_gesv<float, TestDevice>("N");  // No pivoting
  test_gesv<float, TestDevice>("Y");  // Partial pivoting
  Kokkos::Profiling::popRegion();
}

TEST_F(TestCategory, gesv_mrhs_float) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::gesv_mrhs_float");
  test_gesv_mrhs<float, TestDevice>("N");  // No pivoting
  test_gesv_mrhs<float, TestDevice>("Y");  // Partial pivoting
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, gesv_double) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::gesv_double");
  test_gesv<double, TestDevice>("N");  // No pivoting
  test_gesv<double, TestDevice>("Y");  // Partial pivoting
  Kokkos::Profiling::popRegion();
}

TEST_F(TestCategory, gesv_mrhs_double) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::gesv_mrhs_double");
  test_gesv_mrhs<double, TestDevice>("N");  // No pivoting
  test_gesv_mrhs<double, TestDevice>("Y");  // Partial pivoting
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, gesv_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::gesv_complex_double");
  test_gesv<Kokkos::complex<double>, TestDevice>("N");  // No pivoting
  test_gesv<Kokkos::complex<double>, TestDevice>("Y");  // Partial pivoting
  Kokkos::Profiling::popRegion();
}

TEST_F(TestCategory, gesv_mrhs_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::gesv_mrhs_complex_double");
  test_gesv_mrhs<Kokkos::complex<double>, TestDevice>("N");  // No pivoting
  test_gesv_mrhs<Kokkos::complex<double>, TestDevice>("Y");  // Partial pivoting
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, gesv_complex_float) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::gesv_complex_float");
  test_gesv<Kokkos::complex<float>, TestDevice>("N");  // No pivoting
  test_gesv<Kokkos::complex<float>, TestDevice>("Y");  // Partial pivoting
  Kokkos::Profiling::popRegion();
}

TEST_F(TestCategory, gesv_mrhs_complex_float) {
  Kokkos::Profiling::pushRegion("KokkosLapack::Test::gesv_mrhs_complex_float");
  test_gesv_mrhs<Kokkos::complex<float>, TestDevice>("N");  // No pivoting
  test_gesv_mrhs<Kokkos::complex<float>, TestDevice>("Y");  // Partial pivoting
  Kokkos::Profiling::popRegion();
}
#endif

#endif  // CUDA+(MAGMA or CUSOLVER) or HIP+ROCSOLVER or LAPACK+HOST
