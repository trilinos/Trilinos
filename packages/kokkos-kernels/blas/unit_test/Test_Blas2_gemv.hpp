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
#include <KokkosBlas2_gemv.hpp>
#include <KokkosBlas1_dot.hpp>
#include <KokkosKernels_TestUtils.hpp>

namespace Test {
template <class ExecutionSpace, class ViewTypeA, class ViewTypeX, class ViewTypeY, class Device>
void impl_test_gemv_streams(ExecutionSpace& space, const char* mode, int M, int N) {
  typedef typename ViewTypeA::value_type ScalarA;
  typedef typename ViewTypeX::value_type ScalarX;
  typedef typename ViewTypeY::value_type ScalarY;
  typedef Kokkos::ArithTraits<ScalarY> KAT_Y;

  const ScalarA alpha                = 3;
  ScalarY beta                       = 5;
  typename KAT_Y::mag_type const eps = KAT_Y::epsilon();

  int ldx;
  int ldy;
  if (mode[0] == 'N') {
    ldx = N;
    ldy = M;
  } else {
    ldx = M;
    ldy = N;
  }

  view_stride_adapter<ViewTypeA> A("A", M, N);
  view_stride_adapter<ViewTypeX> x("X", ldx);
  view_stride_adapter<ViewTypeY> y("Y", ldy);
  view_stride_adapter<ViewTypeY> org_y("Org_Y", ldy);

  Kokkos::Random_XorShift64_Pool<ExecutionSpace> rand_pool(13718);

  constexpr double max_valX = 1;
  constexpr double max_valY = 1;
  constexpr double max_valA = 1;
  {
    ScalarX randStart, randEnd;
    Test::getRandomBounds(max_valX, randStart, randEnd);
    Kokkos::fill_random(space, x.d_view, rand_pool, randStart, randEnd);
  }
  {
    ScalarY randStart, randEnd;
    Test::getRandomBounds(max_valY, randStart, randEnd);
    Kokkos::fill_random(space, y.d_view, rand_pool, randStart, randEnd);
  }
  {
    ScalarA randStart, randEnd;
    Test::getRandomBounds(max_valA, randStart, randEnd);
    Kokkos::fill_random(space, A.d_view, rand_pool, randStart, randEnd);
  }

  const typename KAT_Y::mag_type max_error = KAT_Y::abs(alpha * max_valA * max_valX * ldx + beta * max_valY);
  const typename KAT_Y::mag_type tol       = max_error * eps * 2;  // adding small fudge factor of 2

  Kokkos::deep_copy(org_y.h_base, y.d_base);
  Kokkos::deep_copy(x.h_base, x.d_base);
  Kokkos::deep_copy(A.h_base, A.d_base);

  Kokkos::View<ScalarY*, Kokkos::HostSpace> expected("expected aAx+by", ldy);
  Kokkos::deep_copy(expected, org_y.h_view);
  vanillaGEMV(mode[0], alpha, A.h_view, x.h_view, beta, expected);

  KokkosBlas::gemv(space, mode, alpha, A.d_view, x.d_view, beta, y.d_view);
  Kokkos::deep_copy(y.h_base, y.d_base);
  int numErrors = 0;
  for (int i = 0; i < ldy; i++) {
    if (KAT_Y::abs(expected(i) - y.h_view(i)) > tol) {
      numErrors++;
      std::cerr << __FILE__ << ":" << __LINE__ << ": expected(i)=" << expected(i) << ", h_y(i)=" << y.h_view(i)
                << std::endl;
    }
  }
  EXPECT_EQ(numErrors, 0) << "Nonconst input, " << M << 'x' << N << ", alpha = " << alpha << ", beta = " << beta
                          << ", mode " << mode << ": gemv incorrect";

  Kokkos::deep_copy(space, y.d_base, org_y.h_base);
  KokkosBlas::gemv(space, mode, alpha, A.d_view, x.d_view_const, beta, y.d_view);
  Kokkos::deep_copy(y.h_base, y.d_base);
  numErrors = 0;
  Kokkos::fence();  // Wait for vanillaGEMV
  for (int i = 0; i < ldy; i++) {
    if (KAT_Y::abs(expected(i) - y.h_view(i)) > tol) numErrors++;
  }
  EXPECT_EQ(numErrors, 0) << "Const vector input, " << M << 'x' << N << ", alpha = " << alpha << ", beta = " << beta
                          << ", mode " << mode << ": gemv incorrect";

  Kokkos::deep_copy(space, y.d_base, org_y.h_base);
  KokkosBlas::gemv(space, mode, alpha, A.d_view_const, x.d_view_const, beta, y.d_view);
  Kokkos::deep_copy(y.h_base, y.d_base);
  numErrors = 0;
  for (int i = 0; i < ldy; i++) {
    if (KAT_Y::abs(expected(i) - y.h_view(i)) > tol) numErrors++;
  }
  EXPECT_EQ(numErrors, 0) << "Const matrix/vector input, " << M << 'x' << N << ", alpha = " << alpha
                          << ", beta = " << beta << ", mode " << mode << ": gemv incorrect";
  // Test once with beta = 0, but with y initially filled with NaN.
  // This should overwrite the NaNs with the correct result.
  beta = KAT_Y::zero();
  // beta changed, so update the correct answer
  vanillaGEMV(mode[0], alpha, A.h_view, x.h_view, beta, expected);
  Kokkos::deep_copy(space, y.d_view, KAT_Y::nan());
  KokkosBlas::gemv(space, mode, alpha, A.d_view, x.d_view, beta, y.d_view);
  Kokkos::deep_copy(y.h_base, y.d_base);

  Kokkos::fence();  // Wait for vanillaGEMV
  numErrors = 0;
  for (int i = 0; i < ldy; i++) {
    if (KAT_Y::isNan(y.h_view(i)) ||
        KAT_Y::abs(expected(i) - y.h_view(i)) > KAT_Y::abs(alpha * max_valA * max_valX * ldx * eps * 2)) {
      numErrors++;
      std::cerr << __FILE__ << ":" << __LINE__ << ": expected(" << i << ")=" << expected(i) << ", h_y(" << i
                << ")=" << y.h_view(i) << ", eps=" << eps << ", 1024*2*eps=" << 1024 * 2 * KAT_Y::epsilon()
                << std::endl;
    }
  }
  EXPECT_EQ(numErrors, 0) << "beta = 0, input contains NaN, A is " << M << 'x' << N << ", mode " << mode
                          << ": gemv incorrect";
}
template <class ViewTypeA, class ViewTypeX, class ViewTypeY, class Device>
void impl_test_gemv(const char* mode, int M, int N) {
  using execution_space = typename Device::execution_space;
  execution_space space;
  impl_test_gemv_streams<execution_space, ViewTypeA, ViewTypeX, ViewTypeY, Device>(space, mode, M, N);
}
}  // namespace Test

template <class ScalarA, class ScalarX, class ScalarY, class Device>
int test_gemv(const char* mode) {
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutLeft, Device> view_type_a_ll;
  typedef Kokkos::View<ScalarX*, Kokkos::LayoutLeft, Device> view_type_b_ll;
  typedef Kokkos::View<ScalarY*, Kokkos::LayoutLeft, Device> view_type_c_ll;
#if 0
  Test::impl_test_gemv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(mode,10,10);
  Test::impl_test_gemv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(mode,100,10);
  Test::impl_test_gemv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(mode,10,150);
  Test::impl_test_gemv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(mode,150,10);
  Test::impl_test_gemv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(mode,10,200);
  Test::impl_test_gemv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(mode,200,10);
#endif
  Test::impl_test_gemv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(mode, 0, 1024);
  Test::impl_test_gemv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(mode, 1024, 0);
  Test::impl_test_gemv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(mode, 13, 13);
  Test::impl_test_gemv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(mode, 13, 1024);
  Test::impl_test_gemv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(mode, 50, 40);
  Test::impl_test_gemv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(mode, 1024, 1024);
  Test::impl_test_gemv<view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(mode, 2131, 2131);
  // Test::impl_test_gemv<view_type_a_ll, view_type_b_ll, view_type_c_ll,
  // Device>(mode,132231,1024);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutRight, Device> view_type_a_lr;
  typedef Kokkos::View<ScalarX*, Kokkos::LayoutRight, Device> view_type_b_lr;
  typedef Kokkos::View<ScalarY*, Kokkos::LayoutRight, Device> view_type_c_lr;
  Test::impl_test_gemv<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(mode, 0, 1024);
  Test::impl_test_gemv<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(mode, 1024, 0);
  Test::impl_test_gemv<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(mode, 13, 13);
  Test::impl_test_gemv<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(mode, 13, 1024);
  Test::impl_test_gemv<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(mode, 50, 40);
  Test::impl_test_gemv<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(mode, 1024, 1024);
  Test::impl_test_gemv<view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(mode, 2131, 2131);
  // Test::impl_test_gemv<view_type_a_lr, view_type_b_lr, view_type_c_lr,
  // Device>(mode,132231,1024);
#endif

#if (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
  typedef Kokkos::View<ScalarA**, Kokkos::LayoutStride, Device> view_type_a_ls;
  typedef Kokkos::View<ScalarX*, Kokkos::LayoutStride, Device> view_type_b_ls;
  typedef Kokkos::View<ScalarY*, Kokkos::LayoutStride, Device> view_type_c_ls;
  Test::impl_test_gemv<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(mode, 0, 1024);
  Test::impl_test_gemv<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(mode, 1024, 0);
  Test::impl_test_gemv<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(mode, 13, 13);
  Test::impl_test_gemv<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(mode, 13, 1024);
  Test::impl_test_gemv<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(mode, 50, 40);
  Test::impl_test_gemv<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(mode, 1024, 1024);
  Test::impl_test_gemv<view_type_a_ls, view_type_b_ls, view_type_c_ls, Device>(mode, 2131, 2131);
  // Test::impl_test_gemv<view_type_a_ls, view_type_b_ls, view_type_c_ls,
  // Device>(mode,132231,1024);
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
  Test::impl_test_gemv<view_type_a_ls, view_type_b_ll, view_type_c_lr, Device>(mode, 1024, 1024);
  Test::impl_test_gemv<view_type_a_ll, view_type_b_ls, view_type_c_lr, Device>(mode, 1024, 1024);
#endif

  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, gemv_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemv_float");
  test_gemv<float, float, float, TestDevice>("N");
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemv_tran_float");
  test_gemv<float, float, float, TestDevice>("T");
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, gemv_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemv_double");
  test_gemv<double, double, double, TestDevice>("N");
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemv_tran_double");
  test_gemv<double, double, double, TestDevice>("T");
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, gemv_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemv_complex_double");
  test_gemv<Kokkos::complex<double>, Kokkos::complex<double>, Kokkos::complex<double>, TestDevice>("N");
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemv_tran_complex_double");
  test_gemv<Kokkos::complex<double>, Kokkos::complex<double>, Kokkos::complex<double>, TestDevice>("T");
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemv_conj_complex_double");
  test_gemv<Kokkos::complex<double>, Kokkos::complex<double>, Kokkos::complex<double>, TestDevice>("C");
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_INT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, gemv_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemv_int");
  test_gemv<int, int, int, TestDevice>("N");
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemv_tran_int");
  test_gemv<int, int, int, TestDevice>("T");
  Kokkos::Profiling::popRegion();
}
#endif

#if !defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS)
TEST_F(TestCategory, gemv_double_int) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemv_double_int");
  test_gemv<double, int, float, TestDevice>("N");
  Kokkos::Profiling::popRegion();

  // Kokkos::Profiling::pushRegion("KokkosBlas::Test::gemvt_double_int");
  //  test_gemv<double,int,float,TestDevice> ("T");
  // Kokkos::Profiling::popRegion();
}
#endif

template <class Scalar, class Ordinal, class Offset, class Device>
int test_gemv_streams(const char* mode) {
  using execution_space = typename Device::execution_space;
  execution_space space;
#if defined(KOKKOSKERNELS_INST_LAYOUTLEFT)
  using view_type_a_ll = Kokkos::View<Scalar**, Kokkos::LayoutLeft, Device>;
  using view_type_b_ll = Kokkos::View<Scalar*, Kokkos::LayoutLeft, Device>;
  using view_type_c_ll = Kokkos::View<Scalar*, Kokkos::LayoutLeft, Device>;
  Test::impl_test_gemv_streams<execution_space, view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(space, mode, 0,
                                                                                                        1024);
  Test::impl_test_gemv_streams<execution_space, view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(space, mode, 13,
                                                                                                        1024);
  Test::impl_test_gemv_streams<execution_space, view_type_a_ll, view_type_b_ll, view_type_c_ll, Device>(space, mode, 50,
                                                                                                        40);
#endif

#if defined(KOKKOSKERNELS_INST_LAYOUTRIGHT)
  using view_type_a_lr = Kokkos::View<Scalar**, Kokkos::LayoutRight, Device>;
  using view_type_b_lr = Kokkos::View<Scalar*, Kokkos::LayoutRight, Device>;
  using view_type_c_lr = Kokkos::View<Scalar*, Kokkos::LayoutRight, Device>;
  Test::impl_test_gemv_streams<execution_space, view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(space, mode, 0,
                                                                                                        1024);
  Test::impl_test_gemv_streams<execution_space, view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(space, mode, 13,
                                                                                                        1024);
  Test::impl_test_gemv_streams<execution_space, view_type_a_lr, view_type_b_lr, view_type_c_lr, Device>(space, mode, 50,
                                                                                                        40);
#endif
  (void)space;
  return 1;
}

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)                          \
  TEST_F(TestCategory, blas##_##gemv_streams##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    test_gemv_streams<SCALAR, ORDINAL, OFFSET, DEVICE>("N");                                 \
    test_gemv_streams<SCALAR, ORDINAL, OFFSET, DEVICE>("T");                                 \
  }

#define NO_TEST_COMPLEX

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST
#undef NO_TEST_COMPLEX