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
#include <KokkosBlas1_rotmg.hpp>

namespace Test {
template <class View0, class PView, class RView>
void test_rotmg_impl(View0& d1, View0& d2, View0& x1, View0& y1, PView& param, RView& ref_vals) {
  using scalar_type = typename View0::non_const_value_type;
  using YView       = typename View0::const_type;

  YView y1_const(y1);

  KokkosBlas::rotmg(d1, d2, x1, y1_const, param);

  const scalar_type eps = Kokkos::ArithTraits<scalar_type>::eps();
  const scalar_type tol =
#if defined(KOKKOSKERNELS_ENABLE_TPL_BLAS) || defined(KOKKOSKERNELS_ENABLE_TPL_MKL)
      100 * eps;  // Guessing MKL implements sin/cos differently so need larger tol
#else
      10 * eps;
#endif
  auto d1_h = Kokkos::create_mirror_view(d1);
  Kokkos::deep_copy(d1_h, d1);
  EXPECT_NEAR_KK_REL(d1_h(), ref_vals(0), tol, "rotmg: d1 is off");

  auto d2_h = Kokkos::create_mirror_view(d2);
  Kokkos::deep_copy(d2_h, d2);
  EXPECT_NEAR_KK_REL(d2_h(), ref_vals(1), tol, "rotmg: d2 is off");

  auto x1_h = Kokkos::create_mirror_view(x1);
  Kokkos::deep_copy(x1_h, x1);
  EXPECT_NEAR_KK_REL(x1_h(), ref_vals(2), tol, "rotmg: x1 is off");

  auto y1_h = Kokkos::create_mirror_view(y1_const);
  Kokkos::deep_copy(y1_h, y1_const);
  EXPECT_NEAR_KK_REL(y1_h(), ref_vals(3), tol, "rotmg: y1 is off");

  auto param_h = Kokkos::create_mirror_view(param);
  Kokkos::deep_copy(param_h, param);
  EXPECT_NEAR_KK_REL(param_h(0), ref_vals(4), tol, "rotmg: param(0) is off");
  EXPECT_NEAR_KK_REL(param_h(1), ref_vals(5), tol, "rotmg: param(1) is off");
  EXPECT_NEAR_KK_REL(param_h(2), ref_vals(6), tol, "rotmg: param(2) is off");
  EXPECT_NEAR_KK_REL(param_h(3), ref_vals(7), tol, "rotmg: param(3) is off");
  EXPECT_NEAR_KK_REL(param_h(4), ref_vals(8), tol, "rotmg: param(4) is off");
}

template <class View0, class PView, class RView>
void set_rotmg_input_ref_vals(const int test_case, View0& d1, View0& d2, View0& x1, View0& y1, PView& param,
                              RView& ref_vals) {
  constexpr double gamma = 4096;
  Kokkos::deep_copy(param, 0.0);
  switch (test_case) {
    case 0:
      Kokkos::deep_copy(d1, 0.1);
      Kokkos::deep_copy(d2, 0.3);
      Kokkos::deep_copy(x1, 1.2);
      Kokkos::deep_copy(y1, 0.2);

      ref_vals(0) = 12.0 / 130.0;
      ref_vals(1) = 36.0 / 130.0;
      ref_vals(2) = 1.3;
      ref_vals(3) = 0.2;
      ref_vals(4) = 0.0;
      ref_vals(5) = 0.0;
      ref_vals(6) = -1.0 / 6.0;
      ref_vals(7) = 0.5;
      ref_vals(8) = 0.0;
      break;
    case 1:
      Kokkos::deep_copy(d1, 0.7);
      Kokkos::deep_copy(d2, 0.2);
      Kokkos::deep_copy(x1, 0.6);
      Kokkos::deep_copy(y1, 4.2);

      ref_vals(0) = 14.0 / 75.0;
      ref_vals(1) = 49.0 / 75.0;
      ref_vals(2) = 4.5;
      ref_vals(3) = 4.2;
      ref_vals(4) = 1.0;
      ref_vals(5) = 0.5;
      ref_vals(6) = 0.0;
      ref_vals(7) = 0.0;
      ref_vals(8) = 1.0 / 7.0;
      break;
    case 2:
      Kokkos::deep_copy(d1, 0.0);
      Kokkos::deep_copy(d2, 0.0);
      Kokkos::deep_copy(x1, 0.0);
      Kokkos::deep_copy(y1, 0.0);

      ref_vals(0) = 0.0;
      ref_vals(1) = 0.0;
      ref_vals(2) = 0.0;
      ref_vals(3) = 0.0;
      ref_vals(4) = -2.0;
      ref_vals(5) = 0.0;
      ref_vals(6) = 0.0;
      ref_vals(7) = 0.0;
      ref_vals(8) = 0.0;
      break;
    case 3:
      Kokkos::deep_copy(d1, 4.0);
      Kokkos::deep_copy(d2, -1.0);
      Kokkos::deep_copy(x1, 2.0);
      Kokkos::deep_copy(y1, 4.0);

      ref_vals(0) = 0.0;
      ref_vals(1) = 0.0;
      ref_vals(2) = 0.0;
      ref_vals(3) = 4.0;
      ref_vals(4) = -1.0;
      ref_vals(5) = 0.0;
      ref_vals(6) = 0.0;
      ref_vals(7) = 0.0;
      ref_vals(8) = 0.0;
      break;
    case 4:
      Kokkos::deep_copy(d1, 6.0e-10);
      Kokkos::deep_copy(d2, 2.0e-2);
      Kokkos::deep_copy(x1, 1.0e5);
      Kokkos::deep_copy(y1, 10.0);

      ref_vals(0) = 45.0e-11 * gamma * gamma;
      ref_vals(1) = 15.0e-3;
      ref_vals(2) = 4.0e5 / (3 * gamma);
      ref_vals(3) = 10.0;
      ref_vals(4) = -1.0;
      ref_vals(5) = 1.0 / gamma;
      ref_vals(6) = -1.0e-4;
      ref_vals(7) = 1.0e4 / (3 * gamma);
      ref_vals(8) = 1.0;
      break;
    case 5:
      Kokkos::deep_copy(d1, 4.0e10);
      Kokkos::deep_copy(d2, 2.0e-2);
      Kokkos::deep_copy(x1, 1.0e-5);
      Kokkos::deep_copy(y1, 10.0);

      ref_vals(0) = 4.0e10 / (1.5 * gamma * gamma);
      ref_vals(1) = 2.0e-2 / 1.5;
      ref_vals(2) = 6144.0e-5;
      ref_vals(3) = 10.0;
      ref_vals(4) = -1.0;
      ref_vals(5) = 4096.0;
      ref_vals(6) = -1.0e6;
      ref_vals(7) = 5.0e-7 * gamma;
      ref_vals(8) = 1.0;
      break;
    case 6:
      Kokkos::deep_copy(d1, 2.0e-10);
      Kokkos::deep_copy(d2, 4.0e-2);
      Kokkos::deep_copy(x1, 1.0e5);
      Kokkos::deep_copy(y1, 10.0);

      ref_vals(0) = 4.0 / 150.0;
      ref_vals(1) = (2.0e-10 / 1.5) * (gamma * gamma);
      ref_vals(2) = 15.0;
      ref_vals(3) = 10.0;
      ref_vals(4) = -1.0;
      ref_vals(5) = 5.0e-5;
      ref_vals(6) = -1.0 / gamma;
      ref_vals(7) = 1.0;
      ref_vals(8) = 1.0e4 / gamma;
      break;
    case 7:
      Kokkos::deep_copy(d1, 2.0e10);
      Kokkos::deep_copy(d2, 4.0e-2);
      Kokkos::deep_copy(x1, 1.0e-5);
      Kokkos::deep_copy(y1, 10.0);

      ref_vals(0) = 4.0 / 150.0;
      ref_vals(1) = 2.0e10 / (1.5 * gamma * gamma);
      ref_vals(2) = 15.0;
      ref_vals(3) = 10.0;
      ref_vals(4) = -1.0;
      ref_vals(5) = 5.0e5;
      ref_vals(6) = -4096.0;
      ref_vals(7) = 1.0;
      ref_vals(8) = 4096.0e-6;
      break;
    case 8:
      Kokkos::deep_copy(d1, 4.0);
      Kokkos::deep_copy(d2, -2.0);
      Kokkos::deep_copy(x1, 8.0);
      Kokkos::deep_copy(y1, 4.0);

      ref_vals(0) = 32.0 / 7.0;
      ref_vals(1) = -16.0 / 7.0;
      ref_vals(2) = 7.0;
      ref_vals(3) = 4.0;
      ref_vals(4) = 0.0;
      ref_vals(5) = 0.0;
      ref_vals(6) = -0.5;
      ref_vals(7) = -0.25;
      ref_vals(8) = 0.0;
      break;
    default: throw std::runtime_error("rotmg test: test case unrecognized!"); break;
  }
}
}  // namespace Test

template <class Scalar, class Device>
int test_rotmg() {
  Kokkos::View<Scalar, Device> d1("d1"), d2("d2"), x1("x1"), y1("y1");
  Kokkos::View<Scalar[5], Device> param("param");
  Kokkos::View<Scalar[9], Kokkos::DefaultHostExecutionSpace> ref_vals("reference values");

  constexpr int num_test_cases = 9;
  for (int test_case = 0; test_case < num_test_cases; ++test_case) {
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
    // There is a bug in the current cuBLAS implementation for this
    // corner case so we do not test it when cuBLAS is enabled.
    if (test_case == 4 || test_case == 5) {
      continue;
    }
#endif
    Test::set_rotmg_input_ref_vals(test_case, d1, d2, x1, y1, param, ref_vals);
    Test::test_rotmg_impl(d1, d2, x1, y1, param, ref_vals);
  }

  return 1;
}

#define KOKKOSKERNELS_EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE)          \
  TEST_F(TestCategory, rotmg##_##SCALAR##_##ORDINAL##_##OFFSET##_##DEVICE) { \
    Kokkos::Profiling::pushRegion("KokkosBlas::Test::rotg");                 \
    test_rotmg<SCALAR, DEVICE>();                                            \
    Kokkos::Profiling::popRegion();                                          \
  }

#define NO_TEST_COMPLEX

#include <Test_Common_Test_All_Type_Combos.hpp>

#undef KOKKOSKERNELS_EXECUTE_TEST
#undef NO_TEST_COMPLEX
