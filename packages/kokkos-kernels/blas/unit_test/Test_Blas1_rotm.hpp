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
#include <KokkosBlas1_rotm.hpp>

namespace Test {

template <class vector_view_type, class param_view_type, class vector_ref_type>
void set_rotm_inputs(const int &test_case, vector_view_type &X, vector_view_type &Y, param_view_type &param,
                     vector_ref_type &Xref, vector_ref_type &Yref) {
  // Initialize X and Y inputs
  typename vector_view_type::HostMirror X_h = Kokkos::create_mirror_view(X);
  typename vector_view_type::HostMirror Y_h = Kokkos::create_mirror_view(Y);

  X_h(0) = 0.6;
  X_h(1) = 0.1;
  X_h(2) = -0.5;
  X_h(3) = 0.8;
  Y_h(0) = 0.5;
  Y_h(1) = -0.9;
  Y_h(2) = 0.3;
  Y_h(3) = 0.7;
  Kokkos::deep_copy(X, X_h);
  Kokkos::deep_copy(Y, Y_h);

  // Initialize Xref, Yref and param (test case dependent)
  typename param_view_type::HostMirror param_h = Kokkos::create_mirror_view(param);
  switch (test_case) {
    case 0:
      param_h(0) = -2.0;
      param_h(1) = 0.0;
      param_h(2) = 0.0;
      param_h(3) = 0.0;
      param_h(4) = 0.0;

      Xref(0) = 0.60;
      Xref(1) = 0.10;
      Xref(2) = -0.50;
      Xref(3) = 0.80;
      Yref(0) = 0.50;
      Yref(1) = -0.90;
      Yref(2) = 0.30;
      Yref(3) = 0.70;
      break;

    case 1:
      param_h(0) = -1.0;
      param_h(1) = 2.0;
      param_h(2) = -3.0;
      param_h(3) = -4.0;
      param_h(4) = 5.0;

      Xref(0) = -0.80;
      Xref(1) = 3.80;
      Xref(2) = -2.20;
      Xref(3) = -1.20;
      Yref(0) = 0.70;
      Yref(1) = -4.80;
      Yref(2) = 3.00;
      Yref(3) = 1.10;
      break;

    case 2:
      param_h(0) = 0.0;
      param_h(1) = 0.0;
      param_h(2) = 2.0;
      param_h(3) = -3.0;
      param_h(4) = 0.0;

      Xref(0) = -0.90;
      Xref(1) = 2.80;
      Xref(2) = -1.40;
      Xref(3) = -1.30;
      Yref(0) = 1.70;
      Yref(1) = -0.70;
      Yref(2) = -0.70;
      Yref(3) = 2.30;
      break;

    case 3:
      param_h(0) = 1.0;
      param_h(1) = 5.0;
      param_h(2) = 2.0;
      param_h(3) = 0.0;
      param_h(4) = -4.0;

      Xref(0) = 3.50;
      Xref(1) = -0.40;
      Xref(2) = -2.20;
      Xref(3) = 4.70;
      Yref(0) = -2.60;
      Yref(1) = 3.50;
      Yref(2) = -0.70;
      Yref(3) = -3.60;
      break;
    default: throw std::runtime_error("rotm: unimplemented test case!");
  }

  Kokkos::deep_copy(param, param_h);

  return;
}

template <class vector_view_type, class vector_ref_type>
void check_results(vector_view_type &X, vector_view_type &Y, vector_ref_type &Xref, vector_ref_type &Yref) {
  using Scalar = typename vector_view_type::value_type;

  typename vector_view_type::HostMirror X_h = Kokkos::create_mirror_view(X);
  typename vector_view_type::HostMirror Y_h = Kokkos::create_mirror_view(Y);
  Kokkos::deep_copy(X_h, X);
  Kokkos::deep_copy(Y_h, Y);

  Scalar const tol = 10 * Kokkos::ArithTraits<Scalar>::eps();
  for (int idx = 0; idx < 4; ++idx) {
    Test::EXPECT_NEAR_KK_REL(X_h(idx), Xref(idx), tol);
    Test::EXPECT_NEAR_KK_REL(Y_h(idx), Yref(idx), tol);
  }

  return;
}

}  // namespace Test

template <class Scalar, class ExecutionSpace>
int test_rotm() {
  using vector_view_type = Kokkos::View<Scalar *, ExecutionSpace>;
  using vector_ref_type  = Kokkos::View<Scalar *, Kokkos::HostSpace>;
  using param_view_type  = Kokkos::View<Scalar[5], ExecutionSpace>;

  vector_view_type X("X", 4), Y("Y", 4);
  vector_ref_type Xref("Xref", 4), Yref("Yref", 4);
  param_view_type param("param");

  for (int test_case = 0; test_case < 4; ++test_case) {
    // Initialize inputs
    Test::set_rotm_inputs(test_case, X, Y, param, Xref, Yref);

    // Compute the rotated vectors
    KokkosBlas::rotm(X, Y, param);
    Kokkos::fence();

    // Check outputs
    Test::check_results(X, Y, Xref, Yref);
  }

  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, rotm_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::rotm");
  test_rotm<float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, rotm_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::rotm");
  test_rotm<double, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif
