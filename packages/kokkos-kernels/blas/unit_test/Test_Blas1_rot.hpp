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
#include <KokkosBlas1_rot.hpp>

template <class Scalar, class ExecutionSpace>
int test_rot() {
  using mag_type        = typename Kokkos::ArithTraits<Scalar>::mag_type;
  using vector_type     = Kokkos::View<Scalar*, ExecutionSpace>;
  using scalar_type     = Kokkos::View<mag_type, ExecutionSpace>;
  using vector_ref_type = Kokkos::View<Scalar*, Kokkos::HostSpace>;

  vector_type X("X", 4), Y("Y", 4);
  vector_ref_type Xref("Xref", 4), Yref("Yref", 4);
  scalar_type c("c"), s("s");

  // Initialize inputs
  typename vector_type::HostMirror X_h = Kokkos::create_mirror_view(X);
  typename vector_type::HostMirror Y_h = Kokkos::create_mirror_view(Y);
  X_h(0)                               = 0.6;
  X_h(1)                               = 0.1;
  X_h(2)                               = -0.5;
  X_h(3)                               = 0.8;
  Y_h(0)                               = 0.5;
  Y_h(1)                               = -0.9;
  Y_h(2)                               = 0.3;
  Y_h(3)                               = 0.7;
  Kokkos::deep_copy(X, X_h);
  Kokkos::deep_copy(Y, Y_h);

  Kokkos::deep_copy(c, 0.8);
  Kokkos::deep_copy(s, 0.6);

  // Compute the rotated vectors
  KokkosBlas::rot(X, Y, c, s);
  Kokkos::fence();

  // Bring solution back to host
  Kokkos::deep_copy(X_h, X);
  Kokkos::deep_copy(Y_h, Y);

  // Check outputs against reference values
  Xref(0) = 0.78;
  Xref(1) = -0.46;
  Xref(2) = -0.22;
  Xref(3) = 1.06;
  Yref(0) = 0.04;
  Yref(1) = -0.78;
  Yref(2) = 0.54;
  Yref(3) = 0.08;

  Scalar const tol = 10 * Kokkos::ArithTraits<Scalar>::eps();
  for (int idx = 0; idx < 4; ++idx) {
    Test::EXPECT_NEAR_KK_REL(X_h(idx), Xref(idx), tol);
    Test::EXPECT_NEAR_KK_REL(Y_h(idx), Yref(idx), tol);
  }

  return 1;
}

#if defined(KOKKOSKERNELS_INST_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, rot_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::rot");
  test_rot<float, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, rot_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::rot");
  test_rot<double, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_FLOAT) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, rot_complex_float) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::rot");
  test_rot<Kokkos::complex<float>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif

#if defined(KOKKOSKERNELS_INST_COMPLEX_DOUBLE) || \
    (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
TEST_F(TestCategory, rot_complex_double) {
  Kokkos::Profiling::pushRegion("KokkosBlas::Test::rot");
  test_rot<Kokkos::complex<double>, TestDevice>();
  Kokkos::Profiling::popRegion();
}
#endif
